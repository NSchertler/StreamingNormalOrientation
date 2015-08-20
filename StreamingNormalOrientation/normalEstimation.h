#pragma once

#include "normalEstimation.cuh"
#include "MortonUtils.cuh"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <assert.h>

#include "common.h"

#include <omp.h>
#include <algorithm>
#include <iostream>

inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
	if (result != cudaSuccess) {
		std::stringstream ss;
		ss << "CUDA Runtime Error: " << cudaGetErrorString(result);
		throw std::runtime_error(ss.str());
	}
#endif
	return result;
}

class NormalEstimation
{
	enum PipelineStream
	{
		First,
		Middle,
		Last
	};
	enum PipelineState
	{
		CopyH2D,
		IssueKernel,
		ProcessResults
	};

	struct Permutation
	{
		unsigned int sortKey;
		unsigned int originalIndex;
	};

public:
	// CUDA implementation of PCA normal estimation with simplified neighbor search
	template <typename TVertexIn, typename TVertexOut>
	std::shared_ptr<PointCloudStreamBinary<TVertexOut>> estimateNormals(std::shared_ptr<PointCloudStreamBinary<TVertexIn>> cloudIn, const char* filenameOut, float searchRadius, unsigned int maxVerticesPerSlice, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{
		const int minVerticesPerSlice = 2048;
		if (maxVerticesPerSlice < minVerticesPerSlice)
			maxVerticesPerSlice = minVerticesPerSlice;
		const int slicesPerRead = 12;

		const float quadTreeCellSize = 2.0f * searchRadius;

		float sliceWidth = 2.0f * searchRadius;
		float searchRadiusSqr = searchRadius * searchRadius;

		assert(slicesPerRead >= 4);

		auto streamIn = cloudIn->OpenStream();
		auto streamOut = OpenPointCloudStreamBinary<TVertexOut>(std::string(filenameOut), "wb");

		CompoundPointCloudBuffer<TVertexOut> vertexStore(4, maxVerticesPerSlice);
				
		unsigned int readBufferSize = maxVerticesPerSlice * slicesPerRead;
		ContiguousPointCloudBuffer<TVertexIn> bufferIn(readBufferSize, streamIn);

		struct LocalData
		{
			VertexPositionMortonIndex* h_transferMemoryUpload;
			float3* h_transferMemoryDownload;

			//Holds the order of incoming vertices.
			//Buffer is valid for all streams. Each stream has its own partial buffer with a stride of verticesPerSlice
			Permutation* permutationBuffer;
			
			VertexPositionMortonIndex* d_vertexBuffer[4];
			float3* d_normalBuffer[4];
			cudaEvent_t copyCompleteEvent[4];
			cudaEvent_t kernelCompleteEvent[4];

			cudaStream_t stream[4];

			~LocalData()
			{
				delete[] permutationBuffer;
				cudaFreeHost(h_transferMemoryUpload);
				cudaFreeHost(h_transferMemoryDownload);

				for (int i = 0; i < 4; ++i)
				{
					checkCuda(cudaStreamDestroy(stream[i]));

					checkCuda(cudaFree(d_vertexBuffer[i]));
					checkCuda(cudaFree(d_normalBuffer[i]));

					checkCuda(cudaEventDestroy(copyCompleteEvent[i]));
					checkCuda(cudaEventDestroy(kernelCompleteEvent[i]));
				}
			}
		} data;
		
				
		checkCuda(cudaMallocHost(&data.h_transferMemoryUpload, maxVerticesPerSlice * 4 * sizeof(VertexPositionMortonIndex)));
		checkCuda(cudaMallocHost(&data.h_transferMemoryDownload, maxVerticesPerSlice * 4 * sizeof(float3)));
					
		data.permutationBuffer = new Permutation[4 * maxVerticesPerSlice];

		//Allocate memory and GPU resources
		for (int i = 0; i < 4; ++i)
		{
			checkCuda(cudaStreamCreate(&data.stream[i]));

			checkCuda(cudaMalloc(&data.d_vertexBuffer[i], sizeof(VertexPositionMortonIndex) * maxVerticesPerSlice));
			checkCuda(cudaMalloc(&data.d_normalBuffer[i], sizeof(float3) * maxVerticesPerSlice));

			checkCuda(cudaEventCreateWithFlags(&data.copyCompleteEvent[i], cudaEventDisableTiming));
			checkCuda(cudaEventCreateWithFlags(&data.kernelCompleteEvent[i], cudaEventDisableTiming));
		}

		PipelineStream currentStream = First;
		PipelineState currentState = CopyH2D;

		int numVertices[4] = { -1, -1, -1, -1 };
		TVertexIn* workingSet[4];

		bool eof = false;
		int finishAfterFinalProcessingOfStream = -1;
		bool finished = false;

		while (!finished)
		{
			int streamFrom, streamTo;
			switch (currentStream)
			{
			case First:
				streamFrom = streamTo = 0; break;
			case Middle:
				streamFrom = 1; streamTo = 2; break;
			case Last:
				streamFrom = streamTo = 3; break;
			};

			switch (currentState)
			{
			case CopyH2D:
			{
				if (eof)
				{
					//only record necessary events
					for (int i = streamFrom; i <= streamTo; ++i)
						checkCuda(cudaEventRecord(data.copyCompleteEvent[i], data.stream[i]));
				}
				else
				{
					//Build working set
					for (int i = streamFrom; i <= streamTo; ++i)
					{
						//Make sure that the read buffer contains more data
						if (bufferIn.NumberOfUnprocessedPoints() == 0)
						{
							bufferIn.Fill();
							if (bufferIn.NumberOfUnprocessedPoints() == 0)
							{
								eof = true;
								finishAfterFinalProcessingOfStream = i == 0 ? 3 : i - 1;
							}
						}
							
						if (!eof)
						{
							float firstX = bufferIn.front().position.x();
							if (bufferIn.back().position.x() - firstX < sliceWidth || bufferIn.NumberOfUnprocessedPoints() < minVerticesPerSlice)
							{
								//we need to read more data.
								bufferIn.Fill();									
							}
							//search for slice size
							TVertexIn* it;
							for (it = bufferIn.begin(); it != bufferIn.end(); ++it)
							{
								if (it + 1 - bufferIn.begin() >= maxVerticesPerSlice || it->position.x() - firstX >= sliceWidth)
								{
									++it;
									break;
								}
							}

							//don't create too small slices
							if (it - bufferIn.begin() < minVerticesPerSlice)
							{
								it = std::min(bufferIn.begin() + minVerticesPerSlice, bufferIn.end());
							}

							numVertices[i] = it - bufferIn.begin();
							workingSet[i] = bufferIn.begin();
							bufferIn.SetProcessedUpTo(it);
						}
						else
						{
							numVertices[i] = 0;
						}


						//Copy working set to transfer and output buffer
						auto vStore = vertexStore[i];

						//calculate Morton index and sort permutation buffer
#pragma omp parallel for						
						for (int j = 0; j < numVertices[i]; ++j)
						{
							data.permutationBuffer[j + maxVerticesPerSlice * i].originalIndex = j;
							data.permutationBuffer[j + maxVerticesPerSlice * i].sortKey
								= getBinnedMortonIndex(workingSet[i][j].position.y(), workingSet[i][j].position.z(), bbxMin.y(), bbxMin.z(), quadTreeCellSize);
						}

						std::sort(data.permutationBuffer + maxVerticesPerSlice * i, data.permutationBuffer + maxVerticesPerSlice * i + numVertices[i],
							[](const Permutation& p1, const Permutation& p2) { return p1.sortKey < p2.sortKey; });

#pragma omp parallel for						
						for (int j = 0; j < numVertices[i]; ++j)
						{
							//copy to vertexStore in the original order
							//TODO: use custom assignment operator for both vertex types
							memcpy(&vStore[j].position, &workingSet[i][j].position, sizeof(float3));

							//copy to transfer memory in custom order
							Permutation& permutation = data.permutationBuffer[j + maxVerticesPerSlice * i];
							memcpy(&data.h_transferMemoryUpload[j + maxVerticesPerSlice * i].position, &workingSet[i][permutation.originalIndex].position, sizeof(float3));
							data.h_transferMemoryUpload[j + maxVerticesPerSlice * i].index = permutation.sortKey;
						}

					} //for each stream

					//Upload to GPU
					for (int i = streamFrom; i <= streamTo; ++i)
					{
						checkCuda(cudaStreamWaitEvent(data.stream[i], data.kernelCompleteEvent[(i + 1) % 4], 0));
						checkCuda(cudaStreamWaitEvent(data.stream[i], data.kernelCompleteEvent[i == 0 ? 3 : i - 1], 0));

						checkCuda(cudaMemcpyAsync(data.d_vertexBuffer[i], data.h_transferMemoryUpload + maxVerticesPerSlice * i, sizeof(VertexPositionMortonIndex) * numVertices[i], cudaMemcpyHostToDevice, data.stream[i]));

						checkCuda(cudaEventRecord(data.copyCompleteEvent[i], data.stream[i]));
					}
				} //if eof, else

				//Update state
				currentState = IssueKernel;
				switch (currentStream)
				{
				case First: currentStream = Last; break;
				case Middle: currentStream = First; break;
				case Last: currentStream = Middle; break;
				}
			} //case CopyH2D
				break;
			case IssueKernel:
			{
				for (int i = streamFrom; i <= streamTo; ++i)
				{
					cudaStreamWaitEvent(data.stream[i], data.copyCompleteEvent[i == 0 ? 3 : i - 1], 0);
					cudaStreamWaitEvent(data.stream[i], data.copyCompleteEvent[(i + 1) % 4], 0);

					if (numVertices[i] > 0)
						calculateNormals(data.d_vertexBuffer[i], data.d_normalBuffer[i], //data buffers
						searchRadiusSqr,
						vertexStore[i].front().position.x(), //xmin, first vertex in original order
						vertexStore[i][numVertices[i] - 1].position.x(), //xmax, last vertex in original order
						bbxMin.y(), bbxMin.z(), quadTreeCellSize, numVertices[i], data.stream[i]);

					checkCuda(cudaEventRecord(data.kernelCompleteEvent[i], data.stream[i]));

					if (numVertices[i] > 0)
						checkCuda(cudaMemcpyAsync(data.h_transferMemoryDownload + i * maxVerticesPerSlice, data.d_normalBuffer[i], sizeof(float3) * numVertices[i], cudaMemcpyDeviceToHost, data.stream[i]));
				}

				currentState = ProcessResults;
				switch (currentStream)
				{
				case First: currentStream = Last; break;
				case Middle: currentStream = First; break;
				case Last: currentStream = Middle; break;
				}

			} // case IssueKernel
				break;

			case ProcessResults:
			{
				//Write to file
				for (int i = streamFrom; i <= streamTo; ++i)
				{
					auto success = cudaStreamSynchronize(data.stream[i]);
					if (success != cudaSuccess)
					{
						std::stringstream ss;
						ss << "Synchronization for stream " << i << " failed: " << success;
						throw std::runtime_error(ss.str());
					}

					if (numVertices[i] > 0)
					{
						auto vStore = vertexStore[i];

						//vertex store is sorted in original order
						//transfer memory is sorted in custom order (according to permutation buffer)

						for (int j = 0; j < numVertices[i]; ++j)
							memcpy(&vStore[data.permutationBuffer[j + i * maxVerticesPerSlice].originalIndex].normal, &data.h_transferMemoryDownload[j + i * maxVerticesPerSlice], sizeof(float3));

						streamOut.Write(vertexStore[i], numVertices[i]);
					}

					if (finishAfterFinalProcessingOfStream == i)
						finished = true;
				}

				currentState = CopyH2D;

			} //case ProcessResults
				break;

			} //switch


		} //while true		

		return std::make_shared<PointCloudStreamBinary<TVertexOut>>(filenameOut);
	}
};