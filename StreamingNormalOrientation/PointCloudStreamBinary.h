#pragma once

#include "OpenPointCloudStreamBinary.h"

#include <string>
#include <iostream>

#include <Eigen/Dense>
#include <omp.h>

using namespace std;

// Represents a point cloud stream in a binary file. This stream is not necessarily sorted.
// TVertex: Data type of vertices. Must have a Eigen::Vector3f position
template <typename TVertex>
class PointCloudStreamBinary
{
	const size_t blockSize = 40960000; // ~40 MB

public:

	PointCloudStreamBinary(const char* filename)
		: filename(std::string(filename))
	{ }

	// Returns the number of points in this file based on the file size.
	long long NumberOfPoints() const { return boost::filesystem::file_size(filename) / sizeof(TVertex); }
	
	// Asserts that the internal file handle exists. Throws a FileReadException if this is not the case.
	void AssertValid() const 
	{ 
		if (!boost::filesystem::exists(filename))
			throw FileReadException(filename);
	}

	// Reads the contents of the file and saves them in the provided array. The array must be pre-allocated.
	void CopyPointsToArray(TVertex* array)
	{
		FILE* f = fopen(filename.c_str(), "rb");
		if (!f)
			throw FileReadException(filename);
		fread(array, sizeof(TVertex), NumberOfPoints(), f);
		fclose(f);
	}

	OpenPointCloudStreamBinary<TVertex> OpenStream()
	{
		OpenPointCloudStreamBinary<TVertex> result(filename);
		return std::move(result);
	}
	
	// Calls a processor method on each vertex in the stream (packed in an array).
	// ProcessFunction: void(TVertex*, size_t vertices, int threadNum, void* parameter)
	// verticesPerBlockPerThread - Specifies the size of the arrays that are passed to the processor function.
	// vertexProcessor - The callback function that is used.
	// parameter - A user defined parameter that is passed to the processor.
	// parallel - Specifies if several processors can be called in parallel.
	template <typename ProcessFunction>
	void Process(int verticesPerBlockPerThread, ProcessFunction vertexProcessor, void* parameter, bool parallel)
	{
		OpenPointCloudStreamBinary<TVertex> stream = OpenStream();
		
		int threads = parallel ? omp_get_num_procs() : 1;

		CompoundPointCloudBuffer<TVertex> readBuffer(threads, verticesPerBlockPerThread);

		bool eof = false;

#pragma omp parallel num_threads(threads)
		{
			int threadNum = omp_get_thread_num();
			PointCloudBuffer<TVertex> partialBuffer = readBuffer[threadNum];

			while (!eof)
			{
				#pragma omp critical (updateEof)
				{
					stream >> partialBuffer;
					eof = partialBuffer.NumberOfPoints() < verticesPerBlockPerThread;
				}

				vertexProcessor(partialBuffer.Ptr(), partialBuffer.NumberOfPoints(), threadNum, parameter);
			}
		}
	}

	// First calls a preparation function for each thread. Then calls a processor method on each vertex in the stream (packed in an array).
	// PrepareFunction: void(int)
	// ProcessFunction: void(TVertex*, size_t vertices, int threadNum, void* parameter)
	// verticesPerBlockPerThread - Specifies the size of the arrays that are passed to the processor function.
	// vertexProcessor - The callback function that is used.
	// parameter - A user defined parameter that is passed to the processor.
	// parallel - Specifies if several processors can be called in parallel.
	template <typename PrepareFunction, typename ProcessFunction>
	void Process(int verticesPerBlockPerThread, PrepareFunction prepare, ProcessFunction vertexProcessor, void* parameter, bool parallel)
	{
		int threads = parallel ? omp_get_num_procs() : 1;
#pragma omp parallel num_threads(threads)
		{
			int threadNum = omp_get_thread_num();
			prepare(threadNum);
		}
		Process(verticesPerBlockPerThread, vertexProcessor, parameter, parallel);
	}
	
	// Prepares the point cloud for streaming, i.e. sorts the points along the x-axis.
	// realign - set to true to allow the point cloud to be rotated to map the longest principal axis to the x-axis.
	void PrepareForStreaming(const Eigen::Vector3f& bbxMin, const Eigen::Vector3f& bbxMax, const Eigen::Vector3f& centroid, bool realign = true)
	{
		int nProcessors = omp_get_num_procs();

		Eigen::Matrix3f rot = Eigen::Matrix3f::Identity();

		if (realign)
		{			
#ifndef ACCURATE_TIMING
			cout << "Calculating covariance..." << endl;
#endif
			CovarianceCalculator covCalc(nProcessors, centroid);

			Process(
				(int)(blockSize / sizeof(TVertex) / nProcessors),
				std::bind(&CovarianceCalculator::PrepareCovariance, &covCalc, placeholders::_1),
				std::bind(&CovarianceCalculator::AccumulateCovariance, &covCalc, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4), nullptr, true);

			//Accumulate the partial covariances
			covCalc.AccumulateCovariances();

#ifndef ACCURATE_TIMING
			cout << "Covariance: " << endl << covCalc.covariance[0] << endl;
#endif

			Eigen::JacobiSVD<Eigen::Matrix3f> svd(covCalc.covariance[0], Eigen::ComputeFullU);
			rot = svd.matrixU().transpose();

#ifndef ACCURATE_TIMING
			cout << "Rotation: " << endl << rot << endl;
#endif
		}

		ExternalMergeSort(centroid, rot);
	}

private:
	std::string filename;

	struct CentroidRotation{ const Eigen::Vector3f& centroid; const Eigen::Matrix3f rotation; };

	int sortChunks;
	const int sortVerticesPerChunk = 65536;
	const std::string sortChunkPrefix = "sort_chunk_";

	inline std::string GetSortChunkFilename(int mergeLevel, int chunkNo)
	{
		return (sortChunkPrefix + std::to_string(mergeLevel) + std::string("_") + std::to_string(chunkNo) + std::string(".bin"));
	}

	// Processor method that re-aligns the vertices with a given centroid and rotation. Then, sorts the vertices and writes them to a temporary output file.
	void BinVertices(TVertex* vertices, size_t nVertices, size_t threadNum, void* parameter)
	{
		CentroidRotation* cr = static_cast<CentroidRotation*>(parameter);

		int chunkNo;
#pragma omp critical (updateChunkCount)
		{
			chunkNo = sortChunks++;
		}

		//Align along longest axis
		for (int i = 0; i < nVertices; ++i)
		{
			vertices[i].Move(-cr->centroid);
			vertices[i].Transform(cr->rotation);
		}

		std::sort(vertices, vertices + nVertices, [](const TVertex& lhs, const TVertex& rhs) { return lhs.position[0] < rhs.position[0]; });

		FILE* file = fopen(GetSortChunkFilename(0, chunkNo).c_str(), "wb");
		fwrite(vertices, sizeof(TVertex), nVertices, file);
		fclose(file);
	}

	struct CompareMergePair
	{
		bool operator()(const std::pair<int, TVertex>& lhs, const std::pair<int, TVertex>& rhs)
		{
			return lhs.second.position[0] > rhs.second.position[0];
		}
	};

	void ExternalMergeSort(const Eigen::Vector3f& centroid, const Eigen::Matrix3f& rot)
	{
#ifndef ACCURATE_TIMING
		cout << "Sorting vertices..." << endl;
#endif

		const int chunksPerThread = 25;
		const int mergeBufferSizePerChunk = 1000;

#ifndef ACCURATE_TIMING
		cout << "Splitting vertices into chunks with " << sortVerticesPerChunk << " vertices..." << endl;
#endif

		CentroidRotation cr = { centroid, rot };

		//Split and sort
		sortChunks = 0;
		Process(
			sortVerticesPerChunk,
			std::bind(&PointCloudStreamBinary<TVertex>::BinVertices, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4), &cr, true);

		//Merge
		int mergeLevel = 1;

		int maxThreads = min((int)ceil((float)sortChunks / chunksPerThread), omp_get_num_procs());		

		int processedChunks = 0;
		int currentLevelChunks = sortChunks;
		int nextOutputNumber = 0;

		//global data buffers
		struct Buffers
		{
			TVertex* buffer;
			FILE** globalChunkInputFiles;
			TVertex** globalChunkBuffer;
			TVertex** globalReadBuffer;
			TVertex** globalChunkBufferEnd;
			bool* globalEof;

			~Buffers()
			{
				if (buffer) delete[] buffer;

				if (globalChunkInputFiles) delete[] globalChunkInputFiles;
				if (globalChunkBuffer) delete[] globalChunkBuffer;
				if (globalReadBuffer) delete[] globalReadBuffer;
				if (globalChunkBufferEnd) delete[] globalChunkBufferEnd;
				if (globalEof) delete[] globalEof;
			}
		} buffers;

		buffers.buffer = new TVertex[chunksPerThread * mergeBufferSizePerChunk * maxThreads];
		buffers.globalChunkInputFiles = new FILE*[maxThreads * chunksPerThread];
		buffers.globalChunkBuffer = new TVertex*[maxThreads * chunksPerThread];
		buffers.globalReadBuffer = new TVertex*[maxThreads * chunksPerThread];
		buffers.globalChunkBufferEnd = new TVertex*[maxThreads * chunksPerThread];
		buffers.globalEof = new bool[maxThreads * chunksPerThread];

		//Merge iterations, as long as data is cluttered across several files
		while (currentLevelChunks > 1)
		{
			int threads = min(maxThreads, (int)ceil((float)currentLevelChunks / chunksPerThread));

#ifndef ACCURATE_TIMING
			cout << "Merging " << currentLevelChunks << " chunks with " << threads << " threads..." << endl;
#endif

#pragma omp parallel num_threads(threads)
			{
				int threadNum = omp_get_thread_num();

				FILE** chunkInputFiles = buffers.globalChunkInputFiles + threadNum * chunksPerThread;
				TVertex** chunkBuffer = buffers.globalChunkBuffer + threadNum * chunksPerThread;
				TVertex** readBuffer = buffers.globalReadBuffer + threadNum * chunksPerThread;
				TVertex** chunkBufferEnd = buffers.globalChunkBufferEnd + threadNum * chunksPerThread;
				bool* eof = buffers.globalEof + threadNum * chunksPerThread;

				//as long as there are unprocessed chunks in the current level...
				while (processedChunks < currentLevelChunks)
				{
					int chunks, chunkStart;

#pragma omp critical (select_chunks)
					{
						chunkStart = processedChunks;
						chunks = min(chunksPerThread, currentLevelChunks - chunkStart);
						processedChunks += chunks;
					}

					if (chunks == 1)
					{
						//we don't need to merge. Just rename the file
#pragma omp critical
					{
						auto newFilename = GetSortChunkFilename(mergeLevel, nextOutputNumber++).c_str();

						//Remove the file if it already exists
						remove(newFilename);

						rename(GetSortChunkFilename(mergeLevel - 1, chunkStart).c_str(), newFilename);
					}
					break;
					}

					std::priority_queue<std::pair<int /*chunkNo*/, TVertex>, std::vector<std::pair<int, TVertex>>, CompareMergePair> queue;

					for (int i = 0; i < chunks; ++i)
					{
						chunkBuffer[i] = buffers.buffer + mergeBufferSizePerChunk * (i + chunksPerThread * threadNum);
						readBuffer[i] = chunkBuffer[i];
						chunkInputFiles[i] = fopen(GetSortChunkFilename(mergeLevel - 1, chunkStart + i).c_str(), "rb");

						size_t read = fread(chunkBuffer[i], sizeof(TVertex), mergeBufferSizePerChunk, chunkInputFiles[i]);
						chunkBufferEnd[i] = chunkBuffer[i] + read;
						eof[i] = read < mergeBufferSizePerChunk;

						if (readBuffer[i] < chunkBufferEnd[i])
							queue.push(std::make_pair(i, *readBuffer[i]));
						readBuffer[i]++;
					}

					FILE* outputFile;
#pragma omp critical
					{
						outputFile = fopen(GetSortChunkFilename(mergeLevel, nextOutputNumber++).c_str(), "wb");
					}

					while (!queue.empty())
					{
						auto nextElement = queue.top();
						queue.pop();

						fwrite(&nextElement.second, sizeof(TVertex), 1, outputFile);

						int chunkNo = nextElement.first;
						bool chunkHasMoreElements = true;
						if (readBuffer[chunkNo] >= chunkBufferEnd[chunkNo])
						{
							//read new data from file
							if (!eof[chunkNo])
							{
								readBuffer[chunkNo] = chunkBuffer[chunkNo];
								int read = fread(chunkBuffer[chunkNo], sizeof(TVertex), mergeBufferSizePerChunk, chunkInputFiles[chunkNo]);
								chunkBufferEnd[chunkNo] = chunkBuffer[chunkNo] + read;
								eof[chunkNo] = read < mergeBufferSizePerChunk;
							}
						}
						if (readBuffer[chunkNo] < chunkBufferEnd[chunkNo])
						{
							queue.push(std::make_pair(chunkNo, *(readBuffer[chunkNo])));
							(readBuffer[chunkNo])++;
						}
					}

					for (int i = 0; i < chunks; ++i)
					{
						fclose(chunkInputFiles[i]);
						remove(GetSortChunkFilename(mergeLevel - 1, chunkStart + i).c_str());
					}
					fclose(outputFile);

				} //while process chunks in current level
			} //parallel block
			mergeLevel++;

			currentLevelChunks = nextOutputNumber;
			processedChunks = 0;
			nextOutputNumber = 0;
		} //for each level

		//Rename the result of the last merge operation
		remove(filename.c_str());
		rename(GetSortChunkFilename(mergeLevel - 1, 0).c_str(), filename.c_str());			
	}

	struct CovarianceCalculator
	{
		Eigen::Matrix3f* covariance;
		double* covarianceNumPoints;
		Eigen::Vector3f centroid;
		int nThreads;

		CovarianceCalculator(int nThreads, Eigen::Vector3f centroid)
			: nThreads(nThreads)
		{
			covariance = new Eigen::Matrix3f[nThreads];
			covarianceNumPoints = new double[nThreads];
			this->centroid = centroid;
		}

		~CovarianceCalculator()
		{
			delete[] covariance;
			delete[] covarianceNumPoints;
		}

		void PrepareCovariance(int threadNum)
		{
			covariance[threadNum] = Eigen::Matrix3f::Zero();
			covarianceNumPoints[threadNum] = 0;
		}

		void AccumulateCovariance(TVertex* vertices, size_t nVertices, int threadNum, void* parameter)
		{
			for (int i = 0; i < nVertices; ++i)
			{
				covarianceNumPoints[threadNum]++;
				Eigen::Vector3f p = vertices[i].position - centroid;
				covariance[threadNum] +=
					(float)(1.0 / covarianceNumPoints[threadNum]) * (p * p.transpose() - covariance[threadNum]); //the numerically stable version
				// p * p.transpose(); //the default version
			}
		}

		void AccumulateCovariances()
		{
			for (int i = 1; i < nThreads; ++i)
				covariance[0] += covariance[i];
		}
	};
};