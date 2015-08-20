#pragma once

#include "PointCloudStreamBinary.h"
#include "exceptions.h"

#include <boost/filesystem.hpp>

// Represents a point cloud stream in a set of ASCII files.
// TVertex: Data type of vertices. Must have a ReadFromString(&readPointer) method.
template <typename TVertex>
class PointCloudStreamAscii
{
	const size_t blockSize = 40960000; // ~40 MB

public:	
	static std::shared_ptr<PointCloudStreamAscii<TVertex>> FromFolder(const char* foldername)
	{
		std::shared_ptr<PointCloudStreamAscii<TVertex>> cloud = std::make_shared < PointCloudStreamAscii<TVertex>>();

		using namespace boost::filesystem;
		for (directory_iterator dirIt(foldername); dirIt != directory_iterator(); ++dirIt)
		{
			if (dirIt->path().extension() == ".xyz")
				cloud->AddFile(dirIt->path().string().c_str());
		}		

		return cloud;
	}

	static std::shared_ptr<PointCloudStreamAscii<TVertex>> FromFile(const char* filename)
	{
		std::shared_ptr<PointCloudStreamAscii<TVertex>> cloud = std::make_shared < PointCloudStreamAscii<TVertex>>();
		cloud->AddFile(filename);		
		return cloud;
	}

	void AddFile(const char* filename)
	{
		//check if files exist
		if (!boost::filesystem::exists(filename))
			throw FileReadException(std::string(filename));
		filenames.push_back(std::string(filename));
	}

	// Reads the referenced ASCII files, parses the entries and writes them to the specified binary file.
	// ProcessFunction: void (const TVertex&)
	// filenameOut - path of the output file
	// vertexProcessor - this callback function is called on every read vertex
	template <typename ProcessFunction>
	std::shared_ptr<PointCloudStreamBinary<TVertex>> ConvertToBinary(const char* filenameOut, ProcessFunction vertexProcessor)
	{
		FILE* fileOut = fopen(filenameOut, "wb");

		int maxFileSize = 0;
		unsigned long long sumFileSize = 0;

		for (auto it = filenames.begin(); it != filenames.end(); ++it)
		{
			uintmax_t fileSize = boost::filesystem::file_size(*it);

			maxFileSize = max<long>(maxFileSize, fileSize);
			sumFileSize = sumFileSize + fileSize;
		}

		int readBufferSize = min<long>(maxFileSize, blockSize);
		char* readBuffer = new char[readBufferSize + 1];
		readBuffer[readBufferSize] = 0;

		unsigned long long processed = 0;
		auto prevPrecision = cout.precision(1);
		cout.setf(std::ios::fixed, std::ios::floatfield);

		for (auto it = filenames.begin(); it != filenames.end(); ++it)
		{
			FILE* fileIn = fopen(it->c_str(), "rb");

			// Get the file size
			uintmax_t fileSize = boost::filesystem::file_size(*it);

			if (fileIn)
			{				
				char* readPointer = readBuffer;
				char* readBufferEnd = readBuffer + readBufferSize;

				long read = 0;

				fread(readBuffer, 1, readBufferSize, fileIn);
				read += readBufferSize;

				TVertex v;

				while (read < fileSize || readPointer < readBufferEnd)
				{
					//look for next newline
					char* newlinePtr = find(readPointer, readBufferEnd, '\n');

					if (newlinePtr == readBufferEnd) //if there is no newline in the read buffer
					{
						if (read < fileSize) //If there is still something left to read.
						{
							//copy the remaining part to the beginning of the read buffer
							size_t unparsedSize = readBufferEnd - readPointer;
							memcpy(readBuffer, readPointer, unparsedSize);
							readPointer = readBuffer;
							//read the next chunk from file
							size_t readSize = fread(readBuffer + unparsedSize, 1, readBufferSize - unparsedSize, fileIn);
							read += (long)readSize;
							readBufferEnd = readBuffer + readSize + unparsedSize;
							if (readBufferEnd - readBuffer < readBufferSize)
								*readBufferEnd = 0; //string end, so strtof does not read previous data
							newlinePtr = find(readPointer, readBufferEnd, '\n');
						}
					}

					v.ReadFromString(&readPointer);
					fwrite(&v, sizeof(TVertex), 1, fileOut);

					vertexProcessor(v);

					readPointer = newlinePtr + 1;
				}

				fclose(fileIn);

				processed += fileSize;
				std::cout << "\r" << (100.0f * processed / sumFileSize) << " % complete...";
			}
			else
			{
				throw FileReadException(*it);
			}
		}

		std::cout << std::endl;

		cout.unsetf(std::ios::floatfield);
		cout.precision(prevPrecision);

		delete[] readBuffer;

		fclose(fileOut);

		return std::make_shared<PointCloudStreamBinary<TVertex>>(filenameOut);
	}

	std::shared_ptr<PointCloudStreamBinary<TVertex>> PrepareForStreaming(const char* binaryFilename, bool realign = true)
	{
		BbxCentroidAccumulator acc;

		std::shared_ptr<PointCloudStreamBinary<TVertex>> binary = ConvertToBinary(binaryFilename,
			std::bind(&BbxCentroidAccumulator::AccumulateBbxAndCentroidAscii, &acc, placeholders::_1));
		
		binary->PrepareForStreaming(acc.bbxMin, acc.bbxMax, acc.centroid, realign);

		return binary;
	}

private:
	std::vector<std::string> filenames;

	struct BbxCentroidAccumulator
	{
		Eigen::Vector3f bbxMin, bbxMax, centroid;
		double nPoints;

		BbxCentroidAccumulator()
			: nPoints(0)
		{
			bbxMin = Eigen::Vector3f(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity());
			bbxMax = Eigen::Vector3f(-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());
			centroid = Eigen::Vector3f::Zero();
		}

		void AccumulateBbxAndCentroidAscii(const TVertex& v)
		{
			++nPoints;
			double nPointsInv = 1.0 / nPoints;
			for (int i = 0; i < 3; ++i)
			{
				if (v.position[i] < bbxMin[i])
					bbxMin[i] = v.position[i];
				if (v.position[i] > bbxMax[i])
					bbxMax[i] = v.position[i];

				centroid[i] += (float)(nPointsInv * (v.position[i] - centroid[i]));
			}
		}
	};
};