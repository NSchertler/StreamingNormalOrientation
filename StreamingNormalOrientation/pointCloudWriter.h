#include <stdio.h>
#include <fstream>
#include <boost/filesystem.hpp>

using namespace std;

class PointCloudWriter
{
	template <typename TVertex>
	void writeToXYZ(TVertex* vertices, size_t nVertices, int threadNum, void*, ofstream& f)
	{
		for (int i = 0; i < nVertices; ++i)
			f << vertices[i] << endl;
	}

public:
	template <typename TVertex>
	void writeXYZ(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut)
	{
		std::ofstream fileStream(filenameOut);

		cloud->Process(1000, std::bind(&PointCloudWriter::writeToXYZ<TVertex>, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4, std::ref(fileStream)), nullptr, false);
		fileStream.close();
	}

	template <typename TVertex>
	void writePLYAscii(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut)
	{
		std::streamoff numElements = cloud->NumberOfPoints();

		std::ofstream fileStream(filenameOut);
		fileStream
			<< "ply" << endl
			<< "format ascii 1.0" << endl
			<< "element vertex " << numElements << endl;
		TVertex::WritePLYProperties(fileStream);
		fileStream
			<< "end_header" << endl;

		cloud->Process(1000, std::bind(&PointCloudWriter::writeToXYZ<TVertex>, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4, std::ref(fileStream)), nullptr, false);
		fileStream.close();
	}
};