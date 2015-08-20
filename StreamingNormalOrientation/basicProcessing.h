#include <mutex>

class BasicProcessing
{
public:
	Eigen::Vector3f bbxMin;
	Eigen::Vector3f bbxMax;
	Eigen::Vector3f centroid;
	double nPoints;

private:

	std::mutex bbxMutex;

	/// TVertex must have a Eigen::Vector3f position
	template <typename TVertex>
	void AccumulateBbxCentroid(TVertex* vertices, size_t nVertices, int threadNum, void*)
	{
		auto localBbxMin = Eigen::Vector3f(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity());
		auto localBbxMax = Eigen::Vector3f(-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());
		auto localCentroidSum = Eigen::Vector3f(0, 0, 0);

		for (int j = 0; j < nVertices; ++j)
		{
			for (int i = 0; i < 3; ++i)
			{
				if (vertices[j].position[i] < localBbxMin[i])
					localBbxMin[i] = vertices[j].position[i];
				if (vertices[j].position[i] > localBbxMax[i])
					localBbxMax[i] = vertices[j].position[i];
				localCentroidSum[i] += vertices[j].position[i];
			}
		}

		std::lock_guard<std::mutex> lock(bbxMutex);
		for (int i = 0; i < 3; ++i)
		{
			if (localBbxMin[i] < bbxMin[i])
				bbxMin[i] = localBbxMin[i];
			if (localBbxMax[i] > bbxMax[i])
				bbxMax[i] = localBbxMax[i];

			nPoints += nVertices;
			centroid += ((localCentroidSum - (float)nVertices * centroid) / nPoints);
		}
	}

	float sliceStart;
	int sliceVertexNumStart;
	int processedVertices;
	float sliceWidth;
	int sliceSize;

	template <typename TVertex>
	void findSliceSize(TVertex* vertices, size_t nVertices, int threadNum, void*)
	{
		if (sliceStart == std::numeric_limits<float>::infinity())
			sliceStart = vertices[0].position.x();
		for (int i = 0; i < nVertices; ++i)
		{
			if (vertices[i].position.x() - sliceStart >= sliceWidth)
			{
				sliceSize = max(sliceSize, i + processedVertices - sliceVertexNumStart);
				sliceStart = vertices[i].position.x();
				sliceVertexNumStart = i + processedVertices;
			}
		}
		processedVertices += nVertices;
	}	

	Eigen::Matrix3f rot;	

public:

	/**
	Updates bbxMin, bbxMax, centroid, nPoints
	TVertex must have a Eigen::Vector3f position
	*/
	template <typename TVertex>
	void calculateBoundingBoxCentroid(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud)
	{
		calculateBoundingBoxCentroid(cloud.get());
	}

	/**
	Updates bbxMin, bbxMax, centroid, nPoints
	TVertex must have a Eigen::Vector3f position
	*/
	template <typename TVertex>
	void calculateBoundingBoxCentroid(PointCloudStreamBinary<TVertex>* cloud)
	{
		bbxMin = Eigen::Vector3f(std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity());
		bbxMax = Eigen::Vector3f(-std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity());
		centroid = Eigen::Vector3f(0, 0, 0);
		nPoints = 0;

		cloud->Process(1024, std::bind(&BasicProcessing::AccumulateBbxCentroid < TVertex >, this, placeholders::_1,
			placeholders::_2, placeholders::_3, placeholders::_4), nullptr, true);
	}
		
	
	/**
	Calculates the number of vertices per slice for a given slice width.
	*/
	template <typename TVertex>
	unsigned int findSliceSizeForWidth(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, float sliceWidth)
	{
		sliceStart = std::numeric_limits<float>::infinity();
		sliceVertexNumStart = 0;
		processedVertices = 0;
		this->sliceWidth = sliceWidth;
		sliceSize = 0;

		cloud->Process(16384, std::bind(&BasicProcessing::findSliceSize < TVertex >, this, placeholders::_1, placeholders::_2, placeholders::_3, placeholders::_4), nullptr, false);

		if (sliceSize == 0)
			return cloud->NumberOfPoints();
		else
			return sliceSize;
	}

	//ref: https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
	template<typename T>
	struct matrix_hash : std::unary_function<T, size_t> {
		std::size_t operator()(T const& matrix) const {
			// Note that it is oblivious to the storage order of Eigen matrix (column- or
			// row-major). It will give you the same hash value for two different matrices if they
			// are the transpose of each other in different storage order.
			size_t seed = 0;
			for (size_t i = 0; i < matrix.size(); ++i) {
				auto elem = *(matrix.data() + i);
				seed ^= std::hash<typename T::Scalar>()(elem)+0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}
	};


	template<typename TVertex1, typename TVertex2>
	void compareNormals(std::shared_ptr<PointCloudStreamBinary<TVertex1>> cloud1, std::shared_ptr<PointCloudStreamBinary<TVertex2>> cloud2, int& correctNormals, int& wrongNormals, int& skippedPoints1, int& skippedPoints2)
	{
		OpenPointCloudStreamBinary<TVertex1> stream1 = cloud1->OpenStream();
		OpenPointCloudStreamBinary<TVertex2> stream2 = cloud2->OpenStream();

		ContiguousPointCloudBuffer<TVertex1> buffer1(1024, stream1);
		ContiguousPointCloudBuffer<TVertex2> buffer2(1024, stream2);

		correctNormals = 0;
		wrongNormals = 0;
		skippedPoints1 = 0;
		skippedPoints2 = 0;

		std::unordered_map<Eigen::Vector3f, TVertex1, matrix_hash<Eigen::Vector3f>> cache1;
		std::unordered_map<Eigen::Vector3f, TVertex2, matrix_hash<Eigen::Vector3f>> cache2;

		while (true)
		{
			if (buffer1.NumberOfUnprocessedPoints() == 0)
				buffer1.Fill();
			if (buffer2.NumberOfUnprocessedPoints() == 0)
				buffer2.Fill();

			if (buffer1.NumberOfUnprocessedPoints() == 0 && buffer2.NumberOfUnprocessedPoints() == 0)
				break;

			if (buffer1.NumberOfUnprocessedPoints() > 0 && buffer2.NumberOfUnprocessedPoints() > 0 && buffer1.front().position == buffer2.front().position)
			{
				if (buffer1.front().normal.dot(buffer2.front().normal) > 0)
					++correctNormals;
				else
					++wrongNormals;
				buffer1.SetProcessedUpTo(buffer1.begin() + 1);
				buffer2.SetProcessedUpTo(buffer2.begin() + 1);
			}
			else
			{
				//Either one buffer has no more points or the front points of the buffers are different
				//Use one point from a buffer and one point from a cache.
				ContiguousPointCloudBuffer<TVertex1>* useBuffer;
				std::unordered_map<Eigen::Vector3f, TVertex1, matrix_hash<Eigen::Vector3f>>* useReadCache;
				std::unordered_map<Eigen::Vector3f, TVertex1, matrix_hash<Eigen::Vector3f>>* useWriteCache;
				if (buffer1.NumberOfUnprocessedPoints() == 0)
				{
					useBuffer = &buffer2;
					useReadCache = &cache1;
					useWriteCache = &cache2;
				}
				else if (buffer2.NumberOfUnprocessedPoints() == 0)
				{
					useBuffer = &buffer1;
					useReadCache = &cache2;
					useWriteCache = &cache1;
				}
				else if (buffer1.front().position.x() < buffer2.front().position.x())
				{
					useBuffer = &buffer1;
					useReadCache = &cache2;
					useWriteCache = &cache1;
				}
				else
				{
					useBuffer = &buffer2;
					useReadCache = &cache1;
					useWriteCache = &cache2;
				}


				auto it = useReadCache->find(useBuffer->front().position);
				if (it == useReadCache->end())
				{
					//The point is not in the cache. Write it to its own cache.
					(*useWriteCache)[useBuffer->front().position] = useBuffer->front();
				}
				else
				{
					if (useBuffer->front().normal.dot(it->second.normal) > 0)
						++correctNormals;
					else
						++wrongNormals;
					useReadCache->erase(it);
				}
				useBuffer->SetProcessedUpTo(useBuffer->begin() + 1);

			}

		}
		skippedPoints1 += (int)cache1.size();
		skippedPoints2 += (int)cache2.size();
	}

	template <typename TVertex>
	std::shared_ptr<PointCloudStreamBinary<TVertex>> filter(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, const std::function<bool(const TVertex&)>& functor)
	{
		FILE* f = fopen(filenameOut, "wb");
		cloud->Process(4096, [&f, &functor](TVertex* v, size_t nVertices, int thread, void* param)
		{
			for (int i = 0; i < nVertices; ++i)
				if (functor(v[i]))
					fwrite(v + i, sizeof(TVertex), 1, f);
		}, nullptr, false);
		fclose(f);
		return std::make_shared<PointCloudStreamBinary<TVertex>>(filenameOut);
	}

	template <typename TVertex>
	std::shared_ptr<PointCloudStreamBinary<TVertex>> map(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, const std::function<void(TVertex&)>& functor)
	{
		FILE* f = fopen(filenameOut, "wb");
		cloud->Process(4096, [&f, &functor](TVertex* v, size_t nVertices, int thread, void* param)
		{
			for (int i = 0; i < nVertices; ++i)
			{
				functor(v[i]);
			}
			fwrite(v, sizeof(TVertex), nVertices, f);
		}, nullptr, false);
		fclose(f);
		return std::make_shared<PointCloudStreamBinary<TVertex>>(filenameOut);
	}

	template <typename TVertexSource, typename TVertexTarget>
	std::shared_ptr<PointCloudStreamBinary<TVertexTarget>> convertDownwards(std::shared_ptr<PointCloudStreamBinary<TVertexSource>> cloud, const char* filenameOut)
	{
		const int blockSize = 4096;
		TVertexTarget* targetBlock = new TVertexTarget[blockSize];

		FILE* f = fopen(filenameOut, "wb");
		cloud->Process(blockSize, [&f, &targetBlock](TVertexSource* v, size_t nVertices, int thread, void* param)
		{
			for (int i = 0; i < nVertices; ++i)
			{
				targetBlock[i] = (TVertexTarget)v[i];
			}
			fwrite(targetBlock, sizeof(TVertexTarget), nVertices, f);
		}, nullptr, false);
		fclose(f);

		delete[] targetBlock;

		return std::make_shared<PointCloudStreamBinary<TVertexTarget>>(filenameOut);
	}
};