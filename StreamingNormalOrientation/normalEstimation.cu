#include "normalEstimation.cuh"
#include "MortonUtils.cuh"
#include <stdio.h>
#include "SVD/SVD.cuh"
#include "nih/priority_queue.h"


__device__ inline float3 operator*(const float a, const float3 &b) {
	return make_float3(a * b.x, a * b.y, a * b.z);
}
__device__ inline float3 operator-(const float3& a, const float3& b) {
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

struct Neighbor
{
	unsigned int index;
	float distanceSqr;

	__host__ __device__ inline bool operator<(const Neighbor& other) const
	{
		return distanceSqr > other.distanceSqr;
	}
};

#define ThreadsPerBlock 256
const unsigned int maxNeighbors = 48 * 1024 / (sizeof(Neighbor) * ThreadsPerBlock); //limited by shared memory
#define NEIGHBOR_COUNT 20

static_assert(NEIGHBOR_COUNT <= maxNeighbors, "Number of neighbors is too large.");

__global__ void kCalculateNormals(VertexPositionMortonIndex* vertices, float3* normals, float searchRadiusSqr, float xmin, float xmax, float ymin, float zmin, float quadTreeCellSize, unsigned int numVertices)
{
	__shared__ Neighbor sharedMemory[NEIGHBOR_COUNT * ThreadsPerBlock];	

	using namespace nih;

	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	Neighbor* neighborStore = sharedMemory + threadIdx.x * NEIGHBOR_COUNT;

	priority_queue<Neighbor, NEIGHBOR_COUNT> neighborQueue(neighborStore);

	if (tid < numVertices)
	{
		VertexPositionMortonIndex p = vertices[tid];

		float3 cov1, cov2, cov3;
		//float sumWeights = 0.0f;

		int i = tid - 1;
		int direction = -1;

		float3 centroid = make_float3(0.0f, 0.0f, 0.0f);

		//Search neighbors
		//TODO: search in neighboring cells
		while (true)
		{
			if (i < 0 || i >= numVertices || vertices[i].index != vertices[tid].index)
			{
				if (direction == 1)
					break;

				i = tid + 1;
				direction = 1;
				
				if (i >= numVertices || vertices[i].index != vertices[tid].index)
					break;
			}
			
			float3 diff = vertices[i].position - p.position;
			float lengthSqr = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
			//TODO: find a way to reduce warp divergence
			if (lengthSqr < searchRadiusSqr && (neighborQueue.empty() || lengthSqr < neighborQueue.top().distanceSqr))
			{				
				if (neighborQueue.size() == NEIGHBOR_COUNT)
					neighborQueue.pop();
				neighborQueue.push({ (unsigned int)i, lengthSqr });
			}

			i += direction;
		}

		for (int i = 0; i < neighborQueue.size(); ++i)
		{
			auto idx = neighborQueue[i].index;
			centroid.x += vertices[idx].position.x;
			centroid.y += vertices[idx].position.y;
			centroid.z += vertices[idx].position.z;
		}
		centroid = (1.0f / (float)neighborQueue.size()) * centroid;

		float3 diffToCentroid = p.position - centroid;
		cov1.x = diffToCentroid.x * diffToCentroid.x;
		cov1.y = diffToCentroid.x * diffToCentroid.y;
		cov1.z = diffToCentroid.x * diffToCentroid.z;

		cov2.y = diffToCentroid.y * diffToCentroid.y;
		cov2.z = diffToCentroid.y * diffToCentroid.z;

		cov3.z = diffToCentroid.z * diffToCentroid.z;

		//Calculate covariance from neighbors
		for (unsigned int i = 0; i < neighborQueue.size(); ++i)
		{
			diffToCentroid = vertices[neighborQueue[i].index].position - centroid;
			float lengthSqrInv = 1.0f / sqrtf(diffToCentroid.x * diffToCentroid.x + diffToCentroid.y * diffToCentroid.y + diffToCentroid.z * diffToCentroid.z);
			diffToCentroid = lengthSqrInv * diffToCentroid;
			//float3 diffToV = vertices[neighbors[i]].position - p.position;
			//float lengthSqr = diffToV.x * diffToV.x + diffToV.y * diffToV.y + diffToV.z * diffToV.z;

			//float weight = fmaxf(0.0f, 1 - lengthSqr / (NeighborSearchRadius * NeighborSearchRadius));
			//sumWeights += weight;

			cov1.x += diffToCentroid.x * diffToCentroid.x;
			cov1.y += diffToCentroid.x * diffToCentroid.y;
			cov1.z += diffToCentroid.x * diffToCentroid.z;

			cov2.y += diffToCentroid.y * diffToCentroid.y;
			cov2.z += diffToCentroid.y * diffToCentroid.z;

			cov3.z += diffToCentroid.z * diffToCentroid.z;
		}

		cov2.x = cov1.y;
		cov3.x = cov1.z;
		cov3.y = cov2.z;
		
		if (neighborQueue.size() < 3)
			normals[tid] = make_float3(0.0, 0.0, 0.0);		
		else
			normals[tid] = NormalFromCovariance(cov1, cov2, cov3);
	}
}

void calculateNormals(VertexPositionMortonIndex* vertices, float3* normals, float searchRadiusSqr, float xmin, float xmax, float ymin, float zmin, float quadTreeCellSize, unsigned int numVertices, cudaStream_t stream)
{	
	unsigned int blocks = (unsigned int)ceil((float)numVertices / ThreadsPerBlock);

	kCalculateNormals <<< blocks, ThreadsPerBlock, 0, stream >>>(vertices, normals, searchRadiusSqr, xmin, xmax, ymin, zmin, quadTreeCellSize, numVertices);
}