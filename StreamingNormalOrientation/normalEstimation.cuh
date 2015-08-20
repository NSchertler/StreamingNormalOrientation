#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

struct VertexPositionMortonIndex
{
	float3 position;
	unsigned int index;
};

extern void calculateNormals(VertexPositionMortonIndex* vertices, float3* normals, float searchRadiusSqr, float xmin, float xmax, float ymin, float zmin, float quadTreeCellSize, unsigned int numVertices, cudaStream_t stream);