#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

extern __host__ __device__ unsigned int getMortonIndex(unsigned int u, unsigned int v);

extern __host__ __device__ unsigned int getBinnedMortonIndex(float x, float y, float minx, float miny, float cellSize);