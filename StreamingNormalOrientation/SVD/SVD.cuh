#pragma once

#include <cuda_runtime.h>

/**
 Computes the normal from a given covariance matrix using SVD.
 c1, c2, c3 ... column vectors of the covariance matrix
 */
__device__ float3 NormalFromCovariance(float3& c1, float3& c2, float3& c3);