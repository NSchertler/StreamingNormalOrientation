#include "SVD.cuh"

#include <stdio.h>

__device__ float3 NormalFromCovariance(float3& c1, float3& c2, float3& c3)
{
#define COMPUTE_U_AS_MATRIX

#include "Kernel_Declarations.cuh"

	Sa11.f = c1.x;
	Sa21.f = c1.y;
	Sa31.f = c1.z;
	Sa12.f = c2.x;
	Sa22.f = c2.y;
	Sa32.f = c2.z;
	Sa13.f = c3.x;
	Sa23.f = c3.y;
	Sa33.f = c3.z;

#include "Main_Kernel_Body.cuh"

	return make_float3(Su13.f, Su23.f, Su33.f);

#undef COMPUTE_U_AS_MATRIX
}