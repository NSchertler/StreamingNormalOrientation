#include "MortonUtils.cuh"

__host__ __device__ unsigned int spreadBits(unsigned int x)
{
	unsigned int r = x;
	r = (r | (r << 8)) & 0x00FF00FF;
	r = (r | (r << 4)) & 0x0F0F0F0F;
	r = (r | (r << 2)) & 0x33333333;
	r = (r | (r << 1)) & 0x55555555;
	return r;
}

__host__ __device__ unsigned int getMortonIndex(unsigned int u, unsigned int v)
{
	return spreadBits(u) | (spreadBits(v) << 1);
}

__host__ __device__ unsigned int getBinnedMortonIndex(float x, float y, float minx, float miny, float cellSize)
{
	unsigned int u = (unsigned int)((x - minx) / cellSize);
	unsigned int v = (unsigned int)((y - miny) / cellSize);

	return getMortonIndex(u, v);
}