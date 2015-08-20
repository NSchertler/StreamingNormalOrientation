#pragma once

#include "exceptions.h"

#include <stdio.h>
#include <vector>
#include <stdexcept>

template <typename TVertex>
class OpenPointCloudStreamBinary;

template <typename TVertex>
class PointCloudBuffer
{
public:
	PointCloudBuffer(size_t points)
		: points(points), selfManaged(true)
	{
		buffer = new TVertex[points];
		bufferEnd = buffer;
	}

	PointCloudBuffer(TVertex* buffer, size_t points)
		: points(points), buffer(buffer), selfManaged(false)
	{
		bufferEnd = buffer;
	}

	~PointCloudBuffer()
	{
		if (buffer && selfManaged)
			delete[] buffer;
	}

	TVertex* begin() const { return buffer; }
	TVertex* end() const { return bufferEnd; }

	TVertex& front() const { return *buffer; }
	TVertex& back() const { return *(bufferEnd - 1); }

	TVertex& operator[](int i)
	{
		return *(buffer + i);
	}

	size_t NumberOfPoints() const { return bufferEnd - buffer; }
	TVertex* Ptr() const { return buffer; }

private:
	bool selfManaged;
	size_t points;
	TVertex* buffer;
	TVertex* bufferEnd;

	friend class OpenPointCloudStreamBinary<TVertex>;
};

template <typename TVertex>
class CompoundPointCloudBuffer
{
public:
	CompoundPointCloudBuffer(size_t nBuffers, size_t pointsPerBuffer)
	{
		buffer = new TVertex[pointsPerBuffer * nBuffers];
		buffers.resize(nBuffers);
		for (auto i = 0; i < nBuffers; ++i)
			buffers[i] = new PointCloudBuffer<TVertex>(buffer + pointsPerBuffer * i, pointsPerBuffer);
	}

	~CompoundPointCloudBuffer()
	{
		for (auto it = buffers.begin(); it != buffers.end(); ++it)
			if (*it)
				delete *it;

		if (buffer)
			delete[] buffer;
	}

	PointCloudBuffer<TVertex>& operator[](int i)
	{
		return *buffers[i];
	}

private:
	TVertex* buffer;
	std::vector<PointCloudBuffer<TVertex>*> buffers;
};

template <typename TVertex>
class ContiguousPointCloudBuffer
{
public:
	ContiguousPointCloudBuffer(size_t points, OpenPointCloudStreamBinary<TVertex>& stream)
		: points(points), stream(stream)
	{
		buffer = new TVertex[points];
		readPointer = buffer;
		bufferEnd = buffer;
	}


	~ContiguousPointCloudBuffer()
	{
		if (buffer)
			delete[] buffer;
	}

	size_t NumberOfUnprocessedPoints() const { return bufferEnd - readPointer; }

	void SetProcessedUpTo(TVertex* ptrExclusive) { readPointer = ptrExclusive; }

	// Moves the unprocessed part of the buffer to the beginning and reads new data
	void Fill()
	{
		size_t unprocessedPoints = NumberOfUnprocessedPoints();
		if (unprocessedPoints == points)
			return;
		//Move existing data to read buffer beginning
		memmove(buffer, readPointer, unprocessedPoints * sizeof(TVertex));
		readPointer = buffer;
		bufferEnd = buffer + unprocessedPoints;
		size_t remainingSpace = points - unprocessedPoints;
		bufferEnd += stream.Read(bufferEnd, remainingSpace);
	}

	TVertex* begin() const { return readPointer; }
	TVertex* end() const { return bufferEnd; }

	TVertex& front() const { return *readPointer; }
	TVertex& back() const { return *(bufferEnd - 1); }
	
private:
	OpenPointCloudStreamBinary<TVertex>& stream;

	size_t points;
	TVertex* buffer;
	TVertex* readPointer; //points to the first unprocessed point
	TVertex* bufferEnd;
};

template <typename TVertex>
class OpenPointCloudStreamBinary
{
public:
	OpenPointCloudStreamBinary(const std::string& filename, const char* fileOpenMode = "rb")
	{
		f = fopen(filename.c_str(), fileOpenMode);
		if (!f)
			throw FileReadException(filename);
	}

	OpenPointCloudStreamBinary(FILE* f)
		: f(f)
	{
		if (!f)
			throw std::runtime_error("Invalid file passed to OpenPointCloudStreamBinary.");
	}

	OpenPointCloudStreamBinary(const OpenPointCloudStreamBinary&) = delete;
	OpenPointCloudStreamBinary& operator=(const OpenPointCloudStreamBinary&) = delete;

	//move ctor
	OpenPointCloudStreamBinary(OpenPointCloudStreamBinary&& other) : f(other.f)
	{
		other.f = nullptr;
	}

	~OpenPointCloudStreamBinary()
	{
		if (f)
			fclose(f);
	}

	void rewind() { ::rewind(f); }

	// Returns the number of read points.
	size_t Read(TVertex* buffer, size_t maxPoints)
	{
		return fread(buffer, sizeof(TVertex), maxPoints, f);
	}

	void Write(const PointCloudBuffer<TVertex>& buffer, size_t points)
	{
		fwrite(buffer.begin(), sizeof(TVertex), points, f);
	}

	// Returns true iff data has been read.
	bool operator>>(PointCloudBuffer<TVertex>& buffer)
	{
		size_t readSize = fread(buffer.buffer, sizeof(TVertex), buffer.points, f);
		buffer.bufferEnd = buffer.buffer + readSize;
		return readSize > 0;
	}

	void operator<<(const PointCloudBuffer<TVertex>& buffer)
	{
		Write(buffer, buffer.points);
	}

	OpenPointCloudStreamBinary<TVertex>& operator<<(const TVertex& vertex)
	{
		fwrite(&vertex, sizeof(TVertex), 1, f);
		return *this;
	}
	

private:

	FILE* f;
};