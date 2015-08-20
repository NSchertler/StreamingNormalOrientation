#pragma once

#include <vector>
#include <string>
#include <iomanip>
#include <chrono>

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif


class BlockTimer
{
	typedef std::chrono::high_resolution_clock Clock;

private:

	class BlockInfo
	{
	private:
#ifdef _WIN32
		LARGE_INTEGER start;
		LARGE_INTEGER duration;
#else
		std::chrono::time_point<Clock> start;
		std::chrono::duration<Clock::rep, Clock::period > duration;
#endif

	public:
		std::string name;

		BlockInfo(std::string name)
			: name(name)
		{
#ifdef _WIN32
			duration.QuadPart = 0;
#endif
		}

		void startNow() 
		{
#ifdef _WIN32
			QueryPerformanceCounter(&start);
#else
			start = Clock::now(); 
#endif
		}

		void stopNow() 
		{ 
#ifdef _WIN32
			LARGE_INTEGER t;
			QueryPerformanceCounter(&t);
			duration.QuadPart += t.QuadPart - start.QuadPart;
#else
			duration += Clock::now() - start;
#endif
		}
		long long microseconds() const 
		{
#ifdef _WIN32
			LARGE_INTEGER freq;
			QueryPerformanceFrequency(&freq);
			return (duration.QuadPart * 1000000 / freq.QuadPart);
#else
			return std::chrono::duration_cast<std::chrono::microseconds>(duration).count(); 
#endif
		}
	};

public:
	typedef size_t BlockId;

	BlockId newBlock(std::string name)
	{
		blocks.emplace(blocks.end(), name);
		return blocks.size() - 1;
	}

	void startBlock(BlockId block)
	{
		blocks.at(block).startNow();
	}

	void stopBlock(BlockId block)
	{
		blocks.at(block).stopNow();		
	}	

	friend std::ostream& operator<<(std::ostream& o, const BlockTimer& timer);

private:
	std::vector<BlockInfo> blocks;

};

extern std::ostream& operator<<(std::ostream& o, const BlockTimer& timer);