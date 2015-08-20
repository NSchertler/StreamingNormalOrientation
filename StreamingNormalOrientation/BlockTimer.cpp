#include "BlockTimer.h"

#include <iostream>

std::ostream& operator<<(std::ostream& o, const BlockTimer& timer)
{
	std::ios init(NULL);
	init.copyfmt(o);

	size_t maxNameSize = 5;
	for (auto it = timer.blocks.begin(); it != timer.blocks.end(); ++it)
	{
		if (it->name.length() > maxNameSize)
			maxNameSize = it->name.length();
	}
	o << "block";
	for (int i = 0; i < maxNameSize - 5; ++i)
		o << " ";
	o << "   time (ms)" << std::endl;
	for (int i = 0; i < maxNameSize + 3 + 12; ++i)
		o << "-";	
	o << std::endl;
	o.setf(std::ios::fixed);
	o.precision(1);
	for (auto it = timer.blocks.begin(); it != timer.blocks.end(); ++it)
	{
		o << std::right << std::setw(maxNameSize) << std::setfill(' ') << it->name;
		o << " " << std::right << std::setw(12) << std::setfill('.');
		o << it->microseconds() / 1000.0f << std::endl;
	}

	o.copyfmt(init);

	return o;
}