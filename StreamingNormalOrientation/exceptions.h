#pragma once

#include <exception>

class FileReadException : public std::runtime_error
{
public:
	FileReadException(const std::string& filename) : std::runtime_error(std::string("Cannot read file \"") + filename + "\".")
	{}
};