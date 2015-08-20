#pragma once

#include <vector>
#include <stdio.h>
#include <boost/filesystem.hpp>

// Represents a union-find data structure where each node is augmented with a sign.
// The template parameter specifies if the sign should be used. Otherwise, it is a
// plain union-find data structure.
template <bool UseSign>
class SignedUnionFindBase
{
	// Node is root iff parent index == index
	// Sign is product of path to root

public:
	typedef unsigned int index_t;

	// Saves the entire structure to a file for later usage
	void saveToFile(const char* filename) const
	{
		FILE* file = fopen(filename, "wb");
		fwrite(&parentIndices[0], sizeof(index_t), parentIndices.size(), file);
		fwrite(&ranks[0], sizeof(unsigned int), ranks.size(), file);
		if (UseSign)
		{
			for (int i = 0; i < signs.size(); ++i)
			{
				fwrite(&signs[i], 1, 1, file);
			}
		}
		fclose(file);
	}

	size_t size() const { return parentIndices.size(); }

	// Loads the entire structure from a file. Existing data in the structure is overridden.
	void loadFromFile(const char* filename)
	{
		size_t sizePerEntry = sizeof(index_t) + sizeof(unsigned int) + (UseSign ? 1 : 0);
		auto entries = boost::filesystem::file_size(filename) / sizePerEntry;
		parentIndices.resize(entries);
		ranks.resize(entries);
		if (UseSign)
			signs.resize(entries);

		FILE* file = fopen(filename, "rb");
		fread(&parentIndices[0], sizeof(index_t), entries, file);
		fread(&ranks[0], sizeof(unsigned int), entries, file);
		if (UseSign)
		{
			for (int i = 0; i < signs.size(); ++i)
			{
				fread(&signs[i], 1, 1, file);
			}
		}
		fclose(file);
	}

	// Adds an item to the structure with the given sign
	void addItem(bool negative)
	{
		parentIndices.push_back((index_t)parentIndices.size());
		ranks.push_back(0);
		signs.push_back(negative);
	}

	// Finds the set representative for a given entry. Two entries are in the same set
	// iff they have the same set representative.
	index_t getRepresentative(index_t index)
	{
		bool signCorrection = false; //pre-multiplied signs of intermediate path (skipped by deflating), start with neutral +

		//Find the root
		index_t current = index;
		while (parentIndices[current] != current)
		{
			current = parentIndices[current];
			if (UseSign)
				signCorrection ^= signs[current];
		}

		//remove root sign
		if (UseSign)
			signCorrection ^= signs[current];

		//Path compression
		index_t root = current;
		current = index;
		while (parentIndices[current] != current)
		{
			index_t i = current;
			current = parentIndices[current];
			parentIndices[i] = root;
			
			if (UseSign)
			{
				bool myCorrection = signCorrection;
				signs[i] = signs[i] ^ myCorrection;

				//remove the next leaf from the sign correction term
				signCorrection ^= signs[current];
			}
		}

		return root;
	}

	// Merges the sets of the two specified entries.
	void merge(index_t i1, index_t i2)
	{
		index_t rep1 = getRepresentative(i1);
		index_t rep2 = getRepresentative(i2);
		if (rep1 == rep2)
			return;

		//Union by rank
		unsigned int rank1 = ranks[rep1];
		unsigned int rank2 = ranks[rep2];

		if (rank1 < rank2)
			concreteMerge(rep2, rep1);
		else if (rank2 < rank1)
			concreteMerge(rep1, rep2);
		else
		{
			concreteMerge(rep1, rep2);
			++ranks[rep1];
		}	
	}

protected:
	void concreteMerge(index_t newRoot, index_t child)
	{
		parentIndices[child] = newRoot;

		//Preserve sign
		if (UseSign)
			signs[child] = signs[child] ^ signs[newRoot];
	}

	std::vector<index_t> parentIndices;
	std::vector<unsigned int> ranks;
	std::vector<bool> signs; //false = +   true = -
};

template <bool UseSign>
class SignedUnionFind : public SignedUnionFindBase<UseSign>
{ };

// Add sign-related functions to the class if UseSign=true
template <>
class SignedUnionFind < true > : public SignedUnionFindBase<true>
{
public:
	// Returns the sign of a given entry.
	// false = +   true = -
	bool getSign(index_t index) const
	{
		bool sign = signs[index];
		// Accumulate the sign along the path to its root
		while (parentIndices[index] != index)
		{
			index = parentIndices[index];
			sign ^= signs[index];
		}
		return sign;

	}

	// Flips the sign of the set that contains the specified entry.
	void flipSign(index_t index)
	{
		index_t rep = getRepresentative(index);
		signs[rep] = !signs[rep];
	}
};