#pragma once

#include "SignedUnionFind.h"
#include <unordered_map>

template <bool UseUnionFindSigns>
struct OrientationProblem
{
	SignedUnionFind<UseUnionFindSigns> _uf;

	//Holds edge information of the neighbor graph	
	//First index:  Larger segment index
	//Second index: Smaller segment index
	//Value:        Result flip criterion. If this is positive, segments are oriented in the same direction.
	std::unordered_map<int, std::unordered_map<int, float>> _orientationEdges;

	std::vector<bool> _solution;


	virtual void addNode();

	template <bool accumulate>
	void addEdgeWeight(int node1, int node2, float weight);

	float getEdgeWeight(int node1, int node2) const;

	float getOrientationFactor(int node) const;

	int   getComponentId(int node);

	struct Edge
	{
		int firstNode;
		int secondNode;
		float weight;
	};

	void getAllEdges(std::vector<Edge>& edges) const;
};