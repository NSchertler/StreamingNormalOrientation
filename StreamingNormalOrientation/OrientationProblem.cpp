#include "OrientationProblem.h"

template <bool UseUnionFindSign>
void  OrientationProblem<UseUnionFindSign>::addNode() { _uf.addItem(false); }

template <bool UseUnionFindSign>
float OrientationProblem<UseUnionFindSign>::getOrientationFactor(int node) const { return _solution[node] ? -1.0f : 1.0f; }

template <bool UseUnionFindSign>
int   OrientationProblem<UseUnionFindSign>::getComponentId(int node) { return _uf.getRepresentative(node); }

template <bool UseUnionFindSign>
template <bool accumulate>
void OrientationProblem<UseUnionFindSign>::addEdgeWeight(int node1, int node2, float weight)
{
	int maxNode = std::max(node1, node2);
	int minNode = std::min(node1, node2);

	auto outgoingEdges = _orientationEdges.find(maxNode);
	if (outgoingEdges == _orientationEdges.end())
	{
		//MaxSegment has no outgoing edges
		auto insertedMap = _orientationEdges.emplace(maxNode, std::unordered_map<int, float>());
		insertedMap.first->second.emplace(minNode, weight);
	}
	else
	{
		//MaxSegment already has outgoing edges
		auto edge = outgoingEdges->second.find(minNode);
		if (edge == outgoingEdges->second.end())
		{
			//there is no outgoing edge to MinSegment
			outgoingEdges->second.emplace(minNode, weight);
		}
		else if (accumulate)
		{
			//there is already an outgoing edge to MinSegment
			edge->second += weight;
		}
	}
}

template void OrientationProblem<false>::addEdgeWeight<false>(int, int, float);
template void OrientationProblem<false>::addEdgeWeight<true>(int, int, float);
template void OrientationProblem<true>::addEdgeWeight<false>(int, int, float);
template void OrientationProblem<true>::addEdgeWeight<true>(int, int, float);

template <bool UseUnionFindSign>
float OrientationProblem<UseUnionFindSign>::getEdgeWeight(int node1, int node2) const
{
	int maxNode = std::max(node1, node2);
	int minNode = std::min(node1, node2);

	auto outgoingEdges = _orientationEdges.find(maxNode);
	if (outgoingEdges == _orientationEdges.end())
	{
		//MaxSegment has no outgoing edges
		return 0;
	}
	else
	{
		//MaxSegment already has outgoing edges
		auto edge = outgoingEdges->second.find(minNode);
		if (edge == outgoingEdges->second.end())
		{
			//there is no outgoing edge to MinSegment
			return 0;
		}
		else
		{
			//there is already an outgoing edge to MinSegment
			return edge->second;
		}
	}
}

template <bool UseUnionFindSign>
void OrientationProblem<UseUnionFindSign>::getAllEdges(std::vector<Edge>& edges) const
{
#ifdef SAVE_MODEL
	FILE* edgeFile = fopen("edges.bin", "wb");
#endif

	//gather all edges of the global connectivity graph
	for (auto node = _orientationEdges.begin(); node != _orientationEdges.end(); ++node)
	{
		auto& outgoingEdges = node->second;
		for (auto connectedNode = outgoingEdges.begin(); connectedNode != outgoingEdges.end(); ++connectedNode)
		{
			edges.push_back({ node->first, connectedNode->first, connectedNode->second });

#ifdef SAVE_MODEL
			if (connectedNode->second != 0)
			{
				Edge e = { node->first, connectedNode->first, connectedNode->second };
				fwrite(&e, sizeof(Edge), 1, edgeFile);
			}
#endif
		}
	}
#ifdef SAVE_MODEL
	fclose(edgeFile);
#endif
}

template struct OrientationProblem < false > ;
template struct OrientationProblem < true >;