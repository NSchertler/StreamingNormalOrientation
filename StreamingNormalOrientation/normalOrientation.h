#pragma once

#include "kdTree.h"
#include "normalOrientationHooks.h"
#include "BlockTimer.h"
#include "utils.h"
#include "PointCloudStreamBinary.h"
#include "exceptions.h"

#include <iostream>
#include <unordered_map>
#include <map>
#include <unordered_set>
#include <set>
#include <vector>
#include <deque>
#include <boost/thread.hpp>
#include <boost/filesystem.hpp>

#include <Eigen/Dense>

#define NORMALIZE_INPUT

struct Neighbor
{
	long long vertexIndex;
	float sqrDistance;

	Neighbor(long long vIndex, float sqrDistance)
		: vertexIndex(vIndex), sqrDistance(sqrDistance)
	{}

	bool operator<(const Neighbor& rhs) const
	{
		return sqrDistance < rhs.sqrDistance;
	}
};

typedef std::priority_queue<Neighbor> NeighborQueue;

/**
TVertex must have Eigen::Vector3f position, Eigen::Vector3f normal
*/
template <typename TVertex>
class NormalOrientation
{
private:
	
	typedef int SegmentIndex;

	struct PointWithSegmentNumber : public TVertex
	{
		PointWithSegmentNumber() {}
		PointWithSegmentNumber(const TVertex& copy) : TVertex(copy), segment(-1) {}

		SegmentIndex segment;
	};

	typedef std::deque<PointWithSegmentNumber> TVertexStore;

	struct GridCell
	{
		std::unordered_set<long long> newPoints;
		std::shared_ptr<KdTree<TVertexStore>> kdTree;
		boost::shared_mutex mutex;
	};

	//Key:   grid index
	//Value: global vertex index
	typedef std::unordered_map<int, GridCell> GridType;

	struct ConnectedComponent
	{
		Eigen::Vector3f bbxMin, bbxMax;
		std::string filename;
		FILE* file;
		float lastEntryX;

		void openFile() { file = fopen(filename.c_str(), "wb"); }
	};

	struct SimpleConnectedComponent
	{
		Eigen::Vector3f bbxMin, bbxMax;
		std::vector<int> vertexIndices;
	};

	struct FinalComponent
	{
		std::string filename;
		unsigned long long fileSize;

		FinalComponent(ConnectedComponent& c)
		{
			filename = c.filename;			
			fileSize = boost::filesystem::file_size(filename);
		}

		bool operator<(const FinalComponent& other) const
		{
			return fileSize < other.fileSize;
		}
	};

	template<typename T>
	struct AccumulatedVote
	{
		T sumPositive;
		T sumNegative;

		AccumulatedVote() : sumPositive(0), sumNegative(0) {}
		AccumulatedVote(T initialVote)
			: sumPositive(0), sumNegative(0)
		{
			addVote(initialVote);
		}

		void addVote(T vote)
		{
			if (vote >= 0)
				sumPositive += vote;
			else
				sumNegative += vote;
		}

		T total() const { return sumPositive + sumNegative; }
	};

	struct SegmentDistanceComparer
	{
		std::map<SegmentIndex, float>& distances;

		SegmentDistanceComparer(std::map<SegmentIndex, float>& distances)
			: distances(distances) {}

		bool operator() (const SegmentIndex lhs, const SegmentIndex rhs) const
		{
			if (lhs == rhs)
				return false;
			return distances.at(lhs) < distances.at(rhs);
		}
	};

	inline int calculateGridIndex(int cellX, int cellY, int cellsX)
	{
		return cellX + cellsX * cellY;
	}

	//calculates a linearized grid index (order: y, z)
	int calculateGridIndex(Eigen::Vector3f& position, Eigen::Vector3f& bbxMin, float gridCellSize, int cellsY)
	{
		int cellY = (int)floor((position.y() - bbxMin.y()) / gridCellSize);
		int cellZ = (int)floor((position.z() - bbxMin.z()) / gridCellSize);

		return calculateGridIndex(cellY, cellZ, cellsY);
	}

	//calculates a linearized grid index (order: y, z)
	int calculateGridIndex(Eigen::Vector3f& position, Eigen::Vector3f& bbxMin, float gridCellSize, int cellsY, float& distanceFromFirstGridEdge, float& distanceFromSecondGridEdge)
	{
		int cellY = (int)floor((position.y() - bbxMin.y()) / gridCellSize);
		int cellZ = (int)floor((position.z() - bbxMin.z()) / gridCellSize);

		distanceFromFirstGridEdge  = position.y() - (bbxMin.y() + cellY * gridCellSize);
		distanceFromSecondGridEdge = position.z() - (bbxMin.z() + cellZ * gridCellSize);

		return calculateGridIndex(cellY, cellZ, cellsY);
	}

	/// Returns the number of performed comparisons with other vertices
	template <bool insertPointInKdTree>
	inline int searchCellForNeighbors(int gridIndex, GridType& grid, NeighborQueue& neighborSet, TVertex v, long long vertexGlobalIndex, unsigned int neighborCount, float searchRadiusSqr, TVertexStore& vertexStore, long long vertexStoreStartIndex)
	{
		int comparisons = 0;
		auto cellEntry = grid.find(gridIndex);
		if (cellEntry != grid.end())
		{
			GridCell& cell = cellEntry->second;
			boost::upgrade_lock<boost::shared_mutex> readLock(cell.mutex);
			auto verticesEnd = cell.newPoints.end();
			for (auto it = cell.newPoints.begin(); it != verticesEnd; ++it)
			{
				long long globalVertexIndex = *it;
				//Look only for points that are located before the current vertex
				if (globalVertexIndex >= vertexGlobalIndex)
					continue;

				TVertex other = vertexStore[(unsigned int)(globalVertexIndex - vertexStoreStartIndex)];
				Eigen::Vector3f diff = other.position - v.position;
				float distSqr = diff.squaredNorm();

				if (distSqr <= searchRadiusSqr && (neighborSet.size() < neighborCount || distSqr < neighborSet.top().sqrDistance))
				{					
					if (neighborSet.size() >= neighborCount)
						neighborSet.pop();
					Neighbor n(*it, distSqr);
					neighborSet.push(n);
				}
			} //for each vertex in cell


#ifdef GATHER_STATISTICS
			comparisons = (int)cell.newPoints.size();
#endif
			if (insertPointInKdTree)
			{
				boost::upgrade_to_unique_lock<boost::shared_mutex> writeLock(readLock);
#ifdef GATHER_STATISTICS
				comparisons +=
#endif 
				cell.kdTree->template findKnn<insertPointInKdTree, Neighbor, true>(v.position, searchRadiusSqr, neighborCount, vertexGlobalIndex, neighborSet);
				cell.newPoints.erase(vertexGlobalIndex);
			}
			else
			{
#ifdef GATHER_STATISTICS
				comparisons +=
#endif
				cell.kdTree->template findKnn<insertPointInKdTree, Neighbor, true>(v.position, searchRadiusSqr, neighborCount, vertexGlobalIndex, neighborSet);
			}
		} //if this grid cell exists
		return comparisons;
	}

	inline int FindNeighbors(int gridIndex, GridType& grid, int cellsY, float gridCellSize, NeighborQueue& neighborSet, TVertex v, long long vertexGlobalIndex, unsigned int neighborCount, float searchRadiusSqr, TVertexStore& vertexStore, long long vertexStoreStartIndex, float distance1, float distance2)
	{
		int comparisons = 0;
		float radiusSqr = searchRadiusSqr;
		//start by searching in the point's cell
		comparisons += searchCellForNeighbors<true>(gridIndex, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
		//update the search radius if necessary
		if (neighborSet.size() >= neighborCount)
			radiusSqr = neighborSet.top().sqrDistance;
		//search in other cells if they might contain closer neighbors
		if (distance1 * distance1 < radiusSqr)
		{
			//look left
			comparisons += searchCellForNeighbors<false>(gridIndex - 1, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
			if (neighborSet.size() >= neighborCount)
				radiusSqr = neighborSet.top().sqrDistance;
			if (distance1 * distance1 + distance2 * distance2 < radiusSqr)
			{
				//look bottom left
				comparisons += searchCellForNeighbors<false>(gridIndex - 1 - cellsY, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
				if (neighborSet.size() >= neighborCount)
					radiusSqr = neighborSet.top().sqrDistance;
			}
			if (distance1 * distance1 + (gridCellSize - distance2) * (gridCellSize - distance2) < radiusSqr)
			{
				//look top left
				comparisons += searchCellForNeighbors<false>(gridIndex - 1 + cellsY, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
				if (neighborSet.size() >= neighborCount)
					radiusSqr = neighborSet.top().sqrDistance;
			}
		}
		if (distance2 * distance2 < radiusSqr)
		{
			//look bottom
			comparisons += searchCellForNeighbors<false>(gridIndex - cellsY, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
			if (neighborSet.size() >= neighborCount)
				radiusSqr = neighborSet.top().sqrDistance;
		}
		if ((gridCellSize - distance1) * (gridCellSize - distance1) < radiusSqr)
		{
			//look right
			comparisons += searchCellForNeighbors<false>(gridIndex + 1, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
			if (neighborSet.size() >= neighborCount)
				radiusSqr = neighborSet.top().sqrDistance;
			if ((gridCellSize - distance1) * (gridCellSize - distance1) + distance2 * distance2 < radiusSqr)
			{
				//look bottom right
				comparisons += searchCellForNeighbors<false>(gridIndex + 1 - cellsY, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
				if (neighborSet.size() >= neighborCount)
					radiusSqr = neighborSet.top().sqrDistance;
			}
			if ((gridCellSize - distance1) * (gridCellSize - distance1) + (gridCellSize - distance2) * (gridCellSize - distance2) < radiusSqr)
			{
				//look top right
				comparisons += searchCellForNeighbors<false>(gridIndex + 1 + cellsY, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
				if (neighborSet.size() >= neighborCount)
					radiusSqr = neighborSet.top().sqrDistance;
			}
		}
		if ((gridCellSize - distance2) * (gridCellSize - distance2) < radiusSqr)
		{
			//look top
			comparisons += searchCellForNeighbors<false>(gridIndex + cellsY, grid, neighborSet, v, vertexGlobalIndex, neighborCount, radiusSqr, vertexStore, vertexStoreStartIndex);
			if (neighborSet.size() >= neighborCount)
				radiusSqr = neighborSet.top().sqrDistance;
		}
		return comparisons;
	}

	inline void closeConnectedComponentFiles(std::unordered_map<SegmentIndex, ConnectedComponent>& connectedComponents, std::vector<FinalComponent>& finalComponents, float maxLastEntryX, float connectedComponentsMinimalSize)
	{
		auto it = connectedComponents.begin();
		while (it != connectedComponents.end())
		{
			if (it->second.lastEntryX < maxLastEntryX)
			{		
				fclose(it->second.file);

				//If the component is smaller than the minimal size...
				if (it->second.bbxMax.x() - it->second.bbxMin.x() < connectedComponentsMinimalSize ||
					it->second.bbxMax.y() - it->second.bbxMin.y() < connectedComponentsMinimalSize ||
					it->second.bbxMax.z() - it->second.bbxMin.z() < connectedComponentsMinimalSize)
				{
					//Remove it altogether
					remove(it->second.filename.c_str());
				}
				else
				{
					finalComponents.push_back(it->second);
				}

				it = connectedComponents.erase(it);
			}
			else
				++it;
		}
	}

public:
	
	//Returns number of generated files
	template <typename FlipCriterion = FlipCriterionHoppe, typename OrientationSolver = OrientationSolverMST>
	int OrientNormalsStreaming(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, float searchRadius, int neighborCount, float connectedComponentsMinimalSize, unsigned int segmentMinSize, float minSegmentVoteCertainty, float minAccumulatedVoteCertainty, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{
		auto cloudStream = cloud->OpenStream();

		int nextConnectedComponent = 0;

		const float gridCellSize = 1.0f * searchRadius;

		const int readElementCount = 4096;		
		const float searchRadiusSqr = searchRadius * searchRadius;		

		//size of the vertex store at which to remove old elements
		unsigned int vertexStoreCleanSize = 32768;

		//number of cells in the Y direction
		int cellsY = (int)floor((bbxMax.y() - bbxMin.y()) / gridCellSize) + 1;

		//Points that leave the vertex store will be written to this file.
		//This file contains the same points as the input file (in the same order), but some normals
		//might be flipped due to the first in-segment orientation pass
		OpenPointCloudStreamBinary<PointWithSegmentNumber> tempCloudStream(tmpfile());		

		BlockTimer timer;
		auto totalTime = timer.newBlock("Total");
		auto neighborSearchTime = timer.newBlock("Neighbor Search");
		auto segmentationTime = timer.newBlock("Segmentation");
		auto optimizationTime = timer.newBlock("Optimization");
		auto finalizationTime = timer.newBlock("Finalization");

		timer.startBlock(totalTime);

		PointCloudBuffer<TVertex> readBuffer(readElementCount);

		GridType grid;

		// The set of active points
		TVertexStore vertexStore;

		//the global index of the vertex store's first point
		long long vertexStoreStartIndex = 0;

		//[vertex index, score, distance]
		std::vector < NeighborQueue > nearestNeighbors;
			
		SegmentIndex nextFreeSegment = 0;

		OrientationProblem<OrientationSolver::NeedsUnionFindSign> orientationProblem;

		//the index of the currently read vertex
		long long currentVertexIndex = 0;

		long long globalComparisons = 0;

		//As long as there is something to read...
		while (cloudStream >> readBuffer)
		{
			//index of the vertex that has been read first in the current read iteration
			long long firstVertexIndex = currentVertexIndex;				

			//number of vertices of the current read iteration
			unsigned int vertices = (unsigned int)readBuffer.NumberOfPoints();

			//for each read vertex
			for (auto it = readBuffer.begin(); it != readBuffer.end(); ++it)
			{
				//sort vertex into grid
				PointWithSegmentNumber vs(*it);
				int gridIndex = calculateGridIndex(vs.position, bbxMin, gridCellSize, cellsY);

#ifdef NORMALIZE_INPUT
				if (vs.normal.norm() > 0)
					vs.normal.normalize();
#endif

				vertexStore.push_back(vs);
				if (grid.find(gridIndex) == grid.end())
					grid[gridIndex].kdTree = std::make_shared<KdTree<TVertexStore>>(vertexStore, vertexStoreStartIndex, 1, 2);
				grid.at(gridIndex).newPoints.insert(currentVertexIndex++);
			}
				
			//Make room for neighbors
			if (nearestNeighbors.size() < vertices)
				nearestNeighbors.resize(vertices);
				
			int localComparisons = 0;
			//Look for neighbors in current and surrounding grid cells
			timer.startBlock(neighborSearchTime);
#pragma omp parallel for
			for (int i = 0; i < (int)vertices; ++i)
			{
				TVertex& v = readBuffer[i];

				//distance to the respective grid cells
				float distance1, distance2;
				int gridIndex = calculateGridIndex(v.position, bbxMin, gridCellSize, cellsY, distance1, distance2);

				NeighborQueue& neighborSet = nearestNeighbors[i];
				neighborSet = NeighborQueue(); //reset the queue

				int comparisons = FindNeighbors(gridIndex, grid, cellsY, gridCellSize, neighborSet, v, i + firstVertexIndex, neighborCount, searchRadiusSqr, vertexStore, vertexStoreStartIndex, distance1, distance2);
				
#ifdef GATHER_STATISTICS
#pragma omp atomic
				localComparisons += comparisons;
#endif
			} //parallel for

#ifdef GATHER_STATISTICS
			globalComparisons += localComparisons;
#endif
				
			timer.stopBlock(neighborSearchTime);

			//records the cumulative votes for a specific point for each surrounding segment
			//Key  : segment index
			//Value: accumulated vote
			std::unordered_map<SegmentIndex, AccumulatedVote<float>> votesForSegment;

			//assign a segment to each vertex
			timer.startBlock(segmentationTime);
			for (unsigned int i = 0; i < vertices; ++i)
			{
				PointWithSegmentNumber& vs = vertexStore[(unsigned int)(i + firstVertexIndex - vertexStoreStartIndex)];
				auto& neighborSet = nearestNeighbors[i];

				vs.segment = -1;
				float bestSegmentVote = minAccumulatedVoteCertainty;
				votesForSegment.clear();

				if (neighborSet.size() > 0)
				{					
					std::map<SegmentIndex, float> segmentDistances;
					SegmentDistanceComparer comparer(segmentDistances);
					std::set<SegmentIndex, SegmentDistanceComparer> neighborSegments(comparer);

					//accumulate orientation votes for each adjacent segment
					while (neighborSet.size() > 0)
					{
						const Neighbor* n = &neighborSet.top();
						PointWithSegmentNumber& neighbor = vertexStore[(unsigned int)(n->vertexIndex - vertexStoreStartIndex)];

						//weight the neighbor's vote by its distance
						float neighborVote = FlipCriterion::calculateVote(
							vs.normal, neighbor.normal,
							vs.position, neighbor.position)
							* (1.0f - n->sqrDistance / searchRadiusSqr);						

						//Record the vote
						SegmentIndex neighborSegment = neighbor.segment;
						auto segmentVote = votesForSegment.find(neighborSegment);
						if (segmentVote == votesForSegment.end())
							votesForSegment.emplace(neighborSegment, neighborVote);
						else
							votesForSegment.at(neighborSegment).addVote(neighborVote);

						//Update the neighbor distance set
						auto entry = segmentDistances.find(neighborSegment);
						if (entry == segmentDistances.end())
						{
							//The segment is not recorded yet
							segmentDistances[neighborSegment] = n->sqrDistance;
							neighborSegments.insert(neighborSegment);
						}
						else
						{
							if (n->sqrDistance < entry->second)
							{
								//Update the order
								neighborSegments.erase(neighborSegment);
								entry->second = n->sqrDistance;
								neighborSegments.insert(neighborSegment);
							}
						}

						neighborSet.pop();
					}

					//Find eligible segments for this vertex

					auto segmentIt = neighborSegments.begin();
					float maxSegmentDistance = 1.69f * segmentDistances.at(*segmentIt); //distances are squared => maxDist = 1.3 nearestDist
					for (; segmentIt != neighborSegments.end(); ++segmentIt)
					{
						SegmentIndex segment = *segmentIt;
						if (segmentDistances.at(segment) > maxSegmentDistance)
							break;

						//Check intra-segment criterion
						//All neighbors in segment must vote for the same orientation
						AccumulatedVote<float> intraSegmentVotes = votesForSegment.at(segment);
						if ((intraSegmentVotes.sumNegative < -minSegmentVoteCertainty && intraSegmentVotes.sumPositive > minSegmentVoteCertainty) || abs(intraSegmentVotes.total()) < abs(bestSegmentVote))
							continue;

						bool eligible = true;

						//Check inter-segment criterion
						for (auto votesIt = votesForSegment.begin(); votesIt != votesForSegment.end(); ++votesIt)
						{
							SegmentIndex forSegment = votesIt->first;
							AccumulatedVote<float> vote = votesIt->second;
							if (forSegment != segment) //consider only inter-segment edges
							{
								float currentSegmentVote = orientationProblem.getEdgeWeight(segment, forSegment);
								float addedSegmentVote = vote.total() * (intraSegmentVotes.total() > 0 ? 1.0f : -1.0f);

								//compare the signs of the edge weights
								if ((currentSegmentVote >= 0) != (addedSegmentVote >= 0) && currentSegmentVote != 0 && addedSegmentVote != 0 && abs(addedSegmentVote) > minSegmentVoteCertainty)
								{
									eligible = false;
									break;
								}
							}
						}
						if (eligible)
						{
							vs.segment = segment;
							bestSegmentVote = intraSegmentVotes.total();
						}
					}
				}

				//orient normal according to the segment's orientation

				// +1.0f if vertex normal is unflipped, -1.0f otherwise
				float vertexNormalFlipped = 1.0f;
						
				if (vs.segment < 0)
				{
					// create a new segment for this vertex
					vs.segment = nextFreeSegment++;
					orientationProblem.addNode();
				}
				else
				{
					if (bestSegmentVote < 0)
					{
						//flip the normal if necessary
						vs.normal *= -1.0f;
						vertexNormalFlipped = -1.0f;
					}
				}

				//vote for segment orientation
				for (auto it = votesForSegment.begin(); it != votesForSegment.end(); ++it)
				{
					SegmentIndex neighboringSegment = it->first;
					if (neighboringSegment != vs.segment)
					{
						orientationProblem.template addEdgeWeight<true>(vs.segment, neighboringSegment, it->second.total() * vertexNormalFlipped);
					} //if neighboringSegment != closestSegment
				} //for each segment vote																
			} //for each vertex - assign segment number
			timer.stopBlock(segmentationTime);

#ifndef ACCURATE_TIMING
			std::cout << "\rProcessed " << currentVertexIndex << " points (" << (localComparisons / vertices) << " comp/vert, " << nextFreeSegment << " segments).    ";
#endif
			if (vertexStore.size() >= vertexStoreCleanSize)
			{
				//It's time to get rid of some old vertices

				unsigned int deleted = 0;
				//Vertices up to this x-coordinate won't be needed anymore.
				float deleteUpToX = vertexStore.back().position.x() - searchRadius;
				auto it = vertexStore.begin();
				//vertex store is sorted in x-direction, so we can just examine point by point
				while (it->position.x() < deleteUpToX)
				{
					PointWithSegmentNumber& vs = *it;
					int gridIndex = calculateGridIndex(vs.position, bbxMin, gridCellSize, cellsY);

					tempCloudStream << vs;

					grid.at(gridIndex).kdTree->remove(vertexStoreStartIndex + deleted);

					deleted++;
					it++;
				}
				vertexStore.erase(vertexStore.begin(), it);
					
				if (deleted <= 100)
					//If very few elements have been removed, increase the clean size.
					//This avoids many few-element clean-ups that might become slow in favor of fewer large clean-ups
					vertexStoreCleanSize *= 2;
				if (deleted > vertexStoreCleanSize / 2)
					//But we want to keep the vertex store as small as possible, so adjust the clean size in the other
					//direction, too
					vertexStoreCleanSize /= 2;

				vertexStoreStartIndex += deleted;
			} //if clean-up necessary
		} //while reading

#ifndef ACCURATE_TIMING
		std::cout << std::endl;

		std::cout << "Finished reading with " << (globalComparisons / currentVertexIndex) << " comparisons per vertex in average." << std::endl;
#endif

		//write remaining vertices in read buffer to file
		while (!vertexStore.empty())
		{
			PointWithSegmentNumber& vs = vertexStore.front();

			tempCloudStream << vs;

			vertexStore.pop_front();
		}

#ifndef ACCURATE_TIMING
		std::cout << "Orienting segments..." << std::endl;
#endif
		{
			timer.startBlock(optimizationTime);
			OrientationSolver solver;
			//orient segments consistently
			solver.solve(orientationProblem, &timer);
			timer.stopBlock(optimizationTime);
		}

		//Now finalize the orientation process		

#ifndef ACCURATE_TIMING
		std::cout << "Finalizing orientation..." << std::endl;
#endif
		tempCloudStream.rewind();
		PointCloudBuffer<PointWithSegmentNumber> readTempBuffer(readElementCount);

		//Calculate set of connected components on the fly
		std::unordered_map<SegmentIndex, ConnectedComponent> connectedComponents;
		//set of components that are large enough to count as a real component
		std::vector<FinalComponent> finalComponents;

		nextConnectedComponent = 0;
			
		timer.startBlock(finalizationTime);
		while (tempCloudStream >> readTempBuffer)
		{			
			for (auto& vs : readTempBuffer)
			{
				//apply the segment's sign to the normal
				vs.normal *= orientationProblem.getOrientationFactor(vs.segment);						

				//colorize vertex according to segment number					
				float r, g, b;
				getSegmentColor(vs.segment, r, g, b);
				vs.color << (unsigned char)(r * 255), (unsigned char)(g * 255), (unsigned char)(b * 255);

				int componentId = orientationProblem.getComponentId(vs.segment);
				auto componentIt = connectedComponents.find(componentId);
				if (componentIt == connectedComponents.end())
				{
					//we have found a new connected component

					float currentX = vs.position.x();

					//create new entry in the map
					auto& cmp = connectedComponents[componentId];
					//generate a temporary filename
					cmp.filename = filenameOut + std::to_string(nextConnectedComponent++) + std::string(".temp.bin");
					cmp.openFile();
					cmp.lastEntryX = currentX;
					using namespace Eigen;
					//initialize bounding box
					cmp.bbxMin << std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity();
					cmp.bbxMax << -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity();

					//if fopen was not successful, we just might have too many open files
					if (cmp.file == nullptr)
					{
						//try to close some files with complete connected components
						closeConnectedComponentFiles(connectedComponents, finalComponents, currentX - searchRadius, connectedComponentsMinimalSize);
						cmp.openFile();
					}
				}

				auto& component = connectedComponents.at(componentId);
				component.lastEntryX = vs.position.x();
				//write the vertex to the component's file
				TVertex v = vs;
				fwrite(&v, sizeof(TVertex), 1, component.file);

				//Update the component's bounding box
				for (int i = 0; i < 3; ++i)
				{
					if (vs.position[i] < component.bbxMin[i])
						component.bbxMin[i] = vs.position[i];
					if (vs.position[i] > component.bbxMax[i])
						component.bbxMax[i] = vs.position[i];
				}
					
			}
		} //while reading and processing points

		closeConnectedComponentFiles(connectedComponents, finalComponents, std::numeric_limits<float>::infinity(), connectedComponentsMinimalSize);

		nextConnectedComponent = 0;
		//sort components by file size
		std::sort(finalComponents.begin(), finalComponents.end());
		for (auto it = finalComponents.rbegin(); it != finalComponents.rend(); ++it)
		{
			//Assign a reasonable filename
			std::string newFilename = filenameOut + std::to_string(nextConnectedComponent++) + std::string(".bin");
			remove(newFilename.c_str());
			rename(it->filename.c_str(), newFilename.c_str());
		}

		timer.stopBlock(finalizationTime);
			
		timer.stopBlock(totalTime);
		std::cout << timer;
		

		return nextConnectedComponent;
	}

private:
	struct VertexStoreData
	{
		TVertex* vertexStore;

		~VertexStoreData()
		{
			if (vertexStore)
				delete[] vertexStore;
		}
	};
public:

	//Returns number of generated files
	template <typename FlipCriterion = FlipCriterionHoppe, typename OrientationSolver = OrientationSolverMST>
	int OrientNormals(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, float searchRadius, float connectedComponentsMinimalSize, const int neighborCount, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{
		cloud->AssertValid();

		const float searchRadiusSqr = searchRadius * searchRadius;
		
		BlockTimer timer;
		BlockTimer::BlockId totalTime = timer.newBlock("Total");
		BlockTimer::BlockId neighborGraphTime = timer.newBlock("Build Neighbor Graph");
		BlockTimer::BlockId solveTime = timer.newBlock("Solve MRF");
		BlockTimer::BlockId finalizationTime = timer.newBlock("Finalization");
			
		int nVertices = cloud->NumberOfPoints();

		VertexStoreData localData;
		localData.vertexStore = new TVertex[nVertices];

		cloud->CopyPointsToArray(localData.vertexStore);

		long long vertexStoreStartIndex = 0;

#ifndef ACCURATE_TIMING
		std::cout << "Shuffle" << std::endl;
#endif
		std::random_shuffle(localData.vertexStore, localData.vertexStore + nVertices);

		OrientationProblem<OrientationSolver::NeedsUnionFindSign> orientationProblem;

		timer.startBlock(totalTime);
		timer.startBlock(neighborGraphTime);

		{
			//kd tree scope
			KdTree<TVertex*> kd(localData.vertexStore, vertexStoreStartIndex, 0, 2);
#ifndef ACCURATE_TIMING
			std::cout << "Building kd tree..." << std::endl;
#endif
			for (int i = 0; i < nVertices; ++i)
				kd.insert(i);

#ifndef ACCURATE_TIMING
			std::cout << "Constructing neighbor graph..." << std::endl;
#endif

#ifdef NORMALIZE_INPUT
#pragma omp parallel for
			for (int i = 0; i < nVertices; ++i)
			{
				if (localData.vertexStore[i].normal.norm() > 0)
					localData.vertexStore[i].normal.normalize();
			}
#endif

#pragma omp parallel for
			for (int i = 0; i < nVertices; ++i)
			{
#pragma omp critical
				{
					orientationProblem.addNode();
				}

				//Look for neighbors in current and surrounding grid cells
				std::priority_queue < Neighbor > neighbors;
				kd.template findKnn<false, Neighbor>(localData.vertexStore[i].position, searchRadiusSqr, neighborCount, i, neighbors);

				while (neighbors.size() > 0)
				{
					auto& n = neighbors.top();
					float neighborVote = FlipCriterion::calculateVote(
						localData.vertexStore[i].normal, localData.vertexStore[n.vertexIndex].normal,
						localData.vertexStore[i].position, localData.vertexStore[n.vertexIndex].position)
						* (1.0f - n.sqrDistance / searchRadiusSqr);
#pragma omp critical
					{
						orientationProblem.template addEdgeWeight<false>(i, (int)n.vertexIndex, neighborVote);
					}
					neighbors.pop();
				}
					
			} //while reading
		} //kd tree scope

			
		timer.stopBlock(neighborGraphTime);
		timer.startBlock(solveTime);

		{
			//orient segments consistently
			OrientationSolver solver;
			solver.solve(orientationProblem, &timer);
		}

		timer.stopBlock(solveTime);
		timer.startBlock(finalizationTime);

		//Now finalize the orientation process
#ifndef ACCURATE_TIMING
		std::cout << "Finalizing orientation..." << std::endl;
#endif
			
		//Calculate set of connected components on the fly
		std::unordered_map<SegmentIndex, SimpleConnectedComponent> connectedComponents;					

		for (int i = 0; i < nVertices; ++i)
		{	
			TVertex* v = localData.vertexStore + i;
			//apply the segment's sign to the normal
			v->normal *= orientationProblem.getOrientationFactor(i);											

			v->color << 255, 255, 255;

			int componentId = orientationProblem.getComponentId(i);
			auto componentIt = connectedComponents.find(componentId);
			if (componentIt == connectedComponents.end())
			{
				//we have found a new connected component

				//create new entry in the map
				auto& cmp = connectedComponents[componentId];
				//generate a temporary filename					
				using namespace Eigen;
				//initialize bounding box
				cmp.bbxMin << std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity(), std::numeric_limits<float>::infinity();
				cmp.bbxMax << -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity(), -std::numeric_limits<float>::infinity();
			}

			auto& component = connectedComponents.at(componentId);
			component.vertexIndices.push_back(i);
				
			//Update the component's bounding box
			for (int j = 0; j < 3; ++j)
			{
				if (v->position[j] < component.bbxMin[j])
					component.bbxMin[j] = v->position[j];
				if (v->position[j] > component.bbxMax[j])
					component.bbxMax[j] = v->position[j];
			}
				
		} //while reading and processing points			

		std::vector<SimpleConnectedComponent> finalComponents;
		for (auto it = connectedComponents.begin(); it != connectedComponents.end(); ++it)
		{
			if (it->second.bbxMax.x() - it->second.bbxMin.x() > connectedComponentsMinimalSize &&
				it->second.bbxMax.y() - it->second.bbxMin.y() > connectedComponentsMinimalSize &&
				it->second.bbxMax.z() - it->second.bbxMin.z() > connectedComponentsMinimalSize)
			{
				finalComponents.push_back(it->second);
			}
		}

		//sort components by size
		std::sort(finalComponents.begin(), finalComponents.end(), [](const SimpleConnectedComponent& c1, const SimpleConnectedComponent& c2){ return c1.vertexIndices.size() > c2.vertexIndices.size(); });
		for (int i = 0; i < finalComponents.size(); ++i)
		{
			//Assign a reasonable filename
			std::string newFilename = filenameOut + std::to_string(i) + std::string(".bin");
			FILE* f = fopen(newFilename.c_str(), "wb");
			for (auto it = finalComponents[i].vertexIndices.begin(); it != finalComponents[i].vertexIndices.end(); ++it)
			{
				fwrite(localData.vertexStore + *it, sizeof(TVertex), 1, f);
			}
			fclose(f);
		}

		timer.stopBlock(finalizationTime);
		timer.stopBlock(totalTime);

		std::cout << timer;

		return (int)finalComponents.size();
	}

	//Calculates the orientation energy for an input file with respect to the chosen flip criterion.
	template <typename FlipCriterion = FlipCriterionHoppe>
	double CalculateEnergy(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, float searchRadius, const int neighborCount, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{		
		cloud->AssertValid();

		const float searchRadiusSqr = searchRadius * searchRadius;
		
		int nVertices = (int)cloud->NumberOfPoints();

		VertexStoreData localData;
		localData.vertexStore = new TVertex[nVertices];

		cloud->CopyPointsToArray(localData.vertexStore);

		long long vertexStoreStartIndex = 0;

		// Shuffling the vertices makes the kd-tree perform better
		std::random_shuffle(localData.vertexStore, localData.vertexStore + nVertices);

		double energy = 0.0f;

		{
			//kd tree scope
			KdTree<TVertex*> kd(localData.vertexStore, vertexStoreStartIndex, 0, 2);
			for (int i = 0; i < nVertices; ++i)
				kd.insert(i);

#ifdef NORMALIZE_INPUT
#pragma omp parallel for
			for (int i = 0; i < nVertices; ++i)
			{
				if (localData.vertexStore[i].normal.norm() > 0)
					localData.vertexStore[i].normal.normalize();
			}
#endif

			std::set<std::pair<int, int>> visitedEdges;

#pragma omp parallel for
			for (int i = 0; i < nVertices; ++i)
			{
				std::priority_queue < Neighbor > neighbors;
				kd.template findKnn<false, Neighbor>(localData.vertexStore[i].position, searchRadiusSqr, neighborCount, i, neighbors);

				while (neighbors.size() > 0)
				{
					auto& n = neighbors.top();
					double neighborVote = FlipCriterion::calculateVote(
						localData.vertexStore[i].normal, localData.vertexStore[n.vertexIndex].normal,
						localData.vertexStore[i].position, localData.vertexStore[n.vertexIndex].position)
						* (1.0f - n.sqrDistance / searchRadiusSqr);
#pragma omp critical
					{
						if (neighborVote < 0)
						{
							std::pair<int, int> p(std::min<int>(i, (int)n.vertexIndex), std::max<int>(i, n.vertexIndex));
							if (visitedEdges.find(p) == visitedEdges.end())
							{
								visitedEdges.insert(p);
								energy -= neighborVote;
							}
						}
					}
					neighbors.pop();
				}

			} //while reading
		} //kd tree scope

		return energy;
	}
};