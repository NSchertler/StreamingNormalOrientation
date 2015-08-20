//#define SAVE_MODEL

#if defined ACCURATE_TIMING || !defined WITH_HDF5
#undef SAVE_MODEL
#endif

#include "SignedUnionFind.h"
#include "OrientationProblem.h"
#include "BlockTimer.h"

#include <unordered_set>
#include <memory>

#include <Eigen/Dense>

#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/graphicalmodel/space/discretespace.hxx>
#include <opengm/inference/external/qpbo.hxx>
#ifdef WITH_LSATR
#include <opengm/inference/lsatr.hxx>
#endif
#include <opengm/functions/function_registration.hxx>
#include <opengm/functions/function_properties_base.hxx>
#ifdef WITH_HDF5
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>
#endif

// ####  Flip Criteria  ####
/* Interface:
static float calculateVote(Eigen::Vector3f& n1, Eigen::Vector3f& n2,
Eigen::Vector3f& p1, Eigen::Vector3f& p2)
*/
struct FlipCriterionHoppe
{
	inline static float calculateVote(Eigen::Vector3f& n1, Eigen::Vector3f& n2,
		Eigen::Vector3f& p1, Eigen::Vector3f& p2)
	{
		return n1.dot(n2);
	}
};

struct FlipCriterionXie
{
	inline static float calculateVote(Eigen::Vector3f& n1, Eigen::Vector3f& n2,
		Eigen::Vector3f& p1, Eigen::Vector3f& p2)
	{
		Eigen::Vector3f e = p1 - p2;
		e.normalize();

		Eigen::Vector3f n1Prime = n1 - 2 * e.dot(n1) * e;
		return n1Prime.dot(n2);
	}
};

// ####  Orientation Solvers  ####
/* Interface:
static const bool NeedsUnionFindSign;
void solve(OrientationProblem<NeedsUnionFindSign>& problem, BlockTimer* timer = nullptr);
*/


// Solves the orientation problem with a maximum spanning tree and the Signed Union Find data structure
class OrientationSolverMST
{
public:
	static const bool NeedsUnionFindSign = true;
	typedef OrientationProblem<NeedsUnionFindSign> TOrientationProblem;

	void solve(TOrientationProblem& problem, BlockTimer* timer = nullptr)
	{
		BlockTimer::BlockId blockId;
		if (timer != nullptr)
		{
			blockId = timer->newBlock("MST Solver");
			timer->startBlock(blockId);
		}

#ifndef ACCURATE_TIMING
		std::cout << "Starting MST solver..." << std::endl;
#endif

		std::vector<TOrientationProblem::Edge> edges;
		problem.getAllEdges(edges);
		
		//sort edges in descending order
		std::sort(edges.begin(), edges.end(), [](const TOrientationProblem::Edge& edge1, const TOrientationProblem::Edge& edge2) {return abs(edge1.weight) > abs(edge2.weight); });

		unsigned int unorientableEdges = 0;
		for (auto edge = edges.begin(); edge != edges.end(); ++edge)
		{
			auto repThis = problem._uf.getRepresentative(edge->firstNode);
			auto repNeighbor = problem._uf.getRepresentative(edge->secondNode);

			bool nodesHaveDifferentSigns = problem._uf.getSign(edge->firstNode) != problem._uf.getSign(edge->secondNode);
			bool nodesShouldHaveDifferentSigns = edge->weight < 0;
			bool nodesOrientedDifferently = nodesHaveDifferentSigns != nodesShouldHaveDifferentSigns;

			if (repThis == repNeighbor)
			{
				if (nodesOrientedDifferently)
				{
					++unorientableEdges;
				}
			}
			else
			{
				//both segments are not merged yet

				if (nodesOrientedDifferently)
					problem._uf.flipSign(repThis > repNeighbor ? edge->firstNode : edge->secondNode);
				problem._uf.merge(edge->firstNode, edge->secondNode);
			}
		}

		if (timer != nullptr)
		{
			timer->stopBlock(blockId);
		}

#ifdef SAVE_MODEL
		problem._uf.saveToFile("nodeUnionFind.bin");
#endif
		
		if (problem._solution.size() < problem._uf.size())
			problem._solution.resize(problem._uf.size());
		for (int i = 0; i < problem._uf.size(); ++i)
			problem._solution[i] = problem._uf.getSign(i);

#ifndef ACCURATE_TIMING
		std::cout << "Found " << unorientableEdges << " unorientable edges (" << 100.0f * unorientableEdges / (edges.size()) << " %)" << std::endl;
#endif
	}
};

// Solves the orientation problem with by propagation along the spanning tree
class OrientationSolverTraditionalMST
{
public:
	static const bool NeedsUnionFindSign = false;
	typedef OrientationProblem<NeedsUnionFindSign> TOrientationProblem;

	void solve(TOrientationProblem& problem, BlockTimer* timer = nullptr)
	{
		BlockTimer::BlockId blockId;
		if (timer != nullptr)
		{
			blockId = timer->newBlock("TMST Solver");
			timer->startBlock(blockId);
		}

#ifndef ACCURATE_TIMING
		std::cout << "Starting TMST solver..." << std::endl;
#endif

		std::vector<TOrientationProblem::Edge> edges;
		problem.getAllEdges(edges);

		//sort edges in descending order
		std::sort(edges.begin(), edges.end(), [](const TOrientationProblem::Edge& edge1, const TOrientationProblem::Edge& edge2) {return abs(edge1.weight) > abs(edge2.weight); });

		//structure for MST
		std::unordered_map<int, std::unordered_map<int, float>> mst;

		//construct the MST
		for (auto edge = edges.begin(); edge != edges.end(); ++edge)
		{
			auto repThis = problem._uf.getRepresentative(edge->firstNode);
			auto repNeighbor = problem._uf.getRepresentative(edge->secondNode);

			if (repThis != repNeighbor)
			{
				problem._uf.merge(edge->firstNode, edge->secondNode);
				mst[edge->firstNode][edge->secondNode] = edge->weight;
				mst[edge->secondNode][edge->firstNode] = edge->weight;
			}			
		}

		problem._solution.resize(problem._uf.size(), false);

		//propagate along the MST
		std::queue<int> openPoints;
		openPoints.emplace(0);

		while (!openPoints.empty())
		{
			int currentPoint = openPoints.front();
			openPoints.pop();

			float currentPointFactor = (problem._solution[currentPoint] ? -1.0f : 1.0f);
			for (auto mstNeighbor : mst[currentPoint])
			{
				float currentNeighborFactor = (problem._solution[mstNeighbor.first] ? -1.0f : 1.0f);

				//flip orientation for negative flip criterion
				if (mstNeighbor.second * currentPointFactor * currentNeighborFactor < 0)
					problem._solution[mstNeighbor.first] = !problem._solution[mstNeighbor.first];

				//break the back link
				mst[mstNeighbor.first].erase(currentPoint);

				openPoints.push(mstNeighbor.first);
			}
		}

		if (timer != nullptr)
		{
			timer->stopBlock(blockId);
		}

#ifdef SAVE_MODEL
		problem._uf.saveToFile("nodeUnionFind.bin");
#endif		
	}
};

//Binary factor function for the orientation problem
template <class T, class I = size_t, class L = size_t>
class OrientationFunction
	: public opengm::FunctionBase < OrientationFunction<T, I, L>, T, I, L >
{
public:
	typedef T ValueType;
	typedef I LabelType;
	typedef L IndexType;

	OrientationFunction(T weight) : weight(weight) {}

	LabelType shape(const size_t j) const { return 2; }; //all nodes have two possible orientations
	template <class ITERATOR>
	ValueType operator()(ITERATOR it) const
	{
		bool isSame = *it == *(it + 1);
		if (isSame)
		{
			if (weight > 0)
				return 0;
			else
				return -weight;
		}
		else
		{
			if (weight > 0)
				return weight;
			else
				return 0;
		}
	}
	
	bool operator==(const OrientationFunction& other) const { return weight == other.weight; }

private:
	T weight;
};

template <class T, class I, class L>
struct opengm::FunctionRegistration < OrientationFunction<T, I, L> > {
	enum ID { Id = 0 };
};

//Base solver for the QPBO solver that does nothing except calculating connected components
class VoidBaseSolver
{
public:
	static const bool NeedsUnionFindSign = false;
	typedef OrientationProblem<NeedsUnionFindSign> TOrientationProblem;

	void solve(TOrientationProblem& problem, BlockTimer* timer = nullptr)
	{
		BlockTimer::BlockId blockId;
		if (timer != nullptr)
		{
			blockId = timer->newBlock("Connected Components");
			timer->startBlock(blockId);
		}

		if (problem._solution.size() < problem._uf.size())
			problem._solution.resize(problem._uf.size());
		for (int i = 0; i < problem._uf.size(); ++i)
			problem._solution[i] = false;

		for (auto node = problem._orientationEdges.begin(); node != problem._orientationEdges.end(); ++node)
		{
			auto& outgoingEdges = node->second;
			for (auto connectedNode = outgoingEdges.begin(); connectedNode != outgoingEdges.end(); ++connectedNode)
			{
				problem._uf.merge(node->first, connectedNode->first);				
			}
		}

		if (timer != nullptr)
		{
			timer->stopBlock(blockId);
		}
	}
};

typedef opengm::GraphicalModel<double, opengm::Adder,
	opengm::meta::TypeList<
#ifndef SAVE_MODEL
	OrientationFunction<double>, opengm::meta::TypeList<
#endif
	opengm::ExplicitFunction<double>, opengm::meta::ListEnd >
#ifndef SAVE_MODEL
	>
#endif
> Model;

struct OpenGMAlgorithmQPBO
{
	static const std::string NAME;

	typedef opengm::external::QPBO<Model> Optimizer;

	static std::unique_ptr<Optimizer> getOptimizer(const Model& gm)
	{	
		Optimizer::Parameter param;
		param.strongPersistency_ = false;
		//param.useProbeing_ = true;
		param.useImproveing_ = true;

		return std::make_unique<Optimizer>(gm, param);		
	}
};

#ifdef WITH_LSATR
struct OpenGMAlgorithmLSATR
{
	static const std::string NAME;

	typedef opengm::LSA_TR<Model, opengm::Minimizer> Optimizer;

	static std::unique_ptr<Optimizer> getOptimizer(const Model& gm)
	{
		return std::make_unique<Optimizer>(gm);
	}
};
#endif

// This class is a wrapper for all OpenGM optimizers. It prepares the model and interprets the optimizer's results.
// TAlgorithm  - Implementation of the actual algorithm that solves the problem
// TBaseSolver - The IOrientationSolver that is executed before the QPBO-I passes
// CleanEdges  - Specifies whether the solver can delete edges in order to free memory
template <class TAlgorithm, class TBaseSolver = VoidBaseSolver, bool CleanEdges = true>
class OrientationSolverOpenGM
{
public:	
	static const bool NeedsUnionFindSign = TBaseSolver::NeedsUnionFindSign;
	typedef OrientationProblem<NeedsUnionFindSign> TOrientationProblem;

	void solve(TOrientationProblem& problem, BlockTimer* timer = nullptr)
	{
		//Run the base solver
		TBaseSolver baseSolver;
		baseSolver.solve(problem, timer);

		BlockTimer::BlockId blockIdBuildModel, blockIdSolveModel;
		if (timer != nullptr)
		{
			blockIdBuildModel = timer->newBlock("OpenGM Build Model");
			blockIdSolveModel = timer->newBlock(TAlgorithm::NAME + " Solve Model");
			timer->startBlock(blockIdBuildModel);
		}

#ifndef ACCURATE_TIMING
		std::cout << "Running OpenGM solver..." << std::endl;
#endif

		opengm::DiscreteSpace<> labelSpace;

		for (unsigned int i = 0; i < problem._uf.size(); ++i)
			labelSpace.addVariable(2);

		Model gm(labelSpace);

		size_t edgeVariables[2];

#ifndef ACCURATE_TIMING
		std::cout << "Building graphical model..." << std::endl;
#endif

#ifdef SAVE_MODEL
		const size_t functionShape[] = { 2, 2 };
#endif

		auto node = problem._orientationEdges.begin();
		int nonSubmodular = 0;

		while (node != problem._orientationEdges.end())
		{
			auto& outgoingEdges = node->second;
			float factor1 = problem.getOrientationFactor(node->first);
			for (auto connectedNode = outgoingEdges.begin(); connectedNode != outgoingEdges.end(); ++connectedNode)
			{
				if (connectedNode->second != 0)
				{
					float adaptedWeight = connectedNode->second * factor1 * problem.getOrientationFactor(connectedNode->first);
#ifdef SAVE_MODEL
					auto edgeFunction = opengm::ExplicitFunction<double>(functionShape, functionShape + 2);
					edgeFunction(0, 0) = adaptedWeight > 0 ? 0 : -adaptedWeight;
					edgeFunction(1, 1) = adaptedWeight > 0 ? 0 : -adaptedWeight;
					edgeFunction(0, 1) = adaptedWeight > 0 ? adaptedWeight : 0;
					edgeFunction(1, 0) = adaptedWeight > 0 ? adaptedWeight : 0;
#else
					auto edgeFunction = OrientationFunction<double>(adaptedWeight);
#endif
					if (adaptedWeight < 0)
						++nonSubmodular;					

					edgeVariables[0] = std::min(node->first, connectedNode->first);
					edgeVariables[1] = std::max(node->first, connectedNode->first);
					gm.addFactor(gm.addFunction(edgeFunction), edgeVariables, edgeVariables + 2);					
				}
			}
			//free memory
			if (CleanEdges)
				node = problem._orientationEdges.erase(node);
			else
				node++;
		}

		const size_t functionShapeUnary[] = { 2 };
		auto unaryFunction = opengm::ExplicitFunction<double>(functionShapeUnary, functionShapeUnary + 1);
		unaryFunction(0) = 0;
		unaryFunction(1) = 1000; //zeta
		int label = 0;
		gm.addFactor(gm.addFunction(unaryFunction), &label, &label + 1);

		if (timer != nullptr)
			timer->stopBlock(blockIdBuildModel);

#ifdef SAVE_MODEL
		opengm::hdf5::save(gm, "model.h5", "gm");
#endif

		std::vector<size_t> argmin;
#ifndef ACCURATE_TIMING
		argmin.resize(gm.numberOfVariables(), 0);
		std::cout << "Energy before optimization: " << gm.evaluate(argmin) << std::endl;

		std::cout << "Solving " << TAlgorithm::NAME << " with " << gm.numberOfVariables() << " nodes and " << gm.numberOfFactors() << " factors (" << (nonSubmodular * 100.0 / gm.numberOfFactors()) << " % non-submodular)." << std::endl;
#endif

		std::unique_ptr<TAlgorithm::Optimizer> optimizer = TAlgorithm::getOptimizer(gm);

		if (timer)
			timer->startBlock(blockIdSolveModel);

		optimizer->infer();

		if (timer)
			timer->stopBlock(blockIdSolveModel);

		optimizer->arg(argmin);

#ifndef ACCURATE_TIMING
		std::cout << "Energy after optimization: " << optimizer->value() << std::endl;
#endif

		if (problem._solution.size() < problem._uf.size())
			problem._solution.resize(problem._uf.size());
		for (int i = 0; i < argmin.size(); ++i)
		{
			problem._solution[i] = problem._solution[i] ^ (argmin[i] == 1);
		}
	}	
};

typedef OrientationSolverOpenGM<OpenGMAlgorithmQPBO, OrientationSolverMST> OrientationSolverMSTQPBO;
#ifdef WITH_LSATR
typedef OrientationSolverOpenGM<OpenGMAlgorithmLSATR, OrientationSolverMST> OrientationSolverMSTLSATR;
#endif