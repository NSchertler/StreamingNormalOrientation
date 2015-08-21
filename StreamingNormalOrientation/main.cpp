#include <iostream>

#include "normalOrientation.h"
#ifdef WITH_CUDA
#include "normalEstimation.h"
#endif
#include "common.h"
#include "basicProcessing.h"
#include "pointCloudWriter.h"
#include "BlockTimer.h"
#include "PointCloudStreamAscii.h"

using namespace std;

template <typename TVertex>
struct OrientationDispatcher
{
	OrientationDispatcher(NormalOrientation<TVertex>& no) : no(no) {}

	virtual int OrientStreaming(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, float searchRadius, int neighborCount, float connectedComponentsMinimalSize, unsigned int segmentMinSize, float minSegmentVoteCertainty, float minAccumulatedVoteCertainty, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax) = 0;
	virtual int Orient(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, float searchRadius, float connectedComponentsMinimalSize, const int neighborCount, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax) = 0;
	virtual double CalculateEnergy(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, float searchRadius, const int neighborCount, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax) = 0;

	std::string name;
protected:
	NormalOrientation<TVertex>& no;
};

template <typename TVertex, typename TFlipCriterion, typename TSolver>
struct OrientationDispatcherTemplate : public OrientationDispatcher<TVertex>
{
	OrientationDispatcherTemplate(NormalOrientation<TVertex>& no) : OrientationDispatcher<TVertex>(no)
	{
		if (std::is_same<TFlipCriterion, FlipCriterionHoppe>::value)
			this->name = "Hoppe";
		else if (std::is_same<TFlipCriterion, FlipCriterionXie>::value)
			this->name = "Xie";
		else
			this->name = "Unknown";
	}

	virtual int OrientStreaming(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, float searchRadius, int neighborCount, float connectedComponentsMinimalSize, unsigned int segmentMinSize, float minSegmentVoteCertainty, float minAccumulatedVoteCertainty, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{
		return this->no.template OrientNormalsStreaming<TFlipCriterion, TSolver>(cloud, filenameOut, searchRadius, neighborCount, connectedComponentsMinimalSize, segmentMinSize, minSegmentVoteCertainty, minAccumulatedVoteCertainty, bbxMin, bbxMax);
	}

	virtual int Orient(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, const char* filenameOut, float searchRadius, float connectedComponentsMinimalSize, const int neighborCount, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{
		return this->no.template OrientNormals<TFlipCriterion, TSolver>(cloud, filenameOut, searchRadius, connectedComponentsMinimalSize, neighborCount, bbxMin, bbxMax);
	}

	virtual double CalculateEnergy(std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud, float searchRadius, const int neighborCount, Eigen::Vector3f& bbxMin, Eigen::Vector3f& bbxMax)
	{
		return this->no.template CalculateEnergy<TFlipCriterion>(cloud, searchRadius, neighborCount, bbxMin, bbxMax);
	}
};


int main(int argc, const char* args[])
{	
	try
	{
#ifdef WIN32
		cout.imbue(locale("en-US"));
#else
		cout.imbue(locale("en_US.UTF-8"));
#endif
	}
	catch (...)
	{
		cout << "Cannot set english local for cout." << endl;
	}

	if (argc < 2)
	{
		cout << "Parameters: filename [options]" << endl;
		cout << "\tfilename        first argument must be the input file name (PLY OR BIN file without extension)!" << endl;
		cout << "options:" << endl;
		cout << "\t-stream         streaming orientation." << endl;
		cout << "\t-prepare        convert the input PLY file to the streaming BIN format." << endl;
		cout << "\t-randomize      randomize normal orientations before the optimization." << endl;
#ifdef WITH_CUDA
		cout << "\t-estimateNormals    run a simple PCA-based normal estimation step." << endl;
#endif
		cout << "\t-calculateEnergies  calculate initial and final energies if file size is less than 50 MB." << endl;
		cout << "\t-k %i           number of neighbors." << endl;
		cout << "\t-r %f           neighborhood search radius." << endl;
		cout << "\t-thetaS %f      theta_S from paper." << endl;
		cout << "\t-thetaAcc %f    theta_acc from paper." << endl;
		cout << "\t-flipcrit [hoppe|xie]  flip criterion." << endl;
		return -1;
	}

	std::string model(args[1]);
	float searchRadius = 1.0f;
	int neighborCount = 6;
	typedef OrientationSolverMSTQPBO TSolver;
	typedef VertexPositionNormalColor TVertex;
	bool stream = false;
	bool randomize = false;
	bool estimateNormals = false;
	bool writeOptInputPly = false;
	bool writeFinalPly = true;
	bool prepare = false;
	bool calculateEnergies = false;
	bool evaluateAccuracy = false;

	float minSegmentVoteCertainty = 0.0; //theta_S from paper
	float minAccumulatedVoteCertainty = 0.5; //theta_acc from paper

	BasicProcessing bp;
	PointCloudWriter pw;
#ifdef WITH_CUDA
	NormalEstimation ne;
#endif
	NormalOrientation<TVertex> no;

	std::shared_ptr<OrientationDispatcher<TVertex>> orientationDispatcher = std::make_shared<OrientationDispatcherTemplate<TVertex, FlipCriterionHoppe, TSolver>>(no);

	//Parse options
	for (int i = 2; i < argc; ++i)
	{
		if (strcmp(args[i], "-stream") == 0)
			stream = true;
		else if (strcmp(args[i], "-prepare") == 0)
			prepare = true;
		else if (strcmp(args[i], "-randomize") == 0)
			randomize = true;
#ifdef WITH_CUDA
		else if (strcmp(args[i], "-estimateNormals") == 0)
			estimateNormals = true;
#endif
		else if (strcmp(args[i], "-calculateEnergies") == 0)
			calculateEnergies = true;
		else if (strcmp(args[i], "-calculateEnergies") == 0)
			calculateEnergies = true;
		else if (strcmp(args[i], "-k") == 0 && i + 1 < argc)
		{
			neighborCount = atoi(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-r") == 0 && i + 1 < argc)
		{
			searchRadius = (float)atof(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-thetaS") == 0 && i + 1 < argc)
		{
			minSegmentVoteCertainty = (float)atof(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-thetaAcc") == 0 && i + 1 < argc)
		{
			minAccumulatedVoteCertainty = (float)atof(args[i + 1]);
			++i;
		}
		else if (strcmp(args[i], "-flipcrit") == 0 && i + 1 < argc)
		{
			if (strcmp(args[i + 1], "hoppe") == 0)
				; // default value is Hoppe
			else if (strcmp(args[i + 1], "xie") == 0)
			{
				orientationDispatcher = std::make_shared<OrientationDispatcherTemplate<TVertex, FlipCriterionXie, TSolver>>(no);
			}
			else
			{
				cerr << "Unknown parameter \"" << args[i + 1] << "\" to option -flipcrit.";
			}
			++i;
		}
		else
		{
			cerr << "Unknown or invalid option \"" << args[i] << ".";
			return -1;
		}
	}

	cout << "Using the following parameters:" << endl;
	cout << "\tPreparation: " << (prepare ? "yes" : "no") << endl;
	cout << "\tStreaming: " << (stream ? "yes" : "no") << endl;
	cout << "\tRandomize: " << (randomize ? "yes" : "no") << endl;
#ifdef WITH_CUDA
	cout << "\tNormal Estimation: " << (estimateNormals ? "yes" : "no") << endl;
#endif
	cout << "\tEnergy Calculation: " << (calculateEnergies ? "yes" : "no") << endl;
	cout << "\tk: " << neighborCount << endl;
	cout << "\tr: " << searchRadius << endl;
	cout << "\ttheta_S: " << minSegmentVoteCertainty << endl;
	cout << "\ttheta_acc: " << minAccumulatedVoteCertainty << endl;
	cout << "\t" << orientationDispatcher->name << "'s flip criterion" << endl;
	cout << "\tMST + QPBO-I solver" << endl;

	Eigen::initParallel();

	
	auto inputName = model + ".ply";
	auto binName = model + ".bin";
	auto randomizedName = model + "RandomizedNormals.bin";
	auto orientedName = model + "WithOrientedNormals";
	auto optInput = model + "OptimizationInput.ply";
	auto outputName = model + "WithOrientedNormals.ply";

	try
	{

		std::shared_ptr<PointCloudStreamBinary<TVertex>> cloud;

		if (prepare)
		{
			cout << "Preparing file \"" << inputName << "\" for streaming..." << endl;
			BlockTimer t;
			auto id = t.newBlock("Preparation");
			t.startBlock(id);

			auto asciiCloud = PointCloudStreamAscii<TVertex>::FromFile(inputName.c_str());
			cloud = asciiCloud->PrepareForStreaming(binName.c_str());
			
			t.stopBlock(id);
			cout << t;
		}
		else
		{
			cout << "Using input file \"" << binName << "\"" << endl;
			cloud = std::make_shared<PointCloudStreamBinary<TVertex>>(binName.c_str());
		}		

		cout << "Calculating bounding box..." << endl;

		bp.calculateBoundingBoxCentroid(cloud);

		cout << "Bounding box: (" << bp.bbxMin[0] << ", " << bp.bbxMin[1] << ", " << bp.bbxMin[2] << ") - (" << bp.bbxMax[0] << ", " << bp.bbxMax[1] << ", " << bp.bbxMax[2] << ")" << endl;

#ifdef WITH_CUDA
		if (estimateNormals)
		{
			cout << "Calculating minimal slice size..." << endl;

			unsigned int sliceSize = bp.findSliceSizeForWidth(cloud, 2 * searchRadius);
			cout << "Minimal size: " << sliceSize << " vertices." << endl;

			cout << "Calculating normals..." << endl;

			auto normalName = model + "Normals.bin";
			auto cleanedName = model + "Cleaned.bin";

			cloud = ne.estimateNormals<TVertex, TVertex>(cloud, normalName.c_str(), 4 * searchRadius, sliceSize, bp.bbxMin, bp.bbxMax);			

			cloud = bp.filter<TVertex>(cloud, cleanedName.c_str(),
				[](const TVertex& v)
				{
					return v.normal.squaredNorm() > 0.5f;
				});
		}
#endif

		if (randomize)
		{
			cout << "Randomizing normal orientations..." << endl;
			FILE* randFile = fopen(randomizedName.c_str(), "wb");
			cloud->Process(4096, [&randFile](TVertex* v, size_t vertices, int threadNum, void*)
			{
				for (int i = 0; i < vertices; ++i)
					v[i].normal *= (rand() & 1) == 1 ? -1.0f : 1.0f;
				fwrite(v, sizeof(TVertex), vertices, randFile);
			}, nullptr, false);
			fclose(randFile);
			cloud = std::make_shared<PointCloudStreamBinary<TVertex>>(randomizedName.c_str());
		}

		if (writeOptInputPly)
		{
			cout << "Writing to PLY file..." << endl;
			pw.writePLYAscii(cloud, optInput.c_str());
		}

		if (calculateEnergies && boost::filesystem::file_size(binName) < 50000000)
		{
			double energy = orientationDispatcher->CalculateEnergy(cloud, searchRadius, neighborCount, bp.bbxMin, bp.bbxMax);
			std::cout << "Initial energy on plain neighbor graph: " << energy << std::endl;
		}

		cout << "Orienting normals..." << endl;

		int components = 0;
		if (stream)
			components = orientationDispatcher->OrientStreaming(cloud, orientedName.c_str(), searchRadius, (int)ceil(neighborCount / 2.0), searchRadius, 100, minSegmentVoteCertainty, minAccumulatedVoteCertainty, bp.bbxMin, bp.bbxMax);			
		else
			components = orientationDispatcher->Orient(cloud, orientedName.c_str(), searchRadius, searchRadius, neighborCount, bp.bbxMin, bp.bbxMax);

		if (components > 0)
		{
			std::string orientedNameComp0 = orientedName + "0.bin";
			auto cloud = std::make_shared<PointCloudStreamBinary<TVertex>>(orientedNameComp0.c_str());
			
			if (evaluateAccuracy)
			{
				auto cloud2 = std::make_shared<PointCloudStreamBinary<VertexPositionNormalColor>>(binName.c_str());
				int skipped1, skipped2, correct, wrong;
				bp.compareNormals(cloud, cloud2, correct, wrong, skipped1, skipped2);
				std::cout << "Correct normals: " << correct << " (" << (100.0 * correct / (wrong + correct)) << " %); wrong: " << wrong << " (" << (100.0 * wrong / (wrong + correct)) << " %); skipped from ground truth: " << skipped2 << "; skipped from result: " << skipped1 << endl;
			}

			if (writeFinalPly)
			{
				cout << "Writing to PLY file..." << endl;
				pw.writePLYAscii(cloud, outputName.c_str());
			}

			if (calculateEnergies && boost::filesystem::file_size(orientedNameComp0) < 50000000)
			{
				double energy = orientationDispatcher->CalculateEnergy(cloud, searchRadius, neighborCount, bp.bbxMin, bp.bbxMax);
				std::cout << "Final energy on plain neighbor graph: " << energy << std::endl;
			}
		}
		else
		{
			std::cout << "No connected components." << std::endl;
		}
	}
	catch (std::exception& e)
	{
		cerr << "Error: " << e.what() << endl;
	}

	return 0;
}

