#include <GlfwApp.h>
#include <SceneGraph.h>

#include <Volume/BasicShapeToVolume.h>
#include <Multiphysics/VolumeBoundary.h>

#include <Module/CalculateNorm.h>
#include <GLRenderEngine.h>
#include <GLPointVisualModule.h>
#include <ColorMapping.h>
#include <ImColorbar.h>

#include "StressAwareSPH/SASPHFluidSystem.h"
#include "ParticleSystem/MakeParticleSystem.h"
#include <BasicShapes/CubeModel.h>
#include <Samplers/ShapeSampler.h>
#include <ParticleSystem/Emitters/SquareEmitter.h>

#include"ABCExporter/ParticleWriterABC.h"





using namespace std;
using namespace dyno;

bool useVTK = false;




std::shared_ptr<SceneGraph> createScene()
{
	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
	scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));

	Real scale = 1.4f;

	auto cube1 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
	cube1->varLocation()->setValue(Vec3f(0.125, 0.125, 0.125));
	cube1->varLength()->setValue(Vec3f(0.15, 0.15, 0.15));
	cube1->graphicsPipeline()->disable();
	auto sampler1 = scn->addNode(std::make_shared<ShapeSampler<DataType3f>>());
	sampler1->varSamplingDistance()->setValue(0.005);
	sampler1->setVisible(false);
	cube1->connect(sampler1->importShape());
	auto initialParticles1 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
	sampler1->statePointSet()->promoteOuput()->connect(initialParticles1->inPoints());

	auto cube2 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
	cube2->varLocation()->setValue(Vec3f(-0.125, 0.125, 0.125));
	cube2->varLength()->setValue(Vec3f(0.15, 0.15, 0.15));
	cube2->graphicsPipeline()->disable();
	auto sampler2 = scn->addNode(std::make_shared<ShapeSampler<DataType3f>>());
	sampler2->varSamplingDistance()->setValue(0.005);
	sampler2->setVisible(false);
	cube2->connect(sampler2->importShape());
	auto initialParticles2 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
	sampler2->statePointSet()->promoteOuput()->connect(initialParticles2->inPoints());

	auto cube3 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
	cube3->varLocation()->setValue(Vec3f(0.125, 0.125, -0.125));
	cube3->varLength()->setValue(Vec3f(0.15, 0.15, 0.15));
	cube3->graphicsPipeline()->disable();
	auto sampler3 = scn->addNode(std::make_shared<ShapeSampler<DataType3f>>());
	sampler3->varSamplingDistance()->setValue(0.005);
	sampler3->setVisible(false);
	cube3->connect(sampler3->importShape());
	auto initialParticles3 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
	sampler3->statePointSet()->promoteOuput()->connect(initialParticles3->inPoints());

	auto cube4 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
	cube4->varLocation()->setValue(Vec3f(-0.125, 0.125, -0.125));
	cube4->varLength()->setValue(Vec3f(0.15, 0.15, 0.15));
	cube4->graphicsPipeline()->disable();
	auto sampler4 = scn->addNode(std::make_shared<ShapeSampler<DataType3f>>());
	sampler4->varSamplingDistance()->setValue(0.005);
	sampler4->setVisible(false);
	cube4->connect(sampler4->importShape());
	auto initialParticles4 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
	sampler4->statePointSet()->promoteOuput()->connect(initialParticles4->inPoints());

	auto fluid = scn->addNode(std::make_shared<SASPHFluidSystem<DataType3f>>());
	fluid->varReshuffleParticles()->setValue(true);
	initialParticles1->connect(fluid->importInitialStates());
	initialParticles2->connect(fluid->importInitialStates());
	initialParticles3->connect(fluid->importInitialStates());
	initialParticles4->connect(fluid->importInitialStates());

	//PBD or VSSPH
	//auto fluid = scn->addNode(std::make_shared<ParticleFluid<DataType3f>>());
	//fluid->varReshuffleParticles()->setValue(true);
	////initialParticles1->connect(fluid->importInitialStates());
	////initialParticles2->connect(fluid->importInitialStates());
	//initialParticles3->connect(fluid->importInitialStates());
	////initialParticles4->connect(fluid->importInitialStates());
	
	//Create a boundary
	auto cubeBoundary = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
	cubeBoundary->varLocation()->setValue(Vec3f(0.0f, 1.0f, 0.0f));
	cubeBoundary->varLength()->setValue(Vec3f(0.5f, 2.0f, 0.5f));
	cubeBoundary->setVisible(false);

	auto cube2vol = scn->addNode(std::make_shared<BasicShapeToVolume<DataType3f>>());
	cube2vol->varGridSpacing()->setValue(0.02f);
	cube2vol->varInerted()->setValue(true);
	cubeBoundary->connect(cube2vol->importShape());

	auto container = scn->addNode(std::make_shared<VolumeBoundary<DataType3f>>());
	cube2vol->connect(container->importVolumes());

	fluid->connect(container->importParticleSystems());

	auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
	fluid->stateVelocity()->connect(calculateNorm->inVec());
	fluid->graphicsPipeline()->pushModule(calculateNorm);

	auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
	colorMapper->varMax()->setValue(5.0f);
	calculateNorm->outNorm()->connect(colorMapper->inScalar());
	fluid->graphicsPipeline()->pushModule(colorMapper);

	auto ptRender = std::make_shared<GLPointVisualModule>();
	ptRender->setColor(Color(1, 0, 0));
	//ptRender->varPointSize()->setValue(0.0035f);
	ptRender->varPointSize()->setValue(0.001f);
	ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
	fluid->statePointSet()->connect(ptRender->inPointSet());
	colorMapper->outColor()->connect(ptRender->inColor());
	fluid->graphicsPipeline()->pushModule(ptRender);

	// A simple color bar widget for node
	auto colorBar = std::make_shared<ImColorbar>();
	colorBar->varMax()->setValue(5.0f);
	colorBar->varFieldName()->setValue("Velocity");
	calculateNorm->outNorm()->connect(colorBar->inScalar());
	// add the widget to app
	fluid->graphicsPipeline()->pushModule(colorBar);

	return scn;
}

int main()
{
	GlfwApp window;
	window.setSceneGraph(createScene());
	window.initialize(1024, 768);
	window.mainLoop();

	return 0;
}

//#include <GlfwApp.h>
//#include "SceneGraph.h"
//#include <Log.h>
//#include "ParticleSystem/StaticBoundary.h"
//#include <Module/CalculateNorm.h>
//#include <GLRenderEngine.h>
//#include <GLPointVisualModule.h>
//#include <ColorMapping.h>
//#include <ImColorbar.h>
//#include  "AKSPH/AKSPHFluidSystem.h"
//#include "ParticleSystem/MakeParticleSystem.h"
//#include <BasicShapes/CubeModel.h>
//#include <BasicShapes/SphereModel.h>
//#include <SemiAnalyticalScheme/TriangularMeshBoundary.h>
//#include <ParticleSystem/CubeSampler.h>
//#include <StaticTriangularMesh.h>
//#include <GLSurfaceVisualModule.h>
//#include <ParticleSystem/SquareEmitter.h>
//#include "PointsLoader.h"
//#include <ParticleSystem/PoissonEmitter.h>
//
//#include"ABCExporter/ParticleWriterABC.h"
////#include "ParticleSystem/GhostParticles.h"
//
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	int i = 3;
//	if (i == 1)
//	{
//		scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//		scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//		//Create a particle emitter
//		auto emitter = scn->addNode(std::make_shared<PoissonEmitter<DataType3f>>());
//		emitter->varSamplingDistance()->setValue(0.008f);
//		emitter->varEmitterShape()->getDataPtr()->setCurrentKey(1);
//		emitter->varWidth()->setValue(0.1f);
//		emitter->varHeight()->setValue(0.1f);
//		emitter->varVelocityMagnitude()->setValue(1.5);
//		//emitter->varLocation()->setValue(Vec3f(0.15f, 0.5f, 0.15f));
//		emitter->varLocation()->setValue(Vec3f(0.0f, 0.5f, 0.0f));
//
//		//Create a particle-based fluid solver
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		//fluid->loadParticles(Vec3f(0.0f), Vec3f(0.2f), 0.005f);
//		emitter->connect(fluid->importParticleEmitters());
//
//		/*auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.25, 0, -0.25), Vec3f(0.25, 2, 0.25), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());*/
//
//		/*auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.25, 0, -0.25), Vec3f(0.25, 2, 0.25), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());*/
//
//		auto ball = scn->addNode(std::make_shared<SphereModel<DataType3f>>());
//		ball->varScale()->setValue(Vec3f(0.38));
//		ball->varLocation()->setValue(Vec3f(0.0, 0.29, 0.0));
//		auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//		sRenderf->setColor(Color(0.8f, 0.52f, 0.25f));
//		sRenderf->setVisible(true);
//		sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//		ball->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//		ball->graphicsPipeline()->pushModule(sRenderf);
//
//		auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//		ball->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//		fluid->connect(pm_collide->importParticleSystems());
//		//fluid->stateVelocity()->connect(pm_collide->inVelocity());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//
//		auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(1, 0, 0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		ptRender->varPointSize()->setValue(0.0035);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//
//		colorMapper->outColor()->connect(ptRender->inColor());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//		fluid->graphicsPipeline()->pushModule(colorMapper);
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//
//		/*auto vpRender = std::make_shared<GLPointVisualModule>();
//		vpRender->setColor(Color(1, 1, 0));
//		vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		fluid->stateVirtualPointSet()->connect(vpRender->inPointSet());
//		vpRender->varPointSize()->setValue(0.0005);
//		fluid->graphicsPipeline()->pushModule(vpRender);*/
//
//		//auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//		//fluid->statePointSet()->connect(ptswriter->inPointSet());
//		//fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//		//dyno::FilePath outpath("E:/StableSPH/file__/play/");
//		//colorBar->inScalar()->connect(ptswriter->inColor());
//		//fluid->stateVelocity()->connect(ptswriter->inVelocity());
//		//ptswriter->varOutputPath()->setValue(outpath);
//		//ptswriter->varPrefix()->setValue("AKSPH");
//		//fluid->animationPipeline()->pushModule(ptswriter);
//
//	}
//	else if (i == 2)
//	{
//
//		scn->setUpperBound(Vec3f(1.5, 1, 1.5));
//		scn->setLowerBound(Vec3f(-0.5, 0, -0.5));
//
//		//Create a cube
//		auto cube = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		cube->varLocation()->setValue(Vec3f(0.6, 0.5, 0.5));
//		cube->varLength()->setValue(Vec3f(0.2, 0.2, 0.2));
//		cube->graphicsPipeline()->disable();
//
//		//Create a sampler
//		auto sampler = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//		sampler->varSamplingDistance()->setValue(0.005);
//		sampler->graphicsPipeline()->disable();
//
//		cube->outCube()->connect(sampler->inCube());
//
//		auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//
//		sampler->statePointSet()->promoteOuput()->connect(initialParticles->inPoints());
//
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		fluid->varReshuffleParticles()->setValue(true);
//		fluid->loadParticles(Vec3f(0.5, 0.2, 0.4), Vec3f(0.7, 1.5, 0.6), 0.005);
//		initialParticles->connect(fluid->importInitialStates());
//
//		//Create a boundary
//		auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>()); ;
//		boundary->loadCube(Vec3f(-0.5, 0, -0.5), Vec3f(1.5, 2, 1.5), 0.02, true);
//		boundary->loadSDF(getAssetPath() + "bowl/bowl.sdf", false);
//		fluid->connect(boundary->importParticleSystems());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//
//		auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());
//		fluid->graphicsPipeline()->pushModule(colorMapper);
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(1, 0, 0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//		colorMapper->outColor()->connect(ptRender->inColor());
//
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		// A simple color bar widget for node
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//
//		/*auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//		fluid->statePointSet()->connect(ptswriter->inPointSet());
//		fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//		dyno::FilePath outpath("E:/StableSPH/file__/play/");
//		colorBar->inScalar()->connect(ptswriter->inColor());
//		fluid->stateVelocity()->connect(ptswriter->inVelocity());
//		ptswriter->varOutputPath()->setValue(outpath);
//		ptswriter->varPrefix()->setValue("AKSPH");
//		fluid->animationPipeline()->pushModule(ptswriter);*/
//	}\
//	else if (i == 3)
//	{
//		scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//		scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//
//		//Create a particle-based fluid solver
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		//fluid->loadParticles(Vec3f(0.0f), Vec3f(0.2f), 0.005f);
//
//		auto cube1 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		//cube1->varLocation()->setValue(Vec3f(0.35, 0.1, 0));
//		//cube1->varLength()->setValue(Vec3f(0.1, 0.19, 0.1));
//		cube1->varLocation()->setValue(Vec3f(0, 0.05, 0));
//		cube1->varLength()->setValue(Vec3f(0.2, 0.09, 0.2));
//		cube1->graphicsPipeline()->disable();
//		auto sampler1 = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//		sampler1->varSamplingDistance()->setValue(0.005);
//		sampler1->setVisible(false);
//		cube1->outCube()->connect(sampler1->inCube());
//		auto initialParticles1 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//		sampler1->statePointSet()->promoteOuput()->connect(initialParticles1->inPoints());
//
//		//auto cube2 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		//cube2->varLocation()->setValue(Vec3f(0,0.0375, 0));
//		//cube2->varLength()->setValue(Vec3f(0.05, 0.075, 0.1));
//		////cube2->graphicsPipeline()->disable();
//		//auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//		//sRenderf->setColor(Color(0.8f, 0.52f, 0.25f));
//		//sRenderf->setVisible(true);
//		//sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//		//cube2->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//		//cube2->graphicsPipeline()->pushModule(sRenderf);
//
//		//auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//		//cube2->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//		//fluid->connect(pm_collide->importParticleSystems());
//		////fluid->stateVelocity()->connect(pm_collide->inVelocity());
//
//		/*auto cube2 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		cube2->varLocation()->setValue(Vec3f(-0.05, 0.05, 0));
//		cube2->varLength()->setValue(Vec3f(0.7, 0.09, 0.1));
//		cube2->graphicsPipeline()->disable();
//		auto sampler2 = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//		sampler2->varSamplingDistance()->setValue(0.005);
//		sampler2->setVisible(false);
//		cube2->outCube()->connect(sampler2->inCube());
//		auto initialParticles2 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//		sampler2->statePointSet()->promoteOuput()->connect(initialParticles2->inPoints());*/
//
//		//auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		////boundary->loadCube(Vec3f(-0.405, 0, -0.055), Vec3f(0.405, 0.24, 0.055), 0.05, true);
//		//boundary->loadCube(Vec3f(-0.125, 0, -0.055), Vec3f(0.125, 0.24, 0.055), 0.05, true);
//		//fluid->connect(boundary->importParticleSystems());
//
//		auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.110, -0.05, -0.110), Vec3f(0.110, 0.2, 0.110), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());
//
//		fluid->varReshuffleParticles()->setValue(true);
//		initialParticles1->connect(fluid->importInitialStates());
//		//initialParticles2->connect(fluid->importInitialStates());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//
//		auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(1, 0, 0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		ptRender->varPointSize()->setValue(0.0035);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//
//		colorMapper->outColor()->connect(ptRender->inColor());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//		fluid->graphicsPipeline()->pushModule(colorMapper);
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//
//		/*auto vpRender = std::make_shared<GLPointVisualModule>();
//		vpRender->setColor(Color(1, 1, 0));
//		vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		fluid->stateVirtualPointSet()->connect(vpRender->inPointSet());
//		vpRender->varPointSize()->setValue(0.0005);
//		fluid->graphicsPipeline()->pushModule(vpRender);*/
//
//		/*auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//		fluid->statePointSet()->connect(ptswriter->inPointSet());
//		fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//		dyno::FilePath outpath("E:/StableSPH/file__/play/");
//		colorBar->inScalar()->connect(ptswriter->inColor());
//		fluid->stateVelocity()->connect(ptswriter->inVelocity());
//		ptswriter->varOutputPath()->setValue(outpath);
//		ptswriter->varPrefix()->setValue("AKSPH");
//		fluid->animationPipeline()->pushModule(ptswriter);*/
//	}
//	return scn;
//}
//
//int main()
//{
//	GlfwApp window;
//	window.setSceneGraph(createScene());
//	window.initialize(1024, 768);
//	window.mainLoop();
//
//	return 0;
//}


//#include <GlfwApp.h>
//#include "SceneGraph.h"
//#include <Log.h>
//#include "ParticleSystem/StaticBoundary.h"
//#include <Module/CalculateNorm.h>
//#include <GLRenderEngine.h>
//#include <GLPointVisualModule.h>
//#include <ColorMapping.h>
//#include <ImColorbar.h>
//#include  "AKSPH/AKSPHFluidSystem.h"
//#include "ParticleSystem/MakeParticleSystem.h"
//#include <BasicShapes/CubeModel.h>
//#include <BasicShapes/SphereModel.h>
//#include <SemiAnalyticalScheme/TriangularMeshBoundary.h>
//#include <ParticleSystem/CubeSampler.h>
//#include <StaticTriangularMesh.h>
//#include <GLSurfaceVisualModule.h>
//#include <ParticleSystem/SquareEmitter.h>
//#include "PointsLoader.h"
//#include <ParticleSystem/PoissonEmitter.h>
//
//#include"ABCExporter/ParticleWriterABC.h"
////#include "ParticleSystem/GhostParticles.h"
//
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	int i = 4;
//	if (i == 1)
//	{
//		scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//		scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//		//Create a particle emitter
//		auto emitter = scn->addNode(std::make_shared<PoissonEmitter<DataType3f>>());
//		emitter->varSamplingDistance()->setValue(0.008f);
//		emitter->varEmitterShape()->getDataPtr()->setCurrentKey(1);
//		emitter->varWidth()->setValue(0.1f);
//		emitter->varHeight()->setValue(0.1f);
//		emitter->varVelocityMagnitude()->setValue(1.5);
//		//emitter->varLocation()->setValue(Vec3f(0.15f, 0.5f, 0.15f));
//		emitter->varLocation()->setValue(Vec3f(0.0f, 0.5f, 0.0f));
//
//		//Create a particle-based fluid solver
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		//fluid->loadParticles(Vec3f(0.0f), Vec3f(0.2f), 0.005f);
//		emitter->connect(fluid->importParticleEmitters());
//
//		/*auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.25, 0, -0.25), Vec3f(0.25, 2, 0.25), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());*/
//
//		/*auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.25, 0, -0.25), Vec3f(0.25, 2, 0.25), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());*/
//
//		auto ball = scn->addNode(std::make_shared<SphereModel<DataType3f>>());
//		ball->varScale()->setValue(Vec3f(0.38));
//		ball->varLocation()->setValue(Vec3f(0.0, 0.29, 0.0));
//		auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//		sRenderf->setColor(Color(0.8f, 0.52f, 0.25f));
//		sRenderf->setVisible(true);
//		sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//		ball->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//		ball->graphicsPipeline()->pushModule(sRenderf);
//
//		auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//		ball->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//		fluid->connect(pm_collide->importParticleSystems());
//		//fluid->stateVelocity()->connect(pm_collide->inVelocity());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//
//		auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(1, 0, 0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		ptRender->varPointSize()->setValue(0.0035);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//
//		colorMapper->outColor()->connect(ptRender->inColor());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//		fluid->graphicsPipeline()->pushModule(colorMapper);
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//
//		/*auto vpRender = std::make_shared<GLPointVisualModule>();
//		vpRender->setColor(Color(1, 1, 0));
//		vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		fluid->stateVirtualPointSet()->connect(vpRender->inPointSet());
//		vpRender->varPointSize()->setValue(0.0005);
//		fluid->graphicsPipeline()->pushModule(vpRender);*/
//
//		//auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//		//fluid->statePointSet()->connect(ptswriter->inPointSet());
//		//fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//		//dyno::FilePath outpath("E:/StableSPH/file__/play/");
//		//colorBar->inScalar()->connect(ptswriter->inColor());
//		//fluid->stateVelocity()->connect(ptswriter->inVelocity());
//		//ptswriter->varOutputPath()->setValue(outpath);
//		//ptswriter->varPrefix()->setValue("AKSPH");
//		//fluid->animationPipeline()->pushModule(ptswriter);
//
//	}
//	else if (i == 2)
//	{
//
//		scn->setUpperBound(Vec3f(1.5, 1, 1.5));
//		scn->setLowerBound(Vec3f(-0.5, 0, -0.5));
//
//		//Create a cube
//		auto cube = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		cube->varLocation()->setValue(Vec3f(0.6, 0.5, 0.5));
//		cube->varLength()->setValue(Vec3f(0.2, 0.2, 0.2));
//		cube->graphicsPipeline()->disable();
//
//		//Create a sampler
//		auto sampler = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//		sampler->varSamplingDistance()->setValue(0.005);
//		sampler->graphicsPipeline()->disable();
//
//		cube->outCube()->connect(sampler->inCube());
//
//		auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//
//		sampler->statePointSet()->promoteOuput()->connect(initialParticles->inPoints());
//
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		fluid->varReshuffleParticles()->setValue(true);
//		fluid->loadParticles(Vec3f(0.5, 0.2, 0.4), Vec3f(0.7, 1.5, 0.6), 0.005);
//		initialParticles->connect(fluid->importInitialStates());
//
//		//Create a boundary
//		auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>()); ;
//		boundary->loadCube(Vec3f(-0.5, 0, -0.5), Vec3f(1.5, 2, 1.5), 0.02, true);
//		//boundary->loadSDF(getAssetPath() + "bowl/bowl.sdf", false);
//		boundary->loadSDF(getAssetPath() + "fountain/mesh.sdf", false);
//		fluid->connect(boundary->importParticleSystems());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//
//		auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());
//		fluid->graphicsPipeline()->pushModule(colorMapper);
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(1, 0, 0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//		colorMapper->outColor()->connect(ptRender->inColor());
//
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		// A simple color bar widget for node
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//
//		/*auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//		fluid->statePointSet()->connect(ptswriter->inPointSet());
//		fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//		dyno::FilePath outpath("E:/StableSPH/file__/play/");
//		colorBar->inScalar()->connect(ptswriter->inColor());
//		fluid->stateVelocity()->connect(ptswriter->inVelocity());
//		ptswriter->varOutputPath()->setValue(outpath);
//		ptswriter->varPrefix()->setValue("AKSPH");
//		fluid->animationPipeline()->pushModule(ptswriter);*/
//	}
//	else if (i == 3)
//	{
//		scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//		scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//
//		//Create a particle-based fluid solver
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		//fluid->loadParticles(Vec3f(0.0f), Vec3f(0.2f), 0.005f);
//
//		auto cube1 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		//cube1->varLocation()->setValue(Vec3f(0.35, 0.1, 0));
//		//cube1->varLength()->setValue(Vec3f(0.1, 0.19, 0.1));
//		cube1->varLocation()->setValue(Vec3f(0, 0.05, 0));
//		cube1->varLength()->setValue(Vec3f(0.2, 0.09, 0.2));
//		cube1->graphicsPipeline()->disable();
//		auto sampler1 = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//		sampler1->varSamplingDistance()->setValue(0.005);
//		sampler1->setVisible(false);
//		cube1->outCube()->connect(sampler1->inCube());
//		auto initialParticles1 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//		sampler1->statePointSet()->promoteOuput()->connect(initialParticles1->inPoints());
//
//		//auto cube2 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		//cube2->varLocation()->setValue(Vec3f(0,0.0375, 0));
//		//cube2->varLength()->setValue(Vec3f(0.05, 0.075, 0.1));
//		////cube2->graphicsPipeline()->disable();
//		//auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//		//sRenderf->setColor(Color(0.8f, 0.52f, 0.25f));
//		//sRenderf->setVisible(true);
//		//sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//		//cube2->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//		//cube2->graphicsPipeline()->pushModule(sRenderf);
//
//		//auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//		//cube2->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//		//fluid->connect(pm_collide->importParticleSystems());
//		////fluid->stateVelocity()->connect(pm_collide->inVelocity());
//
//		/*auto cube2 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//		cube2->varLocation()->setValue(Vec3f(-0.05, 0.05, 0));
//		cube2->varLength()->setValue(Vec3f(0.7, 0.09, 0.1));
//		cube2->graphicsPipeline()->disable();
//		auto sampler2 = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//		sampler2->varSamplingDistance()->setValue(0.005);
//		sampler2->setVisible(false);
//		cube2->outCube()->connect(sampler2->inCube());
//		auto initialParticles2 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//		sampler2->statePointSet()->promoteOuput()->connect(initialParticles2->inPoints());*/
//
//		//auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		////boundary->loadCube(Vec3f(-0.405, 0, -0.055), Vec3f(0.405, 0.24, 0.055), 0.05, true);
//		//boundary->loadCube(Vec3f(-0.125, 0, -0.055), Vec3f(0.125, 0.24, 0.055), 0.05, true);
//		//fluid->connect(boundary->importParticleSystems());
//
//		auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.105, 0.0, -0.105), Vec3f(0.105, 0.2, 0.105), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());
//
//		fluid->varReshuffleParticles()->setValue(true);
//		initialParticles1->connect(fluid->importInitialStates());
//		//initialParticles2->connect(fluid->importInitialStates());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//
//		auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(1, 0, 0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		ptRender->varPointSize()->setValue(0.0035);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//
//		colorMapper->outColor()->connect(ptRender->inColor());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//		fluid->graphicsPipeline()->pushModule(colorMapper);
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//
//		/*auto vpRender = std::make_shared<GLPointVisualModule>();
//		vpRender->setColor(Color(1, 1, 0));
//		vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		fluid->stateVirtualPointSet()->connect(vpRender->inPointSet());
//		vpRender->varPointSize()->setValue(0.0005);
//		fluid->graphicsPipeline()->pushModule(vpRender);*/
//
//		/*auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//		fluid->statePointSet()->connect(ptswriter->inPointSet());
//		fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//		dyno::FilePath outpath("E:/StableSPH/file__/play/");
//		colorBar->inScalar()->connect(ptswriter->inColor());
//		fluid->stateVelocity()->connect(ptswriter->inVelocity());
//		ptswriter->varOutputPath()->setValue(outpath);
//		ptswriter->varPrefix()->setValue("AKSPH");
//		fluid->animationPipeline()->pushModule(ptswriter);*/
//	}
//	else if (i == 4)
//	{
//		scn->setUpperBound(Vec3f(3.0, 30.0, 3.0));
//		scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//		auto ptsLoader = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader->varLocation()->setValue(Vec3f(0.0f, 0.2f, 0.0f));
//		//ptsLoader->varScale()->setValue(Vec3f(0.5f));
//		auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader->outPointSet()->promoteOuput()->connect(initialParticles->inPoints());
//
//		auto ptsLoader2 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader2->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader2->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader2->varLocation()->setValue(Vec3f(0.30f, 0.3f, 0.30f));
//		auto initialParticles2 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles2->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader2->outPointSet()->promoteOuput()->connect(initialParticles2->inPoints());
//
//		auto ptsLoader3 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader3->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader3->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader3->varLocation()->setValue(Vec3f(0.30f, 0.4f, -0.30f));
//		auto initialParticles3 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles3->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader3->outPointSet()->promoteOuput()->connect(initialParticles3->inPoints());
//
//		auto ptsLoader4 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader4->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader4->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader4->varLocation()->setValue(Vec3f(-0.30f, 0.5f, -0.30f));
//		auto initialParticles4 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles4->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader4->outPointSet()->promoteOuput()->connect(initialParticles4->inPoints());
//
//		auto ptsLoader5 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader5->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader5->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader5->varLocation()->setValue(Vec3f(-0.30f, 0.6f, 0.30f));
//		auto initialParticles5 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles5->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader5->outPointSet()->promoteOuput()->connect(initialParticles5->inPoints());
//
//		auto ptsLoader6 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader6->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader6->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader6->varLocation()->setValue(Vec3f(0.10f, 0.7f, -0.20f));
//		auto initialParticles6 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles6->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader6->outPointSet()->promoteOuput()->connect(initialParticles6->inPoints());
//
//		auto ptsLoader7 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader7->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader7->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader7->varLocation()->setValue(Vec3f(-0.15f, 0.8f, 0.25f));
//		auto initialParticles7 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles7->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader7->outPointSet()->promoteOuput()->connect(initialParticles7->inPoints());
//
//		auto ptsLoader8 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//		ptsLoader8->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//		ptsLoader8->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
//		ptsLoader8->varLocation()->setValue(Vec3f(0.2f, 0.9f, 0.05f));
//		auto initialParticles8 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//		initialParticles8->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 0.0f));
//		ptsLoader8->outPointSet()->promoteOuput()->connect(initialParticles8->inPoints());
//
//		//Create a particle-based fluid solver
//		auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//		fluid->varReshuffleParticles()->setValue(true);
//		initialParticles->connect(fluid->importInitialStates());
//		initialParticles2->connect(fluid->importInitialStates());
//		initialParticles3->connect(fluid->importInitialStates());
//		initialParticles4->connect(fluid->importInitialStates());
//		initialParticles5->connect(fluid->importInitialStates());
//		initialParticles6->connect(fluid->importInitialStates());
//		initialParticles7->connect(fluid->importInitialStates());
//		initialParticles8->connect(fluid->importInitialStates());
//
//		auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//		boundary->loadCube(Vec3f(-0.5, 0, -0.5), Vec3f(0.5, 30.0, 0.5), 0.05, true);
//		fluid->connect(boundary->importParticleSystems());
//
//		//Create a foutain boundary
//		//auto fountain = scn->addNode(std::make_shared<StaticTriangularMesh<DataType3f>>());
//		//fountain->varFileName()->setValue(getAssetPath() + "gargoyle/gargoyle1.0.obj");
//		//fountain->varLocation()->setValue(Vec3f(0.0f, 0.0f, 0.0f));//obj 中心位置
//		//fountain->varScale()->setValue(Vec3f(1.0f));//放缩模型大小
//		//auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//		//sRenderf->setColor(Color(0.43f, 0.5f, 0.56f));
//		//sRenderf->setVisible(true);
//		//sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//		//fountain->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//		//fountain->graphicsPipeline()->pushModule(sRenderf);
//		////创建碰撞
//		//auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//		//fountain->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//		//fluid->connect(pm_collide->importParticleSystems());
//
//		auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//		fluid->stateVelocity()->connect(calculateNorm->inVec());
//
//		/*auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//		colorMapper->varMax()->setValue(5.0f);
//		calculateNorm->outNorm()->connect(colorMapper->inScalar());*/
//
//		auto ptRender = std::make_shared<GLPointVisualModule>();
//		ptRender->setColor(Color(0.0, 1.0, 1.0));
//		ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//		ptRender->varPointSize()->setValue(0.0035);
//
//		fluid->statePointSet()->connect(ptRender->inPointSet());
//
//		//colorMapper->outColor()->connect(ptRender->inColor());
//		fluid->graphicsPipeline()->pushModule(calculateNorm);
//		//fluid->graphicsPipeline()->pushModule(colorMapper);
//		fluid->graphicsPipeline()->pushModule(ptRender);
//
//		auto colorBar = std::make_shared<ImColorbar>();
//		colorBar->varMax()->setValue(5.0f);
//		colorBar->varFieldName()->setValue("Velocity");
//		calculateNorm->outNorm()->connect(colorBar->inScalar());
//		// add the widget to app
//		fluid->graphicsPipeline()->pushModule(colorBar);
//	}
//	return scn;
//}
//
//int main()
//{
//	GlfwApp window;
//	window.setSceneGraph(createScene());
//	window.initialize(1024, 768);
//	window.mainLoop();
//
//	return 0;
//}




//#include <GlfwApp.h>
////#include <QtApp.h>
//#include "SceneGraph.h"
//#include <Log.h>
//#include "ParticleSystem/StaticBoundary.h"
//#include <Module/CalculateNorm.h>
//#include <GLRenderEngine.h>
//#include <GLPointVisualModule.h>
//#include <ColorMapping.h>
//#include <ImColorbar.h>
////#include "DualParticleSystem/DualParticleFluidSystem.h"
//#include "AKSPH/AKSPHFluidSystem.h"
//#include "ParticleSystem/MakeParticleSystem.h"
//#include <BasicShapes/SphereModel.h>
//#include <SemiAnalyticalScheme/TriangularMeshBoundary.h>
//#include <StaticTriangularMesh.h>
//#include <GLSurfaceVisualModule.h>
//#include "Auxiliary/DataSource.h"
//#include "PointsLoader.h"
//
//#include "ParticleWriter.h"
//#include"ABCExporter/ParticleWriterABC.h"
//
//#include <ParticleSystem/ParticleFluid.h>
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	//scn->setGravity(Vec3f(0.0, 0.0, 0.0));
//	scn->setUpperBound(Vec3f(3.0, 3, 3.0));
//	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//	const Real R = 0.5f;
//	const Real theta = M_PI / 6;
//	const Real phi = 2 * M_PI / 5;
//	const Real v = 2.0f;
//
//	const Vec3f pos_h = Vec3f(0.0, 0.7, 0.0);
//
//	auto ptsLoader = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	const Real pitch = 0.0f;
//	const Real yaw = theta;
//	const Real roll = 3 * M_PI / 2 - phi;
//	ptsLoader->varRotation()->setValue(Vec3f(pitch, yaw, roll));
//	ptsLoader->varLocation()->setValue(R * Vec3f(cos(theta) * cos(phi), sin(theta), cos(theta) * sin(phi)) + pos_h);
//	auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles->varInitialVelocity()->setValue(-v * Vec3f(cos(theta) * cos(phi), sin(theta), cos(theta) * sin(phi)));
//	ptsLoader->outPointSet()->promoteOuput()->connect(initialParticles->inPoints());
//
//	const Real theta2 = M_PI / 6;
//	const Real phi2 = 4 * M_PI / 5;
//
//	auto ptsLoader2 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader2->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	const Real pitch2 = 0.0f;
//	const Real yaw2 = theta2;
//	const Real roll2 = 3 * M_PI / 2 - phi2;
//	ptsLoader2->varRotation()->setValue(Vec3f(pitch2, yaw2, roll2));
//	ptsLoader2->varLocation()->setValue(R * Vec3f(cos(theta2) * cos(phi2), sin(theta2), cos(theta2) * sin(phi2)) + pos_h);
//	auto initialParticles2 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles2->varInitialVelocity()->setValue(-v * Vec3f(cos(theta2) * cos(phi2), sin(theta2), cos(theta2) * sin(phi2)));
//	ptsLoader2->outPointSet()->promoteOuput()->connect(initialParticles2->inPoints());
//
//	const Real theta3 = M_PI / 6;
//	const Real phi3 = 6 * M_PI / 5;
//
//	auto ptsLoader3 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader3->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	const Real pitch3 = 0.0f;
//	const Real yaw3 = theta3;
//	const Real roll3 = 3 * M_PI / 2 - phi3;
//	ptsLoader3->varRotation()->setValue(Vec3f(pitch3, yaw3, roll3));
//	ptsLoader3->varLocation()->setValue(R * Vec3f(cos(theta3) * cos(phi3), sin(theta3), cos(theta3) * sin(phi3)) + pos_h);
//	auto initialParticles3 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles3->varInitialVelocity()->setValue(-v * Vec3f(cos(theta3) * cos(phi3), sin(theta3), cos(theta3) * sin(phi3)));
//	ptsLoader3->outPointSet()->promoteOuput()->connect(initialParticles3->inPoints());
//
//	const Real theta4 = M_PI / 6;
//	const Real phi4 = 8 * M_PI / 5;
//
//	auto ptsLoader4 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader4->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	const Real pitch4 = 0.0f;
//	const Real yaw4 = theta4;
//	const Real roll4 = 3 * M_PI / 2 - phi4;
//	ptsLoader4->varRotation()->setValue(Vec3f(pitch4, yaw4, roll4));
//	ptsLoader4->varLocation()->setValue(R * Vec3f(cos(theta4) * cos(phi4), sin(theta4), cos(theta4) * sin(phi4)) + pos_h);
//	auto initialParticles4 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles4->varInitialVelocity()->setValue(-v * Vec3f(cos(theta4) * cos(phi4), sin(theta4), cos(theta4) * sin(phi4)));
//	ptsLoader4->outPointSet()->promoteOuput()->connect(initialParticles4->inPoints());
//
//	const Real theta5 = M_PI / 6;
//	const Real phi5 = 0.0f;
//
//	auto ptsLoader5 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader5->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	const Real pitch5 = 0.0f;
//	const Real yaw5 = theta5;
//	const Real roll5 = 3 * M_PI / 2 - phi5;
//	ptsLoader5->varRotation()->setValue(Vec3f(pitch5, yaw5, roll5));
//	ptsLoader5->varLocation()->setValue(R * Vec3f(cos(theta5) * cos(phi5), sin(theta5), cos(theta5) * sin(phi5)) + pos_h);
//	auto initialParticles5 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles5->varInitialVelocity()->setValue(-v * Vec3f(cos(theta5) * cos(phi5), sin(theta5), cos(theta5) * sin(phi5)));
//	ptsLoader5->outPointSet()->promoteOuput()->connect(initialParticles5->inPoints());
//
//	//Create a particle-based fluid solver
//	auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//	fluid->varReshuffleParticles()->setValue(true);
//	initialParticles->connect(fluid->importInitialStates());
//	initialParticles2->connect(fluid->importInitialStates());
//	initialParticles3->connect(fluid->importInitialStates());
//	initialParticles4->connect(fluid->importInitialStates());
//	initialParticles5->connect(fluid->importInitialStates());
//
//	auto armadillo = scn->addNode(std::make_shared<StaticTriangularMesh<DataType3f>>());
//	armadillo->varFileName()->setValue(getAssetPath() + "armadillo/armadillo.obj");
//	armadillo->varLocation()->setValue(Vec3f(-1.0f, -0.3f, -1.0f));//obj 中心位置
//	armadillo->varScale()->setValue(Vec3f(2.0f));//放缩模型大小
//	auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//	sRenderf->setColor(Color(0.43f, 0.5f, 0.56f));
//	sRenderf->setVisible(true);
//	sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//	armadillo->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//	armadillo->graphicsPipeline()->pushModule(sRenderf);
//
//	//创建碰撞
//	auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//	armadillo->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//	fluid->connect(pm_collide->importParticleSystems());
//
//	auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
//	fluid->stateVelocity()->connect(calculateNorm->inVec());
//	fluid->graphicsPipeline()->pushModule(calculateNorm);
//
//	auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
//	colorMapper->varMax()->setValue(5.0f);
//	calculateNorm->outNorm()->connect(colorMapper->inScalar());
//	fluid->graphicsPipeline()->pushModule(colorMapper);
//
//	auto ptRender = std::make_shared<GLPointVisualModule>();
//	ptRender->setColor(Color(1, 0, 0));
//	ptRender->varPointSize()->setValue(0.0035f);
//	ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//
//	fluid->statePointSet()->connect(ptRender->inPointSet());
//	colorMapper->outColor()->connect(ptRender->inColor());
//	fluid->graphicsPipeline()->pushModule(ptRender);
//
//	// A simple color bar widget for node
//	auto colorBar = std::make_shared<ImColorbar>();
//	colorBar->varMax()->setValue(5.0f);
//	colorBar->varFieldName()->setValue("Velocity");
//	calculateNorm->outNorm()->connect(colorBar->inScalar());
//	// add the widget to app
//	fluid->graphicsPipeline()->pushModule(colorBar);
//
//
//	/*auto vpRender = std::make_shared<GLPointVisualModule>();
//	vpRender->setColor(Color(1, 1, 0));
//	vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//	fluid->stateVirtualPointSet()->connect(vpRender->inPointSet());
//	vpRender->varPointSize()->setValue(0.0005);
//	fluid->graphicsPipeline()->pushModule(vpRender);*/
//
//	//auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//	//fluid->statePointSet()->connect(ptswriter->inPointSet());
//	//fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//	//dyno::FilePath outpath("");
//	//ptswriter->varOutputPath()->setValue(outpath);
//	//ptswriter->inColor()->connect(colorMapper->outColor());
//	//ptswriter->varPrefix()->setValue("armadillowater");
//	//fluid->animationPipeline()->pushModule(ptswriter);
//	
//
//	auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//	fluid->statePointSet()->connect(ptswriter->inPointSet());
//	fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//	dyno::FilePath outpath("D:/Coding/output");
//	//calculateNorm->outNorm()->connect(ptswriter->inColor());
//	fluid->stateVelocity()->connect(ptswriter->inVelocity());
//	ptswriter->varOutputPath()->setValue(outpath);
//	ptswriter->varPrefix()->setValue("armadillowater");
//	fluid->animationPipeline()->pushModule(ptswriter);
//
//
//	return scn;
//}
//
//int main()
//{
//
//	GlfwApp window;
//	//QtApp window;
//	window.setSceneGraph(createScene());
//	window.initialize(1024, 768);
//	window.mainLoop();
//
//	return 0;
//}


