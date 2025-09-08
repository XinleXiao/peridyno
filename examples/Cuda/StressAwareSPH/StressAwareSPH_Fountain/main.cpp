
//#include <GlfwApp.h>
//#include "SceneGraph.h"
//#include <Log.h>
//#include "ParticleSystem/StaticBoundary.h"
//#include <Module/CalculateNorm.h>
//#include <GLRenderEngine.h>
//#include <GLPointVisualModule.h>
//#include <ColorMapping.h>
//#include <ImColorbar.h>
////#include "DualParticleSystem/DualParticleFluidSystem.h"
//
//#include"AKSPH/AKSPHFluidSystem.h"
//#include "ParticleSystem/MakeParticleSystem.h"
//#include <BasicShapes/CubeModel.h>
//#include <ParticleSystem/CubeSampler.h>
//#include <ParticleSystem/SquareEmitter.h>
//#include "PointsLoader.h"
//
//#include <ParticleSystem/ParticleFluid.h>
//
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//	scn->setGravity(Vec3f(0.0f));
//
//	auto ptsLoader = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	ptsLoader->varRotation()->setValue(Vec3f(0.0f, 0.0f, 3.1415926f));
//	ptsLoader->varLocation()->setValue(Vec3f(0.0f, 0.0f, 0.23f));
//	auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, -1.5f));
//	ptsLoader->outPointSet()->promoteOuput()->connect(initialParticles->inPoints());
//
//	auto ptsLoader2 = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
//	ptsLoader2->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
//	ptsLoader2->varRotation()->setValue(Vec3f(0.0f, 0.0f, 0.0));
//	ptsLoader2->varLocation()->setValue(Vec3f(0.0f, 0.0f, -0.23f));
//	auto initialParticles2 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
//	initialParticles2->varInitialVelocity()->setValue(Vec3f(0.0f, 0.0f, 1.5f));
//	ptsLoader2->outPointSet()->promoteOuput()->connect(initialParticles2->inPoints());
//
//	auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//	fluid->varReshuffleParticles()->setValue(true);
//	initialParticles->connect(fluid->importInitialStates());
//	initialParticles2->connect(fluid->importInitialStates());
//
//	//PBD
//	//auto fluid = scn->addNode(std::make_shared<ParticleFluid<DataType3f>>());
//	//fluid->varReshuffleParticles()->setValue(true);
//	//initialParticles->connect(fluid->importInitialStates());
//	//initialParticles2->connect(fluid->importInitialStates());
//
//	//Create a boundary
//	auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>()); ;
//	//boundary->loadCube(Vec3f(-0.5, 0, -0.5), Vec3f(1.5, 2, 1.5), 0.02, true);
//	//boundary->loadSDF(getAssetPath() + "bowl/bowl.sdf", false);
//	fluid->connect(boundary->importParticleSystems());
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
//	//auto vpRender = std::make_shared<GLPointVisualModule>();
//	//vpRender->setColor(Color(1, 1, 0));
//	//vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
//	//fluid->stateBoundaryVirtualPosition()->connect(vpRender->inPointSet());
//	//vpRender->varPointSize()->setValue(0.0005);
//	//fluid->graphicsPipeline()->pushModule(vpRender);
//
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
//#include  "AKSPH/AKSPHFluidSystem.h"
//#include"AKSPH/StableFluidSystem.h"
//#include "ParticleSystem/MakeParticleSystem.h"
//#include <BasicShapes/CubeModel.h>
//#include <ParticleSystem/CubeSampler.h>
//#include <ParticleSystem/SquareEmitter.h>
//#include "PointsLoader.h"
//#include <ParticleSystem/PoissonEmitter.h>
//
//#include <ParticleSystem/ParticleFluid.h>
//#include <ParticleSystem/Module/ProjectionBasedFluidModel.h>
//
//
//#include <SemiAnalyticalScheme/TriangularMeshBoundary.h>
//#include <StaticTriangularMesh.h>
//#include <GLSurfaceVisualModule.h>
//
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//	auto emitter = scn->addNode(std::make_shared<PoissonEmitter<DataType3f>>());
//	emitter->varRotation()->setValue(Vec3f(0.0f, 0.0f, -0.0));
//	emitter->varSamplingDistance()->setValue(0.008f);
//	emitter->varEmitterShape()->getDataPtr()->setCurrentKey(1);
//	emitter->varWidth()->setValue(0.1f);
//	emitter->varHeight()->setValue(0.1f);
//	emitter->varVelocityMagnitude()->setValue(1.0f);
//	emitter->varLocation()->setValue(Vec3f(0.0f, 0.5f, 0.0f));
//
//	auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//	//auto fluid = scn->addNode(std::make_shared<DualParticleFluidSystem<DataType3f>>());
//	fluid->varReshuffleParticles()->setValue(true);
//	emitter->connect(fluid->importParticleEmitters());
//
//
//	//Create a foutain boundary
//	auto fountain = scn->addNode(std::make_shared<StaticTriangularMesh<DataType3f>>());
//	fountain->varFileName()->setValue(getAssetPath() + "bowl/bowl.obj");
//	fountain->varLocation()->setValue(Vec3f(-0.25f, 0.01f, -0.25f));//obj 中心位置
//	fountain->varScale()->setValue(Vec3f(0.5f));//放缩模型大小
//
//	auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//	sRenderf->setColor(Color(0.43f, 0.5f, 0.56f));
//	sRenderf->setVisible(true);
//	sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//	fountain->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//	fountain->graphicsPipeline()->pushModule(sRenderf);
//	//创建碰撞
//	auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//	fountain->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//	fluid->connect(pm_collide->importParticleSystems());
//
//	auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//	//boundary->loadCube(Vec3f(0.045, 0.045, 0.045), Vec3f(0.25, 0.25, 0.25), 0.02, true);
//	boundary->loadCube(Vec3f(-0.25, 0, -0.25), Vec3f(0.25, 2, 0.25), 0.02, true);
//	fluid->connect(boundary->importParticleSystems());
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
//	ptRender->varPointSize()->setValue(0.0015f);
//	ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
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
//	return scn;
//}
//
//int main()
//{
//
//	GlfwApp window;
//	window.setSceneGraph(createScene());
//	window.initialize(1024, 768);
//	window.mainLoop();
//
//	return 0;
//}


#include <GlfwApp.h>
//#include <QtApp.h>
#include "SceneGraph.h"
#include <Log.h>
#include <Module/CalculateNorm.h>
#include <GLRenderEngine.h>
#include <GLPointVisualModule.h>
#include <ColorMapping.h>
#include <ImColorbar.h>
//#include "DualParticleSystem/DualParticleFluidSystem.h"
#include "StressAwareSPH/SASPHFluidSystem.h"
#include "ParticleSystem/MakeParticleSystem.h"
#include <BasicShapes/SphereModel.h>
#include <SemiAnalyticalScheme/TriangularMeshBoundary.h>
//#include <StaticTriangularMesh.h>
#include <GLSurfaceVisualModule.h>
#include "Auxiliary/DataSource.h"
#include "PointsLoader.h"

#include <ParticleSystem/ParticleFluid.h>
using namespace std;
using namespace dyno;

bool useVTK = false;

std::shared_ptr<SceneGraph> createScene()
{
	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
	scn->setUpperBound(Vec3f(3.0, 3, 3.0));
	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));

	auto ptsLoader = scn->addNode(std::make_shared<PointsLoader<DataType3f>>());
	ptsLoader->varFileName()->setValue(getAssetPath() + "fish/FishPoints.obj");
	ptsLoader->varRotation()->setValue(Vec3f(0.0f, 3.14 * 2 / 5, 0.0f));
	ptsLoader->varLocation()->setValue(Vec3f(0.0f, 0.4f, 0.30f));
	auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f >>());
	ptsLoader->outPointSet()->promoteOuput()->connect(initialParticles->inPoints());

	auto fluid = scn->addNode(std::make_shared<SASPHFluidSystem<DataType3f>>());
	fluid->varReshuffleParticles()->setValue(true);
	initialParticles->connect(fluid->importInitialStates());

	//PBD
	//auto fluid = scn->addNode(std::make_shared<ParticleFluid<DataType3f>>());
	//fluid->varReshuffleParticles()->setValue(true);
	//initialParticles->connect(fluid->importInitialStates());

	//Create a boundary
	
	auto ball = scn->addNode(std::make_shared<SphereModel<DataType3f>>());
	ball->varScale()->setValue(Vec3f(0.38));
	ball->varLocation()->setValue(Vec3f(0.0, 0.0, 0.3));
	auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
	sRenderf->setColor(Color(0.8f, 0.52f, 0.25f));
	sRenderf->setVisible(true);
	sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
	ball->stateTriangleSet()->connect(sRenderf->inTriangleSet());
	ball->graphicsPipeline()->pushModule(sRenderf);

	auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
	ball->stateTriangleSet()->connect(pm_collide->inTriangleSet());
	fluid->connect(pm_collide->importParticleSystems());
	//fluid->stateVelocity()->connect(pm_collide->inVelocity());

	auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
	fluid->stateVelocity()->connect(calculateNorm->inVec());
	fluid->graphicsPipeline()->pushModule(calculateNorm);

	auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
	colorMapper->varMax()->setValue(5.0f);
	calculateNorm->outNorm()->connect(colorMapper->inScalar());
	fluid->graphicsPipeline()->pushModule(colorMapper);

	auto ptRender = std::make_shared<GLPointVisualModule>();
	ptRender->setColor(Color(1, 0, 0));
	ptRender->varPointSize()->setValue(0.0035f);
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


	/*auto vpRender = std::make_shared<GLPointVisualModule>();
	vpRender->setColor(Color(1, 1, 0));
	vpRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
	fluid->stateVirtualPointSet()->connect(vpRender->inPointSet());
	vpRender->varPointSize()->setValue(0.0005);
	fluid->graphicsPipeline()->pushModule(vpRender);*/

	return scn;
}

int main()
{

	GlfwApp window;
	//QtApp window;
	window.setSceneGraph(createScene());
	window.initialize(1024, 768);
	window.mainLoop();

	return 0;
}



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
//#include <BasicShapes/CubeModel.h>
//#include <ParticleSystem/CubeSampler.h>
//#include <RigidBody/RigidBodySystem.h>
//#include "RigidBody/Module/PJSConstraintSolver.h"
//#include <Mapping/DiscreteElementsToTriangleSet.h>
//
//#include <ParticleSystem/ParticleFluid.h>
//#include <Volume/Volume.h>
//#include <Volume/VolumeGenerator.h>
//
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	scn->setUpperBound(Vec3f(3.0, 3.0, 3.0));
//	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//	auto staticfluid = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//	staticfluid->varLocation()->setValue(Vec3f(0.0, 0.101, 0.0));
//	staticfluid->varLength()->setValue(Vec3f(0.5 - 0.01, 0.20, 0.5 - 0.01));
//	staticfluid->graphicsPipeline()->disable();
//	auto sampler2 = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//	sampler2->varSamplingDistance()->setValue(0.005);
//	sampler2->setVisible(false);
//	staticfluid->outCube()->connect(sampler2->inCube());
//	auto initialParticles = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//	sampler2->statePointSet()->promoteOuput()->connect(initialParticles->inPoints());
//
//	auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//	fluid->varReshuffleParticles()->setValue(true);
//	initialParticles->connect(fluid->importInitialStates());
//
//	//Create a boundary
//	auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//	boundary->loadCube(Vec3f(-0.25, 0, -0.25), Vec3f(0.25, 2, 0.25), 0.02, true);
//	
//
//
//	//设置固体
//	auto cube = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//	cube->varLocation()->setValue(Vec3f(0.0, 0.20, -2.0));
//	cube->varLength()->setValue(Vec3f(0.1, 0.1, 0.1));	
//
//	//auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//	//sRenderf->setColor(Color(0.8f, 0.52f, 0.25f));
//	//sRenderf->setVisible(true);
//	//cube->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//	//cube->graphicsPipeline()->pushModule(sRenderf);
//
//	//auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//	//cube->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//	//fluid->connect(pm_collide->importParticleSystems());
//	////fluid->stateVelocity()->connect(pm_collide->inVelocity());
//
//	fluid->connect(boundary->importParticleSystems());
//
//
//	//auto volGen = scn->addNode(std::make_shared<VolumeGenerator<DataType3f>>());
//	//cube->stateTriangleSet()->connect(volGen->inTriangleSet());
//	
//
//	//auto volTranslator = std::make_shared<MoveVolume<DataType3f>>();
//	//volGen->outGenSDF()->connect(volTranslator->inLevelSet());
//	//volGen->animationPipeline()->pushModule(volTranslator);
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
//	return scn;
//}

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
//
//#include "ParticleSystem/StaticBoundary.h"
//#include <Module/CalculateNorm.h>
//#include <GLRenderEngine.h>
//#include <GLPointVisualModule.h>
//#include <ColorMapping.h>
//#include <ImColorbar.h>
//
//#include "AKSPH/AKSPHFluidSystem.h"
//#include "ParticleSystem/MakeParticleSystem.h"
//#include <BasicShapes/CubeModel.h>
//#include <ParticleSystem/CubeSampler.h>
//#include <SemiAnalyticalScheme/TriangularMeshBoundary.h>
//#include <StaticTriangularMesh.h>
//#include <GLSurfaceVisualModule.h>
//#include "Auxiliary/DataSource.h"
//#include "PointsLoader.h"
//
//#include <ParticleSystem/ParticleFluid.h>
//
//using namespace std;
//using namespace dyno;
//
//bool useVTK = false;
//
//std::shared_ptr<SceneGraph> createScene()
//{
//	std::shared_ptr<SceneGraph> scn = std::make_shared<SceneGraph>();
//	//scn->setGravity(Vec3f(0));
//	scn->setUpperBound(Vec3f(3.0, 3, 3.0));
//	scn->setLowerBound(Vec3f(-3.0, -3.0, -3.0));
//
//	auto cube1 = scn->addNode(std::make_shared<CubeModel<DataType3f>>());
//	cube1->varLocation()->setValue(Vec3f(0, 0.5, 0.0));
//	cube1->varLength()->setValue(Vec3f(0.5, 1.0, 0.5));
//	cube1->graphicsPipeline()->disable();
//
//	auto sampler1 = scn->addNode(std::make_shared<CubeSampler<DataType3f>>());
//	sampler1->varSamplingDistance()->setValue(0.005);
//	sampler1->setVisible(false);
//	cube1->outCube()->connect(sampler1->inCube());
//	auto initialParticles1 = scn->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
//	sampler1->statePointSet()->promoteOuput()->connect(initialParticles1->inPoints());
//
//	//auto fluid = scn->addNode(std::make_shared<AKSPHFluidSystem<DataType3f>>());
//	auto fluid = scn->addNode(std::make_shared<ParticleFluid<DataType3f>>());
//	fluid->varReshuffleParticles()->setValue(true);
//	initialParticles1->connect(fluid->importInitialStates());
//
//	auto fountain = scn->addNode(std::make_shared<StaticTriangularMesh<DataType3f>>());
//	fountain->varFileName()->setValue(getAssetPath() + "bunny/bunny.obj");
//	fountain->varLocation()->setValue(Vec3f(0.0, 0.025, -0.1));//obj 中心位置
//	fountain->varScale()->setValue(Vec3f(0.01f));//放缩模型大小
//	//渲染喷泉
//	auto sRenderf = std::make_shared<GLSurfaceVisualModule>();
//	sRenderf->setColor(Color(0.43f, 0.5f, 0.56f));
//	sRenderf->setVisible(true);
//	sRenderf->varUseVertexNormal()->setValue(true);	// use generated smooth normal
//	fountain->stateTriangleSet()->connect(sRenderf->inTriangleSet());
//	fountain->graphicsPipeline()->pushModule(sRenderf);
//	//创建碰撞
//	auto pm_collide = scn->addNode(std::make_shared <TriangularMeshBoundary<DataType3f >>());
//	fountain->stateTriangleSet()->connect(pm_collide->inTriangleSet());
//	fluid->connect(pm_collide->importParticleSystems());
//
//	auto boundary = scn->addNode(std::make_shared<StaticBoundary<DataType3f>>());
//	boundary->loadCube(Vec3f(-0.25, 0, -0.5), Vec3f(0.25, 2, 1.0), 0.02, true);
//	fluid->connect(boundary->importParticleSystems());
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
//	//auto ptswriter = std::make_shared<ParticleWriterABC<DataType3f>>();
//	//fluid->statePointSet()->connect(ptswriter->inPointSet());
//	//fluid->stateFrameNumber()->connect(ptswriter->inFrameNumber());
//	//dyno::FilePath outpath("E:/StableSPH/file__/Jets/");
//	//colorBar->inScalar()->connect(ptswriter->inColor());
//	//fluid->stateVelocity()->connect(ptswriter->inVelocity());
//	//ptswriter->varOutputPath()->setValue(outpath);
//	//ptswriter->varPrefix()->setValue("AKSPH");
//	//fluid->animationPipeline()->pushModule(ptswriter);
//
//	return scn;
//}
//
//int main()
//{
//
//	GlfwApp window;
//	window.setSceneGraph(createScene());
//	window.initialize(1024, 768);
//	window.mainLoop();
//
//	return 0;
//}



