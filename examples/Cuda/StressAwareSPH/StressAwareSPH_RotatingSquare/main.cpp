#include <GlfwApp.h>
#include "SceneGraph.h"
#include <Log.h>

#include <Module/CalculateNorm.h>
#include <GLRenderEngine.h>
#include <GLPointVisualModule.h>
#include <ColorMapping.h>
#include <ImColorbar.h>
#include"StressAwareSPH/SASPHFluidSystem.h"
#include "ParticleSystem/MakeParticleSystem.h"
#include <BasicShapes/CubeModel.h>
#include <Samplers/ShapeSampler.h>
#include <ParticleSystem/Emitters/SquareEmitter.h>
#include <StaticMeshLoader.h>
#include <GLSurfaceVisualModule.h>
#include "Collision/Attribute.h"
#include "ParticleSystem/Module/ImplicitViscosity.h"
#include "RotatingSquarePatchModule.h"
#include "Auxiliary/DataSource.h"
#include"Collision/NeighborPointQuery.h"

#include"StressAwareSPH/SA_DFSPHsolver.h"
#include "ParticleSystem/Module/IterativeDensitySolver.h"
#include "ParticleSystem/Module/VariationalApproximateProjection.h"



using namespace std;
using namespace dyno;

bool useVTK = false;

std::shared_ptr<SceneGraph> createScene()
{
	std::shared_ptr<SceneGraph> scene = std::make_shared<SceneGraph>();
	scene->setGravity(Vec3f(0));
	scene->setUpperBound(Vec3f(3.0));
	scene->setLowerBound(Vec3f(-3.0));

	//Create a cube
	auto cube = scene->addNode(std::make_shared<CubeModel<DataType3f>>());
	cube->varLocation()->setValue(Vec3f(0.5, 0.5, 0.0));
	cube->varLength()->setValue(Vec3f(0.24, 0.001, 0.24));
	cube->graphicsPipeline()->disable();

	//Create a sampler
	auto sampler = scene->addNode(std::make_shared<ShapeSampler<DataType3f>>());
	sampler->varSamplingDistance()->setValue(0.005);
	sampler->setVisible(false);
	cube->connect(sampler->importShape());
	auto initialParticles = scene->addNode(std::make_shared<MakeParticleSystem<DataType3f>>());
	sampler->statePointSet()->promoteOuput()->connect(initialParticles->inPoints());



	auto fluid = scene->addNode(std::make_shared<SASPHFluidSystem<DataType3f>>());
	initialParticles->connect(fluid->importInitialStates());

	fluid->animationPipeline()->clear();
	{
		fluid->varReshuffleParticles()->setValue(false);

		fluid->varSmoothingLength()->setValue(2.5f);

		auto m_nbrQuery = std::make_shared<NeighborPointQuery<DataType3f>>();
		fluid->stateSmoothingLength()->connect(m_nbrQuery->inRadius());
		fluid->statePosition()->connect(m_nbrQuery->inPosition());
		fluid->animationPipeline()->pushModule(m_nbrQuery);

		auto m_stableISPHsolver = std::make_shared<SA_DFSPHsolver<DataType3f>>();
		fluid->stateSmoothingLength()->connect(m_stableISPHsolver->varSmoothingLength());
		//samplingDistance->outFloating()->connect(m_stableISPHsolver->varSamplingDistance());
		fluid->stateTimeStep()->connect(m_stableISPHsolver->inTimeStep());
		fluid->statePosition()->connect(m_stableISPHsolver->inPosition());
		fluid->stateVelocity()->connect(m_stableISPHsolver->inVelocity());
		//this->stateParticleAttribute()->connect(m_stableISPHsolver->inParticleAttribute());
		//this->stateBoundaryNorm()->connect(m_stableISPHsolver->inBoundaryNorm());
		m_nbrQuery->outNeighborIds()->connect(m_stableISPHsolver->inNeighborIds());
		fluid->animationPipeline()->pushModule(m_stableISPHsolver);

		//PBD
		//auto density = std::make_shared<IterativeDensitySolver<DataType3f>>();
		//smoothingLength->outFloating()->connect(density->inSmoothingLength());
		//samplingDistance->outFloating()->connect(density->inSamplingDistance());
		//fluid->stateTimeStep()->connect(density->inTimeStep());
		//fluid->statePosition()->connect(density->inPosition());
		//fluid->stateVelocity()->connect(density->inVelocity());
		////this->inNormal()->connect(density->inNormal());
		////this->inAttribute()->connect(density->inAttribute());
		//m_nbrQuery->outNeighborIds()->connect(density->inNeighborIds());
		//fluid->animationPipeline()->pushModule(density);

		//VSSPH
		//auto density = std::make_shared<VariationalApproximateProjection<DataType3f>>();
		//smoothingLength->outFloating()->connect(density->inSmoothingLength());
		//samplingDistance->outFloating()->connect(density->inSamplingDistance());
		//fluid->stateTimeStep()->connect(density->inTimeStep());
		//fluid->statePosition()->connect(density->inPosition());
		//fluid->stateVelocity()->connect(density->inVelocity());
		////this->inNormal()->connect(density->inNormal());
		////this->inAttribute()->connect(density->inAttribute());
		//m_nbrQuery->outNeighborIds()->connect(density->inNeighborIds());
		//fluid->animationPipeline()->pushModule(density);

		auto m_integrator = std::make_shared<RotatingSquarePatchModule<DataType3f>>();
		fluid->stateFrameNumber()->connect(m_integrator->inFrameNumber());
		fluid->stateTimeStep()->connect(m_integrator->inTimeStep());
		fluid->statePosition()->connect(m_integrator->inPosition());
		fluid->stateVelocity()->connect(m_integrator->inVelocity());
		fluid->stateParticleAttribute()->connect(m_integrator->inAttribute());
		fluid->animationPipeline()->pushModule(m_integrator);

		auto m_visModule = std::make_shared<ImplicitViscosity<DataType3f>>();
		m_visModule->varViscosity()->setValue(Real(0.1));
		fluid->stateTimeStep()->connect(m_visModule->inTimeStep());
		fluid->stateSmoothingLength()->connect(m_visModule->inSmoothingLength());
		fluid->stateTimeStep()->connect(m_visModule->inTimeStep());
		fluid->statePosition()->connect(m_visModule->inPosition());
		fluid->stateVelocity()->connect(m_visModule->inVelocity());
		m_nbrQuery->outNeighborIds()->connect(m_visModule->inNeighborIds());
		fluid->animationPipeline()->pushModule(m_visModule);
	
	}



	auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
	auto colorMapper = std::make_shared<ColorMapping<DataType3f >>();
	colorMapper->varMax()->setValue(5.0f);

	auto ptRender = std::make_shared<GLPointVisualModule>();
	ptRender->setColor(Color(1, 0, 0));
	ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);
	ptRender->varPointSize()->setValue(0.004);

	fluid->stateVelocity()->connect(calculateNorm->inVec());
	fluid->statePointSet()->connect(ptRender->inPointSet());

	calculateNorm->outNorm()->connect(colorMapper->inScalar());
	colorMapper->outColor()->connect(ptRender->inColor());

	fluid->graphicsPipeline()->pushModule(calculateNorm);
	fluid->graphicsPipeline()->pushModule(colorMapper);
	fluid->graphicsPipeline()->pushModule(ptRender);

	// A simple color bar widget for node
	auto colorBar = std::make_shared<ImColorbar>();
	colorBar->varMax()->setValue(5.0f);
	calculateNorm->outNorm()->connect(colorBar->inScalar());
	// add the widget to app
	fluid->graphicsPipeline()->pushModule(colorBar);

	return scene;
}

int main()
{
	GlfwApp window;
	window.setSceneGraph(createScene());
	window.initialize(1024, 768);
	window.mainLoop();

	return 0;
}




