#include "initializeSASPHSystem.h"
#include"SASPHFluidSystem.h"
#include"SAISPHsolver.h"

#include "ParticleSystem/Emitters/CircularEmitter.h"
#include "ParticleSystem/Emitters/SquareEmitter.h"
#include "ParticleSystem/ParticleFluid.h"
#include "GLPointVisualModule.h"
#include "ColorMapping.h"
#include "Module/CalculateNorm.h"

#include "NodeFactory.h"

namespace dyno
{
	std::atomic<SASPHInitializer*> SASPHInitializer::gInstance;
	std::mutex SASPHInitializer::gMutex;

	PluginEntry* SASPHInitializer::instance()
	{
		SASPHInitializer* ins = gInstance.load(std::memory_order_acquire);
		if (!ins)
		{
			std::lock_guard<std::mutex> tLock(gMutex);
			ins = gInstance.load(std::memory_order_relaxed);
			if (!ins)
			{
				ins = new SASPHInitializer();
				ins->setName("Stress Aware Kernel SPH");
				ins->setVersion("1.0");
				ins->setDescription("An stress adaptive kernel method to remove tensile instability");

				gInstance.store(ins, std::memory_order_release);
			}
		}

		return ins;
	}

	SASPHInitializer::SASPHInitializer()
		:PluginEntry()
	{
	}

	void SASPHInitializer::initializeActions()
	{
		NodeFactory* factory = NodeFactory::instance();

		auto page = factory->addPage(
			"Particle System",
			"ToolBarIco/SASPH/SASPH_v1.png"
		);

		auto group = page->addGroup("AKSPH");

		group->addAction(
			"SASPH Fluid",
			"ToolBarIco/SASPH/SASPH_v1.png",
			[=]()->std::shared_ptr<Node> {
				auto fluid = std::make_shared<SASPHFluidSystem<DataType3f>>();

				auto calculateNorm = std::make_shared<CalculateNorm<DataType3f>>();
				fluid->stateVelocity()->connect(calculateNorm->inVec());
				fluid->graphicsPipeline()->pushModule(calculateNorm);

				auto colorMapper = std::make_shared<ColorMapping<DataType3f>>();
				colorMapper->varMax()->setValue(5.0f);
				calculateNorm->outNorm()->connect(colorMapper->inScalar());
				fluid->graphicsPipeline()->pushModule(colorMapper);

				auto ptRender = std::make_shared<GLPointVisualModule>();
				ptRender->setColor(Color(1, 0, 0));
				ptRender->setColorMapMode(GLPointVisualModule::PER_VERTEX_SHADER);

				fluid->statePointSet()->connect(ptRender->inPointSet());
				colorMapper->outColor()->connect(ptRender->inColor());
				fluid->graphicsPipeline()->pushModule(ptRender);

				return fluid;
			});
	}
}

PERIDYNO_API dyno::PluginEntry* SASPH::initDynoPlugin()
{
	if (dyno::SASPHInitializer::instance()->initialize())
		return dyno::SASPHInitializer::instance();

	return nullptr;
}

dyno::PluginEntry* SASPH::initStaticPlugin()
{
	if (dyno::SASPHInitializer::instance()->initialize())
		return dyno::SASPHInitializer::instance();

	return nullptr;
}