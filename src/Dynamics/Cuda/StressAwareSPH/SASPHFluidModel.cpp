#include "SASPHFluidModel.h"

#include "SAISPHsolver.h"

#include "ParticleSystem/Module/ParticleIntegrator.h"
#include "StressAwareKernelSummation.h"

#include "ParticleSystem/Module/ImplicitViscosity.h"
#include "Collision/NeighborPointQuery.h"

namespace dyno
{
	IMPLEMENT_TCLASS(SASPHFluidModel, TDataType)

		template<typename TDataType>
	SASPHFluidModel<TDataType>::SASPHFluidModel()
		:GroupModule()
	{
		auto integrator = std::make_shared<ParticleIntegrator<TDataType>>();
		this->inTimeStep()->connect(integrator->inTimeStep());
		this->inPosition()->connect(integrator->inPosition());
		this->inVelocity()->connect(integrator->inVelocity());
		this->inAttribute()->connect(integrator->inAttribute());
		this->pushModule(integrator);

		auto nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		this->varSmoothingLength()->connect(nbrQuery->inRadius());
		this->inPosition()->connect(nbrQuery->inPosition());
		this->pushModule(nbrQuery);

		auto density = std::make_shared<SAISPHsolver<TDataType>>();
		this->varSmoothingLength()->connect(density->varSmoothingLength());
		this->varSamplingDistance()->connect(density->varSamplingDistance());
		this->inTimeStep()->connect(density->inTimeStep());
		this->inPosition()->connect(density->inPosition());
		this->inVelocity()->connect(density->inVelocity());
		this->inNormal()->connect(density->inBoundaryNorm());
		this->inAttribute()->connect(density->inParticleAttribute());
		nbrQuery->outNeighborIds()->connect(density->inNeighborIds());
		this->pushModule(density);

		auto viscosity = std::make_shared<ImplicitViscosity<TDataType>>();
		viscosity->varViscosity()->setValue(Real(0.3));
		this->inTimeStep()->connect(viscosity->inTimeStep());
		this->varSmoothingLength()->connect(viscosity->inSmoothingLength());
		this->inPosition()->connect(viscosity->inPosition());
		this->inVelocity()->connect(viscosity->inVelocity());
		nbrQuery->outNeighborIds()->connect(viscosity->inNeighborIds());
		this->pushModule(viscosity);
	}

	DEFINE_CLASS(SASPHFluidModel);
}