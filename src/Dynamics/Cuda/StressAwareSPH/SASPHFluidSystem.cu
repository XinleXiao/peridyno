#include "SASPHFluidSystem.h"
//DataType
#include"Auxiliary/DataSource.h"
//Collision
#include"Collision/NeighborPointQuery.h"

//ParticleSystem
#include "ParticleSystem/Module/ImplicitViscosity.h"
#include "ParticleSystem/Module/ParticleIntegrator.h"

#include "SAISPHsolver.h"
#include "SA_DFSPHsolver.h"

//#define STRESSAWAREISPH
#define STRESSAWAREDFSPH

namespace dyno
{
	__global__ void  SA_init_AttributeReset(
		DArray<Attribute> att
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= att.size()) return;

		att[pId].setFluid();
		att[pId].setDynamic();
	}

	template<typename TDataType>
	SASPHFluidSystem<TDataType>::SASPHFluidSystem()
		:ParticleFluid<TDataType>()
	{
		this->varReshuffleParticles()->setValue(false);

		this->animationPipeline()->clear();

		auto smoothingLength = std::make_shared<FloatingNumber<TDataType>>();
		this->animationPipeline()->pushModule(smoothingLength);
		smoothingLength->varValue()->setValue(Real(0.005f*2.5f));

		auto samplingDistance = std::make_shared<FloatingNumber<TDataType>>();
		this->animationPipeline()->pushModule(smoothingLength);
		samplingDistance->varValue()->setValue(Real(0.005f));

		auto restDensity = std::make_shared<FloatingNumber<TDataType>>();
		this->animationPipeline()->pushModule(restDensity);
		restDensity->varValue()->setValue(Real(1000.0f));

		auto m_nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		smoothingLength->outFloating()->connect(m_nbrQuery->inRadius());
		this->statePosition()->connect(m_nbrQuery->inPosition());
		this->animationPipeline()->pushModule(m_nbrQuery);

#ifdef STRESSAWAREISPH
		auto m_SAISPHsolver = std::make_shared<SAISPHsolver<TDataType>>();
		smoothingLength->outFloating()->connect(m_SAISPHsolver->varSmoothingLength());
		samplingDistance->outFloating()->connect(m_SAISPHsolver->varSamplingDistance());
		this->stateTimeStep()->connect(m_SAISPHsolver->inTimeStep());
		this->statePosition()->connect(m_SAISPHsolver->inPosition());
		this->stateVelocity()->connect(m_SAISPHsolver->inVelocity());
		//this->stateParticleAttribute()->connect(m_SAISPHsolver->inParticleAttribute());
		//this->stateBoundaryNorm()->connect(m_SAISPHsolver->inBoundaryNorm());
		m_nbrQuery->outNeighborIds()->connect(m_SAISPHsolver->inNeighborIds());
		this->animationPipeline()->pushModule(m_SAISPHsolver);
#endif // STRESSAWAREISPH

#ifdef STRESSAWAREDFSPH
		auto m_SADFSPH = std::make_shared<SA_DFSPHsolver<TDataType>>();
		smoothingLength->outFloating()->connect(m_SADFSPH->varSmoothingLength());
		samplingDistance->outFloating()->connect(m_SADFSPH->varSamplingDistance());
		this->stateTimeStep()->connect(m_SADFSPH->inTimeStep());
		this->statePosition()->connect(m_SADFSPH->inPosition());
		this->stateVelocity()->connect(m_SADFSPH->inVelocity());
		m_nbrQuery->outNeighborIds()->connect(m_SADFSPH->inNeighborIds());
		this->animationPipeline()->pushModule(m_SADFSPH);
#endif // STRESSAWAREDFSPH

		auto m_integrator = std::make_shared<ParticleIntegrator<TDataType>>();
		this->stateTimeStep()->connect(m_integrator->inTimeStep());
		this->statePosition()->connect(m_integrator->inPosition());
		this->stateVelocity()->connect(m_integrator->inVelocity());
		this->stateParticleAttribute()->connect(m_integrator->inAttribute());
		this->animationPipeline()->pushModule(m_integrator);

		auto m_visModule = std::make_shared<ImplicitViscosity<TDataType>>();
		m_visModule->varViscosity()->setValue(Real(0.3f));
		this->stateTimeStep()->connect(m_visModule->inTimeStep());
		smoothingLength->outFloating()->connect(m_visModule->inSmoothingLength());
		this->stateTimeStep()->connect(m_visModule->inTimeStep());
		this->statePosition()->connect(m_visModule->inPosition());
		this->stateVelocity()->connect(m_visModule->inVelocity());
		m_nbrQuery->outNeighborIds()->connect(m_visModule->inNeighborIds());
		this->animationPipeline()->pushModule(m_visModule);

		this->setDt(Real(0.001f));
	}

	template<typename TDataType>
	SASPHFluidSystem<TDataType>::~SASPHFluidSystem()
	{

	}

	template<typename TDataType>
	void SASPHFluidSystem<TDataType>::resetStates()
	{
		this->ParticleFluid<TDataType>::resetStates();

		auto ptSet = this->statePointSet()->getDataPtr();
		if (ptSet != nullptr)
		{
			auto pts = ptSet->getPoints();
			this->stateBoundaryNorm()->resize(pts.size());
			this->stateParticleAttribute()->resize(pts.size());

			cuExecute(pts.size(), SA_init_AttributeReset,
				this->stateParticleAttribute()->getData());

			this->stateBoundaryNorm()->getDataPtr()->reset();
		}
	}

	template<typename TDataType>
	void SASPHFluidSystem<TDataType>::preUpdateStates()
	{
		this->varReshuffleParticles()->setValue(false);
		this->ParticleFluid<TDataType>::preUpdateStates();

		this->stateBoundaryNorm()->resize(this->statePosition()->size());
		this->stateBoundaryNorm()->reset();
		this->stateParticleAttribute()->resize(this->statePosition()->size());


		cuExecute(this->statePosition()->size(), SA_init_AttributeReset,
			this->stateParticleAttribute()->getData());
	}

	template<typename TDataType>
	void SASPHFluidSystem<TDataType>::postUpdateStates()
	{
		this->ParticleSystem<TDataType>::postUpdateStates();

		return;
	}

	DEFINE_CLASS(SASPHFluidSystem);
}