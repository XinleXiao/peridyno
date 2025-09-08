//#pragma once
#include "ParticleSystem/ParticleFluid.h"
#include "ParticleSystem/Emitters/ParticleEmitter.h"
#include "Topology/PointSet.h"
#include "Collision/Attribute.h"

//#include "ParticleSystem/GhostParticles.h"

namespace dyno
{
	template<typename TDataType>
	class SASPHFluidSystem : public ParticleFluid<TDataType>
	{
		DECLARE_TCLASS(SASPHFluidSystem, TDataType)

	public:

		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		SASPHFluidSystem();
		~SASPHFluidSystem();

		DEF_ARRAY_STATE(Attribute, ParticleAttribute, DeviceType::GPU, "Particle Attribute");

		DEF_ARRAY_STATE(Coord, BoundaryNorm, DeviceType::GPU, "Boundary Norm");

		DEF_ARRAY_STATE(Coord, BoundaryVirtualPosition, DeviceType::GPU, "");


	protected:

		void resetStates();

		void preUpdateStates();

		void postUpdateStates();

	private:

	};

	IMPLEMENT_TCLASS(SASPHFluidSystem, TDataType)
}