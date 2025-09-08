#pragma once
#include "StressAwareKernelParticleApproximation.h"

namespace dyno {
	/*
	* summation density for adpative kernel
	* kernel has parameter s 
	*/
	template<typename TDataType>
	class StressAwareKernelSummation : public virtual StressAwareKernelParticleApproximation<TDataType>
	{
		DECLARE_TCLASS(StressAwareKernelSummation, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		StressAwareKernelSummation();
		~StressAwareKernelSummation() override {};

		void compute() override;
	public:
		void compute(
			DArray<Real>& rho,
			DArray<Coord>& pos,
			DArrayList<int>& neighbors,
			DArray<Real>& nonpositivepressuredis,
			Real smoothingLength,
			Real mass);

		void compute(
			DArray<Real>& rho,
			DArray<Coord>& pos,
			DArray<Coord>& posQueried,
			DArrayList<int>& neighbors,
			DArray<Real>& nonpositivepressuredis,
			Real smoothingLength,
			Real mass);
	public:
		DEF_VAR(Real, RestDensity, 1000, "Rest Density");

		///Define inputs
		/**
		* @brief Particle positions
		*/
		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "Particle position");

		/**
		 * @brief Particle positions
		 */
		DEF_ARRAY_IN(Coord, Other, DeviceType::GPU, "Particle position");

		/**
		 * @brief Kernel Parameters
		 */
		DEF_ARRAY_IN(Real, KernelParameters, DeviceType::GPU, "Nearest Non-Positive Pressure Particle Distance");
		/**
		 * @brief Neighboring particles
		 *
		 */
		DEF_ARRAYLIST_IN(int, NeighborIds, DeviceType::GPU, "Neighboring particles' ids");

		///Define outputs
		/**
		 * @brief Particle densities
		 */
		DEF_ARRAY_OUT(Real, Density, DeviceType::GPU, "Return the particle density");

		Real getScalingFactor() {
			return this->mScalingFactor;
		}
		Real getParticleMass() {
			return m_particle_mass;
		}

	private:
		void calculateParticleMass();

		Real m_particle_mass;
		//Real m_factor;
	};

	IMPLEMENT_TCLASS(StressAwareKernelSummation, TDataType)
}
