#include "StressAwareKernelSummation.h"

namespace dyno {
	template<typename TDataType>
	StressAwareKernelSummation<TDataType>::StressAwareKernelSummation()
		:StressAwareKernelParticleApproximation<TDataType>()
	{
		this->varRestDensity()->setValue(Real(1000));

		auto callback = std::make_shared<FCallBackFunc>(
			std::bind(&StressAwareKernelSummation<TDataType>::calculateParticleMass, this));

		this->varRestDensity()->attach(callback);
		this->inSamplingDistance()->attach(callback);

		this->inOther()->tagOptional(true);
		this->inKernelParameters()->tagOptional(true);
	}

	__global__ void KernelParameters_Init
	(
		DArray<Real> KernelParameters
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= KernelParameters.size()) return;

		KernelParameters[pId] = 0.0f;
	}

	template<typename TDataType>
	void StressAwareKernelSummation<TDataType>::compute()
	{
		int p_num = this->inPosition()->getDataPtr()->size();
		int n_num = this->inNeighborIds()->getDataPtr()->size();
		if (p_num != n_num) {
			Log::sendMessage(Log::Error, "The input array sizes of DensitySummation are not compatible!");
			return;
		}

		if (this->inKernelParameters()->isEmpty())
		{
			this->inKernelParameters()->allocate();
			this->inKernelParameters()->resize(p_num);
			cuExecute(p_num, KernelParameters_Init,
				this->inKernelParameters()->getData());
		}

		if (this->outDensity()->size() != p_num) {
			this->outDensity()->resize(p_num);
		}

		if (this->inOther()->isEmpty()) {
			compute(
				this->outDensity()->getData(),
				this->inPosition()->getData(),
				this->inNeighborIds()->getData(),
				this->inKernelParameters()->getData(),
				this->inSmoothingLength()->getData(),
				m_particle_mass);
		}
		else {
			compute(this->outDensity()->getData(),
				this->inPosition()->getData(),
				this->inOther()->getData(),
				this->inNeighborIds()->getData(),
				this->inKernelParameters()->getData(),
				this->inSmoothingLength()->getData(),
				m_particle_mass);
		}
	}

	template<typename Real, typename Coord, typename SA_AdaptiveKernel>
	__global__ void SA_ComputeDensity(
		DArray<Real> rhoArr,
		DArray<Coord> posArr,
		DArray<Real> ParameterArr,
		DArrayList<int> neighbors,
		Real smoothingLength,
		Real mass,
		SA_AdaptiveKernel weight,
		Real scale)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= posArr.size()) return;

		Real r;
		Real rho_i = Real(0);
		Coord pos_i = posArr[pId];
		Real para_i = ParameterArr[pId];
		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			r = (pos_i - posArr[j]).norm();
			rho_i += mass * weight(r, smoothingLength, para_i, scale);
		}
		rhoArr[pId] = rho_i;
	}

	template<typename Real, typename Coord, typename SA_AdaptiveKernel>
	__global__ void SA_ComputeDensity(
		DArray<Real> rhoArr,
		DArray<Coord> posArr,
		DArray<Real> ParameterArr,
		DArray<Coord> posQueried,
		DArrayList<int> neighbors,
		Real smoothingLength,
		Real mass,
		SA_AdaptiveKernel weight,
		Real scale)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= posArr.size()) return;

		Real r;
		Real rho_i = Real(0);
		Coord pos_i = posArr[pId];
		Real para_i = ParameterArr[pId];
		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			r = (pos_i - posQueried[j]).norm();
			rho_i += mass * weight(r, smoothingLength, para_i, scale);
		}
		rhoArr[pId] = rho_i;
	}

	template<typename TDataType>
	void StressAwareKernelSummation<TDataType>::compute(
		DArray<Real>& rho,
		DArray<Coord>& pos,
		DArrayList<int>& neighbors,
		DArray<Real>& ParameterArr,
		Real smoothingLength,
		Real mass)
	{
		cuZerothOrder_SASPH(rho.size(), this->varKernelType()->getDataPtr()->currentKey(), this->mScalingFactor,
			SA_ComputeDensity,
			rho,
			pos,
			ParameterArr,
			neighbors,
			smoothingLength,
			mass);
	}

	template<typename TDataType>
	void StressAwareKernelSummation<TDataType>::compute(
		DArray<Real>& rho,
		DArray<Coord>& pos,
		DArray<Coord>& posQueried,
		DArrayList<int>& neighbors,
		DArray<Real>& ParameterArr,
		Real smoothingLength,
		Real mass)
	{
		cuZerothOrder_SASPH(rho.size(), this->varKernelType()->getDataPtr()->currentKey(), this->mScalingFactor,
			SA_ComputeDensity,
			rho,
			pos,
			ParameterArr,
			posQueried,
			neighbors,
			smoothingLength,
			mass);

	}

	template<typename TDataType>
	void StressAwareKernelSummation<TDataType>::calculateParticleMass()
	{
		Real rho_0 = this->varRestDensity()->getData();
		Real d = this->inSamplingDistance()->getData();

		m_particle_mass = d * d * d * rho_0;
	}

	DEFINE_CLASS(StressAwareKernelSummation);
}