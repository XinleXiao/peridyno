#include "StressAwareKernelParticleApproximation.h"

namespace dyno 
{
	IMPLEMENT_TCLASS(StressAwareKernelParticleApproximation, TDataType)

	template<typename TDataType>
	StressAwareKernelParticleApproximation<TDataType>::StressAwareKernelParticleApproximation()
		:ComputeModule()
	{
		auto callback = std::make_shared<FCallBackFunc>(
			std::bind(&StressAwareKernelParticleApproximation<TDataType>::calculateScalingFactor, this));

		this->varKernelType()->attach(callback);
		this->inSmoothingLength()->attach(callback);
		this->inSamplingDistance()->attach(callback);
	}

	template<typename TDataType>
	StressAwareKernelParticleApproximation<TDataType>::~StressAwareKernelParticleApproximation()
	{

	}

	template<typename TDataType>
	void StressAwareKernelParticleApproximation<TDataType>::calculateScalingFactor()
	{
		Real d = this->inSamplingDistance()->getValue();
		Real H = this->inSmoothingLength()->getValue();
		Real V = d * d * d;
		SA_AdaptiveKernel<Real>* kern;
		switch (this->varKernelType()->currentKey())
		{
		case KT_ACubic:
			kern = new SA_CubicKernel<Real>();
			break;
		case KT_ASpiky:
			kern = new SA_SpikyKernel<Real>();
			break;
		case KT_AStressAware:
			kern = new StressAwareKernel<Real>();
			break;
		default:
			break;
		}

		Real total_weight(0);
		int half_res = (int)(H / d + 1);
		for (int i = -half_res; i <= half_res; i++)
			for (int j = -half_res; j <= half_res; j++)
				for (int k = -half_res; k <= half_res; k++)
				{
					Real x = i * d;
					Real y = j * d;
					Real z = k * d;
					Real r = sqrt(x * x + y * y + z * z);
					total_weight += V * kern->Weight(r, H, r);
				}

		mScalingFactor = Real(1) / total_weight;// *0.9995f; //sight reduce the factor to maintain the init sturct of fluids

		delete kern;
	}

	DEFINE_CLASS(StressAwareKernelParticleApproximation);
}