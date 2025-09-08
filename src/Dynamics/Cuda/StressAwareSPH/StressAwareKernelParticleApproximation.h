#pragma once
#include"Module/ComputeModule.h"

#include"StressAwareKernel.h"

namespace dyno
{
#define cuIntegralAdh_SASPH(size, type, scale, Func,...){					\
		uint pDims = cudaGridSize((uint)size, BLOCK_SIZE);				\
		if (type == 0)											\
		{																\
			auto lambdaFunc = [=] __device__(Real r, Real h, Real r_t, Real s) -> Real {		\
				return SA_CubicKernel<Real>::_integral(r, h,  r, s );	\
			};																\
			Func << <pDims, BLOCK_SIZE >> > (__VA_ARGS__, lambdaFunc, scale);	\
		}																\
                                                                        \
		cuSynchronize();												\
	}

#define cuIntegral_SASPH(size, type, scale, Func,...){					\
		uint pDims = cudaGridSize((uint)size, BLOCK_SIZE);				\
		if (type == 0)											\
		{																\
			auto lambdaFunc = [=] __device__(Real r, Real h, Real r_t, Real s) -> Real {		\
				return SA_CubicKernel<Real>::_integral(r, h, r, s );	\
			};																\
			Func << <pDims, BLOCK_SIZE >> > (__VA_ARGS__, lambdaFunc, scale);	\
		}																\
                                                                        \
		cuSynchronize();												\
	}

#define cuZerothOrder_SASPH(size, type, scale, Func,...){					\
		uint pDims = cudaGridSize((uint)size, BLOCK_SIZE);				\
		if (type == 0)											\
		{																\
			auto lambdaFunc = [=] __device__(Real r, Real h, Real r_t, Real s) -> Real {		\
				return SA_CubicKernel<Real>::weight(r, h, r, s);	\
			};																\
			Func << <pDims, BLOCK_SIZE >> > (__VA_ARGS__, lambdaFunc, scale);	\
		}																\
                                                                        \
		cuSynchronize();												\
	}

#define cuFirstOrder_SASPH(size, type, scale, Func,...){					\
		uint pDims = cudaGridSize((uint)size, BLOCK_SIZE);				\
		if (type == 0)											\
		{																\
			auto lambdaFunc = [=] __device__(Real r, Real h, Real r_t, Real s) -> Real {		\
				return SA_CubicKernel<Real>::gradient(r, h,  r, s);	\
			};																\
			Func << <pDims, BLOCK_SIZE >> > (__VA_ARGS__, lambdaFunc, scale);	\
		}																\
                                                                        \
		cuSynchronize();												\
	}

	template<typename TDataType>
	class StressAwareKernelParticleApproximation : public ComputeModule
	{
		DECLARE_TCLASS(StressAwareKernelParticleApproximation, TDataType)
	public:
		typedef typename TDataType::Real Real;

		StressAwareKernelParticleApproximation();
		virtual ~StressAwareKernelParticleApproximation();

		DECLARE_ENUM(EKernelType,
			KT_ACubic = 0,
			KT_ASpiky = 1,
			KT_AStressAware = 2);

		void compute() override {};
	public:
		DEF_VAR_IN(Real, SmoothingLength, "Smoothing Length");
		DEF_VAR_IN(Real, SamplingDistance, "Particle sampling distance");

		DEF_ENUM(EKernelType, KernelType, EKernelType::KT_ACubic, "Kernel type");

	protected:
		Real mScalingFactor = Real(1);

	private:
		void calculateScalingFactor();
	};
}