#pragma once

#include "Module/ConstraintModule.h"
#include"Algorithm/Reduction.h"
#include "Algorithm/Functional.h"
#include "Algorithm/Arithmetic.h"
#include "Collision/Attribute.h"

#include "StressAwareKernel.h"
#include "StressAwareKernelSummation.h"

#include"Collision/NeighborPointQuery.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#define STRESSAWARE
//#define DFSPH_ORIGIN  //Do not apply semi-analytical free surface boundary 

namespace dyno
{
	class Attribute;
	template<typename TDataType> class StressAwareKernelSummation;
	template<typename TDataType> class NeighborPointQuery;

	template<typename TDataType>
	class SA_DFSPHsolver : public ConstraintModule
	{
		DECLARE_TCLASS(SA_DFSPHsolver, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		SA_DFSPHsolver();
		~SA_DFSPHsolver() override;

	public:
		DEF_VAR_IN(Real, TimeStep, "Time Step");

		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "Input particle position");

		DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "Input particle velocity");

		DEF_ARRAYLIST_IN(int, NeighborIds, DeviceType::GPU, "Neighboring particles' ids");

		DEF_ARRAY_OUT(Real, Density, DeviceType::GPU, "Final particle density");

		DEF_VAR(bool, DivergenceSolverDisabled, false, "Disable the Divergence solver in the DFSPH method");

		DEF_VAR(bool, DensitySolverDisabled, false, "Disable the Density solver in the DFSPH method");

		DEF_VAR(Real, RestDensity, 1000, "Reference density");

		DEF_VAR(Real, SamplingDistance, Real(0.005f), "");

		DEF_VAR(Real, SmoothingLength, Real(0.0100f), "");

		DEF_VAR(Real, DivergenceErrorThreshold, 0.001f, "Error Thershold for the Divergence solver");

		DEF_VAR(Real, DensityErrorThreshold, 0.001f, "Error Thershold for the Density solver");

		DEF_VAR(Real, MaxIterationNumber, 50, "Maximum number of iteration of each solver");

	public:
		void constrain() override;

	private:


		DArray<Real> mKappa_r;
		DArray<Real> mKappa_v;
		DArray<Real> mAlpha;
		DArray<Real> mAlpha_c;
		DArray<Real> mPressureEstimation;
		DArray<Real> mPara;
		DArray<Coord> mAccel;
		DArray<Real> mDensityAdv;
		DArray<Real> m_resv;
		DArray<Real> m_resr;

		//DArray<Real> m_Ap;
		//DArray<Real> m_source;
		//DArray<Real> m_p;
		//DArray<Real> m_r;

		//surface tension term
		DArray<Real> m_color;
		DArray<Real> m_gradC2;


	#ifdef STRESSAWARE
		StressAwareKernel<Real> m_kernel;
	#else
		SA_CubicKernel<Real> m_kernel;
	#endif

		unsigned int frag_number = 0;
		Real Alpha_min = 10.0f;

	private:
		std::shared_ptr<StressAwareKernelSummation<TDataType>> mSummation;


	};
}