#pragma once
#include "Module/ConstraintModule.h"
#include "Algorithm/Reduction.h"
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

namespace dyno {

	class Attribute;
	template<typename TDataType> class StressAwareKernelSummation;
	template<typename TDataType> class NeighborPointQuery;

	template<typename TDataType>
	class SAISPHsolver : public ConstraintModule
	{
		DECLARE_TCLASS(SAISPHsolver, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		//typedef typename TDataType::Matrix Matrix;

		SAISPHsolver();
		~SAISPHsolver();

		void constrain() override;

	public:
		DEF_VAR(Real, RestDensity, Real(1000.0f), "Reference density");

		DEF_VAR(Real, SamplingDistance, Real(0.005f), "");

		DEF_VAR(Real, SmoothingLength, Real(0.0125), "Smoothing length in most cases");

		DEF_VAR(int, PBDIterationNumber, 5, "Iteration number of the PBD solver for boundary Particle");

		DEF_VAR_IN(Real, TimeStep, "Time Step");

		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "Input Particle Position");

		//DEF_ARRAY_OUT(Coord, VPosition, DeviceType::GPU, "Input virtual Particle Position");

		DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "");

		DEF_ARRAY_IN(Attribute, ParticleAttribute, DeviceType::GPU, "");

		DEF_ARRAY_IN(Coord, BoundaryNorm, DeviceType::GPU, "Real-solid particle normal");

		DEF_ARRAYLIST_IN(int, NeighborIds, DeviceType::GPU, "Return real neighbor ids of real particles");

		DEF_ARRAY_OUT(Real, NPDistance, DeviceType::GPU, "Farest Non-Positive Pressure Particle Distance");

		DEF_ARRAY_OUT(Coord, BoundaryParticles, DeviceType::GPU, "Vitrual Boundary position");

		DEF_VAR(Real, ResiThreshold, 0.001f, "Convergence threshold for PPE");

		DEF_ARRAYLIST_OUT(int, VRNeighborIds, DeviceType::GPU, "Return virtual particles' real neighbor ids");

		DEF_ARRAYLIST_OUT(int, RVNeighborIds, DeviceType::GPU, "Return real particles' virtual neighbor ids");

		DEF_ARRAYLIST_OUT(int, VVNeighborIds, DeviceType::GPU, "Return virtual particles' virtual neighbor ids");

	protected:
		bool initializeImpl() override;
		bool ArraysResize();
		bool ArraysResize_Boundary(int num);


		StressAwareKernel<Real> kernel;
		SA_CubicKernel<Real> kernel_mass;
		Real m_maxAii;
		Real m_particleMass;
		Real m_airPressure;

		DArray<bool> m_BoundaryFlag;

		DArray<Real> m_source;
		DArray<Real> m_boundarySource;
		DArray<Real> m_Ax;
		DArray<Real> m_boundaryAx;
		DArray<Real> m_Aii;
		DArray<Real> m_boundaryAii;
		DArray<Real> m_r;
		DArray<Real> m_boundary_r;
		DArray<Real> m_p;
		DArray<Real> m_boundary_p;
		DArray<Real> m_pressure;
		DArray<Real> m_pressure_;
		DArray<Real> m_boundary_pressure;
		//DArray<Real> m_residual;
		//DArray<Real> m_boundary_residual;

		DArray<Coord> m_velocity;
		DArray<Coord> m_boundary_velocity;
		DArray<Coord> m_Gp;
		DArray<Coord> m_GpNearSolid;

		DArray<Real> m_Npdistance;
		DArray<Real> m_boundary_Npdistance;

		DArray<Real> m_rho;
		//DArray<Real> m_rho_boundary;

		DArray<Real> m_boundary_rho;
		DArray<Real> m_boundary_rho_projection;


		DArray<bool> m_boundaryAirFlag;
		DArray<bool> m_boundarySolidFlag;

		DArray<bool> m_boundary_deleteFlag;

		Reduction<Real>* m_reduce;
		Arithmetic<Real>* m_arithmetic;

		Reduction<Real>* m_boundary_reduce;
		Arithmetic<Real>* m_boundary_arithmetic;

		//
		// using morton codes for generate boundary points
		//
		DArray<Coord> m_anchorPoint;
		DArray<bool> m_anchorPointFlag;
		DArray<uint32_t> m_anchorPointCodes;
		DArray<uint32_t> m_nonRepeatedCount;
		DArray<uint32_t> m_candidateCodes;

		DArray<Coord> m_boundary_position;

		DArray<Coord> m_boundary_position_filtered;

		unsigned int boundaryVitrualNumber_old = 0;
		unsigned int points_num_old = 0;
		unsigned int frag_number = 0;
		unsigned int iter_sum = 0;

		Real max_Aii = 0.0f;
		Real max_Aii_f = 0.0f;
		Real max_Aii_b = 0.0f;

		//for PBD
		SA_SpikyKernel<Real> m_kernel_PBD;
		DArray<Real> m_lambda_PBD;
		DArray<Coord> m_deltaPos_PBD;

	private:
		std::shared_ptr<StressAwareKernelSummation<TDataType>> m_ak_summation; //Adaptive kernel summation

		std::shared_ptr <NeighborPointQuery<TDataType>> m_vv_nbrQuery;
		std::shared_ptr <NeighborPointQuery<TDataType>> m_rv_nbrQuery;
		std::shared_ptr <NeighborPointQuery<TDataType>> m_vr_nbrQuery;

	};
}