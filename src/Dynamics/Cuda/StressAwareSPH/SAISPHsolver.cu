#include <cuda_runtime.h>
#include "SAISPHsolver.h"
#include "Node.h"
#include "Field.h"
#include "StressAwareKernelSummation.h"
#include "Collision/Attribute.h"
#include "StressAwareKernel.h"
#include "Algorithm/Function2Pt.h"
#include "Algorithm/CudaRand.h"

#include"Collision/NeighborPointQuery.h"

#include<thrust/sort.h>

namespace dyno {

	IMPLEMENT_TCLASS(SAISPHsolver, TDataType)

		template<typename TDataType>
	SAISPHsolver<TDataType>::SAISPHsolver()
		:ConstraintModule()
		,m_airPressure(Real(0))
		,m_reduce(NULL)
		,m_arithmetic(NULL)
		,m_boundary_reduce(NULL)
		,m_boundary_arithmetic(NULL)
	{
		this->inParticleAttribute()->tagOptional(true);
		this->inBoundaryNorm()->tagOptional(true);

		this->varSamplingDistance()->setValue(Real(0.005));
		this->varRestDensity()->setValue(Real(1000));
		this->varSmoothingLength()->setValue(Real(0.010));

		this->outBoundaryParticles()->allocate();

		m_ak_summation = std::make_shared<StressAwareKernelSummation<TDataType>>();
		this->varRestDensity()->connect(m_ak_summation->varRestDensity());
		this->varSmoothingLength()->connect(m_ak_summation->inSmoothingLength());
		this->varSamplingDistance()->connect(m_ak_summation->inSamplingDistance());
		this->inPosition()->connect(m_ak_summation->inPosition());
		this->inNeighborIds()->connect(m_ak_summation->inNeighborIds());

		this->outVVNeighborIds()->allocate();
		m_vv_nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		this->varSmoothingLength()->connect(m_vv_nbrQuery->inRadius());
		this->outBoundaryParticles()->connect(m_vv_nbrQuery->inPosition());
		m_vv_nbrQuery->outNeighborIds()->connect(this->outVVNeighborIds());


		this->outRVNeighborIds()->allocate();
		m_rv_nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		this->varSmoothingLength()->connect(m_rv_nbrQuery->inRadius());
		this->outBoundaryParticles()->connect(m_rv_nbrQuery->inPosition());
		this->inPosition()->connect(m_rv_nbrQuery->inOther());
		m_rv_nbrQuery->outNeighborIds()->connect(this->outRVNeighborIds());

		this->outVRNeighborIds()->allocate();
		m_vr_nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
		this->varSmoothingLength()->connect(m_vr_nbrQuery->inRadius());
		this->inPosition()->connect(m_vr_nbrQuery->inPosition());
		this->outBoundaryParticles()->connect(m_vr_nbrQuery->inOther());
		m_vr_nbrQuery->outNeighborIds()->connect(this->outVRNeighborIds());
	}

	template<typename TDataType>
	SAISPHsolver<TDataType>::~SAISPHsolver()
	{
		m_r.clear();
		m_Aii.clear();
		m_boundary_r.clear();
		m_boundaryAii.clear();
		m_pressure.clear();
		m_boundary_pressure.clear();

		if (m_reduce)
		{
			delete m_reduce;
		}
		if (m_arithmetic)
		{
			delete m_arithmetic;
		}

		if (m_boundary_reduce)
		{
			delete m_boundary_reduce;
		}
		if (m_boundary_arithmetic)
		{
			delete m_boundary_arithmetic;
		}
	}

	template<typename TDataType>
	bool SAISPHsolver<TDataType>::ArraysResize()
	{
		int num = this->inPosition()->size();

		if (m_Ax.size() != num)
			m_Ax.resize(num);
		if (m_Aii.size() != num)
			m_Aii.resize(num);
		if (m_r.size() != num)
			m_r.resize(num);
		if (m_p.size() != num)
			m_p.resize(num);
		if (m_BoundaryFlag.size() != num)
			m_BoundaryFlag.resize(num);
		if (m_velocity.size() != num)
			m_velocity.resize(num);
		if (m_source.size() != num)
			m_source.resize(num);
		if (m_pressure.size() != num)
			m_pressure.resize(num);
		if (m_pressure_.size() != num)
			m_pressure_.resize(num);
		if (m_Gp.size() != num)
			m_Gp.resize(num);
		if (m_GpNearSolid.size() != num)
			m_GpNearSolid.resize(num);
		if (m_rho.size() != num)
			m_rho.resize(num);

		if (m_reduce)
		{
			delete m_reduce;
			m_reduce = Reduction<Real>::Create(num);
		}
		else
		{
			m_reduce = Reduction<Real>::Create(num);
		}

		if (m_arithmetic)
		{
			delete m_arithmetic;
			m_arithmetic = Arithmetic<Real>::Create(num);
		}
		else
		{
			m_arithmetic = Arithmetic<Real>::Create(num);
		}

		return true;
	}

	template<typename TDataType>
	bool SAISPHsolver<TDataType>::ArraysResize_Boundary(int num_b)
	{
		int num = num_b;

		if (m_boundary_rho.size() != num)
			m_boundary_rho.resize(num);
		if (m_boundary_position.size() != num)
			m_boundary_position.resize(num);
		if (m_boundarySource.size() != num)
			m_boundarySource.resize(num);
		if (m_boundaryAx.size() != num)
			m_boundaryAx.resize(num);
		if (m_boundaryAii.size() != num)
			m_boundaryAii.resize(num);
		if (m_boundary_r.size() != num)
			m_boundary_r.resize(num);
		if (m_boundary_p.size() != num)
			m_boundary_p.resize(num);
		if (m_boundary_pressure.size() != num)
			m_boundary_pressure.resize(num);
		if (m_boundaryAirFlag.size() != num)
			m_boundaryAirFlag.resize(num);
		if (m_boundarySolidFlag.size() != num)
			m_boundarySolidFlag.resize(num);
		if (m_boundary_deleteFlag.size() != num)
			m_boundary_deleteFlag.resize(num);
		if (m_boundary_velocity.size() != num)
			m_boundary_velocity.resize(num);
		if (m_deltaPos_PBD.size() != num)
			m_deltaPos_PBD.resize(num);
		if (m_lambda_PBD.size() != num)
			m_lambda_PBD.resize(num);
		if (m_boundary_rho_projection.size() != num)
			m_boundary_rho_projection.resize(num);

		if (m_boundary_reduce)
		{
			delete m_boundary_reduce;
			m_boundary_reduce = Reduction<Real>::Create(num);
		}
		else
		{
			m_boundary_reduce = Reduction<Real>::Create(num);
		}

		if (m_boundary_arithmetic)
		{
			delete m_boundary_arithmetic;
			m_boundary_arithmetic = Arithmetic<Real>::Create(num);
		}
		else
		{
			m_boundary_arithmetic = Arithmetic<Real>::Create(num);
		}

		return true;
	}

	__global__ void SAISPH_AttributeInit
	(
		DArray<Attribute> atts
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= atts.size()) return;

		atts[pId].setFluid();
		atts[pId].setDynamic();
	}

	template<typename Real>
	__global__ void SAISPH_BoundaryParticleDetection
	(
		DArray<bool> boundaryFlag,
		DArray<Attribute> attribute,
		DArrayList<int> neighborIds,
		DArray<Real> rho,
		Real threshold
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= boundaryFlag.size()) return;
		if (!attribute[pId].isFluid()) return;

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		boundaryFlag[pId] = true;

		if ((rho[pId] > threshold) && (nbSize > 0))
		{
			boundaryFlag[pId] = false;

		}
		else
		{
			boundaryFlag[pId] = true;
		}

	}

	__constant__ int diff_v_33[33][3] =
	{
		0, 0, 0,
		0, 0, 1,
		0, 1, 0,
		1, 0, 0,
		0, 0, -1,
		0, -1, 0,
		-1, 0, 0,
		0, 1, 1,
		0, 1, -1,
		0, -1, 1,
		0, -1, -1,
		1, 0, 1,
		1, 0, -1,
		-1, 0, 1,
		-1, 0, -1,
		1, 1, 0,
		1, -1, 0,
		-1, 1, 0,
		-1, -1, 0,
		1, 1, 1,
		1, 1, -1,
		1, -1, 1,
		-1, 1, 1,
		1, -1, -1,
		-1, 1, -1,
		-1, -1, 1,
		-1, -1, -1,
		2, 0, 0,
		-2, 0, 0,
		0, 2, 0,
		0, -2, 0,
		0, 0, 2,
		0, 0, -2

	};

	template<typename Real, typename Coord>
	__global__ void SAISPH_AnchorNeighbor_33
	(
		DArray<Coord> anchorPoint,
		DArray<bool> anchorPointFlag,
		DArray<Coord> pos,
		DArray<bool> boundaryFlag,
		Coord origin,
		Real dh
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size())return;
		if (boundaryFlag[pId] == false)
		{
			for (int i = 0; i < 33; i++)
				anchorPointFlag[pId * 33 + i] = false;
			return;
		}

		Coord pos_ref = pos[pId] - origin + dh / 2;
		Coord a(0.0f);

		a[0] = (Real)((int)(floor(pos_ref[0] / dh))) * dh;
		a[1] = (Real)((int)(floor(pos_ref[1] / dh))) * dh;
		a[2] = (Real)((int)(floor(pos_ref[2] / dh))) * dh;

		for (int i = 0; i < 33; i++)
		{
			anchorPoint[pId * 33 + i] = Coord(a[0] + dh * Real(diff_v_33[i][0]),
				a[1] + dh * Real(diff_v_33[i][1]),
				a[2] + dh * Real(diff_v_33[i][2])
			);
			anchorPointFlag[pId * 33 + i] = true;
		}
	}

	template<typename Coord>
	__global__ void SAISPH_PositionMortonCode
	(
		DArray<uint32_t> mortonCode,
		DArray<bool> anchorPointFlag,
		DArray<Coord> gridPoint,
		Real dh
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId > gridPoint.size())return;
		if (anchorPointFlag[pId] == false)
		{
			mortonCode[pId] = 0;
			return;
		}

		uint32_t  x = (uint32_t)(gridPoint[pId][0] / dh + 100 * EPSILON);
		uint32_t  y = (uint32_t)(gridPoint[pId][1] / dh + 100 * EPSILON);
		uint32_t  z = (uint32_t)(gridPoint[pId][2] / dh + 100 * EPSILON);
		

		uint32_t  xi = x;
		uint32_t  yi = y;
		uint32_t  zi = z;
		
		uint32_t  key = 0;


		xi = (xi | (xi << 16)) & 0x030000FF;
		xi = (xi | (xi << 8)) & 0x0300F00F;
		xi = (xi | (xi << 4)) & 0x030C30C3;
		xi = (xi | (xi << 2)) & 0x09249249;

		yi = (yi | (yi << 16)) & 0x030000FF;
		yi = (yi | (yi << 8)) & 0x0300F00F;
		yi = (yi | (yi << 4)) & 0x030C30C3;
		yi = (yi | (yi << 2)) & 0x09249249;

		zi = (zi | (zi << 16)) & 0x030000FF;
		zi = (zi | (zi << 8)) & 0x0300F00F;
		zi = (zi | (zi << 4)) & 0x030C30C3;
		zi = (zi | (zi << 2)) & 0x09249249;

		key = xi | (yi << 1) | (zi << 2);

		mortonCode[pId] = key;
	}

	__global__ void SAISPH_RepeatedElementSearch
	(
		DArray<uint32_t> morton,
		DArray<uint32_t> counter
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= morton.size()) return;
		if (pId == 0)
		{
			counter[pId] = 0;
			return;
		}

		//if (pId == 0 || morton[pId] != morton[pId - 1])
		if (morton[pId] != morton[pId - 1])
		{
			counter[pId] = 1;
		}
		else
		{
			counter[pId] = 0;
		}
	}

	__global__ void SAISPH_nonRepeatedElementsCal

	(
		DArray<uint32_t>  non_repeated_elements,
		DArray<uint32_t>  post_elements,
		DArray<uint32_t> counter
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= post_elements.size()) return;
		if (pId == 0)return;

		if (post_elements[pId] != post_elements[pId - 1])
		{
			//printf("counter:%d,%d,%d", counter[pId - 1], counter[pId], counter[pId + 1]);
			non_repeated_elements[counter[pId]] = post_elements[pId];
		}
	}

	template<typename Coord>
	__global__ void SAISPH_MortonCodeToPosition
	(
		DArray<Coord> gridPoint,
		DArray<uint32_t>  mortonCode,
		Real dh
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= mortonCode.size()) return;


		uint32_t key = mortonCode[pId];

		uint32_t  xv = key & 0x09249249;
		uint32_t  yv = (key >> 1) & 0x09249249;
		uint32_t  zv = (key >> 2) & 0x09249249;

		xv = ((xv >> 2) | xv) & 0x030C30C3;
		xv = ((xv >> 4) | xv) & 0x0300F00F;
		xv = ((xv >> 8) | xv) & 0x030000FF;
		xv = ((xv >> 16) | xv) & 0x000003FF;

		yv = ((yv >> 2) | yv) & 0x030C30C3;
		yv = ((yv >> 4) | yv) & 0x0300F00F;
		yv = ((yv >> 8) | yv) & 0x030000FF;
		yv = ((yv >> 16) | yv) & 0x000003FF;

		zv = ((zv >> 2) | zv) & 0x030C30C3;
		zv = ((zv >> 4) | zv) & 0x0300F00F;
		zv = ((zv >> 8) | zv) & 0x030000FF;
		zv = ((zv >> 16) | zv) & 0x000003FF;

		Real x = (float)(xv)*dh;
		Real y = (float)(yv)*dh;
		Real z = (float)(zv)*dh;

		gridPoint[pId] = Coord(x, y, z);
	}

	template< typename Coord>
	__global__ void SAISPH_CopyToBoundaryPos(
		DArray<Coord> vpos,
		Coord origin,
		DArray<Coord>  elements

	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= elements.size()) return;
		vpos[pId] = elements[pId] + origin;
		elements[pId] = elements[pId] + origin;
	}

	template<typename Real, typename Coord>
	__global__ void SAISPH_compute_b_rho(
		DArray<Coord> b_pos,
		DArrayList<int> vr_neighbors,
		DArray<Coord> position,
		DArray<bool> boundaryFlag,
		DArray<Real> Rho_b,
		DArray<Real> Rho_b_p,
		DArray<Attribute> attribute,
		DArrayList<int> vv_neighbors,
		SA_CubicKernel<Real> kernel,
		Real mass,
		Real rho_factor,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= b_pos.size())return;

		//Real fluid particle
		List<int>& list_i = vr_neighbors[pId];
		int nbSize = list_i.size();
		Real temp_v(0.0f);
		Real temp_v_p(0.0f);
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real rij = (b_pos[pId] - position[j]).norm();
			Real wij = kernel.Weight(rij, h, rij);
			if (boundaryFlag[j] == true)
			{
				temp_v_p += mass * wij * rho_factor;
			}
			else
			{
				temp_v += mass * wij * rho_factor;
				temp_v_p += mass * wij * rho_factor;
			}

		}
		//virtual boundary particle
		List<int>& list_i_b = vv_neighbors[pId];
		int nbSize_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_b; ne_b++)
		{
			int j_b = list_i_b[ne_b];
			Real rij_b = (b_pos[pId] - b_pos[j_b]).norm();
			Real wij = kernel.Weight(rij_b, h, rij_b);
			temp_v += mass * wij * rho_factor;

		}
		Rho_b[pId] = temp_v;
		Rho_b_p[pId] = temp_v_p;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_DeleteNearParticles
	(
		DArray<Coord> boudaryParticlePos,
		DArray<bool> boundarydeleteflag,
		DArray<Real> Rho_vv,
		Real rho_0,
		Real MaxDensity,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= boudaryParticlePos.size())return;
		boundarydeleteflag[pId] = false;
		//delete particles with large a
		if (Rho_vv[pId] > 1.05f * MaxDensity && Rho_vv[pId] > rho_0 * 1.0f)
		{
			boundarydeleteflag[pId] = true;
		}
		/*if (Rho_vv[pId] < rho_0 * 0.2f)
		{
			boundarydeleteflag[pId] = true;
		}*/
	}

	template<typename Real, typename Coord>
	__global__ void SAISPH_compute_b_rho_filtered(
		DArray<Coord> b_pos,
		DArray<bool> boundarydeleteFlag,
		DArrayList<int> vr_neighbors,
		DArray<Coord> position,
		DArray<bool> boundaryFlag,
		DArray<Real> Rho_b,
		DArray<Real> Rho_b_p,
		DArray<Attribute> attribute,
		DArrayList<int> vv_neighbors,
		SA_CubicKernel<Real> kernel,
		Real mass,
		Real rho_factor,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= b_pos.size())return;
		if (boundarydeleteFlag[pId] == true)
		{
			Rho_b[pId] = 1000;
			Rho_b_p[pId] = 1000;
			return;
		}

		//Real fluid particle
		List<int>& list_i = vr_neighbors[pId];
		int nbSize = list_i.size();
		Real temp_v(0.0f);
		Real temp_v_p(0.0f);
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			if (boundaryFlag[j] == true)
			{
				Real rij = (b_pos[pId] - position[j]).norm();
				Real wij = kernel.Weight(rij, h, rij);
				temp_v_p += mass * wij * rho_factor;
			}
			else
			{
				Real rij = (b_pos[pId] - position[j]).norm();
				Real wij = kernel.Weight(rij, h, rij);
				temp_v += mass * wij * rho_factor;
				temp_v_p += mass * wij * rho_factor;
			}

		}
		//virtual boundary particle
		List<int>& list_i_b = vv_neighbors[pId];
		int nbSize_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_b; ne_b++)
		{
			
			int j_b = list_i_b[ne_b];
			if (boundarydeleteFlag[j_b] == false) {
				Real rij_b = (b_pos[pId] - b_pos[j_b]).norm();
				Real wij = kernel.Weight(rij_b, h, rij_b);
				temp_v += mass * wij * rho_factor;
			}

		}
	
		Rho_b[pId] = temp_v;
		Rho_b_p[pId] = temp_v_p;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_InitBoundaryParticles
	(
		DArray<Coord> boudaryParticlePos,
		DArray<bool> boundarydeleteFlag,
		DArray<Real> Rho_b_p,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArrayList<int> vr_neighbors,
		DArray<Coord> particlePos,
		DArray<Real> particleRho,
		DArray<Attribute> attribute,
		DArray<bool> particleBoundaryFlag,
		SA_CubicKernel<Real> kernel,
		Real mass,
		Real rho_factor,
		Real h,
		Real threshold_air,
		Real threshold_solid
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= boudaryParticlePos.size())return;
		if (boundarydeleteFlag[pId] == true)
		{
			boundaryAirFlag[pId] = true;
			return;
		}

		List<int>& list_i = vr_neighbors[pId];


		int nbSize = list_i.size();
		//Real rho_boundary = 0.0f;
		Real rho_forairflag = 0.0f;
		Real c = 0.0f;
		Real w = 0.0f;
		//Real wsum = 0.0f;

		//Real temp = 0.0f;

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (boudaryParticlePos[pId] - particlePos[j]).norm();
			Real wij = kernel.Weight(r_ij, h, r_ij);

			rho_forairflag += mass * wij * rho_factor;

			if (!attribute[j].isFluid())
			{
				c += wij * mass / particleRho[j];
			}
			w += wij * mass / particleRho[j];
		}

		if (rho_forairflag < threshold_air && nbSize > 0)
		{
			boundaryAirFlag[pId] = true;
		}
		else
		{
			boundaryAirFlag[pId] = false;
		}

		if (w < EPSILON)w = EPSILON;

		if (c / w > threshold_solid)
		{
			boundarySolidFlag[pId] = true;
		}
		else
		{
			boundarySolidFlag[pId] = false;
		}

		if (Rho_b_p[pId] < threshold_air)
		{
			boundaryAirFlag[pId] = true;
		}
	}

	

	template<typename Real, typename Coord>
	__global__ void SAISPH_SphBoundaryVelocity(
		DArray<Coord> boundaryVelocity,
		DArray<Coord> velocity,
		DArray<Coord> boundaryPos,
		DArray<Coord> pos,
		DArrayList<int> vr_neighbors,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArray<Real> rho,
		SA_CubicKernel<Real> kernel,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= boundaryVelocity.size()) return;
		if (boundarySolidFlag[pId] == true)
		{
			return;
		}
		if (boundaryAirFlag[pId] == true)
		{
			boundaryVelocity[pId] = Coord(0.0f);
			return;
		}
		List<int>& list_i = vr_neighbors[pId];
		int nbSize = list_i.size();
		Real total_w(0);
		Coord total_vw(0);
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (boundaryPos[pId] - pos[j]).norm();
			Real wij = kernel.Weight(r_ij, h, r_ij);
			total_w += wij;
			total_vw += velocity[j] * wij;
		}

		if (total_w > EPSILON)
		{
			if (nbSize == 0)printf("1-ERROR: pId  %d \r\n", pId);
			boundaryVelocity[pId] = total_vw / total_w;
		}
		else
		{
			if (nbSize != 0) printf("0-ERROR: pId  %d \r\n", pId);
			boundaryVelocity[pId] = Coord(0.0f);
		}
	}

	template <typename Coord>
	__global__ void SAISPH_SolidVelocityReset
	(
		DArray<Coord> velocity,
		DArray<Attribute> attribute
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= velocity.size()) return;
		if (!attribute[pId].isFluid())
		{
			velocity[pId] = Coord(0.0f);
		}

	}

	template<typename Real>
	__global__ void SAISPH_Init_pressure_
	(
		DArray<Real> pressure,
		Real pressure_max,
		DArray<Real> pressure_
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId > pressure.size()) return;

		pressure_[pId] = pressure[pId] / pressure_max * 1.0f;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_SourceTerm
	(
		DArray<Real> source,
		DArray<Coord> velocity,
		DArray<Coord> position,
		DArray<Real> pressure_,
		DArray<Real> rho,
		DArray<bool> boundaryFlag,
		DArrayList<int> neighborIds,
		DArray<Attribute> attribute,
		DArray<Coord> boundaryPos,
		DArray<Coord> boundaryVelocity,
		DArray<Real> boundaryRho,
		DArrayList<int> RVneighborIds,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArray<Coord> boundaryNorm,
		StressAwareKernel<Real> kernel,
		Real restDensity,
		Real mass,
		Real dt,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= source.size()) return;

		if (boundaryFlag[pId] == true)
		{
			source[pId] = 0.0f;
			return;
		}
		else
		{
			List<int> list_i = neighborIds[pId];
			int nbSize = list_i.size();
			Real value(0.0f);
			for (int ne = 0; ne < nbSize; ne++)
			{
				int j = list_i[ne];
				//if (boundaryFlag[j] == false)
				{
					Real r_ij = (position[pId] - position[j]).norm();

					Coord Gwij(0.0f);
					if (r_ij > EPSILON && attribute[j].isFluid())
					{
						Real s = r_ij * expf(-(pressure_[pId] + pressure_[j]) / 2.0f);
						Gwij = kernel.Gradient(r_ij, h, s) * (position[pId] - position[j]) / r_ij;
					}
					value += (velocity[j] - velocity[pId]).dot(Gwij) * mass / rho[pId];
					//value += (velocity[j] / rho[j] / rho[j] + velocity[pId] / rho[pId] / rho[pId]).dot(Gwij) * mass / rho[pId];
				}
			}
			source[pId] = -value * restDensity / dt;
		}
	}


	template<typename Coord, typename Real>
	__global__ void SAISPH_BoundarySourceTerm(
		DArray<Real> boundarySource,
		DArray<Coord> velocity,
		DArray<Coord> position,
		DArray<Real> pressure_,
		DArray<Real> rho,
		DArray<bool> boundaryFlag,
		DArrayList<int> VRneighborIds,
		DArray<Attribute> attribute,
		DArray<Coord> boundaryPos,
		DArray<Coord> boundaryVelocity,
		DArray<Real> boundaryRho,
		DArrayList<int> VVneighborIds,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArray<Coord> boundaryNorm,
		StressAwareKernel<Real> kernel,
		Real restDensity,
		Real mass,
		Real dt,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= boundarySource.size())return;

		if (boundarySolidFlag[pId] == true)
		{
			boundarySource[pId] = 0.0f;
			return;
		}

		if (boundaryAirFlag[pId] == true)
		{
			boundarySource[pId] = 0.0f;
			return;
		}
		else
		{
			List<int> list_i_vr = VRneighborIds[pId];
			int nbSize_vr = list_i_vr.size();

			Real value(0.0f);
			for (int ne = 0; ne < nbSize_vr; ne++)
			{
				int j = list_i_vr[ne];
				//if (boundaryFlag[j] == false) 
				{
					Real r_ij = (boundaryPos[pId] - position[j]).norm();
					Coord Gwij(0);
					if (r_ij > EPSILON && attribute[j].isFluid())
					{
						Real s = r_ij * expf(-pressure_[j] / 2.0f);


						Gwij = kernel.Gradient(r_ij, h, s) * (boundaryPos[pId] - position[j]) / r_ij;
					}
					value += (velocity[j] - boundaryVelocity[pId]).dot(Gwij) * mass / boundaryRho[pId];// rho[j];
					//value += (velocity[j] / rho[j] / rho[j] + boundaryVelocity[pId] / boundaryRho[pId] / boundaryRho[pId]).dot(Gwij) * mass / boundaryRho[pId];
				}
			}

			boundarySource[pId] = -value * restDensity / dt;
		}
	}

	template<typename Real>
	__global__ void SAISPH_DensityCompensate
	(
		DArray<Real> source,
		DArray<Real> rho,
		DArray<Attribute> attribute,
		DArray<bool> boundaryFlag,
		Real rho_0,
		Real dt,
		Real h
	)
	{

		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= source.size()) return;
		if (!attribute[pId].isFluid()) return;
		if (boundaryFlag[pId] == true) return;
		if (rho[pId] > rho_0)
		{
			source[pId] += 1000000.0f * (rho[pId] - rho_0) / rho_0;
		}
	}

	template<typename Real>
	__global__ void SAISPH_BoundaryDensityCompensate
	(
		DArray<Real> boundarySource,
		DArray<Real> boundaryRho,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		Real rho_0,
		Real dt,
		Real h
	)
	{

		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= boundarySource.size()) return;
		if (boundaryAirFlag[pId] == true)return;
		if (boundarySolidFlag[pId] == true)return;

		if (boundaryRho[pId] > rho_0)
		{
			boundarySource[pId] += 1000000.0f * (boundaryRho[pId] - rho_0) / rho_0;
		}
	}

	template<typename Real, typename Coord>
	__global__ void SAISPH_AiiInLaplacian
	(
		DArray<Real> Aii,
		DArray<Coord> position,
		DArray<Real> pressure_,
		DArray<Real> rho,
		DArray<bool> boundaryFlag,
		DArrayList<int> neighborIds,
		DArray<Attribute> attribute,
		DArray<Coord>    b_pos,
		DArray<Real> b_rho,
		DArrayList<int> rv_neighbors,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		StressAwareKernel<Real> kernel,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Aii.size())return;

		if (!attribute[pId].isFluid())
		{
			Aii[pId] = EPSILON;
			return;
		}
		if (boundaryFlag[pId] == true)
		{
			Aii[pId] = EPSILON;
			return;
		}

		//Real fluid particle
		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();
		Real temp(0.0f);
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real rij = (position[pId] - position[j]).norm();
			if (rij > EPSILON && attribute[j].isFluid() && boundaryFlag[j] == false)
				//if (rij > EPSILON && attribute[j].isFluid())
			{
				Real s = rij * expf(-(pressure_[pId] + pressure_[j]) / 2.0f);
				Real dwij = kernel.Gradient(rij, h, s) / rij;
				temp += 8 * mass * dwij / (rho[pId] + rho[j]) / (rho[pId] + rho[j]);
			}
		}

		//virtual boundary particle
		List<int>& list_i_b = rv_neighbors[pId];
		int nbSize_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_b; ne_b++)
		{
			int j_b = list_i_b[ne_b];
			Real rij_b = (position[pId] - b_pos[j_b]).norm();
			if (rij_b > EPSILON && boundaryAirFlag[j_b] == false && boundarySolidFlag[j_b] == false)
			{
				Real s = rij_b * expf(-(pressure_[pId]) / 2.0f);
				Real dwij_b = kernel.Gradient(rij_b, h, s) / rij_b;
				temp += 8 * mass * dwij_b / (rho[pId] + b_rho[j_b]) / (rho[pId] + b_rho[j_b]);
			}
		}

		Aii[pId] = -temp;
	}
	
	template<typename Real, typename Coord>
	__global__ void SAISPH_BoundaryAiiInLaplacian
	(
		DArray<Real> b_Aii,
		DArray<Coord> position,
		DArray<Real> pressure_,
		DArray<Real> rho,
		DArray<bool> boundaryFlag,
		DArrayList<int> vr_neighbors,
		DArray<Attribute> attribute,
		DArray<Coord> b_pos,
		DArray<Real> b_rho,
		DArrayList<int> vv_neighbors,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		StressAwareKernel<Real> kernel,
		Real mass,
		Real h
	)
	{

		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= b_Aii.size())return;
		if (boundarySolidFlag[pId] == true)
		{
			b_Aii[pId] = EPSILON;
			return;
		}
		if (boundaryAirFlag[pId] == true)
		{
			b_Aii[pId] = EPSILON;
			return;
		}
		//Real fluid particle
		List<int>& list_i = vr_neighbors[pId];
		int nbSize = list_i.size();
		Real temp(0.0f);
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real rij = (b_pos[pId] - position[j]).norm();
			if (rij > EPSILON && attribute[j].isFluid() && boundaryFlag[j] == false)
			//if (rij > EPSILON && boundaryFlag[j] == false)
			{
				Real s = rij * expf(-(pressure_[j]) / 2.0f);
				Real dwij = kernel.Gradient(rij, h, s) / rij;
				temp += 8 * mass * dwij / (b_rho[pId] + rho[j]) / (b_rho[pId] + rho[j]);

			}
		}

		//virtual boundary particle
		List<int>& list_i_b = vv_neighbors[pId];
		int nbSize_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_b; ne_b++)
		{
			int j_b = list_i_b[ne_b];
			Real rij_b = (b_pos[pId] - b_pos[j_b]).norm();
			if (rij_b > EPSILON  && boundaryAirFlag[j_b] == false && boundarySolidFlag[j_b]== false)
			//if (rij_b > EPSILON && boundaryAirFlag[j_b] == false)
			{
				Real dwij_b = kernel.Gradient(rij_b, h, rij_b) / rij_b;
				temp += 8 * mass * dwij_b / (b_rho[pId] + b_rho[j_b]) / (b_rho[pId] + b_rho[j_b]);
			}
		}

		b_Aii[pId] = -temp;
	}

	template <typename Real>
	__global__ void SAISPH_CorrectAii
	(
		DArray<Real> Aii,
		DArray<bool> boundaryflag,
		Real max_Aii,
		Real c
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Aii.size()) return;
		//if (!attribute[pId].isFluid()) return;
		if (boundaryflag[pId] == true)return;
		Real temp = Aii[pId];
		if (temp < max_Aii * c)
		{
			temp = max_Aii;
		}
		Aii[pId] = temp;
	}

	template <typename Real>
	__global__ void  SAISPH_BoundaryCorrectAii
	(
		DArray<Real> b_Aii,
		DArray<bool> boundary_deleteFlag,
		Real max_Aii,
		Real c
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= b_Aii.size()) return;
		if (boundary_deleteFlag[pId] == true)return;
		Real temp = b_Aii[pId];
		if (temp < max_Aii * c)
		{
			temp = max_Aii * c;
		}
		b_Aii[pId] = temp;
	}

	template <typename Real, typename Coord>
	__global__ void SAISPH_AiiNeumannCorrect
	(
		DArray<Real> Aii,
		DArray<Attribute> attribute,
		DArray<bool> boundaryFlag,
		DArray<Coord> pos,
		DArray<Real> pressure_,
		DArray<Real> rho,
		DArrayList<int> neighbors,
		DArray<Coord> b_pos,
		DArray<Real> b_rho,
		DArrayList<int> rv_neighbors,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		StressAwareKernel<Real> kernel,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Aii.size())return;


		if (!attribute[pId].isFluid())return;
		if (boundaryFlag[pId] == true)return;
		List<int>& list_i_b = rv_neighbors[pId];
		int nbSize_i_b = list_i_b.size();

		Real temp(0.0f);
		int solidCounter = 0;
		for (int ne = 0; ne < nbSize_i_b; ne++)
		{
			int j_b = list_i_b[ne];
			Real r_ij_b = (pos[pId] - b_pos[j_b]).norm();
			if (r_ij_b > EPSILON && boundarySolidFlag[j_b] == true)
			{
				Real s = r_ij_b * expf(-(pressure_[pId]) / 2.0f);
				Real dwij = kernel.Gradient(r_ij_b, h, s) / r_ij_b;
				temp += 8 * mass * dwij / (rho[pId] + b_rho[j_b]) / (rho[pId] + b_rho[j_b]);
				solidCounter++;
			}
		}
		Real value = Aii[pId] + temp;
		Aii[pId] = value;
	}

	template <typename Real, typename Coord>
	__global__ void SAISPH_BoundaryAiiNeumannCorrect
	(
		DArray<Real> b_Aii,
		DArray<Attribute> attribute,
		DArray<bool> boundaryFlag,
		DArray<Coord> pos,
		DArray<Real> rho,
		DArrayList<int> vr_neighbors,
		DArray<Coord> b_pos,
		DArray<Real> b_rho,
		DArrayList<int> vv_neighbors,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		StressAwareKernel<Real> kernel,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= b_Aii.size())return;

		if (boundaryAirFlag[pId] == true)return;
		if (boundarySolidFlag[pId] == true)return;
		List<int>& list_i_b = vv_neighbors[pId];
		int nbSize_i_b = list_i_b.size();

		Real temp(0.0f);
		int solidCounter = 0;
		for (int ne = 0; ne < nbSize_i_b; ne++)
		{
			int j_b = list_i_b[ne];
			Real r_ij_b = (b_pos[pId] - b_pos[j_b]).norm();
			if (r_ij_b > EPSILON && boundarySolidFlag[j_b] == true)
			{
				Real dwij = kernel.Gradient(r_ij_b, h, r_ij_b) / r_ij_b;
				temp += 8 * mass * dwij / (b_rho[pId] + b_rho[j_b]) / (b_rho[pId] + b_rho[j_b]);
				solidCounter++;
			}
		}

		Real value = b_Aii[pId] + temp;
		b_Aii[pId] = value;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_LaplacianPressure
	(
		DArray<Real> Ax,
		DArray<Real> Aii,
		DArray<Attribute> attribute,
		DArray<Real> pressure,
		DArray<Coord> pos,
		DArray<Real> pressure_,
		DArrayList<int> neighborIds,
		DArray<bool> boundaryFlag,
		DArray<Real> rho,
		DArray<Real> Aii_b,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArray<Real> pressure_b,
		DArray<Coord> pos_b,
		DArrayList<int> rv_neighbors,
		DArray<Real> rho_b,
		StressAwareKernel<Real> kernel,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Ax.size()) return;

		if (!attribute[pId].isFluid())
		{
			Ax[pId] = 0.0f;
			pressure[pId] = 0.0f;
			return;
		}
		if (boundaryFlag[pId] == true)
		{
			Ax[pId] = 0.0f;
			pressure[pId] = 0.0f;
			return;
		}

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		Real temp = 0.0f;
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real rij = (pos[pId] - pos[j]).norm();

			if (rij > EPSILON && attribute[j].isFluid() && boundaryFlag[j] == false)
			{
				Real s = rij * expf(-(pressure_[pId] + pressure_[j]) / 2.0f);
				Real dwij = kernel.Gradient(rij, h, s) / rij;
				temp += 8 * mass * dwij / (rho[pId] + rho[j]) / (rho[pId] + rho[j]) * pressure[j];
			}
		}

		List<int>& list_i_b = rv_neighbors[pId];
		int nbSize_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_b; ne_b++)
		{
			int j_b = list_i_b[ne_b];
			Real rij_b = (pos[pId] - pos_b[j_b]).norm();
			if (rij_b > EPSILON && boundaryAirFlag[j_b] == false && boundarySolidFlag[j_b] == false)
			{
				Real s = rij_b * expf(-(pressure_[pId]) / 2.0f);
				Real dwij = kernel.Gradient(rij_b, h, s) / rij_b;
				temp += 8 * mass * dwij / (rho[pId] + rho_b[j_b]) / (rho[pId] + rho_b[j_b]) * pressure_b[j_b];
			}
		}
		Ax[pId] = Aii[pId] * pressure[pId] + temp;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_BoundaryLaplacianPressure
	(
		DArray<Real> Ax_b,
		DArray<Real> Aii_b,
		DArray<Attribute> attribute,
		DArray<Real> pressure,
		DArray<Coord> pos,
		DArray<Real> pressure_,
		DArrayList<int> vr_neighbors,
		DArray<bool> boundaryFlag,
		DArray<Real> rho,
		DArray<Real> Aii,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArray<Real> pressure_b,
		DArray<Coord> pos_b,
		DArrayList<int> vv_neighbors,
		DArray<Real> rho_b,
		StressAwareKernel<Real> kernel,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Ax_b.size()) return;
		if (boundaryAirFlag[pId] == true)
		{
			Ax_b[pId] = 0.0f;
			pressure_b[pId] = 0.0f;
			return;
		}
		if (boundarySolidFlag[pId] == true)
		{
			//Ax_b[pId] = 0.0f;
			//pressure_b[pId] = 0.0f;
			return;
		}

		List<int>& list_i = vr_neighbors[pId];
		int nbSize = list_i.size();

		Real temp = 0.0f;
		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real rij = (pos_b[pId] - pos[j]).norm();

			//if (rij > EPSILON && attribute[j].isFluid() && boundaryFlag[j] == false)
			if (rij > EPSILON && boundaryFlag[j] == false)
			{
				Real s = rij * expf(-(pressure_[j]) / 2.0f);
				Real dwij = kernel.Gradient(rij, h, s) / rij;
				temp += 8 * mass * dwij / (rho_b[pId] + rho[j]) / (rho_b[pId] + rho[j]) * pressure[j];
			}
		}

		List<int>& list_i_b = vv_neighbors[pId];
		int nbSize_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_b; ne_b++)
		{
			int j_b = list_i_b[ne_b];
			Real rij_b = (pos_b[pId] - pos_b[j_b]).norm();
			if (rij_b > EPSILON && boundaryAirFlag[j_b] == false && boundarySolidFlag[j_b] == false)
			{
				Real dwij = kernel.Gradient(rij_b, h, rij_b) / rij_b;
				temp += 8 * mass * dwij / (rho_b[pId] + rho_b[j_b]) / (rho_b[pId] + rho_b[j_b]) * pressure_b[j_b];
			}
		}

		Ax_b[pId] = Aii_b[pId] * pressure_b[pId] + temp;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_GradientPressure
	(
		DArray<Coord> gradient,
		DArray<Real> pressure,
		DArray<Coord> velocity,
		DArray<Coord> pos,
		DArray<Real> pressure_,
		DArray<Attribute> attribute,
		DArray<bool> boundaryFlag,
		DArray<Real> b_pressure,
		DArray<Coord> b_pos,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArrayList<int> neighbors,
		DArrayList<int> rv_neighbors,
		DArray<Real> rho,
		DArray<Real> b_rho,
		StressAwareKernel<Real> kernel,
		Real rho_0,
		Real dt,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= gradient.size()) return;
		if (!attribute[pId].isFluid()) return;

		Coord value(0.0f);
		List<int> list_i = neighbors[pId];
		int nbSize_i = list_i.size();
		for (int ne = 0; ne < nbSize_i; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (pos[pId] - pos[j]).norm();
			//if ((r_ij > EPSILON) && (attribute[j].isFluid()) && (boundaryFlag[j] == false))
			if ((r_ij > EPSILON) && (boundaryFlag[j] == false))
			{
				Real s = r_ij * expf(-(pressure_[j] + pressure_[pId]) / 2.0f);
				//value += (mass / rho[j]) * pressure[j] * (pos[pId] - pos[j]) * kernel.Gradient(r_ij, h, s) / r_ij;
				value += (mass / rho_0) * pressure[j] * (pos[pId] - pos[j]) * kernel.Gradient(r_ij, h, s) / r_ij;
			}
		}

		List<int> list_i_b = rv_neighbors[pId];
		int nbSize_i_b = list_i_b.size();
		for (int ne_b = 0; ne_b < nbSize_i_b; ne_b++)
		{
			int j_b = list_i_b[ne_b];
			Real r_ij_b = (pos[pId] - b_pos[j_b]).norm();
			if ((r_ij_b > EPSILON) && (boundaryAirFlag[j_b] == false) && (boundarySolidFlag[j_b] == false))
			{
				Real s = r_ij_b * expf(-(pressure_[pId]) / 2.0f);
				//value += (mass / b_rho[j_b]) * b_pressure[j_b] * (pos[pId] - b_pos[j_b]) * kernel.Gradient(r_ij_b, h, s) / r_ij_b;
				value += (mass / rho_0) * b_pressure[j_b] * (pos[pId] - b_pos[j_b]) * kernel.Gradient(r_ij_b, h, s) / r_ij_b;
			}
		}

		value = value * dt / rho[pId];
		gradient[pId] = value;
		velocity[pId] -= value;
	}

	template<typename Coord, typename Real>
	__global__ void SAISPH_GradientNearSolid
	(
		DArray<Coord> gradientComp,
		DArray<Real> pressure,
		DArray<Coord> velocity,
		DArray<Coord> pos,
		DArray<Attribute> attribute,
		DArray<Coord> norm,
		DArray<bool> boundaryFlag,
		DArray<Real> b_pressure,
		DArray<Coord> b_pos,
		DArray<bool> boundaryAirFlag,
		DArray<bool> boundarySolidFlag,
		DArrayList<int> neighbors,
		DArrayList<int> rv_neighbors,
		DArray<Real> rho,
		DArray<Real> b_rho,
		SA_CubicKernel<Real> kernel,
		Real rho_0,
		Real dt,
		Real mass,
		Real h
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= gradientComp.size()) return;

		gradientComp[pId] = Coord(0.0f);

		List<int> list_i = neighbors[pId];
		int nbSize_i = list_i.size();
		Coord& vel_i = velocity[pId];
		Coord dvij(0.0f);
		for (int ne = 0; ne < nbSize_i; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (pos[pId] - pos[j]).norm();
			Real weight = kernel.Weight(r_ij, h, 0.0f);
			if (!attribute[j].isFluid())
			{
				Coord nij = (pos[pId] - pos[j]);
				if (nij.norm() > EPSILON)
				{
					nij = nij.normalize();
				}
				else
				{
					nij = Coord(1.0f, 0.0f, 0.0f);
				}

				Coord normal_j = norm[j];
				Coord dVel = vel_i - velocity[j];
				Real magNVel = dVel.dot(normal_j);
				Coord nVel = magNVel * normal_j;
				Coord tVel = dVel - nVel;

				if (magNVel < -EPSILON)
				{
					dvij = nij.dot(nVel + 0.01 * dVel) * weight * nij;
				}
				else
				{
					dvij += nij.dot(0.1 * nVel + 0.01 * dVel) * weight * nij;
				}
			}

			gradientComp[pId] = 2 * dt * dvij / rho[pId];
			velocity[pId] -= gradientComp[pId];
		}
	}

	template<typename Real, typename Coord>
	__global__ void SAISPH_artifical_viscosity(
		DArray<Coord> pos,
		DArray<Real> pressure_,
		DArrayList<int> neighborIds,
		DArray<Attribute> attribute,
		DArray<Real> density,
		DArray<Real> pressure,
		DArray<Coord> velocity,
		SA_CubicKernel<Real> kernel,
		Real pressure_max,
		Real pressure_min,
		Real rho_max,
		Real h,
		Real dt,
		Real mass,
		Real factor
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size()) return;
		if (!attribute[pId].isFluid()) return;

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		Coord viscous_term(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			if (attribute[j].isFluid()) {
				//J. J. Monaghan 1997 , SPH and Riemann Solvers
				Coord x_ij = pos[pId] - pos[j];
				Real x_ij_v_ij = x_ij.dot(velocity[pId] - velocity[j]);
				Real r_ij = x_ij.norm();
				if (r_ij > EPSILON) {
					Real s = r_ij * expf(-(pressure_[j] + pressure_[pId]) / 2.0f);
					Coord gradW = kernel.Gradient(r_ij, h, s) * x_ij / r_ij;
					Real mu_ij = h * x_ij_v_ij / (r_ij * r_ij + 0.01 * h * h);
					//viscous_term += mass * (-0.005f * h * 1000.0f* mu_ij + 2.0f * mu_ij * mu_ij) / (density[pId] + density[j]) * gradW;
					if (x_ij_v_ij < 0)viscous_term += mass * (-3000.0f * mu_ij + 2.0f * mu_ij * mu_ij) / (density[pId] + density[j]) * gradW;
				}
			}
		}

		Real fac(0.0f);
		Real p_ratio = ((pressure[pId] - pressure_min) / (pressure_max - pressure_min));
		fac = 0.1f + 0.7f * (1.0f - (p_ratio));

		velocity[pId] -= factor * fac * dt * viscous_term;
	}

	template<typename TDataType>
	void SAISPHsolver<TDataType>::constrain()
	{
		int num = this->inPosition()->size();

		cudaDeviceSynchronize();

		if (m_Gp.size() != num)
		{
			ArraysResize();
		}

		Real dt = this->inTimeStep()->getData();
		//uint pDims = cudaGridSize(inPosition()->size(), BLOCK_SIZE);
		Real h = this->varSmoothingLength()->getData();

		if (this->inParticleAttribute()->isEmpty()
			|| this->inParticleAttribute()->size() != num
			|| this->inBoundaryNorm()->size() != num
			)
		{
			this->inParticleAttribute()->allocate();
			this->inParticleAttribute()->resize(num);
			cuExecute(num, SAISPH_AttributeInit,
				this->inParticleAttribute()->getData());
			this->inBoundaryNorm()->resize(num);
			this->inBoundaryNorm()->reset();
		}

		m_ak_summation->update();
		m_particleMass = m_ak_summation->getParticleMass();
		Real rho_factor = m_ak_summation->getScalingFactor();

		Real restDensity = this->varRestDensity()->getValue();


		Real MaxDensity = m_reduce->maximum(
			m_ak_summation->outDensity()->getData().begin(),
			m_ak_summation->outDensity()->getData().size()
		);

		Real MinDensity = m_reduce->minimum(
			m_ak_summation->outDensity()->getData().begin(),
			m_ak_summation->outDensity()->getData().size()
		);

		std::cout << "MaxDensity:" << MaxDensity << "MinDensity:" << MinDensity << std::endl;

		Real Density_t = restDensity * 0.7f;
		if (MaxDensity < Density_t) Density_t = MaxDensity * 0.7f;
		if (MinDensity > Density_t) Density_t = MinDensity * 1.01f;

		cuExecute(num, SAISPH_BoundaryParticleDetection,
			m_BoundaryFlag,
			this->inParticleAttribute()->getData(),
			this->inNeighborIds()->getData(),
			m_ak_summation->outDensity()->getData(),
			Density_t);


		//generate boundary particles
		Real gridSize = this->varSamplingDistance()->getData();

		int node_num = num * 33;

		if (m_anchorPoint.size() != node_num)
		{
			m_anchorPoint.resize(node_num);
			m_anchorPointFlag.resize(node_num);
			m_anchorPointCodes.resize(node_num);
			m_nonRepeatedCount.resize(node_num);
		}

		Reduction<Coord> reduce_coord;
		Coord hiBound = reduce_coord.maximum(this->inPosition()->getData().begin(), this->inPosition()->getData().size());
		Coord loBound = reduce_coord.minimum(this->inPosition()->getData().begin(), this->inPosition()->getData().size());

		//delete reduce_coord;

		int padding = 2;
		hiBound += Coord((padding + 1) * gridSize);
		loBound -= Coord(padding * gridSize);

		loBound[0] = (Real)((int)(floor(loBound[0] / gridSize))) * gridSize;
		loBound[1] = (Real)((int)(floor(loBound[1] / gridSize))) * gridSize;
		loBound[2] = (Real)((int)(floor(loBound[2] / gridSize))) * gridSize;

		Coord origin = Coord(loBound[0], loBound[1], loBound[2]);

		//get virtual position for boundary points
		cuExecute(num,
			SAISPH_AnchorNeighbor_33,
			m_anchorPoint,
			m_anchorPointFlag,
			this->inPosition()->getData(),
			m_BoundaryFlag,
			origin,
			gridSize
		);

		//get virtual points Morton Codes
		cuExecute(m_anchorPointCodes.size(),
			SAISPH_PositionMortonCode,
			m_anchorPointCodes,
			m_anchorPointFlag,
			m_anchorPoint,
			gridSize
		);

		thrust::sort(thrust::device, m_anchorPointCodes.begin(), m_anchorPointCodes.begin() + m_anchorPointCodes.size());

		cuExecute(m_anchorPointCodes.size(),
			SAISPH_RepeatedElementSearch,
			m_anchorPointCodes,
			m_nonRepeatedCount
		);

		int num_b = thrust::reduce(thrust::device, m_nonRepeatedCount.begin(), m_nonRepeatedCount.begin() + m_nonRepeatedCount.size(), int(0), thrust::plus<uint32_t>()) - 1;

		thrust::exclusive_scan(thrust::device, m_nonRepeatedCount.begin(), m_nonRepeatedCount.begin() + m_nonRepeatedCount.size(), m_nonRepeatedCount.begin());

		m_candidateCodes.resize(num_b);

		cuExecute(m_anchorPointCodes.size(),
			SAISPH_nonRepeatedElementsCal,
			m_candidateCodes,
			m_anchorPointCodes,
			m_nonRepeatedCount
		);

		cudaDeviceSynchronize();

		ArraysResize_Boundary(num_b);

		cuExecute(num_b,
			SAISPH_MortonCodeToPosition,
			m_boundary_position,
			m_candidateCodes,
			gridSize
		);

		if (num_b != this->outBoundaryParticles()->size())
		{
			this->outBoundaryParticles()->resize(num_b);
		}

		cuExecute(num_b,
			SAISPH_CopyToBoundaryPos,
			this->outBoundaryParticles()->getData(),
			origin,
			m_boundary_position
		);

		{
			m_anchorPoint.resize(0);
			m_anchorPointFlag.resize(0);
			m_anchorPointCodes.resize(0);
			m_nonRepeatedCount.resize(0);
		}

		m_vv_nbrQuery->update();
		m_vr_nbrQuery->update();
		m_rv_nbrQuery->update();

		//计算基于虚粒子+实粒子的粒子密度
		cuExecute(num_b, SAISPH_compute_b_rho,
			m_boundary_position,
			this->outVRNeighborIds()->getData(),
			this->inPosition()->getData(),
			m_BoundaryFlag,
			m_boundary_rho,
			m_boundary_rho_projection,
			this->inParticleAttribute()->getData(),
			this->outVVNeighborIds()->getData(),
			kernel_mass,
			m_particleMass,
			rho_factor,
			h
		);

		//基于粒子密度删除过近的粒子
		cuExecute(num_b,
			SAISPH_DeleteNearParticles,
			m_boundary_position,
			m_boundary_deleteFlag,
			m_boundary_rho,
			restDensity,
			MaxDensity,
			h
		);

		std::cout << "格点生成虚粒子数目" << num_b << std::endl;

		//int num_b_filtered = thrust::reduce(thrust::device, m_boundary_deleteFlag.begin(), m_boundary_deleteFlag.begin() + m_boundary_deleteFlag.size(), int(0), thrust::plus<bool>());

		//std::cout << "参与计算的虚粒子数目" << num_b_filtered << std::endl;

		cuExecute(num_b, SAISPH_compute_b_rho_filtered,
			m_boundary_position,
			m_boundary_deleteFlag,
			this->outVRNeighborIds()->getData(),
			this->inPosition()->getData(),
			m_BoundaryFlag,
			m_boundary_rho,
			m_boundary_rho_projection,
			this->inParticleAttribute()->getData(),
			this->outVVNeighborIds()->getData(),
			kernel_mass,
			m_particleMass,
			rho_factor,
			h
		);
		//m_anchorPoint.clear();


		Real MaxBoundaryDensity =
			m_boundary_reduce->maximum(
				m_boundary_rho.begin(),
				m_boundary_rho.size()
			);

		cuExecute(num_b,
			SAISPH_InitBoundaryParticles,
			m_boundary_position,
			m_boundary_deleteFlag,
			m_boundary_rho_projection,
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			this->outVRNeighborIds()->getData(),
			this->inPosition()->getData(),
			m_ak_summation->outDensity()->getData(),
			this->inParticleAttribute()->getData(),
			m_BoundaryFlag,
			kernel_mass,
			m_particleMass,
			rho_factor,
			h,
			MaxBoundaryDensity * 0.1f,
			0.999f
		);

		//get boundary particle velocity by csph projection
		cuExecute(num_b,
			SAISPH_SphBoundaryVelocity,
			m_boundary_velocity,
			this->inVelocity()->getData(),
			m_boundary_position,
			this->inPosition()->getData(),
			this->outVRNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			m_ak_summation->outDensity()->getData(),
			kernel_mass,
			m_particleMass,
			h
		);


		cuExecute(num,
			SAISPH_SolidVelocityReset,
			this->inVelocity()->getData(),
			this->inParticleAttribute()->getData()
		);

		Real pressure_max = m_reduce->maximum(m_pressure.begin(), m_pressure.size());
		if (fabs(pressure_max) < EPSILON) pressure_max = EPSILON;
		Real pressure_min = m_reduce->minimum(m_pressure.begin(), m_pressure.size());
		if (-pressure_min > pressure_max) pressure_max = -pressure_min;

		cuExecute(num,
			SAISPH_Init_pressure_,
			m_pressure,
			pressure_max,
			m_pressure_
		);

		//isph module

		m_source.reset();
		m_boundarySource.reset();

		//sourceTerm
		cuExecute(num,
			SAISPH_SourceTerm,
			m_source,
			this->inVelocity()->getData(),
			this->inPosition()->getData(),
			m_pressure_,
			m_ak_summation->outDensity()->getData(),
			m_BoundaryFlag,
			this->inNeighborIds()->getData(),
			this->inParticleAttribute()->getData(),
			m_boundary_position,
			m_boundary_velocity,
			m_boundary_rho_projection,
			this->outRVNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			this->inBoundaryNorm()->getData(),
			kernel,
			restDensity,
			m_particleMass,
			dt,
			h
		);

		cuExecute(num_b,
			SAISPH_BoundarySourceTerm,
			m_boundarySource,
			this->inVelocity()->getData(),
			this->inPosition()->getData(),
			m_pressure_,
			m_ak_summation->outDensity()->getData(),
			m_BoundaryFlag,
			this->outVRNeighborIds()->getData(),
			this->inParticleAttribute()->getData(),
			m_boundary_position,
			m_boundary_velocity,
			m_boundary_rho_projection,
			this->outVVNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			this->inBoundaryNorm()->getData(),
			kernel,
			restDensity,
			m_particleMass,
			dt,
			h
		);

		//source term correct
		cuExecute(num,
			SAISPH_DensityCompensate,
			m_source,
			m_ak_summation->outDensity()->getData(),
			this->inParticleAttribute()->getData(),
			m_BoundaryFlag,
			restDensity,
			dt,
			h
		);
		cuExecute(num_b,
			SAISPH_BoundaryDensityCompensate,
			m_boundarySource,
			m_boundary_rho_projection,
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			restDensity,
			dt,
			h
		);

		m_r.reset();
		m_Ax.reset();
		m_boundary_r.reset();
		m_boundaryAx.reset();

		if (m_boundary_pressure.size() != points_num_old)
			m_boundary_pressure.reset();
		points_num_old = m_boundary_pressure.size();

		m_boundary_pressure.reset();

		cuExecute(num,
			SAISPH_AiiInLaplacian,
			m_Aii,
			this->inPosition()->getData(),
			m_pressure_,
			m_ak_summation->outDensity()->getData(),
			m_BoundaryFlag,
			this->inNeighborIds()->getData(),
			this->inParticleAttribute()->getData(),
			m_boundary_position,
			m_boundary_rho,
			this->outRVNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			kernel,
			m_particleMass,
			h
		);
		cuExecute(num_b,
			SAISPH_BoundaryAiiInLaplacian,
			m_boundaryAii,
			this->inPosition()->getData(),
			m_pressure_,
			m_ak_summation->outDensity()->getData(),
			m_BoundaryFlag,
			this->outVRNeighborIds()->getData(),
			this->inParticleAttribute()->getData(),
			m_boundary_position,
			m_boundary_rho,
			this->outVVNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			kernel,
			m_particleMass,
			h
		);

		if ((frag_number < 3) || abs(max_Aii) < EPSILON)
		{
			Real max_f_Aii = m_reduce->maximum(m_Aii.begin(), m_Aii.size());
			Real max_b_Aii = m_boundary_reduce->maximum(m_boundaryAii.begin(), m_boundaryAii.size());

			max_Aii = max_f_Aii > max_Aii ? max_f_Aii : max_Aii;
			max_Aii = max_b_Aii > max_Aii ? max_b_Aii : max_Aii;
		}
		frag_number += 1;

		cuExecute(num, SAISPH_CorrectAii,
			m_Aii,
			//this->inParticleAttribute()->getData(),
			m_BoundaryFlag,
			max_Aii,
			1.0f
		);
		cuExecute(num_b, SAISPH_BoundaryCorrectAii,
			m_boundaryAii,
			m_boundary_deleteFlag,
			max_Aii,
			1.0f
		);

		//NeumannCorrect
		cuExecute(num,
			SAISPH_AiiNeumannCorrect,
			m_Aii,
			this->inParticleAttribute()->getData(),
			m_BoundaryFlag,
			this->inPosition()->getData(),
			m_pressure_,
			m_ak_summation->outDensity()->getData(),
			this->inNeighborIds()->getData(),
			m_boundary_position,
			m_boundary_rho,
			this->outRVNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			kernel,
			m_particleMass,
			h
		);

		cuExecute(num_b,
			SAISPH_BoundaryAiiNeumannCorrect,
			m_boundaryAii,
			this->inParticleAttribute()->getData(),
			m_BoundaryFlag,
			this->inPosition()->getData(),
			m_ak_summation->outDensity()->getData(),
			this->outVRNeighborIds()->getData(),
			m_boundary_position,
			m_boundary_rho,
			this->outVVNeighborIds()->getData(),
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			kernel,
			m_particleMass,
			h
		);

		//LaplacianPressure
		cuExecute(num,
			SAISPH_LaplacianPressure,
			m_Ax,
			m_Aii,
			this->inParticleAttribute()->getData(),
			m_pressure,
			this->inPosition()->getData(),
			m_pressure_,
			this->inNeighborIds()->getData(),
			m_BoundaryFlag,
			m_ak_summation->outDensity()->getData(),
			m_boundaryAii,
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			m_boundary_pressure,
			m_boundary_position,
			this->outRVNeighborIds()->getData(),
			m_boundary_rho,
			kernel,
			m_particleMass,
			h
		);

		cuExecute(num_b,
			SAISPH_BoundaryLaplacianPressure,
			m_boundaryAx,
			m_boundaryAii,
			this->inParticleAttribute()->getData(),
			m_pressure,
			this->inPosition()->getData(),
			m_pressure_,
			this->outVRNeighborIds()->getData(),
			m_BoundaryFlag,
			m_ak_summation->outDensity()->getData(),
			m_Aii,
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			m_boundary_pressure,
			m_boundary_position,
			this->outVVNeighborIds()->getData(),
			m_boundary_rho,
			kernel,
			m_particleMass,
			h
		);

		Function2Pt::subtract(m_r, m_source, m_Ax);

		Function2Pt::subtract(m_boundary_r, m_boundarySource, m_boundaryAx);

		m_p.assign(m_r);
		m_boundary_p.assign(m_boundary_r);

		Real rr = m_arithmetic->Dot(m_r, m_r) + m_boundary_arithmetic->Dot(m_boundary_r, m_boundary_r);

		Real err = num + num_b > 0 ? sqrt(rr / (num + num_b)) : 0.0f;
		Real max_err = err;
		if (abs(max_err) < EPSILON) max_err = EPSILON;

		unsigned int iter = 0;
		Real threshold = this->varResiThreshold()->getValue();

		while ((iter < 500) && (err / max_err > threshold) && (err > threshold))
		{
			iter++;
			m_Ax.reset();
			m_boundaryAx.reset();

			cuExecute(num,
				SAISPH_LaplacianPressure,
				m_Ax,
				m_Aii,
				this->inParticleAttribute()->getData(),
				m_p,
				this->inPosition()->getData(),
				m_pressure_,
				this->inNeighborIds()->getData(),
				m_BoundaryFlag,
				m_ak_summation->outDensity()->getData(),
				m_boundaryAii,
				m_boundaryAirFlag,
				m_boundarySolidFlag,
				m_boundary_p,
				m_boundary_position,
				this->outRVNeighborIds()->getData(),
				m_boundary_rho,
				kernel,
				m_particleMass,
				h
			);

			cuExecute(num_b,
				SAISPH_BoundaryLaplacianPressure,
				m_boundaryAx,
				m_boundaryAii,
				this->inParticleAttribute()->getData(),
				m_p,
				this->inPosition()->getData(),
				m_pressure_,
				this->outVRNeighborIds()->getData(),
				m_BoundaryFlag,
				m_ak_summation->outDensity()->getData(),
				m_Aii,
				m_boundaryAirFlag,
				m_boundarySolidFlag,
				m_boundary_p,
				m_boundary_position,
				this->outVVNeighborIds()->getData(),
				m_boundary_rho,
				kernel,
				m_particleMass,
				h
			);

			float alpha = rr / (m_arithmetic->Dot(m_p, m_Ax) + m_boundary_arithmetic->Dot(m_boundary_p, m_boundaryAx));
			Function2Pt::saxpy(m_pressure, m_p, m_pressure, alpha);
			Function2Pt::saxpy(m_boundary_pressure, m_boundary_p, m_boundary_pressure, alpha);

			Function2Pt::saxpy(m_r, m_Ax, m_r, -alpha);
			Function2Pt::saxpy(m_boundary_r, m_boundaryAx, m_boundary_r, -alpha);

			Real rr_old = rr;

			rr = m_arithmetic->Dot(m_r, m_r) + m_boundary_arithmetic->Dot(m_boundary_r, m_boundary_r);

			Real beta = rr / rr_old;
			Function2Pt::saxpy(m_p, m_p, m_r, beta);
			Function2Pt::saxpy(m_boundary_p, m_boundary_p, m_boundary_r, beta);

			err = sqrt(rr / (num + num_b));
		}

		std::cout << "*STABLE-ISPH::Solver::Iteration:" << iter << "||RelativeError:" << err / max_err * 100 << "%" << std::endl;

		iter_sum += iter;

		std::cout << "ave iters" << Real(iter_sum / frag_number) << std::endl;

		m_GpNearSolid.reset();

		cuExecute(num,
			SAISPH_GradientPressure,
			m_Gp,
			m_pressure,
			this->inVelocity()->getData(),
			this->inPosition()->getData(),
			m_pressure_,
			this->inParticleAttribute()->getData(),
			m_BoundaryFlag,
			m_boundary_pressure,
			m_boundary_position,
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			this->inNeighborIds()->getData(),
			this->outRVNeighborIds()->getData(),
			m_ak_summation->outDensity()->getData(),
			m_boundary_rho,
			kernel,
			restDensity,
			dt,
			m_particleMass,
			h
		);



		cuExecute(num,
			SAISPH_GradientNearSolid,
			m_GpNearSolid,
			m_pressure,
			this->inVelocity()->getData(),
			this->inPosition()->getData(),
			this->inParticleAttribute()->getData(),
			this->inBoundaryNorm()->getData(),
			m_BoundaryFlag,
			m_boundary_pressure,
			m_boundary_position,
			m_boundaryAirFlag,
			m_boundarySolidFlag,
			this->inNeighborIds()->getData(),
			this->outRVNeighborIds()->getData(),
			m_ak_summation->outDensity()->getData(),
			m_boundary_rho,
			kernel_mass,
			restDensity,
			dt,
			m_particleMass,
			h
		);

		////artificalViscosity, control zero-energy mode
		pressure_max = m_reduce->maximum(m_pressure.begin(), m_pressure.size());
		if (fabs(pressure_max) < EPSILON) pressure_max = EPSILON;
		pressure_min = m_reduce->minimum(m_pressure.begin(), m_pressure.size());

		cuExecute(num, SAISPH_artifical_viscosity,
			this->inPosition()->getData(),
			m_pressure_,
			this->inNeighborIds()->getData(),
			this->inParticleAttribute()->getData(),
			m_ak_summation->outDensity()->getData(),
			m_pressure,
			this->inVelocity()->getData(),
			kernel_mass,
			pressure_max,
			pressure_min,
			MaxDensity,
			h,
			dt,
			m_particleMass,
			Real(1.0f)
		);

	}

	template<typename TDataType>
	bool SAISPHsolver<TDataType>::initializeImpl()
	{
		cuSynchronize();
		int num = this->inPosition()->size();

		m_reduce = Reduction<float>::Create(num);
		m_arithmetic = Arithmetic<float>::Create(num);

		return true;
	}

	DEFINE_CLASS(SAISPHsolver);
}