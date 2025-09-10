#include <cuda_runtime.h>
#include "SA_DFSPHsolver.h"
#include "Node.h"
#include "Field.h"
#include "StressAwareKernelSummation.h"
#include "Collision/Attribute.h"
#include "StressAwareKernel.h"
#include "Algorithm/Function2Pt.h"
#include "Algorithm/CudaRand.h"

#include"Collision/NeighborPointQuery.h"

#ifdef STRESSAWARE
#define Kernel StressAwareKernel
#else
#define Kernel SA_CubicKernel
#endif // STRESSAWARE


namespace dyno
{
	IMPLEMENT_TCLASS(SA_DFSPHsolver, TDataType)

		template<typename TDataType>
	SA_DFSPHsolver<TDataType>::SA_DFSPHsolver()
		: ConstraintModule()
	{
		//this->varSamplingDistance()->setValue(Real(0.005));
		//this->varRestDensity()->setValue(Real(1000));
		//this->varSmoothingLength()->setValue(Real(0.010));

		mSummation = std::make_shared<StressAwareKernelSummation<TDataType>>();
		this->varSmoothingLength()->connect(mSummation->inSmoothingLength());
		this->varRestDensity()->connect(mSummation->varRestDensity());
		this->varSamplingDistance()->connect(mSummation->inSamplingDistance());
		this->inPosition()->connect(mSummation->inPosition());
		this->inNeighborIds()->connect(mSummation->inNeighborIds());

	}

	template<typename TDataType>
	SA_DFSPHsolver<TDataType>::~SA_DFSPHsolver()
	{
		mKappa_r.clear();
		mKappa_v.clear();
		mAlpha.clear();
		mPara.clear();
		mPressureEstimation.clear();
		mAccel.clear();
		mDensityAdv.clear();
		m_resv.clear();
		m_resr.clear();
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_AlphaCompute(
		DArray<Real> alpha,
		DArray<Coord> pos,
		DArray<Real> density,
		DArrayList<int> neighbors,
		Real density_0,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= alpha.size()) return;

		Coord pos_i = pos[pId];

		Real inv_alpha_i = Real(0);
		Coord grad_ci(0);

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				grad_ci += g;
				inv_alpha_i += g.dot(g);
			}
		}

		inv_alpha_i += grad_ci.dot(grad_ci);

		Real density_i = density[pId];// > density_0 ? density[pId] : density_0;

		if (inv_alpha_i < EPSILON)// || density_i < 0.95 * density_0 || density_i > 1.05 * density_0)
#ifdef DFSPH_ORIGIN
			alpha[pId] = 0.0f;
#else
			alpha[pId] = 1.0f / EPSILON;
#endif
		else
			alpha[pId] = density_i / (inv_alpha_i);

		//printf("density:%f,alpha:%f\n", density[pId], alpha[pId]);
	}

	template<typename Real>
	__global__ void SADFSPH_AlphaCorrect
	(
		DArray<Real> alpha_c,
		Real alpha_min,
		Real c
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= alpha_c.size()) return;


		alpha_c[pId] = c * alpha_min;

	}

	template<typename Real>
	__global__ void SADFSPH_PressureEstimation
	(
		DArray<Real> Pressure,
		DArray<Real> density,
		Real density_0
	)
	{
		//Estimate with Tait’s equation, 
		//Markus Becker and Matthias Teschner, Weakly compressible SPH for free surface flows, 2007
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Pressure.size()) return;

		Pressure[pId] = 1119 * (pow((density[pId] / density_0), 7) - 1);
	}

	template<typename Real>
	__global__ void SADFSPH_Paramaters_Init
	(
		DArray<Real> Para,
		DArray<Real> PressureEst,
		Real pressure_max
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId > Para.size()) return;

		Para[pId] = PressureEst[pId] / pressure_max * 10.0f;
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_PredictDensityAdv_Div(
		DArray<Real> densityAdv,
		DArray<Coord> pos,
		DArray<Real> density,
		DArray<Coord> vel,
		DArrayList<int> neighbors,
		Real density_0,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= densityAdv.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real densityAdv_i(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				densityAdv_i += (vel[pId] - vel[j]).dot(g);
			}
		}

#ifdef DFSPH_ORIGIN
		densityAdv_i = densityAdv_i > 0.0f ? densityAdv_i : 0.0f;

		if (nbSize < 5) densityAdv_i = 0.0f;
#endif 

		densityAdv[pId] = densityAdv_i;

	}

	template<typename Real>
	__global__ void SADFSPH_DivergenceKappa(
		DArray<Real> kappa_v,
		DArray<Real> densityAdv,
		DArrayList<int> neighbors,
		DArray<Real> alpha,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= kappa_v.size()) return;

		kappa_v[pId] = densityAdv[pId] * alpha[pId] / dt;
		//if (fabs(densityAdv[pId]) > EPSILON) printf("densityAdv_i:%f,kappa_i:%f,alpha_i:%f,1/dt:1/dt\n", densityAdv[pId], kappa_v[pId], alpha[pId], 1 / dt);

		//List<int>& list_i = neighbors[pId];
		/*int nbSize = list_i.size();
		if (nbSize < 20)
			kappa_v[pId] = 0.0f;*/
	}

	template<typename Real, typename Coord>
	__global__ void SADFSAPH_KappaAccel(
		DArray<Coord> accelArr,
		DArray<Coord> pos,
		DArray<Real> kappaArr,
		DArray<Real> density,
		DArrayList<int> neighbors,
		Real density_0,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= accelArr.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real Density_i = density[pId];// > density_0 ? density[pId] : density_0;
		Coord Accel(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				Real Density_j = density[j];// > density_0 ? density[j] : density_0;
				Accel -= (kappaArr[pId] / Density_i + kappaArr[j] / Density_j) * g;
				//Accel -= (kappaArr[j] / Density_j) * g;
				//Accel -= (kappaArr[j] + kappaArr[pId]) / density_0 * g;
			}

		}

		accelArr[pId] = Accel;
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_UpdateDivKappa(
		DArray<Real> mKappa_v,
		DArray<Real> res_v,
		DArray<Coord> pos,
		DArray<Coord> accelArr,
		DArray<Real> densityAdv,
		DArray<Real> Alpha,
		Real Alpha_min,
		DArrayList<int> neighbors,
		Real density_0,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= mKappa_v.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real aij_pj(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				aij_pj += (accelArr[pId] - accelArr[j]).dot(g);
			}
		}

		aij_pj *= dt;

#ifndef DFSPH_ORIGIN
		aij_pj += mKappa_v[pId] * 2.0f * dt * (1.0f / Alpha[pId] - 1.0f / Alpha_min);
#endif

		Real s_i = -densityAdv[pId];

		Real res = s_i - aij_pj;

#ifdef DFSPH_ORIGIN
		res = res < 0.0f ? res : 0.0f;
		if (nbSize < 20)
			res = 0.0f;
#endif

		Real mKappa_v_i = mKappa_v[pId];

#ifdef DFSPH_ORIGIN
		mKappa_v_i = mKappa_v_i - 0.5f * (s_i - aij_pj) * Alpha[pId] / dt;
		mKappa_v_i = mKappa_v_i > 0.0f ? mKappa_v_i : 0.0f;
#else
		mKappa_v_i = mKappa_v_i - 0.5f * (s_i - aij_pj) * Alpha_min / dt;
#endif
		
		mKappa_v[pId] = mKappa_v_i;
		res_v[pId] = fabs(res);
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_updateVelocity
	(
		DArray<Coord> velArr,
		DArray<Coord> accelArr,
		Real density_0,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= velArr.size()) return;

		velArr[pId] += accelArr[pId] * dt;

	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_PredictDensityAdv_Den(
		DArray<Real> densityAdv,
		DArray<Coord> pos,
		DArray<Real> density,
		DArray<Coord> vel,
		DArrayList<int> neighbors,
		Real density_0,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= densityAdv.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real delta(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				delta += (vel[pId] - vel[j]).dot(g);
			}
		}

		//prevent particle deficiency 
		//semi-anlytical free boundary
		//Real ep = 0.0f;
		Real density_i = density[pId];// > density_0 ? density[pId] : density_0;
		

		density_i = density_i + dt * delta;

#ifndef DFSPH_ORIGIN
		Real ep = 0.0f; //this parameter can control surface tension?
		density_i = density_i > density_0 ? density_i : ep * density_i + (1.0f - ep) * density_0;
#endif // DFSPH_JACOBI

		densityAdv[pId] = density_i;
	}


	template<typename Real>
	__global__ void SADFSPH_DensityKappa(
		DArray<Real> kappa_r,
		DArray<Real> densityAdv,
		DArray<Real> alpha,
		Real density_0,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= kappa_r.size()) return;

		Real res = density_0 - densityAdv[pId];
		//res = res < 0.0f ? res : 0.0f;
		kappa_r[pId] = -res * alpha[pId] / dt / dt / density_0;
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_UpdateDenKappa(
		DArray<Real> mKappa_r,
		DArray<Real> res_r,
		DArray<Coord> pos,
		DArray<Coord> accelArr,
		DArray<Real> densityAdv,
		DArray<Real> Alpha,
		Real Alpha_min,
		DArrayList<int> neighbors,
		Real density_0,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= mKappa_r.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real aij_pj(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				aij_pj += (accelArr[pId] - accelArr[j]).dot(g);
			}
		}

		aij_pj = aij_pj * dt * dt;

#ifndef DFSPH_ORIGIN
		aij_pj += mKappa_r[pId] * 2.0f * dt * dt * (1.0f / Alpha[pId] - 1.0f / Alpha_min);
#endif // !DFSPH_ORIGIN

		Real s_i = (density_0 - densityAdv[pId]) / density_0;

#ifndef DFSPH_ORIGIN
		Real res = s_i - aij_pj;

		Real mKappa_r_i = mKappa_r[pId];
		mKappa_r_i = mKappa_r_i - 0.5f * (s_i - aij_pj) * Alpha_min / dt / dt;

		mKappa_r[pId] = mKappa_r_i;
		res_r[pId] = fabs(res) / density_0;
#else
		Real res = s_i - aij_pj;
		res = res < 0.0f ? res : 0.0f;


		Real mKappa_r_i = mKappa_r[pId];
		mKappa_r_i = mKappa_r_i - 0.5f * (s_i - aij_pj) * Alpha[pId] / dt / dt;
		mKappa_r_i = mKappa_r_i > 0.0f ? mKappa_r_i : 0.0f;

		mKappa_r[pId] = mKappa_r_i;
		res_r[pId] = fabs(res) / density_0;
#endif // !DFSPH_ORIGIN

		
	}

	template<typename Real>
	__global__ void SADFSPH_SourceTerm(
		DArray<Real> source,
		DArray<Real> densityAdv,
		Real density_0
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= source.size()) return;

		source[pId] = (density_0 - densityAdv[pId]) / density_0;
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_LaplacianKappa(
		DArray<Real> Ap,
		DArray<Coord> pos,
		DArray<Coord> accelArr,
		DArrayList<int> neighbors,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Ap.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real aij_pj(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				aij_pj += (accelArr[pId] - accelArr[j]).dot(g);
			}
		}

		aij_pj = aij_pj * dt * dt;

		Ap[pId] = aij_pj;
	}

	template<typename Real>
	__global__ void SADFSPH_LaplacianKappaCorrect(
		DArray<Real> Ap,
		DArray<Real> mKappa_r,
		DArray<Real> Alpha,
		Real Alpha_min,
		Real dt
	)
	{
		// Aii = -2* dt*dt/Alpha_i 
		// clamp Aii to -2* dt*dt/Alpha_min 

		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Ap.size()) return;

		Ap[pId] = Ap[pId] + mKappa_r[pId] * 2.0f * dt * dt *
			(1.0f / Alpha[pId] - 1.0f / Alpha_min);
		//printf("the value:%f\n", Ap[pId]);

	}

	template<typename Real>
	__global__ void  SADFSPH_SourceTerm_Div(
		DArray<Real> source,
		DArray<Real> densityAdv
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= source.size()) return;

		source[pId] = -densityAdv[pId];
	}

	template<typename Real>
	__global__ void  SADFSPH_SourceTerm_Div_Correct(
		DArray<Real> source,
		DArray<Real> density,
		Real density_0,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= source.size()) return;

		if (density[pId] > density_0)
		{
			source[pId] += 100000.0f * (density[pId] - density_0) * dt;
		}
	}

	template<typename Real, typename Coord>
	__global__ void   SADFSPH_LaplacianKappa_Div(
		DArray<Real> Ap,
		DArray<Coord> pos,
		DArray<Coord> accelArr,
		DArrayList<int> neighbors,
		Real mass,
		Real smoothingLength,
		DArray<Real> para,
		Kernel<Real> kernel,
		Real scale,
		Real dt
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Ap.size()) return;

		Coord pos_i = pos[pId];

		List<int>& list_i = neighbors[pId];
		int nbSize = list_i.size();

		Real aij_pj(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r = (pos_i - pos[j]).norm();

			if (r > EPSILON)
			{
				Real s = r * expf(-(para[pId] + para[j]) / 2.0f);
				Coord g = mass * kernel.Gradient(r, smoothingLength, s) * (pos_i - pos[j]) * (1.0f / r) * scale;
				aij_pj += (accelArr[pId] - accelArr[j]).dot(g);
			}
		}

		aij_pj = aij_pj * dt;

		Ap[pId] = aij_pj;
	}

	template<typename Real>
	__global__ void SADFSPH_LaplacianKappaCorrect_Div(
		DArray<Real> Ap,
		DArray<Real> mKappa_v,
		DArray<Real> Alpha,
		Real Alpha_min,
		Real dt
	)
	{
		// Aii = -2* dt/Alpha_i 
		// clamp Aii to -2* dt/Alpha_min 

		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= Ap.size()) return;

		Ap[pId] = Ap[pId] + mKappa_v[pId] * 2.0f * dt *
			(1.0f / Alpha[pId] - 1.0f / Alpha_min);
		//printf("the value:%f\n", Ap[pId]);
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_artifical_viscosity(
		DArray<Coord> pos,
		DArray<Real> para,
		DArrayList<int> neighborIds,
		DArray<Real> density,
		DArray<Coord> velocity,
		DArray<Real> pressure,
		Real pressure_max,
		Real pressure_min,
		Kernel<Real> kernel,
		Real h,
		Real dt,
		Real mass,
		Real factor
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size()) return;

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		Coord viscous_term(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			//J. J. Monaghan 1997 , SPH and Riemann Solvers
			Coord x_ij = pos[pId] - pos[j];
			Real x_ij_v_ij = x_ij.dot(velocity[pId] - velocity[j]);
			Real r_ij = x_ij.norm();
			if (r_ij > EPSILON) {
				Real s = r_ij * expf(-(para[j] + para[pId]) / 2.0f);
				Coord gradW = kernel.Gradient(r_ij, h, s) * x_ij / r_ij;
				Real mu_ij = h * x_ij_v_ij / (r_ij * r_ij + 0.01 * h * h);
				//viscous_term += mass * (-0.005f * h * 1000.0f* mu_ij + 2.0f * mu_ij * mu_ij) / (density[pId] + density[j]) * gradW;
				viscous_term += mass * (-3000.0f * mu_ij + 2.0f * mu_ij * mu_ij) / (density[pId] + density[j]) * gradW;

			}
		}

		Real fac(0.0f);
		Real p_ratio = ((pressure[pId] - pressure_min) / (pressure_max - pressure_min));
		fac = 0.05f + 0.75f * (1.0f - (p_ratio));

		velocity[pId] -= factor * fac * dt * viscous_term;
	}

	//compute color 
	template<typename Real, typename Coord>
	__global__ void SADFSPH_TensioncomputeColor(
		DArray<Real> color,
		DArray<Coord> pos,
		DArray<Real> para,
		DArrayList<int> neighborIds,
		Real density_0,
		DArray<Real> density,
		Kernel<Real> kernel,
		Real h,
		Real dt,
		Real mass,
		Real scale
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size()) return;

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		Real c_i = mass / density[pId] * kernel.Weight(0.0f, h, 0.0f) * scale;

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (pos[pId] - pos[j]).norm();
			if (r_ij > EPSILON) {
				Real s = r_ij * expf(-(para[j] + para[pId]) / 2.0f);
				c_i += mass / density[j] * kernel.Weight(r_ij, h, s) * scale;
			}
		}

		color[pId] = c_i;

		//if (c_i > EPSILON)printf("PId:%d,c_i:%f", pId, c_i);
	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_TensioncomputeGradColor(
		DArray<Real> gradColor2,
		DArray<Real> color,
		DArray<Coord> pos,
		DArray<Real> para,
		DArrayList<int> neighborIds,
		Real density_0,
		DArray<Real> density,
		Kernel<Real> kernel,
		Real h,
		Real dt,
		Real mass,
		Real scale
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size()) return;

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		Coord gradC_i(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (pos[pId] - pos[j]).norm();
			if (r_ij > EPSILON) {
				Real s = r_ij * expf(-(para[j] + para[pId]) / 2.0f);
				Coord g =  kernel.Gradient(r_ij, h, s) * (pos[pId] - pos[j]) * (1.0f / r_ij) * scale;
				gradC_i += mass / density[j] * color[j] * g;
			}
		}

		gradC_i = gradC_i / color[pId];
		gradColor2[pId] = gradC_i.dot(gradC_i);
		//if (gradC_i.dot(gradC_i) > EPSILON) printf("PId:%d,gradc_i2:%f\n", pId, gradC_i.dot(gradC_i));

	}

	template<typename Real, typename Coord>
	__global__ void SADFSPH_TensionForceUpdateVel(
		DArray<Coord> vel,
		DArray<Real> gradColor2,
		DArray<Real> color,
		DArray<Coord> pos,
		DArray<Real> para,
		DArrayList<int> neighborIds,
		Real density_0,
		DArray<Real> density,
		Kernel<Real> kernel,
		Real h,
		Real dt,
		Real mass,
		Real scale,
		Real tensionfactor
	)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size()) return;

		List<int>& list_i = neighborIds[pId];
		int nbSize = list_i.size();

		Real factor = 0.25f * tensionfactor;

		Coord ai(0.0f);

		for (int ne = 0; ne < nbSize; ne++)
		{
			int j = list_i[ne];
			Real r_ij = (pos[pId] - pos[j]).norm();
			if (r_ij > EPSILON) {
				Real s = r_ij * expf(-(para[j] + para[pId]) / 2.0f);
				Coord g = kernel.Gradient(r_ij, h, s) * (pos[pId] - pos[j]) * (1.0f / r_ij) * scale;
				ai += factor * mass / density[j] * (gradColor2[pId] + gradColor2[j]) * g;
			}
		}

		vel[pId] += ai * dt;
		//if (pId==1000)printf("pId:%d,a:%f,%f,%f\n", pId, ai[0], ai[1], ai[2]);
		//if (ai.norm() > EPSILON)printf("pId:%d,a:%f,%f,%f\n", pId, ai[0], ai[1], ai[2]);
	}

	template<typename TDataType>
	void SA_DFSPHsolver<TDataType>::constrain()
	{
		int num = this->inPosition()->size();

		if (this->outDensity()->size() != this->inPosition()->size());
		this->outDensity()->resize(this->inPosition()->size());

		if (mKappa_r.size() != this->inPosition()->size())
			mKappa_r.resize(this->inPosition()->size());
		if (mKappa_v.size() != this->inPosition()->size())
			mKappa_v.resize(this->inPosition()->size());
		if (mAlpha.size() != this->inPosition()->size())
			mAlpha.resize(this->inPosition()->size());
		if (mAlpha_c.size() != this->inPosition()->size())
			mAlpha_c.resize(this->inPosition()->size());
		if (mPressureEstimation.size() != this->inPosition()->size())
			mPressureEstimation.resize(this->inPosition()->size());
		if (mPara.size() != this->inPosition()->size())
			mPara.resize(this->inPosition()->size());
		if (mAccel.size() != this->inPosition()->size())
			mAccel.resize(this->inPosition()->size());
		if (mDensityAdv.size() != this->inPosition()->size())
			mDensityAdv.resize(this->inPosition()->size());
		if (m_resv.size() != this->inPosition()->size())
			m_resv.resize(this->inPosition()->size());
		if (m_resr.size() != this->inPosition()->size())
			m_resr.resize(this->inPosition()->size());
		//if (m_Ap.size() != this->inPosition()->size())
		//	m_Ap.resize(this->inPosition()->size());
		//if (m_source.size() != this->inPosition()->size())
		//	m_source.resize(this->inPosition()->size());
		//if (m_p.size() != this->inPosition()->size())
		//	m_p.resize(this->inPosition()->size());
		//if (m_r.size() != this->inPosition()->size())
		//	m_r.resize(this->inPosition()->size());

		if (m_color.size() != this->inPosition()->size())
			m_color.resize(this->inPosition()->size());
		if (m_gradC2.size() != this->inPosition()->size())
			m_gradC2.resize(this->inPosition()->size());


		mSummation->update();
		int MaxItNum = this->varMaxIterationNumber()->getData();

		Real h = this->varSmoothingLength()->getValue();
		Real mass = mSummation->getParticleMass();
		Real dt = this->inTimeStep()->getValue();
		Real rho_0 = this->varRestDensity()->getData();
		Real scalingfactor = mSummation->getScalingFactor(); //scaling factor for deensity compute

		auto m_reduce_p = Reduction<Real>::Create(num);
		Real pressure_max = m_reduce_p->maximum(mKappa_r.begin(), mKappa_r.size());
		if (fabs(pressure_max) < EPSILON) pressure_max = EPSILON;
		Real pressure_min = m_reduce_p->minimum(mKappa_r.begin(), mKappa_r.size());
		if (-pressure_min > pressure_max) pressure_max = -pressure_min;
		

		//估计核函数参数值
		//cuExecute(num, SADFSPH_PressureEstimation,
		//	mPressureEstimation,
		//	mSummation->outDensity()->getData(),
		//	rho_0
		//);
		//auto m_reduce_p = Reduction<Real>::Create(num);
		//Real pressure_max = m_reduce_p->maximum(mPressureEstimation.begin(), mPressureEstimation.size());
		//if (fabs(pressure_max) < EPSILON) pressure_max = EPSILON;
		//Real pressure_min = m_reduce_p->minimum(mPressureEstimation.begin(), mPressureEstimation.size());
		//if (-pressure_min > pressure_max) pressure_max = -pressure_min;
		//Kappa_r better than Kappa_v, why?
		Real density_max = m_reduce_p->maximum(mSummation->outDensity()->getData().begin(), mSummation->outDensity()->getData().size());
		Real density_min = m_reduce_p->minimum(mSummation->outDensity()->getData().begin(), mSummation->outDensity()->getData().size());
		std::cout << "density range::" << density_min << "-" << density_max << std::endl;
		
		cuExecute(num, SADFSPH_Paramaters_Init,
			mPara,
			mKappa_r,
			pressure_max
		);

		std::cout << "Paramax:" << m_reduce_p->maximum(mPara.begin(), mPara.size()) << "Paramin" << m_reduce_p->minimum(mPara.begin(), mPara.size()) << std::endl;
		
		//DFSPH solver

		//Compute Alpha
		cuExecute(num, SADFSPH_AlphaCompute,
			mAlpha,
			this->inPosition()->getData(),
			mSummation->outDensity()->getData(),
			this->inNeighborIds()->getData(),
			rho_0,
			mass,
			h,
			mPara,
			m_kernel,
			scalingfactor
		);

#ifndef DFSPH_ORIGIN
		//semi-implicit boundary condition
		//alpha_i = min(alpha)
		//if (frag_number < 10)
		{
			Real alpha_min = m_reduce_p->minimum(mAlpha.begin(), mAlpha.size());
			std::cout << "alpha_min:" << alpha_min << std::endl;
			Alpha_min = alpha_min;// Alpha_min > alpha_min ? alpha_min : Alpha_min;
			std::cout << "Alpha_min:" << Alpha_min << std::endl;
		}

		cuExecute(num, SADFSPH_AlphaCorrect,
			mAlpha_c,
			Alpha_min,
			1.0f
		);
#endif
		int it1 = 0;
		int it2 = 0;
		//Divergence Solver	
		if (this->varDivergenceSolverDisabled()->getValue() == false)
		{
			cuExecute(num, SADFSPH_PredictDensityAdv_Div,
				mDensityAdv,
				this->inPosition()->getData(),
				mSummation->outDensity()->getData(),
				this->inVelocity()->getData(),
				this->inNeighborIds()->getData(),
				rho_0,
				mass,
				h,
				mPara,
				m_kernel,
				scalingfactor
			);

			cuExecute(num, SADFSPH_DivergenceKappa,
				mKappa_v,
				mDensityAdv,
				this->inNeighborIds()->getData(),
				mAlpha_c,
				dt
			);
		
			//Jacobi iteration
			Real v_err = 10000.0f;
			while (v_err > this->varDivergenceErrorThreshold()->getValue() && (it1 <= MaxItNum))
			{
				cuExecute(num, SADFSAPH_KappaAccel,
					mAccel,
					this->inPosition()->getData(),
					mKappa_v,
					mSummation->outDensity()->getData(),
					this->inNeighborIds()->getData(),
					rho_0,
					mass,
					h,
					mPara,
					m_kernel,
					scalingfactor
				);

				cuExecute(num, SADFSPH_UpdateDivKappa,
					mKappa_v,
					m_resv,
					this->inPosition()->getData(),
					mAccel,
					mDensityAdv,
					mAlpha,
					Alpha_min,
					this->inNeighborIds()->getData(),
					rho_0,
					mass,
					h,
					mPara,
					m_kernel,
					scalingfactor,
					dt
				);

				auto m_reduce = Reduction<Real>::Create(num);
				Real ave_diverr = m_reduce->average(m_resv.begin(), num);
				v_err = fabs(ave_diverr) * dt / rho_0;
				delete m_reduce;

				std::cout << "divergence error:" << v_err << std::endl;

				it1++;
			}
/*
			////using CG solver
			//{
			//	cuExecute(num, SADFSPH_SourceTerm_Div,
			//		m_source,
			//		mDensityAdv
			//	);

			//	cuExecute(num, SADFSAPH_KappaAccel,
			//		mAccel,
			//		this->inPosition()->getData(),
			//		mKappa_v,
			//		mSummation->outDensity()->getData(),
			//		this->inNeighborIds()->getData(),
			//		rho_0,
			//		mass,
			//		h,
			//		mPara,
			//		m_kernel,
			//		scalingfactor
			//	);

			//	cuExecute(num, SADFSPH_LaplacianKappa_Div,
			//		m_Ap,
			//		this->inPosition()->getData(),
			//		mAccel,
			//		this->inNeighborIds()->getData(),
			//		mass,
			//		h,
			//		mPara,
			//		m_kernel,
			//		scalingfactor,
			//		dt
			//	);

			//	cuExecute(num, SADFSPH_LaplacianKappaCorrect_Div,
			//		m_Ap,
			//		mKappa_v,
			//		mAlpha,
			//		Alpha_min,
			//		dt
			//	);

			//	Function2Pt::subtract(m_r, m_source, m_Ap);
			//	m_p.assign(m_r);

			//	auto m_arithmetic = Arithmetic<Real>::Create(num);

			//	Real rr = m_arithmetic->Dot(m_r, m_r);
			//	Real err = num > 0 ? sqrt(rr / num) : 0.0f;

			//	std::cout << "0:err" << err << std::endl;
			//	Real max_err = err;
			//	if (abs(max_err) < EPSILON) max_err = EPSILON;
			//	Real threshold = this->varDivergenceErrorThreshold()->getValue();

			//	while ((err / max_err > threshold) && (it1 <= MaxItNum) && (err > threshold))
			//	{
			//		it1++;
			//		m_Ap.reset();

			//		cuExecute(num, SADFSAPH_KappaAccel,
			//			mAccel,
			//			this->inPosition()->getData(),
			//			m_p,
			//			mSummation->outDensity()->getData(),
			//			this->inNeighborIds()->getData(),
			//			rho_0,
			//			mass,
			//			h,
			//			mPara,
			//			m_kernel,
			//			scalingfactor
			//		);

			//		cuExecute(num, SADFSPH_LaplacianKappa_Div,
			//			m_Ap,
			//			this->inPosition()->getData(),
			//			mAccel,
			//			this->inNeighborIds()->getData(),
			//			mass,
			//			h,
			//			mPara,
			//			m_kernel,
			//			scalingfactor,
			//			dt
			//		);
			//		//semi-implicit boundary here
			//		cuExecute(num, SADFSPH_LaplacianKappaCorrect_Div,
			//			m_Ap,
			//			m_p,
			//			mAlpha,
			//			Alpha_min,
			//			dt
			//		);

			//		float alpha = rr / (m_arithmetic->Dot(m_p, m_Ap));

			//		//std::cout <<"111" << m_arithmetic->Dot(m_p, m_Ap) << std::endl;
			//		Function2Pt::saxpy(mKappa_v, m_p, mKappa_v, alpha);

			//		Function2Pt::saxpy(m_r, m_Ap, m_r, -alpha);

			//		Real rr_old = rr;
			//		rr = m_arithmetic->Dot(m_r, m_r);

			//		Real beta = rr / rr_old;

			//		Function2Pt::saxpy(m_p, m_p, m_r, beta);
			//		err = sqrt(rr / num);

			//		std::cout << " Divergence iter:" << it1 << "||RelativeError:" << err / max_err * 100 << " % " << std::endl;
			//	}

			//	delete m_arithmetic;
			//}
*/

			//Update Velocity
			cuExecute(num, SADFSAPH_KappaAccel,
				mAccel,
				this->inPosition()->getData(),
				mKappa_v,
				mSummation->outDensity()->getData(),
				this->inNeighborIds()->getData(),
				rho_0,
				mass,
				h,
				mPara,
				m_kernel,
				scalingfactor
			);

			cuExecute(num, SADFSPH_updateVelocity,
				this->inVelocity()->getData(),
				mAccel,
				rho_0,
				dt
			);

		}
		
		//Apply non-pressure force here
		//Adaptive AV here, scaling with Kappa_v
		pressure_max = m_reduce_p->maximum(mKappa_v.begin(), mKappa_v.size());
		if (fabs(pressure_max) < EPSILON) pressure_max = EPSILON;
		pressure_min = m_reduce_p->minimum(mKappa_v.begin(), mKappa_v.size());
		delete m_reduce_p;

		//adaptive viscosity
		cuExecute(num, SADFSPH_artifical_viscosity,
			this->inPosition()->getData(),
			mPara,
			this->inNeighborIds()->getData(),
			mSummation->outDensity()->getData(),
			this->inVelocity()->getData(),
			mKappa_v,
			pressure_max,
			pressure_min,
			m_kernel,
			h,
			dt,
			mass,
			Real(1.0f)
		);

		
		//{
		//	//caculate explicit surface tension here
		//	/*
		//	*References:
		//	*-[HWZ + 14] Xiaowei He, Huamin Wang, Fengjun Zhang, Hongan Wang, Guoping Wang, and Kun Zhou.
		//	* Robust simulation of sparsely sampled thin features in SPH - based free surface flows.ACM Trans.Graph., 34(1) : 7 : 1 - 7 : 9,
		//	* December 2014. URL : http ://doi.acm.org/10.1145/2682630
		//	*/

		//	//compute color feild
		//	cuExecute(num, SADFSPH_TensioncomputeColor,
		//		m_color,
		//		this->inPosition()->getData(),
		//		mPara,
		//		this->inNeighborIds()->getData(),
		//		rho_0,
		//		mSummation->outDensity()->getData(),
		//		m_kernel,
		//		h,
		//		dt,
		//		mass,
		//		scalingfactor
		//	);

		//	cuExecute(num, SADFSPH_TensioncomputeGradColor,
		//		m_gradC2,
		//		m_color,
		//		this->inPosition()->getData(),
		//		mPara,
		//		this->inNeighborIds()->getData(),
		//		rho_0,
		//		mSummation->outDensity()->getData(),
		//		m_kernel,
		//		h,
		//		dt,
		//		mass,
		//		scalingfactor
		//	);
	
		//	cuExecute(num, SADFSPH_TensionForceUpdateVel,
		//		this->inVelocity()->getData(),
		//		m_gradC2,
		//		m_color,
		//		this->inPosition()->getData(),
		//		mPara,
		//		this->inNeighborIds()->getData(),
		//		rho_0,
		//		mSummation->outDensity()->getData(),
		//		m_kernel,
		//		h,
		//		dt,
		//		mass,
		//		scalingfactor,
		//		Real(500.0f)
		//	);	
		//}
		
		//Density Solver
		mAccel.reset();
		if (this->varDensitySolverDisabled()->getValue() == false)
		{
			cuExecute(num, SADFSPH_PredictDensityAdv_Den,
				mDensityAdv,
				this->inPosition()->getData(),
				mSummation->outDensity()->getData(),
				this->inVelocity()->getData(),
				this->inNeighborIds()->getData(),
				rho_0,
				mass,
				h,
				mPara,
				m_kernel,
				scalingfactor,
				dt
			);

			cuExecute(num, SADFSPH_DensityKappa,
				mKappa_r,
				mDensityAdv,
				mAlpha_c,
				rho_0,
				dt
			);

			auto m_reduce2 = Reduction<Real>::Create(num);
			Real densityP_ave = m_reduce2->average(mDensityAdv.begin(), num);
			Real mKappa_r_ave = m_reduce2->average(mKappa_r.begin(), num);
			delete m_reduce2;

			std::cout << "densityP_ave:" << densityP_ave << std::endl;


			Real d_err = 10000.0f;
			while(d_err > this->varDensityErrorThreshold()->getValue() && (it2 < MaxItNum))
			{
				cuExecute(num, SADFSAPH_KappaAccel,
					mAccel,
					this->inPosition()->getData(),
					mKappa_r,
					mSummation->outDensity()->getData(),
					this->inNeighborIds()->getData(),
					rho_0,
					mass,
					h,
					mPara,
					m_kernel,
					scalingfactor
				);

				cuExecute(num, SADFSPH_UpdateDenKappa,
					mKappa_r,
					m_resr,
					this->inPosition()->getData(),
					mAccel,
					mDensityAdv,
					mAlpha,
					Alpha_min,
					this->inNeighborIds()->getData(),
					rho_0,
					mass,
					h,
					mPara,
					m_kernel,
					scalingfactor,
					dt
				);

				auto m_reduce2 = Reduction<Real>::Create(num);
				Real ave_denerr = m_reduce2->average(m_resr.begin(), num);
				d_err = fabs(ave_denerr);
				delete m_reduce2;
				it2++;

				std::cout << "d_err:" << d_err << std::endl;
			}

/*
			////Jacobi迭代不收敛，修正为CG迭代如下：
			//cuExecute(num, SADFSPH_SourceTerm,
			//	m_source,
			//	mDensityAdv,
			//	rho_0
			//);

			//cuExecute(num, SADFSAPH_KappaAccel,
			//	mAccel,
			//	this->inPosition()->getData(),
			//	mKappa_r,
			//	mSummation->outDensity()->getData(),
			//	this->inNeighborIds()->getData(),
			//	rho_0,
			//	mass,
			//	h,
			//	mPara,
			//	m_kernel,
			//	scalingfactor
			//);

			//cuExecute(num, SADFSPH_LaplacianKappa,
			//	m_Ap,
			//	this->inPosition()->getData(),
			//	mAccel,
			//	this->inNeighborIds()->getData(),
			//	mass,
			//	h,
			//	mPara,
			//	m_kernel,
			//	scalingfactor,
			//	dt
			//);
			////semi-implicit boundary
			//cuExecute(num, SADFSPH_LaplacianKappaCorrect,
			//	m_Ap,
			//	mKappa_r,
			//	mAlpha,
			//	Alpha_min,
			//	dt
			//);

			//Function2Pt::subtract(m_r, m_source, m_Ap);

			//m_p.assign(m_r);

			//auto m_arithmetic = Arithmetic<Real>::Create(num);

			//Real rr = m_arithmetic->Dot(m_r, m_r);
			//Real err = num > 0 ? sqrt(rr / num) : 0.0f;

			//Real max_err = err;
			//if (abs(max_err) < EPSILON) max_err = EPSILON;
			//Real threshold = this->varDensityErrorThreshold()->getValue();

			//while ((err / max_err > threshold) && (it2 <= MaxItNum) && (err > threshold))
			//{
			//	it2++;
			//	m_Ap.reset();

			//	cuExecute(num, SADFSAPH_KappaAccel,
			//		mAccel,
			//		this->inPosition()->getData(),
			//		m_p,
			//		mSummation->outDensity()->getData(),
			//		this->inNeighborIds()->getData(),
			//		rho_0,
			//		mass,
			//		h,
			//		mPara,
			//		m_kernel,
			//		scalingfactor
			//	);

			//	cuExecute(num, SADFSPH_LaplacianKappa,
			//		m_Ap,
			//		this->inPosition()->getData(),
			//		mAccel,
			//		this->inNeighborIds()->getData(),
			//		mass,
			//		h,
			//		mPara,
			//		m_kernel,
			//		scalingfactor,
			//		dt
			//	);
			//	//semi-implicit boundary here
			//	cuExecute(num, SADFSPH_LaplacianKappaCorrect,
			//		m_Ap,
			//		m_p,
			//		mAlpha,
			//		Alpha_min,
			//		dt
			//	);

			//	float alpha = rr / (m_arithmetic->Dot(m_p, m_Ap));

			//	//std::cout <<"111" << m_arithmetic->Dot(m_p, m_Ap) << std::endl;
			//	Function2Pt::saxpy(mKappa_r, m_p, mKappa_r, alpha);

			//	Function2Pt::saxpy(m_r, m_Ap, m_r, -alpha);

			//	Real rr_old = rr;
			//	rr = m_arithmetic->Dot(m_r, m_r);

			//	Real beta = rr / rr_old;

			//	Function2Pt::saxpy(m_p, m_p, m_r, beta);
			//	err = sqrt(rr / num);

			//	std::cout << "Density iter:" << it2 << "||RelativeError:" << err / max_err * 100 << " % " << std::endl;
			//}

			//delete m_arithmetic;
*/

			//Update Velocity
			cuExecute(num, SADFSAPH_KappaAccel,
				mAccel,
				this->inPosition()->getData(),
				mKappa_r,
				mSummation->outDensity()->getData(),
				this->inNeighborIds()->getData(),
				rho_0,
				mass,
				h,
				mPara,
				m_kernel,
				scalingfactor
			);

			cuExecute(num, SADFSPH_updateVelocity,
				this->inVelocity()->getData(),
				mAccel,
				rho_0,
				dt
			);

			//std::cout << "*DFSPH::Density Solver::Iteration:" << it2 << "||RelativeError:" << err / max_err * 100 << "%" << std::endl;
		}

		frag_number += 1;
		
		std::cout << "DFSPH with StressAwareKernel || Divergence iters::" << it1 << "||Density iters::" << it2 << std::endl;
	}

	DEFINE_CLASS(SA_DFSPHsolver)
}
