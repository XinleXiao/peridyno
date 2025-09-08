#pragma once
#include "Platform.h"
#include "DeclareEnum.h"

#define ADAPTIVE

namespace dyno {
	template<typename Real>
	class SA_AdaptiveKernel 
	{
	public:
		DYN_FUNC SA_AdaptiveKernel() {};
		DYN_FUNC ~SA_AdaptiveKernel() {};
		/*
		* r:particle distance
		* h:influnce radius
		* r_t:farest negative pressure particles distance
		*/
		DYN_FUNC inline  virtual Real Weight(const Real r, const Real h, const Real r_t)
		{
			return Real(0);
		}

		DYN_FUNC inline virtual Real Gradient(const Real r, const Real h, const Real r_t)
		{
			return Real(0);
		}

		DYN_FUNC inline virtual Real integral(const Real r, const Real h, const Real r_t)
		{
			return Real(0);
		}

	public:
		Real m_scale = Real(1.0f);
	};

	
	//cubic kernel
	template<typename Real>
	class SA_CubicKernel :public SA_AdaptiveKernel<Real>
	{
	public:
		DYN_FUNC SA_CubicKernel() :SA_AdaptiveKernel<Real>() {};
		DYN_FUNC ~SA_CubicKernel() {};
	
		DYN_FUNC inline Real Weight(const Real r, const Real h, const Real r_t) override
		{
			return SA_CubicKernel<Real>::weight(r, h, r_t, this->m_scale);
		}
	
		DYN_FUNC inline Real Gradient(const Real r, const Real h, const Real r_t) override
		{
			return SA_CubicKernel<Real>::gradient(r, h, r_t, this->m_scale);
		}
	
		DYN_FUNC inline Real integral(const Real r, const Real h, const Real r_t) override
		{
			return SA_CubicKernel<Real>::_integral(r, h, r_t, this->m_scale);
		}
	
		DYN_FUNC static inline Real weight(const Real r, const Real h, const Real r_t , Real scale)
		{
			Real q = 2.0f * r / h;
			const Real hh = h * h;
			const Real alpha = 3.0f / (2.0f * (Real)M_PI * hh * h);

			if (q > 2.0f) return 0.0f;
			else if (q >= 1.0f)
			{
				//1/6*(2-q)*(2-q)*(2-q)
				const Real d = 2.0f - q;
				return alpha / 6.0f * d * d * d * scale;
			}
			else
			{
				//(2/3)-q*q+0.5f*q*q*q
				const Real qq = q * q;
				const Real qqq = qq * q;
				return alpha * (2.0f / 3.0f - qq + 0.5f * qqq) * scale;
			}
		}
		
		DYN_FUNC static inline Real gradient(const Real r, const Real h, const Real r_t, Real scale)
		{
			Real q = 2.0f * r / h;
			const Real hh = h * h;
			const Real alpha = 3.0f / (2.0f * (Real)M_PI * hh * h);

			if (q > 2.0f) return 0.0f;
			else if (q >= 1.0f)
			{
				//-0.5*(2.0-q)*(2.0-q)
				const Real d = 2.0f - q;
				return -0.5f * alpha * d * d * scale;
			}
			else
			{
				//-2q+1.5*q*q
				const Real qq = q * q;
				return alpha * (-2.0f * q + 1.5f * qq) * scale;
			}
		}
	
		DYN_FUNC static inline Real _integral(const Real r, const Real h, const Real r_t, Real scale)
		{
			return Real(1);
		}

	};
	
	//adaptive kernel
	template<typename Real>
	class StressAwareKernel :public SA_AdaptiveKernel<Real>
	{
	public:
		DYN_FUNC StressAwareKernel() :SA_AdaptiveKernel<Real>() {};
		DYN_FUNC ~StressAwareKernel() {};

		DYN_FUNC inline Real Weight(const Real r, const Real h, const Real s) override
		{
			return StressAwareKernel<Real>::weight(r, h, s, this->m_scale);
		}

		DYN_FUNC inline Real Gradient(const Real r, const Real h, const Real s) override
		{
			return StressAwareKernel<Real>::gradient(r, h, s, this->m_scale);
		}

		DYN_FUNC inline Real integral(const Real r, const Real h, const Real s) override
		{
			return StressAwareKernel<Real>::_integral(r, h, s, this->m_scale);
		}

		DYN_FUNC static inline Real weight(const Real r, const Real h, const Real r_t, Real scale)
		{

#ifdef ADAPTIVE
			const Real q = r / h;
			const Real qq = q * q;

			//Real s = 0.9f;
			Real s = r_t / h;
			if (s > 0.98f)s = 0.98f;
			if (s < 0.02f)s = 0.02f;
			const Real ss = s * s;
			const Real hh = h * h;
			const Real alpha = 105.0f / (8 * (Real)M_PI * ss * s * (5.0f + 8.0f * s + 9.0f * ss + 8.0f * ss * s + 2.0f * ss * ss) * hh * h) / 16.0f;
			//const Real alpha2 = 3.0f / (4.0 * (Real)M_PI * ss * s * hh * h) / 16.0f;

			//const Real k = 0.0f;
			//const Real alpha = (1 - k) * alpha1 + k * alpha2;


			if (q > 1.0f) return 0.0f;
			else if (q > s)
			{
				Real d = 1.0f - q;
				Real ds = 1.0f - s;
				return alpha * 4.0f * ss * s * d * d * d * (1.0f + q - 2.0f * s) / ds / ds / ds * scale;
			}
			else
			{
				return alpha * (qq * qq - 6.0f * ss * qq + 4.0f * ss * s + ss * ss) * scale;
			}
#else
			Real q = 2.0f * r / h;
			const Real hh = h * h;
			const Real alpha = 3.0f / (2.0f * (Real)M_PI * hh * h);

			if (q > 2.0f) return 0.0f;
			else if (q >= 1.0f)
			{
				//1/6*(2-q)*(2-q)*(2-q)
				const Real d = 2.0f - q;
				return alpha / 6.0f * d * d * d * scale;
			}
			else
			{
				//(2/3)-q*q+0.5f*q*q*q
				const Real qq = q * q;
				const Real qqq = qq * q;
				return alpha * (2.0f / 3.0f - qq + 0.5f * qqq) * scale;
			}
#endif // ADAPTIVE		
		}

		DYN_FUNC static inline Real gradient(const Real r, const Real h, const Real r_t, Real scale)
		{
#ifdef ADAPTIVE
			const Real q = r / h;
			const Real qq = q * q;
			//Real s = 0.9f;
			Real s = r_t / h;
			if (s > 0.98f)s = 0.98f;
			if (s < 0.02f)s = 0.02f;
			const Real ss = s * s;
			const Real hh = h * h;
			const Real alpha = 105.0f / (8 * (Real)M_PI * ss * s * (5.0f + 8.0f * s + 9.0f * ss + 8.0f * ss * s + 2.0f * ss * ss) * hh * h) / 16.0f;
			//const Real alpha2 = 3.0f / (4.0 * (Real)M_PI * ss * s * hh * h) / 16.0f;

			//const Real k = 0.0f;
			//const Real alpha = (1 - k) * alpha1 + k * alpha2;
			

			if (q > 1.0f) return 0.0f;
			else if (q > s)
			{
				Real d = 1.0f - q;
				Real ds = 1.0f - s;
				return alpha * (d * d * 8.0f * ss * s * (3.0f * s - 1.0f - 2.0f * q) / ds / ds / ds) * scale;
			}
			else
			{
				return alpha * (4.0f * qq * q - 12.0f * ss * q) * scale;
			}
#else 
			Real q = 2.0f * r / h;
			const Real hh = h * h;
			const Real alpha = 3.0f / (2.0f * (Real)M_PI * hh * h);

			if (q > 2.0f) return 0.0f;
			else if (q >= 1.0f)
			{
				//-0.5*(2.0-q)*(2.0-q)
				const Real d = 2.0f - q;
				return -0.5f * alpha * d * d * scale;
			}
			else
			{
				//-2q+1.5*q*q
				const Real qq = q * q;
				return alpha * (-2.0f * q + 1.5f * qq) * scale;
			}
#endif
		}

		DYN_FUNC static inline Real _integral(const Real r, const Real h, const Real r_t, Real scale)
		{
			return 1.0f;
		}
	};
	
	//spiky kernel
	template<typename Real>
	class SA_SpikyKernel : public SA_AdaptiveKernel<Real>
	{
	public:
		DYN_FUNC SA_SpikyKernel() : SA_AdaptiveKernel<Real>() {};
		DYN_FUNC ~SA_SpikyKernel() {};

		DYN_FUNC inline Real Weight(const Real r, const Real h, const Real s) override
		{
			return SA_SpikyKernel<Real>::weight(r, h, s , this->m_scale);
		}

		DYN_FUNC inline Real Gradient(const Real r, const Real h, const Real s) override
		{
			return SA_SpikyKernel<Real>::gradient(r, h, s , this->m_scale);
		}

		DYN_FUNC inline Real integral(const Real r, const Real h, const Real s) override
		{
			return SA_SpikyKernel<Real>::_integral(r, h, s , this->m_scale);
		}

		DYN_FUNC static inline Real weight(const Real r, const Real h, const Real r_t, Real scale)
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0f;
			else {
				const Real d = Real(1) - q;
				const Real hh = h * h;
				return 15.0f / ((Real)M_PI * hh * h) * d * d * d * scale;
			}
		}

		DYN_FUNC static inline Real gradient(const Real r, const Real h, const Real r_t, Real scale)
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0;
			//else if (r==0.0f) return 0.0f;
			else {
				const Real d = Real(1) - q;
				const Real hh = h * h;
				return -45.0f / ((Real)M_PI * hh * hh) * d * d * scale;
			}
		}

		DYN_FUNC static inline Real _integral(const Real r, const Real h, const Real r_t, Real scale)
		{
			const Real q = r / h;
			if (q > 1.0f) return 0.0f;
			else {
				const Real qq = q * q;
				const Real hh = h * h;
				return -15.0f / ((Real)M_PI * hh) * (q - Real(1.5) * qq + q * qq - Real(0.25) * qq * qq - Real(0.25)) * scale;
			}
		}


	};

	//cohesion kernel [Akinci et.al 2013, Versatile Surface Tension and Adhesion for SPH Fluids]

	/* template<typename Real>
	 class SA_CohesionKernel : public SA_AdaptiveKernel<Real>
	 {
		 DYN_FUNC SA_CohesionKernel() : SA_AdaptiveKernel<Real>() {};
		 DYN_FUNC ~SA_CohesionKernel() {};

	 };*/
}
