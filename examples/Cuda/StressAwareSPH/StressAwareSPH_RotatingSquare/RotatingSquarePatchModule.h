#pragma once
#include "Module/ComputeModule.h"
#include "Collision/Attribute.h"

namespace dyno {
	template<typename TDataType>
	class RotatingSquarePatchModule : public ComputeModule
	{
		DECLARE_TCLASS(RotatingSquarePatchModule, TDataType)

	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		RotatingSquarePatchModule();
		~RotatingSquarePatchModule() override {};

	public:


		DEF_VAR_IN(Real, TimeStep, "Time step size");
		/**
		* @brief Position
		* Particle position
		*/
		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "Particle position");

		/**
		* @brief Velocity
		* Particle velocity
		*/
		DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "Particle velocity");

		/**
		* @brief Attribute
		* Particle attribute
		*/
		DEF_ARRAY_IN(Attribute, Attribute, DeviceType::GPU, "Particle attribute");

		DEF_VAR_IN(uint, FrameNumber, "Frame number");

		DEF_VAR(Real, InitialAngularVelocity, 0.667f, "");

	protected:
		void compute() override;

	private:
		void begin();
		void end();

		bool integrate();

		bool updateVelocity();
		bool updatePosition();

		unsigned int fragNum = 0;
	private:
		DArray<Coord> m_prePosition;
		DArray<Coord> m_preVelocity;
	};

	IMPLEMENT_TCLASS(RotatingSquarePatchModule, TDataType)
}