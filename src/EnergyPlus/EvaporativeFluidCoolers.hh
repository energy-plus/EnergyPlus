// EnergyPlus, Copyright (c) 1996-2016, The Board of Trustees of the University of Illinois and
// The Regents of the University of California, through Lawrence Berkeley National Laboratory
// (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights
// reserved.
//
// If you have questions about your rights to use or distribute this software, please contact
// Berkeley Lab's Innovation & Partnerships Office at IPO@lbl.gov.
//
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the
// U.S. Government consequently retains certain rights. As such, the U.S. Government has been
// granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable,
// worldwide license in the Software to reproduce, distribute copies to the public, prepare
// derivative works, and perform publicly and display publicly, and to permit others to do so.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// (1) Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory,
//     the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be
//     used to endorse or promote products derived from this software without specific prior
//     written permission.
//
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form
//     without changes from the version obtained under this License, or (ii) Licensee makes a
//     reference solely to the software portion of its product, Licensee must refer to the
//     software as "EnergyPlus version X" software, where "X" is the version number Licensee
//     obtained under this License and may not use a different name for the software. Except as
//     specifically required in this Section (4), Licensee shall not use in a company name, a
//     product name, in advertising, publicity, or other promotional activities any name, trade
//     name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly
//     similar designation, without Lawrence Berkeley National Laboratory's prior written consent.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
// AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the
// features, functionality or performance of the source code ("Enhancements") to anyone; however,
// if you choose to make your Enhancements available either publicly, or directly to Lawrence
// Berkeley National Laboratory, without imposing a separate written license agreement for such
// Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free
// perpetual license to install, use, modify, prepare derivative works, incorporate into other
// computer software, distribute, and sublicense such enhancements or derivative works thereof,
// in binary and source code form.

#ifndef EvaporativeFluidCoolers_hh_INCLUDED
#define EvaporativeFluidCoolers_hh_INCLUDED

#include <memory>
#include <vector>

// EnergyPlus Headers
#include <EnergyPlus.hh>
#include <PlantComponent.hh>
#include <PlantLocation.hh>

namespace EnergyPlus {

	class EvaporativeFluidCooler : public PlantComponent
	{

	public:

		virtual ~EvaporativeFluidCooler() {};

		static PlantComponent * factory( int objectType, std::string objectName );

		virtual void simulate( const PlantLocation & calledFromLocation, bool const FirstHVACIteration, Real64 & CurLoad ) override;

		virtual void getDesignCapacities( const PlantLocation & calledFromLocation, Real64 & MaxLoad, Real64 & MinLoad, Real64 & OptLoad ) override;

		virtual void onInitLoopEquip( const PlantLocation & calledFromLocation ) override;

		virtual void getSizingFactor( Real64 & SizFac ) override;

		static void clear_state();

		enum class PIM : int {
			None = 0,
			StandardDesignCapacity,
			UFactor,
			UserSpecifiedDesignCapacity
		};

		enum class EvaporativeFluidCoolerType : int {
			None = 0,
			SingleSpeed,
			TwoSpeed
		};

		enum class EvaporativeLossBy : int {
			None = 0,
			UserFactor = 80,
			MoistTheory = 81
		};

		enum class BlowdownBy : int {
			None = 0,
			Concentration = 90,
			Schedule = 91
		};

	private:

		// Default Constructor
		EvaporativeFluidCooler() = default;

		// Copy Constructor
		EvaporativeFluidCooler( EvaporativeFluidCooler const & ) = default;

		// Move Constructor
		#if !defined(_MSC_VER) || defined(__INTEL_COMPILER) || (_MSC_VER>=1900)
			EvaporativeFluidCooler( EvaporativeFluidCooler && ) = default;
		#endif

		// Copy Assignment
		EvaporativeFluidCooler &
		operator =( EvaporativeFluidCooler const & ) = default;

		// Move Assignment
		#if !defined(_MSC_VER) || defined(__INTEL_COMPILER) || (_MSC_VER>=1900)
			EvaporativeFluidCooler &
			operator =( EvaporativeFluidCooler && ) = default;
		#endif

		static void getInput();

		void init();

		void size();

		void update();

		void report( bool const RunFlag );

		void calculateWaterUseage();

		void calcSingleSpeedEvapFluidCooler();

		void calcTwoSpeedEvapFluidCooler();

		void simSimpleEvapFluidCooler(
				Real64 const WaterMassFlowRate,
				Real64 const AirFlowRate,
				Real64 const UAdesign,
				Real64 & OutletWaterTemp
		);

		Real64 simpleEvapFluidCoolerUAResidual(
				Real64 const UA, // UA of fluid cooler
				Array1< Real64 > const & Par // par(1) = design fluid cooler load [W]
		);

		// Members
		static std::vector< std::unique_ptr< EvaporativeFluidCooler > > instances;
		static bool getInputFlag;
		static std::string const cEvapFluidCooler_SingleSpeed;
		static std::string const cEvapFluidCooler_TwoSpeed;

	public:
		std::string Name; // User identifier
	private:
		PlantLocation location; // connection location structure
		std::string EvapFluidCoolerType; // Type of evaporative fluid cooler
		EvaporativeFluidCoolerType EvapFluidCoolerType_Num = EvaporativeFluidCoolerType::None;
		PIM PerformanceInputMethod_Num = PIM::None;
		int PlantType_Num = 0;
		Real64 DesignWaterFlowRate = 0.0; // Design water flow rate through the evaporative fluid cooler [m3/s]
		bool DesignWaterFlowRateWasAutoSized = false; // true if design water rate was autosize on input
		Real64 DesignSprayWaterFlowRate = 0.0; // Design spray water flow rate through the evaporative fluid cooler [m3/s]
		Real64 DesWaterMassFlowRate = 0.0; // Design water flow rate through the evaporative fluid cooler [kg/s]
		Real64 HighSpeedAirFlowRate = 0.0; // Air flow rate through evaporative fluid cooler at high speed [m3/s]
		bool HighSpeedAirFlowRateWasAutoSized = false; //true if high speed air rate was autosized
		Real64 HighSpeedFanPower = 0.0; // Fan power at high fan speed [W]
		bool HighSpeedFanPowerWasAutoSized = false; // true if high fan power was autosize on input
		Real64 HighSpeedEvapFluidCoolerUA = 0.0; // UA of evaporative fluid cooler at high fan speed [W/C]
		bool HighSpeedEvapFluidCoolerUAWasAutoSized = false; // true if high speed UA was autosized on input
		Real64 LowSpeedAirFlowRate = 0.0; // Air flow rate through evaporative fluid cooler at low speed [m3/s]
		bool LowSpeedAirFlowRateWasAutoSized = false; // true if low speed air rate was autosize on input
		Real64 LowSpeedAirFlowRateSizingFactor = 0.0; // sizing factor for low speed air flow rate []
		Real64 LowSpeedFanPower = 0.0; // Fan power at low fan speed [W]
		bool LowSpeedFanPowerWasAutoSized = false; // true if low speed fan power set to autosize on input
		Real64 LowSpeedFanPowerSizingFactor = 0.0; // Sizing factor for low speed fan power []
		Real64 LowSpeedEvapFluidCoolerUA = 0.0; // UA of evaporative fluid cooler at low fan speed [W/C]
		bool LowSpeedEvapFluidCoolerUAWasAutoSized = false; //true if low speed UA set to autosize on input
		Real64 LowSpeedEvapFluidCoolerUASizingFactor = 0.0; // sizing factor for low speed UA []
		Real64 DesignEnteringWaterTemp = 0.0; // Entering water temperature at design conditions
		Real64 DesignEnteringAirTemp = 0.0; // Design inlet air dry-bulb temperature (C)
		Real64 DesignEnteringAirWetBulbTemp = 0.0; // Design inlet air wet-bulb temperature (C)
		Real64 EvapFluidCoolerMassFlowRateMultiplier = 0.0; // Maximum evaporative fluid cooler flow rate is
		// this multiplier times design flow rate
		Real64 HeatRejectCapNomCapSizingRatio = 0.0; // ratio of actual cap to nominal capacity []
		Real64 HighSpeedStandardDesignCapacity = 0.0; // Standard Design Capacity of the evaporative fluid cooler [W]
		// with entering water at 35C (95F),
		//  leaving water at 29.44C (85F), entering air at 25.56C (78F) wet-bulb
		//  temp and 35C (95F) dry-bulb temp, and water flow
		//  rate of 5.382E-8 m3/s per watt (3 gpm/ton)
		Real64 LowSpeedStandardDesignCapacity = 0.0; // Standard Design Capacity of the evaporative fluid cooler [W]
		// with entering water at 35C (95F),
		//  leaving water at 29.44C (85F), entering air at 25.56C (78F) wet-bulb
		//  temp and 35C (95F) dry-bulb temp, and water flow
		//  rate of 5.382E-8 m3/s per watt (3 gpm/ton)
		Real64 LowSpeedStandardDesignCapacitySizingFactor = 0.0; // sizing factor for low speed capacity []
		Real64 HighSpeedUserSpecifiedDesignCapacity = 0.0; // User specified design capacity [W]
		Real64 LowSpeedUserSpecifiedDesignCapacity = 0.0; // User specified design capacity for at low speed for
		// two speed fluid cooler[W]
		Real64 LowSpeedUserSpecifiedDesignCapacitySizingFactor = 0.0; // sizing factor for low speed user capacity []
		int FluidIndex = 0; // Index to Property arrays
		Real64 SizFac = 0.0; // sizing factor
		int WaterInletNodeNum = 0; // Node number on the water inlet side of the evaporative fluid cooler
		int WaterOutletNodeNum = 0; // Node number on the water outlet side of the evaporative fluid cooler
		int OutdoorAirInletNodeNum = 0; // Node number of outdoor air inlet for the evaporative fluid cooler
		int HighMassFlowErrorCount = 0; // Counter when mass flow rate is >
		// Design*EvapFluidCoolerMassFlowRateMultiplier
		int HighMassFlowErrorIndex = 0; // Index for high mass flow recurring error message
		int OutletWaterTempErrorCount = 0; // Counter when outlet water temperature is < minimum allowed temperature
		int OutletWaterTempErrorIndex = 0; // Index for outlet water temperature recurring error message
		int SmallWaterMassFlowErrorCount = 0; // Counter when water mass flow rate is very small
		int SmallWaterMassFlowErrorIndex = 0; // Index for very small water mass flow rate recurring error message
		//fluid bypass
		int CapacityControl = 0; // Type of capacity control for single speed cooling tower:
		//  0 - FanCycling, 1 - FluidBypass
		Real64 BypassFraction = 0.0; // Fraction of fluid bypass as a ratio of total fluid flow
		//  through the tower sump
		//begin water system interactions
		EvaporativeLossBy EvapLossMode = EvaporativeLossBy::MoistTheory; // sets how evaporative fluid cooler water evaporation is modeled
		BlowdownBy BlowdownMode = BlowdownBy::Concentration; // sets how evaporative fluid cooler water blowdown is modeled
		int SchedIDBlowdown = 0; // index "pointer" to schedule of blowdown in [m3/s]
		int WaterTankID = 0; // index "pointer" to WaterStorage structure
		int WaterTankDemandARRID = 0; // index "pointer" to demand array inside WaterStorage structure
		Real64 UserEvapLossFactor = 0.0; // simple model [%/Delt C]
		Real64 DriftLossFraction = 0.0;
		Real64 ConcentrationRatio = 0.0; // ratio of solids in blowdown vs make up water
		bool SuppliedByWaterSystem = false;
		//end water system variables

		Real64 WaterUsage = 0; // Evaporative fluid cooler water usage (m3/s)
		bool envrnFlag = true;
		bool oneTimeFlag = true;

		/////// EvaporativeFluidCoolerInletConds ///////
		Real64 WaterTemp = 0.0; // Evaporative fluid cooler water inlet temperature (C)
		Real64 AirTemp = 0.0; // Evaporative fluid cooler air inlet dry-bulb temperature (C)
		Real64 AirWetBulb = 0.0; // Evaporative fluid cooler air inlet wet-bulb temperature (C)
		Real64 AirPress = 0.0; // Evaporative fluid cooler air barometric pressure
		Real64 AirHumRat = 0.0; // Evaporative fluid cooler air inlet humidity ratio (kg/kg)
		////////////////////////////////////

		/////// ReportVars ///////
		Real64 InletWaterTemp = 0.0; // Evaporative fluid cooler inlet water temperature (C)
		Real64 OutletWaterTemp = 0.0; // Evaporative fluid cooler outlet water temperature (C)
		Real64 WaterMassFlowRate = 0.0; // Evaporative fluid cooler water mass flow rate (m3/s)
		Real64 Qactual = 0.0; // Evaporative fluid cooler heat rejection rate (W)
		Real64 FanPower = 0.0; // Evaporative fluid cooler fan power (W)
		Real64 FanEnergy = 0.0; // Evaporative fluid cooler fan energy consumption (J)
		Real64 AirFlowRatio = 0.0; // Air flow ratio through variable speed evaporative fluid cooler
		Real64 WaterAmountUsed = 0.0; // Evaporative fluid cooler make up water usage (m3)
		Real64 EvaporationVdot = 0.0;
		Real64 EvaporationVol = 0.0;
		Real64 DriftVdot = 0.0;
		Real64 DriftVol = 0.0;
		Real64 BlowdownVdot = 0.0;
		Real64 BlowdownVol = 0.0;
		Real64 MakeUpVdot = 0.0;
		Real64 MakeUpVol = 0.0;
		Real64 TankSupplyVdot = 0.0;
		Real64 TankSupplyVol = 0.0;
		Real64 StarvedMakeUpVdot = 0.0;
		Real64 StarvedMakeUpVol = 0.0;
		///// duplicate... ////
		// Real64 BypassFraction = 0.0; // Added for fluid bypass
		//////////////////////////

	};

} // EnergyPlus

#endif
