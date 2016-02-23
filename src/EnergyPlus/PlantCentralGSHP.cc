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

// C++ Headers
#include <cmath>
#include <string>

// ObjexxFCL Headers
#include <ObjexxFCL/Array.functions.hh>
#include <ObjexxFCL/Fmath.hh>
#include <ObjexxFCL/gio.hh>

// EnergyPlus Headers
#include <PlantCentralGSHP.hh>
#include <BranchNodeConnections.hh>
#include <CurveManager.hh>
#include <DataBranchAirLoopPlant.hh>
#include <DataEnvironment.hh>
#include <DataHVACGlobals.hh>
#include <DataIPShortCuts.hh>
#include <DataLoopNode.hh>
#include <DataPlant.hh>
#include <DataPrecisionGlobals.hh>
#include <DataSizing.hh>
#include <EMSManager.hh>
#include <FluidProperties.hh>
#include <General.hh>
#include <InputProcessor.hh>
#include <NodeInputManager.hh>
#include <OutputProcessor.hh>
#include <OutputReportPredefined.hh>
#include <PlantComponent.hh>
#include <PlantLocation.hh>
#include <PlantUtilities.hh>
#include <Psychrometrics.hh>
#include <ReportSizingManager.hh>
#include <ScheduleManager.hh>
#include <UtilityRoutines.hh>

namespace EnergyPlus {

// Contents:
// CentralHeatPumpSystem (CGSHP) System
// ChillerHeaterPerformance:Electric:EIR

namespace PlantCentralGSHP {

	// MODULE INFORMATION:
	//       AUTHOR         PNNL
	//       DATE WRITTEN   Feb 2013
	//       MODIFIED       na
	//       RE-ENGINEERED  na
	// PURPOSE OF THIS MODULE:
	// This module simulates the performance of the Central Plant GSHP systems
	// It currently includes one object: ChillerHeaterPerformance:Electric:EIR.
	// The other object available for this central CGSHP system such as HeatPumpPerformance:WaterToWater:EIR
	//      will be impletemented later.
	// METHODOLOGY EMPLOYED:
	//  Once the PlantLoopManager determines that the Central Plant GSHP
	//  is available to meet a loop cooling and heating demands, it calls SimCentralGroundSourceHeatPump
	//  which in turn calls the electric PlantCentralGSHP model. The PlantCentralGSHP model is based on
	//  polynomial fits of chiller/heater or heat pump performance data.
	// REFERENCES:
	// USE STATEMENTS:

	// Using/Aliasing
	using namespace DataPrecisionGlobals;
	using namespace DataLoopNode;
	using Psychrometrics::PsyCpAirFnWTdb;
	using Psychrometrics::PsyRhoAirFnPbTdbW;
	using FluidProperties::GetDensityGlycol;
	using FluidProperties::GetSpecificHeatGlycol;
	using DataPlant::PlantLoop;

	// Data
	// MODULE PARAMETER DEFINITIONS:
	// Chiller type parameters: Only water cooled chiller is allowed
	int const WaterCooled( 2 );
	int const SmartMixing( 1 );
	int const FullyMixed( 2 );
	bool SimulClgDominant( false );
	bool SimulHtgDominant( false );

	// MODULE VARIABLE DECLARATIONS:
	bool GetInputWrapper( true ); // When TRUE, calls subroutine to read input file.
	int NumWrappers( 0 ); // Number of Wrappers specified in input
	int NumChillerHeaters( 0 ); // Number of Chiller/heaters specified in input
	Real64 CondenserFanPower( 0.0 ); // Condenser Fan Power (fan cycles with compressor) [W]
	Real64 nsvChillerCapFT( 0.0 ); // Chiller/heater capacity fraction (evaluated as a function of temperature)
	Real64 nsvChillerEIRFT( 0.0 ); // Chiller/heater electric input ratio (EIR = 1 / COP) as a function of temperature
	Real64 nsvChillerEIRFPLR( 0.0 ); // Chiller/heater EIR as a function of part-load ratio (PLR)
	Real64 nsvChillerPartLoadRatio( 0.0 ); // Chiller/heater part-load ratio (PLR)
	Real64 nsvChillerCyclingRatio( 0.0 ); // Chiller/heater cycling ratio
	Real64 nsvChillerFalseLoadRate( 0.0 ); // Chiller/heater false load over and above the water-side load [W]

	// Type defining the component specifications

	Array1D_bool CheckEquipName;
	Array1D_bool CHCheckEquipName;
	Array1D_bool HPCheckEquipName;

	// SUBROUTINE SPECIFICATIONS FOR MODULE ChillerElectricEIR

	// Object Data
	Array1D< WrapperSpecs > Wrapper;
	Array1D< ChillerHeaterSpecs > ChillerHeater;
	Array1D< CHReportVars > ChillerHeaterReport;
	Array1D< WrapperReportVars > WrapperReport;

	// MODULE SUBROUTINES:

	// Beginning of Chiller/Heater Module Driver Subroutine
	//*************************************************************************

	// Functions

	PlantComponent * 
	WrapperSpecs::factory( int const EP_UNUSED(objectType), std::string objectName ) {

		// Get user input values
		if ( GetInputWrapper ) {
			GetWrapperInput();
			GetInputWrapper = false;
		}
		
		// Now look for this particular component in the list
		for ( auto & WrapperItem : Wrapper ) {
			if ( WrapperItem.Name == objectName ) {
				return &WrapperItem;
			}
		}
		// If we didn't find it, fatal
		ShowFatalError( "WrapperSpecs::factory : Error getting inputs for Wrapper named: " + objectName );
		// Shut up the compiler
		return nullptr;
	}
	
	void 
	WrapperSpecs::onInitLoopEquip( const PlantLocation & calledFromLocation )
	{
		this -> InitWrapper( 0.0, calledFromLocation.loopNum );
		
		if ( calledFromLocation.loopNum == this->CWLoopNum ) { // Chilled water loop
			this->SizeWrapper();
		}

	}
	
	void 
	WrapperSpecs::getDesignCapacities( const PlantLocation & calledFromLocation, Real64 & MaxCap, Real64 & MinCap, Real64 & OptCap )
	{
		MinCap = 0.0;
		MaxCap = 0.0;
		OptCap = 0.0;
		
		if ( calledFromLocation.loopNum == this->CWLoopNum ) { // Chilled water loop
			if ( this->ControlMode == SmartMixing ) { // control mode is SmartMixing
				for ( int NumChillerHeater = 1; NumChillerHeater <= this->ChillerHeaterNums; ++NumChillerHeater ) {
					MaxCap += this->ChillerHeater( NumChillerHeater ).RefCapCooling * this->ChillerHeater( NumChillerHeater ).MaxPartLoadRatCooling;

					OptCap += this->ChillerHeater( NumChillerHeater ).RefCapCooling * this->ChillerHeater( NumChillerHeater ).OptPartLoadRatCooling;

					MinCap += this->ChillerHeater( NumChillerHeater ).RefCapCooling * this->ChillerHeater( NumChillerHeater ).MinPartLoadRatCooling;
				}
			}
		} else if ( calledFromLocation.loopNum == this->HWLoopNum ) { // Hot water loop
			if ( this->ControlMode == SmartMixing ) { // control mode is SmartMixing
				for ( int NumChillerHeater = 1; NumChillerHeater <= this->ChillerHeaterNums; ++NumChillerHeater ) {
					MaxCap += this->ChillerHeater( NumChillerHeater ).RefCapClgHtg * this->ChillerHeater( NumChillerHeater ).MaxPartLoadRatClgHtg;

					OptCap += this->ChillerHeater( NumChillerHeater ).RefCapClgHtg * this->ChillerHeater( NumChillerHeater ).OptPartLoadRatClgHtg;

					MinCap += this->ChillerHeater( NumChillerHeater ).RefCapClgHtg * this->ChillerHeater( NumChillerHeater ).MinPartLoadRatClgHtg;
				}
			} // End of control mode determination
		} // End of loop determination
	}
	
	void 
	WrapperSpecs::getSizingFactor( Real64 & SizFac )
	{
		SizFac = 1.0; // Always equal to one now. The conponent may have its own sizing factor
	}
	
	void
	WrapperSpecs::simulate( 
		const PlantLocation & calledFromLocation, 
		bool const FirstHVACIteration, 
		Real64 & CurLoad )
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Yunzhi Huang, PNNL
		//       DATE WRITTEN   Feb 2013
		//       MODIFIED       Nov. 2013, Daeho Kang, add component sizing table entries
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		//  This subroutine manages the simulation of the wrapper.

		// METHODOLOGY EMPLOYED:
		//  na

		// REFERENCES:
		//  na

		// Using/Aliasing
		using InputProcessor::GetNumObjectsFound;
		using InputProcessor::GetObjectItem;
		using InputProcessor::VerifyName;
		using InputProcessor::SameString;
		using InputProcessor::FindItemInList;
		using namespace DataIPShortCuts;
		using CurveManager::GetCurveIndex;
		using CurveManager::CurveValue;
		using DataPlant::TypeOf_CentralGroundSourceHeatPump;
		using General::TrimSigDigits;
		using General::RoundSigDigits;
		using PlantUtilities::UpdateChillerComponentCondenserSide;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:
		int WrapperNum; // Wrapper number pointer
		int NumChillerHeater; // Chiller heater number pointer
		int LoopSide; // Plant loop side
		int LoopNum; // plant loop index pointer
		Real64 SimulLoadRatio; // Cooling/heating ratio to determine a load domination
		
		// Not really used. Should be deleted later
		int const EquipFlowCtrl = true; // Flow control mode for the equipment
		bool const RunFlag = true; // Simulate chiller when TRUE

		// Initialization
		LoopNum = calledFromLocation.loopNum;
		WrapperNum = InputProcessor::FindItemInList( this->Name, Wrapper );
		
		if ( LoopNum != this->GLHELoopNum ) {

			this->InitWrapper( CurLoad, LoopNum );
			this->CalcWrapperModel( WrapperNum, CurLoad, RunFlag, FirstHVACIteration, EquipFlowCtrl, LoopNum );

		} else if ( LoopNum == this->GLHELoopNum ) {
			LoopSide = this->GLHELoopSideNum;
			UpdateChillerComponentCondenserSide( LoopNum, LoopSide, TypeOf_CentralGroundSourceHeatPump, this->GLHEInletNodeNum, this->GLHEOutletNodeNum, WrapperReport( WrapperNum ).GLHERate, WrapperReport( WrapperNum ).GLHEInletTemp, WrapperReport( WrapperNum ).GLHEOutletTemp, WrapperReport( WrapperNum ).GLHEmdot, FirstHVACIteration );

			// Use the first chiller heater's evaporator capacity ratio to determine dominant load
			SimulClgDominant = false;
			SimulHtgDominant = false;
			if ( this->WrapperCoolingLoad > 0 && this->WrapperHeatingLoad > 0 ) {
				SimulLoadRatio = this->WrapperCoolingLoad / this->WrapperHeatingLoad;
				if ( SimulLoadRatio > this->ChillerHeater( 1 ).ClgHtgToCoolingCapRatio ) {
					SimulClgDominant = true;
					SimulHtgDominant = false;
				} else {
					SimulHtgDominant = true;
					SimulClgDominant = false;
				}
			}
		}
	}
 
	void
	WrapperSpecs::SizeWrapper()
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Yunzhi Huang, PNNL
		//       DATE WRITTEN   Feb 2013
		//       MODIFIED       Nov. 2013, Daeho Kang, add component sizing table entries
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		//  This subroutine is for sizing all the components under each 'CentralHeatPumpSystem' object,
		//  for which capacities and flow rates have not been specified in the input.

		// METHODOLOGY EMPLOYED:
		//  Obtains evaporator flow rate from the plant sizing array. Calculates reference capacity from
		//  the evaporator (or load side) flow rate and the chilled water loop design delta T. The condenser
		//  flow (or sourse side) rate is calculated from the reference capacity, the COP, and the condenser
		//  loop design delta T.

		// REFERENCES:
		//  na

		// Using/Aliasing
		using namespace DataSizing;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using DataPlant::PlantFirstSizesOkayToReport;
		using DataPlant::PlantFinalSizesOkayToReport;
		using PlantUtilities::RegisterPlantCompDesignFlow;
		using ReportSizingManager::ReportSizingOutput;
		using DataHVACGlobals::SmallWaterVolFlow;
		using DataGlobals::InitConvTemp;
		using DataGlobals::DisplayExtraWarnings;
		using namespace OutputReportPredefined;
		using General::RoundSigDigits;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "SizeCGSHPChillerHeater" );

		// INTERFACE BLOCK SPECIFICATIONS:
		//  na

		// DERIVED TYPE DEFINITIONS:
		//  na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int PltSizNum; // Plant Sizing index corresponding to CurLoopNum
		int PltSizCondNum; // Plant Sizing index for condenser loop
		bool ErrorsFound; // If errors detected in input
		std::string equipName;
		Real64 rho;
		Real64 Cp;
		Real64 tmpNomCap; // local nominal capacity cooling power
		Real64 tmpEvapVolFlowRate; // local evaporator design volume flow rate
		Real64 tmpCondVolFlowRate; // local condenser design volume flow rate
		int NumChillerHeater; // Number of Chiller heater pointer
		int CHWInletNodeNum; // Chilled water inlet node index number
		int CHWOutletNodeNum; // Chilled water outlet node index number
		int GLHEInletNodeNum; // Geo-field water inlet node index number
		int GLHEOutletNodeNum; // Geo-field water outlet node index number
		int HWInletNodeNum; // Hot water inlet node index number
		int HWOutletNodeNum; // Hot water outlet node index number
		Real64 TotalEvapVolFlowRate;
		Real64 TotalCondVolFlowRate;
		Real64 TotalHotWaterVolFlowRate;
		Real64 NomCapUser; // Hardsized nominal capacity cooling power for reporting
		Real64 EvapVolFlowRateUser; // Hardsized evaporator design volume flow rate for reporting
		Real64 CondVolFlowRateUser; // Hardsized condenser design volume flow rate for reporting

		// get all the nodes' indices
		CHWInletNodeNum = this->CHWInletNodeNum;
		CHWOutletNodeNum = this->CHWOutletNodeNum;
		GLHEInletNodeNum = this->GLHEInletNodeNum;
		GLHEOutletNodeNum = this->GLHEOutletNodeNum;
		HWInletNodeNum = this->HWInletNodeNum;
		HWOutletNodeNum = this->HWOutletNodeNum;

		// auto-size the chiller heater components
		if ( this->ControlMode == SmartMixing ) {
			for ( NumChillerHeater = 1; NumChillerHeater <= this->ChillerHeaterNums; ++NumChillerHeater ) {
				ErrorsFound = false;

				// find the appropriate Plant Sizing object
				PltSizNum = PlantLoop( this->CWLoopNum ).PlantSizNum;

				//if ( this->ChillerHeater( NumChillerHeater ).CondVolFlowRate == AutoSize ) {
				PltSizCondNum = PlantLoop( this->GLHELoopNum ).PlantSizNum;
				//}

				tmpNomCap = this->ChillerHeater( NumChillerHeater ).RefCapCooling;
				tmpEvapVolFlowRate = this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate;
				tmpCondVolFlowRate = this->ChillerHeater( NumChillerHeater ).CondVolFlowRate;
				NomCapUser = 0.0;
				EvapVolFlowRateUser = 0.0;
				CondVolFlowRateUser = 0.0;

				// auto-size the Evaporator Flow Rate
				if ( PltSizNum > 0 ) {
					if ( PlantSizData( PltSizNum ).DesVolFlowRate >= SmallWaterVolFlow ) {
						tmpEvapVolFlowRate = PlantSizData( PltSizNum ).DesVolFlowRate * this->ChillerHeater( NumChillerHeater ).SizFac;
						this->ChillerHeater( NumChillerHeater ).tmpEvapVolFlowRate = tmpEvapVolFlowRate;
						if ( ! this->ChillerHeater( NumChillerHeater ).EvapVolFlowRateWasAutoSized ) tmpEvapVolFlowRate = this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate;

					} else {
						if ( this->ChillerHeater( NumChillerHeater ).EvapVolFlowRateWasAutoSized ) tmpEvapVolFlowRate = 0.0;
						this->ChillerHeater( NumChillerHeater ).tmpEvapVolFlowRate = tmpEvapVolFlowRate;

					}
					if ( PlantFirstSizesOkayToFinalize ) {
						if ( this->ChillerHeater( NumChillerHeater ).EvapVolFlowRateWasAutoSized ) {
							this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate = tmpEvapVolFlowRate;
							if ( PlantFinalSizesOkayToReport ) {
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Design Size Reference Chilled Water Flow Rate [m3/s]", tmpEvapVolFlowRate );
							}
							if ( PlantFirstSizesOkayToReport ) {
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Initial Design Size Reference Chilled Water Flow Rate [m3/s]", tmpEvapVolFlowRate );
							}
						} else {
							if ( this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate > 0.0 && tmpEvapVolFlowRate > 0.0
									&& PlantFinalSizesOkayToReport) {
								EvapVolFlowRateUser = this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate;
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Design Size Reference Chilled Water Flow Rate [m3/s]", tmpEvapVolFlowRate,
									"User-Specified Reference Chilled Water Flow Rate [m3/s]", EvapVolFlowRateUser );
								tmpEvapVolFlowRate = EvapVolFlowRateUser;
								if ( DisplayExtraWarnings ) {
									if ( ( std::abs( tmpEvapVolFlowRate - EvapVolFlowRateUser ) / EvapVolFlowRateUser ) > AutoVsHardSizingThreshold ) {
										ShowMessage( "SizeChillerHeaterPerformanceElectricEIR: Potential issue with equipment sizing for " + this->ChillerHeater( NumChillerHeater ).Name );
										ShowContinueError( "User-Specified Reference Chilled Water Flow Rate of " + RoundSigDigits( EvapVolFlowRateUser, 5 ) + " [m3/s]" );
										ShowContinueError( "differs from Design Size Reference Chilled Water Flow Rate of " + RoundSigDigits( tmpEvapVolFlowRate, 5 ) + " [m3/s]" );
										ShowContinueError( "This may, or may not, indicate mismatched component sizes." );
										ShowContinueError( "Verify that the value entered is intended and is consistent with other components." );
									}
								}
							}
						}
					}
				} else {
					if ( this->ChillerHeater( NumChillerHeater ).EvapVolFlowRateWasAutoSized ) {
						if ( PlantFirstSizesOkayToFinalize ) {
							ShowSevereError( "Autosizing of CGSHP Chiller Heater evap flow rate requires a loop Sizing:Plant object" );
							ShowContinueError( "Occurs in CGSHP Chiller Heater Performance object=" + this->ChillerHeater( NumChillerHeater ).Name );
							ErrorsFound = true;
						}
					} else {
						if ( this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate > 0.0 && PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
								"User-Specified Reference Chilled Water Flow Rate [m3/s]", this->ChillerHeater( NumChillerHeater ).EvapVolFlowRate );
						}
					}
				}

				// auto-size the Reference Cooling Capacity
				// each individual chiller heater module is sized to be capable of supporting the total load on the wrapper
				if ( PltSizNum > 0 ) {
					if ( PlantSizData( PltSizNum ).DesVolFlowRate >= SmallWaterVolFlow && tmpEvapVolFlowRate > 0.0 ) {
						Cp = GetSpecificHeatGlycol( PlantLoop( this->CWLoopNum ).FluidName, InitConvTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );

						rho = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, InitConvTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );
						tmpNomCap = Cp * rho * PlantSizData( PltSizNum ).DeltaT * tmpEvapVolFlowRate;
						if ( ! this->ChillerHeater( NumChillerHeater ).RefCapCoolingWasAutoSized ) tmpNomCap = this->ChillerHeater( NumChillerHeater ).RefCapCooling;
					} else {
						if ( this->ChillerHeater( NumChillerHeater ).RefCapCoolingWasAutoSized ) tmpNomCap = 0.0;
					}
					if ( PlantFirstSizesOkayToFinalize ) {
						if ( this->ChillerHeater( NumChillerHeater ).RefCapCoolingWasAutoSized ) {
							this->ChillerHeater( NumChillerHeater ).RefCapCooling = tmpNomCap;

							this->ChillerHeater( NumChillerHeater ).RefCapClgHtg = this->ChillerHeater( NumChillerHeater ).RefCapCooling
								* this->ChillerHeater( NumChillerHeater ).ClgHtgToCoolingCapRatio;
							if ( PlantFinalSizesOkayToReport ) {
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
								"Design Size Reference Capacity [W]", tmpNomCap );
							}
							if ( PlantFirstSizesOkayToReport ) {
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
								"Initial Design Size Reference Capacity [W]", tmpNomCap );
							}
						} else {
							if ( this->ChillerHeater( NumChillerHeater ).RefCapCooling > 0.0 && tmpNomCap > 0.0
									&& PlantFinalSizesOkayToReport ) {
								NomCapUser = this->ChillerHeater( NumChillerHeater ).RefCapCooling;
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Design Size Reference Capacity [W]", tmpNomCap,
									"User-Specified Reference Capacity [W]", NomCapUser );
								tmpNomCap = NomCapUser;
								if ( DisplayExtraWarnings ) {
									if ( ( std::abs( tmpNomCap - NomCapUser ) / NomCapUser ) > AutoVsHardSizingThreshold ) {
										ShowMessage( "SizeChillerHeaterPerformanceElectricEIR: Potential issue with equipment sizing for " + this->ChillerHeater( NumChillerHeater ).Name );
										ShowContinueError( "User-Specified Reference Capacity of " + RoundSigDigits( NomCapUser, 2 ) + " [W]" );
										ShowContinueError( "differs from Design Size Reference Capacity of " + RoundSigDigits( tmpNomCap, 2 ) + " [W]" );
										ShowContinueError( "This may, or may not, indicate mismatched component sizes." );
										ShowContinueError( "Verify that the value entered is intended and is consistent with other components." );
									}
								}
							}
						}
					}
				} else {
					if ( this->ChillerHeater( NumChillerHeater ).RefCapCoolingWasAutoSized ) {
						if ( PlantFirstSizesOkayToFinalize ) {
							ShowSevereError( "Size ChillerHeaterPerformance:Electric:EIR=\"" + this->ChillerHeater( NumChillerHeater ).Name + "\", autosize error." );
							ShowContinueError( "Autosizing of CGSHP Chiller Heater reference capacity requires" );
							ShowContinueError( "a cooling loop Sizing:Plant object." );
							ErrorsFound = true;
						}
					} else {
						if ( this->ChillerHeater( NumChillerHeater ).RefCapCooling > 0.0 && PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
								"User-Specified Reference Capacity [W]", this->ChillerHeater( NumChillerHeater ).RefCapCooling );
						}
					}
				}

				// auto-size the condenser volume flow rate
				// each individual chiller heater module is sized to be capable of supporting the total load on the wrapper
				if ( PltSizCondNum > 0 ) {
					if ( PlantSizData( PltSizNum ).DesVolFlowRate >= SmallWaterVolFlow ) {
						rho = GetDensityGlycol( PlantLoop( this->GLHELoopNum ).FluidName, InitConvTemp, PlantLoop( this->GLHELoopNum ).FluidIndex, RoutineName );
						Cp = GetSpecificHeatGlycol( PlantLoop( this->GLHELoopNum ).FluidName, this->ChillerHeater( NumChillerHeater ).TempRefCondIn, PlantLoop( this->GLHELoopNum ).FluidIndex, RoutineName );
						tmpCondVolFlowRate = tmpNomCap * ( 1.0 + ( 1.0 / this->ChillerHeater( NumChillerHeater ).RefCOPCooling ) * this->ChillerHeater( NumChillerHeater ).OpenMotorEff ) / ( PlantSizData( PltSizCondNum ).DeltaT * Cp * rho );
						this->ChillerHeater( NumChillerHeater ).tmpCondVolFlowRate = tmpCondVolFlowRate;
						if ( ! this->ChillerHeater( NumChillerHeater ).CondVolFlowRateWasAutoSized ) tmpCondVolFlowRate = this->ChillerHeater( NumChillerHeater ).CondVolFlowRate;

					} else {
						if ( this->ChillerHeater( NumChillerHeater ).CondVolFlowRateWasAutoSized ) tmpCondVolFlowRate = 0.0;
						this->ChillerHeater( NumChillerHeater ).tmpCondVolFlowRate = tmpCondVolFlowRate;

					}
					if ( PlantFirstSizesOkayToFinalize ) {
						if ( this->ChillerHeater( NumChillerHeater ).CondVolFlowRateWasAutoSized ) {
							this->ChillerHeater( NumChillerHeater ).CondVolFlowRate = tmpCondVolFlowRate;
							if ( PlantFinalSizesOkayToReport ) {
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Design Size Reference Condenser Water Flow Rate [m3/s]", tmpCondVolFlowRate );
							}
							if ( PlantFirstSizesOkayToReport ) {
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Initial Design Size Reference Condenser Water Flow Rate [m3/s]", tmpCondVolFlowRate );
							}
						} else {
							if ( this->ChillerHeater( NumChillerHeater ).CondVolFlowRate > 0.0 && tmpCondVolFlowRate > 0.0
									&& PlantFinalSizesOkayToReport ) {
								CondVolFlowRateUser = this->ChillerHeater( NumChillerHeater ).CondVolFlowRate;
								ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
									"Design Size Reference Condenser Water Flow Rate [m3/s]", tmpCondVolFlowRate,
									"User-Specified Reference Condenser Water Flow Rate [m3/s]", CondVolFlowRateUser );
								if ( DisplayExtraWarnings ) {
									if ( ( std::abs( tmpCondVolFlowRate - CondVolFlowRateUser ) / CondVolFlowRateUser ) > AutoVsHardSizingThreshold ) {
										ShowMessage( "SizeChillerHeaterPerformanceElectricEIR: Potential issue with equipment sizing for " + this->ChillerHeater( NumChillerHeater ).Name );
										ShowContinueError( "User-Specified Reference Condenser Water Flow Rate of " + RoundSigDigits( CondVolFlowRateUser, 5 ) + " [m3/s]" );
										ShowContinueError( "differs from Design Size Reference Condenser Water Flow Rate of " + RoundSigDigits( tmpCondVolFlowRate, 5 ) + " [m3/s]" );
										ShowContinueError( "This may, or may not, indicate mismatched component sizes." );
										ShowContinueError( "Verify that the value entered is intended and is consistent with other components." );
									}
								}
							}
						}
					}
				} else {
					if ( this->ChillerHeater( NumChillerHeater ).CondVolFlowRateWasAutoSized ) {
						if ( PlantFirstSizesOkayToFinalize ) {
							ShowSevereError( "Size ChillerHeaterPerformance:Electric:EIR=\"" + this->ChillerHeater( NumChillerHeater ).Name + "\", autosize error." );
							ShowContinueError( "Autosizing of CGSHP Chiller Heater condenser flow rate requires" );
							ShowContinueError( "a condenser loop Sizing:Plant object." );
							ErrorsFound = true;
						}
					} else {
						if ( this->ChillerHeater( NumChillerHeater ).CondVolFlowRate > 0.0 && PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( "ChillerHeaterPerformance:Electric:EIR", this->ChillerHeater( NumChillerHeater ).Name,
								"User-Specified Reference Condenser Water Flow Rate [m3/s]", this->ChillerHeater( NumChillerHeater ).CondVolFlowRate );
						}
					}
				}

				if ( PlantFinalSizesOkayToReport ) {
					//create predefined report
					equipName = this->ChillerHeater( NumChillerHeater ).Name;
					PreDefTableEntry( pdchMechType, equipName, "ChillerHeaterPerformance:Electric:EIR" );
					PreDefTableEntry( pdchMechNomEff, equipName, this->ChillerHeater( NumChillerHeater ).RefCOPCooling );
					PreDefTableEntry( pdchMechNomCap, equipName, this->ChillerHeater( NumChillerHeater ).RefCapCooling );
				}

				if ( ErrorsFound ) {
					ShowFatalError( "Preceding sizing errors cause program termination" );
				}

			}

			// sum individual volume flows and register wrapper inlets
			TotalEvapVolFlowRate = 0.0;
			TotalCondVolFlowRate = 0.0;
			TotalHotWaterVolFlowRate = 0.0;
			for ( NumChillerHeater = 1; NumChillerHeater <= this->ChillerHeaterNums; ++NumChillerHeater ) {
				TotalEvapVolFlowRate += this->ChillerHeater( NumChillerHeater ).tmpEvapVolFlowRate;
				TotalCondVolFlowRate += this->ChillerHeater( NumChillerHeater ).tmpCondVolFlowRate;
				TotalHotWaterVolFlowRate += this->ChillerHeater( NumChillerHeater ).DesignHotWaterVolFlowRate;
			}

			RegisterPlantCompDesignFlow( this->CHWInletNodeNum, TotalEvapVolFlowRate );
			RegisterPlantCompDesignFlow( this->HWInletNodeNum, TotalHotWaterVolFlowRate );
			// save the reference condenser water volumetric flow rate for use by the condenser water loop sizing algorithms
			RegisterPlantCompDesignFlow( this->GLHEInletNodeNum, TotalCondVolFlowRate );

			return;
		}
	}

	void
	GetWrapperInput()
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Yunzhi Huang and Daeho Kang, PNNL
		//       DATE WRITTEN:    Feb 2013

		// PURPOSE OF THIS SUBROUTINE:
		//  This routine will get the input required by the Wrapper model.

		// METHODOLOGY EMPLOYED:

		// REFERENCES: na

		// Using/Aliasing
		using InputProcessor::GetNumObjectsFound;
		using InputProcessor::GetObjectItem;
		using InputProcessor::VerifyName;
		using InputProcessor::SameString;
		using InputProcessor::FindItemInList;
		using namespace DataIPShortCuts;
		using BranchNodeConnections::TestCompSet;
		using BranchNodeConnections::SetUpCompSets;
		using NodeInputManager::GetOnlySingleNode;
		using CurveManager::GetCurveIndex;
		using DataGlobals::ScheduleAlwaysOn;
		using CurveManager::CurveValue;
		using ScheduleManager::GetScheduleIndex;
		using General::TrimSigDigits;
		using General::RoundSigDigits;

		// Locals
		static int NumChillerHeaters( 0 ); // total number of chiller heater (without identical multiplier)

		// PARAMETERS
		// na

		// LOCAL VARIABLES
		static std::string CompName; // component name
		static bool ErrorsFound( false ); // True when input errors are found
		bool IsNotOK; // Flag to verify name
		bool IsBlank; // Flag for blank name
		static bool AllocatedFlag( false ); // True when arrays are allocated
		static bool CHAllocatedFlag( false ); // True when arrays are allocated
		int NumAlphas; // Number of elements in the alpha array
		int NumNums; // Number of elements in the numeric array
		int IOStat; // IO Status when calling get input subroutine
		int i_CH; // chiller heater index pointer
		static int WrapperNum( 0 ); // wrapper number
		static int NumberOfComp( 0 ); // number of components under each wrapper
		static int Comp( 0 ); // an index number for input all the components
		static int loop( 0 ); // an index number for read in all the parameters of a component
		static int CompIndex( 0 ); // component index in the sequence of internal input array
		static int ChillerHeaterNum( 1 ); // chiller heater index pointer for current wrapper object
		static int NumChHtrPerWrapper( 0 ); // total number of chiller heaters (including identical units) per wrapper

		if ( AllocatedFlag ) return;
		cCurrentModuleObject = "CentralHeatPumpSystem";
		NumWrappers = GetNumObjectsFound( cCurrentModuleObject );

		if ( NumWrappers <= 0 ) {
			ShowSevereError( "No " + cCurrentModuleObject + " equipment specified in input file" );
			ErrorsFound = true;
		}

		// ALLOCATE ARRAYS
		Wrapper.allocate( NumWrappers );
		WrapperReport.allocate( NumWrappers );
		CheckEquipName.dimension( NumWrappers, true );
		AllocatedFlag = true;

		// Load arrays with electric EIR chiller data
		for ( WrapperNum = 1; WrapperNum <= NumWrappers; ++WrapperNum ) {
			GetObjectItem( cCurrentModuleObject, WrapperNum, cAlphaArgs, NumAlphas, rNumericArgs, NumNums, IOStat, _, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

			Wrapper( WrapperNum ).Name = cAlphaArgs( 1 );

			// intialize nth chiller heater index (including identical units) for current wrapper
			NumChHtrPerWrapper = 0;

			IsNotOK = false;
			IsBlank = false;
			VerifyName( cAlphaArgs( 1 ), Wrapper, WrapperNum - 1, IsNotOK, IsBlank, cCurrentModuleObject + " Name" );

			if ( IsNotOK ) {
				ErrorsFound = true;
				if ( IsBlank ) cAlphaArgs( 1 ) = "xxxxx";
				continue;
			}

			if ( cAlphaArgs( 2 ) == "SMARTMIXING" ) {
				Wrapper( WrapperNum ).ControlMode = SmartMixing;
			}

			Wrapper( WrapperNum ).CHWInletNodeNum = GetOnlySingleNode( cAlphaArgs( 3 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Inlet, 1, ObjectIsNotParent ); //node name : connection should be careful!
			Wrapper( WrapperNum ).CHWOutletNodeNum = GetOnlySingleNode( cAlphaArgs( 4 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Outlet, 1, ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, cAlphaArgs( 1 ), cAlphaArgs( 3 ), cAlphaArgs( 4 ), "Chilled Water Nodes" );

			Wrapper( WrapperNum ).GLHEInletNodeNum = GetOnlySingleNode( cAlphaArgs( 5 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Inlet, 2, ObjectIsNotParent ); //node name : connection should be careful!
			Wrapper( WrapperNum ).GLHEOutletNodeNum = GetOnlySingleNode( cAlphaArgs( 6 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Outlet, 2, ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, cAlphaArgs( 1 ), cAlphaArgs( 5 ), cAlphaArgs( 6 ), "GLHE Nodes" );

			Wrapper( WrapperNum ).HWInletNodeNum = GetOnlySingleNode( cAlphaArgs( 7 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Inlet, 3, ObjectIsNotParent ); //node name : connection should be careful!
			Wrapper( WrapperNum ).HWOutletNodeNum = GetOnlySingleNode( cAlphaArgs( 8 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Outlet, 3, ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, cAlphaArgs( 1 ), cAlphaArgs( 7 ), cAlphaArgs( 8 ), "Hot Water Nodes" );

			Wrapper( WrapperNum ).AncilliaryPower = rNumericArgs( 1 );
			Wrapper( WrapperNum ).AncilliaryPwSchedule = cAlphaArgs( 9 );
			if ( lAlphaFieldBlanks( 9 ) ) {
				Wrapper( WrapperNum ).SchedPtr = 0;
			} else {
				Wrapper( WrapperNum ).SchedPtr = GetScheduleIndex( cAlphaArgs( 9 ) );
			}

			NumberOfComp = ( NumAlphas - 9 ) / 3;
			Wrapper( WrapperNum ).NumOfComp = NumberOfComp;
			Wrapper( WrapperNum ).WrapperComp.allocate( NumberOfComp );

			if ( Wrapper( WrapperNum ).NumOfComp == 0 ) {
				ShowSevereError( "GetWrapperInput: No component names on " + cCurrentModuleObject + '=' + Wrapper( WrapperNum ).Name );
				ErrorsFound = true;
			} else {
				Comp = 0;
				for ( loop = 10; loop <= NumAlphas; loop += 3 ) {
					++Comp;
					Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperPerformanceObjectType = cAlphaArgs( loop );
					Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperComponentName = cAlphaArgs( loop + 1 );
					Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperPerformanceObjectSch = cAlphaArgs( loop + 2 );
					if ( lAlphaFieldBlanks( loop + 2 ) ) {
						Wrapper( WrapperNum ).WrapperComp( Comp ).CHSchedPtr = ScheduleAlwaysOn;
					} else {
						Wrapper( WrapperNum ).WrapperComp( Comp ).CHSchedPtr = GetScheduleIndex( cAlphaArgs( loop + 2 ) );
					}
					Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperIdenticalObjectNum = rNumericArgs( 1 + Comp );
					if ( Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperPerformanceObjectType == "CHILLERHEATERPERFORMANCE:ELECTRIC:EIR" ) {

						// count number of chiller heaters (including identical units) for current wrapper
						if ( Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperIdenticalObjectNum > 1 ) {
							NumChHtrPerWrapper += Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperIdenticalObjectNum;
						} else {
							++NumChHtrPerWrapper;
						}

						// count total number of chiller heaters (not including identical units) for ALL wrappers
						++NumChillerHeaters;

					}
				}

				Wrapper( WrapperNum ).ChillerHeaterNums = NumChHtrPerWrapper;
			}

			if ( ErrorsFound ) {
				ShowFatalError( "GetWrapperInput: Invalid " + cCurrentModuleObject + " Input, preceding condition(s) cause termination." );
			}

			// ALLOCATE ARRAYS
			if ( ( NumChillerHeaters == 0 ) && ( Wrapper( WrapperNum ).ControlMode == SmartMixing ) ) {
				ShowFatalError( "SmartMixing Control Mode in object " + cCurrentModuleObject + " : " + Wrapper( WrapperNum ).Name + " need to apply to ChillerHeaterPerformance:Electric:EIR object(s)." );
			}

		}

		if ( NumChillerHeaters > 0 ) {
			if ( CHAllocatedFlag ) return;

			for ( WrapperNum = 1; WrapperNum <= NumWrappers; ++WrapperNum ) {
				Wrapper( WrapperNum ).ChillerHeater.allocate( Wrapper( WrapperNum ).ChillerHeaterNums );
				Wrapper( WrapperNum ).ChillerHeaterReport.allocate( Wrapper( WrapperNum ).ChillerHeaterNums );
			}
			GetChillerHeaterInput();
			CHAllocatedFlag = true;
		}

		for ( WrapperNum = 1; WrapperNum <= NumWrappers; ++WrapperNum ) {
			ChillerHeaterNum = 0; // intialize nth chiller heater index (including identical units) for current wrapper
			for ( Comp = 1; Comp <= Wrapper( WrapperNum ).NumOfComp; ++Comp ) {
				if ( Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperPerformanceObjectType == "CHILLERHEATERPERFORMANCE:ELECTRIC:EIR" ) {
					CompName = Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperComponentName;
					CompIndex = FindItemInList( CompName, ChillerHeater );
					// User may enter invalid name rather than selecting one from the object list
					if ( CompIndex <= 0 ) {
						ShowSevereError( "GetWrapperInput: Invalid Chiller Heater Modules Performance Component Name =" + CompName );
						ShowContinueError( "Select the name of ChillerHeaterPerformance:Electric:EIR object(s) from the object list." );
						ShowFatalError( "Program terminates due to preceding condition." );
					}
					Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperPerformanceObjectIndex = CompIndex;
					if ( ChillerHeater( CompIndex ).VariableFlow ) {
						Wrapper( WrapperNum ).VariableFlowCH = true;
					}
					for ( i_CH = 1; i_CH <= Wrapper( WrapperNum ).WrapperComp( Comp ).WrapperIdenticalObjectNum; ++i_CH ) {
						// increment nth chiller heater index (including identical units) for current wrapper
						++ChillerHeaterNum;
						Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ) = ChillerHeater( CompIndex );
						Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ) = ChillerHeaterReport( CompIndex );
					}
				}
			}
		}

		//Release memory from temporary arrays; values now copied into their associated Wrapper in above loop
		if ( allocated( ChillerHeater ) ) ChillerHeater.deallocate();
		if ( allocated( ChillerHeaterReport ) ) ChillerHeaterReport.deallocate();

		//Set up output variables
		for ( WrapperNum = 1; WrapperNum <= NumWrappers; ++WrapperNum ) {
			SetupOutputVariable( "Chiller Heater System Cooling Electric Energy [J]", WrapperReport( WrapperNum ).TotElecCooling, "System", "Sum", Wrapper( WrapperNum ).Name, _, "ELECTRICITY", "Cooling", _, "Plant" );

			SetupOutputVariable( "Chiller Heater System Heating Electric Energy [J]", WrapperReport( WrapperNum ).TotElecHeating, "System", "Sum", Wrapper( WrapperNum ).Name, _, "ELECTRICITY", "Heating", _, "Plant" );

			SetupOutputVariable( "Chiller Heater System Cooling Electric Power [W]", WrapperReport( WrapperNum ).TotElecCoolingPwr, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Heating Electric Power [W]", WrapperReport( WrapperNum ).TotElecHeatingPwr, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Cooling Energy [J]", WrapperReport( WrapperNum ).CoolingEnergy, "System", "Sum", Wrapper( WrapperNum ).Name, _, "ENERGYTRANSFER", "CHILLERS", _, "Plant" );

			SetupOutputVariable( "Chiller Heater System Heating Energy [J]", WrapperReport( WrapperNum ).HeatingEnergy, "System", "Sum", Wrapper( WrapperNum ).Name, _, "ENERGYTRANSFER", "BOILER", _, "Plant" );

			SetupOutputVariable( "Chiller Heater System Source Heat Transfer Energy [J]", WrapperReport( WrapperNum ).GLHEEnergy, "System", "Sum", Wrapper( WrapperNum ).Name, _, "ENERGYTRANSFER", "HEATREJECTION", _, "Plant" );

			SetupOutputVariable( "Chiller Heater System Cooling Rate [W]", WrapperReport( WrapperNum ).CoolingRate, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Heating Rate [W]", WrapperReport( WrapperNum ).HeatingRate, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Source Heat Transfer Rate [W]", WrapperReport( WrapperNum ).GLHERate, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Cooling Mass Flow Rate [kg/s]", WrapperReport( WrapperNum ).CHWmdot, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Heating Mass Flow Rate [kg/s]", WrapperReport( WrapperNum ).HWmdot, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Source Mass Flow Rate [kg/s]", WrapperReport( WrapperNum ).GLHEmdot, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Cooling Inlet Temperature [C]", WrapperReport( WrapperNum ).CHWInletTemp, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Heating Inlet Temperature [C]", WrapperReport( WrapperNum ).HWInletTemp, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Source Inlet Temperature [C]", WrapperReport( WrapperNum ).GLHEInletTemp, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Cooling Outlet Temperature [C]", WrapperReport( WrapperNum ).CHWOutletTemp, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Heating Outlet Temperature [C]", WrapperReport( WrapperNum ).HWOutletTemp, "System", "Average", Wrapper( WrapperNum ).Name );

			SetupOutputVariable( "Chiller Heater System Source Outlet Temperature [C]", WrapperReport( WrapperNum ).GLHEOutletTemp, "System", "Average", Wrapper( WrapperNum ).Name );

			if ( Wrapper( WrapperNum ).ChillerHeaterNums > 0 ) {

				for ( ChillerHeaterNum = 1; ChillerHeaterNum <= Wrapper( WrapperNum ).ChillerHeaterNums; ++ChillerHeaterNum ) {

					SetupOutputVariable( "Chiller Heater Operation Mode Unit " + TrimSigDigits( ChillerHeaterNum ) + " []", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CurrentMode, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Part Load Ratio Unit " + TrimSigDigits( ChillerHeaterNum ) + " []", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatio, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Cycling Ratio Unit " + TrimSigDigits( ChillerHeaterNum ) + " []", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatio, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Cooling Electric Power Unit " + TrimSigDigits( ChillerHeaterNum ) + " [W]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingPower, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Heating Electric Power Unit " + TrimSigDigits( ChillerHeaterNum ) + " [W]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).HeatingPower, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Cooling Electric Energy Unit " + TrimSigDigits( ChillerHeaterNum ) + " [J]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy, "System", "Sum", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Heating Electric Energy Unit " + TrimSigDigits( ChillerHeaterNum ) + " [J]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy, "System", "Sum", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Cooling Rate Unit " + TrimSigDigits( ChillerHeaterNum ) + " [W]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).QEvap, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Cooling Energy Unit " + TrimSigDigits( ChillerHeaterNum ) + " [J]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy, "System", "Sum", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater False Load Heat Transfer Rate Unit " + TrimSigDigits( ChillerHeaterNum ) + " [W]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater False Load Heat Transfer Energy Unit " + TrimSigDigits( ChillerHeaterNum ) + " [J]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoad, "System", "Sum", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Evaporator Inlet Temperature Unit " + TrimSigDigits( ChillerHeaterNum ) + " [C]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapInletTemp, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Evaporator Outlet Temperature Unit " + TrimSigDigits( ChillerHeaterNum ) + " [C]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Evaporator Mass Flow Rate Unit " + TrimSigDigits( ChillerHeaterNum ) + " [kg/s]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).Evapmdot, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Condenser Heat Transfer Rate Unit " + TrimSigDigits( ChillerHeaterNum ) + " [W]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).QCond, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Condenser Heat Transfer Energy Unit " + TrimSigDigits( ChillerHeaterNum ) + " [J]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondEnergy, "System", "Sum", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater COP Unit " + TrimSigDigits( ChillerHeaterNum ) + " [W/W]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ActualCOP, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Capacity Temperature Modifier Multiplier Unit " + TrimSigDigits( ChillerHeaterNum ) + " []", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFT, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater EIR Temperature Modifier Multiplier Unit " + TrimSigDigits( ChillerHeaterNum ) + " []", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFT, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater EIR Part Load Modifier Multiplier Unit " + TrimSigDigits( ChillerHeaterNum ) + " []", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLR, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Condenser Inlet Temperature Unit " + TrimSigDigits( ChillerHeaterNum ) + " [C]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondInletTemp, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Condenser Outlet Temperature Unit " + TrimSigDigits( ChillerHeaterNum ) + " [C]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );

					SetupOutputVariable( "Chiller Heater Condenser Mass Flow Rate Unit " + TrimSigDigits( ChillerHeaterNum ) + " [kg/s]", Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).Condmdot, "System", "Average", Wrapper( WrapperNum ).ChillerHeater( ChillerHeaterNum ).Name );
				} // End of individual chiller heater count for current wrapper

			} // End of individual chiller heater output

		} // End of wrapper count

	}

	void
	GetChillerHeaterInput()
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Kyung Tae Yun, Mississippi State University
		//       DATE WRITTEN:    Feb 2013

		// PURPOSE OF THIS SUBROUTINE:
		//  This routine will get the input required by the ChillerHeaterPerformance:Electric:EIR model.

		// METHODOLOGY EMPLOYED:

		// REFERENCES: na

		// Using/Aliasing
		using InputProcessor::GetNumObjectsFound;
		using InputProcessor::GetObjectItem;
		using InputProcessor::VerifyName;
		using InputProcessor::SameString;
		using namespace DataIPShortCuts;
		using BranchNodeConnections::TestCompSet;
		using NodeInputManager::GetOnlySingleNode;
		using CurveManager::GetCurveIndex;
		using CurveManager::GetCurveMinMaxValues;
		using CurveManager::CurveValue;
		using PlantUtilities::RegisterPlantCompDesignFlow;
		using ScheduleManager::GetScheduleIndex;
		using General::TrimSigDigits;
		using General::RoundSigDigits;
		using DataSizing::AutoSize;

		// Locals
		// PARAMETERS
		// na

		// LOCAL VARIABLES
		std::string StringVar; // Used for EIRFPLR warning messages
		static bool CHErrorsFound( false ); // True when input errors are found
		bool IsNotOK; // Flag to verify name
		bool IsBlank; // Flag for blank name
		static bool FoundNegValue( false ); // Used to evaluate PLFFPLR curve objects
		int CurveValPtr; // Index to EIRFPLR curve output
		static int CurveCheck( 0 ); // Used to evaluate PLFFPLR curve objects
		int ChillerHeaterNum; // Chiller counter
		int NumAlphas; // Number of elements in the alpha array
		int NumNums; // Number of elements in the numeric array
		int IOStat; // IO Status when calling get input subroutine
		Real64 CurveVal; // Used to verify EIR-FT and CAP-FT curves
		Array1D< Real64 > CurveValArray( 11 ); // Used to evaluate PLFFPLR curve objects
		Real64 CurveValTmp; // Used to evaluate PLFFPLR curve objects

		// Formats
		static gio::Fmt Format_530( "('Curve Output = ',11(F7.2))" );
		static gio::Fmt Format_550( "('Curve Output = ',11(F7.2))" );

		cCurrentModuleObject = "ChillerHeaterPerformance:Electric:EIR";
		NumChillerHeaters = GetNumObjectsFound( cCurrentModuleObject );

		if ( NumChillerHeaters <= 0 ) {
			ShowSevereError( "No " + cCurrentModuleObject + " equipment specified in input file" );
			CHErrorsFound = true;
		}

		// Allocate temporary ChillerHeater and ChillerHeaterReport arrays
		if ( allocated( ChillerHeater ) ) ChillerHeater.deallocate();
		if ( allocated( ChillerHeaterReport ) ) ChillerHeaterReport.deallocate();
		ChillerHeater.allocate( NumChillerHeaters );
		ChillerHeaterReport.allocate( NumChillerHeaters );

		// Load arrays with electric EIR chiller data
		for ( ChillerHeaterNum = 1; ChillerHeaterNum <= NumChillerHeaters; ++ChillerHeaterNum ) {
			GetObjectItem( cCurrentModuleObject, ChillerHeaterNum, cAlphaArgs, NumAlphas, rNumericArgs, NumNums, IOStat, _, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

			ChillerHeater( ChillerHeaterNum ).Name = cAlphaArgs( 1 );

			IsNotOK = false;
			IsBlank = false;
			VerifyName( cAlphaArgs( 1 ), ChillerHeater, ChillerHeaterNum - 1, IsNotOK, IsBlank, cCurrentModuleObject + " Name" );
			if ( IsNotOK ) {
				CHErrorsFound = true;
				if ( IsBlank ) cAlphaArgs( 1 ) = "xxxxx";
			}

			ChillerHeater( ChillerHeaterNum ).CondModeCooling = cAlphaArgs( 4 );

			// Performance curves
			ChillerHeater( ChillerHeaterNum ).ChillerCapFTCooling = GetCurveIndex( cAlphaArgs( 5 ) );
			if ( ChillerHeater( ChillerHeaterNum ).ChillerCapFTCooling == 0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 5 ) + '=' + cAlphaArgs( 5 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).ChillerEIRFTCooling = GetCurveIndex( cAlphaArgs( 6 ) );
			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFTCooling == 0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 6 ) + '=' + cAlphaArgs( 6 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling = GetCurveIndex( cAlphaArgs( 7 ) );
			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling == 0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 7 ) + '=' + cAlphaArgs( 7 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).CondModeHeating = cAlphaArgs( 8 );

			// Performance curves
			ChillerHeater( ChillerHeaterNum ).ChillerCapFTHeating = GetCurveIndex( cAlphaArgs( 9 ) );
			if ( ChillerHeater( ChillerHeaterNum ).ChillerCapFTHeating == 0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 9 ) + '=' + cAlphaArgs( 9 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).ChillerEIRFTHeating = GetCurveIndex( cAlphaArgs( 10 ) );
			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFTHeating == 0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 10 ) + '=' + cAlphaArgs( 10 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating = GetCurveIndex( cAlphaArgs( 11 ) );
			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating == 0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 11 ) + '=' + cAlphaArgs( 11 ) );
				CHErrorsFound = true;
			}

			if ( cAlphaArgs( 2 ) == "CONSTANTFLOW" ) {
				ChillerHeater( ChillerHeaterNum ).ConstantFlow = true;
				ChillerHeater( ChillerHeaterNum ).VariableFlow = false;
			} else if ( cAlphaArgs( 2 ) == "VARIABLEFLOW" ) {
				ChillerHeater( ChillerHeaterNum ).ConstantFlow = false;
				ChillerHeater( ChillerHeaterNum ).VariableFlow = true;
			} else { // Assume a constant flow chiller if none is specified
				ChillerHeater( ChillerHeaterNum ).ConstantFlow = true;
				ChillerHeater( ChillerHeaterNum ).VariableFlow = false;
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 2 ) + '=' + cAlphaArgs( 2 ) );
				ShowContinueError( "simulation assumes CONSTANTFLOW and continues.." );
			}

			if ( ChillerHeaterNum > 1 ) {
				if ( ChillerHeater( ChillerHeaterNum ).ConstantFlow != ChillerHeater( ChillerHeaterNum - 1 ).ConstantFlow ) {
					ChillerHeater( ChillerHeaterNum ).ConstantFlow = true;
					ShowWarningError( "Water flow mode is different from the other chiller heater(s) " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
					ShowContinueError( "Entered in " + cAlphaFieldNames( 2 ) + '=' + cAlphaArgs( 2 ) );
					ShowContinueError( "Simulation assumes CONSTANTFLOW and continues.." );
				}
			}

			if ( SameString( cAlphaArgs( 3 ), "WaterCooled" ) ) {
				ChillerHeater( ChillerHeaterNum ).CondenserType = WaterCooled;
			} else {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cAlphaFieldNames( 3 ) + '=' + cAlphaArgs( 3 ) );
				ShowContinueError( "Valid entries is WaterCooled" );
				CHErrorsFound = true;
			}

			// Chiller rated performance data
			ChillerHeater( ChillerHeaterNum ).RefCapCooling = rNumericArgs( 1 );
			if ( ChillerHeater( ChillerHeaterNum ).RefCapCooling == AutoSize ) {
				ChillerHeater( ChillerHeaterNum ).RefCapCoolingWasAutoSized = true;
			}
			if ( rNumericArgs( 1 ) == 0.0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cNumericFieldNames( 1 ) + '=' + RoundSigDigits( rNumericArgs( 1 ), 2 ) );
				CHErrorsFound = true;
			}
			ChillerHeater( ChillerHeaterNum ).RefCOPCooling = rNumericArgs( 2 );
			if ( rNumericArgs( 2 ) == 0.0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cNumericFieldNames( 2 ) + '=' + RoundSigDigits( rNumericArgs( 2 ), 2 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).TempRefEvapOutCooling = rNumericArgs( 3 );
			ChillerHeater( ChillerHeaterNum ).TempRefCondInCooling = rNumericArgs( 4 );
			ChillerHeater( ChillerHeaterNum ).TempRefCondOutCooling = rNumericArgs( 5 );
			ChillerHeater( ChillerHeaterNum ).ClgHtgToCoolingCapRatio = rNumericArgs( 6 );
			if ( ! ChillerHeater( ChillerHeaterNum ).RefCapCoolingWasAutoSized ) {
				ChillerHeater( ChillerHeaterNum ).RefCapClgHtg = rNumericArgs( 6 ) * ChillerHeater( ChillerHeaterNum ).RefCapCooling;
			}

			if ( rNumericArgs( 6 ) == 0.0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cNumericFieldNames( 6 ) + '=' + RoundSigDigits( rNumericArgs( 6 ), 2 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).ClgHtgtoCogPowerRatio = rNumericArgs( 7 );
			ChillerHeater( ChillerHeaterNum ).RefPowerClgHtg = ChillerHeater( ChillerHeaterNum ).RefCapCooling / ChillerHeater( ChillerHeaterNum ).RefCOPCooling * rNumericArgs( 7 );
			ChillerHeater( ChillerHeaterNum ).RefCOPClgHtg = ChillerHeater( ChillerHeaterNum ).RefCapClgHtg / ChillerHeater( ChillerHeaterNum ).RefPowerClgHtg;

			if ( rNumericArgs( 7 ) == 0.0 ) {
				ShowSevereError( "Invalid " + cCurrentModuleObject + '=' + cAlphaArgs( 1 ) );
				ShowContinueError( "Entered in " + cNumericFieldNames( 7 ) + '=' + RoundSigDigits( rNumericArgs( 7 ), 2 ) );
				CHErrorsFound = true;
			}

			ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg = rNumericArgs( 8 );
			ChillerHeater( ChillerHeaterNum ).TempRefCondOutClgHtg = rNumericArgs( 9 );
			ChillerHeater( ChillerHeaterNum ).TempRefCondInClgHtg = rNumericArgs( 10 );
			ChillerHeater( ChillerHeaterNum ).TempLowLimitEvapOut = rNumericArgs( 11 );
			ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate = rNumericArgs( 12 );
			if ( ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate == AutoSize ) {
				ChillerHeater( ChillerHeaterNum ).EvapVolFlowRateWasAutoSized = true;
			}
			ChillerHeater( ChillerHeaterNum ).CondVolFlowRate = rNumericArgs( 13 );
			if ( ChillerHeater( ChillerHeaterNum ).CondVolFlowRate == AutoSize ) {
				ChillerHeater( ChillerHeaterNum ).CondVolFlowRateWasAutoSized = true;
			}
			ChillerHeater( ChillerHeaterNum ).DesignHotWaterVolFlowRate = rNumericArgs( 14 );
			ChillerHeater( ChillerHeaterNum ).OpenMotorEff = rNumericArgs( 15 );
			ChillerHeater( ChillerHeaterNum ).OptPartLoadRatCooling = rNumericArgs( 16 );
			ChillerHeater( ChillerHeaterNum ).OptPartLoadRatClgHtg = rNumericArgs( 17 );
			ChillerHeater( ChillerHeaterNum ).SizFac = rNumericArgs( 18 );

			if ( ChillerHeater( ChillerHeaterNum ).SizFac <= 0.0 ) ChillerHeater( ChillerHeaterNum ).SizFac = 1.0;

			if ( ChillerHeater( ChillerHeaterNum ).OpenMotorEff < 0.0 || ChillerHeater( ChillerHeaterNum ).OpenMotorEff > 1.0 ) {
				ShowSevereError( "GetCurveInput: For " + cCurrentModuleObject + ": " + cAlphaArgs( 1 ) );
				ShowContinueError( cNumericFieldNames( 14 ) + " = " + RoundSigDigits( rNumericArgs( 14 ), 3 ) );
				ShowContinueError( cNumericFieldNames( 14 ) + " must be greater than or equal to zero" );
				ShowContinueError( cNumericFieldNames( 14 ) + " must be less than or equal to one" );
				CHErrorsFound = true;
			}

			// Check the CAP-FT, EIR-FT, and PLR curves and warn user if different from 1.0 by more than +-10%
			if ( ChillerHeater( ChillerHeaterNum ).ChillerCapFTCooling > 0 ) {
				CurveVal = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerCapFTCooling, ChillerHeater( ChillerHeaterNum ).TempRefEvapOutCooling, ChillerHeater( ChillerHeaterNum ).TempRefCondInCooling );
				if ( CurveVal > 1.10 || CurveVal < 0.90 ) {
					ShowWarningError( "Capacity ratio as a function of temperature curve output is not equal to 1.0" );
					ShowContinueError( "(+ or - 10%) at reference conditions for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "Curve output at reference conditions = " + TrimSigDigits( CurveVal, 3 ) );
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFTCooling > 0 ) {
				CurveVal = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerEIRFTCooling, ChillerHeater( ChillerHeaterNum ).TempRefEvapOutCooling, ChillerHeater( ChillerHeaterNum ).TempRefCondInCooling );
				if ( CurveVal > 1.10 || CurveVal < 0.90 ) {
					ShowWarningError( "Energy input ratio as a function of temperature curve output is not equal to 1.0" );
					ShowContinueError( "(+ or - 10%) at reference conditions for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "Curve output at reference conditions = " + TrimSigDigits( CurveVal, 3 ) );
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling > 0 ) {
				CurveVal = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling, 1.0 );

				if ( CurveVal > 1.10 || CurveVal < 0.90 ) {
					ShowWarningError( "Energy input ratio as a function of part-load ratio curve output is not equal to 1.0" );
					ShowContinueError( "(+ or - 10%) at reference conditions for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "Curve output at reference conditions = " + TrimSigDigits( CurveVal, 3 ) );
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling > 0 ) {
				FoundNegValue = false;
				for ( CurveCheck = 0; CurveCheck <= 10; ++CurveCheck ) {
					CurveValTmp = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling, double( CurveCheck / 10.0 ) );
					if ( CurveValTmp < 0.0 ) FoundNegValue = true;
					CurveValArray( CurveCheck + 1 ) = int( CurveValTmp * 100.0 ) / 100.0;
				}
				if ( FoundNegValue ) {
					ShowWarningError( "Energy input ratio as a function of part-load ratio curve shows negative values " );
					ShowContinueError( "for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "EIR as a function of PLR curve output at various part-load ratios shown below:" );
					ShowContinueError( "PLR   =  0.00   0.10   0.20   0.30   0.40   0.50   0.60   0.70   0.80   0.90   1.00" );
					gio::write( StringVar, "'Curve Output = '" );
					for ( CurveValPtr = 1; CurveValPtr <= 11; ++CurveValPtr ) {
						gio::write( StringVar, "(F7.2,$)" )
							<< CurveValArray( CurveValPtr );
					}
					gio::write( StringVar );
					ShowContinueError( StringVar );
					CHErrorsFound = true;
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerCapFTHeating > 0 ) {
				CurveVal = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerCapFTHeating, ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg, ChillerHeater( ChillerHeaterNum ).TempRefCondInClgHtg );
				if ( CurveVal > 1.10 || CurveVal < 0.90 ) {
					ShowWarningError( "Capacity ratio as a function of temperature curve output is not equal to 1.0" );
					ShowContinueError( "(+ or - 10%) at reference conditions for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "Curve output at reference conditions = " + TrimSigDigits( CurveVal, 3 ) );
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFTHeating > 0 ) {
				CurveVal = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerEIRFTHeating, ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg, ChillerHeater( ChillerHeaterNum ).TempRefCondInClgHtg );
				if ( CurveVal > 1.10 || CurveVal < 0.90 ) {
					ShowWarningError( "Energy input ratio as a function of temperature curve output is not equal to 1.0" );
					ShowContinueError( "(+ or - 10%) at reference conditions for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "Curve output at reference conditions = " + TrimSigDigits( CurveVal, 3 ) );
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating > 0 ) {
				CurveVal = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating, 1.0 );

				if ( CurveVal > 1.10 || CurveVal < 0.90 ) {
					ShowWarningError( "Energy input ratio as a function of part-load ratio curve output is not equal to 1.0" );
					ShowContinueError( "(+ or - 10%) at reference conditions for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "Curve output at reference conditions = " + TrimSigDigits( CurveVal, 3 ) );
				}
			}

			if ( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating > 0 ) {
				FoundNegValue = false;
				for ( CurveCheck = 0; CurveCheck <= 10; ++CurveCheck ) {
					CurveValTmp = CurveValue( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating, double( CurveCheck / 10.0 ) );
					if ( CurveValTmp < 0.0 ) FoundNegValue = true;
					CurveValArray( CurveCheck + 1 ) = int( CurveValTmp * 100.0 ) / 100.0;
				}
				if ( FoundNegValue ) {
					ShowWarningError( "Energy input ratio as a function of part-load ratio curve shows negative values " );
					ShowContinueError( "for " + cCurrentModuleObject + "= " + cAlphaArgs( 1 ) );
					ShowContinueError( "EIR as a function of PLR curve output at various part-load ratios shown below:" );
					ShowContinueError( "PLR          =    0.00   0.10   0.20   0.30   0.40   0.50   0.60   0.70   0.80   0.90   1.00" );
					gio::write( StringVar, "'Curve Output = '" );
					for ( CurveValPtr = 1; CurveValPtr <= 11; ++CurveValPtr ) {
						gio::write( StringVar, "(F7.2,$)" ) << CurveValArray( CurveValPtr );
					}
					gio::write( StringVar );
					ShowContinueError( StringVar );
					CHErrorsFound = true;
				}
			}

			GetCurveMinMaxValues( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating, ChillerHeater( ChillerHeaterNum ).MinPartLoadRatClgHtg, ChillerHeater( ChillerHeaterNum ).MaxPartLoadRatClgHtg );

			GetCurveMinMaxValues( ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling, ChillerHeater( ChillerHeaterNum ).MinPartLoadRatCooling, ChillerHeater( ChillerHeaterNum ).MaxPartLoadRatCooling );

		}

		if ( CHErrorsFound ) {
			ShowFatalError( "Errors found in processing input for " + cCurrentModuleObject );
		}

	}

	void
	WrapperSpecs::InitWrapper(
		Real64 const MyLoad, // Demand Load
		int const LoopNum // Loop Number Index
	)
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Daeho Kang, PNNL
		//       DATE WRITTEN   Feb. 2013
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		//  This subroutine is for initializations of the CentralHeatPumpSystem variables

		// METHODOLOGY EMPLOYED:
		//  Uses the status flags to trigger initializations.

		// REFERENCES:
		//  na

		// Using/Aliasing
		using DataGlobals::BeginEnvrnFlag;
		using DataGlobals::AnyEnergyManagementSystemInModel;
		using DataGlobals::InitConvTemp;
		using DataPlant::PlantLoop;
		using DataPlant::TypeOf_CentralGroundSourceHeatPump;
		using DataPlant::ScanPlantLoopsForObject;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using DataPlant::LoopFlowStatus_NeedyIfLoopOn;
		using InputProcessor::SameString;
		using Psychrometrics::PsyRhoAirFnPbTdbW;
		using CurveManager::GetCurveMinMaxValues;
		using PlantUtilities::InterConnectTwoPlantLoopSides;
		using PlantUtilities::InitComponentNodes;
		using PlantUtilities::SetComponentFlowRate;
		using EMSManager::iTemperatureSetPoint;
		using EMSManager::CheckIfNodeSetPointManagedByEMS;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "InitCGSHPHeatPump" );

		// INTERFACE BLOCK SPECIFICATIONS:
		//  na

		// DERIVED TYPE DEFINITIONS:
		//  na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int ChillerHeaterNum; // Chiller Heater index
		bool errFlag; // Err flag
		bool FatalError; // Fatal error indicator
		int CHWInletNodeNum; // Chilled water inlet node number
		int CHWOutletNodeNum; // Chilled water outlet node number
		int HWInletNodeNum; // Hot water inlet node number
		int HWOutletNodeNum; // Hot water outlet node number
		int GLHEInletNodeNum; // Condenser water inlet node number
		int GLHEOutletNodeNum; // Condenser water outlet node number
		Real64 rho; // local fluid density
		Real64 mdotCHW; // Chilled water mass flow rate
		Real64 mdotHW; // Hot water mass flow rate
		Real64 mdotGLHE; // Condenser water mass flow rate


		if ( this->MyWrapperFlag ) {
			// Locate the chillers on the plant loops for later usage
			errFlag = false;
			ScanPlantLoopsForObject( this->Name, TypeOf_CentralGroundSourceHeatPump, this->CWLoopNum, this->CWLoopSideNum, this->CWBranchNum, this->CWCompNum, _, _, _, this->CHWInletNodeNum, _, errFlag );

			ScanPlantLoopsForObject( this->Name, TypeOf_CentralGroundSourceHeatPump, this->HWLoopNum, this->HWLoopSideNum, this->HWBranchNum, this->HWCompNum, _, _, _, this->HWInletNodeNum, _, errFlag );

			ScanPlantLoopsForObject( this->Name, TypeOf_CentralGroundSourceHeatPump, this->GLHELoopNum, this->GLHELoopSideNum, this->GLHEBranchNum, this->GLHECompNum, _, _, _, this->GLHEInletNodeNum, _, errFlag );

			InterConnectTwoPlantLoopSides( this->CWLoopNum, this->CWLoopSideNum, this->GLHELoopNum, this->GLHELoopSideNum, TypeOf_CentralGroundSourceHeatPump, true );

			InterConnectTwoPlantLoopSides( this->HWLoopNum, this->HWLoopSideNum, this->GLHELoopNum, this->GLHELoopSideNum, TypeOf_CentralGroundSourceHeatPump, true );

			InterConnectTwoPlantLoopSides( this->CWLoopNum, this->CWLoopSideNum, this->HWLoopNum, this->HWLoopSideNum, TypeOf_CentralGroundSourceHeatPump, true );

			if ( this->VariableFlowCH ) {
				// Reset flow priority
				if ( LoopNum == this->CWLoopNum ) {
					PlantLoop( this->CWLoopNum ).LoopSide( this->CWLoopSideNum ).Branch( this->CWBranchNum ).Comp( this->CWCompNum ).FlowPriority = LoopFlowStatus_NeedyIfLoopOn;
				} else if ( LoopNum == this->HWLoopNum ) {
					PlantLoop( this->HWLoopNum ).LoopSide( this->HWLoopSideNum ).Branch( this->HWBranchNum ).Comp( this->HWCompNum ).FlowPriority = LoopFlowStatus_NeedyIfLoopOn;
				}

				// check if setpoint on outlet node - chilled water loop
				if ( Node( this->CHWOutletNodeNum ).TempSetPoint == SensedNodeFlagValue ) {
					if ( ! AnyEnergyManagementSystemInModel ) {
						if ( ! this->CoolSetPointErrDone ) {
							ShowWarningError( "Missing temperature setpoint on cooling side for CentralHeatPumpSystem named " + this->Name );
							ShowContinueError( "  A temperature setpoint is needed at the outlet node of a CentralHeatPumpSystem, use a SetpointManager" );
							ShowContinueError( "  The overall loop setpoint will be assumed for CentralHeatPumpSystem. The simulation continues ... " );
							this->CoolSetPointErrDone = true;
						}
					} else {
						// need call to EMS to check node
						FatalError = false; // but not really fatal yet, but should be.
						CheckIfNodeSetPointManagedByEMS( this->CHWOutletNodeNum, iTemperatureSetPoint, FatalError );
						if ( FatalError ) {
							if ( ! this->CoolSetPointErrDone ) {
								ShowWarningError( "Missing temperature setpoint on cooling side for CentralHeatPumpSystem named " + this->Name );
								ShowContinueError( "A temperature setpoint is needed at the outlet node of a CentralHeatPumpSystem " );
								ShowContinueError( "use a Setpoint Manager to establish a setpoint at the chiller side outlet node " );
								ShowContinueError( "or use an EMS actuator to establish a setpoint at the outlet node " );
								ShowContinueError( "The overall loop setpoint will be assumed for chiller side. The simulation continues ... " );
								this->CoolSetPointErrDone = true;
							}
						}
					}
					this->CoolSetPointSetToLoop = true;
					Node( this->CHWOutletNodeNum ).TempSetPoint = Node( PlantLoop( this->CWLoopNum ).TempSetPointNodeNum ).TempSetPoint;
				}

				if ( Node( this->HWOutletNodeNum ).TempSetPoint == SensedNodeFlagValue ) {
					if ( ! AnyEnergyManagementSystemInModel ) {
						if ( ! this->HeatSetPointErrDone ) {
							ShowWarningError( "Missing temperature setpoint on heating side for CentralHeatPumpSystem named " + this->Name );
							ShowContinueError( "  A temperature setpoint is needed at the outlet node of a CentralHeatPumpSystem, use a SetpointManager" );
							ShowContinueError( "  The overall loop setpoint will be assumed for CentralHeatPumpSystem. The simulation continues ... " );
							this->HeatSetPointErrDone = true;
						}
					} else {
						// need call to EMS to check node
						FatalError = false; // but not really fatal yet, but should be.
						CheckIfNodeSetPointManagedByEMS( this->HWOutletNodeNum, iTemperatureSetPoint, FatalError );
						if ( FatalError ) {
							if ( ! this->HeatSetPointErrDone ) {
								ShowWarningError( "Missing temperature setpoint on heating side for CentralHeatPumpSystem named " + this->Name );
								ShowContinueError( "A temperature setpoint is needed at the outlet node of a CentralHeatPumpSystem " );
								ShowContinueError( "use a Setpoint Manager to establish a setpoint at the chiller side outlet node " );
								ShowContinueError( "or use an EMS actuator to establish a setpoint at the outlet node " );
								ShowContinueError( "The overall loop setpoint will be assumed for chiller side. The simulation continues ... " );
								this->HeatSetPointErrDone = true;
							}
						}
					}
					this->HeatSetPointSetToLoop = true;
					Node( this->HWOutletNodeNum ).TempSetPoint = Node( PlantLoop( this->HWLoopNum ).TempSetPointNodeNum ).TempSetPoint;
				}
			}
			this->MyWrapperFlag = false;
		}

		CHWInletNodeNum = this->CHWInletNodeNum;
		CHWOutletNodeNum = this->CHWOutletNodeNum;
		HWInletNodeNum = this->HWInletNodeNum;
		HWOutletNodeNum = this->HWOutletNodeNum;
		GLHEInletNodeNum = this->GLHEInletNodeNum;
		GLHEOutletNodeNum = this->GLHEOutletNodeNum;

		if ( this->MyWrapperEnvrnFlag && BeginEnvrnFlag && ( PlantFirstSizesOkayToFinalize ) ) {

			if ( this->ControlMode == SmartMixing ) {

				this->CHWVolFlowRate = 0.0;
				this->HWVolFlowRate = 0.0;
				this->GLHEVolFlowRate = 0.0;

				for ( ChillerHeaterNum = 1; ChillerHeaterNum <= this->ChillerHeaterNums; ++ChillerHeaterNum ) {
					this->CHWVolFlowRate += this->ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate;
					this->HWVolFlowRate += this->ChillerHeater( ChillerHeaterNum ).DesignHotWaterVolFlowRate;
					this->GLHEVolFlowRate += this->ChillerHeater( ChillerHeaterNum ).CondVolFlowRate;
				}

				rho = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, InitConvTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );

				this->CHWMassFlowRateMax = this->CHWVolFlowRate * rho;
				this->HWMassFlowRateMax = this->HWVolFlowRate * rho;
				this->GLHEMassFlowRateMax = this->GLHEVolFlowRate * rho;

				InitComponentNodes( 0.0, this->CHWMassFlowRateMax, CHWInletNodeNum, CHWOutletNodeNum, this->CWLoopNum, this->CWLoopSideNum, this->CWBranchNum, this->CWCompNum );
				InitComponentNodes( 0.0, this->HWMassFlowRateMax, HWInletNodeNum, HWOutletNodeNum, this->HWLoopNum, this->HWLoopSideNum, this->HWBranchNum, this->HWCompNum );
				InitComponentNodes( 0.0, this->GLHEMassFlowRateMax, GLHEInletNodeNum, GLHEOutletNodeNum, this->GLHELoopNum, this->GLHELoopSideNum, this->GLHEBranchNum, this->GLHECompNum );

				// Initialize nodes for individual chiller heaters
				for ( ChillerHeaterNum = 1; ChillerHeaterNum <= this->ChillerHeaterNums; ++ChillerHeaterNum ) {
					this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRateMin = 0.0;
					this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRateMinAvail = 0.0;
					this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRateMax = rho * this->ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate;
					this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRateMaxAvail = rho * this->ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate;
					this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRate = 0.0;
					this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateMin = 0.0;
					this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateMinAvail = 0.0;
					this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateMax = rho * this->ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate;
					this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateMaxAvail = rho * this->ChillerHeater( ChillerHeaterNum ).EvapVolFlowRate;
					this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRate = 0.0;
					this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateRequest = 0.0;

				}

			}
			this->MyWrapperEnvrnFlag = false;
		}

		if ( ! BeginEnvrnFlag ) {
			this->MyWrapperEnvrnFlag = true;
		}

		if ( this->CoolSetPointSetToLoop ) {
			//IF (CurCoolingLoad > 0.0d0) THEN
			Node( this->CHWOutletNodeNum ).TempSetPoint = Node( PlantLoop( this->CWLoopNum ).TempSetPointNodeNum ).TempSetPoint;
		}
		//IF (CurHeatingLoad > 0.0d0) THEN
		if ( this->HeatSetPointSetToLoop ) {
			Node( this->HWOutletNodeNum ).TempSetPoint = Node( PlantLoop( this->HWLoopNum ).TempSetPointNodeNum ).TempSetPoint;
			//ENDIF
		}

		// Switch over the mass flow rate to the condenser loop, i.e., ground heat exchanger
		if ( LoopNum == this->CWLoopNum ) { // called for on cooling loop
			if ( MyLoad < -1.0 ) { // calling for cooling
				mdotCHW = Node( CHWInletNodeNum ).MassFlowRateMax;
			} else {
				mdotCHW = 0.0;
			}
			if ( this->WrapperHeatingLoad > 1.0 ) {
				mdotHW = Node( HWInletNodeNum ).MassFlowRateMax;
			} else {
				mdotHW = 0.0;
			}
			if ( ( MyLoad < -1.0 ) || ( this->WrapperHeatingLoad > 1.0 ) ) {
				mdotGLHE = Node( GLHEInletNodeNum ).MassFlowRateMax;
			} else {
				mdotGLHE = 0.0;
			}

		} else if ( LoopNum == this->HWLoopNum ) {
			if ( MyLoad > 1.0 ) {
				mdotHW = Node( HWInletNodeNum ).MassFlowRateMax;
			} else {
				mdotHW = 0.0;
			}
			if ( this->WrapperCoolingLoad > 1.0 ) {
				mdotCHW = Node( CHWInletNodeNum ).MassFlowRateMax;
			} else {
				mdotCHW = 0.0;
			}
			if ( ( MyLoad > 1.0 ) || ( this->WrapperCoolingLoad > 1.0 ) ) {
				mdotGLHE = Node( GLHEInletNodeNum ).MassFlowRateMax;
			} else {
				mdotGLHE = 0.0;
			}

		} else if ( LoopNum == this->GLHELoopNum ) {
			if ( this->WrapperCoolingLoad > 1.0 ) {
				mdotCHW = Node( CHWInletNodeNum ).MassFlowRateMax;
			} else {
				mdotCHW = 0.0;
			}
			if ( this->WrapperHeatingLoad > 1.0 ) {
				mdotHW = Node( HWInletNodeNum ).MassFlowRateMax;
			} else {
				mdotHW = 0.0;
			}
			if ( ( this->WrapperHeatingLoad > 1.0 ) || ( this->WrapperCoolingLoad > 1.0 ) ) {
				mdotGLHE = Node( GLHEInletNodeNum ).MassFlowRateMax;
			} else {
				mdotGLHE = 0.0;
			}
		}

		SetComponentFlowRate( mdotCHW, CHWInletNodeNum, CHWOutletNodeNum, this->CWLoopNum, this->CWLoopSideNum, this->CWBranchNum, this->CWCompNum );

		SetComponentFlowRate( mdotHW, HWInletNodeNum, HWOutletNodeNum, this->HWLoopNum, this->HWLoopSideNum, this->HWBranchNum, this->HWCompNum );

		SetComponentFlowRate( mdotGLHE, GLHEInletNodeNum, GLHEOutletNodeNum, this->GLHELoopNum, this->GLHELoopSideNum, this->GLHEBranchNum, this->GLHECompNum );

	}

	void
	WrapperSpecs::CalcChillerModel(
		int const EP_UNUSED( OpMode ), // Operation mode
		Real64 & EP_UNUSED( MyLoad ), // Operating load
		bool const EP_UNUSED( RunFlag ), // TRUE when chiller operating
		bool const EP_UNUSED( FirstIteration ), // TRUE when first iteration of timestep
		int const EP_UNUSED( EquipFlowCtrl ), // Flow control mode for the equipment
		int const EP_UNUSED( LoopNum ) // Plant loop number
	)
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Daeho Kang, PNNL
		//       DATE WRITTEN   Feb 2013
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		//  Simulate a ChillerHeaterPerformance:Electric:EIR using curve fit

		// METHODOLOGY EMPLOYED:
		//  Use empirical curve fits to model performance at off-reference conditions

		// REFERENCES:
		// 1. DOE-2 Engineers Manual, Version 2.1A, November 1982, LBL-11353

		// Using/Aliasing
		using DataGlobals::WarmupFlag;
		using DataGlobals::InitConvTemp;
		using DataHVACGlobals::SmallLoad;
		using CurveManager::CurveValue;
		using CurveManager::GetCurveMinMaxValues;
		using DataPlant::DeltaTempTol;
		using DataBranchAirLoopPlant::MassFlowTolerance;
		using ScheduleManager::GetCurrentScheduleValue;
		using General::TrimSigDigits;
		using General::RoundSigDigits;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "CalcChillerHeaterModel" );
		static std::string const RoutineNameElecEIRChiller( "CalcElectricEIRChillerModel" );

		//CHARACTER(len=*), PARAMETER :: OutputFormat  = '(F6.2)'

		// INTERFACE BLOCK SPECIFICATIONS
		//  na

		// DERIVED TYPE DEFINITIONS
		//  na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		static bool IsLoadCoolRemaining( true );
		static bool NextCompIndicator( false ); // Component indicator when identical chiller heaters exist
		int LoopSideNum; // Plant loop side which contains the current chiller (usually supply side)
		int CompNum( 0 ); // Component number in the loop  REAL(r64) :: FRAC
		int ChillerHeaterNum; // Chiller heater number
		int CurrentMode; // Current operational mode, cooling or simultaneous cooling and heating mode
		int IdenticalUnitCounter; // Pointer to count number of identical unit passed
		int IdenticalUnitRemaining; // Pointer to count number of identical unit available for a component
		Real64 FRAC; // Chiller cycling ratio
		Real64 MinPartLoadRat; // Min allowed operating fraction of full load
		Real64 MaxPartLoadRat; // Max allowed operating fraction of full load
		Real64 EvapInletTemp; // Evaporator inlet temperature [C]
		Real64 CondInletTemp; // Condenser inlet temperature [C]
		Real64 EvapOutletTempSetPoint; // Evaporator outlet temperature setpoint [C]
		Real64 AvailChillerCap; // Chiller available capacity at current operating conditions [W]
		Real64 ChillerRefCap; // Chiller reference capacity
		Real64 EvapDeltaTemp; // Evaporator temperature difference [C]
		Real64 ReferenceCOP; // Reference coefficient of performance, from user input
		Real64 PartLoadRat; // Operating part load ratio
		Real64 TempLowLimitEout; // Evaporator low temp. limit cut off [C]
		Real64 Cp; // Local fluid specific heat
		Real64 CondTempforCurve; // Condenser temp used for performance curve
		Real64 RemainingEvapMassPrevCH; // Bypass water from the previous variable chiller heater
		Real64 CoolingLoadToMeet; // Remaining cooling load the other chiller heaters should meet
		Real64 GLHEDensityRatio; // Fraction between starndarized density and local density in the condenser side
		Real64 CHWDensityRatio; // Fraction between starndarized density and local density in the chilled water side
		Real64 EvaporatorCapMin; // Minimum capacity of the evaporator
		Real64 EvaporatorLoad( 0.0 ); // Cooling load evaporator should meet
		Real64 HeatingPower; // Electric power use for heating
		Real64 CHWInletMassFlowRate; // Chilled water inlet mass flow rate
		Real64 CurAvailCHWMassFlowRate( 0.0 ); // Maximum available mass flow rate for current chiller heater
		Real64 EvapMassFlowRateCalc; // Evaporator mass flow rate calculated
		Real64 EvapDeltaTempCalc; // Evaporator temperature difference calculated
		Real64 EvapOutletTempCalc; // Evaporator outlet temperature calculated
		Real64 EvapMassFlowRate; // Actual evaporator mass flow rate
		Real64 CondMassFlowRate; // Condenser mass flow rate
		Real64 EvapOutletTemp; // Evaporator outlet temperature
		Real64 CondOutletTemp; // Condenser outlet temperature
		Real64 QCondenser; // Condenser heat transfer rate
		Real64 QEvaporator; // Evaporator heat transfer rate
		Real64 CHPower; // Evaporator power rate
		Real64 InitDensity; // Water density at the initial temperature
		Real64 EvapDensity; // Evaporator water density
		Real64 CondDensity; // Condenser water density
		Real64 ActualCOP; // Actual performance of individual chiller heater

		EvaporatorLoad = this->WrapperCoolingLoad;
		LoopSideNum = this->CWLoopSideNum;
		CHWInletMassFlowRate = Node( this->CHWInletNodeNum ).MassFlowRate;

		for ( ChillerHeaterNum = 1; ChillerHeaterNum <= this->ChillerHeaterNums; ++ChillerHeaterNum ) {

			// Initialize local variables for each chiller heater
			CurrentMode = 0;
			nsvChillerCapFT = 0.0;
			nsvChillerEIRFT = 0.0;
			nsvChillerEIRFPLR = 0.0;
			CoolingLoadToMeet = 0.0;
			nsvChillerPartLoadRatio = 0.0;
			nsvChillerCyclingRatio = 0.0;
			nsvChillerFalseLoadRate = 0.0;
			EvapMassFlowRate = 0.0;
			CondMassFlowRate = 0.0;
			CHPower = 0.0;
			HeatingPower = 0.0;
			QCondenser = 0.0;
			QEvaporator = 0.0;
			CondenserFanPower = 0.0;
			FRAC = 1.0;
			EvapDeltaTemp = 0.0;
			ActualCOP = 0.0;
			RemainingEvapMassPrevCH = 0.0;
			EvapInletTemp = Node( this->CHWInletNodeNum ).Temp;
			CondInletTemp = Node( this->GLHEInletNodeNum ).Temp;
			EvapOutletTemp = EvapInletTemp;
			CondOutletTemp = CondInletTemp;
			this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = 0;

			// Find proper schedule values
			if ( this->NumOfComp != this->ChillerHeaterNums ) { // Identical units exist
				if ( ChillerHeaterNum == 1 ) {
					IdenticalUnitCounter = 0;
					IdenticalUnitRemaining = 0;
					NextCompIndicator = false;
					CompNum = ChillerHeaterNum;
				}
				if ( NextCompIndicator ) {
					++CompNum;
				}
				if ( CompNum == 1 ) {
					if ( ChillerHeaterNum != this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = false;
					} else if ( ChillerHeaterNum == this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = true;
					}
				} else if ( CompNum > 1 ) {
					if ( ( ChillerHeaterNum - ( ( ChillerHeaterNum - 1 ) - IdenticalUnitCounter ) ) != this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = false;
					} else if ( ( ChillerHeaterNum - ( ( ChillerHeaterNum - 1 ) - IdenticalUnitCounter ) ) == this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = true;
					}
				}
				++IdenticalUnitCounter;
				IdenticalUnitRemaining = this->WrapperComp( CompNum ).WrapperIdenticalObjectNum - IdenticalUnitCounter;
				if ( IdenticalUnitRemaining == 0 ) IdenticalUnitCounter = 0;
			} else if ( this->NumOfComp == this->ChillerHeaterNums ) {
				++CompNum;
			}

			if ( CompNum > this->NumOfComp ) {
				ShowSevereError( "CalcChillerModel: ChillerHeater=\"" + this->Name + "\", calculated component number too big." );
				ShowContinueError( "Max number of components=[" + RoundSigDigits( this->NumOfComp ) + "], indicated component number=[" + RoundSigDigits( CompNum ) + "]." );
				ShowFatalError( "Program terminates due to preceding condition." );
			}

			// Check whether this chiller heater needs to run
			if ( EvaporatorLoad > 0.0 && ( GetCurrentScheduleValue( this->WrapperComp( CompNum ).CHSchedPtr ) > 0.0 ) ) {
				IsLoadCoolRemaining = true;

				// Calculate density ratios to adjust mass flow rates from initialized ones
				// Hot water temperature is known, but evaporator mass flow rates will be adjusted in the following "Do" loop
				InitDensity = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, InitConvTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );
				EvapDensity = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, EvapInletTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );
				CondDensity = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, CondInletTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );

				// Calculate density ratios to adjust mass flow rates from initialized ones
				CHWDensityRatio = EvapDensity / InitDensity;
				GLHEDensityRatio = CondDensity / InitDensity;
				CondMassFlowRate = this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateMaxAvail;
				EvapMassFlowRate = this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRateMaxAvail;
				EvapMassFlowRate *= CHWDensityRatio;
				CondMassFlowRate *= GLHEDensityRatio;

				// Check available flows from plant and then adjust as necessary
				if ( CurAvailCHWMassFlowRate == 0 ) { // The very first chiller heater to operate
					CurAvailCHWMassFlowRate = CHWInletMassFlowRate;
				} else if ( ChillerHeaterNum > 1 ) {
					CurAvailCHWMassFlowRate -= this->ChillerHeater( ChillerHeaterNum - 1 ).EvapOutletNode.MassFlowRate;
				}
				EvapMassFlowRate = min( CurAvailCHWMassFlowRate, EvapMassFlowRate );
			} else {
				IsLoadCoolRemaining = false;
				EvapMassFlowRate = 0.0;
				CondMassFlowRate = 0.0;
				CurrentMode = 0;
			}

			// Chiller heater is on when cooling load for this chiller heater remains and chilled water available
			if ( IsLoadCoolRemaining && ( EvapMassFlowRate > 0 ) && ( GetCurrentScheduleValue( this->WrapperComp( CompNum ).CHSchedPtr ) > 0 ) ) {
				// Indicate current mode is cooling-only mode. Simulataneous clg/htg mode will be set later
				CurrentMode = 1;

				// Assign proper performance curve information depending on the control mode
				// Cooling curve is used only for cooling-only mode, and the others (Simulataneous and heating) read the heating curve
				if ( SimulClgDominant || SimulHtgDominant ) {
					this->ChillerHeater( ChillerHeaterNum ).RefCap = this->ChillerHeater( ChillerHeaterNum ).RefCapClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).RefCOP = this->ChillerHeater( ChillerHeaterNum ).RefCOPClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOut = this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).TempRefCondOut = this->ChillerHeater( ChillerHeaterNum ).TempRefCondOutClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).OptPartLoadRat = this->ChillerHeater( ChillerHeaterNum ).OptPartLoadRatClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).CondMode = this->ChillerHeater( ChillerHeaterNum ).CondModeHeating;
					this->ChillerHeater( ChillerHeaterNum ).ChillerCapFT = this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTHeating;
					this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFT = this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFTHeating;
					this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR = this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating;
				} else {
					this->ChillerHeater( ChillerHeaterNum ).RefCap = this->ChillerHeater( ChillerHeaterNum ).RefCapCooling;
					this->ChillerHeater( ChillerHeaterNum ).RefCOP = this->ChillerHeater( ChillerHeaterNum ).RefCOPCooling;
					this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOut = this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOutCooling;
					this->ChillerHeater( ChillerHeaterNum ).TempRefCondIn = this->ChillerHeater( ChillerHeaterNum ).TempRefCondInCooling;
					this->ChillerHeater( ChillerHeaterNum ).TempRefCondOut = this->ChillerHeater( ChillerHeaterNum ).TempRefCondOutCooling;
					this->ChillerHeater( ChillerHeaterNum ).OptPartLoadRat = this->ChillerHeater( ChillerHeaterNum ).OptPartLoadRatCooling;
					this->ChillerHeater( ChillerHeaterNum ).CondMode = this->ChillerHeater( ChillerHeaterNum ).CondModeCooling;
					this->ChillerHeater( ChillerHeaterNum ).ChillerCapFT = this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTCooling;
					this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFT = this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFTCooling;
					this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR = this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRCooling;
				}

				// Only used to read curve values
				CondOutletTemp = this->ChillerHeater( ChillerHeaterNum ).TempRefCondOutCooling;
				if ( this->ChillerHeater( ChillerHeaterNum ).CondMode == "ENTERINGCONDENSER" ) {
					CondTempforCurve = CondInletTemp;
				} else if ( this->ChillerHeater( ChillerHeaterNum ).CondMode == "LEAVINGCONDENSER" ) {
					CondTempforCurve = CondOutletTemp;
				} else {
					ShowWarningError( "ChillerHeaterPerformance:Electric:EIR \"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\":" );
					ShowContinueError( "Chiller condensor temperature for curve fit are not decided, defalt value= cond_leaving (" + RoundSigDigits( nsvChillerCapFT, 3 ) + ")." );
					CondTempforCurve = CondOutletTemp;
				}

				// Bind local variables from the curve
				GetCurveMinMaxValues( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR, MinPartLoadRat, MaxPartLoadRat );
				ChillerRefCap = this->ChillerHeater( ChillerHeaterNum ).RefCap;
				ReferenceCOP = this->ChillerHeater( ChillerHeaterNum ).RefCOP;
				EvapOutletTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.Temp;
				TempLowLimitEout = this->ChillerHeater( ChillerHeaterNum ).TempLowLimitEvapOut;
				EvapOutletTempSetPoint = this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOutCooling;
				nsvChillerCapFT = CurveValue( this->ChillerHeater( ChillerHeaterNum ).ChillerCapFT, EvapOutletTempSetPoint, CondTempforCurve );

				if ( nsvChillerCapFT < 0 ) {
					if ( this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTError < 1 && ! WarmupFlag ) {
						++this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTError;
						ShowWarningError( "ChillerHeaterPerformance:Electric:EIR \"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\":" );
						ShowContinueError( " ChillerHeater Capacity as a Function of Temperature curve output is negative (" + RoundSigDigits( nsvChillerCapFT, 3 ) + ")." );
						ShowContinueError( " Negative value occurs using an Evaporator Outlet Temp of " + RoundSigDigits( EvapOutletTempSetPoint, 1 ) + " and a Condenser Inlet Temp of " + RoundSigDigits( CondInletTemp, 1 ) + '.' );
						ShowContinueErrorTimeStamp( " Resetting curve output to zero and continuing simulation." );
					} else if ( ! WarmupFlag ) {
						++this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTError;
						ShowRecurringWarningErrorAtEnd( "ChillerHeaterPerformance:Electric:EIR \"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\": ChillerHeater Capacity as a Function of Temperature curve output is negative warning continues...", this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTErrorIndex, nsvChillerCapFT, nsvChillerCapFT );
					}
					nsvChillerCapFT = 0.0;
				}

				// Calculate the specific heat of chilled water
				Cp = GetSpecificHeatGlycol( PlantLoop( this->CWLoopNum ).FluidName, EvapInletTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );

				// Calculate cooling load this chiller should meet and the other chillers are demanded
				EvapOutletTempSetPoint = Node( PlantLoop( this->CWLoopNum ).TempSetPointNodeNum ).TempSetPoint;
				EvaporatorCapMin = this->ChillerHeater( ChillerHeaterNum ).MinPartLoadRatCooling * this->ChillerHeater( ChillerHeaterNum ).RefCapCooling;
				CoolingLoadToMeet = min( this->ChillerHeater( ChillerHeaterNum ).RefCapCooling, max( std::abs( EvaporatorLoad ), EvaporatorCapMin ) );

				// Available chiller capacity as a function of temperature
				AvailChillerCap = ChillerRefCap * nsvChillerCapFT;

				// Part load ratio based on load and available chiller capacity, cap at max part load ratio
				if ( AvailChillerCap > 0 ) {
					PartLoadRat = max( 0.0, min( CoolingLoadToMeet / AvailChillerCap, MaxPartLoadRat ) );
				} else {
					PartLoadRat = 0.0;
				}

				if ( this->ChillerHeater( ChillerHeaterNum ).PossibleSubcooling ) {
					QEvaporator = CoolingLoadToMeet;
					EvapDeltaTemp = QEvaporator / EvapMassFlowRate / Cp;
					EvapOutletTemp = EvapInletTemp - EvapDeltaTemp;
				}

				// Set load this chiller heater should meet
				QEvaporator = min( CoolingLoadToMeet, ( AvailChillerCap * MaxPartLoadRat ) );
				EvapOutletTemp = EvapOutletTempSetPoint;
				EvapDeltaTemp = EvapInletTemp - EvapOutletTemp;

				// Calculate temperatures for constant flow and mass flow rates for variable flow
				if ( EvapMassFlowRate > MassFlowTolerance ) {
					if ( SimulHtgDominant ) { // Evaporator operates at full capacity for heating
						PartLoadRat = max( 0.0, min( ( ChillerRefCap / AvailChillerCap ), MaxPartLoadRat ) );
						QEvaporator = AvailChillerCap * PartLoadRat;
						EvapDeltaTemp = QEvaporator / EvapMassFlowRate / Cp;
						EvapOutletTemp = EvapInletTemp - EvapDeltaTemp;
					} else { // Cooling only mode or cooling dominant simultaneous htg/clg mode
						if ( this->VariableFlowCH ) { // Variable flow
							EvapMassFlowRateCalc = QEvaporator / EvapDeltaTemp / Cp;
							if ( EvapMassFlowRateCalc > EvapMassFlowRate ) {
								EvapMassFlowRateCalc = EvapMassFlowRate;
								EvapDeltaTempCalc = QEvaporator / EvapMassFlowRate / Cp;
								EvapOutletTemp = EvapInletTemp - EvapDeltaTempCalc;
								if ( EvapDeltaTempCalc > EvapDeltaTemp ) {
									EvapDeltaTempCalc = EvapDeltaTemp;
									QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
								}
							}
							EvapMassFlowRate = EvapMassFlowRateCalc;
						} else { // Constant Flow
							EvapDeltaTempCalc = QEvaporator / EvapMassFlowRate / Cp;
							EvapOutletTempCalc = EvapInletTemp - EvapDeltaTemp;
							if ( EvapOutletTempCalc > EvapOutletTemp ) { // Load to meet should be adjusted
								EvapOutletTempCalc = EvapOutletTemp;
								QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
							}
							EvapOutletTemp = EvapOutletTempCalc;
						} // End of flow control decision
					} // End of operation mode
				} else {
					QEvaporator = 0.0;
					EvapOutletTemp = EvapInletTemp;
				}

				// Check evaporator temperature low limit and adjust capacity if needed
				if ( EvapOutletTemp < TempLowLimitEout ) {
					if ( ( EvapInletTemp - TempLowLimitEout ) > DeltaTempTol ) {
						EvapOutletTemp = TempLowLimitEout;
						EvapDeltaTemp = EvapInletTemp - EvapOutletTemp;
						QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
					} else {
						QEvaporator = 0.0;
						EvapOutletTemp = EvapInletTemp;
					}
				}

				// Check if the outlet temperature exceeds the node minimum temperature and adjust capacity if needed
				if ( EvapOutletTemp < this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin ) {
					if ( ( this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp - this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin ) > DeltaTempTol ) {
						EvapOutletTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin;
						EvapDeltaTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin - EvapOutletTemp;
						QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
					} else {
						QEvaporator = 0.0;
						EvapOutletTemp = EvapInletTemp;
					}
				}

				// Calculate part load once more since evaporator capacity might be modified
				if ( AvailChillerCap > 0.0 ) {
					PartLoadRat = max( 0.0, min( ( QEvaporator / AvailChillerCap ), MaxPartLoadRat ) );
				} else {
					PartLoadRat = 0.0;
				}

				// Chiller cycles below minimum part load ratio, FRAC = amount of time chiller is ON during this time step
				if ( PartLoadRat < MinPartLoadRat ) FRAC = min( 1.0, ( PartLoadRat / MinPartLoadRat ) );

				// set the module level variable used for reporting FRAC
				nsvChillerCyclingRatio = FRAC;

				// Chiller is false loading below PLR = minimum unloading ratio, find PLR used for energy calculation
				if ( AvailChillerCap > 0.0 ) {
					PartLoadRat = max( PartLoadRat, MinPartLoadRat );
				} else {
					PartLoadRat = 0.0;
				}

				// set the module level variable used for reporting PLR
				nsvChillerPartLoadRatio = PartLoadRat;

				// calculate the load due to false loading on chiller over and above water side load
				nsvChillerFalseLoadRate = ( AvailChillerCap * PartLoadRat * FRAC ) - QEvaporator;
				if ( nsvChillerFalseLoadRate < SmallLoad ) {
					nsvChillerFalseLoadRate = 0.0;
				}

				// Determine chiller compressor power and transfer heat calculation
				nsvChillerEIRFT = max( 0.0, CurveValue( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFT, EvapOutletTemp, CondTempforCurve ) );
				nsvChillerEIRFPLR = max( 0.0, CurveValue( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR, PartLoadRat ) );
				CHPower = ( AvailChillerCap / ReferenceCOP ) * nsvChillerEIRFPLR * nsvChillerEIRFT * FRAC;
				QCondenser = CHPower * this->ChillerHeater( ChillerHeaterNum ).OpenMotorEff + QEvaporator + nsvChillerFalseLoadRate;
				ActualCOP = ( QEvaporator + nsvChillerFalseLoadRate ) / CHPower;

				if ( CondMassFlowRate > MassFlowTolerance ) {
					Cp = GetSpecificHeatGlycol( PlantLoop( this->GLHELoopNum ).FluidName, CondInletTemp, PlantLoop( this->GLHELoopNum ).FluidIndex, RoutineNameElecEIRChiller );
					CondOutletTemp = QCondenser / CondMassFlowRate / Cp + CondInletTemp;
				} else {
					ShowSevereError( "CalcChillerheaterModel: Condenser flow = 0, for Chillerheater=" + this->ChillerHeater( ChillerHeaterNum ).Name );
					ShowContinueErrorTimeStamp( "" );
				}

				// Determine load next chillers should meet
				if ( EvaporatorLoad < QEvaporator ) {
					EvaporatorLoad = 0.0; // No remaining load so the rest will be off
				} else {
					EvaporatorLoad -= QEvaporator;
				}

				// Initialize reporting variable when this chiller doesn't need to operate
				if ( QEvaporator == 0.0 ) {
					CurrentMode = 0;
					nsvChillerPartLoadRatio = 0.0;
					nsvChillerCyclingRatio = 0.0;
					nsvChillerFalseLoadRate = 0.0;
					EvapMassFlowRate = 0.0;
					CondMassFlowRate = 0.0;
					CHPower = 0.0;
					QCondenser = 0.0;
					CondenserFanPower = 0.0;
					EvapOutletTemp = EvapInletTemp;
					CondOutletTemp = CondInletTemp;
					EvaporatorLoad = 0.0;
				}

			} // End of calculation for cooling

			// Set variables to the arrays
			this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.MassFlowRate = EvapMassFlowRate;
			this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.MassFlowRate = CondMassFlowRate;
			this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.Temp = EvapOutletTemp;
			this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp = EvapInletTemp;
			this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.Temp = CondOutletTemp;
			this->ChillerHeater( ChillerHeaterNum ).CondInletNode.Temp = CondInletTemp;
			this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = CurrentMode;
			this->ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatio = nsvChillerPartLoadRatio;
			this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatio = nsvChillerCyclingRatio;
			this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate = nsvChillerFalseLoadRate;
			this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFT = nsvChillerCapFT;
			this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFT = nsvChillerEIRFT;
			this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLR = nsvChillerEIRFPLR;
			this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPower = CHPower;
			this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower = HeatingPower;
			this->ChillerHeaterReport( ChillerHeaterNum ).QEvap = QEvaporator;
			this->ChillerHeaterReport( ChillerHeaterNum ).QCond = QCondenser;
			this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp = EvapOutletTemp;
			this->ChillerHeaterReport( ChillerHeaterNum ).EvapInletTemp = EvapInletTemp;
			this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp = CondOutletTemp;
			this->ChillerHeaterReport( ChillerHeaterNum ).CondInletTemp = CondInletTemp;
			this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot = EvapMassFlowRate;
			this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot = CondMassFlowRate;
			this->ChillerHeaterReport( ChillerHeaterNum ).ActualCOP = ActualCOP;

			if ( SimulClgDominant || SimulHtgDominant ) { // Store for using these cooling side data in the hot water loop
				this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = CurrentMode;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatioSimul = nsvChillerPartLoadRatio;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatioSimul = nsvChillerCyclingRatio;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRateSimul = nsvChillerFalseLoadRate;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFTSimul = nsvChillerCapFT;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFTSimul = nsvChillerEIRFT;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLRSimul = nsvChillerEIRFPLR;
				this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPowerSimul = CHPower;
				this->ChillerHeaterReport( ChillerHeaterNum ).QEvapSimul = QEvaporator;
				this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTempSimul = EvapOutletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).EvapInletTempSimul = EvapInletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul = EvapMassFlowRate;
				if ( SimulClgDominant ) {
					this->ChillerHeaterReport( ChillerHeaterNum ).QCondSimul = QCondenser;
					this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTempSimul = CondOutletTemp;
					this->ChillerHeaterReport( ChillerHeaterNum ).CondInletTempSimul = CondInletTemp;
					this->ChillerHeaterReport( ChillerHeaterNum ).CondmdotSimul = CondMassFlowRate;
				}
			}
		}

	}

	void
	WrapperSpecs::CalcChillerHeaterModel(
		int const EP_UNUSED( OpMode ), // Operation mode
		Real64 & EP_UNUSED( MyLoad ), // Heating load plant should meet
		bool const EP_UNUSED( RunFlag ), // TRUE when chiller operating
		bool const EP_UNUSED( FirstIteration ), // TRUE when first iteration of timestep
		int const EP_UNUSED( EquipFlowCtrl ), // Flow control mode for the equipment
		int const EP_UNUSED( LoopNum ) // Loop number
	)
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Daeho Kang, PNNL
		//       DATE WRITTEN   Feb 2013
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		//  Simulate a ChillerHeaterPerformance:Electric:EIR using curve fit

		// METHODOLOGY EMPLOYED:
		//  Use empirical curve fits to model performance at off-reference conditions

		// REFERENCES:
		// 1. DOE-2 Engineers Manual, Version 2.1A, November 1982, LBL-11353

		// Using/Aliasing
		using DataGlobals::WarmupFlag;
		using DataGlobals::InitConvTemp;
		using DataHVACGlobals::SmallLoad;
		using CurveManager::CurveValue;
		using CurveManager::GetCurveMinMaxValues;
		using DataPlant::DeltaTempTol;
		using ScheduleManager::GetCurrentScheduleValue;
		using General::TrimSigDigits;
		using General::RoundSigDigits;
		using DataBranchAirLoopPlant::MassFlowTolerance;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "CalcChillerHeaterModel" );
		static std::string const RoutineNameElecEIRChiller( "CalcElectricEIRChillerModel" );

		//CHARACTER(len=*), PARAMETER :: OutputFormat  = '(F6.2)'

		// INTERFACE BLOCK SPECIFICATIONS
		//  na

		// DERIVED TYPE DEFINITIONS
		//  na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		static bool IsLoadHeatRemaining( true ); // Ture if heating load remains for this chiller heater
		static bool NextCompIndicator( false ); // Component indicator when identical chiller heaters exist
		int LoopSideNum; // Plant loop side which contains the current chiller (usually supply side)
		int CompNum( 0 ); // Component number
		int ChillerHeaterNum; // Chiller heater number
		int CurrentMode; // Current operational mode, heating or simultaneous cooling and heating mode
		int IdenticalUnitCounter; // Pointer to count number of identical unit passed
		int IdenticalUnitRemaining; // Pointer to count number of identical unit available for a component
		Real64 Cp; // Local fluid specific heat
		Real64 CondTempforCurve; // Reference condenser temperature for the performance curve reading
		Real64 FRAC; // Chiller cycling ratio
		Real64 MinPartLoadRat; // Min allowed operating fraction of full load
		Real64 MaxPartLoadRat; // Max allowed operating fraction of full load
		Real64 EvapInletTemp; // Evaporator inlet temperature [C]
		Real64 CondInletTemp; // Condenser inlet temperature [C]
		Real64 EvapOutletTempSetPoint; // Condenser outlet temperature setpoint [C]
		Real64 AvailChillerCap; // Chiller available capacity at current operating conditions [W]
		Real64 ChillerRefCap; // Chiller reference capacity
		Real64 EvapDeltaTemp; // Evaporator temperature difference [C]
		Real64 CondDeltaTemp; // Condenser temperature difference [C]
		Real64 ReferenceCOP; // Reference coefficient of performance, from user input
		Real64 PartLoadRat; // Operating part load ratio
		Real64 TempLowLimitEout; // Evaporator low temp. limit cut off [C]
		Real64 CondenserLoad( 0.0 ); // Remaining heating load that this wrapper should meet
		Real64 HeatingLoadToMeet; // Heating load that this chiller heater should meet
		Real64 GLHEDensityRatio; // The density ratio of source water to the initialized source water
		Real64 HWDensityRatio; // The density ratio of hot water to the initialized hot water
		Real64 CondenserCapMin; // Minimum condenser capacity
		Real64 CoolingPower; // Evaporator cooling power to produce heat for heating
		Real64 HWInletMassFlowRate; // Hot water inlet mass flow rate
		Real64 CurAvailHWMassFlowRate( 0.0 ); // Maximum available hot water mass within the wrapper bank
		Real64 CondDeltaTempCalc; // Temperature differnece between condenser inlet and outlet calculated
		Real64 CondOutletTempCalc; // Condenser outlet temperature calculated
		Real64 CondMassFlowRateCalc; // Condenser mass flow rate calculated
		Real64 EvapMassFlowRate; // Evaporator mass flow rate through this chiller heater
		Real64 CondMassFlowRate; // Condenser mass flow rate through this chiller heater
		Real64 EvapOutletTemp; // Evaporator outlet temperature
		Real64 CondOutletTemp; // Condenser outlet temperature
		Real64 QCondenser; // Condenser heat transfer rate
		Real64 QEvaporator; // Evaporator heat transfer rate
		Real64 CHPower; // Evaporator compressor power added to heating power
		Real64 InitDensity; // Water density at the initial temperature
		Real64 EvapDensity; // Evaporator water density
		Real64 CondDensity; // Condenser water density
		Real64 ActualCOP; // Actual performance of individual chiller heater

		CondenserLoad = this->WrapperHeatingLoad;
		LoopSideNum = this->HWLoopSideNum;
		HWInletMassFlowRate = Node( this->HWInletNodeNum ).MassFlowRate;

		// Flow
		for ( ChillerHeaterNum = 1; ChillerHeaterNum <= this->ChillerHeaterNums; ++ChillerHeaterNum ) {

			// Set module level inlet and outlet nodes and initialize other local variables
			CurrentMode = 0;
			HeatingLoadToMeet = 0.0;
			nsvChillerPartLoadRatio = 0.0;
			nsvChillerCyclingRatio = 0.0;
			nsvChillerFalseLoadRate = 0.0;
			EvapMassFlowRate = 0.0;
			CondMassFlowRate = 0.0;
			CHPower = 0.0;
			QCondenser = 0.0;
			QEvaporator = 0.0;
			CondenserFanPower = 0.0;
			FRAC = 1.0;
			CondDeltaTemp = 0.0;
			EvapDeltaTemp = 0.0;
			CoolingPower = 0.0;
			ActualCOP = 0.0;
			EvapInletTemp = Node( this->GLHEInletNodeNum ).Temp;
			CondInletTemp = Node( this->HWInletNodeNum ).Temp;
			EvapOutletTemp = EvapInletTemp;
			CondOutletTemp = CondInletTemp;

			// Find proper schedule values
			if ( this->NumOfComp != this->ChillerHeaterNums ) { // Identical units exist
				if ( ChillerHeaterNum == 1 ) {
					IdenticalUnitCounter = 0;
					IdenticalUnitRemaining = 0;
					NextCompIndicator = false;
					CompNum = ChillerHeaterNum;
				}
				if ( NextCompIndicator ) {
					++CompNum;
				}
				if ( CompNum == 1 ) {
					if ( ChillerHeaterNum != this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = false;
					} else if ( ChillerHeaterNum == this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = true;
					}
				} else if ( CompNum > 1 ) {
					if ( ( ChillerHeaterNum - ( ( ChillerHeaterNum - 1 ) - IdenticalUnitCounter ) ) != this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = false;
					} else if ( ( ChillerHeaterNum - ( ( ChillerHeaterNum - 1 ) - IdenticalUnitCounter ) ) == this->WrapperComp( CompNum ).WrapperIdenticalObjectNum ) {
						NextCompIndicator = true;
					}
				}
				++IdenticalUnitCounter;
				IdenticalUnitRemaining = this->WrapperComp( CompNum ).WrapperIdenticalObjectNum - IdenticalUnitCounter;
				if ( IdenticalUnitRemaining == 0 ) IdenticalUnitCounter = 0;
			} else if ( this->NumOfComp == this->ChillerHeaterNums ) {
				++CompNum;
			}

			// Check to see if this chiiller heater needs to run
			if ( CondenserLoad > 0.0 && ( GetCurrentScheduleValue( this->WrapperComp( CompNum ).CHSchedPtr ) > 0 ) ) {
				IsLoadHeatRemaining = true;

				// Calculate density ratios to adjust mass flow rates from initialized ones
				// Hot water temperature is known, but condenser mass flow rates will be adjusted in the following "Do" loop
				InitDensity = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, InitConvTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );
				EvapDensity = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, EvapInletTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );
				CondDensity = GetDensityGlycol( PlantLoop( this->CWLoopNum ).FluidName, CondInletTemp, PlantLoop( this->CWLoopNum ).FluidIndex, RoutineName );

				// Calculate density ratios to adjust mass flow rates from initialized ones
				HWDensityRatio = CondDensity / InitDensity;
				GLHEDensityRatio = EvapDensity / InitDensity;
				EvapMassFlowRate = this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.MassFlowRateMaxAvail;
				CondMassFlowRate = this->ChillerHeater( ChillerHeaterNum ).CondInletNode.MassFlowRateMaxAvail;
				EvapMassFlowRate *= GLHEDensityRatio;
				CondMassFlowRate *= HWDensityRatio;

				// Check flows from plant to adjust as necessary
				if ( CurAvailHWMassFlowRate == 0 ) { // First chiller heater which is on
					CurAvailHWMassFlowRate = HWInletMassFlowRate;
				} else if ( ChillerHeaterNum > 1 ) {
					CurAvailHWMassFlowRate -= this->ChillerHeater( ChillerHeaterNum - 1 ).CondOutletNode.MassFlowRate;
				}
				CondMassFlowRate = min( CurAvailHWMassFlowRate, CondMassFlowRate );

				// It is not enforced to be the smaller of CH max temperature and plant temp setpoint.
				// Hot water temperatures at the individual CHs' outlet may be greater than plant setpoint temp,
				// but should be lower than the CHs max temp
				CondOutletTemp = this->ChillerHeater( ChillerHeaterNum ).TempRefCondOutClgHtg;
				CondDeltaTemp = CondOutletTemp - CondInletTemp;

				if ( CondDeltaTemp < 0.0 ) { // Hot water temperature is greater than the maximum
					if ( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRRefTempErrorIndex == 0 ) {
						ShowSevereMessage( "CalcChillerHeaterModel: ChillerHeaterPerformance:Electric:EIR=\"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\", DeltaTemp < 0" );
						ShowContinueError( " Reference Simultaneous Cooling-Heating Mode Leaving Condenser Water Temperature [" + RoundSigDigits( CondOutletTemp, 1 ) + ']' );
						ShowContinueError( "is below condenser inlet temperature of [" + RoundSigDigits( CondInletTemp, 1 ) + "]." );
						ShowContinueErrorTimeStamp( "" );
						ShowContinueError( " Reset reference temperature to one greater than the inlet temperature " );
					}
					ShowRecurringSevereErrorAtEnd( "ChillerHeaterPerformance:Electric:EIR=\"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\": Reference temperature problems continue.", this->ChillerHeater( ChillerHeaterNum ).ChillerEIRRefTempErrorIndex, CondDeltaTemp, CondDeltaTemp, _, "deltaC", "deltaC" );
					QCondenser = 0.0;
					IsLoadHeatRemaining = false;
				}

				if ( ChillerHeaterNum > 1 ) {
					// Operation mode needs to be set in a simultaneous clg/htg mode
					// Always off even heating load remains if this CH is assumed to be off in the loop 1
					if ( SimulClgDominant ) {
						if ( this->ChillerHeaterReport( ChillerHeaterNum ).QEvapSimul == 0.0 ) {
							CurrentMode = 0;
							IsLoadHeatRemaining = false;
						} else { // Heat recovery
							CurrentMode = 3;
						}
					}
				} // End of simulataneous clg/htg mode detemination

			} else { // chiller heater is off
				IsLoadHeatRemaining = false;
				CondMassFlowRate = 0.0;
				EvapMassFlowRate = 0.0;
				CurrentMode = 0;
				if ( SimulClgDominant ) {
					if ( this->ChillerHeaterReport( ChillerHeaterNum ).QEvapSimul > 0.0 ) {
						CurrentMode = 4; // Simultaneous cooling dominant mode: 4
					}
				} // End of mode determination
			} // End of system operation determinatoin

			if ( IsLoadHeatRemaining && CondMassFlowRate > 0.0 && ( GetCurrentScheduleValue( this->WrapperComp( CompNum ).CHSchedPtr ) > 0 ) ) { // System is on
				// Operation mode
				if ( SimulHtgDominant ) {
					if ( this->ChillerHeaterReport( ChillerHeaterNum ).QEvapSimul == 0.0 ) {
						CurrentMode = 5; // No cooling necessary
					} else { // Heat recovery mode. Both chilled water and hot water loops are connected. No condenser flow.
						CurrentMode = 3;
					}
				}

				// Mode 3 and 5 use cooling side data stored from the chilled water loop
				// Mode 4 uses all data from the chilled water loop due to no heating demand
				if ( SimulClgDominant || CurrentMode == 3 ) {
					CurrentMode = 3;
					Cp = GetSpecificHeatGlycol( PlantLoop( this->HWLoopNum ).FluidName, CondInletTemp, PlantLoop( this->HWLoopNum ).FluidIndex, RoutineName );

					QCondenser = this->ChillerHeaterReport( ChillerHeaterNum ).QCondSimul;

					if ( this->VariableFlowCH ) { // Variable flow
						CondMassFlowRateCalc = QCondenser / CondDeltaTemp / Cp;
						if ( CondMassFlowRateCalc > CondMassFlowRate ) {
							CondMassFlowRateCalc = CondMassFlowRate;
							CondDeltaTempCalc = QCondenser / CondMassFlowRate / Cp;
							if ( CondDeltaTempCalc > CondDeltaTemp ) { // Load to meet should be adjusted
								CondDeltaTempCalc = CondDeltaTemp;
								QCondenser = CondMassFlowRate * Cp * CondDeltaTemp;
							}
						}
						CondMassFlowRate = CondMassFlowRateCalc;
					} else { // Constant flow control
						CondDeltaTempCalc = QCondenser / CondMassFlowRate / Cp;
						CondOutletTempCalc = CondDeltaTempCalc + CondInletTemp;
						if ( CondOutletTempCalc > CondOutletTemp ) {
							CondOutletTempCalc = CondOutletTemp;
							QCondenser = CondMassFlowRate * Cp * CondDeltaTemp;
						}
						CondOutletTemp = CondOutletTempCalc;
					}

				} else { // Either Mode 2 or 3 or 5
					if ( SimulHtgDominant ) {
						CurrentMode = 5;
					} else {
						CurrentMode = 2;
					}

					nsvChillerCapFT = 0.0;
					nsvChillerEIRFT = 0.0;
					nsvChillerEIRFPLR = 0.0;

					// Assign curve values to local data array
					this->ChillerHeater( ChillerHeaterNum ).RefCap = this->ChillerHeater( ChillerHeaterNum ).RefCapClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).RefCOP = this->ChillerHeater( ChillerHeaterNum ).RefCOPClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOut = this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).TempRefCondOut = this->ChillerHeater( ChillerHeaterNum ).TempRefCondOutClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).OptPartLoadRat = this->ChillerHeater( ChillerHeaterNum ).OptPartLoadRatClgHtg;
					this->ChillerHeater( ChillerHeaterNum ).CondMode = this->ChillerHeater( ChillerHeaterNum ).CondModeHeating;
					this->ChillerHeater( ChillerHeaterNum ).ChillerCapFT = this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTHeating;
					this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFT = this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFTHeating;
					this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR = this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLRHeating;

					if ( this->ChillerHeater( ChillerHeaterNum ).CondMode == "ENTERINGCONDENSER" ) {
						CondTempforCurve = CondInletTemp;
					} else if ( this->ChillerHeater( ChillerHeaterNum ).CondMode == "LEAVINGCONDENSER" ) {
						CondTempforCurve = this->ChillerHeater( ChillerHeaterNum ).TempRefCondOutClgHtg; //!CondOutletTemp
					} else {
						ShowWarningError( "ChillerHeaterPerformance:Electric:EIR \"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\":" );
						ShowContinueError( "Chiller condensor temperature for curve fit are not decided, defalt value= cond_leaving (" + RoundSigDigits( nsvChillerCapFT, 3 ) + ")." );
						CondTempforCurve = Node( PlantLoop( this->HWLoopNum ).TempSetPointNodeNum ).TempSetPoint;
					}

					GetCurveMinMaxValues( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR, MinPartLoadRat, MaxPartLoadRat );
					ChillerRefCap = this->ChillerHeater( ChillerHeaterNum ).RefCap;
					ReferenceCOP = this->ChillerHeater( ChillerHeaterNum ).RefCOP;
					EvapOutletTemp = this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg;
					TempLowLimitEout = this->ChillerHeater( ChillerHeaterNum ).TempLowLimitEvapOut;
					EvapOutletTempSetPoint = this->ChillerHeater( ChillerHeaterNum ).TempRefEvapOutClgHtg;
					nsvChillerCapFT = CurveValue( this->ChillerHeater( ChillerHeaterNum ).ChillerCapFT, EvapOutletTempSetPoint, CondTempforCurve );

					if ( nsvChillerCapFT < 0 ) {
						if ( this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTError < 1 && ! WarmupFlag ) {
							++this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTError;
							ShowWarningError( "ChillerHeaterPerformance:Electric:EIR \"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\":" );
							ShowContinueError( " ChillerHeater Capacity as a Function of Temperature curve output is negative (" + RoundSigDigits( nsvChillerCapFT, 3 ) + ")." );
							ShowContinueError( " Negative value occurs using an Evaporator Outlet Temp of " + RoundSigDigits( EvapOutletTempSetPoint, 1 ) + " and a Condenser Inlet Temp of " + RoundSigDigits( CondInletTemp, 1 ) + '.' );
							ShowContinueErrorTimeStamp( " Resetting curve output to zero and continuing simulation." );
						} else if ( ! WarmupFlag ) {
							++this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTError;
							ShowRecurringWarningErrorAtEnd( "ChillerHeaterPerformance:Electric:EIR \"" + this->ChillerHeater( ChillerHeaterNum ).Name + "\": ChillerHeater Capacity as a Function of Temperature curve output is negative warning continues...", this->ChillerHeater( ChillerHeaterNum ).ChillerCapFTErrorIndex, nsvChillerCapFT, nsvChillerCapFT );
						}
						nsvChillerCapFT = 0.0;
					}

					// Available chiller capacity as a function of temperature
					AvailChillerCap = ChillerRefCap * nsvChillerCapFT;

					// Part load ratio based on reference capacity and available chiller capacity
					if ( AvailChillerCap > 0 ) {
						PartLoadRat = max( 0.0, min( ( ChillerRefCap / AvailChillerCap ), MaxPartLoadRat ) );
					} else {
						PartLoadRat = 0.0;
					}

					Cp = GetSpecificHeatGlycol( PlantLoop( this->HWLoopNum ).FluidName, this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp, PlantLoop( this->HWLoopNum ).FluidIndex, RoutineName );

					// Calculate evaporator heat transfer
					if ( EvapMassFlowRate > MassFlowTolerance ) {
						QEvaporator = AvailChillerCap * PartLoadRat;
						EvapDeltaTemp = QEvaporator / EvapMassFlowRate / Cp;
						EvapOutletTemp = EvapInletTemp - EvapDeltaTemp;
					}

					// Check that the evaporator outlet temp honors both plant loop temp low limit and also the chiller low limit
					if ( EvapOutletTemp < TempLowLimitEout ) {
						if ( ( this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp - TempLowLimitEout ) > DeltaTempTol ) {
							EvapOutletTemp = TempLowLimitEout;
							EvapDeltaTemp = this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp - EvapOutletTemp;
							QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
						} else {
							EvapOutletTemp = this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp;
							EvapDeltaTemp = this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp - EvapOutletTemp;
							QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
						}
					}

					if ( EvapOutletTemp < this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin ) {
						if ( ( this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp - this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin ) > DeltaTempTol ) {
							EvapOutletTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin;
							EvapDeltaTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin - EvapOutletTemp;
							QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
						} else {
							EvapOutletTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin;
							EvapDeltaTemp = this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.TempMin - EvapOutletTemp;
							QEvaporator = EvapMassFlowRate * Cp * EvapDeltaTemp;
						}
					}

					// Evaporator operates at full load
					if ( AvailChillerCap > 0.0 ) {
						PartLoadRat = max( 0.0, min( ( QEvaporator / AvailChillerCap ), MaxPartLoadRat ) );
					} else {
						PartLoadRat = 0.0;
					}

					// Chiller cycles below minimum part load ratio, FRAC = amount of time chiller is ON during this time step
					if ( PartLoadRat < MinPartLoadRat ) FRAC = min( 1.0, ( PartLoadRat / MinPartLoadRat ) );
					if ( FRAC <= 0.0 ) FRAC = 1.0; // CR 9303 COP reporting issue, it should be greater than zero in this routine
					nsvChillerCyclingRatio = FRAC;

					// Chiller is false loading below PLR = minimum unloading ratio, find PLR used for energy calculation
					if ( AvailChillerCap > 0.0 ) {
						PartLoadRat = max( PartLoadRat, MinPartLoadRat );
					} else {
						PartLoadRat = 0.0;
					}
					// Evaporator part load ratio
					nsvChillerPartLoadRatio = PartLoadRat;

					// calculate the load due to false loading on chiller over and above water side load
					nsvChillerFalseLoadRate = ( AvailChillerCap * PartLoadRat * FRAC ) - QEvaporator;
					if ( nsvChillerFalseLoadRate < SmallLoad ) {
						nsvChillerFalseLoadRate = 0.0;
					}

					nsvChillerEIRFT = max( 0.0, CurveValue( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFT, EvapOutletTemp, CondTempforCurve ) );
					nsvChillerEIRFPLR = max( 0.0, CurveValue( this->ChillerHeater( ChillerHeaterNum ).ChillerEIRFPLR, PartLoadRat ) );
					CHPower = ( AvailChillerCap / ReferenceCOP ) * nsvChillerEIRFPLR * nsvChillerEIRFT * FRAC;
					ActualCOP = ( QEvaporator + nsvChillerFalseLoadRate ) / CHPower;
					QCondenser = CHPower * this->ChillerHeater( ChillerHeaterNum ).OpenMotorEff + QEvaporator + nsvChillerFalseLoadRate;

					// Determine heating load for this heater and pass the remaining load to the next chiller heater
					CondenserCapMin = QCondenser * MinPartLoadRat;
					HeatingLoadToMeet = min( QCondenser, max( std::abs( CondenserLoad ), CondenserCapMin ) );

					// Set load this chiller heater should meet and temperatures given
					QCondenser = min( HeatingLoadToMeet, QCondenser );

					Cp = GetSpecificHeatGlycol( PlantLoop( this->HWLoopNum ).FluidName, CondInletTemp, PlantLoop( this->HWLoopNum ).FluidIndex, RoutineNameElecEIRChiller );

					// Calculate temperatures for constant flow and mass flow rate for variable flow
					// Limit mass for this chiller heater to the available mass at given temperature conditions
					// when mass calculated to meet the load is greater than the maximum available
					// then recalculate heating load this chiller heater can meet
					if ( CurrentMode == 2 || SimulHtgDominant ) {
						if ( CondMassFlowRate > MassFlowTolerance && CondDeltaTemp > 0.0 ) {
							if ( this->VariableFlowCH ) { // Variable flow
								CondMassFlowRateCalc = QCondenser / CondDeltaTemp / Cp;
								if ( CondMassFlowRateCalc > CondMassFlowRate ) {
									CondMassFlowRateCalc = CondMassFlowRate;
									CondDeltaTempCalc = QCondenser / CondMassFlowRate / Cp;
									if ( CondDeltaTempCalc > CondDeltaTemp ) { // Load to meet should be adjusted
										CondDeltaTempCalc = CondDeltaTemp;
										QCondenser = CondMassFlowRate * Cp * CondDeltaTemp;
									}
								}
								CondMassFlowRate = CondMassFlowRateCalc;
							} else { // Constant Flow at a fixed flow rate and capacity
								CondDeltaTempCalc = QCondenser / CondMassFlowRate / Cp;
								CondOutletTempCalc = CondDeltaTempCalc + CondInletTemp;
								if ( CondOutletTempCalc > CondOutletTemp ) { // Load to meet should be adjusted
									CondOutletTempCalc = CondOutletTemp;
									QCondenser = CondMassFlowRate * Cp * CondDeltaTemp;
								}
								CondOutletTemp = CondOutletTempCalc;
							}
						} else {
							QCondenser = 0.0;
							CondOutletTemp = CondInletTemp;
						}
					}

				} // End of calculaton dependiong on the modes

				// Determine load next chiller heater meets
				if ( CondenserLoad < QCondenser ) { // Heating load is met by this chiller heater
					CondenserLoad = 0.0;
				} else {
					CondenserLoad -= QCondenser;
				}

				if ( QCondenser == 0.0 ) {
					CurrentMode = 0;
					nsvChillerPartLoadRatio = 0.0;
					nsvChillerCyclingRatio = 0.0;
					nsvChillerFalseLoadRate = 0.0;
					EvapMassFlowRate = 0.0;
					CondMassFlowRate = 0.0;
					CHPower = 0.0;
					QEvaporator = 0.0;
					CondenserFanPower = 0.0;
					EvapOutletTemp = EvapInletTemp;
					CondOutletTemp = CondInletTemp;
					CondenserLoad = 0.0;
				}

				// Heat recovery or cooling dominant modes need to use the evaporator side information
				if ( CurrentMode == 3 || CurrentMode == 4 ) {
					nsvChillerPartLoadRatio = this->ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatioSimul;
					nsvChillerCyclingRatio = this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatioSimul;
					nsvChillerFalseLoadRate = this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRateSimul;
					nsvChillerCapFT = this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFTSimul;
					nsvChillerEIRFT = this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFTSimul;
					nsvChillerEIRFPLR = this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLRSimul;
					QEvaporator = this->ChillerHeaterReport( ChillerHeaterNum ).QEvapSimul;
					EvapOutletTemp = this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTempSimul;
					EvapInletTemp = this->ChillerHeaterReport( ChillerHeaterNum ).EvapInletTempSimul;
					EvapMassFlowRate = this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul;
					if ( SimulClgDominant ) {
						CHPower = this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPowerSimul;
						this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower = 0.0;
					}
				}
			}

			// Check if it is mode 4, then skip binding local variables
			if ( CurrentMode == 4 ) {
				this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = CurrentMode;
			} else {
				this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.MassFlowRate = EvapMassFlowRate;
				this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.MassFlowRate = CondMassFlowRate;
				this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.Temp = EvapOutletTemp;
				this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp = EvapInletTemp;
				this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.Temp = CondOutletTemp;
				this->ChillerHeater( ChillerHeaterNum ).CondInletNode.Temp = CondInletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = CurrentMode;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatio = nsvChillerPartLoadRatio;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatio = nsvChillerCyclingRatio;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate = nsvChillerFalseLoadRate;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFT = nsvChillerCapFT;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFT = nsvChillerEIRFT;
				this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLR = nsvChillerEIRFPLR;
				this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPower = CoolingPower;
				this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower = CHPower;
				this->ChillerHeaterReport( ChillerHeaterNum ).QEvap = QEvaporator;
				this->ChillerHeaterReport( ChillerHeaterNum ).QCond = QCondenser;
				this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp = EvapOutletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).EvapInletTemp = EvapInletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp = CondOutletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).CondInletTemp = CondInletTemp;
				this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot = EvapMassFlowRate;
				this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot = CondMassFlowRate;
				this->ChillerHeaterReport( ChillerHeaterNum ).ActualCOP = ActualCOP;
			}

		}

	}

	void
	WrapperSpecs::CalcWrapperModel(
		int const WrapperNum,
		Real64 & MyLoad,
		bool const RunFlag,
		bool const FirstIteration,
		int const EquipFlowCtrl,
		int const LoopNum
	)
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Daeho Kang, PNNL
		//       DATE WRITTEN   Feb 2013
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		//  Calculate node information connected to plnat & condenser loop

		// METHODOLOGY EMPLOYED:
		//  Use empirical curve fits to model performance at off-reference conditions

		// REFERENCES:

		// Using/Aliasing
		using DataGlobals::WarmupFlag;
		using DataHVACGlobals::SmallLoad;
		using CurveManager::CurveValue;
		using DataPlant::DeltaTempTol;
		using DataPlant::TypeOf_CentralGroundSourceHeatPump;
		using DataBranchAirLoopPlant::MassFlowTolerance;
		using PlantUtilities::SetComponentFlowRate;
		using ScheduleManager::GetCurrentScheduleValue;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// LOCAL VARIABLES
		int ChillerHeaterNum; // Chiller heater number
		int CHWInletNodeNum; // Chiller heater bank chilled water inlet node number
		int CHWOutletNodeNum; // Chiller heater bank chilled water Outlet node number
		int GLHEInletNodeNum; // Chiller heater bank condenser water inlet node number
		int GLHEOutletNodeNum; // Chiller heater bank condenser water outlet node number
		int HWInletNodeNum; // Chiller heater bank hot water inlet node number
		int HWOutletNodeNum; // Chiller heater bank hot water outlet node number
		int LoopSideNum; // Loop side number
		int LoopSide; // Loop side
		int OpMode; // Operation mode
		int ChillerHeaterNums; // Total number of chiller heaters
		Real64 CurCoolingLoad; // Total cooling load chiller heater bank (wrapper) meets
		Real64 CurHeatingLoad; // Total heating load chiller heater bank (wrapper) meets
		Real64 CHWInletTemp; // Chiller heater bank chilled water inlet temperature
		Real64 CHWOutletTemp; // Chiller heater bank chilled water outlet temperature
		Real64 CHWInletMassFlowRate; // Chiller heater bank chilled water inlet mass flow rate
		Real64 CHWOutletMassFlowRate; // Chiller heater bank chilled water outlet mass flow rate
		Real64 CHWBypassMassFlowRate; // Chiller heater bank chilled water bypass mass flow rate
		Real64 HWInletTemp; // Chiller heater bank hot water inlet temperature
		Real64 HWOutletTemp; // Chiller heater bank hot water outlet temperature
		Real64 HWInletMassFlowRate; // Chiller heater bank hot water inlet mass flow rate
		Real64 HWOutletMassFlowRate; // Chiller heater bank hot water outlet mass flow rate
		Real64 HWBypassMassFlowRate; // Chiller heater bank hot water bypass mass flow rate
		Real64 GLHEInletTemp; // Chiller heater bank condenser loop inlet temperature
		Real64 GLHEOutletTemp; // Chiller heater bank condenser loop outlet temperature
		Real64 GLHEInletMassFlowRate; // Chiller heater bank condenser loop intlet mass flow rate
		Real64 GLHEOutletMassFlowRate; // Chiller heater bank condenser loop outlet mass flow rate
		Real64 GLHEBypassMassFlowRate; // Chiller heater bank condenser loop bypass mass flow rate
		Real64 WrapperElecPowerCool( 0.0 ); // Chiller heater bank total cooling electricity [W]
		Real64 WrapperElecPowerHeat( 0.0 ); // Chiller heater bank total heating electricity [W]
		Real64 WrapperCoolRate( 0.0 ); // Chiller heater bank total cooling rate [W]
		Real64 WrapperHeatRate( 0.0 ); // Chiller heater bank total heating rate [W]
		Real64 WrapperGLHERate( 0.0 ); // Chiller heater bank total condenser heat transfer rate [W]
		Real64 WrapperElecEnergyCool( 0.0 ); // Chiller heater bank total electric cooling energy [J]
		Real64 WrapperElecEnergyHeat( 0.0 ); // Chiller heater bank total electric heating energy [J]
		Real64 WrapperCoolEnergy( 0.0 ); // Chiller heater bank total cooling energy [J]
		Real64 WrapperHeatEnergy( 0.0 ); // Chiller heater bank total heating energy [J]
		Real64 WrapperGLHEEnergy( 0.0 ); // Chiller heater bank total condenser heat transfer energy [J]
		int CurrentMode; // Current operation mode indicator

		//Autodesk:Uninit Initialize variables used uninitialized
		OpMode = 0; //Autodesk:Uninit Force default initialization: This didn't cause problems because it's not actually used by functions it is passed to

		// Read note information
		CHWInletNodeNum = this->CHWInletNodeNum;
		CHWOutletNodeNum = this->CHWOutletNodeNum;
		HWInletNodeNum = this->HWInletNodeNum;
		HWOutletNodeNum = this->HWOutletNodeNum;
		GLHEInletNodeNum = this->GLHEInletNodeNum;
		GLHEOutletNodeNum = this->GLHEOutletNodeNum;

		CHWInletMassFlowRate = 0.0;
		HWInletMassFlowRate = 0.0;
		GLHEInletMassFlowRate = 0.0;
		CHWInletTemp = Node( CHWInletNodeNum ).Temp;
		HWInletTemp = Node( HWInletNodeNum ).Temp;
		GLHEInletTemp = Node( GLHEInletNodeNum ).Temp;

		ChillerHeaterNums = this->ChillerHeaterNums;

		// Initiate loads and inlet temperatures each loop
		if ( LoopNum == this->CWLoopNum ) {
			CHWInletMassFlowRate = Node( CHWInletNodeNum ).MassFlowRateMaxAvail;
			HWInletMassFlowRate = Node( HWInletNodeNum ).MassFlowRate;
			GLHEInletMassFlowRate = Node( GLHEInletNodeNum ).MassFlowRateMaxAvail;
			LoopSideNum = this->CWLoopSideNum;
			LoopSide = this->CWLoopSideNum;
			this->WrapperCoolingLoad = 0.0;
			CurCoolingLoad = std::abs( MyLoad );
			this->WrapperCoolingLoad = CurCoolingLoad;
			// Set actual mass flow rate at the nodes when it's locked
			if ( PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 1 ) {
				CHWInletMassFlowRate = Node( CHWInletNodeNum ).MassFlowRate;
			}
			if ( CHWInletMassFlowRate == 0.0 ) GLHEInletMassFlowRate = 0.0;

		} else if ( LoopNum == this->HWLoopNum ) {
			CHWInletMassFlowRate = Node( CHWInletNodeNum ).MassFlowRate;
			HWInletMassFlowRate = Node( HWInletNodeNum ).MassFlowRateMaxAvail;
			GLHEInletMassFlowRate = Node( GLHEInletNodeNum ).MassFlowRateMaxAvail;
			LoopSideNum = this->HWLoopSideNum;
			this->WrapperHeatingLoad = 0.0;
			CurHeatingLoad = MyLoad;
			this->WrapperHeatingLoad = CurHeatingLoad;
			// Set actual mass flow rate at the nodes when it's locked
			if ( PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 1 ) {
				HWInletMassFlowRate = Node( HWInletNodeNum ).MassFlowRate;
			}
			if ( HWInletMassFlowRate == 0.0 ) GLHEInletMassFlowRate = 0.0;
		}

		if ( LoopNum == this->CWLoopNum ) {
			if ( this->ControlMode == SmartMixing ) {
				if ( CurCoolingLoad > 0.0 && CHWInletMassFlowRate > 0.0 && GLHEInletMassFlowRate > 0 ) {

					this->CalcChillerModel( OpMode, MyLoad, RunFlag, FirstIteration, EquipFlowCtrl, LoopNum ); //Autodesk:Uninit OpMode was uninitialized
					UpdateChillerRecords( WrapperNum );

					// Initialize local variables only for calculating mass-weighed temperatures
					CHWOutletTemp = 0.0;
					HWOutletTemp = 0.0;
					GLHEOutletTemp = 0.0;
					CHWOutletMassFlowRate = 0.0;
					HWOutletMassFlowRate = 0.0;
					GLHEOutletMassFlowRate = 0.0;

					for ( ChillerHeaterNum = 1; ChillerHeaterNum <= ChillerHeaterNums; ++ChillerHeaterNum ) {

						// Calculated mass flow rate used by individual chiller heater and bypasses
						CHWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot;
						CHWOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp * ( this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot / CHWInletMassFlowRate );
						WrapperElecPowerCool += this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPower;
						WrapperCoolRate += this->ChillerHeaterReport( ChillerHeaterNum ).QEvap;
						WrapperElecEnergyCool += this->ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy;
						WrapperCoolEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy;
						if ( GLHEInletMassFlowRate > 0.0 ) {
							GLHEOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot;
							if ( GLHEOutletMassFlowRate > GLHEInletMassFlowRate ) GLHEOutletMassFlowRate = GLHEInletMassFlowRate;
							GLHEOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp * ( this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot / GLHEInletMassFlowRate );
							WrapperGLHERate += this->ChillerHeaterReport( ChillerHeaterNum ).QCond;
							WrapperGLHEEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).CondEnergy;
						} else {
							GLHEInletMassFlowRate = 0.0;
							GLHEOutletMassFlowRate = 0.0;
							GLHEOutletTemp = GLHEInletTemp;
							WrapperGLHERate = 0.0;
							WrapperGLHEEnergy = 0.0;
						}
					} // End of summation of mass flow rates and mass weighted temperatrue

					// Calculate temperatures for the mixed flows in the chiller bank
					CHWBypassMassFlowRate = CHWInletMassFlowRate - CHWOutletMassFlowRate;
					if ( CHWBypassMassFlowRate > 0.0 ) {
						CHWOutletTemp += CHWInletTemp * CHWBypassMassFlowRate / CHWInletMassFlowRate;
					} else {
						//CHWOutletTemp = CHWOutletTemp; // Self-assignment commented out
					}

					if ( GLHEInletMassFlowRate > 0.0 ) {
						GLHEBypassMassFlowRate = GLHEInletMassFlowRate - GLHEOutletMassFlowRate;
						if ( GLHEBypassMassFlowRate > 0.0 ) {
							GLHEOutletTemp += GLHEInletTemp * GLHEBypassMassFlowRate / GLHEInletMassFlowRate;
						} else {
							//GLHEOutletTemp = GLHEOutletTemp; // Self-assignment commented out
						}
					} else {
						GLHEOutletTemp = GLHEInletTemp;
					}

					HWOutletTemp = HWInletTemp;

					if ( GetCurrentScheduleValue( this->SchedPtr ) > 0 ) {
						WrapperElecPowerCool += ( this->AncilliaryPower * GetCurrentScheduleValue( this->SchedPtr ) );
					}

					Node( CHWOutletNodeNum ).Temp = CHWOutletTemp;
					Node( HWOutletNodeNum ).Temp = HWOutletTemp;
					Node( GLHEOutletNodeNum ).Temp = GLHEOutletTemp;

				} else {

					// Initialize local variables
					CHWOutletTemp = CHWInletTemp;
					HWOutletTemp = HWInletTemp;
					GLHEOutletTemp = GLHEInletTemp;

					for ( ChillerHeaterNum = 1; ChillerHeaterNum <= ChillerHeaterNums; ++ChillerHeaterNum ) {
						this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.MassFlowRate = 0.0;
						this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.MassFlowRate = 0.0;
						this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.Temp = CHWInletTemp;
						this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp = CHWInletTemp;
						this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.Temp = GLHEInletTemp;
						this->ChillerHeater( ChillerHeaterNum ).CondInletNode.Temp = GLHEInletTemp;
						this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = 0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatio = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatio = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFT = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFT = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLR = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPower = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).QEvap = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).QCond = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp = CHWOutletTemp;
						this->ChillerHeaterReport( ChillerHeaterNum ).EvapInletTemp = CHWInletTemp;
						this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp = GLHEOutletTemp;
						this->ChillerHeaterReport( ChillerHeaterNum ).CondInletTemp = GLHEInletTemp;
						this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoad = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).CondEnergy = 0.0;
						this->ChillerHeaterReport( ChillerHeaterNum ).ActualCOP = 0.0;
					}

				}

				if ( SimulHtgDominant || SimulClgDominant ) {
					Node( CHWOutletNodeNum ).Temp = CHWOutletTemp;
					WrapperReport( WrapperNum ).CHWInletTempSimul = CHWInletTemp;
					WrapperReport( WrapperNum ).CHWOutletTempSimul = CHWOutletTemp;
					WrapperReport( WrapperNum ).CHWmdotSimul = CHWInletMassFlowRate;
					WrapperReport( WrapperNum ).GLHEInletTempSimul = GLHEInletTemp;
					WrapperReport( WrapperNum ).GLHEOutletTempSimul = GLHEOutletTemp;
					WrapperReport( WrapperNum ).GLHEmdotSimul = GLHEInletMassFlowRate;
					WrapperReport( WrapperNum ).TotElecCoolingSimul = WrapperElecEnergyCool;
					WrapperReport( WrapperNum ).CoolingEnergySimul = WrapperCoolEnergy;
					WrapperReport( WrapperNum ).TotElecCoolingPwrSimul = WrapperElecPowerCool;
					WrapperReport( WrapperNum ).CoolingRateSimul = WrapperCoolRate;

				} else {

					Node( CHWOutletNodeNum ).Temp = CHWOutletTemp;
					Node( HWOutletNodeNum ).Temp = HWOutletTemp;
					Node( GLHEOutletNodeNum ).Temp = GLHEOutletTemp;
					WrapperReport( WrapperNum ).CHWInletTemp = CHWInletTemp;
					WrapperReport( WrapperNum ).CHWOutletTemp = CHWOutletTemp;
					WrapperReport( WrapperNum ).HWInletTemp = HWInletTemp;
					WrapperReport( WrapperNum ).HWOutletTemp = HWOutletTemp;
					WrapperReport( WrapperNum ).GLHEInletTemp = GLHEInletTemp;
					WrapperReport( WrapperNum ).GLHEOutletTemp = GLHEOutletTemp;
					WrapperReport( WrapperNum ).CHWmdot = CHWInletMassFlowRate;
					WrapperReport( WrapperNum ).HWmdot = HWInletMassFlowRate;
					WrapperReport( WrapperNum ).GLHEmdot = GLHEInletMassFlowRate;
					WrapperReport( WrapperNum ).TotElecCooling = WrapperElecEnergyCool;
					WrapperReport( WrapperNum ).TotElecHeating = WrapperElecEnergyHeat;
					WrapperReport( WrapperNum ).CoolingEnergy = WrapperCoolEnergy;
					WrapperReport( WrapperNum ).HeatingEnergy = WrapperHeatEnergy;
					WrapperReport( WrapperNum ).GLHEEnergy = WrapperGLHEEnergy;
					WrapperReport( WrapperNum ).TotElecCoolingPwr = WrapperElecPowerCool;
					WrapperReport( WrapperNum ).TotElecHeatingPwr = WrapperElecPowerHeat;
					WrapperReport( WrapperNum ).CoolingRate = WrapperCoolRate;
					WrapperReport( WrapperNum ).HeatingRate = WrapperHeatRate;
					WrapperReport( WrapperNum ).GLHERate = WrapperGLHERate;

				}
				SetComponentFlowRate( CHWInletMassFlowRate, CHWInletNodeNum, CHWOutletNodeNum, this->CWLoopNum, this->CWLoopSideNum, this->CWBranchNum, this->CWCompNum );

				SetComponentFlowRate( HWInletMassFlowRate, HWInletNodeNum, HWOutletNodeNum, this->HWLoopNum, this->HWLoopSideNum, this->HWBranchNum, this->HWCompNum );

				SetComponentFlowRate( GLHEInletMassFlowRate, GLHEInletNodeNum, GLHEOutletNodeNum, this->GLHELoopNum, this->GLHELoopSideNum, this->GLHEBranchNum, this->GLHECompNum );

			} // End of cooling

		} else if ( LoopNum == this->HWLoopNum ) { // Hot water loop
			if ( this->ControlMode == SmartMixing ) { // Chiller heater component
				if ( CurHeatingLoad > 0.0 && HWInletMassFlowRate > 0.0 ) {

					this->CalcChillerHeaterModel( OpMode, MyLoad, RunFlag, FirstIteration, EquipFlowCtrl, LoopNum ); //Autodesk:Uninit OpMode was uninitialized
					UpdateChillerHeaterRecords( WrapperNum );

					// Calculate individual CH units's temperatures and mass flow rates
					CHWOutletTemp = 0.0;
					HWOutletTemp = 0.0;
					GLHEOutletTemp = 0.0;
					CHWOutletMassFlowRate = 0.0;
					HWOutletMassFlowRate = 0.0;
					GLHEOutletMassFlowRate = 0.0;

					if ( SimulHtgDominant || SimulClgDominant ) {
						if ( SimulClgDominant ) {
							for ( ChillerHeaterNum = 1; ChillerHeaterNum <= ChillerHeaterNums; ++ChillerHeaterNum ) {
								CurrentMode = this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode;
								CHWInletTemp = WrapperReport( WrapperNum ).CHWInletTempSimul;
								GLHEInletTemp = WrapperReport( WrapperNum ).GLHEInletTempSimul;
								CHWInletMassFlowRate = WrapperReport( WrapperNum ).CHWmdotSimul;
								GLHEInletMassFlowRate = WrapperReport( WrapperNum ).GLHEmdotSimul;

								if ( CurrentMode != 0 ) { // This chiller heater unit is on
									if ( CurrentMode == 3 ) { // Heat recovery mode. Both chilled water and hot water connections
										CHWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul; // Wrapper evaporator side to plant chilled water loop
										HWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot; // Wrapper condenser side to plant hot water loop
										if ( HWInletMassFlowRate > 0.0 ) {
											HWOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp * ( this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot / HWInletMassFlowRate ); // Only calculate in the heat recovery mode
										} else {
											HWOutletTemp = HWInletTemp;
										}
									} else { // Mode 4. Cooling-only mode with other heat recovery units. Condenser flows.
										CHWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul; // Wrapper evaporator side to plant chilled water loop
										// Sum condenser node mass flow rates and mass weighed temperatures
										if ( GLHEInletMassFlowRate > 0.0 ) {
											GLHEOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).CondmdotSimul;
											if ( GLHEOutletMassFlowRate > GLHEInletMassFlowRate ) GLHEOutletMassFlowRate = GLHEInletMassFlowRate;
											GLHEOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTempSimul * ( this->ChillerHeaterReport( ChillerHeaterNum ).CondmdotSimul / GLHEInletMassFlowRate );
											WrapperGLHERate += this->ChillerHeaterReport( ChillerHeaterNum ).QCondSimul;
											WrapperGLHEEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).CondEnergySimul;
										} else {
											GLHEInletMassFlowRate = 0.0;
											GLHEOutletMassFlowRate = 0.0;
											GLHEOutletTemp = GLHEInletTemp;
											WrapperGLHERate = 0.0;
											WrapperGLHEEnergy = 0.0;
										}
									}
								} else { // This chiller heater is off
									// Check if any unit is cooling only mode
									if ( ChillerHeaterNum == ChillerHeaterNums ) { // All units are heat revocery mode. No condenser flow
										GLHEOutletMassFlowRate = 0.0;
										GLHEInletMassFlowRate = 0.0;
										GLHEOutletTemp = GLHEInletTemp;
									} else { // At leaset, one of chiller heater units is cooling-only mode
										//GLHEOutletMassFlowRate = GLHEOutletMassFlowRate; // Self-assignment commented out
										//GLHEOutletTemp = GLHEOutletTemp; // Self-assignment commented out
									}
								}
								// Calculate mass weighed chilled water temperatures
								if ( CHWInletMassFlowRate > 0.0 ) {
									CHWOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTempSimul * ( this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul / CHWInletMassFlowRate );
								} else {
									CHWOutletTemp = CHWInletTemp;
								}

								WrapperElecPowerCool += this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPowerSimul; // Cooling electricity
								WrapperCoolRate += this->ChillerHeaterReport( ChillerHeaterNum ).QEvapSimul;
								WrapperElecEnergyCool += this->ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergySimul;
								WrapperCoolEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).EvapEnergySimul;
								// Avoid double counting wrapper energy use
								WrapperElecPowerHeat = 0.0;
								WrapperHeatRate = 0.0;
								WrapperElecEnergyHeat = 0.0;
								WrapperHeatEnergy = 0.0;

							}

							// Calculate chilled water temperature
							if ( CHWInletMassFlowRate > 0.0 ) {
								CHWBypassMassFlowRate = CHWInletMassFlowRate - CHWOutletMassFlowRate;
								if ( CHWBypassMassFlowRate > 0.0 ) {
									CHWOutletTemp += CHWInletTemp * CHWBypassMassFlowRate / CHWInletMassFlowRate;
								} else { // No bypass withnin a wrapper
									//CHWOutletTemp = CHWOutletTemp; // Self-assignment commented out
								}
							} else {
								CHWOutletTemp = CHWInletTemp;
							}
							// Calculate hot water outlet temperature
							if ( HWInletMassFlowRate > 0.0 ) {
								HWBypassMassFlowRate = HWInletMassFlowRate - HWOutletMassFlowRate;
								if ( HWBypassMassFlowRate > 0.0 ) {
									HWOutletTemp += HWInletTemp * HWBypassMassFlowRate / HWInletMassFlowRate;
								} else {
									//HWOutletTemp = HWOutletTemp; // Self-assignment commented out
								}
							} else {
								HWOutletTemp = HWInletTemp;
							}
							// Calculate condenser outlet temperature
							if ( GLHEInletMassFlowRate > 0.0 ) {
								GLHEBypassMassFlowRate = GLHEInletMassFlowRate - GLHEOutletMassFlowRate;
								if ( GLHEBypassMassFlowRate > 0.0 ) {
									GLHEOutletTemp += GLHEInletTemp * GLHEBypassMassFlowRate / GLHEInletMassFlowRate;
								} else {
									//GLHEOutletTemp = GLHEOutletTemp; // Self-assignment commented out
								}
							} else {
								GLHEOutletTemp = GLHEInletTemp;
							}

							// Add ancilliary power if scheduled
							if ( GetCurrentScheduleValue( this->SchedPtr ) > 0 ) {
								WrapperElecPowerCool += ( this->AncilliaryPower * GetCurrentScheduleValue( this->SchedPtr ) );
							}

							// Electricity should be counted once for cooling in this mode
							WrapperElecEnergyHeat = 0.0;

						} else if ( SimulHtgDominant ) { // Heating dominant simultaneous clg/htg mode

							for ( ChillerHeaterNum = 1; ChillerHeaterNum <= ChillerHeaterNums; ++ChillerHeaterNum ) {
								// Set temperatures and mass flow rates for the cooling side
								CurrentMode = this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode;
								CHWInletTemp = WrapperReport( WrapperNum ).CHWInletTempSimul;
								CHWInletMassFlowRate = WrapperReport( WrapperNum ).CHWmdotSimul;

								if ( CurrentMode != 0 ) { // This chiller heater unit is on
									if ( CurrentMode == 3 ) { // Heat recovery mode. Both chilled water and hot water connections
										CHWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul; // Wrapper evaporator side to plant chilled water loop
										HWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot; // Wrapper condenser side to plant hot water loop
										if ( CHWInletMassFlowRate > 0.0 ) {
											CHWOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTempSimul * ( this->ChillerHeaterReport( ChillerHeaterNum ).EvapmdotSimul / CHWInletMassFlowRate ); // Only need to calculate in the heat recovery mode
										} else {
											CHWOutletTemp = CHWInletTemp;
										}
									} else { // Mode 5. Heating only mode with other heat recovery units
										HWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot; // Wrapper condenser side to plant hot water loop
										if ( GLHEInletMassFlowRate > 0.0 ) {
											GLHEOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot; // Wrapper evaporator side to plant condenser loop
											if ( GLHEOutletMassFlowRate > GLHEInletMassFlowRate ) GLHEOutletMassFlowRate = GLHEInletMassFlowRate;
											GLHEOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp * ( this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot / GLHEInletMassFlowRate );
											WrapperGLHERate += this->ChillerHeaterReport( ChillerHeaterNum ).QEvap;
											WrapperGLHEEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy;
										} else {
											GLHEInletMassFlowRate = 0.0;
											GLHEOutletMassFlowRate = 0.0;
											GLHEOutletTemp = GLHEInletTemp;
											WrapperGLHERate = 0.0;
											WrapperGLHEEnergy = 0.0;
										}
									} // End of heat recovery mode

								} else { // This chiller heater is off

									// Check if any unit is heating only mode
									if ( ChillerHeaterNum == ChillerHeaterNums ) { // All are heat revocery mode. No condenser flow
										GLHEOutletMassFlowRate = 0.0;
										GLHEInletMassFlowRate = 0.0;
										GLHEOutletTemp = GLHEInletTemp;
									} else { // At leaset, one of chiller heater units is heating only mode
										//GLHEOutletMassFlowRate = GLHEOutletMassFlowRate; // Self-assignment commented out
										//GLHEOutletTemp = GLHEOutletTemp; // Self-assignment commented out
									}
								}

								// Calculate mass weighed hot water temperatures
								if ( HWInletMassFlowRate > 0.0 ) {
									HWOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp * ( this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot / HWInletMassFlowRate ); // Always heating as long as heating load remains
								} else {
									HWOutletTemp = HWInletTemp;
								}

								WrapperElecPowerHeat += this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower;
								WrapperHeatRate += this->ChillerHeaterReport( ChillerHeaterNum ).QCond;
								WrapperElecEnergyHeat += this->ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy;
								WrapperHeatEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).CondEnergy;

								// Avoid double counting wrapper energy use
								WrapperElecPowerCool = 0.0;
								WrapperCoolRate = 0.0;
								WrapperElecEnergyCool = 0.0;
							}
							// Calculate chilled water outlet temperature
							if ( CHWInletMassFlowRate > 0.0 ) {
								CHWBypassMassFlowRate = CHWInletMassFlowRate - CHWOutletMassFlowRate;
								if ( CHWBypassMassFlowRate > 0.0 ) {
									CHWOutletTemp += CHWInletTemp * CHWBypassMassFlowRate / CHWInletMassFlowRate;
								} else { // No bypass withnin a wrapper
									//CHWOutletTemp = CHWOutletTemp; // Self-assignment commented out
								}
							} else {
								CHWOutletTemp = CHWInletTemp;
							}
							// Calculate hot water outlet temperature
							if ( HWInletMassFlowRate > 0.0 ) {
								HWBypassMassFlowRate = HWInletMassFlowRate - HWOutletMassFlowRate;
								if ( HWBypassMassFlowRate > 0.0 ) {
									HWOutletTemp += HWInletTemp * HWBypassMassFlowRate / HWInletMassFlowRate;
								} else {
									//HWOutletTemp = HWOutletTemp; // Self-assignment commented out
								}
							} else {
								HWOutletTemp = HWInletTemp;
							}
							// Calculate condenser outlet temperature
							if ( GLHEInletMassFlowRate > 0.0 ) {
								GLHEBypassMassFlowRate = GLHEInletMassFlowRate - GLHEOutletMassFlowRate;
								if ( GLHEBypassMassFlowRate > 0.0 ) {
									GLHEOutletTemp += GLHEInletTemp * GLHEBypassMassFlowRate / GLHEInletMassFlowRate;
								} else {
									//GLHEOutletTemp = GLHEOutletTemp; // Self-assignment commented out
								}
							} else {
								GLHEOutletTemp = GLHEInletTemp;
							}

							// Check if ancilliary power is used
							if ( GetCurrentScheduleValue( this->SchedPtr ) > 0 ) {
								WrapperElecPowerHeat += ( this->AncilliaryPower * GetCurrentScheduleValue( this->SchedPtr ) );
							}

							// Electricity should be counted once
							WrapperElecEnergyCool = 0.0;

						} // End of simultaneous clg/htg mode calculations

					} else { // Heating only mode (mode 2)

						for ( ChillerHeaterNum = 1; ChillerHeaterNum <= ChillerHeaterNums; ++ChillerHeaterNum ) {
							HWOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot;
							HWOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp * this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot / HWInletMassFlowRate;
							WrapperElecPowerHeat += this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower;
							WrapperHeatRate += this->ChillerHeaterReport( ChillerHeaterNum ).QCond;
							WrapperElecEnergyHeat += this->ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy;
							WrapperHeatEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).CondEnergy;

							if ( GLHEInletMassFlowRate > 0.0 ) {
								GLHEOutletMassFlowRate += this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot;
								if ( GLHEOutletMassFlowRate > GLHEInletMassFlowRate ) GLHEOutletMassFlowRate = GLHEInletMassFlowRate;
								GLHEOutletTemp += this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp * ( this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot / GLHEInletMassFlowRate );
								WrapperGLHERate += this->ChillerHeaterReport( ChillerHeaterNum ).QEvap;
								WrapperGLHEEnergy += this->ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy;
							} else { // No source water flow
								GLHEOutletMassFlowRate = 0.0;
								GLHEInletMassFlowRate = 0.0;
								GLHEOutletTemp = GLHEInletTemp;
								WrapperGLHERate = 0.0;
								WrapperGLHEEnergy = 0.0;
							}
						}

						// Calculate hot water outlet temperature
						if ( HWInletMassFlowRate > 0.0 ) {
							HWBypassMassFlowRate = HWInletMassFlowRate - HWOutletMassFlowRate;
							if ( HWBypassMassFlowRate > 0.0 ) {
								HWOutletTemp += HWInletTemp * HWBypassMassFlowRate / HWInletMassFlowRate;
							} else {
								//HWOutletTemp = HWOutletTemp; // Self-assignment commented out
								if ( HWOutletTemp > HWInletTemp ) HWOutletTemp = HWInletTemp;
							}
						} else {
							HWOutletTemp = HWInletTemp;
						}

						// Calculate condenser outlet temperature
						if ( GLHEInletMassFlowRate > 0.0 ) {
							GLHEBypassMassFlowRate = GLHEInletMassFlowRate - GLHEOutletMassFlowRate;
							if ( GLHEBypassMassFlowRate > 0.0 ) {
								GLHEOutletTemp += GLHEInletTemp * GLHEBypassMassFlowRate / GLHEInletMassFlowRate;
							} else {
								//GLHEOutletTemp = GLHEOutletTemp; // Self-assignment commented out
							}
						} else {
							GLHEOutletTemp = GLHEInletTemp;
						}

						CHWOutletTemp = CHWInletTemp;

						// Add ancilliary power if necessary
						if ( GetCurrentScheduleValue( this->SchedPtr ) > 0 ) {
							WrapperElecPowerHeat += ( this->AncilliaryPower * GetCurrentScheduleValue( this->SchedPtr ) );
						}

					} // End of calculations

					SetComponentFlowRate( CHWInletMassFlowRate, CHWInletNodeNum, CHWOutletNodeNum, this->CWLoopNum, this->CWLoopSideNum, this->CWBranchNum, this->CWCompNum );

					SetComponentFlowRate( HWInletMassFlowRate, HWInletNodeNum, HWOutletNodeNum, this->HWLoopNum, this->HWLoopSideNum, this->HWBranchNum, this->HWCompNum );

					SetComponentFlowRate( GLHEInletMassFlowRate, GLHEInletNodeNum, GLHEOutletNodeNum, this->GLHELoopNum, this->GLHELoopSideNum, this->GLHEBranchNum, this->GLHECompNum );

					// Local variables
					WrapperReport( WrapperNum ).CHWInletTemp = CHWInletTemp;
					WrapperReport( WrapperNum ).CHWOutletTemp = CHWOutletTemp;
					WrapperReport( WrapperNum ).HWInletTemp = HWInletTemp;
					WrapperReport( WrapperNum ).HWOutletTemp = HWOutletTemp;
					WrapperReport( WrapperNum ).GLHEInletTemp = GLHEInletTemp;
					WrapperReport( WrapperNum ).GLHEOutletTemp = GLHEOutletTemp;
					WrapperReport( WrapperNum ).CHWmdot = CHWInletMassFlowRate;
					WrapperReport( WrapperNum ).HWmdot = HWInletMassFlowRate;
					WrapperReport( WrapperNum ).GLHEmdot = GLHEInletMassFlowRate;
					WrapperReport( WrapperNum ).TotElecCooling = WrapperElecEnergyCool;
					WrapperReport( WrapperNum ).TotElecHeating = WrapperElecEnergyHeat;
					WrapperReport( WrapperNum ).CoolingEnergy = WrapperCoolEnergy;
					WrapperReport( WrapperNum ).HeatingEnergy = WrapperHeatEnergy;
					WrapperReport( WrapperNum ).GLHEEnergy = WrapperGLHEEnergy;
					WrapperReport( WrapperNum ).TotElecCoolingPwr = WrapperElecPowerCool;
					WrapperReport( WrapperNum ).TotElecHeatingPwr = WrapperElecPowerHeat;
					WrapperReport( WrapperNum ).CoolingRate = WrapperCoolRate;
					WrapperReport( WrapperNum ).HeatingRate = WrapperHeatRate;
					WrapperReport( WrapperNum ).GLHERate = WrapperGLHERate;

					Node( CHWOutletNodeNum ).Temp = CHWOutletTemp;
					Node( HWOutletNodeNum ).Temp = HWOutletTemp;
					Node( GLHEOutletNodeNum ).Temp = GLHEOutletTemp;

				} else { // Central chiller heater system is off

					CHWOutletTemp = CHWInletTemp;
					HWOutletTemp = HWInletTemp;
					GLHEOutletTemp = GLHEInletTemp;
					Node( CHWOutletNodeNum ).Temp = CHWOutletTemp;
					Node( HWOutletNodeNum ).Temp = HWOutletTemp;
					Node( GLHEOutletNodeNum ).Temp = GLHEOutletTemp;

					if ( this->WrapperCoolingLoad == 0.0 && ! SimulHtgDominant ) {

						for ( ChillerHeaterNum = 1; ChillerHeaterNum <= ChillerHeaterNums; ++ChillerHeaterNum ) {
							this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.MassFlowRate = 0.0;
							this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.MassFlowRate = 0.0;
							this->ChillerHeater( ChillerHeaterNum ).EvapOutletNode.Temp = CHWInletTemp;
							this->ChillerHeater( ChillerHeaterNum ).EvapInletNode.Temp = CHWInletTemp;
							this->ChillerHeater( ChillerHeaterNum ).CondOutletNode.Temp = GLHEInletTemp;
							this->ChillerHeater( ChillerHeaterNum ).CondInletNode.Temp = GLHEInletTemp;
							this->ChillerHeaterReport( ChillerHeaterNum ).CurrentMode = 0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerPartLoadRatio = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCyclingRatio = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerCapFT = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFT = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerEIRFPLR = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).CoolingPower = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).HeatingPower = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).QEvap = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).QCond = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).EvapOutletTemp = CHWOutletTemp;
							this->ChillerHeaterReport( ChillerHeaterNum ).EvapInletTemp = CHWInletTemp;
							this->ChillerHeaterReport( ChillerHeaterNum ).CondOutletTemp = GLHEOutletTemp;
							this->ChillerHeaterReport( ChillerHeaterNum ).CondInletTemp = GLHEInletTemp;
							this->ChillerHeaterReport( ChillerHeaterNum ).Evapmdot = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).Condmdot = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoad = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).CondEnergy = 0.0;
							this->ChillerHeaterReport( ChillerHeaterNum ).ActualCOP = 0.0;
						}

						WrapperReport( WrapperNum ).CHWInletTemp = CHWInletTemp;
						WrapperReport( WrapperNum ).CHWOutletTemp = CHWOutletTemp;
						WrapperReport( WrapperNum ).HWInletTemp = HWInletTemp;
						WrapperReport( WrapperNum ).HWOutletTemp = HWOutletTemp;
						WrapperReport( WrapperNum ).GLHEInletTemp = GLHEInletTemp;
						WrapperReport( WrapperNum ).GLHEOutletTemp = GLHEOutletTemp;
						WrapperReport( WrapperNum ).CHWmdot = CHWInletMassFlowRate;
						WrapperReport( WrapperNum ).HWmdot = HWInletMassFlowRate;
						WrapperReport( WrapperNum ).GLHEmdot = GLHEInletMassFlowRate;
						WrapperReport( WrapperNum ).TotElecCooling = WrapperElecEnergyCool;
						WrapperReport( WrapperNum ).TotElecHeating = WrapperElecEnergyHeat;
						WrapperReport( WrapperNum ).CoolingEnergy = WrapperCoolEnergy;
						WrapperReport( WrapperNum ).HeatingEnergy = WrapperHeatEnergy;
						WrapperReport( WrapperNum ).GLHEEnergy = WrapperGLHEEnergy;
						WrapperReport( WrapperNum ).TotElecCoolingPwr = WrapperElecPowerCool;
						WrapperReport( WrapperNum ).TotElecHeatingPwr = WrapperElecPowerHeat;
						WrapperReport( WrapperNum ).CoolingRate = WrapperCoolRate;
						WrapperReport( WrapperNum ).HeatingRate = WrapperHeatRate;
						WrapperReport( WrapperNum ).GLHERate = WrapperGLHERate;

						SetComponentFlowRate( CHWInletMassFlowRate, CHWInletNodeNum, CHWOutletNodeNum, this->CWLoopNum, this->CWLoopSideNum, this->CWBranchNum, this->CWCompNum );

						SetComponentFlowRate( HWInletMassFlowRate, HWInletNodeNum, HWOutletNodeNum, this->HWLoopNum, this->HWLoopSideNum, this->HWBranchNum, this->HWCompNum );

						SetComponentFlowRate( GLHEInletMassFlowRate, GLHEInletNodeNum, GLHEOutletNodeNum, this->GLHELoopNum, this->GLHELoopSideNum, this->GLHEBranchNum, this->GLHECompNum );
					}

				} // Heating loop calculation

				//!Node(CHWOutletNodeNum)%Temp  = CHWOutletTemp
				//!Node(HWOutletNodeNum)%Temp   = HWOutletTemp
				//!Node(GLHEOutletNodeNum)%Temp = GLHEOutletTemp

			}

		}

	}

	void
	UpdateChillerRecords( int const WrapperNum ) // Wrapper number
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Daeho Kang, PNNL
		//       DATE WRITTEN:    Feb 2013

		// PURPOSE OF THIS SUBROUTINE:
		//  Update cihller heater variables

		// METHODOLOGY EMPLOYED:
		//  na

		// REFERENCES:
		//  na

		// Using/Aliasing
		using DataGlobals::SecInHour;
		using DataHVACGlobals::TimeStepSys;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		// na

		// DERIVED TYPE DEFINITIONS:
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 SecInTimeStep; // Number of seconds per HVAC system time step, to convert from W (J/s) to J
		int ChillerHeaterNum; // Chiller heater number

		SecInTimeStep = TimeStepSys * SecInHour;

		for ( ChillerHeaterNum = 1; ChillerHeaterNum <= Wrapper( WrapperNum ).ChillerHeaterNums; ++ChillerHeaterNum ) {
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoad = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingPower * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).HeatingPower * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).QEvap * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).QCond * SecInTimeStep;
			if ( SimulClgDominant || SimulHtgDominant ) {
				Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadSimul = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoad;
				Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergySimul = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy;
				Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapEnergySimul = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy;
				Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondEnergySimul = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondEnergy;
			}
		}

	}

	void
	UpdateChillerHeaterRecords( int const WrapperNum ) // Wrapper number
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Daeho Kang, PNNL
		//       DATE WRITTEN:    Feb 2013

		// PURPOSE OF THIS SUBROUTINE:
		//  Reporting

		// METHODOLOGY EMPLOYED:
		//  na

		// REFERENCES:
		//  na

		// Using/Aliasing
		using DataGlobals::SecInHour;
		using DataHVACGlobals::TimeStepSys;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		// na

		// DERIVED TYPE DEFINITIONS:
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 SecInTimeStep; // Number of seconds per HVAC system time step, to convert from W (J/s) to J
		int ChillerHeaterNum; // Chiller heater number

		SecInTimeStep = TimeStepSys * SecInHour;

		for ( ChillerHeaterNum = 1; ChillerHeaterNum <= Wrapper( WrapperNum ).ChillerHeaterNums; ++ChillerHeaterNum ) {
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoad = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).ChillerFalseLoadRate * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CoolingPower * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).HeatingEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).HeatingPower * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).EvapEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).QEvap * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondEnergy = Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).QCond * SecInTimeStep;
			Wrapper( WrapperNum ).ChillerHeaterReport( ChillerHeaterNum ).CondenserFanEnergyConsumption = CondenserFanPower * SecInTimeStep;
		}

	}

} // PlantCentralGSHP

} // EnergyPlus
