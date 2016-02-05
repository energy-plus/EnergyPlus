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
#include <cassert>
#include <cmath>

// ObjexxFCL Headers
#include <ObjexxFCL/Array.functions.hh>
#include <ObjexxFCL/gio.hh>

// EnergyPlus Headers
#include <EvaporativeFluidCoolers.hh>
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
#include <DataWater.hh>
#include <FluidProperties.hh>
#include <General.hh>
#include <InputProcessor.hh>
#include <NodeInputManager.hh>
#include <OutAirNodeManager.hh>
#include <OutputProcessor.hh>
#include <OutputReportPredefined.hh>
#include <PlantUtilities.hh>
#include <Psychrometrics.hh>
#include <ReportSizingManager.hh>
#include <ScheduleManager.hh>
#include <UtilityRoutines.hh>
#include <WaterManager.hh>

namespace EnergyPlus {

	// Module containing the routines dealing with the objects EvaporativeFluidCooler:SingleSpeed and
	// EvaporativeFluidCooler:TwoSpeed

	// MODULE INFORMATION:
	//       AUTHOR         Chandan Sharma
	//       DATE WRITTEN   May 2009
	//       MODIFIED       na
	//       RE-ENGINEERED  na

	// PURPOSE OF THIS MODULE:
	// Model the performance of evaporative fluid coolers

	// METHODOLOGY EMPLOYED:
	// Based on cooling tower by Shirey, Raustad: Dec 2000; Shirey, Sept 2002

	using namespace DataPrecisionGlobals;
	using DataGlobals::KelvinConv;
	using DataGlobals::SecInHour;
	using DataGlobals::InitConvTemp;
	using namespace DataHVACGlobals;
	using namespace DataLoopNode;
	using DataEnvironment::StdBaroPress;
	using DataEnvironment::OutDryBulbTemp;
	using General::TrimSigDigits;
	using DataPlant::PlantLoop;
	using DataBranchAirLoopPlant::MassFlowTolerance;
	using Psychrometrics::PsyWFnTdbTwbPb;
	using Psychrometrics::PsyRhoAirFnPbTdbW;
	using Psychrometrics::PsyHFnTdbRhPb;
	using FluidProperties::GetDensityGlycol;
	using FluidProperties::GetSpecificHeatGlycol;

	// Data
	// MODULE PARAMETER DEFINITIONS

	std::string const EvaporativeFluidCooler::cEvapFluidCooler_SingleSpeed( "EvaporativeFluidCooler:SingleSpeed" );
	std::string const EvaporativeFluidCooler::cEvapFluidCooler_TwoSpeed( "EvaporativeFluidCooler:TwoSpeed" );
	std::vector< std::unique_ptr< EvaporativeFluidCooler > > EvaporativeFluidCooler::instances;
	bool EvaporativeFluidCooler::getInputFlag = true;

	void EvaporativeFluidCooler::clear_state()
	{
		EvaporativeFluidCooler::instances.clear();
		EvaporativeFluidCooler::getInputFlag = true;
	}

	PlantComponent * EvaporativeFluidCooler::factory( int objectType, std::string objectName )
	{
		// Process the input data for pipes if it hasn't been done already
		if ( getInputFlag ) {
			EvaporativeFluidCooler::getInput();
			getInputFlag = false;
		}
		// Now look for this particular pipe in the list
		for ( auto const & evaporativeFluidCooler : instances ) {
			if ( evaporativeFluidCooler->PlantType_Num == objectType && evaporativeFluidCooler->Name == objectName ) {
				return evaporativeFluidCooler.get();
			}
		}
		// If we didn't find it, fatal
		ShowFatalError( "EvaporativeFluidCoolerFactory: Error getting inputs for evaporative fluid cooler named: " + objectName );
		// Shut up the compiler
		return nullptr;
	}

	void EvaporativeFluidCooler::onInitLoopEquip( PlantLocation const & EP_UNUSED( calledFromLocation ) )
	{
		this->init();
		this->size();
	}

	void EvaporativeFluidCooler::getDesignCapacities( PlantLocation const & EP_UNUSED( calledFromLocation ), Real64 & MaxLoad, Real64 & MinLoad, Real64 & OptLoad )
	{
		MinLoad = 0.0;
		if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
			MaxLoad = this->HighSpeedStandardDesignCapacity * this->HeatRejectCapNomCapSizingRatio;
			OptLoad = this->HighSpeedStandardDesignCapacity;
		}
		else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
			MaxLoad = 0.0;
			OptLoad = 0.0;
		}
	}

	void EvaporativeFluidCooler::getSizingFactor( Real64 & SizingFactor )
	{
		SizingFactor = this->SizFac;
	}

	void EvaporativeFluidCooler::simulate(
			PlantLocation const & EP_UNUSED( calledFromLocation ),
			bool const EP_UNUSED( FirstHVACIteration ),
			Real64 & CurLoad
	)
	{
		if ( ! ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ||
				this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) ) {
			ShowFatalError( "SimEvapFluidCoolers: Invalid evaporative fluid cooler Type Requested = " + this->EvapFluidCoolerType );
		}

		this->init();
		if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
			this->calcSingleSpeedEvapFluidCooler();
		}
		if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
			this->calcTwoSpeedEvapFluidCooler();
		}
		calculateWaterUseage();
		update();

		bool runFlag = CurLoad > DataPlant::LoopDemandTol || CurLoad < ( -DataPlant::LoopDemandTol );
		report( runFlag );
	}

	void
	EvaporativeFluidCooler::getInput()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    May 2009
		//       MODIFIED         Chandan Sharma, April 2010
		//       RE-ENGINEERED    na

		// PURPOSE OF THIS SUBROUTINE:
		// Obtains input data for evaporative fluid coolers and stores it in SimpleEvapFluidCooler data structure.

		// METHODOLOGY EMPLOYED:
		// Uses "Get" routines to read in the data.

		// REFERENCES:
		// Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002
		// B.A. Qureshi and S.M. Zubair , Prediction of evaporation losses in evaporative fluid coolers
		// Applied thermal engineering 27 (2007) 520-527

		// Using/Aliasing
		using namespace DataSizing;
		using namespace DataLoopNode;
		//  USE DataPlant,          ONLY: PlantLoop
		using InputProcessor::GetNumObjectsFound;
		using InputProcessor::GetObjectItem;
		using InputProcessor::VerifyName;
		using InputProcessor::SameString;
		using InputProcessor::MakeUPPERCase;
		using namespace DataIPShortCuts; // Data for field names, blank numerics
		using NodeInputManager::GetOnlySingleNode;
		using BranchNodeConnections::TestCompSet;
		using CurveManager::GetCurveIndex;
		using ScheduleManager::GetScheduleIndex;
		using WaterManager::SetupTankDemandComponent;
		using OutAirNodeManager::CheckOutAirNodeNumber;
		using General::TrimSigDigits;
		using FluidProperties::CheckFluidPropertyName;
		using FluidProperties::FindGlycol;
		using FluidProperties::GetGlycolNameByIndex;
		using DataEnvironment::OutDryBulbTemp;
		using DataEnvironment::OutRelHumValue;

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int EvapFluidCoolerNum; // Evaporative fluid cooler number,
		// reference counter for SimpleEvapFluidCooler data array
		int NumSingleSpeedEvapFluidCoolers; // Total number of single-speed evaporative fluid coolers
		int SingleSpeedEvapFluidCoolerNumber; // Specific single-speed evaporative fluid cooler of interest
		int NumTwoSpeedEvapFluidCoolers; // Number of two-speed evaporative fluid coolers
		int TwoSpeedEvapFluidCoolerNumber; // Specific two-speed evaporative fluid cooler of interest
		int NumAlphas; // Number of elements in the alpha array
		int NumNums; // Number of elements in the numeric array
		int IOStat; // IO Status when calling get input subroutine
		bool IsNotOK; // Flag to verify name
		bool IsBlank; // Flag for blank name
		bool ErrorsFound( false ); // Logical flag set .TRUE. if errors found while getting input data
		Array1D< Real64 > NumArray( 25 ); // Numeric input data array
		Array1D_string AlphArray( 13 ); // Character string input data array
		std::string FluidName;

		// Get number of all evaporative fluid coolers specified in the input data file (idf)
		NumSingleSpeedEvapFluidCoolers = GetNumObjectsFound( cEvapFluidCooler_SingleSpeed );
		NumTwoSpeedEvapFluidCoolers = GetNumObjectsFound( cEvapFluidCooler_TwoSpeed );
		int NumSimpleEvapFluidCoolers = NumSingleSpeedEvapFluidCoolers + NumTwoSpeedEvapFluidCoolers;

		if ( NumSimpleEvapFluidCoolers <= 0 ) ShowFatalError( "No evaporative fluid cooler objects found in input, however, a branch object has specified an evaporative fluid cooler. Search the input for evaporative fluid cooler to determine the cause for this error." );

		// See if load distribution manager has already gotten the input
		if ( instances.size() > 0 ) return;

		// Allocate data structures to hold evaporative fluid cooler input data,
		// report data and evaporative fluid cooler inlet conditions
		for ( int i = 0; i < NumSimpleEvapFluidCoolers; ++i ) {
			instances.emplace_back( new EvaporativeFluidCooler() );
		}

		// Load data structures with evaporative fluid cooler input data
		cCurrentModuleObject = cEvapFluidCooler_SingleSpeed;
		for ( SingleSpeedEvapFluidCoolerNumber = 1; SingleSpeedEvapFluidCoolerNumber <= NumSingleSpeedEvapFluidCoolers; ++SingleSpeedEvapFluidCoolerNumber ) {
			EvapFluidCoolerNum = SingleSpeedEvapFluidCoolerNumber;
			GetObjectItem( cCurrentModuleObject, SingleSpeedEvapFluidCoolerNumber, AlphArray, NumAlphas, NumArray, NumNums, IOStat, _, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );
			IsNotOK = false;
			IsBlank = false;
			VerifyName( instances.begin(), instances.begin() + EvapFluidCoolerNum - 1, AlphArray( 1 ), IsNotOK, IsBlank, cCurrentModuleObject + " Name" );
			if ( IsNotOK ) {
				ErrorsFound = true;
				if ( IsBlank ) AlphArray( 1 ) = "xxxxx";
			}
			auto & singleSpeedEvapFluidCooler( instances[ EvapFluidCoolerNum - 1 ] );
			singleSpeedEvapFluidCooler->Name = AlphArray( 1 );
			singleSpeedEvapFluidCooler->EvapFluidCoolerType = cCurrentModuleObject;
			singleSpeedEvapFluidCooler->EvapFluidCoolerType_Num = EvaporativeFluidCoolerType::SingleSpeed;
			singleSpeedEvapFluidCooler->PlantType_Num = DataPlant::TypeOf_EvapFluidCooler_SingleSpd;
			singleSpeedEvapFluidCooler->EvapFluidCoolerMassFlowRateMultiplier = 2.5;
			singleSpeedEvapFluidCooler->WaterInletNodeNum = GetOnlySingleNode( AlphArray( 2 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), NodeType_Water, NodeConnectionType_Inlet, 1, ObjectIsNotParent );
			singleSpeedEvapFluidCooler->WaterOutletNodeNum = GetOnlySingleNode( AlphArray( 3 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), NodeType_Water, NodeConnectionType_Outlet, 1, ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, AlphArray( 1 ), AlphArray( 2 ), AlphArray( 3 ), "Chilled Water Nodes" );
			singleSpeedEvapFluidCooler->HighSpeedAirFlowRate = NumArray( 1 );
			if ( singleSpeedEvapFluidCooler->HighSpeedAirFlowRate == AutoSize ) {
				singleSpeedEvapFluidCooler->HighSpeedAirFlowRateWasAutoSized = true;
			}
			singleSpeedEvapFluidCooler->HighSpeedFanPower = NumArray( 2 );
			if ( singleSpeedEvapFluidCooler->HighSpeedFanPower == AutoSize ) {
				singleSpeedEvapFluidCooler->HighSpeedFanPowerWasAutoSized = true;
			}
			singleSpeedEvapFluidCooler->DesignSprayWaterFlowRate = NumArray( 3 );
			singleSpeedEvapFluidCooler->HeatRejectCapNomCapSizingRatio = NumArray( 4 );
			singleSpeedEvapFluidCooler->HighSpeedStandardDesignCapacity = NumArray( 5 );
			singleSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA = NumArray( 6 );
			if ( singleSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA == AutoSize ) {
				singleSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUAWasAutoSized = true;
			}
			singleSpeedEvapFluidCooler->DesignWaterFlowRate = NumArray( 7 );
			if ( singleSpeedEvapFluidCooler->DesignWaterFlowRate == AutoSize ) {
				singleSpeedEvapFluidCooler->DesignWaterFlowRateWasAutoSized = true;
			}
			singleSpeedEvapFluidCooler->HighSpeedUserSpecifiedDesignCapacity = NumArray( 8 );
			singleSpeedEvapFluidCooler->DesignEnteringWaterTemp = NumArray( 9 );
			singleSpeedEvapFluidCooler->DesignEnteringAirTemp = NumArray( 10 );
			singleSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp = NumArray( 11 );

			if ( lAlphaFieldBlanks( 4 ) || AlphArray( 4 ).empty() ) {
				ShowSevereError( cCurrentModuleObject + ", \"" + singleSpeedEvapFluidCooler->Name + "\" Performance input method is not specified. " );
				ErrorsFound = true;
			}
			if ( SameString( AlphArray( 4 ), "STANDARDDESIGNCAPACITY" ) ) {
				singleSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::StandardDesignCapacity;
			}

			//outdoor air inlet node
			if ( lAlphaFieldBlanks( 5 ) ) {
				singleSpeedEvapFluidCooler->OutdoorAirInletNodeNum = 0;
			} else {
				singleSpeedEvapFluidCooler->OutdoorAirInletNodeNum = GetOnlySingleNode( AlphArray( 5 ), ErrorsFound, cCurrentModuleObject, singleSpeedEvapFluidCooler->Name, NodeType_Air, NodeConnectionType_OutsideAirReference, 1, ObjectIsNotParent );
				if ( ! CheckOutAirNodeNumber( singleSpeedEvapFluidCooler->OutdoorAirInletNodeNum ) ) {
					ShowSevereError( cCurrentModuleObject + ", \"" + singleSpeedEvapFluidCooler->Name + "\" Outdoor Air Inlet Node Name not valid Outdoor Air Node= " + AlphArray( 5 ) );
					ShowContinueError( "...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node." );
					ErrorsFound = true;
				}
			}

			//   fluid bypass for single speed evaporative fluid cooler
			if ( lAlphaFieldBlanks( 6 ) || AlphArray( 6 ).empty() ) {
				singleSpeedEvapFluidCooler->CapacityControl = 0; // FanCycling
			} else {
				{ auto const SELECT_CASE_var( MakeUPPERCase( AlphArray( 6 ) ) );
				if ( SELECT_CASE_var == "FANCYCLING" ) {
					singleSpeedEvapFluidCooler->CapacityControl = 0;
				} else if ( SELECT_CASE_var == "FLUIDBYPASS" ) {
					singleSpeedEvapFluidCooler->CapacityControl = 1;
				} else {
					singleSpeedEvapFluidCooler->CapacityControl = 0;
					ShowWarningError( cCurrentModuleObject + ", \"" + singleSpeedEvapFluidCooler->Name + "\" The Capacity Control is not specified correctly. The default Fan Cycling is used." );
				}}
			}

			singleSpeedEvapFluidCooler->SizFac = NumArray( 12 ); //  N11  \field Sizing Factor
			if ( singleSpeedEvapFluidCooler->SizFac <= 0.0 ) singleSpeedEvapFluidCooler->SizFac = 1.0;

			// begin water use and systems get input
			if ( SameString( AlphArray( 7 ), "LossFactor" ) ) {
				singleSpeedEvapFluidCooler->EvapLossMode = EvaporativeLossBy::UserFactor;
			} else if ( SameString( AlphArray( 7 ), "SaturatedExit" ) ) {
				singleSpeedEvapFluidCooler->EvapLossMode = EvaporativeLossBy::MoistTheory;
			} else if ( AlphArray( 7 ).empty() ) {
				singleSpeedEvapFluidCooler->EvapLossMode = EvaporativeLossBy::MoistTheory;
			} else {
				ShowSevereError( "Invalid, " + cAlphaFieldNames( 7 ) + " = " + AlphArray( 7 ) );
				ShowContinueError( "Entered in " + cCurrentModuleObject + " = " + AlphArray( 1 ) );
				ErrorsFound = true;
			}

			singleSpeedEvapFluidCooler->UserEvapLossFactor = NumArray( 13 ); //  N13 , \field Evaporation Loss Factor
			if ( ( NumNums < 13 ) && ( singleSpeedEvapFluidCooler->UserEvapLossFactor == 0.0 ) ) {
				// assume Evaporation loss factor not entered and should be calculated
				if ( ( OutRelHumValue >= 0.1 ) && ( OutRelHumValue <= 0.7 ) ) {
					//Use correlation by B.A. Qureshi and S.M. Zubair if within these limits
					singleSpeedEvapFluidCooler->UserEvapLossFactor = ( 113.0 - 8.417 * OutRelHumValue + 1.6147 * OutDryBulbTemp ) * 1.0e-5;
				} else { // Inlet conditions are out of the limit of correlation; An approximate default value of loss factor is used
					singleSpeedEvapFluidCooler->UserEvapLossFactor = 0.2;
				}
			}

			singleSpeedEvapFluidCooler->DriftLossFraction = NumArray( 14 ) / 100.0; //  N14, \field Drift Loss Percent

			if ( ( NumNums < 13 ) && ( singleSpeedEvapFluidCooler->DriftLossFraction == 0.0 ) ) {
				// assume Drift loss not entered and should be defaulted
				singleSpeedEvapFluidCooler->DriftLossFraction = 0.008 / 100.0;
			}
			singleSpeedEvapFluidCooler->ConcentrationRatio = NumArray( 15 ); //  N15, \field Blowdown Concentration Ratio

			if ( SameString( AlphArray( 8 ), "ScheduledRate" ) ) {
				singleSpeedEvapFluidCooler->BlowdownMode = BlowdownBy::Schedule;
			} else if ( SameString( AlphArray( 8 ), "ConcentrationRatio" ) ) {
				singleSpeedEvapFluidCooler->BlowdownMode = BlowdownBy::Concentration;
			} else if ( AlphArray( 8 ).empty() ) {
				singleSpeedEvapFluidCooler->BlowdownMode = BlowdownBy::Concentration;
				if ( ( NumNums < 15 ) && ( singleSpeedEvapFluidCooler->ConcentrationRatio == 0.0 ) ) {
					// assume Concetration ratio was omitted and should be defaulted
					singleSpeedEvapFluidCooler->ConcentrationRatio = 3.0;
				}
			} else {
				ShowSevereError( "Invalid, " + cAlphaFieldNames( 8 ) + " = " + AlphArray( 8 ) );
				ShowContinueError( "Entered in " + cCurrentModuleObject + " =" + AlphArray( 1 ) );
				ErrorsFound = true;
			}

			singleSpeedEvapFluidCooler->SchedIDBlowdown = GetScheduleIndex( AlphArray( 9 ) );
			if ( ( singleSpeedEvapFluidCooler->SchedIDBlowdown == 0 ) && ( singleSpeedEvapFluidCooler->BlowdownMode == BlowdownBy::Schedule ) ) {
				ShowSevereError( "Invalid, " + cAlphaFieldNames( 9 ) + " = " + AlphArray( 9 ) );
				ShowContinueError( "Entered in " + cCurrentModuleObject + " =" + AlphArray( 1 ) );
				ErrorsFound = true;
			}

			if ( AlphArray( 10 ).empty() ) {
				singleSpeedEvapFluidCooler->SuppliedByWaterSystem = false;
			} else { // water from storage tank
				SetupTankDemandComponent( AlphArray( 1 ), cCurrentModuleObject, AlphArray( 10 ), ErrorsFound, singleSpeedEvapFluidCooler->WaterTankID, singleSpeedEvapFluidCooler->WaterTankDemandARRID );
				singleSpeedEvapFluidCooler->SuppliedByWaterSystem = true;
			}

			//   Check various inputs to ensure that all the required variables are specified.

			if ( singleSpeedEvapFluidCooler->DesignSprayWaterFlowRate <= 0.0 ) {
				ShowSevereError( cCurrentModuleObject + " \"" + singleSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler input requires a design spray water flow rate greater than zero for all performance input methods." );
				ErrorsFound = true;
			}
			if ( singleSpeedEvapFluidCooler->HighSpeedAirFlowRate <= 0.0 && singleSpeedEvapFluidCooler->HighSpeedAirFlowRate != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 1 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( singleSpeedEvapFluidCooler->HighSpeedFanPower <= 0.0 && singleSpeedEvapFluidCooler->HighSpeedFanPower != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 2 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}

			if ( SameString( AlphArray( 4 ), "UFACTORTIMESAREAANDDESIGNWATERFLOWRATE" ) ) {
				singleSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::UFactor;
				if ( singleSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA <= 0.0 && singleSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 6 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->DesignWaterFlowRate <= 0.0 && singleSpeedEvapFluidCooler->DesignWaterFlowRate != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 7 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
			} else if ( SameString( AlphArray( 4 ), "STANDARDDESIGNCAPACITY" ) ) {
				singleSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::StandardDesignCapacity;
				if ( singleSpeedEvapFluidCooler->HighSpeedStandardDesignCapacity <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 5 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
			} else if ( SameString( AlphArray( 4 ), "USERSPECIFIEDDESIGNCAPACITY" ) ) {
				singleSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::UserSpecifiedDesignCapacity;
				if ( singleSpeedEvapFluidCooler->DesignWaterFlowRate <= 0.0 && singleSpeedEvapFluidCooler->DesignWaterFlowRate != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 7 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->HighSpeedUserSpecifiedDesignCapacity <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 8 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->DesignEnteringWaterTemp <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 9 ) + "\", entered value <= 0.0, but must be >0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->DesignEnteringAirTemp <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 10 ) + "\", entered value <= 0.0, but must be >0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 11 ) + "\", entered value <= 0.0, but must be >0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->DesignEnteringWaterTemp <= singleSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", " + cNumericFieldNames( 9 ) + " must be greater than " + cNumericFieldNames( 11 ) + '.' );
					ErrorsFound = true;
				}
				if ( singleSpeedEvapFluidCooler->DesignEnteringAirTemp <= singleSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", " + cNumericFieldNames( 10 ) + " must be greater than " + cNumericFieldNames( 11 ) + '.' );
					ErrorsFound = true;
				}
			} else { // Evaporative fluid cooler performance input method is not specified as a valid "choice"
				ShowSevereError( cCurrentModuleObject + " = \"" + singleSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler Performance Input Method must be \"UFactorTimesAreaAndDesignWaterFlowRate\" or \"StandardDesignCapacity\" or \"UserSpecifiedDesignCapacity\"." );
				ShowContinueError( "Evaporative fluid cooler Performance Input Method currently specified as: " + AlphArray( 4 ) );
				ErrorsFound = true;
			}
		} // End Single-Speed Evaporative Fluid Cooler Loop

		cCurrentModuleObject = cEvapFluidCooler_TwoSpeed;
		for ( TwoSpeedEvapFluidCoolerNumber = 1; TwoSpeedEvapFluidCoolerNumber <= NumTwoSpeedEvapFluidCoolers; ++TwoSpeedEvapFluidCoolerNumber ) {
			EvapFluidCoolerNum = NumSingleSpeedEvapFluidCoolers + TwoSpeedEvapFluidCoolerNumber;
			GetObjectItem( cCurrentModuleObject, TwoSpeedEvapFluidCoolerNumber, AlphArray, NumAlphas, NumArray, NumNums, IOStat, _, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

			IsNotOK = false;
			IsBlank = false;
			VerifyName( instances.begin(), instances.begin() + EvapFluidCoolerNum - 1, AlphArray( 1 ), IsNotOK, IsBlank, cCurrentModuleObject + " Name" );
			if ( IsNotOK ) {
				ErrorsFound = true;
				if ( IsBlank ) AlphArray( 1 ) = "xxxxx";
			}
			auto & twoSpeedEvapFluidCooler( instances[ EvapFluidCoolerNum - 1 ] );
			twoSpeedEvapFluidCooler->Name = AlphArray( 1 );
			twoSpeedEvapFluidCooler->EvapFluidCoolerType = cCurrentModuleObject;
			twoSpeedEvapFluidCooler->EvapFluidCoolerType_Num = EvaporativeFluidCoolerType::TwoSpeed;
			twoSpeedEvapFluidCooler->PlantType_Num = DataPlant::TypeOf_EvapFluidCooler_TwoSpd;
			twoSpeedEvapFluidCooler->EvapFluidCoolerMassFlowRateMultiplier = 2.5;
			twoSpeedEvapFluidCooler->WaterInletNodeNum = GetOnlySingleNode( AlphArray( 2 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), NodeType_Water, NodeConnectionType_Inlet, 1, ObjectIsNotParent );
			twoSpeedEvapFluidCooler->WaterOutletNodeNum = GetOnlySingleNode( AlphArray( 3 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), NodeType_Water, NodeConnectionType_Outlet, 1, ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, AlphArray( 1 ), AlphArray( 2 ), AlphArray( 3 ), "Chilled Water Nodes" );

			twoSpeedEvapFluidCooler->HighSpeedAirFlowRate = NumArray( 1 );
			if ( twoSpeedEvapFluidCooler->HighSpeedAirFlowRate == AutoSize ) {
				twoSpeedEvapFluidCooler->HighSpeedAirFlowRateWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->HighSpeedFanPower = NumArray( 2 );
			if ( twoSpeedEvapFluidCooler->HighSpeedFanPower == AutoSize ) {
				twoSpeedEvapFluidCooler->HighSpeedFanPowerWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->LowSpeedAirFlowRate = NumArray( 3 );
			if ( twoSpeedEvapFluidCooler->LowSpeedAirFlowRate == AutoSize ) {
				twoSpeedEvapFluidCooler->LowSpeedAirFlowRateWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->LowSpeedAirFlowRateSizingFactor = NumArray( 4 );
			twoSpeedEvapFluidCooler->LowSpeedFanPower = NumArray( 5 );
			if ( twoSpeedEvapFluidCooler->LowSpeedFanPower == AutoSize ) {
				twoSpeedEvapFluidCooler->LowSpeedFanPowerWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->LowSpeedFanPowerSizingFactor = NumArray( 6 );
			twoSpeedEvapFluidCooler->DesignSprayWaterFlowRate = NumArray( 7 );
			twoSpeedEvapFluidCooler->HeatRejectCapNomCapSizingRatio = NumArray( 8 );
			twoSpeedEvapFluidCooler->HighSpeedStandardDesignCapacity = NumArray( 9 );
			twoSpeedEvapFluidCooler->LowSpeedStandardDesignCapacity = NumArray( 10 );
			twoSpeedEvapFluidCooler->LowSpeedStandardDesignCapacitySizingFactor = NumArray( 11 );
			twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA = NumArray( 12 );
			if ( twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA == AutoSize ) {
				twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUAWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA = NumArray( 13 );
			if ( twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA == AutoSize ) {
				twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUAWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUASizingFactor = NumArray( 14 );
			twoSpeedEvapFluidCooler->DesignWaterFlowRate = NumArray( 15 );
			if ( twoSpeedEvapFluidCooler->DesignWaterFlowRate == AutoSize ) {
				twoSpeedEvapFluidCooler->DesignWaterFlowRateWasAutoSized = true;
			}
			twoSpeedEvapFluidCooler->HighSpeedUserSpecifiedDesignCapacity = NumArray( 16 );
			twoSpeedEvapFluidCooler->LowSpeedUserSpecifiedDesignCapacity = NumArray( 17 );
			twoSpeedEvapFluidCooler->LowSpeedUserSpecifiedDesignCapacitySizingFactor = NumArray( 18 );
			twoSpeedEvapFluidCooler->DesignEnteringWaterTemp = NumArray( 19 );
			twoSpeedEvapFluidCooler->DesignEnteringAirTemp = NumArray( 20 );
			twoSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp = NumArray( 21 );

			if ( lAlphaFieldBlanks( 4 ) ) {
				ShowSevereError( cCurrentModuleObject + ", \"" + twoSpeedEvapFluidCooler->Name + "\" Performance input method is not specified. " );
				ErrorsFound = true;
			}

			if ( SameString( AlphArray( 4 ), "STANDARDDESIGNCAPACITY" ) ) {
				twoSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::StandardDesignCapacity;
			}

			// outdoor air inlet node
			if ( lAlphaFieldBlanks( 5 ) ) {
				twoSpeedEvapFluidCooler->OutdoorAirInletNodeNum = 0;
			} else {
				twoSpeedEvapFluidCooler->OutdoorAirInletNodeNum = GetOnlySingleNode( AlphArray( 5 ), ErrorsFound, cCurrentModuleObject, twoSpeedEvapFluidCooler->Name, NodeType_Air, NodeConnectionType_OutsideAirReference, 1, ObjectIsNotParent );
				if ( ! CheckOutAirNodeNumber( twoSpeedEvapFluidCooler->OutdoorAirInletNodeNum ) ) {
					ShowSevereError( cCurrentModuleObject + ", \"" + twoSpeedEvapFluidCooler->Name + "\" Outdoor Air Inlet Node Name not valid Outdoor Air Node= " + AlphArray( 5 ) );
					ShowContinueError( "...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node." );
					ErrorsFound = true;
				}
			}

			twoSpeedEvapFluidCooler->SizFac = NumArray( 22 ); //  N16  \field Sizing Factor
			if ( twoSpeedEvapFluidCooler->SizFac <= 0.0 ) twoSpeedEvapFluidCooler->SizFac = 1.0;

			// begin water use and systems get input
			if ( SameString( AlphArray( 6 ), "LossFactor" ) ) {
				twoSpeedEvapFluidCooler->EvapLossMode = EvaporativeLossBy::UserFactor;
			} else if ( SameString( AlphArray( 6 ), "SaturatedExit" ) ) {
				twoSpeedEvapFluidCooler->EvapLossMode = EvaporativeLossBy::MoistTheory;
			} else if ( lAlphaFieldBlanks( 6 ) ) {
				twoSpeedEvapFluidCooler->EvapLossMode = EvaporativeLossBy::MoistTheory;
			} else {
				ShowSevereError( "Invalid " + cAlphaFieldNames( 6 ) + " = " + AlphArray( 6 ) );
				ShowContinueError( "Entered in " + cCurrentModuleObject + " = " + AlphArray( 1 ) );
				ErrorsFound = true;
			}

			twoSpeedEvapFluidCooler->UserEvapLossFactor = NumArray( 23 ); //  N23 , \field Evaporation Loss Factor
			if ( ( NumNums < 23 ) && ( twoSpeedEvapFluidCooler->UserEvapLossFactor == 0.0 ) ) {
				// assume Evaporation loss factor not entered and should be calculated
				if ( ( OutRelHumValue >= 0.1 ) && ( OutRelHumValue <= 0.7 ) ) {
					//Use correlation by B.A. Qureshi and S.M. Zubair if within these limits
					twoSpeedEvapFluidCooler->UserEvapLossFactor = ( 113.0 - 8.417 * OutRelHumValue + 1.6147 * OutDryBulbTemp ) * 1.0e-5;
				} else { // Inlet conditions are out of the limit of correlation; An approximate default value of loss factor is used
					twoSpeedEvapFluidCooler->UserEvapLossFactor = 0.2;
				}
			}
			twoSpeedEvapFluidCooler->DriftLossFraction = NumArray( 24 ) / 100.0; //  N24, \field Drift Loss Percent
			if ( ( NumNums < 24 ) && ( twoSpeedEvapFluidCooler->DriftLossFraction == 0.0 ) ) {
				// assume Drift loss not entered and should be defaulted
				twoSpeedEvapFluidCooler->DriftLossFraction = 0.008 / 100.0;
			}

			twoSpeedEvapFluidCooler->ConcentrationRatio = NumArray( 25 ); //  N25, \field Blowdown Concentration Ratio

			if ( SameString( AlphArray( 7 ), "ScheduledRate" ) ) {
				twoSpeedEvapFluidCooler->BlowdownMode = BlowdownBy::Schedule;
			} else if ( SameString( AlphArray( 7 ), "ConcentrationRatio" ) ) {
				twoSpeedEvapFluidCooler->BlowdownMode = BlowdownBy::Concentration;
			} else if ( lAlphaFieldBlanks( 7 ) ) {
				twoSpeedEvapFluidCooler->BlowdownMode = BlowdownBy::Concentration;
				if ( ( NumNums < 25 ) && ( twoSpeedEvapFluidCooler->ConcentrationRatio == 0.0 ) ) {
					// assume Concetration ratio was omitted and should be defaulted
					twoSpeedEvapFluidCooler->ConcentrationRatio = 3.0;
				}
			} else {
				ShowSevereError( "Invalid " + cAlphaFieldNames( 7 ) + " = " + AlphArray( 7 ) );
				ShowContinueError( "Entered in " + cCurrentModuleObject + " = " + AlphArray( 1 ) );
				ErrorsFound = true;
			}
			twoSpeedEvapFluidCooler->SchedIDBlowdown = GetScheduleIndex( AlphArray( 8 ) );
			if ( ( twoSpeedEvapFluidCooler->SchedIDBlowdown == 0 ) && ( twoSpeedEvapFluidCooler->BlowdownMode == BlowdownBy::Schedule ) ) {
				ShowSevereError( "Invalid " + cAlphaFieldNames( 8 ) + " = " + AlphArray( 8 ) );
				ShowContinueError( "Entered in " + cCurrentModuleObject + " = " + AlphArray( 1 ) );
				ErrorsFound = true;
			}

			if ( lAlphaFieldBlanks( 9 ) ) {
				twoSpeedEvapFluidCooler->SuppliedByWaterSystem = false;
			} else { // water from storage tank
				SetupTankDemandComponent( AlphArray( 1 ), cCurrentModuleObject, AlphArray( 9 ), ErrorsFound, twoSpeedEvapFluidCooler->WaterTankID, twoSpeedEvapFluidCooler->WaterTankDemandARRID );
				twoSpeedEvapFluidCooler->SuppliedByWaterSystem = true;
			}

			//   Check various inputs to ensure that all the required variables are specified.

			if ( twoSpeedEvapFluidCooler->DesignSprayWaterFlowRate <= 0.0 ) {
				ShowSevereError( cCurrentModuleObject + " \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler input requires a design spray water flow rate greater than zero for all performance input methods." );
				ErrorsFound = true;
			}
			if ( twoSpeedEvapFluidCooler->HighSpeedAirFlowRate <= 0.0 && twoSpeedEvapFluidCooler->HighSpeedAirFlowRate != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler input requires design air flow rate at high fan speed to be greater than zero for all performance input methods." );
				ErrorsFound = true;
			}
			if ( twoSpeedEvapFluidCooler->LowSpeedAirFlowRate <= 0.0 && twoSpeedEvapFluidCooler->LowSpeedAirFlowRate != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler input requires design air flow rate at low fan speed to be greater than zero for all performance input methods." );
				ErrorsFound = true;
			}
			//   High speed air flow rate must be greater than low speed air flow rate.
			//   Can't tell yet if autosized, check later in InitEvapFluidCooler.
			if ( twoSpeedEvapFluidCooler->HighSpeedAirFlowRate <= twoSpeedEvapFluidCooler->LowSpeedAirFlowRate && twoSpeedEvapFluidCooler->HighSpeedAirFlowRate != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler air flow rate at low fan speed must be less than the air flow rate at high fan speed." );
				ErrorsFound = true;
			}
			if ( twoSpeedEvapFluidCooler->HighSpeedFanPower <= 0.0 && twoSpeedEvapFluidCooler->HighSpeedFanPower != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 2 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( twoSpeedEvapFluidCooler->LowSpeedFanPower <= 0.0 && twoSpeedEvapFluidCooler->LowSpeedFanPower != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 5 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( twoSpeedEvapFluidCooler->HighSpeedFanPower <= twoSpeedEvapFluidCooler->LowSpeedFanPower && twoSpeedEvapFluidCooler->HighSpeedFanPower != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler low speed fan power must be less than the high speed fan power ." );
				ErrorsFound = true;
			}

			if ( SameString( AlphArray( 4 ), "UFACTORTIMESAREAANDDESIGNWATERFLOWRATE" ) ) {
				twoSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::UFactor;
				if ( twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA <= 0.0 && twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 12 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA <= 0.0 && twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 13 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA <= twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA && twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler U-factor Times Area Value at Low Fan Speed must be less than the U-factor Times Area Value at High Fan Speed." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->DesignWaterFlowRate <= 0.0 && twoSpeedEvapFluidCooler->DesignWaterFlowRate != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 15 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
			} else if ( SameString( AlphArray( 4 ), "STANDARDDESIGNCAPACITY" ) ) {
				twoSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::StandardDesignCapacity;
				if ( twoSpeedEvapFluidCooler->HighSpeedStandardDesignCapacity <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 9 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->LowSpeedStandardDesignCapacity <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 10 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->LowSpeedStandardDesignCapacity >= twoSpeedEvapFluidCooler->HighSpeedStandardDesignCapacity ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". Low-Speed Standard Design Capacity must be less than the High-Speed Standard Design Capacity." );
					ErrorsFound = true;
				}
			} else if ( SameString( AlphArray( 4 ), "USERSPECIFIEDDESIGNCAPACITY" ) ) {
				twoSpeedEvapFluidCooler->PerformanceInputMethod_Num = PIM::UserSpecifiedDesignCapacity;
				if ( twoSpeedEvapFluidCooler->DesignWaterFlowRate <= 0.0 && twoSpeedEvapFluidCooler->DesignWaterFlowRate != AutoSize ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 15 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->HighSpeedUserSpecifiedDesignCapacity <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 16 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->LowSpeedUserSpecifiedDesignCapacity <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 17 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA != 0.0 ) {
					if ( twoSpeedEvapFluidCooler->HighSpeedEvapFluidCoolerUA > 0.0 ) {
						ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". UserSpecifiedDesignCapacity performance input method and evaporative fluid cooler UA at high fan speed have been specified." );
					} else {
						ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". UserSpecifiedDesignCapacity performance input method has been specified and evaporative fluid cooler UA at high fan speed is being autosized." );
					}
					ShowContinueError( "Evaporative fluid cooler UA at high fan speed must be left blank when UserSpecifiedDesignCapacity performance input method is used." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA != 0.0 ) {
					if ( twoSpeedEvapFluidCooler->LowSpeedEvapFluidCoolerUA > 0.0 ) {
						ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". UserSpecifiedDesignCapacity performance input method and evaporative fluid cooler UA at low fan speed have been specified." );
					} else {
						ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". UserSpecifiedDesignCapacity performance input method has been specified and evaporative fluid cooler UA at low fan speed is being autosized." );
					}
					ShowContinueError( "Evaporative fluid cooler UA at low fan speed must be left blank when UserSpecifiedDesignCapacity performance input method is used." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->LowSpeedUserSpecifiedDesignCapacity >= twoSpeedEvapFluidCooler->HighSpeedUserSpecifiedDesignCapacity ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". Low-Speed User Specified Design Capacity must be less than the High-Speed User Specified Design Dapacity." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->DesignEnteringWaterTemp <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 19 ) + "\", entered value <= 0.0, but must be >0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->DesignEnteringAirTemp <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 20 ) + "\", entered value <= 0.0, buy must be >0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp <= 0.0 ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 21 ) + "\", entered value <= 0.0, but must be >0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->DesignEnteringWaterTemp <= twoSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", " + cNumericFieldNames( 19 ) + " must be greater than " + cNumericFieldNames( 15 ) + '.' );
					ErrorsFound = true;
				}
				if ( twoSpeedEvapFluidCooler->DesignEnteringAirTemp <= twoSpeedEvapFluidCooler->DesignEnteringAirWetBulbTemp ) {
					ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", " + cNumericFieldNames( 20 ) + " must be greater than " + cNumericFieldNames( 15 ) + '.' );
					ErrorsFound = true;
				}
			} else { // Evaporative fluid cooler performance input method is not specified as a valid "choice"
				ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedEvapFluidCooler->Name + "\". Evaporative fluid cooler Performance Input Method must be \"UFactorTimesAreaAndDesignWaterFlowRate\" or \"StandardDesignCapacity\" or \"UserSpecifiedDesignCapacity\"." );
				ShowContinueError( "Evaporative fluid cooler Performanace Input Method currently specified as: " + AlphArray( 4 ) );
				ErrorsFound = true;
			}
		} // End Two-Speed Evaporative Fluid Cooler Loop

		if ( ErrorsFound ) {
			ShowFatalError( "Errors found in getting evaporative fluid cooler input." );
		}

		// Set up output variables
		// CurrentModuleObject='EvaporativeFluidCooler:SingleSpeed'
		for ( EvapFluidCoolerNum = 1; EvapFluidCoolerNum <= NumSingleSpeedEvapFluidCoolers; ++EvapFluidCoolerNum ) {
			auto & singleSpeedEvapFluidCooler( instances[ EvapFluidCoolerNum - 1 ] );
			SetupOutputVariable( "Cooling Tower Inlet Temperature [C]", singleSpeedEvapFluidCooler->InletWaterTemp, "System", "Average", singleSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Outlet Temperature [C]", singleSpeedEvapFluidCooler->OutletWaterTemp, "System", "Average", singleSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Mass Flow Rate [kg/s]", singleSpeedEvapFluidCooler->WaterMassFlowRate, "System", "Average", singleSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Heat Transfer Rate [W]", singleSpeedEvapFluidCooler->Qactual, "System", "Average", singleSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Power [W]", singleSpeedEvapFluidCooler->FanPower, "System", "Average", singleSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Energy [J]", singleSpeedEvapFluidCooler->FanEnergy, "System", "Sum", singleSpeedEvapFluidCooler->Name, _, "Electric", "HeatRejection", _, "Plant" );
			// Added for fluid bypass
			SetupOutputVariable( "Cooling Tower Bypass Fraction []", singleSpeedEvapFluidCooler->BypassFraction, "System", "Average", singleSpeedEvapFluidCooler->Name );
		}

		// CurrentModuleObject='EvaporativeFluidCooler:TwoSpeed'
		for ( EvapFluidCoolerNum = NumSingleSpeedEvapFluidCoolers + 1; EvapFluidCoolerNum <= NumSingleSpeedEvapFluidCoolers + NumTwoSpeedEvapFluidCoolers; ++EvapFluidCoolerNum ) {
			auto & twoSpeedEvapFluidCooler( instances[ EvapFluidCoolerNum - 1 ] );
			SetupOutputVariable( "Cooling Tower Inlet Temperature [C]", twoSpeedEvapFluidCooler->InletWaterTemp, "System", "Average", twoSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Outlet Temperature [C]", twoSpeedEvapFluidCooler->OutletWaterTemp, "System", "Average", twoSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Mass Flow Rate [kg/s]", twoSpeedEvapFluidCooler->WaterMassFlowRate, "System", "Average", twoSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Heat Transfer Rate [W]", twoSpeedEvapFluidCooler->Qactual, "System", "Average", twoSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Power [W]", twoSpeedEvapFluidCooler->FanPower, "System", "Average", twoSpeedEvapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Energy [J]", twoSpeedEvapFluidCooler->FanEnergy, "System", "Sum", twoSpeedEvapFluidCooler->Name, _, "Electric", "HeatRejection", _, "Plant" );

		}

		// setup common water reporting for all types of evaporative fluid coolers.
		// CurrentModuleObject='EvaporativeFluidCooler:*'
		for ( EvapFluidCoolerNum = 1; EvapFluidCoolerNum <= NumSingleSpeedEvapFluidCoolers + NumTwoSpeedEvapFluidCoolers; ++EvapFluidCoolerNum ) {
			auto & evapFluidCooler( instances[ EvapFluidCoolerNum - 1 ] );
			if ( evapFluidCooler->SuppliedByWaterSystem ) {
				SetupOutputVariable( "Cooling Tower Make Up Water Volume Flow Rate [m3/s]", evapFluidCooler->MakeUpVdot, "System", "Average", evapFluidCooler->Name );
				SetupOutputVariable( "Cooling Tower Make Up Water Volume [m3]", evapFluidCooler->MakeUpVol, "System", "Sum", evapFluidCooler->Name );
				SetupOutputVariable( "Cooling Tower Storage Tank Water Volume Flow Rate [m3/s]", evapFluidCooler->TankSupplyVdot, "System", "Average", evapFluidCooler->Name );
				SetupOutputVariable( "Cooling Tower Storage Tank Water Volume [m3]", evapFluidCooler->TankSupplyVol, "System", "Sum", evapFluidCooler->Name, _, "Water", "HeatRejection", _, "Plant" );
				SetupOutputVariable( "Cooling Tower Starved Storage Tank Water Volume Flow Rate [m3/s]", evapFluidCooler->StarvedMakeUpVdot, "System", "Average", evapFluidCooler->Name );
				SetupOutputVariable( "Cooling Tower Starved Storage Tank Water Volume [m3]", evapFluidCooler->StarvedMakeUpVol, "System", "Sum", evapFluidCooler->Name, _, "Water", "HeatRejection", _, "Plant" );
				SetupOutputVariable( "Cooling Tower Make Up Mains Water Volume [m3]", evapFluidCooler->StarvedMakeUpVol, "System", "Sum", evapFluidCooler->Name, _, "MainsWater", "HeatRejection", _, "Plant" );
			} else { // Evaporative fluid cooler water from mains and gets metered
				SetupOutputVariable( "Cooling Tower Make Up Water Volume Flow Rate [m3/s]", evapFluidCooler->MakeUpVdot, "System", "Average", evapFluidCooler->Name );
				SetupOutputVariable( "Cooling Tower Make Up Water Volume [m3]", evapFluidCooler->MakeUpVol, "System", "Sum", evapFluidCooler->Name, _, "Water", "HeatRejection", _, "Plant" );
				SetupOutputVariable( "Cooling Tower Make Up Mains Water Volume [m3]", evapFluidCooler->MakeUpVol, "System", "Sum", evapFluidCooler->Name, _, "MainsWater", "HeatRejection", _, "Plant" );
			}

			SetupOutputVariable( "Cooling Tower Water Evaporation Volume Flow Rate [m3/s]", evapFluidCooler->EvaporationVdot, "System", "Average", evapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Water Evaporation Volume [m3]", evapFluidCooler->EvaporationVol, "System", "Sum", evapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Water Drift Volume Flow Rate [m3/s]", evapFluidCooler->DriftVdot, "System", "Average", evapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Water Drift Volume [m3]", evapFluidCooler->DriftVol, "System", "Sum", evapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Water Blowdown Volume Flow Rate [m3/s]", evapFluidCooler->BlowdownVdot, "System", "Average", evapFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Water Blowdown Volume [m3]", evapFluidCooler->BlowdownVol, "System", "Sum", evapFluidCooler->Name );
		} // loop all evaporative fluid coolers

	}

	void
	EvaporativeFluidCooler::init()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for initializations of the evaporative fluid cooler components and for
		// final checking of evaporative fluid cooler inputs (post autosizing)

		// METHODOLOGY EMPLOYED:
		// Uses the status flags to trigger initializations.

		// REFERENCES:
		// Based on InitTower subroutine by Don Shirey Sept/Oct 2002, F Buhl Oct 2002

		// Using/Aliasing
		using DataGlobals::BeginEnvrnFlag;
		using Psychrometrics::PsyTwbFnTdbWPb;
		using InputProcessor::SameString;
		//  USE FluidProperties, ONLY : GetDensityGlycol
		using DataPlant::TypeOf_EvapFluidCooler_SingleSpd;
		using DataPlant::TypeOf_EvapFluidCooler_TwoSpd;
		using DataPlant::ScanPlantLoopsForObject;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using PlantUtilities::InitComponentNodes;
		using PlantUtilities::SetComponentFlowRate;
		using PlantUtilities::RegulateCondenserCompFlowReqOp;

		static std::string const RoutineName( "InitEvapFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		bool ErrorsFound( false ); // Flag if input data errors are found
		int TypeOf_Num( 0 );
		Real64 rho; // local density of fluid

		// from InitSimVars()
		this->InletWaterTemp = 0.0; // CW temperature at evaporative fluid cooler inlet
		this->OutletWaterTemp = 0.0; // CW temperature at evaporative fluid cooler outlet
		this->WaterMassFlowRate = 0.0; // WaterMassFlowRate through evaporative fluid cooler
		this->Qactual = 0.0; // Evaporative fluid cooler heat transfer
		this->FanPower = 0.0; // Evaporative fluid cooler fan power used
		this->AirFlowRatio = 0.0; // Ratio of air flow rate through VS Evaporative fluid cooler to design air flow rate
		this->WaterUsage = 0.0; // Evaporative fluid cooler water usage (m3/s)

		if ( oneTimeFlag ) {

			if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
				TypeOf_Num = TypeOf_EvapFluidCooler_SingleSpd;
			} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
				TypeOf_Num = TypeOf_EvapFluidCooler_TwoSpd;
			} else {
				assert( false );
			}
			ErrorsFound = false;
			// Locate the tower on the plant loops for later usage
			ScanPlantLoopsForObject( this->Name, TypeOf_Num, this->location.loopNum, this->location.loopSideNum, this->location.branchNum, this->location.compNum, _, _, _, _, _, ErrorsFound );

			if ( ErrorsFound ) {
				ShowFatalError( "InitEvapFluidCooler: Program terminated due to previous condition(s)." );
			}

			if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
				if ( this->DesignWaterFlowRate > 0.0 ) {
					if ( this->HighSpeedAirFlowRate <= this->LowSpeedAirFlowRate ) {
						ShowSevereError( "EvaporativeFluidCooler:TwoSpeed \"" + this->Name + "\". Low speed air flow rate must be less than the high speed air flow rate." );
						ErrorsFound = true;
					}
					if ( this->HighSpeedEvapFluidCoolerUA <= this->LowSpeedEvapFluidCoolerUA ) {
						ShowSevereError( "EvaporativeFluidCooler:TwoSpeed \"" + this->Name + "\". Evaporative fluid cooler UA at low fan speed must be less than the evaporative fluid cooler UA at high fan speed." );
						ErrorsFound = true;
					}
				}
			}

			this->FluidIndex = PlantLoop( this->location.loopNum ).FluidIndex;
			if ( this->PerformanceInputMethod_Num == PIM::StandardDesignCapacity ) {
				std::string fluidName = FluidProperties::GetGlycolNameByIndex( this->FluidIndex );
				if ( fluidName != "WATER" ) {
					ShowSevereError( this->EvapFluidCoolerType + " = \"" + this->Name + "\". StandardDesignCapacity performance input method is only valid for fluid type = \"Water\"." );
					ShowContinueError( "Currently, Fluid Type = " + fluidName + " in CondenserLoop = " + PlantLoop( this->location.loopNum ).Name );
					ErrorsFound = true;
				}
			}

			if ( ErrorsFound ) {
				ShowFatalError( "InitEvapFluidCooler: Program terminated due to previous condition(s)." );
			}

			oneTimeFlag = false;

		}

		// Begin environment initializations
		if ( envrnFlag && BeginEnvrnFlag && ( PlantFirstSizesOkayToFinalize ) ) {

			rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
			this->DesWaterMassFlowRate = this->DesignWaterFlowRate * rho;
			InitComponentNodes( 0.0, this->DesWaterMassFlowRate, this->WaterInletNodeNum, this->WaterOutletNodeNum, this->location.loopNum, this->location.loopSideNum, this->location.branchNum, this->location.compNum );
			envrnFlag = false;
		}

		if ( ! BeginEnvrnFlag ) {
			envrnFlag = true;
		}

		// Each time initializations
		this->WaterTemp = Node( this->WaterInletNodeNum ).Temp;

		if ( this->OutdoorAirInletNodeNum != 0 ) {
			this->AirTemp = Node( this->OutdoorAirInletNodeNum ).Temp;
			this->AirHumRat = Node( this->OutdoorAirInletNodeNum ).HumRat;
			this->AirPress = Node( this->OutdoorAirInletNodeNum ).Press;
			this->AirWetBulb = Node( this->OutdoorAirInletNodeNum ).OutAirWetBulb;
		} else {
			this->AirTemp = DataEnvironment::OutDryBulbTemp;
			this->AirHumRat = DataEnvironment::OutHumRat;
			this->AirPress = DataEnvironment::OutBaroPress;
			this->AirWetBulb = DataEnvironment::OutWetBulbTemp;
		}

		WaterMassFlowRate = RegulateCondenserCompFlowReqOp( this->location.loopNum, this->location.loopSideNum, this->location.branchNum, this->location.compNum, this->DesWaterMassFlowRate * this->EvapFluidCoolerMassFlowRateMultiplier );

		SetComponentFlowRate( WaterMassFlowRate, this->WaterInletNodeNum, this->WaterOutletNodeNum, this->location.loopNum, this->location.loopSideNum, this->location.branchNum, this->location.compNum );
	}

	void
	EvaporativeFluidCooler::size()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       Chandan Sharma, April 2010
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for sizing evaporative fluid cooler Components for which capacities and flow rates
		// have not been specified in the input. This subroutine also calculates evaporative fluid cooler UA if the user
		// has specified evaporative fluid cooler performance via the "Standard Design Capacity" method.

		// METHODOLOGY EMPLOYED:
		// Obtains condenser flow rate from the plant sizing array. If evaporative fluid cooler performance is specified
		// via the "Standard Design Capacity" method, the water flow rate is directly proportional to capacity.

		// REFERENCES:
		// Based on SizeTower by Don Shirey, Sept/Oct 2002; Richard Raustad, Feb 2005

		// Using/Aliasing
		using namespace DataSizing;
		using General::SolveRegulaFalsi;
		using General::RoundSigDigits;
		using PlantUtilities::RegisterPlantCompDesignFlow;
		using ReportSizingManager::ReportSizingOutput;
		using namespace OutputReportPredefined;
		using InputProcessor::SameString;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using DataPlant::PlantFirstSizesOkayToReport;
		using DataPlant::PlantFinalSizesOkayToReport;

		int const MaxIte( 500 ); // Maximum number of iterations
		Real64 const Acc( 0.0001 ); // Accuracy of result
		static std::string const CalledFrom( "SizeEvapFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int PltSizCondNum; // Plant Sizing index for condenser loop
		int SolFla; // Flag of solver
		Real64 DesEvapFluidCoolerLoad; // Design evaporative fluid cooler load [W]
		Real64 UA0; // Lower bound for UA [W/C]
		Real64 UA1; // Upper bound for UA [W/C]
		Real64 UA; // Calculated UA value [W/C]
		Real64 OutWaterTempAtUA0; // Water outlet temperature at UA0
		Real64 OutWaterTempAtUA1; // Water outlet temperature at UA1
		Real64 DesignEnteringAirWetBulb; // Intermediate variable to check that design exit
		// temperature specified in the plant:sizing object
		// is higher than the design entering air wet-bulb temp
		// when autosize feature is used
		Array1D< Real64 > Par( 4 ); // Parameter array need for RegulaFalsi routine
		std::string equipName;
		Real64 Cp; // local specific heat for fluid
		Real64 rho; // local density for fluid
		Real64 tmpDesignWaterFlowRate; // local temporary for water volume flow rate
		Real64 tmpHighSpeedFanPower; // local temporary for high speed fan power
		Real64 tmpHighSpeedAirFlowRate; // local temporary for high speed air flow rate
		Real64 tmpHighSpeedEvapFluidCoolerUA; // local temporary for high speed cooler UA

		DesEvapFluidCoolerLoad = 0.0;
		tmpDesignWaterFlowRate = this->DesignWaterFlowRate;
		tmpHighSpeedFanPower = this->HighSpeedFanPower;
		tmpHighSpeedAirFlowRate = this->HighSpeedAirFlowRate;
		tmpHighSpeedEvapFluidCoolerUA = this->HighSpeedEvapFluidCoolerUA;

		PltSizCondNum = PlantLoop( this->location.loopNum ).PlantSizNum;

		if ( this->DesignWaterFlowRateWasAutoSized && this->PerformanceInputMethod_Num != PIM::StandardDesignCapacity ) {
			if ( PltSizCondNum > 0 ) {
				if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= SmallWaterVolFlow ) {
					tmpDesignWaterFlowRate = PlantSizData( PltSizCondNum ).DesVolFlowRate * this->SizFac;
					if ( PlantFirstSizesOkayToFinalize ) this->DesignWaterFlowRate = tmpDesignWaterFlowRate;

				} else {
					tmpDesignWaterFlowRate = 0.0;
					if ( PlantFirstSizesOkayToFinalize ) this->DesignWaterFlowRate = tmpDesignWaterFlowRate;
				}
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
							"Design Water Flow Rate [m3/s]", this->DesignWaterFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
							"Initial Design Water Flow Rate [m3/s]", this->DesignWaterFlowRate );
					}
				}
			} else {
				if ( PlantFirstSizesOkayToFinalize ) {
					ShowSevereError( "Autosizing error for evaporative fluid cooler object = " + this->Name );
					ShowFatalError( "Autosizing of evaporative fluid cooler condenser flow rate requires a loop Sizing:Plant object." );
				}
			}
			// Check when the user specified Condenser/Evaporative Fluid Cooler water design setpoint
			// temperature is less than design inlet air wet bulb temperature
			if ( this->PerformanceInputMethod_Num == PIM::UFactor ) {
				DesignEnteringAirWetBulb = 25.6;
			} else {
				DesignEnteringAirWetBulb = this->DesignEnteringAirWetBulbTemp;
			}
			if ( PlantSizData( PltSizCondNum ).ExitTemp <= DesignEnteringAirWetBulb ) {
				ShowSevereError( "Error when autosizing the UA value for Evaporative Fluid Cooler = " + this->Name + '.' );
				ShowContinueError( "Design Loop Exit Temperature (" + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) + " C) must be greater than design entering air wet-bulb temperature (" + RoundSigDigits( DesignEnteringAirWetBulb, 2 ) + " C) when autosizing the Evaporative Fluid Cooler UA." );
				ShowContinueError( "It is recommended that the Design Loop Exit Temperature = Design Entering Air Wet-bulb Temp plus the Evaporative Fluid Cooler design approach temperature (e.g., 4 C)." );
				ShowContinueError( "If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be > Design Entering Air Wet-bulb Temp if autosizing the Evaporative Fluid Cooler." );
				ShowFatalError( "Review and revise design input values as appropriate." );
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::UFactor
				&& ! this->HighSpeedEvapFluidCoolerUAWasAutoSized ) {
			if ( PltSizCondNum > 0 ) {
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesEvapFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
				this->HighSpeedStandardDesignCapacity = DesEvapFluidCoolerLoad / this->HeatRejectCapNomCapSizingRatio;
			} else {
				this->HighSpeedStandardDesignCapacity = 0.0;
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::StandardDesignCapacity ) {
			// Design water flow rate is assumed to be 3 gpm per ton (SI equivalent 5.382E-8 m3/s per watt)
			tmpDesignWaterFlowRate = 5.382e-8 * this->HighSpeedStandardDesignCapacity;
			if ( PlantFirstSizesOkayToFinalize ) {
				this->DesignWaterFlowRate = tmpDesignWaterFlowRate;
				if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Design Water Flow Rate based on evaporative fluid cooler Standard Design Capacity [m3/s]", this->DesignWaterFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Initial Design Water Flow Rate based on evaporative fluid cooler Standard Design Capacity [m3/s]", this->DesignWaterFlowRate );
					}
				} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Design Water Flow Rate based on evaporative fluid cooler high-speed Standard Design Capacity [m3/s]", this->DesignWaterFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Initial Design Water Flow Rate based on evaporative fluid cooler high-speed Standard Design Capacity [m3/s]", this->DesignWaterFlowRate );
					}
				}
			}
		}

		RegisterPlantCompDesignFlow( this->WaterInletNodeNum, tmpDesignWaterFlowRate );

		if ( this->HighSpeedFanPowerWasAutoSized ) {
			// We assume the nominal fan power is 0.0105 times the design load
			if ( this->PerformanceInputMethod_Num == PIM::StandardDesignCapacity ) {
				tmpHighSpeedFanPower = 0.0105 * this->HighSpeedStandardDesignCapacity;
				if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
			} else if ( this->PerformanceInputMethod_Num == PIM::UserSpecifiedDesignCapacity ) {
				tmpHighSpeedFanPower = 0.0105 * this->HighSpeedUserSpecifiedDesignCapacity;
				if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
			} else {
				if ( DesEvapFluidCoolerLoad > 0 ) {
					tmpHighSpeedFanPower = 0.0105 * DesEvapFluidCoolerLoad;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
				} else if ( PltSizCondNum > 0 ) {
					if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= SmallWaterVolFlow ) {
						rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
						Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
						DesEvapFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
						tmpHighSpeedFanPower = 0.0105 * DesEvapFluidCoolerLoad;
						if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
					} else {
						tmpHighSpeedFanPower = 0.0;
						if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
					}
				} else {
					if ( PlantFirstSizesOkayToFinalize ) {
						ShowSevereError( "Autosizing of evaporative fluid cooler fan power requires a loop Sizing:Plant object." );
						ShowFatalError( " Occurs in evaporative fluid cooler object= " + this->Name );
					}
				}
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Fan Power at Design Air Flow Rate [W]", this->HighSpeedFanPower );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Initial Fan Power at Design Air Flow Rate [W]", this->HighSpeedFanPower );
					}
				} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Fan Power at High Fan Speed [W]", this->HighSpeedFanPower );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Initial Fan Power at High Fan Speed [W]", this->HighSpeedFanPower );
					}
				}
			}
		}

		if ( this->HighSpeedAirFlowRateWasAutoSized ) {
			// Plant Sizing Object is not required to AUTOSIZE this field since its simply a multiple of another field.

			tmpHighSpeedAirFlowRate = tmpHighSpeedFanPower * 0.5 * ( 101325.0 / StdBaroPress ) / 190.0;
			if ( PlantFirstSizesOkayToFinalize ) {
				this->HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;

				if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Design Air Flow Rate [m3/s]", this->HighSpeedAirFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Initial Design Air Flow Rate [m3/s]", this->HighSpeedAirFlowRate );
					}
				} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Air Flow Rate at High Fan Speed [m3/s]", this->HighSpeedAirFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Initial Air Flow Rate at High Fan Speed [m3/s]", this->HighSpeedAirFlowRate );
					}
				}
			}
		}

		if ( this->HighSpeedEvapFluidCoolerUAWasAutoSized && this->PerformanceInputMethod_Num == PIM::UFactor ) {
			if ( PltSizCondNum > 0 ) {
				if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= SmallWaterVolFlow ) {
					// This conditional statement is to trap when the user specified Condenser/Evaporative Fluid Cooler water design setpoint
					// temperature is less than design inlet air wet bulb temperature of 25.6 C
					if ( PlantSizData( PltSizCondNum ).ExitTemp <= 25.6 ) {
						ShowSevereError( "Error when autosizing the UA value for Evaporative Fluid Cooler = " + this->Name + '.' );
						ShowContinueError( "Design Loop Exit Temperature (" + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) + " C) must be greater than 25.6 C when autosizing the Evaporative Fluid Cooler UA." );
						ShowContinueError( "The Design Loop Exit Temperature specified in Sizing:Plant object = " + PlantSizData( PltSizCondNum ).PlantLoopName );
						ShowContinueError( "It is recommended that the Design Loop Exit Temperature = 25.6 C plus the Evaporative Fluid Cooler design approach temperature (e.g., 4 C)." );
						ShowContinueError( "If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be > 25.6 C if autosizing the Evaporative Fluid Cooler." );
						ShowFatalError( "Review and revise design input values as appropriate." );
					}
					rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
					Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
					DesEvapFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
					Par( 1 ) = DesEvapFluidCoolerLoad;
					Par( 2 ) = rho * tmpDesignWaterFlowRate; // Design water mass flow rate
					Par( 3 ) = tmpHighSpeedAirFlowRate; // Design air volume flow rate
					Par( 4 ) = Cp;
					UA0 = 0.0001 * DesEvapFluidCoolerLoad; // Assume deltaT = 10000K (limit)
					UA1 = DesEvapFluidCoolerLoad; // Assume deltaT = 1K
					this->WaterTemp = PlantSizData( PltSizCondNum ).ExitTemp + PlantSizData( PltSizCondNum ).DeltaT;
					this->AirTemp = 35.0;
					this->AirWetBulb = 25.6;
					this->AirPress = StdBaroPress;
					this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress );
					SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleEvapFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
					if ( SolFla == -1 ) {
						ShowWarningError( "Iteration limit exceeded in calculating evaporative fluid cooler UA." );
						ShowContinueError( "Autosizing of fluid cooler UA failed for evaporative fluid cooler = " + this->Name );
						ShowContinueError( "The final UA value = " + RoundSigDigits( UA, 2 ) + "W/C, and the simulation continues..." );
					} else if ( SolFla == -2 ) {
						this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA0, OutWaterTempAtUA0 );
						this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA1, OutWaterTempAtUA1 );
						ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
						ShowContinueError( "reasonable UA value. Review and revise design input values as appropriate. Specifying hard" );
						ShowContinueError( "sizes for some \"autosizable\" fields while autosizing other \"autosizable\" fields may be contributing to this problem." );
						ShowContinueError( "This model iterates on UA to find the heat transfer required to provide the design outlet " );
						ShowContinueError( "water temperature. Initially, the outlet water temperatures at high and low UA values are " );
						ShowContinueError( "calculated. The Design Exit Water Temperature should be between the outlet water " );
						ShowContinueError( "temperatures calculated at high and low UA values. If the Design Exit Water Temperature is " );
						ShowContinueError( "out of this range, the solution will not converge and UA will not be calculated. " );
						ShowContinueError( "The possible solutions could be to manually input adjusted water and/or air flow rates " );
						ShowContinueError( "based on the autosized values shown below or to adjust design evaporative fluid cooler air inlet wet-bulb temperature." );
						ShowContinueError( "Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp)." );
						ShowContinueError( "Inputs to the evaporative fluid cooler object:" );
						ShowContinueError( "Design Evaporative Fluid Cooler Load [W]                      = " + RoundSigDigits( Par( 1 ), 2 ) );
						ShowContinueError( "Design Evaporative Fluid Cooler Water Volume Flow Rate [m3/s] = " + RoundSigDigits( this->DesignWaterFlowRate, 6 ) );
						ShowContinueError( "Design Evaporative Fluid Cooler Air Volume Flow Rate [m3/s]   = " + RoundSigDigits( Par( 3 ), 2 ) );
						ShowContinueError( "Design Evaporative Fluid Cooler Air Inlet Wet-bulb Temp [C]   = " + RoundSigDigits( this->AirWetBulb, 2 ) );
						ShowContinueError( "Design Evaporative Fluid Cooler Water Inlet Temp [C]          = " + RoundSigDigits( this->WaterTemp, 2 ) );
						ShowContinueError( "Inputs to the plant sizing object:" );
						ShowContinueError( "Design Exit Water Temp [C]                                    = " + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) );
						ShowContinueError( "Loop Design Temperature Difference [C]                        = " + RoundSigDigits( PlantSizData( PltSizCondNum ).DeltaT, 2 ) );
						ShowContinueError( "Design Evaporative Fluid Cooler Water Inlet Temp [C]          = " + RoundSigDigits( this->WaterTemp, 2 ) );
						ShowContinueError( "Calculated water outlet temperature at low UA [C](UA = " + RoundSigDigits( UA0, 2 ) + " W/C)  = " + RoundSigDigits( OutWaterTempAtUA0, 2 ) );
						ShowContinueError( "Calculated water outlet temperature at high UA [C](UA = " + RoundSigDigits( UA1, 2 ) + " W/C)  = " + RoundSigDigits( OutWaterTempAtUA1, 2 ) );
						ShowFatalError( "Autosizing of Evaporative Fluid Cooler UA failed for Evaporative Fluid Cooler = " + this->Name );
					}
					tmpHighSpeedEvapFluidCoolerUA = UA;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedEvapFluidCoolerUA = tmpHighSpeedEvapFluidCoolerUA;
					this->HighSpeedStandardDesignCapacity = DesEvapFluidCoolerLoad / this->HeatRejectCapNomCapSizingRatio;
				} else {
					tmpHighSpeedEvapFluidCoolerUA = 0.0;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedEvapFluidCoolerUA = tmpHighSpeedEvapFluidCoolerUA;
				}
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
						if ( PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
								"U-Factor Times Area Value at Design Air Flow Rate [W/C]", this->HighSpeedEvapFluidCoolerUA );
						}
						if ( PlantFirstSizesOkayToReport ) {
							ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
								"Initial U-Factor Times Area Value at Design Air Flow Rate [W/C]", this->HighSpeedEvapFluidCoolerUA );
						}
					} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
						if ( PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
								"U-Factor Times Area Value at High Fan Speed [W/C]", this->HighSpeedEvapFluidCoolerUA );
						}
						if ( PlantFirstSizesOkayToReport ) {
							ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
								"Initial U-Factor Times Area Value at High Fan Speed [W/C]", this->HighSpeedEvapFluidCoolerUA );
						}
					}
				}
			} else {
				if ( PlantFirstSizesOkayToFinalize ) {
					ShowSevereError( "Autosizing error for evaporative fluid cooler object = " + this->Name );
					ShowFatalError( "Autosizing of evaporative fluid cooler UA requires a loop Sizing:Plant object." );
				}
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::StandardDesignCapacity ) {
			if ( this->DesignWaterFlowRate >= SmallWaterVolFlow ) {
				// Standard Design Capacity doesn't include compressor heat;
				// predefined factor was 1.25 W heat rejection per W of delivered cooling, now a user input with 1.25 default
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, 35.0, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesEvapFluidCoolerLoad = this->HighSpeedStandardDesignCapacity * this->HeatRejectCapNomCapSizingRatio;
				Par( 1 ) = DesEvapFluidCoolerLoad;
				Par( 2 ) = rho * this->DesignWaterFlowRate; // Design water mass flow rate
				Par( 3 ) = this->HighSpeedAirFlowRate; // Design air volume flow rate
				Par( 4 ) = Cp;
				UA0 = 0.0001 * DesEvapFluidCoolerLoad; // Assume deltaT = 10000K (limit)
				UA1 = DesEvapFluidCoolerLoad; // Assume deltaT = 1K
				this->WaterTemp = 35.0; // 95F design inlet water temperature
				this->AirTemp = 35.0; // 95F design inlet air dry-bulb temp
				this->AirWetBulb = 25.6; // 78F design inlet air wet-bulb temp
				this->AirPress = StdBaroPress;
				this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress );
				SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleEvapFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
				if ( SolFla == -1 ) {
					ShowWarningError( "Iteration limit exceeded in calculating evaporative fluid cooler UA." );
					ShowContinueError( "Autosizing of fluid cooler UA failed for evaporative fluid cooler = " + this->Name );
					ShowContinueError( "The final UA value = " + RoundSigDigits( UA, 2 ) + "W/C, and the simulation continues..." );
				} else if ( SolFla == -2 ) {
					ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
					ShowContinueError( "reasonable UA value. Review and revise design input values as appropriate. " );
					ShowFatalError( "Autosizing of Evaporative Fluid Cooler UA failed for Evaporative Fluid Cooler = " + this->Name );
				}
				this->HighSpeedEvapFluidCoolerUA = UA;
			} else {
				this->HighSpeedEvapFluidCoolerUA = 0.0;
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"U-Factor Times Area Value at Design Air Flow Rate [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Initial U-Factor Times Area Value at Design Air Flow Rate [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
				} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"U-Factor Times Area Value at High Fan Speed [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Initial U-Factor Times Area Value at High Fan Speed [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
				}
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::UserSpecifiedDesignCapacity ) {
			if ( this->DesignWaterFlowRate >= SmallWaterVolFlow ) {
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, this->DesignEnteringWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesEvapFluidCoolerLoad = this->HighSpeedUserSpecifiedDesignCapacity;
				Par( 1 ) = DesEvapFluidCoolerLoad;
				Par( 2 ) = rho * tmpDesignWaterFlowRate; // Design water mass flow rate
				Par( 3 ) = tmpHighSpeedAirFlowRate; // Design air volume flow rate
				Par( 4 ) = Cp;
				UA0 = 0.0001 * DesEvapFluidCoolerLoad; // Assume deltaT = 10000K (limit)
				UA1 = DesEvapFluidCoolerLoad; // Assume deltaT = 1K

				this->WaterTemp = this->DesignEnteringWaterTemp;
				this->AirTemp = this->DesignEnteringAirTemp;
				this->AirWetBulb = this->DesignEnteringAirWetBulbTemp;
				this->AirPress = StdBaroPress;
				this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress );
				SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleEvapFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
				if ( SolFla == -1 ) {
					ShowWarningError( "Iteration limit exceeded in calculating evaporative fluid cooler UA." );
					ShowContinueError( "Autosizing of fluid cooler UA failed for evaporative fluid cooler = " + this->Name );
					ShowContinueError( "The final UA value = " + RoundSigDigits( UA, 2 ) + "W/C, and the simulation continues..." );
				} else if ( SolFla == -2 ) {
					this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA0, OutWaterTempAtUA0 );
					this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA1, OutWaterTempAtUA1 );
					ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
					ShowContinueError( "reasonable UA value. Review and revise design input values as appropriate. Specifying hard" );
					ShowContinueError( "sizes for some \"autosizable\" fields while autosizing other \"autosizable\" fields may be contributing to this problem." );
					ShowContinueError( "This model iterates on UA to find the heat transfer required to provide the design outlet " );
					ShowContinueError( "water temperature. Initially, the outlet water temperatures at high and low UA values are " );
					ShowContinueError( "calculated. The Design Exit Water Temperature should be between the outlet water " );
					ShowContinueError( "temperatures calculated at high and low UA values. If the Design Exit Water Temperature is " );
					ShowContinueError( "out of this range, the solution will not converge and UA will not be calculated. " );
					ShowContinueError( "The possible solutions could be to manually input adjusted water and/or air flow rates " );
					ShowContinueError( "based on the autosized values shown below or to adjust design evaporative fluid cooler air inlet wet-bulb temperature." );
					ShowContinueError( "Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp)." );
					ShowContinueError( "Inputs to the evaporative fluid cooler object:" );
					ShowContinueError( "Design Evaporative Fluid Cooler Load [W]                      = " + RoundSigDigits( Par( 1 ), 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Water Volume Flow Rate [m3/s] = " + RoundSigDigits( this->DesignWaterFlowRate, 6 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Air Volume Flow Rate [m3/s]   = " + RoundSigDigits( Par( 3 ), 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Air Inlet Wet-bulb Temp [C]   = " + RoundSigDigits( this->AirWetBulb, 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Water Inlet Temp [C]          = " + RoundSigDigits( this->WaterTemp, 2 ) );
					ShowContinueError( "Inputs to the plant sizing object:" );
					ShowContinueError( "Design Exit Water Temp [C]                                    = " + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) );
					ShowContinueError( "Loop Design Temperature Difference [C]                        = " + RoundSigDigits( PlantSizData( PltSizCondNum ).DeltaT, 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Water Inlet Temp [C]          = " + RoundSigDigits( this->WaterTemp, 2 ) );
					ShowContinueError( "Calculated water outlet temperature at low UA [C](UA = " + RoundSigDigits( UA0, 2 ) + " W/C)  = " + RoundSigDigits( OutWaterTempAtUA0, 2 ) );
					ShowContinueError( "Calculated water outlet temperature at high UA [C](UA = " + RoundSigDigits( UA1, 2 ) + " W/C)  = " + RoundSigDigits( OutWaterTempAtUA1, 2 ) );
					ShowFatalError( "Autosizing of Evaporative Fluid Cooler UA failed for Evaporative Fluid Cooler = " + this->Name );
				}
				this->HighSpeedEvapFluidCoolerUA = UA;
			} else {
				this->HighSpeedEvapFluidCoolerUA = 0.0;
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::SingleSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"U-Factor Times Area Value at Design Air Flow Rate [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_SingleSpeed, this->Name,
							"Initial U-Factor Times Area Value at Design Air Flow Rate [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
				} else if ( this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"U-Factor Times Area Value at High Fan Speed [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( cEvapFluidCooler_TwoSpeed, this->Name,
							"Initial U-Factor Times Area Value at High Fan Speed [W/C]", this->HighSpeedEvapFluidCoolerUA );
					}
				}
			}
		}

		if ( this->LowSpeedAirFlowRateWasAutoSized && PlantFirstSizesOkayToFinalize ) {
			this->LowSpeedAirFlowRate = this->LowSpeedAirFlowRateSizingFactor * this->HighSpeedAirFlowRate;
			if ( PlantFinalSizesOkayToReport ) {
				ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
					"Air Flow Rate at Low Fan Speed [m3/s]", this->LowSpeedAirFlowRate );
			}
			if ( PlantFirstSizesOkayToReport ) {
				ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
					"Initial Air Flow Rate at Low Fan Speed [m3/s]", this->LowSpeedAirFlowRate );
			}
		}

		if ( this->LowSpeedFanPowerWasAutoSized && PlantFirstSizesOkayToFinalize ) {
			this->LowSpeedFanPower = this->LowSpeedFanPowerSizingFactor * this->HighSpeedFanPower;
			if ( PlantFinalSizesOkayToReport ) {
				ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
					"Fan Power at Low Fan Speed [W]", this->LowSpeedFanPower );
			}
			if ( PlantFirstSizesOkayToReport ) {
				ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
					"Initial Fan Power at Low Fan Speed [W]", this->LowSpeedFanPower );
			}
		}

		if ( this->LowSpeedEvapFluidCoolerUAWasAutoSized && PlantFirstSizesOkayToFinalize ) {
			this->LowSpeedEvapFluidCoolerUA = this->LowSpeedEvapFluidCoolerUASizingFactor * this->HighSpeedEvapFluidCoolerUA;
			if ( PlantFinalSizesOkayToReport ) {
				ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
					"U-Factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedEvapFluidCoolerUA );
			}
			if ( PlantFirstSizesOkayToReport ) {
				ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
					"Initial U-Factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedEvapFluidCoolerUA );
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::StandardDesignCapacity
				&& this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
			if ( this->DesignWaterFlowRate >= SmallWaterVolFlow
					&& this->LowSpeedStandardDesignCapacity > 0.0 ) {
				// Standard design capacity doesn't include compressor heat;
				// predefined factor was 1.25 W heat rejection per W of delivered cooling, now user input with default 1.25
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, this->DesignEnteringWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesEvapFluidCoolerLoad = this->LowSpeedStandardDesignCapacity * this->HeatRejectCapNomCapSizingRatio;
				Par( 1 ) = DesEvapFluidCoolerLoad;
				Par( 2 ) = rho * tmpDesignWaterFlowRate; // Design water mass flow rate
				Par( 3 ) = this->LowSpeedAirFlowRate; // Air volume flow rate at low fan speed
				Par( 4 ) = Cp;
				UA0 = 0.0001 * DesEvapFluidCoolerLoad; // Assume deltaT = 10000K (limit)
				UA1 = DesEvapFluidCoolerLoad; // Assume deltaT = 1K
				this->WaterTemp = 35.0; // 95F design inlet water temperature
				this->AirTemp = 35.0; // 95F design inlet air dry-bulb temp
				this->AirWetBulb = 25.6; // 78F design inlet air wet-bulb temp
				this->AirPress = StdBaroPress;
				this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress );
				SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleEvapFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
				if ( SolFla == -1 ) {
					ShowWarningError( "Iteration limit exceeded in calculating evaporative fluid cooler UA." );
					ShowContinueError( "Autosizing of fluid cooler UA failed for evaporative fluid cooler = " + this->Name );
					ShowContinueError( "The final UA value = " + RoundSigDigits( UA, 2 ) + "W/C, and the simulation continues..." );
				} else if ( SolFla == -2 ) {
					ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
					ShowContinueError( "reasonable low-speed UA value. Review and revise design input values as appropriate. " );
					ShowFatalError( "Autosizing of Evaporative Fluid Cooler UA failed for Evaporative Fluid Cooler = " + this->Name );
				}
				this->LowSpeedEvapFluidCoolerUA = UA;
			} else {
				this->LowSpeedEvapFluidCoolerUA = 0.0;
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( PlantFinalSizesOkayToReport ) {
					ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
						"U-Factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedEvapFluidCoolerUA );
				}
				if ( PlantFirstSizesOkayToReport ) {
					ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
						"Initial U-Factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedEvapFluidCoolerUA );
				}
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::UserSpecifiedDesignCapacity
				&& this->EvapFluidCoolerType_Num == EvaporativeFluidCoolerType::TwoSpeed ) {
			if ( this->DesignWaterFlowRate >= SmallWaterVolFlow
					&& this->LowSpeedUserSpecifiedDesignCapacity > 0.0 ) {
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, this->DesignEnteringWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesEvapFluidCoolerLoad = this->LowSpeedUserSpecifiedDesignCapacity;
				Par( 1 ) = DesEvapFluidCoolerLoad;
				Par( 2 ) = rho * tmpDesignWaterFlowRate; // Design water mass flow rate
				Par( 3 ) = this->LowSpeedAirFlowRate; // Air volume flow rate at low fan speed
				Par( 4 ) = Cp;
				UA0 = 0.0001 * DesEvapFluidCoolerLoad; // Assume deltaT = 10000K (limit)
				UA1 = DesEvapFluidCoolerLoad; // Assume deltaT = 1K
				this->WaterTemp = this->DesignEnteringWaterTemp;
				this->AirTemp = this->DesignEnteringAirTemp;
				this->AirWetBulb = this->DesignEnteringAirWetBulbTemp;
				this->AirPress = StdBaroPress;
				this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress );
				SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleEvapFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
				if ( SolFla == -1 ) {
					ShowSevereError( "Iteration limit exceeded in calculating EvaporativeFluidCooler UA" );
					ShowFatalError( "Autosizing of EvaporativeFluidCooler UA failed for EvaporativeFluidCooler " + this->Name );
				} else if ( SolFla == -2 ) {
					this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA0, OutWaterTempAtUA0 );
					this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA1, OutWaterTempAtUA1 );
					ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
					ShowContinueError( "reasonable UA value. Review and revise design input values as appropriate. Specifying hard" );
					ShowContinueError( "sizes for some \"autosizable\" fields while autosizing other \"autosizable\" fields may be contributing to this problem." );
					ShowContinueError( "This model iterates on UA to find the heat transfer required to provide the design outlet " );
					ShowContinueError( "water temperature. Initially, the outlet water temperatures at high and low UA values are " );
					ShowContinueError( "calculated. The Design Exit Water Temperature should be between the outlet water " );
					ShowContinueError( "temperatures calculated at high and low UA values. If the Design Exit Water Temperature is " );
					ShowContinueError( "out of this range, the solution will not converge and UA will not be calculated. " );
					ShowContinueError( "Inputs to the Evaporative Fluid Cooler model are:" );
					ShowContinueError( "Design Evaporative Fluid Cooler Load                    = " + RoundSigDigits( Par( 1 ), 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Water Volume Flow Rate  = " + RoundSigDigits( Par( 3 ), 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Air Volume Flow Rate    = " + RoundSigDigits( Par( 3 ), 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Air Inlet Wet-bulb Temp = " + RoundSigDigits( this->AirWetBulb, 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Water Inlet Temp        = " + RoundSigDigits( this->WaterTemp, 2 ) );
					ShowContinueError( "Design Exit Water Temp                                  = " + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) );
					ShowContinueError( "Design Evaporative Fluid Cooler Water Inlet Temp [C]    = " + RoundSigDigits( this->WaterTemp, 2 ) );
					ShowContinueError( "Calculated water outlet temperature at low UA(" + RoundSigDigits( UA0, 2 ) + ")  = " + RoundSigDigits( OutWaterTempAtUA0, 2 ) );
					ShowContinueError( "Calculated water outlet temperature at high UA(" + RoundSigDigits( UA1, 2 ) + ")  = " + RoundSigDigits( OutWaterTempAtUA1, 2 ) );
					ShowFatalError( "Autosizing of Evaporative Fluid Cooler UA failed for Evaporative Fluid Cooler = " + this->Name );
				}
				this->LowSpeedEvapFluidCoolerUA = UA;
			} else {
				this->LowSpeedEvapFluidCoolerUA = 0.0;
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( PlantFinalSizesOkayToReport ) {
					ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
						"U-Factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedEvapFluidCoolerUA );
				}
				if ( PlantFirstSizesOkayToReport ) {
					ReportSizingOutput( this->EvapFluidCoolerType, this->Name,
						"Initial U-Factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedEvapFluidCoolerUA );
				}

			}
		}

		if ( PlantFinalSizesOkayToReport ) {
			//create predefined report
			equipName = this->Name;
			PreDefTableEntry( pdchMechType, equipName, this->EvapFluidCoolerType );
			PreDefTableEntry( pdchMechNomCap, equipName, this->HighSpeedStandardDesignCapacity );
		}

	}

	void
	EvaporativeFluidCooler::calcSingleSpeedEvapFluidCooler()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// To simulate the operation of a single-speed fan evaporative fluid cooler.

		// METHODOLOGY EMPLOYED:
		// The evaporative fluid cooler is modeled using effectiveness-NTU relationships for
		// counterflow heat exchangers based on Merkel's theory.
		// The subroutine calculates the period of time required to meet a
		// leaving water temperature setpoint. It assumes that part-load
		// operation represents a linear interpolation of two steady-state regimes.
		// Cyclic losses are neglected. The period of time required to meet the
		// leaving water temperature setpoint is used to determine the required
		// fan power and energy.
		// A RunFlag is passed by the upper level manager to indicate the ON/OFF status,
		// or schedule, of the evaporative fluid cooler. If the evaporative fluid cooler is OFF, outlet water
		// temperature and flow rate are passed through the model from inlet node to
		// outlet node without intervention. Reports are also updated with fan power and energy being zero.
		// When the RunFlag indicates an ON condition for the evaporative fluid cooler, the
		// mass flow rate and water temperature are read from the inlet node of the
		// evaporative fluid cooler (water-side). The outdoor air wet-bulb temperature is used
		// as the entering condition to the evaporative fluid cooler (air-side).
		// The evaporative fluid cooler fan is turned on and design parameters are used
		// to calculate the leaving water temperature.
		// If the calculated leaving water temperature is below the setpoint, a fan
		// run-time fraction is calculated and used to determine fan power. The leaving
		// water temperature setpoint is placed on the outlet node. If the calculated
		// leaving water temperature is at or above the setpoint, the calculated
		// leaving water temperature is placed on the outlet node and the fan runs at
		// full power. Water mass flow rate is passed from inlet node to outlet node
		// with no intervention.
		// REFERENCES:
		// ASHRAE HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculation. 1999.

		// Based on SingleSpeedTower subroutine by Dan Fisher ,Sept 1998
		// Dec. 2008. BG. added RunFlag logic per original methodology

		// USE STATEMENTS:
		//  USE FluidProperties, ONLY : GetSpecificHeatGlycol
		// Using/Aliasing
		using DataPlant::PlantLoop;
		using DataPlant::SingleSetPoint;
		using DataPlant::DualSetPointDeadBand;

		// Locals
		Real64 InletWaterTemp;

		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "CalcSingleSpeedEvapFluidCooler" );
		int const MaxIteration( 100 ); // Maximum fluid bypass iteration calculations
		static std::string const MaxItChar( "100" );
		Real64 const BypassFractionThreshold( 0.01 ); // Threshold to stop bypass iteration
		Real64 const OWTLowerLimit( 0.0 ); // The limit of evaporative fluid cooler exit fluid temperature used
		// in the fluid bypass calculation to avoid fluid freezing. For water,
		// it is 0 degreeC and for glycols, it can be much lower. The fluid type
		// is stored at the loop. Current choices are Water and Steam,
		// needs to expand for glycols

		Real64 AirFlowRate;
		Real64 UAdesign; // UA value at design conditions (entered by user or calculated)
		Real64 FanModeFrac;
		Real64 FanPowerOn;
		Real64 CpWater;
		Real64 TempSetPoint;

		//Added variables for fluid bypass
		int NumIteration;
		int CapacityControl; // Capacity Control (0 - FanCycling, 1 - FluidBypass)
		int BypassFlag; // Flag indicator for fluid bypass (0 - no bypass, 1 - bypass)
		Real64 BypassFraction; // Fluid bypass fraction
		Real64 BypassFraction2; // Fluid bypass fraction
		Real64 BypassFractionPrev;
		Real64 OutletWaterTempPrev;
		int LoopNum;
		int LoopSideNum;

		//set inlet and outlet nodes
		Qactual = 0.0;
		FanPower = 0.0;
		InletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
		OutletWaterTemp = InletWaterTemp;
		LoopNum = this->location.loopNum;
		LoopSideNum = this->location.loopSideNum;
		AirFlowRate = 0.0;
		{ auto const SELECT_CASE_var( PlantLoop( LoopNum ).LoopDemandCalcScheme );
		if ( SELECT_CASE_var == SingleSetPoint ) {
			TempSetPoint = PlantLoop( LoopNum ).LoopSide( LoopSideNum ).TempSetPoint;
		} else if ( SELECT_CASE_var == DualSetPointDeadBand ) {
			TempSetPoint = PlantLoop( LoopNum ).LoopSide( LoopSideNum ).TempSetPointHi;
		}}

		// Added for fluid bypass. First assume no fluid bypass
		BypassFlag = 0;
		BypassFraction = 0.0;
		BypassFraction2 = 0.0;
		this->BypassFraction = 0.0;
		CapacityControl = this->CapacityControl;

		//   MassFlowTol is a parameter to indicate a no flow condition
		if ( WaterMassFlowRate <= MassFlowTolerance || PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 0 ) return;

		if ( InletWaterTemp > TempSetPoint ) {
			//     Turn on evaporative fluid cooler fan
			UAdesign = this->HighSpeedEvapFluidCoolerUA;
			AirFlowRate = this->HighSpeedAirFlowRate;
			FanPowerOn = this->HighSpeedFanPower;

			this->simSimpleEvapFluidCooler( WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp );

			if ( OutletWaterTemp <= TempSetPoint ) {
				if ( CapacityControl == 0 || OutletWaterTemp <= OWTLowerLimit ) {
					//         Setpoint was met with pump ON and fan ON, calculate run-time fraction
					FanModeFrac = ( TempSetPoint - InletWaterTemp ) / ( OutletWaterTemp - InletWaterTemp );
					FanPower = FanModeFrac * FanPowerOn;
					OutletWaterTemp = TempSetPoint;
				} else {
					//FluidBypass, fan runs at full speed for the entire time step
					FanModeFrac = 1.0;
					FanPower = FanPowerOn;
					BypassFlag = 1;
				}
			} else {
				//       Setpoint was not met, evaporative fluid cooler ran at full capacity
				FanModeFrac = 1.0;
				FanPower = FanPowerOn;
			}
		} else if ( InletWaterTemp <= TempSetPoint ) {
			//Inlet water temperature lower than setpoint, assume 100% bypass, evaporative fluid cooler fan off
			if ( CapacityControl == 1 ) {
				if ( InletWaterTemp > OWTLowerLimit ) {
					FanPower = 0.0;
					BypassFraction = 1.0;
					this->BypassFraction = 1.0;
					OutletWaterTemp = InletWaterTemp;
				}
			}
		}

		// Calculate bypass fraction since OWTLowerLimit < OutletWaterTemp < TempSetPoint.
		// The iteration ends when the numer of iteration exceeds the limit or the difference
		//  between the new and old bypass fractions is less than the threshold.
		if ( BypassFlag == 1 ) {
			BypassFraction = ( TempSetPoint - OutletWaterTemp ) / ( InletWaterTemp - OutletWaterTemp );
			if ( BypassFraction > 1.0 || BypassFraction < 0.0 ) {
				// Bypass cannot meet setpoint, assume no bypass
				BypassFlag = 0;
				BypassFraction = 0.0;
				this->BypassFraction = 0.0;
				AirFlowRate = 0.0;
			} else {
				NumIteration = 0;
				BypassFractionPrev = BypassFraction;
				OutletWaterTempPrev = OutletWaterTemp;
				while ( NumIteration < MaxIteration ) {
					++NumIteration;
					// need to iterate for the new OutletWaterTemp while bypassing evaporative fluid cooler water
					this->simSimpleEvapFluidCooler( WaterMassFlowRate * ( 1.0 - BypassFraction ), AirFlowRate, UAdesign, OutletWaterTemp );
					// Calc new BypassFraction based on the new OutletWaterTemp
					if ( std::abs( OutletWaterTemp - OWTLowerLimit ) <= 0.01 ) {
						BypassFraction2 = BypassFraction;
						break;
					} else if ( OutletWaterTemp < OWTLowerLimit ) {
						// Set OutletWaterTemp = OWTLowerLimit, and use linear interpolation to calculate the bypassFraction
						BypassFraction2 = BypassFractionPrev - ( BypassFractionPrev - BypassFraction ) * ( OutletWaterTempPrev - OWTLowerLimit ) / ( OutletWaterTempPrev - OutletWaterTemp );
						this->simSimpleEvapFluidCooler( WaterMassFlowRate * ( 1.0 - BypassFraction2 ), AirFlowRate, UAdesign, OutletWaterTemp );
						if ( OutletWaterTemp < OWTLowerLimit ) {
							//Use previous iteraction values
							BypassFraction2 = BypassFractionPrev;
							OutletWaterTemp = OutletWaterTempPrev;
						}
						break;
					} else {
						BypassFraction2 = ( TempSetPoint - OutletWaterTemp ) / ( InletWaterTemp - OutletWaterTemp );
					}
					// Compare two BypassFraction to determine when to stop
					if ( std::abs( BypassFraction2 - BypassFraction ) <= BypassFractionThreshold ) break;
					BypassFractionPrev = BypassFraction;
					OutletWaterTempPrev = OutletWaterTemp;
					BypassFraction = BypassFraction2;
				}
				if ( NumIteration > MaxIteration ) {
					ShowWarningError( "Evaporative fluid cooler fluid bypass iteration exceeds maximum limit of " + MaxItChar + " for " + this->Name );
				}
				this->BypassFraction = BypassFraction2;
				// may not meet TempSetPoint due to limit of evaporative fluid cooler outlet temp to OWTLowerLimit
				OutletWaterTemp = ( 1.0 - BypassFraction2 ) * OutletWaterTemp + BypassFraction2 * InletWaterTemp;
			}
		}

		//Should this be water inlet node num?????
		CpWater = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, Node( this->WaterInletNodeNum ).Temp, PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
		Qactual = WaterMassFlowRate * CpWater * ( Node( this->WaterInletNodeNum ).Temp - OutletWaterTemp );
		this->AirFlowRatio = AirFlowRate / this->HighSpeedAirFlowRate;

	}

	void
	EvaporativeFluidCooler::calcTwoSpeedEvapFluidCooler()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// To simulate the operation of a evaporative fluid cooler with a two-speed fan.

		// METHODOLOGY EMPLOYED:
		// The evaporative fluid cooler is modeled using effectiveness-NTU relationships for
		// counterflow heat exchangers based on Merkel's theory.
		// The subroutine calculates the period of time required to meet a
		// leaving water temperature setpoint. It assumes that part-load
		// operation represents a linear interpolation of three steady-state regimes
		// (high-speed fan operation and low-speed fan operation ).
		// Cyclic losses are neglected. The period of time required to meet the
		// leaving water temperature setpoint is used to determine the required
		// fan power and energy. When the leaving water temperature is at or above the setpoint
		// the evaporative fluid cooler fan is turned on,
		// .
		// A RunFlag is passed by the upper level manager to indicate the ON/OFF status,
		// or schedule, of the evaporative fluid cooler. If the evaporative fluid cooler is OFF, outlet water
		// temperature and flow rate are passed through the model from inlet node to
		// outlet node without intervention. Reports are also updated with fan power and fan energy being zero.
		// When the RunFlag indicates an ON condition for the evaporative fluid cooler, the
		// mass flow rate and water temperature are read from the inlet node of the
		// evaporative fluid cooler (water-side). The outdoor air wet-bulb temperature is used
		// as the entering condition to the evaporative fluid cooler (air-side). If the incoming
		// water temperature is above the setpoint, the evaporative fluid cooler fan is turned on
		// and parameters for low fan speed are used to again calculate the leaving
		// water temperature. If the calculated leaving water temperature is
		// below the setpoint, a fan run-time fraction (FanModeFrac) is calculated and
		// used to determine fan power. The leaving water temperature setpoint is placed
		// on the outlet node. If the calculated leaving water temperature is at or above
		// the setpoint, the evaporative fluid cooler fan is turned on 'high speed' and the routine is
		// repeated. If the calculated leaving water temperature is below the setpoint,
		// a fan run-time fraction is calculated for the second stage fan and fan power
		// is calculated as FanModeFrac*HighSpeedFanPower+(1-FanModeFrac)*LowSpeedFanPower.
		// If the calculated leaving water temperature is above the leaving water temp.
		// setpoint, the calculated leaving water temperature is placed on the outlet
		// node and the fan runs at full power (High Speed Fan Power). Water mass flow
		// rate is passed from inlet node to outlet node with no intervention.
		// REFERENCES:
		// ASHRAE HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculation. 1999.
		// Based on TwoSpeedTower by Dan Fisher ,Sept. 1998
		// Dec. 2008. BG. added RunFlag logic per original methodology

		using DataPlant::SingleSetPoint;
		using DataPlant::DualSetPointDeadBand;

		// Locals
		Real64 InletWaterTemp;

		static std::string const RoutineName( "CalcTwoSpeedEvapFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 AirFlowRate;
		Real64 UAdesign; // UA value at design conditions (entered by user) [W/C]
		Real64 OutletWaterTemp1stStage;
		Real64 OutletWaterTemp2ndStage;
		Real64 FanModeFrac;
		Real64 FanPowerLow;
		Real64 FanPowerHigh;
		Real64 CpWater;
		Real64 TempSetPoint;
		int LoopNum;
		int LoopSideNum;

		Qactual = 0.0;
		FanPower = 0.0;
		InletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
		OutletWaterTemp = InletWaterTemp;

		OutletWaterTemp1stStage = OutletWaterTemp;
		OutletWaterTemp2ndStage = OutletWaterTemp;
		FanModeFrac = 0.0;
		AirFlowRate = 0.0;
		LoopNum = this->location.loopNum;
		LoopSideNum = this->location.loopSideNum;
		{ auto const SELECT_CASE_var( PlantLoop( LoopNum ).LoopDemandCalcScheme );
		if ( SELECT_CASE_var == SingleSetPoint ) {
			TempSetPoint = PlantLoop( LoopNum ).LoopSide( LoopSideNum ).TempSetPoint;
		} else if ( SELECT_CASE_var == DualSetPointDeadBand ) {
			TempSetPoint = PlantLoop( LoopNum ).LoopSide( LoopSideNum ).TempSetPointHi;
		}}

		//   MassFlowTol is a parameter to indicate a no flow condition
		if ( WaterMassFlowRate <= MassFlowTolerance || PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 0 ) return;

		if ( InletWaterTemp > TempSetPoint ) {
			//     Setpoint was not met ,turn on evaporative fluid cooler 1st stage fan
			UAdesign = this->LowSpeedEvapFluidCoolerUA;
			AirFlowRate = this->LowSpeedAirFlowRate;
			FanPowerLow = this->LowSpeedFanPower;

			this->simSimpleEvapFluidCooler( WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp1stStage );

			if ( OutletWaterTemp1stStage <= TempSetPoint ) {
				//         Setpoint was met with pump ON and fan ON 1st stage, calculate fan mode fraction
				FanModeFrac = ( TempSetPoint - InletWaterTemp ) / ( OutletWaterTemp1stStage - InletWaterTemp );
				FanPower = FanModeFrac * FanPowerLow;
				OutletWaterTemp = TempSetPoint;
				Qactual *= FanModeFrac;
			} else {
				//         Setpoint was not met, turn on evaporative fluid cooler 2nd stage fan
				UAdesign = this->HighSpeedEvapFluidCoolerUA;
				AirFlowRate = this->HighSpeedAirFlowRate;
				FanPowerHigh = this->HighSpeedFanPower;

				this->simSimpleEvapFluidCooler( WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp2ndStage );

				if ( ( OutletWaterTemp2ndStage <= TempSetPoint ) && UAdesign > 0.0 ) {
					//           Setpoint was met with pump ON and fan ON 2nd stage, calculate fan mode fraction
					FanModeFrac = ( TempSetPoint - OutletWaterTemp1stStage ) / ( OutletWaterTemp2ndStage - OutletWaterTemp1stStage );
					FanPower = ( FanModeFrac * FanPowerHigh ) + ( 1.0 - FanModeFrac ) * FanPowerLow;
					OutletWaterTemp = TempSetPoint;
				} else {
					//           Setpoint was not met, evaporative fluid cooler ran at full capacity
					OutletWaterTemp = OutletWaterTemp2ndStage;
					FanPower = FanPowerHigh;
				}

			}

		}

		//Should this be water inlet node num??
		CpWater = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, Node( this->WaterInletNodeNum ).Temp, PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
		Qactual = WaterMassFlowRate * CpWater * ( Node( this->WaterInletNodeNum ).Temp - OutletWaterTemp );
		this->AirFlowRatio = AirFlowRate / this->HighSpeedAirFlowRate;

	}

	void
	EvaporativeFluidCooler::simSimpleEvapFluidCooler(
		Real64 const WaterMassFlowRate,
		Real64 const AirFlowRate,
		Real64 const UAdesign,
		Real64 & OutletWaterTemp
	)
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// See purpose for single speed or Two speed evaporative fluid cooler model

		// METHODOLOGY EMPLOYED:
		// See methodology for single speed or two speed evaporative fluid cooler model

		// REFERENCES:
		// Based on SimTower subroutine by Dan Fisher Sept. 1998
		// Merkel, F. 1925.  Verduftungskuhlung. VDI Forschungsarbeiten, Nr 275, Berlin.
		// ASHRAE     1999.  HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculations.

		// Locals
		Real64 Qactual; // Actual heat transfer rate between evaporative fluid cooler water and air [W]

		int const IterMax( 50 ); // Maximum number of iterations allowed
		Real64 const WetBulbTolerance( 0.00001 ); // Maximum error for exiting wet-bulb temperature between iterations
		// [delta K/K]
		Real64 const DeltaTwbTolerance( 0.001 ); // Maximum error (tolerance) in DeltaTwb for iteration convergence [C]
		static std::string const RoutineName( "SimSimpleEvapFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int Iter; // Number of iterations completed
		Real64 MdotCpWater; // Water mass flow rate times the heat capacity [W/K]
		Real64 InletAirTemp; // Dry-bulb temperature of air entering the evaporative fluid cooler [C]
		Real64 CpWater; // Heat capacity of water [J/kg/K]
		Real64 CpAir; // Heat capacity of air [J/kg/K]
		Real64 AirDensity; // Density of air [kg/m3]
		Real64 AirMassFlowRate; // Mass flow rate of air [kg/s]
		Real64 effectiveness; // Effectiveness of the heat exchanger [-]
		Real64 UAactual; // UA value at actual conditions [W/C]
		Real64 InletAirEnthalpy; // Enthalpy of entering moist air [J/kg]
		Real64 InletAirWetBulb; // Wetbulb temp of entering moist air [C]
		Real64 OutletAirEnthalpy; // Enthalpy of exiting moist air [J/kg]
		Real64 OutletAirWetBulb; // Wetbulb temp of exiting moist air [C]
		Real64 OutletAirWetBulbLast; // temporary Wetbulb temp of exiting moist air [C]
		Real64 AirCapacity; // MdotCp of air through the evaporative fluid cooler
		Real64 CapacityRatioMin; // Minimum capacity of airside and waterside
		Real64 CapacityRatioMax; // Maximum capacity of airside and waterside
		Real64 CapacityRatio; // Ratio of minimum to maximum capacity
		Real64 NumTransferUnits; // Number of transfer Units [NTU]
		Real64 WetBulbError; // Calculated error for exiting wet-bulb temperature between iterations [delta K/K]
		Real64 CpAirside; // Delta enthalpy of the evaporative fluid cooler air /
		// delta air wet-bulb temp [J/kg/K]
		Real64 DeltaTwb; // Absolute value of difference between inlet and outlet air wet-bulb temp [C]

		// set inlet and outlet node numbers, and initialize some local variables

		Qactual = 0.0;
		WetBulbError = 1.0;
		DeltaTwb = 1.0;

		// set local evaporative fluid cooler inlet and outlet temperature variables
		InletWaterTemp = this->WaterTemp;
		OutletWaterTemp = InletWaterTemp;
		InletAirTemp = this->AirTemp;
		InletAirWetBulb = this->AirWetBulb;

		if ( UAdesign == 0.0 ) return;

		// set water and air properties
		AirDensity = PsyRhoAirFnPbTdbW( this->AirPress, InletAirTemp, this->AirHumRat );
		AirMassFlowRate = AirFlowRate * AirDensity;
		CpAir = Psychrometrics::PsyCpAirFnWTdb( this->AirHumRat, InletAirTemp );
		CpWater = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, InletWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
		InletAirEnthalpy = PsyHFnTdbRhPb( InletAirWetBulb, 1.0, this->AirPress );

		// initialize exiting wet bulb temperature before iterating on final solution
		OutletAirWetBulb = InletAirWetBulb + 6.0;

		// Calcluate mass flow rates
		MdotCpWater = WaterMassFlowRate * CpWater;
		Iter = 0;
		while ( ( WetBulbError > WetBulbTolerance ) && ( Iter <= IterMax ) && ( DeltaTwb > DeltaTwbTolerance ) ) {
			++Iter;
			OutletAirEnthalpy = PsyHFnTdbRhPb( OutletAirWetBulb, 1.0, this->AirPress );
			// calculate the airside specific heat and capacity
			CpAirside = ( OutletAirEnthalpy - InletAirEnthalpy ) / ( OutletAirWetBulb - InletAirWetBulb );
			AirCapacity = AirMassFlowRate * CpAirside;
			// calculate the minimum to maximum capacity ratios of airside and waterside
			CapacityRatioMin = min( AirCapacity, MdotCpWater );
			CapacityRatioMax = max( AirCapacity, MdotCpWater );
			CapacityRatio = CapacityRatioMin / CapacityRatioMax;
			// Calculate heat transfer coefficient and number of transfer units (NTU)
			UAactual = UAdesign * CpAirside / CpAir;
			NumTransferUnits = UAactual / CapacityRatioMin;
			// calculate heat exchanger effectiveness
			if ( CapacityRatio <= 0.995 ) {
				effectiveness = ( 1.0 - std::exp( -1.0 * NumTransferUnits * ( 1.0 - CapacityRatio ) ) ) / ( 1.0 - CapacityRatio * std::exp( -1.0 * NumTransferUnits * ( 1.0 - CapacityRatio ) ) );
			} else {
				effectiveness = NumTransferUnits / ( 1.0 + NumTransferUnits );
			}
			// calculate water to air heat transfer and store last exiting WB temp of air
			Qactual = effectiveness * CapacityRatioMin * ( InletWaterTemp - InletAirWetBulb );
			OutletAirWetBulbLast = OutletAirWetBulb;
			// calculate new exiting wet bulb temperature of airstream
			OutletAirWetBulb = InletAirWetBulb + Qactual / AirCapacity;
			// Check error tolerance and exit if satisfied
			DeltaTwb = std::abs( OutletAirWetBulb - InletAirWetBulb );
			// Add KelvinConv to denominator below convert OutletAirWetBulbLast to Kelvin to avoid divide by zero.
			// Wet bulb error units are delta K/K
			WetBulbError = std::abs( ( OutletAirWetBulb - OutletAirWetBulbLast ) / ( OutletAirWetBulbLast + KelvinConv ) );
		}

		if ( Qactual >= 0.0 ) {
			OutletWaterTemp = InletWaterTemp - Qactual / MdotCpWater;
		} else {
			OutletWaterTemp = InletWaterTemp;
		}

	}

	Real64
	EvaporativeFluidCooler::simpleEvapFluidCoolerUAResidual(
		Real64 const UA, // UA of evaporative fluid cooler
		Array1< Real64 > const & Par // par(1) = design evaporative fluid cooler load [W]
	)
	{

		// FUNCTION INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// Calculates residual function (Design evaporative fluid cooler load - evaporative fluid cooler cooling output)
		//                                    / Design evaporative fluid cooler load.
		// Evaporative fluid cooler Cooling Output depends on the UA which is being varied to zero the residual.

		// METHODOLOGY EMPLOYED:
		// Puts UA into the evaporative fluid cooler data structure, calls SimSimpleEvapFluidCooler, and calculates
		// the residual as defined above.

		// REFERENCES:
		// Based on SimpleTowerUAResidual by Fred Buhl, May 2002

		// Return value
		Real64 Residuum; // residual to be minimized to zero

		// Argument array dimensioning

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:
		// par(2) = Evaporative fluid cooler number
		// par(3) = design water mass flow rate [kg/s]
		// par(4) = design air volume flow rate [m3/s]
		// par(5) = water specific heat [J/(kg*C)]

		// FUNCTION LOCAL VARIABLE DECLARATIONS:
		Real64 OutWaterTemp; // outlet water temperature [C]
		Real64 CoolingOutput; // Evaporative fluid cooler cooling output [W]

		this->simSimpleEvapFluidCooler( Par( 2 ), Par( 3 ), UA, OutWaterTemp );
		CoolingOutput = Par( 4 ) * Par( 2 ) * ( this->WaterTemp - OutWaterTemp );
		Residuum = ( Par( 1 ) - CoolingOutput ) / Par( 1 );
		return Residuum;
	}

	// End of the EvaporativeFluidCoolers Module Simulation Subroutines
	// *****************************************************************************

	void
	EvaporativeFluidCooler::calculateWaterUseage()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   May 2009
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// Collect evaporative fluid cooler water useage calculations for
		// reuse by all the evaporative fluid cooler models.

		// METHODOLOGY EMPLOYED:
		// <description>

		// REFERENCES:
		// Based on CalculateWaterUseage subroutine for cooling tower by B. Griffith, August 2006

		using DataGlobals::SecInHour;
		using DataHVACGlobals::TimeStepSys;
		using ScheduleManager::GetCurrentScheduleValue;
		using DataWater::WaterStorage;

		static std::string const RoutineName( "CalculateWaterUseage" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 AirDensity;
		Real64 AirMassFlowRate;
		Real64 AvailTankVdot;
		Real64 BlowDownVdot( 0.0 );
		Real64 DriftVdot( 0.0 );
		Real64 EvapVdot( 0.0 );
		Real64 InletAirEnthalpy;
		Real64 InSpecificHumRat;
		Real64 OutSpecificHumRat;
		Real64 TairAvg;
		Real64 MakeUpVdot;
		Real64 OutletAirEnthalpy;
		Real64 OutletAirHumRatSat;
		Real64 OutletAirTSat;
		Real64 StarvedVdot;
		Real64 TankSupplyVdot;
		Real64 rho;
		Real64 AverageWaterTemp;

		AverageWaterTemp = ( InletWaterTemp + OutletWaterTemp ) / 2.0;

		// Set water and air properties
		if ( this->EvapLossMode == EvaporativeLossBy::MoistTheory ) {

			AirDensity = PsyRhoAirFnPbTdbW( this->AirPress, this->AirTemp, this->AirHumRat );
			AirMassFlowRate = this->AirFlowRatio * this->HighSpeedAirFlowRate * AirDensity;
			InletAirEnthalpy = PsyHFnTdbRhPb( this->AirWetBulb, 1.0, this->AirPress );

			if ( AirMassFlowRate > 0.0 ) {
				// Calculate outlet air conditions for determining water usage

				OutletAirEnthalpy = InletAirEnthalpy + Qactual / AirMassFlowRate;
				OutletAirTSat = Psychrometrics::PsyTsatFnHPb( OutletAirEnthalpy, this->AirPress );
				OutletAirHumRatSat = Psychrometrics::PsyWFnTdbH( OutletAirTSat, OutletAirEnthalpy );

				// calculate specific humidity ratios (HUMRAT to mass of moist air not dry air)
				InSpecificHumRat = this->AirHumRat / ( 1 + this->AirHumRat );
				OutSpecificHumRat = OutletAirHumRatSat / ( 1 + OutletAirHumRatSat );

				// calculate average air temp for density call
				TairAvg = ( this->AirTemp + OutletAirTSat ) / 2.0;

				// Amount of water evaporated
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, TairAvg, PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
				EvapVdot = ( AirMassFlowRate * ( OutSpecificHumRat - InSpecificHumRat ) ) / rho; // [m3/s]
				if ( EvapVdot < 0.0 ) EvapVdot = 0.0;
			} else {
				EvapVdot = 0.0;
			}

		} else if ( this->EvapLossMode == EvaporativeLossBy::UserFactor ) {
			rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, AverageWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
			EvapVdot = this->UserEvapLossFactor * ( InletWaterTemp - OutletWaterTemp ) * ( WaterMassFlowRate / rho );
			if ( EvapVdot < 0.0 ) EvapVdot = 0.0;
		} else {
			// should never come here
			assert( false );
		}

		//   amount of water lost due to drift
		DriftVdot = this->DesignSprayWaterFlowRate * this->DriftLossFraction * this->AirFlowRatio;

		if ( this->BlowdownMode == BlowdownBy::Schedule ) {
			// Amount of water lost due to blow down (purging contaminants from evaporative fluid cooler basin)
			if ( this->SchedIDBlowdown > 0 ) {
				BlowDownVdot = GetCurrentScheduleValue( this->SchedIDBlowdown );
			} else {
				BlowDownVdot = 0.0;
			}
		} else if ( this->BlowdownMode == BlowdownBy::Concentration ) {
			if ( this->ConcentrationRatio > 2.0 ) { // protect divide by zero
				BlowDownVdot = EvapVdot / ( this->ConcentrationRatio - 1 ) - DriftVdot;
			} else {
				BlowDownVdot = EvapVdot - DriftVdot;
			}
			if ( BlowDownVdot < 0.0 ) BlowDownVdot = 0.0;
		} else {
			//should never come here
			assert( false );
		}

		// Added for fluid bypass
		if ( this->CapacityControl == 1 ) {
			if ( this->EvapLossMode == EvaporativeLossBy::UserFactor ) EvapVdot *= ( 1 - this->BypassFraction );
			DriftVdot *= ( 1 - this->BypassFraction );
			BlowDownVdot *= ( 1 - this->BypassFraction );
		}

		MakeUpVdot = EvapVdot + DriftVdot + BlowDownVdot;

		// set demand request in Water STorage if needed
		StarvedVdot = 0.0;
		TankSupplyVdot = 0.0;
		if ( this->SuppliedByWaterSystem ) {

			// set demand request
			WaterStorage( this->WaterTankID ).VdotRequestDemand( this->WaterTankDemandARRID ) = MakeUpVdot;

			AvailTankVdot = WaterStorage( this->WaterTankID ).VdotAvailDemand( this->WaterTankDemandARRID ); // check what tank can currently provide

			TankSupplyVdot = MakeUpVdot; // init
			if ( AvailTankVdot < MakeUpVdot ) { // calculate starved flow
				StarvedVdot = MakeUpVdot - AvailTankVdot;
				TankSupplyVdot = AvailTankVdot;
			}
		} else { // supplied by mains

		}

		//   total water usage
		// update report variables
		this->EvaporationVdot = EvapVdot;
		this->EvaporationVol = EvapVdot * ( TimeStepSys * SecInHour );
		this->DriftVdot = DriftVdot;
		this->DriftVol = DriftVdot * ( TimeStepSys * SecInHour );
		this->BlowdownVdot = BlowDownVdot;
		this->BlowdownVol = BlowDownVdot * ( TimeStepSys * SecInHour );
		this->MakeUpVdot = MakeUpVdot;
		this->MakeUpVol = MakeUpVdot * ( TimeStepSys * SecInHour );
		this->TankSupplyVdot = TankSupplyVdot;
		this->TankSupplyVol = TankSupplyVdot * ( TimeStepSys * SecInHour );
		this->StarvedMakeUpVdot = StarvedVdot;
		this->StarvedMakeUpVol = StarvedVdot * ( TimeStepSys * SecInHour );

	}

	void
	EvaporativeFluidCooler::update()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    May 2009
		//       MODIFIED         na
		//       RE-ENGINEERED    na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for passing results to the outlet water node.

		// USE STATEMENTS:
		//unused0909  USE DataEnvironment, ONLY: EnvironmentName, CurMnDy
		// Using/Aliasing
		using General::TrimSigDigits;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static gio::Fmt LowTempFmt( "(' ',F6.2)" );
		Real64 const TempAllowance( 0.02 ); // Minimum difference b/w fluid cooler water outlet temp and
		// minimum condenser loop temp [C]

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		std::string CharErrOut;
		std::string CharLowOutletTemp;
		Real64 TempDifference;
		int LoopNum;
		int LoopSideNum;
		Real64 LoopMinTemp;

		// set node information

		Node( this->WaterOutletNodeNum ).Temp = OutletWaterTemp;

		LoopNum = this->location.loopNum;
		LoopSideNum = this->location.loopSideNum;
		if ( PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 0 || DataGlobals::WarmupFlag ) return;

		// Check flow rate through evaporative fluid cooler and compare to design flow rate,
		// show warning if greater than Design * Mulitplier
		if ( Node( this->WaterOutletNodeNum ).MassFlowRate > this->DesWaterMassFlowRate * this->EvapFluidCoolerMassFlowRateMultiplier ) {
			++this->HighMassFlowErrorCount;
			if ( this->HighMassFlowErrorCount < 2 ) {
				ShowWarningError( this->EvapFluidCoolerType + " \"" + this->Name + "\"" );
				ShowContinueError( " Condenser Loop Mass Flow Rate is much greater than the evaporative fluid coolers design mass flow rate." );
				ShowContinueError( " Condenser Loop Mass Flow Rate = " + TrimSigDigits( Node( this->WaterOutletNodeNum ).MassFlowRate, 6 ) );
				ShowContinueError( " Evaporative Fluid Cooler Design Mass Flow Rate   = " + TrimSigDigits( this->DesWaterMassFlowRate, 6 ) );
				ShowContinueErrorTimeStamp( "" );
			} else {
				ShowRecurringWarningErrorAtEnd( this->EvapFluidCoolerType + " \"" + this->Name + "\"  Condenser Loop Mass Flow Rate is much greater than the evaporative fluid coolers design mass flow rate error", this->HighMassFlowErrorIndex, Node( this->WaterOutletNodeNum ).MassFlowRate, Node( this->WaterOutletNodeNum ).MassFlowRate );
			}
		}

		// Check if OutletWaterTemp is below the minimum condenser loop temp and warn user
		LoopMinTemp = PlantLoop( LoopNum ).MinTemp;
		TempDifference = PlantLoop( LoopNum ).MinTemp - OutletWaterTemp;
		if ( TempDifference > TempAllowance && WaterMassFlowRate > 0.0 ) {
			++this->OutletWaterTempErrorCount;
			gio::write( CharLowOutletTemp, LowTempFmt ) << LoopMinTemp;
			gio::write( CharErrOut, LowTempFmt ) << OutletWaterTemp;
			strip( CharErrOut );
			if ( this->OutletWaterTempErrorCount < 2 ) {
				ShowWarningError( this->EvapFluidCoolerType + " \"" + this->Name + "\"" );
				ShowContinueError( "Evaporative fluid cooler water outlet temperature (" + CharErrOut + " C) is below the specified minimum condenser loop temp of " + stripped( CharLowOutletTemp ) + " C" );
				ShowContinueErrorTimeStamp( "" );
			} else {
				ShowRecurringWarningErrorAtEnd( this->EvapFluidCoolerType + " \"" + this->Name + "\" Evaporative fluid cooler water outlet temperature is below the specified minimum condenser loop temp error", this->OutletWaterTempErrorIndex, OutletWaterTemp, OutletWaterTemp );
			}
		}

		// Check if water mass flow rate is small (e.g. no flow) and warn user
		if ( WaterMassFlowRate > 0.0 && WaterMassFlowRate <= MassFlowTolerance ) {
			++this->SmallWaterMassFlowErrorCount;
			if ( this->SmallWaterMassFlowErrorCount < 2 ) {
				ShowWarningError( this->EvapFluidCoolerType + " \"" + this->Name + "\"" );
				ShowContinueError( "Evaporative fluid cooler water mass flow rate near zero." );
				ShowContinueErrorTimeStamp( "" );
				ShowContinueError( "Actual Mass flow = " + TrimSigDigits( WaterMassFlowRate, 2 ) );
			} else {
				ShowRecurringWarningErrorAtEnd( this->EvapFluidCoolerType + " \"" + this->Name + "\" Evaporative fluid cooler water mass flow rate near zero error continues...", this->SmallWaterMassFlowErrorIndex, WaterMassFlowRate, WaterMassFlowRate );
			}
		}

	}

	void
	EvaporativeFluidCooler::report(
		bool const RunFlag
	)
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    May 2009
		//       MODIFIED         na
		//       RE-ENGINEERED    na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine updates the report variables for the evaporative fluid cooler.

		if ( ! RunFlag ) {
			this->InletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
			this->OutletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
			this->Qactual = 0.0;
			this->FanPower = 0.0;
			this->FanEnergy = 0.0;
			this->AirFlowRatio = 0.0;
			this->WaterAmountUsed = 0.0;
			this->BypassFraction = 0.0; // added for fluid bypass
		} else {
			Real64 ReportingConstant = TimeStepSys * SecInHour;
			this->InletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
			this->FanEnergy = this->FanPower * ReportingConstant;
			this->WaterAmountUsed = this->WaterUsage * ReportingConstant;
		}

	}

} // EnergyPlus
