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
#include <FluidCoolers.hh>
#include <BranchNodeConnections.hh>
#include <CurveManager.hh>
#include <DataBranchAirLoopPlant.hh>
#include <DataEnvironment.hh>
#include <DataHVACGlobals.hh>
#include <DataIPShortCuts.hh>
#include <DataLoopNode.hh>
#include <DataPlant.hh>
#include <DataSizing.hh>
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

namespace EnergyPlus {

	// Module containing the routines dealing with the objects FluidCooler:SingleSpeed and
	// FluidCooler:TwoSpeed

	// MODULE INFORMATION:
	//       AUTHOR         Chandan Sharma
	//       DATE WRITTEN   August 2008
	//       MODIFIED       April 2010, Chandan Sharma, FSEC
	//       RE-ENGINEERED  na

	// PURPOSE OF THIS MODULE:
	// Model the performance of fluid coolers

	// METHODOLOGY EMPLOYED:

	// REFERENCES:
	// Based on cooling tower by Shirey, Raustad: Dec 2000; Shirey, Sept 2002

	using DataHVACGlobals::TimeStepSys;
	using DataGlobals::SecInHour;
	using DataPlant::PlantLoop;
	using DataLoopNode::Node;
	using Psychrometrics::PsyWFnTdbTwbPb;
	using FluidProperties::GetDensityGlycol;
	using FluidProperties::GetSpecificHeatGlycol;

	std::vector< std::unique_ptr< FluidCooler > > FluidCooler::instances;
	bool FluidCooler::GetInputFlag = true;
	std::string const FluidCooler::cFluidCooler_SingleSpeed( "FluidCooler:SingleSpeed" );
	std::string const FluidCooler::cFluidCooler_TwoSpeed( "FluidCooler:TwoSpeed" );

	void FluidCooler::clear_state()
	{
		FluidCooler::instances.clear();
		FluidCooler::GetInputFlag = true;
	}

	PlantComponent * FluidCooler::factory( int objectType, std::string objectName )
	{
		// Process the input data for pipes if it hasn't been done already
		if ( GetInputFlag ) {
			FluidCooler::getInput();
			GetInputFlag = false;
		}
		// Now look for this particular pipe in the list
		for ( auto const & fluidCooler : instances ) {
			if ( fluidCooler->PlantType_Num == objectType && fluidCooler->Name == objectName ) {
				return fluidCooler.get();
			}
		}
		// If we didn't find it, fatal
		ShowFatalError( "FluidCoolerSpecsFactory: Error getting inputs for fluid cooler named: " + objectName );
		// Shut up the compiler
		return nullptr;
	}

	void FluidCooler::onInitLoopEquip( PlantLocation const & EP_UNUSED( calledFromLocation ) )
	{
		this->init();
		this->size();
	}

	void FluidCooler::getDesignCapacities( PlantLocation const & EP_UNUSED( calledFromLocation ), Real64 & MaxLoad, Real64 & MinLoad, Real64 & OptLoad )
	{
		MinLoad = 0.0;
		MaxLoad = this->FluidCoolerNominalCapacity;
		OptLoad = this->FluidCoolerNominalCapacity;
	}

	void FluidCooler::simulate(
			PlantLocation const & EP_UNUSED( calledFromLocation ),
			bool const EP_UNUSED( FirstHVACIteration ),
			Real64 & CurLoad
	)
	{
		if ( ! ( this->FluidCoolerType_Num == FluidCoolerEnum::SingleSpeed ||
				this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed ) ) {
			ShowFatalError( "SimFluidCoolers: Invalid Fluid Cooler Type Requested = " + this->FluidCoolerType );
		}

		this->init();
		if ( this->FluidCoolerType_Num == FluidCoolerEnum::SingleSpeed ) {
			this->singleSpeedFluidCooler();
		}
		if ( this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed ) {
			this->twoSpeedFluidCooler();
		}
		update();
		bool runFlag = CurLoad > DataPlant::LoopDemandTol || CurLoad < ( -DataPlant::LoopDemandTol );
		report( runFlag );
	}

	void
	FluidCooler::getInput()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    August 2008
		//       MODIFIED         Chandan Sharma, FSEC, April 2010
		//       RE-ENGINEERED    na

		// PURPOSE OF THIS SUBROUTINE:
		// Obtains input data for fluid coolers and stores it in SimpleFluidCooler data structure.

		// METHODOLOGY EMPLOYED:
		// Uses "Get" routines to read in the data.

		// REFERENCES:
		// Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002;

		// Using/Aliasing
		using InputProcessor::GetNumObjectsFound;
		using InputProcessor::GetObjectItem;
		using InputProcessor::VerifyName;
		using InputProcessor::SameString;
		using namespace DataIPShortCuts; // Data for field names, blank numerics
		using NodeInputManager::GetOnlySingleNode;
		using BranchNodeConnections::TestCompSet;
		using DataSizing::AutoSize;
		using CurveManager::GetCurveIndex;
		using ScheduleManager::GetScheduleIndex;
		using OutAirNodeManager::CheckOutAirNodeNumber;
		using General::TrimSigDigits;
		using FluidProperties::CheckFluidPropertyName;
		using FluidProperties::FindGlycol;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:
		// na

		// SUBROUTINE PARAMETER DEFINITIONS:
		// na

		// INTERFACE BLOCK SPECIFICATIONS
		// na

		// DERIVED TYPE DEFINITIONS
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int FluidCoolerNum; // Fluid cooler number, reference counter for SimpleFluidCooler data array
		int NumSingleSpeedFluidCoolers; // Total number of single-speed Fluid Coolers
		int SingleSpeedFluidCoolerNumber; // Specific single-speed fluid cooler of interest
		int NumTwoSpeedFluidCoolers; // Number of two-speed Fluid Coolers
		int TwoSpeedFluidCoolerNumber; // Specific two-speed fluid cooler of interest
		int NumSimpleFluidCoolers;
		int NumAlphas; // Number of elements in the alpha array
		int NumNums; // Number of elements in the numeric array
		int IOStat; // IO Status when calling get input subroutine
		bool IsNotOK; // Flag to verify name
		bool IsBlank; // Flag for blank name
		bool ErrorsFound( false ); // Logical flag set .TRUE. if errors found while getting input data
		Array1D< Real64 > NumArray( 16 ); // Numeric input data array
		Array1D_string AlphArray( 5 ); // Character string input data array

		//! LKL - still more renaming stuff to go.

		// Get number of all Fluid Coolers specified in the input data file (idf)
		NumSingleSpeedFluidCoolers = GetNumObjectsFound( "FluidCooler:SingleSpeed" );
		NumTwoSpeedFluidCoolers = GetNumObjectsFound( "FluidCooler:TwoSpeed" );
		NumSimpleFluidCoolers = NumSingleSpeedFluidCoolers + NumTwoSpeedFluidCoolers;

		if ( NumSimpleFluidCoolers <= 0 ) ShowFatalError( "No fluid cooler objects found in input, however, a branch object has specified a fluid cooler. Search the input for fluid cooler to determine the cause for this error." );

		// See if load distribution manager has already gotten the input
		if ( instances.size() > 0 ) return;

		// Allocate data structures to hold fluid cooler input data, report data and fluid cooler inlet conditions
		for ( int i = 0; i < NumSimpleFluidCoolers; ++i ) {
			instances.emplace_back( new FluidCooler() );
		}

		// Load data structures with fluid cooler input data
		cCurrentModuleObject = cFluidCooler_SingleSpeed;
		for ( SingleSpeedFluidCoolerNumber = 1; SingleSpeedFluidCoolerNumber <= NumSingleSpeedFluidCoolers; ++SingleSpeedFluidCoolerNumber ) {
			FluidCoolerNum = SingleSpeedFluidCoolerNumber;
			GetObjectItem( cCurrentModuleObject, SingleSpeedFluidCoolerNumber, AlphArray, NumAlphas, NumArray, NumNums, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );
			IsNotOK = false;
			IsBlank = false;
			VerifyName( instances.begin(), instances.begin() + FluidCoolerNum, AlphArray( 1 ), IsNotOK, IsBlank, cCurrentModuleObject + " Name" );
			if ( IsNotOK ) {
				ErrorsFound = true;
				if ( IsBlank ) AlphArray( 1 ) = "xxxxx";
			}
			auto & singleSpeedFluidCooler( instances[ FluidCoolerNum - 1 ] );
			singleSpeedFluidCooler->Name = AlphArray( 1 );
			singleSpeedFluidCooler->FluidCoolerType = cCurrentModuleObject;
			singleSpeedFluidCooler->FluidCoolerType_Num = FluidCoolerEnum::SingleSpeed;
			singleSpeedFluidCooler->PlantType_Num = DataPlant::TypeOf_FluidCooler_SingleSpd;
			singleSpeedFluidCooler->FluidCoolerMassFlowRateMultiplier = 2.5;
			singleSpeedFluidCooler->WaterInletNodeNum = GetOnlySingleNode( AlphArray( 2 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Inlet, 1, DataLoopNode::ObjectIsNotParent );
			singleSpeedFluidCooler->WaterOutletNodeNum = GetOnlySingleNode( AlphArray( 3 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Outlet, 1, DataLoopNode::ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, AlphArray( 1 ), AlphArray( 2 ), AlphArray( 3 ), "Chilled Water Nodes" );
			singleSpeedFluidCooler->HighSpeedFluidCoolerUA = NumArray( 1 );
			singleSpeedFluidCooler->FluidCoolerNominalCapacity = NumArray( 2 );
			singleSpeedFluidCooler->DesignEnteringWaterTemp = NumArray( 3 );
			singleSpeedFluidCooler->DesignEnteringAirTemp = NumArray( 4 );
			singleSpeedFluidCooler->DesignEnteringAirWetBulbTemp = NumArray( 5 );
			singleSpeedFluidCooler->DesignWaterFlowRate = NumArray( 6 );
			if ( singleSpeedFluidCooler->DesignWaterFlowRate == AutoSize ) {
				singleSpeedFluidCooler->DesignWaterFlowRateWasAutoSized = true;
			}
			singleSpeedFluidCooler->HighSpeedAirFlowRate = NumArray( 7 );
			if ( singleSpeedFluidCooler->HighSpeedAirFlowRate == AutoSize ) {
				singleSpeedFluidCooler->HighSpeedAirFlowRateWasAutoSized = true;
			}
			singleSpeedFluidCooler->HighSpeedFanPower = NumArray( 8 );
			if ( singleSpeedFluidCooler->HighSpeedFanPower == AutoSize ) {
				singleSpeedFluidCooler->HighSpeedFanPowerWasAutoSized = true;
			}

			//   outdoor air inlet node
			if ( AlphArray( 5 ).empty() ) {
				singleSpeedFluidCooler->OutdoorAirInletNodeNum = 0;
			} else {
				singleSpeedFluidCooler->OutdoorAirInletNodeNum = GetOnlySingleNode( AlphArray( 5 ), ErrorsFound, cCurrentModuleObject, singleSpeedFluidCooler->Name, DataLoopNode::NodeType_Air, DataLoopNode::NodeConnectionType_OutsideAirReference, 1, DataLoopNode::ObjectIsNotParent );
				if ( ! CheckOutAirNodeNumber( singleSpeedFluidCooler->OutdoorAirInletNodeNum ) ) {
					ShowSevereError( cCurrentModuleObject + "= \"" + singleSpeedFluidCooler->Name + "\" " + cAlphaFieldNames( 5 ) + "= \"" + AlphArray( 5 ) + "\" not valid." );
					ShowContinueError( "...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node." );
					ErrorsFound = true;
				}
			}

			ErrorsFound = ErrorsFound || FluidCooler::verifySingleSpeedDesignInputs( singleSpeedFluidCooler, cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );

		} // End Single-Speed fluid cooler Loop

		cCurrentModuleObject = cFluidCooler_TwoSpeed;
		for ( TwoSpeedFluidCoolerNumber = 1; TwoSpeedFluidCoolerNumber <= NumTwoSpeedFluidCoolers; ++TwoSpeedFluidCoolerNumber ) {
			FluidCoolerNum = NumSingleSpeedFluidCoolers + TwoSpeedFluidCoolerNumber;
			GetObjectItem( cCurrentModuleObject, TwoSpeedFluidCoolerNumber, AlphArray, NumAlphas, NumArray, NumNums, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );
			IsNotOK = false;
			IsBlank = false;
			VerifyName( instances.begin(), instances.begin() + FluidCoolerNum, AlphArray( 1 ), IsNotOK, IsBlank, cCurrentModuleObject + " Name" );
			if ( IsNotOK ) {
				ErrorsFound = true;
				if ( IsBlank ) AlphArray( 1 ) = "xxxxx";
			}
			auto & twoSpeedFluidCooler( instances[ FluidCoolerNum - 1 ] );
			twoSpeedFluidCooler->Name = AlphArray( 1 );
			twoSpeedFluidCooler->FluidCoolerType = cCurrentModuleObject;
			twoSpeedFluidCooler->FluidCoolerType_Num = FluidCoolerEnum::TwoSpeed;
			twoSpeedFluidCooler->PlantType_Num = DataPlant::TypeOf_FluidCooler_TwoSpd;
			twoSpeedFluidCooler->FluidCoolerMassFlowRateMultiplier = 2.5;
			twoSpeedFluidCooler->WaterInletNodeNum = GetOnlySingleNode( AlphArray( 2 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Inlet, 1, DataLoopNode::ObjectIsNotParent );
			twoSpeedFluidCooler->WaterOutletNodeNum = GetOnlySingleNode( AlphArray( 3 ), ErrorsFound, cCurrentModuleObject, AlphArray( 1 ), DataLoopNode::NodeType_Water, DataLoopNode::NodeConnectionType_Outlet, 1, DataLoopNode::ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, AlphArray( 1 ), AlphArray( 2 ), AlphArray( 3 ), "Chilled Water Nodes" );

			twoSpeedFluidCooler->HighSpeedFluidCoolerUA = NumArray( 1 );
			if ( twoSpeedFluidCooler->HighSpeedFluidCoolerUA == AutoSize ) {
				twoSpeedFluidCooler->HighSpeedFluidCoolerUAWasAutoSized = true;
			}
			twoSpeedFluidCooler->LowSpeedFluidCoolerUA = NumArray( 2 );
			if ( twoSpeedFluidCooler->LowSpeedFluidCoolerUA == AutoSize ) {
				twoSpeedFluidCooler->LowSpeedFluidCoolerUAWasAutoSized = true;
			}
			twoSpeedFluidCooler->LowSpeedFluidCoolerUASizingFactor = NumArray( 3 );
			twoSpeedFluidCooler->FluidCoolerNominalCapacity = NumArray( 4 );
			twoSpeedFluidCooler->FluidCoolerLowSpeedNomCap = NumArray( 5 );
			if ( twoSpeedFluidCooler->FluidCoolerLowSpeedNomCap == AutoSize ) {
				twoSpeedFluidCooler->FluidCoolerLowSpeedNomCapWasAutoSized = true;
			}
			twoSpeedFluidCooler->FluidCoolerLowSpeedNomCapSizingFactor = NumArray( 6 );
			twoSpeedFluidCooler->DesignEnteringWaterTemp = NumArray( 7 );
			twoSpeedFluidCooler->DesignEnteringAirTemp = NumArray( 8 );
			twoSpeedFluidCooler->DesignEnteringAirWetBulbTemp = NumArray( 9 );
			twoSpeedFluidCooler->DesignWaterFlowRate = NumArray( 10 );
			if ( twoSpeedFluidCooler->DesignWaterFlowRate == AutoSize ) {
				twoSpeedFluidCooler->DesignWaterFlowRateWasAutoSized = true;
			}
			twoSpeedFluidCooler->HighSpeedAirFlowRate = NumArray( 11 );
			if ( twoSpeedFluidCooler->HighSpeedAirFlowRate == AutoSize ) {
				twoSpeedFluidCooler->HighSpeedAirFlowRateWasAutoSized = true;
			}
			twoSpeedFluidCooler->HighSpeedFanPower = NumArray( 12 );
			if ( twoSpeedFluidCooler->HighSpeedFanPower == AutoSize ) {
				twoSpeedFluidCooler->HighSpeedFanPowerWasAutoSized = true;
			}
			twoSpeedFluidCooler->LowSpeedAirFlowRate = NumArray( 13 );
			if ( twoSpeedFluidCooler->LowSpeedAirFlowRate == AutoSize ) {
				twoSpeedFluidCooler->LowSpeedAirFlowRateWasAutoSized = true;
			}
			twoSpeedFluidCooler->LowSpeedAirFlowRateSizingFactor = NumArray( 14 );
			twoSpeedFluidCooler->LowSpeedFanPower = NumArray( 15 );
			if ( twoSpeedFluidCooler->LowSpeedFanPower == AutoSize ) {
				twoSpeedFluidCooler->LowSpeedFanPowerWasAutoSized = true;
			}
			twoSpeedFluidCooler->LowSpeedFanPowerSizingFactor = NumArray( 16 );

			//   outdoor air inlet node
			if ( AlphArray( 5 ).empty() ) {
				twoSpeedFluidCooler->OutdoorAirInletNodeNum = 0;
			} else {
				twoSpeedFluidCooler->OutdoorAirInletNodeNum = GetOnlySingleNode( AlphArray( 5 ), ErrorsFound, cCurrentModuleObject, twoSpeedFluidCooler->Name, DataLoopNode::NodeType_Air, DataLoopNode::NodeConnectionType_OutsideAirReference, 1, DataLoopNode::ObjectIsNotParent );
				if ( !CheckOutAirNodeNumber( twoSpeedFluidCooler->OutdoorAirInletNodeNum ) ) {
					ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\" " + cAlphaFieldNames( 5 ) + "= \"" + AlphArray( 5 ) + "\" not valid." );
					ShowContinueError( "...does not appear in an OutdoorAir:NodeList or as an OutdoorAir:Node." );
					ErrorsFound = true;
				}
			}

			ErrorsFound = ErrorsFound || FluidCooler::verifyTwoSpeedDesignInputs( twoSpeedFluidCooler, cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
		}

		if ( ErrorsFound ) {
			ShowFatalError( "Errors found in getting fluid cooler input." );
		}

		// Set up output variables, CurrentModuleObject='FluidCooler:SingleSpeed'
		for ( FluidCoolerNum = 1; FluidCoolerNum <= NumSingleSpeedFluidCoolers; ++FluidCoolerNum ) {
			auto & singleSpeedFluidCooler( instances[ FluidCoolerNum - 1 ] );
			SetupOutputVariable( "Cooling Tower Inlet Temperature [C]", singleSpeedFluidCooler->InletWaterTemp, "System", "Average", singleSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Outlet Temperature [C]", singleSpeedFluidCooler->OutletWaterTemp, "System", "Average", singleSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Mass Flow Rate [kg/s]", singleSpeedFluidCooler->WaterMassFlowRate, "System", "Average", singleSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Heat Transfer Rate [W]", singleSpeedFluidCooler->Qactual, "System", "Average", singleSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Power [W]", singleSpeedFluidCooler->FanPower, "System", "Average", singleSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Energy [J]", singleSpeedFluidCooler->FanEnergy, "System", "Sum", singleSpeedFluidCooler->Name, _, "Electric", "HeatRejection", _, "Plant" );
		}

		// CurrentModuleObject='FluidCooler:TwoSpeed'
		for ( FluidCoolerNum = NumSingleSpeedFluidCoolers + 1; FluidCoolerNum <= NumSingleSpeedFluidCoolers + NumTwoSpeedFluidCoolers; ++FluidCoolerNum ) {
			auto & twoSpeedFluidCooler( instances[ FluidCoolerNum - 1 ] );
			SetupOutputVariable( "Cooling Tower Inlet Temperature [C]", twoSpeedFluidCooler->InletWaterTemp, "System", "Average", twoSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Outlet Temperature [C]", twoSpeedFluidCooler->OutletWaterTemp, "System", "Average", twoSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Mass Flow Rate [kg/s]", twoSpeedFluidCooler->WaterMassFlowRate, "System", "Average", twoSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Heat Transfer Rate [W]", twoSpeedFluidCooler->Qactual, "System", "Average", twoSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Power [W]", twoSpeedFluidCooler->FanPower, "System", "Average", twoSpeedFluidCooler->Name );
			SetupOutputVariable( "Cooling Tower Fan Electric Energy [J]", twoSpeedFluidCooler->FanEnergy, "System", "Sum", twoSpeedFluidCooler->Name, _, "Electric", "HeatRejection", _, "Plant" );
		}

	}

	bool
	FluidCooler::verifySingleSpeedDesignInputs(
		std::unique_ptr< FluidCooler > & singleSpeedFluidCooler,
		std::string const & cCurrentModuleObject,
		Array1D<std::string> const & AlphArray,
		Array1D<std::string> const & cNumericFieldNames,
		Array1D<std::string> const & cAlphaFieldNames
	)
	{
		// FUNCTION INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    August 2008
		//       MODIFIED         Chandan Sharma, FSEC, April 2010
		//       RE-ENGINEERED    Jason Glazer, GARD Analytics, February 2015, refactor into a separate function

		// PURPOSE OF THIS FUNCTION:
		// Separate the testing of inputs related to design so that it could be called from the unit tests

		// METHODOLOGY EMPLOYED:
		// na

		// REFERENCES:
		// Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002;

		// Using/Aliasing
		using DataSizing::AutoSize;
		using InputProcessor::SameString;

		// FUNCTION LOCAL VARIABLE DECLARATIONS:
		bool ErrorsFound = false;

		//   Design entering water temperature, design entering air temperature and design entering air
		//   wetbulb temperature must be specified for the both the performance input methods
		if ( singleSpeedFluidCooler->DesignEnteringWaterTemp <= 0.0 ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 3 ) + "\", entered value <= 0.0, but must be > 0 " );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->DesignEnteringAirTemp <= 0.0 ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 4 ) + "\", entered value <= 0.0, but must be > 0 " );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->DesignEnteringAirWetBulbTemp <= 0.0 ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 5 ) + "\", entered value <= 0.0, but must be > 0 " );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->DesignEnteringWaterTemp <= singleSpeedFluidCooler->DesignEnteringAirTemp ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\"," + cNumericFieldNames( 3 ) + " must be greater than " + cNumericFieldNames( 4 ) + '.' );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->DesignEnteringAirTemp <= singleSpeedFluidCooler->DesignEnteringAirWetBulbTemp ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\"," + cNumericFieldNames( 4 ) + " must be greater than " + cNumericFieldNames( 5 ) + '.' );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->HighSpeedAirFlowRate <= 0.0 && singleSpeedFluidCooler->HighSpeedAirFlowRate != AutoSize ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 7 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->DesignWaterFlowRate <= 0.0 && !singleSpeedFluidCooler->DesignWaterFlowRateWasAutoSized ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 6 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		if ( singleSpeedFluidCooler->HighSpeedFanPower <= 0.0 && singleSpeedFluidCooler->HighSpeedFanPower != AutoSize ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 8 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}

		//   Check various inputs for both the performance input methods
		if ( SameString( AlphArray( 4 ), "UFactorTimesAreaAndDesignWaterFlowRate" ) ) {
			singleSpeedFluidCooler->PerformanceInputMethod_Num = FluidCooler::PIM::UFactor;
			if ( singleSpeedFluidCooler->HighSpeedFluidCoolerUA <= 0.0 && singleSpeedFluidCooler->HighSpeedFluidCoolerUA != AutoSize ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 1 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
		}
		else if ( SameString( AlphArray( 4 ), "NominalCapacity" ) ) {
			singleSpeedFluidCooler->PerformanceInputMethod_Num = FluidCooler::PIM::NominalCapacity;
			if ( singleSpeedFluidCooler->FluidCoolerNominalCapacity <= 0.0 ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 2 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( singleSpeedFluidCooler->HighSpeedFluidCoolerUA != 0.0 ) {
				if ( singleSpeedFluidCooler->HighSpeedFluidCoolerUA > 0.0 ) {
					ShowSevereError( cCurrentModuleObject + "= \"" + singleSpeedFluidCooler->Name + "\". Nominal fluid cooler capacity and design fluid cooler UA have been specified." );
				}
				else {
					ShowSevereError( cCurrentModuleObject + "= \"" + singleSpeedFluidCooler->Name + "\". Nominal fluid cooler capacity has been specified and design fluid cooler UA is being autosized." );
				}
				ShowContinueError( "Design fluid cooler UA field must be left blank when nominal fluid cooler capacity performance input method is used." );
				ErrorsFound = true;
			}
		}
		else { // Fluid cooler performance input method is not specified as a valid "choice"
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\", invalid " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
			ShowContinueError( "... must be \"UFactorTimesAreaAndDesignWaterFlowRate\" or \"NominalCapacity\"." );
			ErrorsFound = true;
		}
		return ErrorsFound;
	}

	bool
	FluidCooler::verifyTwoSpeedDesignInputs(
		std::unique_ptr< FluidCooler > & twoSpeedFluidCooler,
		std::string const & cCurrentModuleObject,
		Array1D<std::string> const & AlphArray,
		Array1D<std::string> const & cNumericFieldNames,
		Array1D<std::string> const & cAlphaFieldNames
	)
	{
		// FUNCTION INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    August 2008
		//       MODIFIED         Chandan Sharma, FSEC, April 2010
		//       RE-ENGINEERED    Jason Glazer, GARD Analytics, February 2015, refactor into a separate function

		// PURPOSE OF THIS FUNCTION:
		// Separate the testing of inputs related to design so that it could be called from the unit tests

		// METHODOLOGY EMPLOYED:
		// na

		// REFERENCES:
		// Based on GetTowerInput subroutine from Don Shirey, Jan 2001 and Sept/Oct 2002;

		// Using/Aliasing
		using DataSizing::AutoSize;
		using InputProcessor::SameString;

		// FUNCTION LOCAL VARIABLE DECLARATIONS:
		bool ErrorsFound = false;

		//   Design entering water temperature, design entering air temperature and design entering air
		//   wetbulb temperature must be specified for the both the performance input methods
		if ( twoSpeedFluidCooler->DesignEnteringWaterTemp <= 0.0 ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 7 ) + "\", entered value <= 0.0, but must be > 0 " );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->DesignEnteringAirTemp <= 0.0 ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 8 ) + "\", entered value <= 0.0, but must be > 0 " );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->DesignEnteringAirWetBulbTemp <= 0.0 ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 9 ) + "\", entered value <= 0.0, but must be > 0 " );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->DesignEnteringWaterTemp <= twoSpeedFluidCooler->DesignEnteringAirTemp ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", " + cNumericFieldNames( 7 ) + " must be greater than " + cNumericFieldNames( 8 ) + '.' );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->DesignEnteringAirTemp <= twoSpeedFluidCooler->DesignEnteringAirWetBulbTemp ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", " + cNumericFieldNames( 8 ) + " must be greater than " + cNumericFieldNames( 9 ) + '.' );
			ErrorsFound = true;
		}

		//   Check various inputs for both the performance input methods
		if ( twoSpeedFluidCooler->DesignWaterFlowRate <= 0.0 && ! twoSpeedFluidCooler->DesignWaterFlowRateWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 10 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + "= \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->HighSpeedAirFlowRate <= 0.0 && ! twoSpeedFluidCooler->HighSpeedAirFlowRateWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 11 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + "= \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->LowSpeedAirFlowRate <= 0.0 && ! twoSpeedFluidCooler->LowSpeedAirFlowRateWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 13 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + "= \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		//   High speed air flow rate must be greater than low speed air flow rate.
		//   Can't tell yet if autosized, check later in InitFluidCooler.
		if ( twoSpeedFluidCooler->HighSpeedAirFlowRate <= twoSpeedFluidCooler->LowSpeedAirFlowRate && ! twoSpeedFluidCooler->HighSpeedAirFlowRateWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Fluid cooler air flow rate at low fan speed must be less than the air flow rate at high fan speed." );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->HighSpeedFanPower <= 0.0 && ! twoSpeedFluidCooler->HighSpeedFanPowerWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 12 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->LowSpeedFanPower <= 0.0 && ! twoSpeedFluidCooler->LowSpeedFanPowerWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 15 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
			ErrorsFound = true;
		}
		if ( twoSpeedFluidCooler->HighSpeedFanPower <= twoSpeedFluidCooler->LowSpeedFanPower && ! twoSpeedFluidCooler->HighSpeedFanPowerWasAutoSized ) {
			ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Fluid cooler low speed fan power must be less than high speed fan power." );
			ErrorsFound = true;
		}

		if ( SameString( AlphArray( 4 ), "UFactorTimesAreaAndDesignWaterFlowRate" ) ) {
			twoSpeedFluidCooler->PerformanceInputMethod_Num = PIM::UFactor;
			if ( twoSpeedFluidCooler->HighSpeedFluidCoolerUA <= 0.0 && ! twoSpeedFluidCooler->HighSpeedFluidCoolerUAWasAutoSized ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 1 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( twoSpeedFluidCooler->LowSpeedFluidCoolerUA <= 0.0 && ! twoSpeedFluidCooler->LowSpeedFluidCoolerUAWasAutoSized ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 2 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + " = \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( twoSpeedFluidCooler->HighSpeedFluidCoolerUA <= twoSpeedFluidCooler->LowSpeedFluidCoolerUA && ! twoSpeedFluidCooler->HighSpeedFluidCoolerUAWasAutoSized ) {
				ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Fluid cooler UA at low fan speed must be less than the fluid cooler UA at high fan speed." );
				ErrorsFound = true;
			}
		} else if ( SameString( AlphArray( 4 ), "NominalCapacity" ) ) {
			twoSpeedFluidCooler->PerformanceInputMethod_Num = PIM::NominalCapacity;
			if ( twoSpeedFluidCooler->FluidCoolerNominalCapacity <= 0.0 ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 4 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + "= \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( twoSpeedFluidCooler->FluidCoolerLowSpeedNomCap <= 0.0 && ! twoSpeedFluidCooler->FluidCoolerLowSpeedNomCapWasAutoSized ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + AlphArray( 1 ) + "\", invalid data for \"" + cNumericFieldNames( 5 ) + "\", entered value <= 0.0, but must be > 0 for " + cAlphaFieldNames( 4 ) + "= \"" + AlphArray( 4 ) + "\"." );
				ErrorsFound = true;
			}
			if ( twoSpeedFluidCooler->HighSpeedFluidCoolerUA != 0.0 ) {
				if ( twoSpeedFluidCooler->HighSpeedFluidCoolerUA > 0.0 ) {
					ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Nominal capacity input method and fluid cooler UA at high fan speed have been specified." );
				} else {
					ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Nominal capacity input method has been specified and fluid cooler UA at high fan speed is being autosized." );
				}
				ShowContinueError( "Fluid cooler UA at high fan speed must be left blank when nominal fluid cooler capacity performance input method is used." );
				ErrorsFound = true;
			}
			if ( twoSpeedFluidCooler->LowSpeedFluidCoolerUA != 0.0 ) {
				if ( twoSpeedFluidCooler->LowSpeedFluidCoolerUA > 0.0 ) {
					ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Nominal capacity input method and fluid cooler UA at low fan speed have been specified." );
				} else {
					ShowSevereError( cCurrentModuleObject + "= \"" + twoSpeedFluidCooler->Name + "\". Nominal capacity input method has been specified and fluid cooler UA at low fan speed is being autosized." );
				}
				ShowContinueError( "Fluid cooler UA at low fan speed must be left blank when nominal fluid cooler capacity performance input method is used." );
				ErrorsFound = true;
			}
			if ( twoSpeedFluidCooler->FluidCoolerLowSpeedNomCap >= twoSpeedFluidCooler->FluidCoolerNominalCapacity ) {
				ShowSevereError( cCurrentModuleObject + " = \"" + twoSpeedFluidCooler->Name + "\". Low-speed nominal capacity must be less than the high-speed nominal capacity." );
				ErrorsFound = true;
			}
		} else { // Fluid cooler performance input method is not specified as a valid "choice"
			ShowSevereError( cCurrentModuleObject + "= \"" + AlphArray( 1 ) + "\", invalid " + cAlphaFieldNames( 4 ) + "= \"" + AlphArray( 4 ) + "\"." );
			ShowContinueError( "... must be \"UFactorTimesAreaAndDesignWaterFlowRate\" or \"NominalCapacity\"." );
			ErrorsFound = true;
		}
		return ErrorsFound;
	}

	void
	FluidCooler::init()
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   August 2008
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for initializations of the fluid cooler components and for
		// final checking of fluid cooler inputs (post autosizing)

		// METHODOLOGY EMPLOYED:
		// Uses the status flags to trigger initializations.

		// REFERENCES:
		// Based on InitTower subroutine by Don Shirey Sept/Oct 2002, F Buhl Oct 2002

		static std::string const RoutineName( "InitFluidCooler" );

		if ( this->oneTimeFlag ) {
			bool ErrorsFound( false );
			// Locate the tower on the plant loops for later usage
			DataPlant::ScanPlantLoopsForObject( this->Name, this->PlantType_Num, this->location.loopNum,
												this->location.loopSideNum, this->location.branchNum, this->location.compNum,
												_, _, _, _, _, ErrorsFound );
			if ( ErrorsFound ) {
				ShowFatalError( "InitFluidCooler: Program terminated due to previous condition(s)." );
			}
		}

		// Begin environment initializations
		if ( this->envrnFlag && DataGlobals::BeginEnvrnFlag && DataPlant::PlantFirstSizesOkayToFinalize ) {

			Real64 rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp,
																		 PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
			this->DesWaterMassFlowRate = this->DesignWaterFlowRate * rho;
			PlantUtilities::InitComponentNodes( 0.0, this->DesWaterMassFlowRate, this->WaterInletNodeNum,
												this->WaterOutletNodeNum, this->location.loopNum, this->location.loopSideNum,
												this->location.branchNum, this->location.compNum );
			this->envrnFlag = false;
		}

		if ( ! DataGlobals::BeginEnvrnFlag ) {
			this->envrnFlag = true;
		}

		this->InletWaterTemp = 0.0; // CW temperature at fluid cooler inlet
		this->OutletWaterTemp = 0.0; // CW temperature at fluid cooler outlet
		this->WaterMassFlowRate = 0.0; // WaterMassFlowRate through fluid cooler
		this->Qactual = 0.0; // Fluid cooler heat transfer
		this->FanPower = 0.0; // Fluid cooler fan power used

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

		WaterMassFlowRate = PlantUtilities::RegulateCondenserCompFlowReqOp(
				this->location.loopNum, this->location.loopSideNum, this->location.branchNum, this->location.compNum,
				this->DesWaterMassFlowRate * this->FluidCoolerMassFlowRateMultiplier );

		PlantUtilities::SetComponentFlowRate( WaterMassFlowRate, this->WaterInletNodeNum, this->WaterOutletNodeNum,
																					this->location.loopNum, this->location.loopSideNum, this->location.branchNum,
																					this->location.compNum );
	}

	void
	FluidCooler::size()
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   August 2008
		//       MODIFIED       April 2010, Chandan Sharma, FSEC
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for sizing fluid cooler Components for which capacities and flow rates
		// have not been specified in the input. This subroutine also calculates fluid cooler UA if the user
		// has specified fluid cooler performance via the "Nominal Capacity" method.

		// METHODOLOGY EMPLOYED:
		// Obtains condenser flow rate from the plant sizing array. If fluid cooler performance is specified
		// via the "Nominal Capacity" method, the water flow rate is directly proportional to capacity.

		// REFERENCES:
		// Based on SizeTower by Don Shirey, Sept/Oct 2002; Richard Raustad, Feb 2005

		// Using/Aliasing
		using namespace DataSizing;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using DataPlant::PlantFirstSizesOkayToReport;
		using DataPlant::PlantFinalSizesOkayToReport;
		using namespace DataIPShortCuts; // Data for field names, blank numerics
		using General::SolveRegulaFalsi;
		using General::RoundSigDigits;
		using PlantUtilities::RegisterPlantCompDesignFlow;
		using ReportSizingManager::ReportSizingOutput;
		using namespace OutputReportPredefined;
		using InputProcessor::SameString;

		// SUBROUTINE PARAMETER DEFINITIONS:
		int const MaxIte( 500 ); // Maximum number of iterations
		Real64 const Acc( 0.0001 ); // Accuracy of result
		static std::string const CalledFrom( "SizeFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int PltSizCondNum( 0 ); // Plant Sizing index for condenser loop
		int SolFla; // Flag of solver
		Real64 DesFluidCoolerLoad( 0.0 ); // Design fluid cooler load [W]
		Real64 UA0; // Lower bound for UA [W/C]
		Real64 UA1; // Upper bound for UA [W/C]
		Real64 UA; // Calculated UA value
		Real64 OutWaterTempAtUA0; // Water outlet temperature at UA0
		Real64 OutWaterTempAtUA1; // Water outlet temperature at UA1
		Array1D< Real64 > Par( 4 ); // Parameter array need for RegulaFalsi routine
		Real64 Cp; // local specific heat for fluid
		Real64 rho; // local density for fluid
		Real64 tmpDesignWaterFlowRate; // local temporary for water volume flow rate
		Real64 tmpHighSpeedFanPower; // local temporary for high speed fan power
		Real64 tmpHighSpeedAirFlowRate; // local temporary for high speed air flow rate
		Real64 tmpHighSpeedEvapFluidCoolerUA; // local temporary for high speed cooler UA
		bool ErrorsFound;

		tmpDesignWaterFlowRate = this->DesignWaterFlowRate;
		tmpHighSpeedFanPower = this->HighSpeedFanPower;
		tmpHighSpeedAirFlowRate = this->HighSpeedAirFlowRate;
		tmpHighSpeedEvapFluidCoolerUA = this->HighSpeedFluidCoolerUA;
		// Find the appropriate Plant Sizing object
		PltSizCondNum = PlantLoop( this->location.loopNum ).PlantSizNum;

		if ( this->DesignWaterFlowRateWasAutoSized ) {
			if ( PltSizCondNum > 0 ) {
				if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow ) {
					tmpDesignWaterFlowRate = PlantSizData( PltSizCondNum ).DesVolFlowRate;
					if ( PlantFirstSizesOkayToFinalize ) this->DesignWaterFlowRate = tmpDesignWaterFlowRate;
				} else {
					tmpDesignWaterFlowRate = 0.0;
					if ( PlantFirstSizesOkayToFinalize ) this->DesignWaterFlowRate = tmpDesignWaterFlowRate;
				}
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Design Water Flow Rate [m3/s]", this->DesignWaterFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Design Water Flow Rate [m3/s]", this->DesignWaterFlowRate );
					}
				}
			} else {
				if ( PlantFirstSizesOkayToFinalize ) {
					ShowSevereError( "Autosizing error for fluid cooler object = " + this->Name );
					ShowFatalError( "Autosizing of fluid cooler condenser flow rate requires a loop Sizing:Plant object." );
				}
			}
			// This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
			// temperature is less than design inlet air dry bulb temperature
			if ( PlantSizData( PltSizCondNum ).ExitTemp <= this->DesignEnteringAirTemp && PlantFirstSizesOkayToFinalize ) {
				ShowSevereError( "Error when autosizing the UA value for fluid cooler = " + this->Name + '.' );
				ShowContinueError( "Design Loop Exit Temperature (" + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) + " C) must be greater than design entering air dry-bulb temperature (" + RoundSigDigits( this->DesignEnteringAirTemp, 2 ) + " C) when autosizing the fluid cooler UA." );
				ShowContinueError( "It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the Fluid Cooler design approach temperature (e.g., 4 C)." );
				ShowContinueError( "If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler." );
				ShowFatalError( "Review and revise design input values as appropriate." );
			}
		}

		RegisterPlantCompDesignFlow( this->WaterInletNodeNum, tmpDesignWaterFlowRate );

		if ( this->PerformanceInputMethod_Num == PIM::UFactor && this->HighSpeedFluidCoolerUAWasAutoSized ) {
			if ( PltSizCondNum > 0 ) {
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
				if ( PlantFirstSizesOkayToFinalize ) this->FluidCoolerNominalCapacity = DesFluidCoolerLoad;
			} else {
				if ( PlantFirstSizesOkayToFinalize ) this->FluidCoolerNominalCapacity = 0.0;
			}
		}

		if ( this->HighSpeedFanPowerWasAutoSized ) {
			// We assume the nominal fan power is 0.0105 times the design load
			if ( this->PerformanceInputMethod_Num == PIM::NominalCapacity ) {
				tmpHighSpeedFanPower = 0.0105 * this->FluidCoolerNominalCapacity;
				if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
			} else {
				if ( DesFluidCoolerLoad > 0.0 ) {
					tmpHighSpeedFanPower = 0.0105 * DesFluidCoolerLoad;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
				} else if ( PltSizCondNum > 0 ) {
					if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow ) {
						// This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
						// temperature is less than design inlet air dry bulb temperature
						if ( PlantSizData( PltSizCondNum ).ExitTemp <= this->DesignEnteringAirTemp && PlantFirstSizesOkayToFinalize ) {
							ShowSevereError( "Error when autosizing the UA value for fluid cooler = " + this->Name + '.' );
							ShowContinueError( "Design Loop Exit Temperature (" + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) + " C) must be greater than design entering air dry-bulb temperature (" + RoundSigDigits( this->DesignEnteringAirTemp, 2 ) + " C) when autosizing the fluid cooler UA." );
							ShowContinueError( "It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the Fluid Cooler design approach temperature (e.g., 4 C)." );
							ShowContinueError( "If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler." );
							ShowFatalError( "Review and revise design input values as appropriate." );
						}
						rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
						Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
						DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
						tmpHighSpeedFanPower = 0.0105 * DesFluidCoolerLoad;
						if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
					} else {
						tmpHighSpeedFanPower = 0.0;
						if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFanPower = tmpHighSpeedFanPower;
					}
				} else {
					if ( PlantFirstSizesOkayToFinalize ) {
						ShowSevereError( "Autosizing of fluid cooler fan power requires a loop Sizing:Plant object." );
						ShowFatalError( " Occurs in fluid cooler object = " + this->Name );
					}
				}
			}
			if ( this->FluidCoolerType_Num == FluidCoolerEnum::SingleSpeed ) {
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Fan Power at Design Air Flow Rate [W]", this->HighSpeedFanPower );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Fan Power at Design Air Flow Rate [W]", this->HighSpeedFanPower );
					}
				}
			} else if ( this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed ) {
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Fan Power at High Fan Speed [W]", this->HighSpeedFanPower );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Fan Power at High Fan Speed [W]", this->HighSpeedFanPower );
					}
				}
			}
		}

		if ( this->HighSpeedAirFlowRateWasAutoSized ) {
			if ( this->PerformanceInputMethod_Num == PIM::NominalCapacity ) {
				tmpHighSpeedAirFlowRate = this->FluidCoolerNominalCapacity / ( this->DesignEnteringWaterTemp - this->DesignEnteringAirTemp ) * 4.0;
				if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
			} else {
				if ( DesFluidCoolerLoad > 0.0 ) {
					tmpHighSpeedAirFlowRate = DesFluidCoolerLoad / ( this->DesignEnteringWaterTemp - this->DesignEnteringAirTemp ) * 4.0;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
				} else if ( PltSizCondNum > 0 ) {
					if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow ) {
						// This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
						// temperature is less than design inlet air dry bulb temperature
						if ( PlantSizData( PltSizCondNum ).ExitTemp <= this->DesignEnteringAirTemp && PlantFirstSizesOkayToFinalize ) {
							ShowSevereError( "Error when autosizing the UA value for fluid cooler = " + this->Name + '.' );
							ShowContinueError( "Design Loop Exit Temperature (" + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) + " C) must be greater than design entering air dry-bulb temperature (" + RoundSigDigits( this->DesignEnteringAirTemp, 2 ) + " C) when autosizing the fluid cooler UA." );
							ShowContinueError( "It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the Fluid Cooler design approach temperature (e.g., 4 C)." );
							ShowContinueError( "If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler." );
							ShowFatalError( "Review and revise design input values as appropriate." );
						}
						rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
						Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
						DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
						tmpHighSpeedAirFlowRate = DesFluidCoolerLoad / ( this->DesignEnteringWaterTemp - this->DesignEnteringAirTemp ) * 4.0;
						if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
					} else {
						tmpHighSpeedAirFlowRate = 0.0;
						if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedAirFlowRate = tmpHighSpeedAirFlowRate;
					}
				} else {
					if ( PlantFirstSizesOkayToFinalize ) {
						ShowSevereError( "Autosizing of fluid cooler air flow rate requires a loop Sizing:Plant object" );
						ShowFatalError( " Occurs in fluid cooler object = " + this->Name );
					}
				}
			}
			if ( this->FluidCoolerType_Num == FluidCoolerEnum::SingleSpeed ) {
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Design Air Flow Rate [m3/s]", this->HighSpeedAirFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Design Air Flow Rate [m3/s]", this->HighSpeedAirFlowRate );
					}
				}
			} else if ( this->FluidCoolerType == "FluidCooler:TwoSpeed" ) {
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Air Flow Rate at High Fan Speed [m3/s]", this->HighSpeedAirFlowRate );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Air Flow Rate at High Fan Speed [m3/s]", this->HighSpeedAirFlowRate );
					}
				}
			}
		}

		if ( this->HighSpeedFluidCoolerUAWasAutoSized ) {
			if ( PltSizCondNum > 0 ) {
				if ( PlantSizData( PltSizCondNum ).DesVolFlowRate >= DataHVACGlobals::SmallWaterVolFlow ) {
					// This conditional statement is to trap when the user specified Condenser/Fluid Cooler water design setpoint
					// temperature is less than design inlet air dry bulb temperature
					if ( PlantSizData( PltSizCondNum ).ExitTemp <= this->DesignEnteringAirTemp && PlantFirstSizesOkayToFinalize ) {
						ShowSevereError( "Error when autosizing the UA value for fluid cooler = " + this->Name + '.' );
						ShowContinueError( "Design Loop Exit Temperature (" + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) + " C) must be greater than design entering air dry-bulb temperature (" + RoundSigDigits( this->DesignEnteringAirTemp, 2 ) + " C) when autosizing the fluid cooler UA." );
						ShowContinueError( "It is recommended that the Design Loop Exit Temperature = design inlet air dry-bulb temp plus the Fluid Cooler design approach temperature (e.g., 4 C)." );
						ShowContinueError( "If using HVACTemplate:Plant:ChilledWaterLoop, then check that input field Condenser Water Design Setpoint must be > design inlet air dry-bulb temp if autosizing the Fluid Cooler." );
						ShowFatalError( "Review and revise design input values as appropriate." );
					}
					rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
					Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, PlantSizData( PltSizCondNum ).ExitTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
					DesFluidCoolerLoad = rho * Cp * tmpDesignWaterFlowRate * PlantSizData( PltSizCondNum ).DeltaT;
					Par( 1 ) = DesFluidCoolerLoad;
					Par( 2 ) = rho * tmpDesignWaterFlowRate; // design water mass flow rate
					Par( 3 ) = tmpHighSpeedAirFlowRate; // design air volume flow rate
					Par( 4 ) = Cp;
					UA0 = 0.0001 * DesFluidCoolerLoad; // Assume deltaT = 10000K (limit)
					UA1 = DesFluidCoolerLoad; // Assume deltaT = 1K
					this->WaterTemp = PlantSizData( PltSizCondNum ).ExitTemp + PlantSizData( PltSizCondNum ).DeltaT;
					this->AirTemp = this->DesignEnteringAirTemp;
					this->AirWetBulb = this->DesignEnteringAirWetBulbTemp;
					this->AirPress = DataEnvironment::StdBaroPress;
					this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress, CalledFrom );
					SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
					if ( SolFla == -1 ) {
						ShowWarningError( "Iteration limit exceeded in calculating fluid cooler UA." );
						ShowContinueError( "Autosizing of fluid cooler UA failed for fluid cooler = " + this->Name );
						ShowContinueError( "The final UA value =" + RoundSigDigits( UA, 2 ) + " W/K, and the simulation continues..." );
					} else if ( SolFla == -2 ) {
						this->simSimpleFluidCooler( Par( 3 ), Par( 4 ), UA0, OutWaterTempAtUA0 );
						this->simSimpleFluidCooler( Par( 3 ), Par( 4 ), UA1, OutWaterTempAtUA1 );
						ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
						ShowContinueError( "reasonable UA value. Review and revise design input values as appropriate. Specifying hard" );
						ShowContinueError( "sizes for some \"autosizable\" fields while autosizing other \"autosizable\" fields may be " );
						ShowContinueError( "contributing to this problem." );
						ShowContinueError( "This model iterates on UA to find the heat transfer required to provide the design outlet " );
						ShowContinueError( "water temperature. Initially, the outlet water temperatures at high and low UA values are " );
						ShowContinueError( "calculated. The Design Exit Water Temperature should be between the outlet water " );
						ShowContinueError( "temperatures calculated at high and low UA values. If the Design Exit Water Temperature is " );
						ShowContinueError( "out of this range, the solution will not converge and UA will not be calculated. " );
						ShowContinueError( "The possible solutions could be to manually input adjusted water and/or air flow rates based " );
						ShowContinueError( "on the autosized values shown below or to adjust design fluid cooler air inlet dry-bulb temperature." );
						ShowContinueError( "Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp)." );
						ShowContinueError( "Inputs to the fluid cooler object:" );
						ShowContinueError( "Design Fluid Cooler Load [W]                       = " + RoundSigDigits( Par( 1 ), 2 ) );
						ShowContinueError( "Design Fluid Cooler Water Volume Flow Rate [m3/s]  = " + RoundSigDigits( this->DesignWaterFlowRate, 6 ) );
						ShowContinueError( "Design Fluid Cooler Air Volume Flow Rate [m3/s]    = " + RoundSigDigits( Par( 4 ), 2 ) );
						ShowContinueError( "Design Fluid Cooler Air Inlet Dry-bulb Temp [C]    = " + RoundSigDigits( this->AirTemp, 2 ) );
						ShowContinueError( "Inputs to the plant sizing object:" );
						ShowContinueError( "Design Exit Water Temp [C]                         = " + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) );
						ShowContinueError( "Loop Design Temperature Difference [C]             = " + RoundSigDigits( PlantSizData( PltSizCondNum ).DeltaT, 2 ) );
						ShowContinueError( "Design Fluid Cooler Water Inlet Temp [C]           = " + RoundSigDigits( this->WaterTemp, 2 ) );
						ShowContinueError( "Calculated water outlet temp at low UA [C] (UA = " + RoundSigDigits( UA0, 2 ) + " W/K) = " + RoundSigDigits( OutWaterTempAtUA0, 2 ) );
						ShowContinueError( "Calculated water outlet temp at high UA [C](UA = " + RoundSigDigits( UA1, 2 ) + " W/K) = " + RoundSigDigits( OutWaterTempAtUA1, 2 ) );
						ShowFatalError( "Autosizing of Fluid Cooler UA failed for fluid cooler = " + this->Name );
					}
					tmpHighSpeedEvapFluidCoolerUA = UA;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFluidCoolerUA = tmpHighSpeedEvapFluidCoolerUA;
					this->FluidCoolerNominalCapacity = DesFluidCoolerLoad;
				} else {
					tmpHighSpeedEvapFluidCoolerUA = 0.0;
					if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFluidCoolerUA = tmpHighSpeedEvapFluidCoolerUA;
				}
				if ( this->FluidCoolerType_Num == FluidCoolerEnum::SingleSpeed ) {
					if ( PlantFirstSizesOkayToFinalize ) {
						if ( PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( this->FluidCoolerType, this->Name,
																	"U-factor Times Area Value at Design Air Flow Rate [W/K]", this->HighSpeedFluidCoolerUA );
						}
						if ( PlantFirstSizesOkayToReport ) {
							ReportSizingOutput( this->FluidCoolerType, this->Name,
																	"Initial U-factor Times Area Value at Design Air Flow Rate [W/K]", this->HighSpeedFluidCoolerUA );
						}
					}
				} else if ( this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed ) {
					if ( PlantFirstSizesOkayToFinalize ) {
						if ( PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( this->FluidCoolerType, this->Name,
																	"U-factor Times Area Value at High Fan Speed [W/K]", this->HighSpeedFluidCoolerUA );
						}
						if ( PlantFirstSizesOkayToReport ) {
							ReportSizingOutput( this->FluidCoolerType, this->Name,
																	"Initial U-factor Times Area Value at High Fan Speed [W/K]", this->HighSpeedFluidCoolerUA );
						}
					}
				}
			} else {
				if ( PlantFirstSizesOkayToFinalize ) {
					ShowSevereError( "Autosizing error for fluid cooler object = " + this->Name );
					ShowFatalError( "Autosizing of fluid cooler UA requires a loop Sizing:Plant object." );
				}
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::NominalCapacity ) {
			if ( this->DesignWaterFlowRate >= DataHVACGlobals::SmallWaterVolFlow ) {
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, this->DesignEnteringWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesFluidCoolerLoad = this->FluidCoolerNominalCapacity;
				Par( 1 ) = DesFluidCoolerLoad;
				Par( 2 ) = rho * tmpDesignWaterFlowRate; // design water mass flow rate
				Par( 3 ) = tmpHighSpeedAirFlowRate; // design air volume flow rate
				Par( 4 ) = Cp;
				UA0 = 0.0001 * DesFluidCoolerLoad; // Assume deltaT = 10000K (limit)
				UA1 = DesFluidCoolerLoad; // Assume deltaT = 1K
				this->WaterTemp = this->DesignEnteringWaterTemp; // design inlet water temperature
				this->AirTemp = this->DesignEnteringAirTemp; // design inlet air dry-bulb temp
				this->AirWetBulb = this->DesignEnteringAirWetBulbTemp; // design inlet air wet-bulb temp
				this->AirPress = DataEnvironment::StdBaroPress;
				this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress );
				SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
				if ( SolFla == -1 ) {
					ShowWarningError( "Iteration limit exceeded in calculating fluid cooler UA." );
					ShowContinueError( "Autosizing of fluid cooler UA failed for fluid cooler = " + this->Name );
					ShowContinueError( "The final UA value =" + RoundSigDigits( UA, 2 ) + " W/K, and the simulation continues..." );
				} else if ( SolFla == -2 ) {
					this->simSimpleFluidCooler( Par( 3 ), Par( 4 ), UA0, OutWaterTempAtUA0 );
					this->simSimpleFluidCooler( Par( 3 ), Par( 4 ), UA1, OutWaterTempAtUA1 );
					ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
					ShowContinueError( "reasonable UA value. Review and revise design input values as appropriate. Specifying hard" );
					ShowContinueError( "sizes for some \"autosizable\" fields while autosizing other \"autosizable\" fields may be " );
					ShowContinueError( "contributing to this problem." );
					ShowContinueError( "This model iterates on UA to find the heat transfer required to provide the design outlet " );
					ShowContinueError( "water temperature. Initially, the outlet water temperatures at high and low UA values are " );
					ShowContinueError( "calculated. The Design Exit Water Temperature should be between the outlet water " );
					ShowContinueError( "temperatures calculated at high and low UA values. If the Design Exit Water Temperature is " );
					ShowContinueError( "out of this range, the solution will not converge and UA will not be calculated. " );
					ShowContinueError( "The possible solutions could be to manually input adjusted water and/or air flow rates based " );
					ShowContinueError( "on the autosized values shown below or to adjust design fluid cooler air inlet dry-bulb temperature." );
					ShowContinueError( "Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp)." );
					ShowContinueError( "Inputs to the fluid cooler object:" );
					ShowContinueError( "Design Fluid Cooler Load [W]                       = " + RoundSigDigits( Par( 1 ), 2 ) );
					ShowContinueError( "Design Fluid Cooler Water Volume Flow Rate [m3/s]  = " + RoundSigDigits( this->DesignWaterFlowRate, 6 ) );
					ShowContinueError( "Design Fluid Cooler Air Volume Flow Rate [m3/s]    = " + RoundSigDigits( Par( 4 ), 2 ) );
					ShowContinueError( "Design Fluid Cooler Air Inlet Dry-bulb Temp [C]    = " + RoundSigDigits( this->AirTemp, 2 ) );
					ShowContinueError( "Inputs to the plant sizing object:" );
					ShowContinueError( "Design Exit Water Temp [C]                         = " + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) );
					ShowContinueError( "Loop Design Temperature Difference [C]             = " + RoundSigDigits( PlantSizData( PltSizCondNum ).DeltaT, 2 ) );
					ShowContinueError( "Design Fluid Cooler Water Inlet Temp [C]           = " + RoundSigDigits( this->WaterTemp, 2 ) );
					ShowContinueError( "Calculated water outlet temp at low UA [C] (UA = " + RoundSigDigits( UA0, 2 ) + " W/K) = " + RoundSigDigits( OutWaterTempAtUA0, 2 ) );
					ShowContinueError( "Calculated water outlet temp at high UA [C] (UA = " + RoundSigDigits( UA1, 2 ) + " W/K) = " + RoundSigDigits( OutWaterTempAtUA1, 2 ) );
					ShowFatalError( "Autosizing of Fluid Cooler UA failed for fluid cooler = " + this->Name );
				}
				if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFluidCoolerUA = UA;
			} else {
				if ( PlantFirstSizesOkayToFinalize ) this->HighSpeedFluidCoolerUA = 0.0;
			}
			if ( this->FluidCoolerType_Num == FluidCoolerEnum::SingleSpeed ) {
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Fluid cooler UA value at design air flow rate based on nominal capacity input [W/K]", this->HighSpeedFluidCoolerUA );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Fluid cooler UA value at design air flow rate based on nominal capacity input [W/K]", this->HighSpeedFluidCoolerUA );
					}
				}
			} else if ( this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed ) {
				if ( PlantFirstSizesOkayToFinalize ) {
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Fluid cooler UA value at high fan speed based on nominal capacity input [W/K]", this->HighSpeedFluidCoolerUA );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( this->FluidCoolerType, this->Name,
																"Initial Fluid cooler UA value at high fan speed based on nominal capacity input [W/K]", this->HighSpeedFluidCoolerUA );
					}
				}
			}
		}

		if ( this->LowSpeedAirFlowRateWasAutoSized && PlantFirstSizesOkayToFinalize ) {
			this->LowSpeedAirFlowRate = this->LowSpeedAirFlowRateSizingFactor * this->HighSpeedAirFlowRate;
			if ( PlantFinalSizesOkayToReport ) {
				ReportSizingOutput( this->FluidCoolerType, this->Name,
														"Air Flow Rate at Low Fan Speed [m3/s]", this->LowSpeedAirFlowRate );
			}
			if ( PlantFirstSizesOkayToReport ) {
				ReportSizingOutput( this->FluidCoolerType, this->Name,
														"Initial Air Flow Rate at Low Fan Speed [m3/s]", this->LowSpeedAirFlowRate );
			}
		}

		if ( this->LowSpeedFanPowerWasAutoSized && PlantFirstSizesOkayToFinalize ) {
			this->LowSpeedFanPower = this->LowSpeedFanPowerSizingFactor * this->HighSpeedFanPower;
			if ( PlantFinalSizesOkayToReport ) {
				ReportSizingOutput( this->FluidCoolerType, this->Name,
														"Fan Power at Low Fan Speed [W]", this->LowSpeedFanPower );
			}
			if ( PlantFirstSizesOkayToReport ) {
				ReportSizingOutput( this->FluidCoolerType, this->Name,
														"Initial Fan Power at Low Fan Speed [W]", this->LowSpeedFanPower );
			}
		}

		if ( this->LowSpeedFluidCoolerUAWasAutoSized && PlantFirstSizesOkayToFinalize ) {
			this->LowSpeedFluidCoolerUA = this->LowSpeedFluidCoolerUASizingFactor * this->HighSpeedFluidCoolerUA;
			if ( PlantFinalSizesOkayToReport ) {
				ReportSizingOutput( this->FluidCoolerType, this->Name,
														"U-factor Times Area Value at Low Fan Speed [W/K]", this->LowSpeedFluidCoolerUA );
			}
			if ( PlantFirstSizesOkayToReport ) {
				ReportSizingOutput( this->FluidCoolerType, this->Name,
														"Initial U-factor Times Area Value at Low Fan Speed [W/K]", this->LowSpeedFluidCoolerUA );
			}
		}

		if ( this->PerformanceInputMethod_Num == PIM::NominalCapacity && this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed ) {
			if ( this->FluidCoolerLowSpeedNomCapWasAutoSized && PlantFirstSizesOkayToFinalize ) {
				this->FluidCoolerLowSpeedNomCap = this->FluidCoolerLowSpeedNomCapSizingFactor * this->FluidCoolerNominalCapacity;
				if ( PlantFinalSizesOkayToReport ) {
					ReportSizingOutput( this->FluidCoolerType, this->Name,
															"Low Fan Speed Nominal Capacity [W]", this->FluidCoolerLowSpeedNomCap );
				}
				if ( PlantFirstSizesOkayToReport ) {
					ReportSizingOutput( this->FluidCoolerType, this->Name,
															"Initial Low Fan Speed Nominal Capacity [W]", this->FluidCoolerLowSpeedNomCap );
				}
			}

			if ( this->DesignWaterFlowRate >= DataHVACGlobals::SmallWaterVolFlow && this->FluidCoolerLowSpeedNomCap > 0.0 ) {
				rho = GetDensityGlycol( PlantLoop( this->location.loopNum ).FluidName, DataGlobals::InitConvTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, this->DesignEnteringWaterTemp, PlantLoop( this->location.loopNum ).FluidIndex, CalledFrom );
				DesFluidCoolerLoad = this->FluidCoolerLowSpeedNomCap;
				Par( 1 ) = DesFluidCoolerLoad;
				Par( 2 ) = rho * tmpDesignWaterFlowRate; // design water mass flow rate
				Par( 3 ) = this->LowSpeedAirFlowRate; // Air volume flow rate at low fan speed
				Par( 4 ) = Cp;
				UA0 = 0.0001 * DesFluidCoolerLoad; // Assume deltaT = 10000K (limit)
				UA1 = DesFluidCoolerLoad; // Assume deltaT = 1K
				this->WaterTemp = this->DesignEnteringWaterTemp; // design inlet water temperature
				this->AirTemp = this->DesignEnteringAirTemp; // design inlet air dry-bulb temp
				this->AirWetBulb = this->DesignEnteringAirWetBulbTemp; // design inlet air wet-bulb temp
				this->AirPress = DataEnvironment::StdBaroPress;
				this->AirHumRat = PsyWFnTdbTwbPb( this->AirTemp, this->AirWetBulb, this->AirPress, CalledFrom );
				SolveRegulaFalsi( Acc, MaxIte, SolFla, UA, [ this ]( Real64 const X, Array1< Real64 > const & Par ) -> Real64 { return simpleFluidCoolerUAResidual( X, Par ); }, UA0, UA1, Par );
				if ( SolFla == -1 ) {
					ShowWarningError( "Iteration limit exceeded in calculating fluid cooler UA." );
					ShowContinueError( "Autosizing of fluid cooler UA failed for fluid cooler = " + this->Name );
					ShowContinueError( "The final UA value at low fan speed =" + RoundSigDigits( UA, 2 ) + " W/C, and the simulation continues..." );
				} else if ( SolFla == -2 ) {
					this->simSimpleFluidCooler( Par( 3 ), Par( 4 ), UA0, OutWaterTempAtUA0 );
					this->simSimpleFluidCooler( Par( 3 ), Par( 4 ), UA1, OutWaterTempAtUA1 );
					ShowSevereError( CalledFrom + ": The combination of design input values did not allow the calculation of a " );
					ShowContinueError( "reasonable low-speed UA value. Review and revise design input values as appropriate. " );
					ShowContinueError( "Specifying hard sizes for some \"autosizable\" fields while autosizing other \"autosizable\" " );
					ShowContinueError( "fields may be contributing to this problem." );
					ShowContinueError( "This model iterates on UA to find the heat transfer required to provide the design outlet " );
					ShowContinueError( "water temperature. Initially, the outlet water temperatures at high and low UA values are " );
					ShowContinueError( "calculated. The Design Exit Water Temperature should be between the outlet water " );
					ShowContinueError( "temperatures calculated at high and low UA values. If the Design Exit Water Temperature is " );
					ShowContinueError( "out of this range, the solution will not converge and UA will not be calculated. " );
					ShowContinueError( "The possible solutions could be to manually input adjusted water and/or air flow rates based " );
					ShowContinueError( "on the autosized values shown below or to adjust design fluid cooler air inlet dry-bulb temperature." );
					ShowContinueError( "Plant:Sizing object inputs also influence these results (e.g. DeltaT and ExitTemp)." );
					ShowContinueError( "Inputs to the fluid cooler object:" );
					ShowContinueError( "Design Fluid Cooler Load [W]                         = " + RoundSigDigits( Par( 1 ), 2 ) );
					ShowContinueError( "Design Fluid Cooler Water Volume Flow Rate [m3/s]    = " + RoundSigDigits( this->DesignWaterFlowRate, 6 ) );
					ShowContinueError( "Design Fluid Cooler Air Volume Flow Rate [m3/s]      = " + RoundSigDigits( Par( 4 ), 2 ) );
					ShowContinueError( "Design Fluid Cooler Air Inlet Dry-bulb Temp [C]      = " + RoundSigDigits( this->AirTemp, 2 ) );
					ShowContinueError( "Inputs to the plant sizing object:" );
					ShowContinueError( "Design Exit Water Temp [C]                           = " + RoundSigDigits( PlantSizData( PltSizCondNum ).ExitTemp, 2 ) );
					ShowContinueError( "Loop Design Temperature Difference [C]               = " + RoundSigDigits( PlantSizData( PltSizCondNum ).DeltaT, 2 ) );
					ShowContinueError( "Design Fluid Cooler Water Inlet Temp [C]             = " + RoundSigDigits( this->WaterTemp, 2 ) );
					ShowContinueError( "Calculated water outlet temp at low UA [C](UA = " + RoundSigDigits( UA0, 2 ) + " W/C) = " + RoundSigDigits( OutWaterTempAtUA0, 2 ) );
					ShowContinueError( "Calculated water outlet temp at high UA [C](UA = " + RoundSigDigits( UA1, 2 ) + " W/C) = " + RoundSigDigits( OutWaterTempAtUA1, 2 ) );
					ShowFatalError( "Autosizing of Fluid Cooler UA failed for fluid cooler = " + this->Name );
				}
				if ( PlantFirstSizesOkayToFinalize ) this->LowSpeedFluidCoolerUA = UA;
			} else {
				if ( PlantFirstSizesOkayToFinalize ) this->LowSpeedFluidCoolerUA = 0.0;
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( PlantFinalSizesOkayToReport ) {
					ReportSizingOutput( this->FluidCoolerType, this->Name,
															"U-factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedFluidCoolerUA );
				}
				if ( PlantFirstSizesOkayToReport ) {
					ReportSizingOutput( this->FluidCoolerType, this->Name,
															"Initial U-factor Times Area Value at Low Fan Speed [W/C]", this->LowSpeedFluidCoolerUA );
				}
			}
		}

		ErrorsFound = false;

		if ( PlantFinalSizesOkayToReport ) {
			//create predefined report
			PreDefTableEntry( pdchMechType, this->Name, this->FluidCoolerType );
			PreDefTableEntry( pdchMechNomCap, this->Name, this->FluidCoolerNominalCapacity );
		}

		if ( this->FluidCoolerType_Num == FluidCoolerEnum::TwoSpeed && PlantFirstSizesOkayToFinalize ) {
			if ( this->DesignWaterFlowRate > 0.0 ) {
				if ( this->HighSpeedAirFlowRate <= this->LowSpeedAirFlowRate ) {
					ShowSevereError( "FluidCooler:TwoSpeed  \"" + this->Name + "\". Low speed air flow rate must be less than high speed air flow rate." );
					ErrorsFound = true;
				}
				if ( this->HighSpeedFluidCoolerUA <= this->LowSpeedFluidCoolerUA ) {
					ShowSevereError( "FluidCooler:TwoSpeed  \"" + this->Name + "\". Fluid cooler UA at low fan speed must be less than the fluid cooler UA at high fan speed." );
					ErrorsFound = true;
				}
			}
		}

		if ( ErrorsFound ) {
			ShowFatalError( "SizeFluidCooler: Program terminated due to previous condition(s)." );
		}
	}

	// End Initialization Section for the CondenserLoopFluidCoolers Module
	//******************************************************************************

	// Beginning of the CondenserLoopFluidCoolers Module Simulation Subroutines
	// *****************************************************************************

	void
	FluidCooler::singleSpeedFluidCooler()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   August 2008
		//       MODIFIED       Dec. 2008. BG. added RunFlag logic per original methodology
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// To simulate the operation of a single-speed fan fluid cooler.

		// METHODOLOGY EMPLOYED:
		// The fluid cooler is modeled using effectiveness-NTU relationships for
		// cross flow heat exchangers (both stream unmixed)based on cooling tower model.
		// The subroutine calculates the period of time required to meet a
		// leaving water temperature setpoint. It assumes that part-load
		// operation represents a linear interpolation of two steady-state regimes.
		// Cyclic losses are neglected. The period of time required to meet the
		// leaving water temperature setpoint is used to determine the required
		// fan power and energy.
		// A RunFlag is passed by the upper level manager to indicate the ON/OFF status,
		// or schedule, of the fluid cooler. If the fluid cooler is OFF, outlet water
		// temperature and flow rate are passed through the model from inlet node to
		// outlet node without intervention. Reports are also updated with fan power
		// and energy being zero.
		// When the RunFlag indicates an ON condition for thefluid cooler, the
		// mass flow rate and water temperature are read from the inlet node of the
		// fluid cooler (water-side). The outdoor air dry-bulb temperature is used
		// as the entering condition to thefluid cooler (air-side).Thefluid cooler
		// fan is turned on and design parameters are used to calculate the leaving
		// water temperature.If the calculated leaving water temperature is below the setpoint,
		// a fan run-time fraction is calculated and used to determine fan power. The leaving
		// water temperature setpoint is placed on the outlet node. If the calculated
		// leaving water temperature is at or above the setpoint, the calculated
		// leaving water temperature is placed on the outlet node and the fan runs at
		// full power. Water mass flow rate is passed from inlet node to outlet node
		// with no intervention.

		// REFERENCES:
		// ASHRAE HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculation. 1999.
		// Based on SingleSpeedTower subroutine by Dan Fisher ,Sept 1998.

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "SingleSpeedFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 AirFlowRate;
		Real64 UAdesign; // UA value at design conditions (entered by user or calculated)
		Real64 OutletWaterTempOFF;
		Real64 FanModeFrac( 0.0 );
		Real64 FanPowerOn;
		Real64 CpWater;
		Real64 TempSetPoint;

		//set inlet and outlet nodes
		int WaterInletNode = this->WaterInletNodeNum;
		Qactual = 0.0;
		FanPower = 0.0;
		OutletWaterTemp = Node( WaterInletNode ).Temp;
		{ auto const SELECT_CASE_var( PlantLoop( this->location.loopNum ).LoopDemandCalcScheme );
		if ( SELECT_CASE_var == DataPlant::SingleSetPoint ) {
			TempSetPoint = PlantLoop( this->location.loopNum ).LoopSide( this->location.loopSideNum ).TempSetPoint;
		} else if ( SELECT_CASE_var == DataPlant::DualSetPointDeadBand ) {
			TempSetPoint = PlantLoop( this->location.loopNum ).LoopSide( this->location.loopSideNum ).TempSetPointHi;
		}}

		//   MassFlowTol is a parameter to indicate a no flow condition
		if ( WaterMassFlowRate <= DataBranchAirLoopPlant::MassFlowTolerance ) return;

		if ( OutletWaterTemp < TempSetPoint ) { //already there don't need to run the cooler
			return;
		}

		//   Initialize local variables
		OutletWaterTempOFF = Node( WaterInletNode ).Temp;
		OutletWaterTemp = OutletWaterTempOFF;

		UAdesign = this->HighSpeedFluidCoolerUA;
		AirFlowRate = this->HighSpeedAirFlowRate;
		FanPowerOn = this->HighSpeedFanPower;

		this->simSimpleFluidCooler( WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp );

		if ( OutletWaterTemp <= TempSetPoint ) {
			//   Setpoint was met with pump ON and fan ON, calculate run-time fraction or just wasn't needed at all
			if ( OutletWaterTemp != OutletWaterTempOFF ) { // don't divide by zero
				FanModeFrac = ( TempSetPoint - OutletWaterTempOFF ) / ( OutletWaterTemp - OutletWaterTempOFF );
			}
			FanPower = max( FanModeFrac * FanPowerOn, 0.0 ); // BG change
			OutletWaterTemp = TempSetPoint;
		} else {
			//    Setpoint was not met, fluid cooler ran at full capacity
			FanPower = FanPowerOn;
		}
		CpWater = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, Node( WaterInletNode ).Temp,
																		 PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
		Qactual = WaterMassFlowRate * CpWater * ( Node( WaterInletNode ).Temp - OutletWaterTemp );

	}

	void
	FluidCooler::twoSpeedFluidCooler()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   August 2008
		//       MODIFIED       Dec. 2008. BG. added RunFlag logic per original methodology
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// To simulate the operation of a fluid cooler with a two-speed fan.

		// METHODOLOGY EMPLOYED:
		// The fluid cooler is modeled using effectiveness-NTU relationships for
		// cross flow heat exchangers (both stream unmixed)based on cooling tower model.
		// The subroutine calculates the period of time required to meet a
		// leaving water temperature setpoint. It assumes that part-load
		// operation represents a linear interpolation of two steady-state regimes
		// (high-speed fan operation and low-speed fan operation).
		// Cyclic losses are neglected. The period of time required to meet the
		// leaving water temperature setpoint is used to determine the required
		// fan power and energy.
		// A RunFlag is passed by the upper level manager to indicate the ON/OFF status,
		// or schedule, of the fluid cooler. If the fluid cooler is OFF, outlet water
		// temperature and flow rate are passed through the model from inlet node to
		// outlet node without intervention.Reports are also updated with fan power
		// and fan energy being zero.
		// When the RunFlag indicates an ON condition for the fluid cooler, the
		// mass flow rate and water temperature are read from the inlet node of the
		// fluid cooler (water-side). The outdoor air dry-bulb temperature is used
		// as the entering condition to the fluid cooler (air-side). Input deck
		// parameters are read for the low fan speed and a leaving water temperature
		// is calculated.
		// If the calculated leaving water temperature is below the setpoint,
		// a fan run-time fraction (FanModeFrac) is calculated and used to determine fan power.
		// The leaving water temperature setpoint is placed on the outlet node.
		// If the calculated leaving water temperature is at or above
		// the setpoint, the fluid cooler fan is turned on 'high speed' and the routine is
		// repeated. If the calculated leaving water temperature is below the setpoint,
		// a fan run-time fraction is calculated for the second stage fan and fan power
		// is calculated as FanModeFrac*HighSpeedFanPower+(1-FanModeFrac)*LowSpeedFanPower.
		// If the calculated leaving water temperature is above the leaving water temp.
		// setpoint, the calculated leaving water temperature is placed on the outlet
		// node and the fan runs at full power (High Speed Fan Power). Water mass flow
		// rate is passed from inlet node to outlet node with no intervention.

		// REFERENCES:
		// ASHRAE HVAC1KIT: A Toolkit for Primary HVAC System Energy Calculation. 1999.
		// Based on TwoSpeedTower by Dan Fisher ,Sept. 1998.

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "TwoSpeedFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 AirFlowRate;
		Real64 UAdesign; // UA value at design conditions (entered by user) [W/C]
		Real64 OutletWaterTempOFF;
		Real64 OutletWaterTemp1stStage;
		Real64 OutletWaterTemp2ndStage;
		Real64 FanModeFrac( 0.0 );
		Real64 FanPowerLow;
		Real64 FanPowerHigh;
		Real64 CpWater;
		Real64 TempSetPoint;
		int LoopNum;
		int LoopSideNum;

		int WaterInletNode = this->WaterInletNodeNum;
		Qactual = 0.0;
		FanPower = 0.0;
		OutletWaterTemp = Node( WaterInletNode ).Temp;
		LoopNum = this->location.loopNum;
		LoopSideNum = this->location.loopSideNum;
		{ auto const SELECT_CASE_var( PlantLoop( LoopNum ).LoopDemandCalcScheme );
		if ( SELECT_CASE_var == DataPlant::SingleSetPoint ) {
			TempSetPoint = PlantLoop( LoopNum ).LoopSide( LoopSideNum ).TempSetPoint;
		} else if ( SELECT_CASE_var == DataPlant::DualSetPointDeadBand ) {
			TempSetPoint = PlantLoop( LoopNum ).LoopSide( LoopSideNum ).TempSetPointHi;
		}}

		// MassFlowTol is a parameter to indicate a no flow condition
		if ( WaterMassFlowRate <= DataBranchAirLoopPlant::MassFlowTolerance || PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 0 ) return;

		// set local variable for fluid cooler
		WaterMassFlowRate = Node( WaterInletNode ).MassFlowRate;
		OutletWaterTempOFF = Node( WaterInletNode ).Temp;
		OutletWaterTemp1stStage = OutletWaterTempOFF;
		OutletWaterTemp2ndStage = OutletWaterTempOFF;
		FanModeFrac = 0.0;

		if ( OutletWaterTempOFF < TempSetPoint ) { //already there don't need to run the cooler
			return;
		}

		UAdesign = this->LowSpeedFluidCoolerUA;
		AirFlowRate = this->LowSpeedAirFlowRate;
		FanPowerLow = this->LowSpeedFanPower;

		this->simSimpleFluidCooler( WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp1stStage );

		if ( OutletWaterTemp1stStage <= TempSetPoint ) {
			// Setpoint was met with pump ON and fan ON 1st stage, calculate fan mode fraction
			if ( OutletWaterTemp1stStage != OutletWaterTempOFF ) { // don't divide by zero
				FanModeFrac = ( TempSetPoint - OutletWaterTempOFF ) / ( OutletWaterTemp1stStage - OutletWaterTempOFF );
			}
			FanPower = FanModeFrac * FanPowerLow;
			OutletWaterTemp = TempSetPoint;
			Qactual *= FanModeFrac;
		} else {
			// Setpoint was not met, turn on fluid cooler 2nd stage fan
			UAdesign = this->HighSpeedFluidCoolerUA;
			AirFlowRate = this->HighSpeedAirFlowRate;
			FanPowerHigh = this->HighSpeedFanPower;

			this->simSimpleFluidCooler( WaterMassFlowRate, AirFlowRate, UAdesign, OutletWaterTemp2ndStage );

			if ( ( OutletWaterTemp2ndStage <= TempSetPoint ) && UAdesign > 0.0 ) {
				// Setpoint was met with pump ON and fan ON 2nd stage, calculate fan mode fraction
				FanModeFrac = ( TempSetPoint - OutletWaterTemp1stStage ) / ( OutletWaterTemp2ndStage - OutletWaterTemp1stStage );
				FanPower = max( ( FanModeFrac * FanPowerHigh ) + ( 1.0 - FanModeFrac ) * FanPowerLow, 0.0 );
				OutletWaterTemp = TempSetPoint;
			} else {
				// Setpoint was not met, fluid cooler ran at full capacity
				OutletWaterTemp = OutletWaterTemp2ndStage;
				FanPower = FanPowerHigh;
			}

		}
		CpWater = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, Node( WaterInletNode ).Temp,
																		 PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );
		Qactual = WaterMassFlowRate * CpWater * ( Node( WaterInletNode ).Temp - OutletWaterTemp );

	}

	void
	FluidCooler::simSimpleFluidCooler(
		Real64 const WaterMassFlowRate,
		Real64 const AirFlowRate,
		Real64 const UAdesign,
		Real64 & OutletWaterTemp
	)
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   August 2008
		//       MODIFIED       April 2010, Chandan Sharma, FSEC
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// See purpose for Single Speed or Two Speed Fluid Cooler model

		// METHODOLOGY EMPLOYED:
		// See methodology for Single Speed or Two Speed Fluid Cooler model

		// Locals
		Real64 InletWaterTemp; // Water inlet temperature
		Real64 Qactual( 0.0 ); // Actual heat transfer rate between fluid cooler water and air [W]

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "SimSimpleFluidCooler" );

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 MdotCpWater; // Water mass flow rate times the heat capacity [W/K]
		Real64 InletAirTemp; // Dry-bulb temperature of air entering the fluid cooler [C]
		Real64 CpWater; // Heat capacity of water [J/kg/K]
		Real64 CpAir; // Heat capacity of air [J/kg/K]
		Real64 AirDensity; // Density of air [kg/m3]
		Real64 AirMassFlowRate; // Mass flow rate of air [kg/s]
		Real64 effectiveness; // Effectiveness of the heat exchanger [-]
//		Real64 OutletAirTemp; // Drybulb temp of exiting moist air [C]
		Real64 AirCapacity; // MdotCp of air through the fluid cooler
		Real64 CapacityRatioMin; // Minimum capacity of airside and waterside
		Real64 CapacityRatioMax; // Maximum capacity of airside and waterside
		Real64 CapacityRatio; // Ratio of minimum to maximum capacity
		Real64 NumTransferUnits; // Number of transfer Units [NTU]
		Real64 ETA; // initialize some local variables
		Real64 A; // initialize some local variables

		// set local fluid cooler inlet and outlet temperature variables
		InletWaterTemp = this->WaterTemp;
		OutletWaterTemp = InletWaterTemp;
		InletAirTemp = this->AirTemp;

		if ( UAdesign == 0.0 ) return;

		// set water and air properties
		AirDensity = Psychrometrics::PsyRhoAirFnPbTdbW( this->AirPress, InletAirTemp, this->AirHumRat );
		AirMassFlowRate = AirFlowRate * AirDensity;
		CpAir = Psychrometrics::PsyCpAirFnWTdb( this->AirHumRat, InletAirTemp );
		CpWater = GetSpecificHeatGlycol( PlantLoop( this->location.loopNum ).FluidName, InletWaterTemp,
																		 PlantLoop( this->location.loopNum ).FluidIndex, RoutineName );

		// Calcluate mass flow rates
		MdotCpWater = WaterMassFlowRate * CpWater;
		AirCapacity = AirMassFlowRate * CpAir;

		// calculate the minimum to maximum capacity ratios of airside and waterside
		CapacityRatioMin = min( AirCapacity, MdotCpWater );
		CapacityRatioMax = max( AirCapacity, MdotCpWater );
		CapacityRatio = CapacityRatioMin / CapacityRatioMax;

		// Calculate number of transfer units (NTU)
		NumTransferUnits = UAdesign / CapacityRatioMin;
		ETA = std::pow( NumTransferUnits, 0.22 );
		A = CapacityRatio * NumTransferUnits / ETA;
		effectiveness = 1.0 - std::exp( ( std::exp( -A ) - 1.0 ) / ( CapacityRatio / ETA ) );

		// calculate water to air heat transfer
		Qactual = effectiveness * CapacityRatioMin * ( InletWaterTemp - InletAirTemp );

		// calculate new exiting dry bulb temperature of airstream
//		OutletAirTemp = InletAirTemp + Qactual / AirCapacity;

		if ( Qactual >= 0.0 ) {
			OutletWaterTemp = InletWaterTemp - Qactual / MdotCpWater;
		} else {
			OutletWaterTemp = InletWaterTemp;
		}

	}

	Real64
	FluidCooler::simpleFluidCoolerUAResidual(
		Real64 const UA, // UA of fluid cooler
		Array1< Real64 > const & Par // par(1) = design fluid cooler load [W]
	)
	{

		// FUNCTION INFORMATION:
		//       AUTHOR         Chandan Sharma
		//       DATE WRITTEN   August 2008
		//       MODIFIED       na
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS FUNCTION:
		// Calculates residual function (Design fluid cooler load - fluid cooler Output) / Design fluid cooler load.
		// Fluid cooler output depends on the UA which is being varied to zero the residual.

		// METHODOLOGY EMPLOYED:
		// Puts UA into the fluid cooler data structure, calls SimSimpleFluidCooler, and calculates
		// the residual as defined above.

		// REFERENCES:
		// Based on SimpleTowerUAResidual by Fred Buhl, May 2002

		// Return value
		Real64 Residuum; // residual to be minimized to zero

		// FUNCTION LOCAL VARIABLE DECLARATIONS:
		Real64 OutWaterTemp; // outlet water temperature [C]
		Real64 Output; // Fluid cooler  output [W]

		this->simSimpleFluidCooler( Par( 2 ), Par( 3 ), UA, OutWaterTemp );
		Output = Par( 4 ) * Par( 2 ) * ( this->WaterTemp - OutWaterTemp );
		Residuum = ( Par( 1 ) - Output ) / Par( 1 );
		return Residuum;
	}

	void
	FluidCooler::update()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    August 2008
		//       MODIFIED         na
		//       RE-ENGINEERED    na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for passing results to the outlet water node.

		using General::TrimSigDigits;

		static gio::Fmt LowTempFmt( "(' ',F6.2)" );

		std::string CharErrOut;
		std::string CharLowOutletTemp;
		int LoopNum;
		int LoopSideNum;
		Real64 LoopMinTemp;

		// set node information

		Node( this->WaterOutletNodeNum ).Temp = OutletWaterTemp;

		LoopNum = this->location.loopNum;
		LoopSideNum = this->location.loopSideNum;
		if ( PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 0 || DataGlobals::WarmupFlag ) return;

		//Check flow rate through fluid cooler and compare to design flow rate, show warning if greater than Design * Mulitplier
		if ( Node( this->WaterOutletNodeNum ).MassFlowRate > this->DesWaterMassFlowRate * this->FluidCoolerMassFlowRateMultiplier ) {
			++this->HighMassFlowErrorCount;
			if ( this->HighMassFlowErrorCount < 2 ) {
				ShowWarningError( this->FluidCoolerType + " \"" + this->Name + "\"" );
				ShowContinueError( " Condenser Loop Mass Flow Rate is much greater than the fluid coolers design mass flow rate." );
				ShowContinueError( " Condenser Loop Mass Flow Rate = " + TrimSigDigits( Node( this->WaterOutletNodeNum ).MassFlowRate, 6 ) );
				ShowContinueError( " Fluid Cooler Design Mass Flow Rate   = " + TrimSigDigits( this->DesWaterMassFlowRate, 6 ) );
				ShowContinueErrorTimeStamp( "" );
			} else {
				ShowRecurringWarningErrorAtEnd( this->FluidCoolerType + " \"" + this->Name + "\"  Condenser Loop Mass Flow Rate is much greater than the fluid coolers design mass flow rate error continues...", this->HighMassFlowErrorIndex, Node( this->WaterOutletNodeNum ).MassFlowRate, Node( this->WaterOutletNodeNum ).MassFlowRate );
			}
		}

		// Check if OutletWaterTemp is below the minimum condenser loop temp and warn user
		LoopMinTemp = PlantLoop( LoopNum ).MinTemp;
		if ( OutletWaterTemp < LoopMinTemp && WaterMassFlowRate > 0.0 ) {
			++this->OutletWaterTempErrorCount;
			gio::write( CharLowOutletTemp, LowTempFmt ) << LoopMinTemp;
			gio::write( CharErrOut, LowTempFmt ) << OutletWaterTemp;
			strip( CharErrOut );
			if ( this->OutletWaterTempErrorCount < 2 ) {
				ShowWarningError( this->FluidCoolerType + " \"" + this->Name + "\"" );
				ShowContinueError( " Fluid cooler water outlet temperature (" + CharErrOut + " C) is below the specified minimum condenser loop temp of " + stripped( CharLowOutletTemp ) + " C" );
				ShowContinueErrorTimeStamp( "" );
			} else {
				ShowRecurringWarningErrorAtEnd( this->FluidCoolerType + " \"" + this->Name + "\" Fluid cooler water outlet temperature is below the specified minimum condenser loop temp error continues...", this->OutletWaterTempErrorIndex, OutletWaterTemp, OutletWaterTemp );
			}
		}

		// Check if water mass flow rate is small (e.g. no flow) and warn user
		if ( WaterMassFlowRate > 0.0 && WaterMassFlowRate <= DataBranchAirLoopPlant::MassFlowTolerance ) {
			++this->SmallWaterMassFlowErrorCount;
			if ( this->SmallWaterMassFlowErrorCount < 2 ) {
				ShowWarningError( this->FluidCoolerType + " \"" + this->Name + "\"" );
				ShowContinueError( " Fluid cooler water mass flow rate near zero." );
				ShowContinueErrorTimeStamp( "" );
				ShowContinueError( "Actual Mass flow = " + TrimSigDigits( WaterMassFlowRate, 2 ) );
			} else {
				ShowRecurringWarningErrorAtEnd( this->FluidCoolerType + " \"" + this->Name + "\" Fluid cooler water mass flow rate near zero error continues...", this->SmallWaterMassFlowErrorIndex, WaterMassFlowRate, WaterMassFlowRate );
			}
		}

	}

	void
	FluidCooler::report(
		bool const RunFlag
	)
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Chandan Sharma
		//       DATE WRITTEN:    August 2008
		//       MODIFIED         na
		//       RE-ENGINEERED    na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine updates the report variables for the fluid cooler.

		if ( ! RunFlag ) {
			this->InletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
			this->OutletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
			this->Qactual = 0.0;
			this->FanPower = 0.0;
			this->FanEnergy = 0.0;
		} else {
			Real64 ReportingConstant = TimeStepSys * SecInHour;
			this->InletWaterTemp = Node( this->WaterInletNodeNum ).Temp;
			this->FanEnergy = FanPower * ReportingConstant;
		}

	}

} // EnergyPlus
