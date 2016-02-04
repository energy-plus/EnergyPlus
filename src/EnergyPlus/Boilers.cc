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
#include <ObjexxFCL/Fmath.hh>

// EnergyPlus Headers
#include <Boilers.hh>
#include <BranchNodeConnections.hh>
#include <CurveManager.hh>
#include <DataBranchAirLoopPlant.hh>
#include <DataGlobalConstants.hh>
#include <DataHVACGlobals.hh>
#include <DataIPShortCuts.hh>
#include <DataLoopNode.hh>
#include <DataPlant.hh>
#include <DataPrecisionGlobals.hh>
#include <DataSizing.hh>
#include <EMSManager.hh>
#include <FluidProperties.hh>
#include <General.hh>
#include <GlobalNames.hh>
#include <InputProcessor.hh>
#include <NodeInputManager.hh>
#include <OutputProcessor.hh>
#include <OutputReportPredefined.hh>
#include <PlantUtilities.hh>
#include <PlantComponent.hh>
#include <PlantLocation.hh>
#include <ReportSizingManager.hh>
#include <UtilityRoutines.hh>

namespace EnergyPlus {

namespace Boilers {

	// Module containing the routines dealing with the Boilers

	// MODULE INFORMATION:
	//       AUTHOR         Dan Fisher, Taecheol Kim
	//       DATE WRITTEN   1998, 2000
	//       MODIFIED       na
	//       RE-ENGINEERED  na

	// PURPOSE OF THIS MODULE:
	// Perform boiler simulation for plant simulation

	// METHODOLOGY EMPLOYED:
	// The BLAST/DOE-2 empirical model based on mfg. data

	// REFERENCES: none

	// OTHER NOTES: none

	// USE STATEMENTS:
	// Use statements for data only modules
	// Using/Aliasing
	using namespace DataLoopNode;
	using namespace DataHVACGlobals;
	using namespace DataPrecisionGlobals;
	using DataGlobals::InitConvTemp;
	using DataGlobals::SecInHour;
	using DataGlobals::DisplayExtraWarnings;
	using DataPlant::PlantLoop;
	using DataPlant::TypeOf_Boiler_Simple;
	using DataPlant::ScanPlantLoopsForObject;
	using DataBranchAirLoopPlant::ControlType_SeriesActive;
	using General::TrimSigDigits;
	using General::RoundSigDigits;

	//USE FunctionFluidProperties
	// Use statements for access to subroutines in other modules

	// Data
	// MODULE PARAMETER DEFINITIONS

	// Boiler normalized efficiency curve types
	int const Linear( 1 );
	int const BiLinear( 2 );
	int const Quadratic( 3 );
	int const BiQuadratic( 4 );
	int const Cubic( 5 );
	int const QuadraticLinear( 6 );
	int const BiCubic( 7 );
	int const TriQuadratic( 8 );

	// water temperature evaluation method
	int const BoilerTempModeNotSet( 100 );
	int const EnteringBoilerTemp( 101 );
	int const LeavingBoilerTemp( 102 );

	//Boiler flow modes
	int const FlowModeNotSet( 200 );
	int const ConstantFlow( 201 );
	int const NotModulated( 202 );
	int const LeavingSetPointModulated( 203 );

	// DERIVED TYPE DEFINITIONS

	// MODULE VARIABLE DECLARATIONS:
	int NumBoilers( 0 ); // Number of boilers
	Real64 nsvFuelUsed( 0.0 ); // W - Boiler fuel used
	Real64 nsvParasiticElecPower( 0.0 ); // W - Parasitic electrical power (e.g. forced draft fan)
	Real64 nsvBoilerLoad( 0.0 ); // W - Boiler Load
	Real64 BoilerMassFlowRate( 0.0 ); // kg/s - Boiler mass flow rate
	Real64 nsvBoilerOutletTemp( 0.0 ); // W - Boiler outlet temperature
	Real64 nsvBoilerPLR( 0.0 ); // Boiler operating part-load ratio

	Array1D_bool CheckEquipName;

	// SUBROUTINE SPECIFICATIONS FOR MODULE Boilers

	// Object Data
	Array1D< BoilerSpecs > Boiler; // boiler data - dimension to number of machines

	// MODULE SUBROUTINES:

	// Beginning of Boiler Module Driver Subroutines
	//*************************************************************************

	// Functions

	void
	clear_state()
	{
		NumBoilers = 0;
		nsvFuelUsed = 0.0;
		nsvParasiticElecPower = 0.0;
		nsvBoilerLoad = 0.0;
		BoilerMassFlowRate = 0.0;
		nsvBoilerOutletTemp = 0.0;
		nsvBoilerPLR = 0.0;
		CheckEquipName.deallocate();
		Boiler.deallocate();
	}

	PlantComponent * BoilerSpecs::factory( int const EP_UNUSED(objectType), std::string objectName ) {

		static bool GetInput( true ); // if TRUE read user input
		
		//Get Input
		if ( GetInput ) {
			GetBoilerInput();
			GetInput = false;
		}
		
		// Now look for this particular component in the list
		for ( auto & BoilerItem : Boiler ) {
			if ( BoilerItem.Name == objectName ) {
				return &BoilerItem;
			}
		}
		// If we didn't find it, fatal
		ShowFatalError( "BoilerSpecs::factory : Error getting inputs for Boiler named: " + objectName );
		// Shut up the compiler
		return nullptr;
	}
	
	void
	BoilerSpecs::simulate( 
		const PlantLocation & calledFromLocation, 
		bool const EP_UNUSED(FirstHVACIteration), 
		Real64 & CurLoad )
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         DAN FISHER
		//       DATE WRITTEN   Sep. 1998
		//       MODIFIED       May. 2000, T. Kim
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subrountine controls the boiler component simulation

		// METHODOLOGY EMPLOYED: na

		// REFERENCES: na

		// Using/Aliasing
		using InputProcessor::FindItemInList;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:
		// na
		
		// SUBROUTINE PARAMETER DEFINITIONS:
		// na

		// DERIVED TYPE DEFINITIONS
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int EquipFlowCtrl; // Flow control mode for the equipment
		bool RunFlag; // if TRUE run boiler simulation--boiler is ON

		//FLOW
		
		//Calculate Load

		//Select boiler type and call boiler model
		this->InitBoiler();
		
		EquipFlowCtrl = DataPlant::PlantLoop( calledFromLocation.loopNum ).LoopSide( calledFromLocation.loopSideNum ).Branch( calledFromLocation.branchNum ).Comp( calledFromLocation.compNum ).FlowCtrl;
		
		if ( std::abs( CurLoad ) < DataPlant::LoopDemandTol ){
			RunFlag = false;
		} else {
			RunFlag = true;
		}
		
		this->CalcBoilerModel( CurLoad, RunFlag, EquipFlowCtrl );
		this->UpdateBoilerRecords( CurLoad, RunFlag );

	}
	
	void 
	BoilerSpecs::onInitLoopEquip(
		const PlantLocation & EP_UNUSED(calledFromLocation) 
	){
		// SUBROUTINE INFORMATION:
		//       AUTHOR         DAN FISHER
		//       DATE WRITTEN   Sep. 1998
		//       MODIFIED       May. 2000, T. Kim
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// User Defined plant generic component

		// METHODOLOGY EMPLOYED:
		// This routine to be called from PlantLoopEquipment.

		// REFERENCES:
		// na

		// Using/Aliasing
		using General::TrimSigDigits;
		using EMSManager::ManageEMS;
		using PlantUtilities::RegisterPlantCompDesignFlow;
		using PlantUtilities::InitComponentNodes;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		// na

		// INTERFACE BLOCK SPECIFICATIONS:
		// na

		// DERIVED TYPE DEFINITIONS:
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		
		// this->IsThisSized = false;
		// this->IsThisSized = true;
		this->InitBoiler();
		this->SizeBoiler();	
	};
	
	void 
	BoilerSpecs::getSizingFactor( 
		Real64 & SizFac 
	){
		// SUBROUTINE INFORMATION:
		//       AUTHOR         DAN FISHER
		//       DATE WRITTEN   Sep. 1998
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// User Defined plant generic component

		// METHODOLOGY EMPLOYED:
		// This routine to be called from PlantLoopEquipment.

		// REFERENCES:
		// na
		
		SizFac = this->SizFac;
	}
			
	void 
	BoilerSpecs::getDesignCapacities( 
		const PlantLocation & EP_UNUSED(calledFromLocation), 
		Real64 & MaxLoad, 
		Real64 & MinLoad, 
		Real64 & OptLoad )
	{
		
		// SUBROUTINE INFORMATION:
		//       AUTHOR         DAN FISHER
		//       DATE WRITTEN   Sep. 1998
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// User Defined plant generic component

		// METHODOLOGY EMPLOYED:
		// This routine to be called from PlantLoopEquipment.

		// REFERENCES:
		// na
		
		// now interface sizing related values with rest of E+
		MinLoad = this->NomCap * this->MinPartLoadRat;
		MaxLoad = this->NomCap * this->MaxPartLoadRat;
		OptLoad = this->NomCap * this->OptPartLoadRat;
		
	}
	
	void
	GetBoilerInput()
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Dan Fisher
		//       DATE WRITTEN:    April 1998
		//       MODIFIED:        R. Raustad - FSEC, June 2008: added boiler efficiency curve object
		//       RE-ENGINEERED:   na

		// PURPOSE OF THIS SUBROUTINE:
		// get all boiler data from input file

		// METHODOLOGY EMPLOYED:
		// standard EnergyPlus input retrieval using input Processor

		// REFERENCES: na

		// Using/Aliasing
		using DataGlobals::AnyEnergyManagementSystemInModel;
		using namespace DataGlobalConstants;
		using InputProcessor::GetNumObjectsFound;
		using InputProcessor::GetObjectItem;
		using InputProcessor::VerifyName;
		using InputProcessor::SameString;
		using namespace DataIPShortCuts; // Data for field names, blank numerics
		using BranchNodeConnections::TestCompSet;
		using NodeInputManager::GetOnlySingleNode;
		using GlobalNames::VerifyUniqueBoilerName;
		using CurveManager::GetCurveIndex;
		using CurveManager::GetCurveType;
		using General::RoundSigDigits;
		using DataSizing::AutoSize;

		// Locals
		// PARAMETERS
		static std::string const RoutineName( "GetBoilerInput: " );

		//LOCAL VARIABLES
		int BoilerNum; // boiler identifier
		int NumAlphas; // Number of elements in the alpha array
		int NumNums; // Number of elements in the numeric array
		int IOStat; // IO Status when calling get input subroutine
		static bool ErrorsFound( false ); // Flag to show errors were found during GetInput
		bool IsNotOK; // Flag to verify name
		bool IsBlank; // Flag for blank name
		bool errFlag; // Flag to show errors were found during function call
		Array1D_string BoilerFuelTypeForOutputVariable; // used to set up report variables

		//GET NUMBER OF ALL EQUIPMENT
		cCurrentModuleObject = "Boiler:HotWater";
		NumBoilers = GetNumObjectsFound( cCurrentModuleObject );

		if ( NumBoilers <= 0 ) {
			ShowSevereError( "No " + cCurrentModuleObject + " Equipment specified in input file" );
			ErrorsFound = true;
		}

		//See if load distribution manager has already gotten the input
		if ( allocated( Boiler ) ) return;

		Boiler.allocate( NumBoilers );
		CheckEquipName.allocate( NumBoilers );
		BoilerFuelTypeForOutputVariable.allocate( NumBoilers );
		CheckEquipName = true;
		BoilerFuelTypeForOutputVariable = "";

		//LOAD ARRAYS WITH CURVE FIT Boiler DATA

		for ( BoilerNum = 1; BoilerNum <= NumBoilers; ++BoilerNum ) {
			GetObjectItem( cCurrentModuleObject, BoilerNum, cAlphaArgs, NumAlphas, rNumericArgs, NumNums, IOStat, lNumericFieldBlanks, lAlphaFieldBlanks, cAlphaFieldNames, cNumericFieldNames );

			IsNotOK = false;
			IsBlank = false;
			VerifyName( cAlphaArgs( 1 ), Boiler, BoilerNum - 1, IsNotOK, IsBlank, cCurrentModuleObject + " Name" );
			if ( IsNotOK ) {
				ErrorsFound = true;
				if ( IsBlank ) cAlphaArgs( 1 ) = "xxxxx";
			}
			VerifyUniqueBoilerName( cCurrentModuleObject, cAlphaArgs( 1 ), errFlag, cCurrentModuleObject + " Name" );
			if ( errFlag ) {
				ErrorsFound = true;
			}
			Boiler( BoilerNum ).Name = cAlphaArgs( 1 );
			Boiler( BoilerNum ).TypeNum = TypeOf_Boiler_Simple;

			{ auto const SELECT_CASE_var( cAlphaArgs( 2 ) );

			if ( ( SELECT_CASE_var == "ELECTRICITY" ) || ( SELECT_CASE_var == "ELECTRIC" ) || ( SELECT_CASE_var == "ELEC" ) ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Electric";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "ELECTRICITY" );

			} else if ( ( SELECT_CASE_var == "GAS" ) || ( SELECT_CASE_var == "NATURALGAS" ) || ( SELECT_CASE_var == "NATURAL GAS" ) ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Gas";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "NATURALGAS" );

			} else if ( SELECT_CASE_var == "DIESEL" ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Diesel";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "DIESEL" );

			} else if ( SELECT_CASE_var == "GASOLINE" ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Gasoline";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "GASOLINE" );

			} else if ( SELECT_CASE_var == "COAL" ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Coal";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "COAL" );

			} else if ( ( SELECT_CASE_var == "FUEL OIL #1" ) || ( SELECT_CASE_var == "FUELOIL#1" ) || ( SELECT_CASE_var == "FUEL OIL" ) || ( SELECT_CASE_var == "DISTILLATE OIL" ) ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "FuelOil#1";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "DISTILLATE OIL" );

			} else if ( ( SELECT_CASE_var == "FUEL OIL #2" ) || ( SELECT_CASE_var == "FUELOIL#2" ) || ( SELECT_CASE_var == "RESIDUAL OIL" ) ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "FuelOil#2";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "RESIDUAL OIL" );

			} else if ( ( SELECT_CASE_var == "PROPANE" ) || ( SELECT_CASE_var == "LPG" ) || ( SELECT_CASE_var == "PROPANEGAS" ) || ( SELECT_CASE_var == "PROPANE GAS" ) ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Propane";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "PROPANE" );

			} else if ( SELECT_CASE_var == "OTHERFUEL1" ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "OtherFuel1";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "OTHERFUEL1" );

			} else if ( SELECT_CASE_var == "OTHERFUEL2" ) {
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "OtherFuel2";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "OTHERFUEL2" );

			} else {
				ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
				ShowContinueError( "Invalid " + cAlphaFieldNames( 2 ) + '=' + cAlphaArgs( 2 ) );
				// Set to Electric to avoid errors when setting up output variables
				BoilerFuelTypeForOutputVariable( BoilerNum ) = "Electric";
				Boiler( BoilerNum ).FuelType = AssignResourceTypeNum( "ELECTRICITY" );
				ErrorsFound = true;
			}}

			Boiler( BoilerNum ).NomCap = rNumericArgs( 1 );
			if ( rNumericArgs( 1 ) == 0.0 ) {
				ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
				ShowContinueError( "Invalid " + cNumericFieldNames( 1 ) + '=' + RoundSigDigits( rNumericArgs( 1 ), 2 ) );
				ShowContinueError( "..." + cNumericFieldNames( 1 ) + " must be greater than 0.0" );
				ErrorsFound = true;
			}
			if ( Boiler( BoilerNum ).NomCap == AutoSize ) {
				Boiler( BoilerNum ).NomCapWasAutoSized = true;
			}

			Boiler( BoilerNum ).Effic = rNumericArgs( 2 );
			if ( rNumericArgs( 2 ) == 0.0 ) {
				ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
				ShowContinueError( "Invalid " + cNumericFieldNames( 2 ) + '=' + RoundSigDigits( rNumericArgs( 2 ), 3 ) );
				ShowSevereError( "..." + cNumericFieldNames( 2 ) + " must be greater than 0.0" );
				ErrorsFound = true;
			}

			{ auto const SELECT_CASE_var( cAlphaArgs( 3 ) );

			if ( SELECT_CASE_var == "ENTERINGBOILER" ) {
				Boiler( BoilerNum ).CurveTempMode = EnteringBoilerTemp;
			} else if ( SELECT_CASE_var == "LEAVINGBOILER" ) {
				Boiler( BoilerNum ).CurveTempMode = LeavingBoilerTemp;
			} else {
				Boiler( BoilerNum ).CurveTempMode = BoilerTempModeNotSet;
			}}

			Boiler( BoilerNum ).EfficiencyCurvePtr = GetCurveIndex( cAlphaArgs( 4 ) );
			if ( Boiler( BoilerNum ).EfficiencyCurvePtr > 0 ) {
				{ auto const SELECT_CASE_var( GetCurveType( Boiler( BoilerNum ).EfficiencyCurvePtr ) );
				if ( SELECT_CASE_var == "LINEAR" ) {
					Boiler( BoilerNum ).EfficiencyCurveType = Linear;
				} else if ( SELECT_CASE_var == "QUADRATIC" ) {
					Boiler( BoilerNum ).EfficiencyCurveType = Quadratic;
				} else if ( SELECT_CASE_var == "QUADRATICLINEAR" ) {
					Boiler( BoilerNum ).EfficiencyCurveType = QuadraticLinear;
				} else if ( SELECT_CASE_var == "CUBIC" ) {
					Boiler( BoilerNum ).EfficiencyCurveType = Cubic;
				} else if ( SELECT_CASE_var == "BICUBIC" ) {
					Boiler( BoilerNum ).EfficiencyCurveType = BiCubic;
				} else if ( SELECT_CASE_var == "BIQUADRATIC" ) {
					Boiler( BoilerNum ).EfficiencyCurveType = BiQuadratic;
				} else {
					ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
					ShowContinueError( "Invalid " + cAlphaFieldNames( 4 ) + '=' + cAlphaArgs( 4 ) );
					ShowContinueError( "...Curve type for " + cAlphaFieldNames( 4 ) + "  = " + GetCurveType( Boiler( BoilerNum ).EfficiencyCurvePtr ) );
					ErrorsFound = true;
				}}
			} else if ( ! lAlphaFieldBlanks( 4 ) ) {
				ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
				ShowContinueError( "Invalid " + cAlphaFieldNames( 4 ) + '=' + cAlphaArgs( 4 ) );
				ShowSevereError( "..." + cAlphaFieldNames( 4 ) + " not found." );
				ErrorsFound = true;
			}

			//if curve uses temperature, make sure water temp mode has been set
			{ auto const SELECT_CASE_var( Boiler( BoilerNum ).EfficiencyCurveType );
			if ( ( SELECT_CASE_var == BiQuadratic ) || ( SELECT_CASE_var == QuadraticLinear ) || ( SELECT_CASE_var == BiCubic ) ) { //curve uses water temperature
				if ( Boiler( BoilerNum ).CurveTempMode == BoilerTempModeNotSet ) { // throw error
					if ( ! lAlphaFieldBlanks( 3 ) ) {
						ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
						ShowContinueError( "Invalid " + cAlphaFieldNames( 3 ) + '=' + cAlphaArgs( 3 ) );
						ShowContinueError( "Boiler using curve type of " + GetCurveType( Boiler( BoilerNum ).EfficiencyCurvePtr ) + " must specify " + cAlphaFieldNames( 3 ) );
						ShowContinueError( "Available choices are EnteringBoiler or LeavingBoiler" );
					} else {
						ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
						ShowContinueError( "Field " + cAlphaFieldNames( 3 ) + " is blank" );
						ShowContinueError( "Boiler using curve type of " + GetCurveType( Boiler( BoilerNum ).EfficiencyCurvePtr ) + " must specify either EnteringBoiler or LeavingBoiler" );
					}
					ErrorsFound = true;
				}
			}}

			Boiler( BoilerNum ).TempDesBoilerOut = rNumericArgs( 3 );
			Boiler( BoilerNum ).VolFlowRate = rNumericArgs( 4 );
			if ( Boiler( BoilerNum ).VolFlowRate ==  AutoSize ) {
				Boiler( BoilerNum ).VolFlowRateWasAutoSized = true;
			}
			Boiler( BoilerNum ).MinPartLoadRat = rNumericArgs( 5 );
			Boiler( BoilerNum ).MaxPartLoadRat = rNumericArgs( 6 );
			Boiler( BoilerNum ).OptPartLoadRat = rNumericArgs( 7 );

			Boiler( BoilerNum ).TempUpLimitBoilerOut = rNumericArgs( 8 );
			// default to 99.9C if upper temperature limit is left blank.
			if ( Boiler( BoilerNum ).TempUpLimitBoilerOut <= 0.0 ) {
				Boiler( BoilerNum ).TempUpLimitBoilerOut = 99.9;
			}

			Boiler( BoilerNum ).ParasiticElecLoad = rNumericArgs( 9 );
			Boiler( BoilerNum ).SizFac = rNumericArgs( 10 );
			if ( Boiler( BoilerNum ).SizFac == 0.0 ) Boiler( BoilerNum ).SizFac = 1.0;

			Boiler( BoilerNum ).BoilerInletNodeNum = GetOnlySingleNode( cAlphaArgs( 5 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Inlet, 1, ObjectIsNotParent );
			Boiler( BoilerNum ).BoilerOutletNodeNum = GetOnlySingleNode( cAlphaArgs( 6 ), ErrorsFound, cCurrentModuleObject, cAlphaArgs( 1 ), NodeType_Water, NodeConnectionType_Outlet, 1, ObjectIsNotParent );
			TestCompSet( cCurrentModuleObject, cAlphaArgs( 1 ), cAlphaArgs( 5 ), cAlphaArgs( 6 ), "Hot Water Nodes" );

			{ auto const SELECT_CASE_var( cAlphaArgs( 7 ) );
			if ( SELECT_CASE_var == "CONSTANTFLOW" ) {
				Boiler( BoilerNum ).FlowMode = ConstantFlow;
			} else if ( SELECT_CASE_var == "VARIABLEFLOW" ) { // backward compatible, clean out eventually
				Boiler( BoilerNum ).FlowMode = LeavingSetPointModulated;
				ShowWarningError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
				ShowContinueError( "Invalid " + cAlphaFieldNames( 7 ) + '=' + cAlphaArgs( 7 ) );
				ShowContinueError( "Key choice is now called \"LeavingSetpointModulated\" and the simulation continues" );
			} else if ( SELECT_CASE_var == "LEAVINGSETPOINTMODULATED" ) {
				Boiler( BoilerNum ).FlowMode = LeavingSetPointModulated;
			} else if ( SELECT_CASE_var == "NOTMODULATED" ) {
				Boiler( BoilerNum ).FlowMode = NotModulated;
			} else {
				ShowSevereError( RoutineName + cCurrentModuleObject + "=\"" + cAlphaArgs( 1 ) + "\"," );
				ShowContinueError( "Invalid " + cAlphaFieldNames( 7 ) + '=' + cAlphaArgs( 7 ) );
				ShowContinueError( "Available choices are ConstantFlow, NotModulated, or LeavingSetpointModulated" );
				ShowContinueError( "Flow mode NotModulated is assumed and the simulation continues." );
				// We will assume variable flow if not specified
				Boiler( BoilerNum ).FlowMode = NotModulated;
			}}

		}

		if ( ErrorsFound ) {
			ShowFatalError( RoutineName + "Errors found in processing " + cCurrentModuleObject + " input." );
		}

		for ( BoilerNum = 1; BoilerNum <= NumBoilers; ++BoilerNum ) {
			SetupOutputVariable( "Boiler Heating Rate [W]", Boiler( BoilerNum ).BoilerLoad, "System", "Average", Boiler( BoilerNum ).Name );
			SetupOutputVariable( "Boiler Heating Energy [J]", Boiler( BoilerNum ).BoilerEnergy, "System", "Sum", Boiler( BoilerNum ).Name, _, "ENERGYTRANSFER", "BOILERS", _, "Plant" );
			if ( SameString( BoilerFuelTypeForOutputVariable( BoilerNum ), "Electric" ) ) {
				SetupOutputVariable( "Boiler " + BoilerFuelTypeForOutputVariable( BoilerNum ) + " Power [W]", Boiler( BoilerNum ).FuelUsed, "System", "Average", Boiler( BoilerNum ).Name );
			} else {
				SetupOutputVariable( "Boiler " + BoilerFuelTypeForOutputVariable( BoilerNum ) + " Rate [W]", Boiler( BoilerNum ).FuelUsed, "System", "Average", Boiler( BoilerNum ).Name );
			}
			SetupOutputVariable( "Boiler " + BoilerFuelTypeForOutputVariable( BoilerNum ) + " Energy [J]", Boiler( BoilerNum ).FuelConsumed, "System", "Sum", Boiler( BoilerNum ).Name, _, BoilerFuelTypeForOutputVariable( BoilerNum ), "Heating", "Boiler", "Plant" );
			SetupOutputVariable( "Boiler Inlet Temperature [C]", Boiler( BoilerNum ).BoilerInletTemp, "System", "Average", Boiler( BoilerNum ).Name );
			SetupOutputVariable( "Boiler Outlet Temperature [C]", Boiler( BoilerNum ).BoilerOutletTemp, "System", "Average", Boiler( BoilerNum ).Name );
			SetupOutputVariable( "Boiler Mass Flow Rate [kg/s]", Boiler( BoilerNum ).Mdot, "System", "Average", Boiler( BoilerNum ).Name );
			SetupOutputVariable( "Boiler Ancillary Electric Power [W]", Boiler( BoilerNum ).ParasiticElecPower, "System", "Average", Boiler( BoilerNum ).Name );
			SetupOutputVariable( "Boiler Ancillary Electric Energy [J]", Boiler( BoilerNum ).ParasiticElecConsumption, "System", "Sum", Boiler( BoilerNum ).Name, _, "ELECTRICITY", "Heating", "Boiler Parasitic", "Plant" );
			SetupOutputVariable( "Boiler Part Load Ratio []", Boiler( BoilerNum ).BoilerPLR, "System", "Average", Boiler( BoilerNum ).Name );
			if ( AnyEnergyManagementSystemInModel ) {
				SetupEMSInternalVariable( "Boiler Nominal Capacity", Boiler( BoilerNum ).Name, "[W]", Boiler( BoilerNum ).NomCap );
			}

		}

		BoilerFuelTypeForOutputVariable.deallocate();

	}

	void
	BoilerSpecs::InitBoiler()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Fred Buhl
		//       DATE WRITTEN   April 2002
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  Brent Griffith, rework for plant upgrade

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for initializations of the Boiler components

		// METHODOLOGY EMPLOYED:
		// Uses the status flags to trigger initializations.

		// REFERENCES:
		// na

		// Using/Aliasing
		using DataGlobals::BeginEnvrnFlag;
		using DataGlobals::AnyEnergyManagementSystemInModel;
		using FluidProperties::GetDensityGlycol;
		using PlantUtilities::InitComponentNodes;
		using DataPlant::TypeOf_Boiler_Simple;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using DataPlant::LoopFlowStatus_NeedyIfLoopOn;
		using DataPlant::SingleSetPoint;
		using DataPlant::DualSetPointDeadBand;
		using EMSManager::iTemperatureSetPoint;
		using EMSManager::CheckIfNodeSetPointManagedByEMS;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "InitBoiler" );

		// INTERFACE BLOCK SPECIFICATIONS
		// na

		// DERIVED TYPE DEFINITIONS
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 rho;
		bool FatalError;
		bool errFlag;
		// FLOW:

		// Init more variables
		if ( this->MyFlag ) {
			// Locate the boilers on the plant loops for later usage
			errFlag = false;
			ScanPlantLoopsForObject( this->Name, TypeOf_Boiler_Simple, this->LoopNum, this->LoopSideNum, this->BranchNum, this->CompNum, _, this->TempUpLimitBoilerOut, _, _, _, errFlag );
			if ( errFlag ) {
				ShowFatalError( "InitBoiler: Program terminated due to previous condition(s)." );
			}

			if ( ( this->FlowMode == LeavingSetPointModulated ) || ( this->FlowMode == ConstantFlow ) ) {
				// reset flow priority
				PlantLoop( this->LoopNum ).LoopSide( this->LoopSideNum ).Branch( this->BranchNum ).Comp( this->CompNum ).FlowPriority = LoopFlowStatus_NeedyIfLoopOn;
			}

			this->MyFlag = false;
		}

		if ( this->MyEnvrnFlag && BeginEnvrnFlag && ( PlantFirstSizesOkayToFinalize ) ) {
			//if ( ! PlantFirstSizeCompleted ) SizeBoiler( BoilerNum );
			rho = GetDensityGlycol( PlantLoop( this->LoopNum ).FluidName, InitConvTemp, PlantLoop( this->LoopNum ).FluidIndex, RoutineName );
			this->DesMassFlowRate = this->VolFlowRate * rho;

			InitComponentNodes( 0.0, this->DesMassFlowRate, this->BoilerInletNodeNum, this->BoilerOutletNodeNum, this->LoopNum, this->LoopSideNum, this->BranchNum, this->CompNum );

			if ( this->FlowMode == LeavingSetPointModulated ) { // check if setpoint on outlet node
				if ( ( Node( this->BoilerOutletNodeNum ).TempSetPoint == SensedNodeFlagValue ) && ( Node( this->BoilerOutletNodeNum ).TempSetPointLo == SensedNodeFlagValue ) ) {
					if ( ! AnyEnergyManagementSystemInModel ) {
						if ( ! this->ModulatedFlowErrDone ) {
							ShowWarningError( "Missing temperature setpoint for LeavingSetpointModulated mode Boiler named " + this->Name );
							ShowContinueError( "  A temperature setpoint is needed at the outlet node of a boiler in variable flow mode, use a SetpointManager" );
							ShowContinueError( "  The overall loop setpoint will be assumed for Boiler. The simulation continues ... " );
							this->ModulatedFlowErrDone = true;
						}
					} else {
						// need call to EMS to check node
						FatalError = false; // but not really fatal yet, but should be.
						CheckIfNodeSetPointManagedByEMS( this->BoilerOutletNodeNum, iTemperatureSetPoint, FatalError );
						if ( FatalError ) {
							if ( ! this->ModulatedFlowErrDone ) {
								ShowWarningError( "Missing temperature setpoint for LeavingSetpointModulated mode Boiler named " + this->Name );
								ShowContinueError( "  A temperature setpoint is needed at the outlet node of a boiler in variable flow mode" );
								ShowContinueError( "  use a Setpoint Manager to establish a setpoint at the boiler outlet node " );
								ShowContinueError( "  or use an EMS actuator to establish a setpoint at the boiler outlet node " );
								ShowContinueError( "  The overall loop setpoint will be assumed for Boiler. The simulation continues ... " );
								this->ModulatedFlowErrDone = true;
							}
						}
					}
					this->ModulatedFlowSetToLoop = true; // this is for backward compatibility and could be removed
				}
			}

			this->MyEnvrnFlag = false;
		}

		if ( ! BeginEnvrnFlag ) {
			this->MyEnvrnFlag = true;
		}

		// every iteration inits.  (most in calc routine)

		if ( ( this->FlowMode == LeavingSetPointModulated ) && this->ModulatedFlowSetToLoop ) {
			// fix for clumsy old input that worked because loop setpoint was spread.
			//  could be removed with transition, testing , model change, period of being obsolete.
			{ auto const SELECT_CASE_var( PlantLoop( this->LoopNum ).LoopDemandCalcScheme );
			if ( SELECT_CASE_var == SingleSetPoint ) {
				Node( this->BoilerOutletNodeNum ).TempSetPoint = Node( PlantLoop( this->LoopNum ).TempSetPointNodeNum ).TempSetPoint;
			} else if ( SELECT_CASE_var == DualSetPointDeadBand ) {
				Node( this->BoilerOutletNodeNum ).TempSetPointLo = Node( PlantLoop( this->LoopNum ).TempSetPointNodeNum ).TempSetPointLo;
			}}
		}

	}
	
	void
	BoilerSpecs::SizeBoiler()
	{

		// SUBROUTINE INFORMATION:
		//       AUTHOR         Fred Buhl
		//       DATE WRITTEN   April 2002
		//       MODIFIED       Nov. 2013, D. Kang, Add component sizing table entries
		//       MODIFIED       Feb. 2016, R. Zhang, Refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine is for sizing Boiler Components for which capacities and flow rates
		// have not been specified in the input.

		// METHODOLOGY EMPLOYED:
		// Obtains hot water flow rate from the plant sizing array. Calculates nominal capacity from
		// the hot water flow rate and the hot water loop design delta T.

		// REFERENCES:
		// na

		// Using/Aliasing
		using namespace DataSizing;
		using DataPlant::PlantLoop;
		using DataPlant::PlantFirstSizesOkayToFinalize;
		using DataPlant::PlantFirstSizesOkayToReport;
		using DataPlant::PlantFinalSizesOkayToReport;
		using FluidProperties::GetDensityGlycol;
		using FluidProperties::GetSpecificHeatGlycol;
		using PlantUtilities::RegisterPlantCompDesignFlow;
		using ReportSizingManager::ReportSizingOutput;
		using namespace OutputReportPredefined;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "SizeBoiler" );

		// INTERFACE BLOCK SPECIFICATIONS
		// na

		// DERIVED TYPE DEFINITIONS
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int PltSizNum( 0 ); // Plant Sizing index corresponding to CurLoopNum
		bool ErrorsFound( false ); // If errors detected in input
		std::string equipName; // Name of boiler object
		Real64 rho;
		Real64 Cp;
		Real64 tmpNomCap; // local nominal capacity cooling power
		Real64 tmpBoilerVolFlowRate; // local boiler design volume flow rate
		Real64 NomCapUser( 0.0 ); // Hardsized nominal capacity for reporting
		Real64 VolFlowRateUser( 0.0 ); // Hardsized volume flow for reporting

		tmpNomCap = this->NomCap;
		tmpBoilerVolFlowRate = this->VolFlowRate;

		PltSizNum = PlantLoop( this->LoopNum ).PlantSizNum;

		if ( PltSizNum > 0 ) {
			if ( PlantSizData( PltSizNum ).DesVolFlowRate >= SmallWaterVolFlow ) {

				rho = GetDensityGlycol( PlantLoop( this->LoopNum ).FluidName, InitConvTemp, PlantLoop( this->LoopNum ).FluidIndex, RoutineName );
				Cp = GetSpecificHeatGlycol( PlantLoop( this->LoopNum ).FluidName, this->TempDesBoilerOut, PlantLoop( this->LoopNum ).FluidIndex, RoutineName );
				tmpNomCap = Cp * rho * this->SizFac * PlantSizData( PltSizNum ).DeltaT * PlantSizData( PltSizNum ).DesVolFlowRate;
				if ( ! this->NomCapWasAutoSized ) tmpNomCap = this->NomCap;

			} else {
				if ( this->NomCapWasAutoSized ) tmpNomCap = 0.0;

			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( this->NomCapWasAutoSized ) {
					this->NomCap = tmpNomCap;
					if ( PlantFinalSizesOkayToReport ) {
						ReportSizingOutput( "Boiler:HotWater", this->Name,
							"Design Size Nominal Capacity [W]", tmpNomCap );
					}
					if ( PlantFirstSizesOkayToReport ) {
						ReportSizingOutput( "Boiler:HotWater", this->Name,
							"Initial Design Size Nominal Capacity [W]", tmpNomCap );
					}
				} else { // Hard-sized with sizing data
					if ( this->NomCap > 0.0 && tmpNomCap > 0.0 ) {
						NomCapUser = this->NomCap;
						if ( PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( "Boiler:HotWater", this->Name, "Design Size Nominal Capacity [W]", tmpNomCap, "User-Specified Nominal Capacity [W]", NomCapUser );
							if ( DisplayExtraWarnings ) {
								if ( ( std::abs( tmpNomCap - NomCapUser ) / NomCapUser ) > AutoVsHardSizingThreshold ) {
									ShowMessage( "SizeBoilerHotWater: Potential issue with equipment sizing for " + this->Name );
									ShowContinueError( "User-Specified Nominal Capacity of " + RoundSigDigits( NomCapUser, 2 ) + " [W]" );
									ShowContinueError( "differs from Design Size Nominal Capacity of " + RoundSigDigits( tmpNomCap, 2 ) + " [W]" );
									ShowContinueError( "This may, or may not, indicate mismatched component sizes." );
									ShowContinueError( "Verify that the value entered is intended and is consistent with other components." );
								}
							}
						}
						tmpNomCap = NomCapUser;
					}
				}
			}
		} else {
			if ( this->NomCapWasAutoSized && PlantFirstSizesOkayToFinalize ) {
				ShowSevereError( "Autosizing of Boiler nominal capacity requires a loop Sizing:Plant object" );
				ShowContinueError( "Occurs in Boiler object=" + this->Name );
				ErrorsFound = true;
			}
			if ( ! this->NomCapWasAutoSized && PlantFinalSizesOkayToReport
					&& ( this->NomCap > 0.0 ) ) { // Hard-sized with no sizing data
					ReportSizingOutput( "Boiler:HotWater", this->Name,
						"User-Specified Nominal Capacity [W]", this->NomCap );
			}
		}

		if ( PltSizNum > 0 ) {
			if ( PlantSizData( PltSizNum ).DesVolFlowRate >= SmallWaterVolFlow ) {
				tmpBoilerVolFlowRate = PlantSizData( PltSizNum ).DesVolFlowRate * this->SizFac;
				if ( ! this->VolFlowRateWasAutoSized ) tmpBoilerVolFlowRate = this->VolFlowRate;
			} else {
				if ( this->VolFlowRateWasAutoSized ) tmpBoilerVolFlowRate = 0.0;
			}
			if ( PlantFirstSizesOkayToFinalize ) {
				if ( this->VolFlowRateWasAutoSized ) {
					this->VolFlowRate = tmpBoilerVolFlowRate;
					if ( PlantFinalSizesOkayToReport) {
						ReportSizingOutput( "Boiler:HotWater", this->Name,
							"Design Size Design Water Flow Rate [m3/s]", tmpBoilerVolFlowRate );
					}
					if ( PlantFirstSizesOkayToReport) {
						ReportSizingOutput( "Boiler:HotWater", this->Name,
							"Initial Design Size Design Water Flow Rate [m3/s]", tmpBoilerVolFlowRate );
					}
				} else {
					if ( this->VolFlowRate > 0.0 && tmpBoilerVolFlowRate > 0.0 ) {
						VolFlowRateUser = this->VolFlowRate;
						if ( PlantFinalSizesOkayToReport ) {
							ReportSizingOutput( "Boiler:HotWater", this->Name, "Design Size Design Water Flow Rate [m3/s]", tmpBoilerVolFlowRate, "User-Specified Design Water Flow Rate [m3/s]", VolFlowRateUser );
							if ( DisplayExtraWarnings ) {
								if ( ( std::abs( tmpBoilerVolFlowRate - VolFlowRateUser ) / VolFlowRateUser ) > AutoVsHardSizingThreshold ) {
									ShowMessage( "SizeBoilerHotWater: Potential issue with equipment sizing for " + this->Name );
									ShowContinueError( "User-Specified Design Water Flow Rate of " + RoundSigDigits( VolFlowRateUser, 2 ) + " [m3/s]" );
									ShowContinueError( "differs from Design Size Design Water Flow Rate of " + RoundSigDigits( tmpBoilerVolFlowRate, 2 ) + " [m3/s]" );
									ShowContinueError( "This may, or may not, indicate mismatched component sizes." );
									ShowContinueError( "Verify that the value entered is intended and is consistent with other components." );
								}
							}
						}
						tmpBoilerVolFlowRate = VolFlowRateUser;
					}
				}
			}
		} else {
			if ( this->VolFlowRateWasAutoSized && PlantFirstSizesOkayToFinalize ) {
				ShowSevereError( "Autosizing of Boiler design flow rate requires a loop Sizing:Plant object" );
				ShowContinueError( "Occurs in Boiler object=" + this->Name );
				ErrorsFound = true;
			}
			if ( ! this->VolFlowRateWasAutoSized && PlantFinalSizesOkayToReport
					&& ( this->VolFlowRate > 0.0 ) ) { // Hard-sized with no sizing data
					ReportSizingOutput( "Boiler:HotWater", this->Name,
						"User-Specified Design Water Flow Rate [m3/s]", this->VolFlowRate );
			}
		}

		RegisterPlantCompDesignFlow( this->BoilerInletNodeNum, tmpBoilerVolFlowRate );

		if ( PlantFinalSizesOkayToReport ) {
			//create predefined report
			equipName = this->Name;
			PreDefTableEntry( pdchMechType, equipName, "Boiler:HotWater" );
			PreDefTableEntry( pdchMechNomEff, equipName, this->Effic );
			PreDefTableEntry( pdchMechNomCap, equipName, this->NomCap );
		}

		if ( ErrorsFound ) {
			ShowFatalError( "Preceding sizing errors cause program termination" );
		}

	}
	
	void
	BoilerSpecs::CalcBoilerModel(
		Real64 const MyLoad, // W - hot water demand to be met by boiler
		bool const RunFlag, // TRUE if boiler operating
		int const EquipFlowCtrl // Flow control mode for the equipment
	)
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR         Dan Fisher
		//       DATE WRITTEN   April 1999
		//       MODIFIED       May. 2000, Taecheol Kim 
		//                      Jun. 2008, R. Raustad - FSEC, added boiler efficiency curve object
		//                      Aug. 2011, B. Griffith - NREL, added switch for temperature to use in curve
		//                      Feb. 2016, R. Zhang - LBNL, refactor plant component
		//       RE-ENGINEERED  na

		// PURPOSE OF THIS SUBROUTINE:
		// This subroutine calculates the boiler fuel consumption and the associated
		// hot water demand met by the boiler

		// METHODOLOGY EMPLOYED:
		// The model is based on a single combustion efficiency (=1 for electric)
		// and a second order polynomial fit of performance data to obtain part
		// load performance

		// REFERENCES:

		// Using/Aliasing
		using DataGlobals::BeginEnvrnFlag;
		using DataGlobals::WarmupFlag;

		using FluidProperties::GetSpecificHeatGlycol;
		using DataBranchAirLoopPlant::ControlType_SeriesActive;
		using CurveManager::CurveValue;
		using General::TrimSigDigits;
		using PlantUtilities::SetComponentFlowRate;
		using DataPlant::SingleSetPoint;
		using DataPlant::DualSetPointDeadBand;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		static std::string const RoutineName( "CalcBoilerModel" );

		// DERIVED TYPE DEFINITIONS
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		Real64 BoilerEff; // boiler efficiency
		Real64 BoilerNomCap; // W - boiler nominal capacity
		Real64 BoilerMaxPLR; // boiler maximum part load ratio
		Real64 BoilerMinPLR; // boiler minimum part load ratio
		Real64 TheorFuelUse; // Theoretical (stoichiometric) fuel use
		Real64 OperPLR; // operating part load ratio
		Real64 BoilerDeltaTemp( 0.0 ); // C - boiler inlet to outlet temperature difference
		Real64 TempUpLimitBout; // C - boiler high temperature limit
		int BoilerInletNode; // Boiler inlet node number
		int BoilerOutletNode; // Boiler outlet node number
		int LoopNum; // Plant loop with boiler
		int LoopSideNum; // Plant loop side with boiler (supply, demand)
		Real64 BoilerMassFlowRateMax; // Max Design Boiler Mass Flow Rate converted from Volume Flow Rate
		Real64 ParasiticElecLoad; // Boiler parasitic electric power at full load
		Real64 EffCurveOutput; // Output of boiler efficiency curve
		Real64 Cp;

		//FLOW

		nsvBoilerLoad = 0.0;
		nsvParasiticElecPower = 0.0;
		BoilerMassFlowRate = 0.0;
		BoilerInletNode = this->BoilerInletNodeNum;
		BoilerOutletNode = this->BoilerOutletNodeNum;
		BoilerNomCap = this->NomCap;
		BoilerMaxPLR = this->MaxPartLoadRat;
		BoilerMinPLR = this->MinPartLoadRat;
		BoilerEff = this->Effic;
		TempUpLimitBout = this->TempUpLimitBoilerOut;
		BoilerMassFlowRateMax = this->DesMassFlowRate;
		ParasiticElecLoad = this->ParasiticElecLoad;
		LoopNum = this->LoopNum;
		LoopSideNum = this->LoopSideNum;

		Cp = GetSpecificHeatGlycol( PlantLoop( this->LoopNum ).FluidName, Node( BoilerInletNode ).Temp, PlantLoop( this->LoopNum ).FluidIndex, RoutineName );

		//If the specified load is 0.0 or the boiler should not run then we leave this subroutine. Before leaving
		//if the component control is SERIESACTIVE we set the component flow to inlet flow so that flow resolver
		//will not shut down the branch
		if ( MyLoad <= 0.0 || ! RunFlag ) {
			if ( EquipFlowCtrl == ControlType_SeriesActive ) BoilerMassFlowRate = Node( BoilerInletNode ).MassFlowRate;
			return;
		}

		//Set the current load equal to the boiler load
		nsvBoilerLoad = MyLoad;

		if ( PlantLoop( LoopNum ).LoopSide( LoopSideNum ).FlowLock == 0 ) {
			// Either set the flow to the Constant value or caluclate the flow for the variable volume
			if ( ( this->FlowMode == ConstantFlow ) || ( this->FlowMode == NotModulated ) ) {
				// Then find the flow rate and outlet temp
				BoilerMassFlowRate = BoilerMassFlowRateMax;
				SetComponentFlowRate( BoilerMassFlowRate, BoilerInletNode, BoilerOutletNode, this->LoopNum, this->LoopSideNum, this->BranchNum, this->CompNum );

				if ( ( BoilerMassFlowRate != 0.0 ) && ( MyLoad > 0.0 ) ) {
					BoilerDeltaTemp = nsvBoilerLoad / BoilerMassFlowRate / Cp;
				} else {
					BoilerDeltaTemp = 0.0;
				}

				nsvBoilerOutletTemp = BoilerDeltaTemp + Node( BoilerInletNode ).Temp;

			} else if ( this->FlowMode == LeavingSetPointModulated ) {
				// Calculate the Delta Temp from the inlet temp to the boiler outlet setpoint
				// Then find the flow rate and outlet temp

				{ auto const SELECT_CASE_var( PlantLoop( this->LoopNum ).LoopDemandCalcScheme );
				if ( SELECT_CASE_var == SingleSetPoint ) {
					BoilerDeltaTemp = Node( BoilerOutletNode ).TempSetPoint - Node( BoilerInletNode ).Temp;
				} else if ( SELECT_CASE_var == DualSetPointDeadBand ) {
					BoilerDeltaTemp = Node( BoilerOutletNode ).TempSetPointLo - Node( BoilerInletNode ).Temp;
				} else {
					assert( false );
				}}

				nsvBoilerOutletTemp = BoilerDeltaTemp + Node( BoilerInletNode ).Temp;

				if ( ( BoilerDeltaTemp > 0.0 ) && ( nsvBoilerLoad > 0.0 ) ) {
					BoilerMassFlowRate = nsvBoilerLoad / Cp / BoilerDeltaTemp;

					BoilerMassFlowRate = min( BoilerMassFlowRateMax, BoilerMassFlowRate );

				} else {
					BoilerMassFlowRate = 0.0;
				}
				SetComponentFlowRate( BoilerMassFlowRate, BoilerInletNode, BoilerOutletNode, this->LoopNum, this->LoopSideNum, this->BranchNum, this->CompNum );

			} //End of Constant/Variable Flow If Block

		} else { // If FlowLock is True
			// Set the boiler flow rate from inlet node and then check performance
			BoilerMassFlowRate = Node( BoilerInletNode ).MassFlowRate;

			if ( ( MyLoad > 0.0 ) && ( BoilerMassFlowRate > 0.0 ) ) { // this boiler has a heat load
				nsvBoilerLoad = MyLoad;
				if ( nsvBoilerLoad > BoilerNomCap * BoilerMaxPLR ) nsvBoilerLoad = BoilerNomCap * BoilerMaxPLR;
				if ( nsvBoilerLoad < BoilerNomCap * BoilerMinPLR ) nsvBoilerLoad = BoilerNomCap * BoilerMinPLR;
				nsvBoilerOutletTemp = Node( BoilerInletNode ).Temp + nsvBoilerLoad / ( BoilerMassFlowRate * Cp );
				BoilerDeltaTemp = nsvBoilerOutletTemp - Node( BoilerInletNode ).Temp;
			} else {
				nsvBoilerLoad = 0.0;
				nsvBoilerOutletTemp = Node( BoilerInletNode ).Temp;
			}

		} //End of the FlowLock If block

		// Limit nsvBoilerOutletTemp.  If > max temp, trip boiler off
		if ( nsvBoilerOutletTemp > TempUpLimitBout ) {
			BoilerDeltaTemp = 0.0;
			nsvBoilerLoad = 0.0;
			nsvBoilerOutletTemp = Node( BoilerInletNode ).Temp;
		}

		OperPLR = nsvBoilerLoad / BoilerNomCap;
		OperPLR = min( OperPLR, BoilerMaxPLR );
		OperPLR = max( OperPLR, BoilerMinPLR );

		// set report variable
		nsvBoilerPLR = OperPLR;

		// calculate theoretical fuel use based on nominal thermal efficiency
		TheorFuelUse = nsvBoilerLoad / BoilerEff;

		// calculate normalized efficiency based on curve object type
		if ( this->EfficiencyCurvePtr > 0 ) {
			if ( this->EfficiencyCurveType == BiQuadratic || this->EfficiencyCurveType == QuadraticLinear || this->EfficiencyCurveType == BiCubic ) {

				if ( this->CurveTempMode == EnteringBoilerTemp ) {
					EffCurveOutput = CurveValue( this->EfficiencyCurvePtr, OperPLR, Node( BoilerInletNode ).Temp );
				} else if ( this->CurveTempMode == LeavingBoilerTemp ) {
					EffCurveOutput = CurveValue( this->EfficiencyCurvePtr, OperPLR, nsvBoilerOutletTemp );
				}

			} else {
				EffCurveOutput = CurveValue( this->EfficiencyCurvePtr, OperPLR );
			}
		} else {
			EffCurveOutput = 1.0;
		}

		// warn if efficiency curve produces zero or negative results
		if ( ! WarmupFlag && EffCurveOutput <= 0.0 ) {
			if ( nsvBoilerLoad > 0.0 ) {
				if ( this->EffCurveOutputError < 1 ) {
					++this->EffCurveOutputError;
					ShowWarningError( "Boiler:HotWater \"" + this->Name + "\"" );
					ShowContinueError( "...Normalized Boiler Efficiency Curve output is less than or equal to 0." );
					ShowContinueError( "...Curve input x value (PLR)     = " + TrimSigDigits( OperPLR, 5 ) );
					if ( this->EfficiencyCurveType == BiQuadratic || this->EfficiencyCurveType == QuadraticLinear || this->EfficiencyCurveType == BiCubic ) {
						if ( this->CurveTempMode == EnteringBoilerTemp ) {
							ShowContinueError( "...Curve input y value (Tinlet) = " + TrimSigDigits( Node( BoilerInletNode ).Temp, 2 ) );
						} else if ( this->CurveTempMode == LeavingBoilerTemp ) {
							ShowContinueError( "...Curve input y value (Toutlet) = " + TrimSigDigits( nsvBoilerOutletTemp, 2 ) );
						}
					}
					ShowContinueError( "...Curve output (normalized eff) = " + TrimSigDigits( EffCurveOutput, 5 ) );
					ShowContinueError( "...Calculated Boiler efficiency  = " + TrimSigDigits( EffCurveOutput * BoilerEff, 5 ) + " (Boiler efficiency = Nominal Thermal Efficiency * Normalized Boiler Efficiency Curve output)" );
					ShowContinueErrorTimeStamp( "...Curve output reset to 0.01 and simulation continues." );
				} else {
					ShowRecurringWarningErrorAtEnd( "Boiler:HotWater \"" + this->Name + "\": Boiler Efficiency Curve output is less than or equal to 0 warning continues...", this->EffCurveOutputIndex, EffCurveOutput, EffCurveOutput );
				}
			}
			EffCurveOutput = 0.01;
		}

		// warn if overall efficiency greater than 1.1
		if ( ! WarmupFlag && EffCurveOutput * BoilerEff > 1.1 ) {
			if ( nsvBoilerLoad > 0.0 && this->EfficiencyCurvePtr > 0 ) {
				if ( this->CalculatedEffError < 1 ) {
					++this->CalculatedEffError;
					ShowWarningError( "Boiler:HotWater \"" + this->Name + "\"" );
					ShowContinueError( "...Calculated Boiler Efficiency is greater than 1.1." );
					ShowContinueError( "...Boiler Efficiency calculations shown below." );
					ShowContinueError( "...Curve input x value (PLR)     = " + TrimSigDigits( OperPLR, 5 ) );
					if ( this->EfficiencyCurveType == BiQuadratic || this->EfficiencyCurveType == QuadraticLinear || this->EfficiencyCurveType == BiCubic ) {
						if ( this->CurveTempMode == EnteringBoilerTemp ) {
							ShowContinueError( "...Curve input y value (Tinlet) = " + TrimSigDigits( Node( BoilerInletNode ).Temp, 2 ) );
						} else if ( this->CurveTempMode == LeavingBoilerTemp ) {
							ShowContinueError( "...Curve input y value (Toutlet) = " + TrimSigDigits( nsvBoilerOutletTemp, 2 ) );
						}
					}
					ShowContinueError( "...Curve output (normalized eff) = " + TrimSigDigits( EffCurveOutput, 5 ) );
					ShowContinueError( "...Calculated Boiler efficiency  = " + TrimSigDigits( EffCurveOutput * BoilerEff, 5 ) + " (Boiler efficiency = Nominal Thermal Efficiency * Normalized Boiler Efficiency Curve output)" );
					ShowContinueErrorTimeStamp( "...Curve output reset to 1.1 and simulation continues." );
				} else {
					ShowRecurringWarningErrorAtEnd( "Boiler:HotWater \"" + this->Name + "\": Calculated Boiler Efficiency is greater than 1.1 warning continues...", this->CalculatedEffIndex, EffCurveOutput * BoilerEff, EffCurveOutput * BoilerEff );
				}
			}
			EffCurveOutput = 1.1;
		}

		// calculate fuel used based on normalized boiler efficiency curve (=1 when no curve used)
		nsvFuelUsed = TheorFuelUse / EffCurveOutput;
		if ( nsvBoilerLoad > 0.0 ) nsvParasiticElecPower = ParasiticElecLoad * OperPLR;

	}

	// Beginning of Record Keeping subroutines for the BOILER:HOTWATER Module
	// *****************************************************************************

	void
	BoilerSpecs::UpdateBoilerRecords(
		Real64 const MyLoad, // boiler operating load
		bool const RunFlag // boiler on when TRUE
	)
	{
		// SUBROUTINE INFORMATION:
		//       AUTHOR:          Dan Fisher
		//       DATE WRITTEN:    October 1998

		// PURPOSE OF THIS SUBROUTINE:
		// boiler simulation reporting

		// METHODOLOGY EMPLOYED:na

		// REFERENCES: na

		// USE STATEMENTS: na

		// Using/Aliasing
		using PlantUtilities::SafeCopyPlantNode;

		// Locals
		// SUBROUTINE ARGUMENT DEFINITIONS:

		// SUBROUTINE PARAMETER DEFINITIONS:
		// na

		// INTERFACE BLOCK SPECIFICATIONS
		// na

		// DERIVED TYPE DEFINITIONS
		// na

		// SUBROUTINE LOCAL VARIABLE DECLARATIONS:
		int BoilerInletNode; // Boiler inlet node number
		int BoilerOutletNode; // Boiler outlet node number
		Real64 ReportingConstant; // constant for converting power to energy

		ReportingConstant = TimeStepSys * SecInHour;

		BoilerInletNode = this->BoilerInletNodeNum;
		BoilerOutletNode = this->BoilerOutletNodeNum;

		if ( MyLoad <= 0 || ! RunFlag ) {
			//set node temperatures
			SafeCopyPlantNode( BoilerInletNode, BoilerOutletNode );
			Node( BoilerOutletNode ).Temp = Node( BoilerInletNode ).Temp;
			this->BoilerOutletTemp = Node( BoilerInletNode ).Temp;
			this->BoilerLoad = 0.0;
			this->FuelUsed = 0.0;
			this->ParasiticElecPower = 0.0;
			this->BoilerPLR = 0.0;

		} else {
			//set node temperatures
			SafeCopyPlantNode( BoilerInletNode, BoilerOutletNode );
			Node( BoilerOutletNode ).Temp = nsvBoilerOutletTemp;
			this->BoilerOutletTemp = nsvBoilerOutletTemp;
			this->BoilerLoad = nsvBoilerLoad;
			this->FuelUsed = nsvFuelUsed;
			this->ParasiticElecPower = nsvParasiticElecPower;
			this->BoilerPLR = nsvBoilerPLR;

		}

		this->BoilerInletTemp = Node( BoilerInletNode ).Temp;
		this->Mdot = Node( BoilerOutletNode ).MassFlowRate;

		this->BoilerEnergy = this->BoilerLoad * ReportingConstant;
		this->FuelConsumed = this->FuelUsed * ReportingConstant;
		this->ParasiticElecConsumption = this->ParasiticElecPower * ReportingConstant;

	}

	// End of Record Keeping subroutines for the BOILER:HOTWATER Module
	// *****************************************************************************

} // Boilers

} // EnergyPlus
