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

// EnergyPlus::FluidCoolers Unit Tests

// Google Test Headers
#include <gtest/gtest.h>

#include "Fixtures/EnergyPlusFixture.hh"

// EnergyPlus Headers
#include <EnergyPlus/FluidCoolers.hh>
#include <DataSizing.hh>
#include <EnergyPlus/UtilityRoutines.hh>


using namespace EnergyPlus;
using namespace EnergyPlus::FluidCoolers;
using namespace DataGlobals;
using namespace EnergyPlus::DataSizing;
using namespace ObjexxFCL;

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test1 )
{

	using DataSizing::AutoSize;
	int StringArraySize = 20;
	Array1D_string cNumericFieldNames;
	cNumericFieldNames.allocate( StringArraySize );
	Array1D_string cAlphaFieldNames;
	cAlphaFieldNames.allocate( StringArraySize );
	Array1D_string AlphArray;
	AlphArray.allocate( StringArraySize );
	for ( int i = 1; i <= StringArraySize; ++i ) {
		cAlphaFieldNames( i ) = "AlphaField";
		cNumericFieldNames( i ) = "NumerField";
		AlphArray( i ) = "FieldValues";
	}
	std::string const cCurrentModuleObject( "FluidCooler:TwoSpeed" );
	int FluidCoolerNum( 1 );
	SimpleFluidCoolers.allocate( FluidCoolerNum );

	SimpleFluidCoolers( FluidCoolerNum ).Name = "Test";
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerMassFlowRateMultiplier = 2.5;
	SimpleFluidCoolers( FluidCoolerNum ).PerformanceInputMethod_Num = PIM::NominalCapacity;
	SimpleFluidCoolers( FluidCoolerNum ).WaterInletNodeNum = 1;
	SimpleFluidCoolers( FluidCoolerNum ).WaterOutletNodeNum = 1;
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerNominalCapacity = 50000;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringWaterTemp = 52;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringAirTemp = 35;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringAirWetBulbTemp = 25;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedAirFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedAirFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFanPower = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFanPowerWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedAirFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedAirFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFanPower = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFanPowerWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerLowSpeedNomCap = 30000;


	AlphArray( 4 ) = "NominalCapacity";
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUA = 0;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFluidCoolerUA = 0;
	SimpleFluidCoolers( 1 ).DesignEnteringWaterTemp = 50;
	bool testResult = FluidCooler::verifyTwoSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_FALSE( testResult ); // no error message triggered

	SimpleFluidCoolers( 1 ).DesignEnteringWaterTemp = -10;
	testResult = FluidCooler::verifyTwoSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_TRUE( testResult ); // error message triggered

	SimpleFluidCoolers( 1 ).DesignEnteringWaterTemp = 50;
	SimpleFluidCoolers( 1 ).FluidCoolerLowSpeedNomCap = AutoSize;
	SimpleFluidCoolers( 1 ).FluidCoolerLowSpeedNomCapWasAutoSized = true;
	testResult = FluidCooler::verifyTwoSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_FALSE( testResult ); // no error message triggered

	SimpleFluidCoolers( 1 ).FluidCoolerLowSpeedNomCap = 0; // this should trigger the original error condition
	SimpleFluidCoolers( 1 ).FluidCoolerLowSpeedNomCapWasAutoSized = false;
	testResult = FluidCooler::verifyTwoSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_TRUE( testResult ); // error message triggered

	SimpleFluidCoolers.deallocate();
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test2 ) {

	using DataSizing::AutoSize;
	int StringArraySize = 20;
	Array1D_string cNumericFieldNames;
	cNumericFieldNames.allocate( StringArraySize );
	Array1D_string cAlphaFieldNames;
	cAlphaFieldNames.allocate( StringArraySize );
	Array1D_string AlphArray;
	AlphArray.allocate( StringArraySize );
	for ( int i = 1; i <= StringArraySize; ++i ) {
		cAlphaFieldNames( i ) = "AlphaField";
		cNumericFieldNames( i ) = "NumerField";
		AlphArray( i ) = "FieldValues";
	}
	std::string const cCurrentModuleObject( "FluidCooler:TwoSpeed" );
	int FluidCoolerNum( 1 );
	bool ErrrorsFound( false );
	SimpleFluidCoolers.allocate( FluidCoolerNum );

	SimpleFluidCoolers( FluidCoolerNum ).Name = "Test";
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerMassFlowRateMultiplier = 1.0;
	SimpleFluidCoolers( FluidCoolerNum ).PerformanceInputMethod_Num = PIM::UFactor;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringWaterTemp = 52;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringAirTemp = 35;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringAirWetBulbTemp = 25;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedAirFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedAirFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFanPower = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFanPowerWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedAirFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedAirFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFanPower = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFanPowerWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerLowSpeedNomCap = 30000;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFluidCoolerUA = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).LowSpeedFluidCoolerUAWasAutoSized = true;

	AlphArray( 4 ) = "UFactorTimesAreaAndDesignWaterFlowRate";
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUA = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUAWasAutoSized = false;
	bool testResult = FluidCooler::verifyTwoSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_TRUE( testResult ); // error message triggered

	ErrrorsFound = false;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUA = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUAWasAutoSized = true;
	testResult = FluidCooler::verifyTwoSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_FALSE( testResult ); // no error message triggered

	SimpleFluidCoolers.deallocate();
	cNumericFieldNames.deallocate();
	cAlphaFieldNames.deallocate();
	AlphArray.deallocate();
}


TEST_F( EnergyPlusFixture, SingleSpeedFluidCoolerInput_Test3 )
{
	using DataSizing::AutoSize;
	int StringArraySize = 20;
	Array1D_string cNumericFieldNames;
	cNumericFieldNames.allocate( StringArraySize );
	Array1D_string cAlphaFieldNames;
	cAlphaFieldNames.allocate( StringArraySize );
	Array1D_string AlphArray;
	AlphArray.allocate( StringArraySize );
	for ( int i = 1; i <= StringArraySize; ++i ) {
		cAlphaFieldNames( i ) = "AlphaField";
		cNumericFieldNames( i ) = "NumerField";
		AlphArray( i ) = "FieldValues";
	}
	std::string const cCurrentModuleObject( "FluidCooler:SingleSpeed" );
	int FluidCoolerNum( 1 );
	SimpleFluidCoolers.allocate( FluidCoolerNum );

	SimpleFluidCoolers( FluidCoolerNum ).Name = "Test";
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerMassFlowRateMultiplier = 2.5;
	SimpleFluidCoolers( FluidCoolerNum ).PerformanceInputMethod_Num = PIM::UFactor;
	SimpleFluidCoolers( FluidCoolerNum ).WaterInletNodeNum = 1;
	SimpleFluidCoolers( FluidCoolerNum ).WaterOutletNodeNum = 1;
	SimpleFluidCoolers( FluidCoolerNum ).FluidCoolerNominalCapacity = 50000;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringWaterTemp = 52;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringAirTemp = 35;
	SimpleFluidCoolers( FluidCoolerNum ).DesignEnteringAirWetBulbTemp = 25;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedAirFlowRate = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedAirFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFanPower = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFanPowerWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUA = AutoSize;
	SimpleFluidCoolers( FluidCoolerNum ).HighSpeedFluidCoolerUAWasAutoSized = true;

	AlphArray( 4 ) = "UFactorTimesAreaAndDesignWaterFlowRate";
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRate = 1;
	bool testResult = FluidCooler::verifySingleSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_FALSE( testResult ); // no error message triggered

	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRateWasAutoSized = true;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRate = 0;
	testResult = FluidCooler::verifySingleSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_FALSE( testResult ); // no error message triggered

	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRateWasAutoSized = false;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRate = 1;
	testResult = FluidCooler::verifySingleSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_FALSE( testResult ); // no error message triggered

	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRateWasAutoSized = false;
	SimpleFluidCoolers( FluidCoolerNum ).DesignWaterFlowRate = 0;
	testResult = FluidCooler::verifySingleSpeedDesignInputs( SimpleFluidCoolers( FluidCoolerNum ), cCurrentModuleObject, AlphArray, cNumericFieldNames, cAlphaFieldNames );
	EXPECT_TRUE( testResult ); // error message triggered

}
