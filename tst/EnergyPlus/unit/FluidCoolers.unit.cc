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
#include <EnergyPlus/DataSizing.hh>
#include <EnergyPlus/DataPlant.hh>
#include <EnergyPlus/UtilityRoutines.hh>


using namespace EnergyPlus;
using namespace DataGlobals;
using namespace EnergyPlus::DataSizing;
using namespace ObjexxFCL;

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test1 )
{
	std::string const idf_objects = 
		R"(FluidCooler:TwoSpeed,
		Big FluidCooler,         !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		NominalCapacity,         !- Performance Input Method
		,                        !- High Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		28601.,                  !- Low Speed Nominal Capacity {W}
		,                        !- Low Speed Nominal Capacity Sizing Factor
		50,                      !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor)";

	ASSERT_FALSE( process_idf( idf_objects ) );
	EXPECT_NO_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test2 )
{
	std::string const idf_objects = 
		R"(FluidCooler:TwoSpeed,
		Big FluidCooler,        !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		NominalCapacity,         !- Performance Input Method
		,                        !- High Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		28601.,                  !- Low Speed Nominal Capacity {W}
		,                        !- Low Speed Nominal Capacity Sizing Factor
		-10,                     !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor)";

	EXPECT_TRUE( process_idf( idf_objects, false ) );
	EXPECT_ANY_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test3 )
{
	std::string const idf_objects = 
		R"(FluidCooler:TwoSpeed,
		Big FluidCooler,       !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		NominalCapacity,         !- Performance Input Method
		,                        !- High Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		autocalculate,           !- Low Speed Nominal Capacity {W}
		,                        !- Low Speed Nominal Capacity Sizing Factor
		50,                      !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor)";

	ASSERT_FALSE( process_idf( idf_objects ) );
	EXPECT_NO_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test4 )
{
	std::string const idf_objects = 
		R"(FluidCooler:TwoSpeed,
		Big FluidCooler,       !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		NominalCapacity,         !- Performance Input Method
		,                        !- High Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		0,                       !- Low Speed Nominal Capacity {W}
		,                        !- Low Speed Nominal Capacity Sizing Factor
		50,                      !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor)";

	EXPECT_TRUE( process_idf( idf_objects, false ) );
	EXPECT_ANY_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test5 )
{
	std::string const idf_objects = 
		R"(FluidCooler:TwoSpeed,
		Big FluidCooler,       !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		UFactorTimesAreaAndDesignWaterFlowRate,         !- Performance Input Method
		Autosize,                !- High Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		28601.,                  !- Low Speed Nominal Capacity {W}
		1.0,                     !- Low Speed Nominal Capacity Sizing Factor
		50,                      !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor)";

	EXPECT_FALSE( process_idf( idf_objects ) );
	EXPECT_ANY_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test6 )
{
	std::string const idf_objects = 
	R"(
		FluidCooler:TwoSpeed,
		Big FluidCooler,         !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		NominalCapacity,         !- Performance Input Method
		,                        !- High Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		autocalculate,           !- Low Speed Nominal Capacity {W}
		0.5,                     !- Low Speed Nominal Capacity Sizing Factor
		51.67,                   !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor
	)";

	ASSERT_FALSE( process_idf( idf_objects ) );
	EXPECT_NO_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, TwoSpeedFluidCoolerInput_Test7 )
{
	std::string const idf_objects = 
	R"(
		FluidCooler:TwoSpeed,
		Big FluidCooler,         !- Name
		Condenser FluidCooler Inlet Node,  !- Water Inlet Node Name
		Condenser FluidCooler Outlet Node,  !- Water Outlet Node Name
		UFactorTimesAreaAndDesignWaterFlowRate,  !- Performance Input Method
		autosize,                !- High Fan Speed U-factor Times Area Value {W/K}
		autosize,                !- Low Fan Speed U-factor Times Area Value {W/K}
		,                        !- Low Fan Speed U-Factor Times Area Sizing Factor
		58601.,                  !- High Speed Nominal Capacity {W}
		28601.,                  !- Low Speed Nominal Capacity {W}
		,                        !- Low Speed Nominal Capacity Sizing Factor
		51.67,                   !- Design Entering Water Temperature {C}
		35,                      !- Design Entering Air Temperature {C}
		25.6,                    !- Design Entering Air Wet-bulb Temperature {C}
		Autosize,                !- Design Water Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Air Flow Rate {m3/s}
		Autosize,                !- High Fan Speed Fan Power {W}
		autocalculate,           !- Low Fan Speed Air Flow Rate {m3/s}
		,                        !- Low Fan Speed Air Flow Rate Sizing Factor
		autocalculate,           !- Low Fan Speed Fan Power {W}
		;                        !- Low Fan Speed Fan Power Sizing Factor
	)";

	ASSERT_FALSE( process_idf( idf_objects ) );
	EXPECT_NO_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_TwoSpd, "BIG FLUIDCOOLER" ) );
}

TEST_F( EnergyPlusFixture, SingleSpeedFluidCoolerInput_Test1 )
{
	std::string const idf_objects = 
	R"(
		FluidCooler:SingleSpeed,
		Fluid Cooler,            !- Name
		Fluid Cooler Water Inlet Node,  !- Water Inlet Node Name
		Fluid Cooler Water Outlet Node,  !- Water Outlet Node Name
		UFactorTimesAreaAndDesignWaterFlowRate,  !- Performance Input Method
		autosize,                !- Design Air Flow Rate U-factor Times Area Value {W/K}
		,                        !- Nominal Capacity {W}
		52.000000,               !- Design Entering Water Temperature {C}
		24.000000,               !- Design Entering Air Temperature {C}
		16.000000,               !- Design Entering Air Wetbulb Temperature {C}
		autosize,                !- Design Water Flow Rate {m3/s}
		autosize,                !- Design Air Flow Rate {m3/s}
		autosize;                !- Design Air Flow Rate Fan Power {W}
	)";

	ASSERT_FALSE( process_idf( idf_objects ) );
	EXPECT_NO_THROW( FluidCooler::factory( DataPlant::TypeOf_FluidCooler_SingleSpd, "FLUID COOLER" ) );
}
