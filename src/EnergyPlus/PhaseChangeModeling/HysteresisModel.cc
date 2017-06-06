// EnergyPlus, Copyright (c) 1996-2017, The Board of Trustees of the University of Illinois and The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.
// 
// NOTICE: This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit others to do so.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 
// (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// 
// (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, the University of Illinois, U.S. Dept. of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 
// (4) Use of EnergyPlus(TM) Name. If Licensee (i) distributes the software in stand-alone form without changes from the version obtained under this License, or (ii) Licensee makes a reference solely to the software portion of its product, Licensee must refer to the software as "EnergyPlus version X" software, where "X" is the version number Licensee obtained under this License and may not use a different name for the software. Except as specifically required in this Section (4), Licensee shall not use in a company name, a product name, in advertising, publicity, or other promotional activities any name, trade name, trademark, logo, or other designation of "EnergyPlus", "E+", "e+" or confusingly similar designation, without the U.S. Department of Energy's prior written consent.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <ObjexxFCL/Array1D.hh>

#include <DataHeatBalance.hh>
#include <EnergyPlus.hh>
#include <PhaseChangeModeling/HysteresisModel.hh>
#include <DataIPShortCuts.hh>
#include <UtilityRoutines.hh>
#include <InputProcessor.hh>

namespace EnergyPlus {

namespace HysteresisPhaseChange {

    struct PhaseChangeStates {
        // keeping these as ints to allow output variable reporting; could refine later into enum class
        static const int LIQUID = -2;
		static const int MELTING = -1;
		static const int TRANSITION = 0;
		static const int FREEZING = 1;
		static const int CRYSTALLIZED = 2;
    };

	bool getHysteresisModels( true );
	int numHysteresisModels = 0;
	std::vector< HysteresisPhaseChange > hysteresisPhaseChangeModels;

	PhaseChangeModel * HysteresisPhaseChange::factory( const std::string objectName ) {
		if ( getHysteresisModels ) {
			readAllHysteresisModels();
			getHysteresisModels = false;
		}
		for ( auto & hm : hysteresisPhaseChangeModels ) {
			if ( hm.name == objectName ) {
				return &hm;
			}
		}
		// because of the passive linking between materials and material property objects,
		// we don't know ahead of time for sure whether we will have a material property
		// so we can't return fatal here if it isn't found, just leave it null
		return nullptr; // just for the compiler warning
	}

	Real64 HysteresisPhaseChange::getEnthalpy( Real64 T, Real64 Tc, Real64 tau1, Real64 tau2, Real64 deltaH, Real64 CpSolid, Real64 CpLiquid ) {
		Real64 eta1 = ( deltaH / 2 ) * exp( -2 * abs( T - Tc ) / tau1 );
		Real64 eta2 = ( deltaH / 2 ) * exp( -2 * abs( T - Tc ) / tau2 );
		if ( T <= Tc ) {
			return ( CpSolid * T ) + eta1;
		} else {
			return ( CpSolid * Tc ) + deltaH + CpLiquid * ( T - Tc ) - eta2;
		}
	}

	Real64 HysteresisPhaseChange::getCurrentSpecificHeat( Real64 prevTempTD, Real64 updatedTempTDT, int prevPhaseChangeState, int & phaseChangeState ) {

		Real64 TempLowPCM = this->peakTempMelting - this->deltaTempMeltingLow;
		Real64 TempHighPCM = this->peakTempMelting + this->deltaTempMeltingHigh;
		Real64 Tau1;  // assigned later
		Real64 Tau2;  // assigned later
		Real64 TempLowPCF = this->peakTempFreezing - this->deltaTempFreezingLow;
		Real64 TempHighPCF = this->peakTempFreezing + this->deltaTempFreezingHigh;
		Real64 DeltaH = this->totalLatentHeat;
		Real64 Cp, Tc;
	        Real64 phaseChangeDeltaT = prevTempTD - updatedTempTDT;

		// this is pulled directly from a chunk of the Fortran PCM code changes
		if ( phaseChangeDeltaT <= 0 ) {
			if ( updatedTempTDT < TempLowPCM ) {
				phaseChangeState = PhaseChangeStates::CRYSTALLIZED;
				Tc = this->peakTempMelting;
				Tau1 = this->deltaTempMeltingLow;
				Tau2 = this->deltaTempMeltingHigh;
			} else if ( updatedTempTDT >= TempLowPCM && updatedTempTDT <= TempHighPCM ) {
				phaseChangeState = PhaseChangeStates::MELTING;
				Tc = this->peakTempMelting;
				Tau1 = this->deltaTempMeltingLow;
				Tau2 = this->deltaTempMeltingHigh;
				if ( ( prevPhaseChangeState == PhaseChangeStates::FREEZING && phaseChangeState == PhaseChangeStates::MELTING ) || ( prevPhaseChangeState == PhaseChangeStates::TRANSITION && phaseChangeState == PhaseChangeStates::MELTING ) )
				{
					phaseChangeState = PhaseChangeStates::TRANSITION;
				}
			} else if ( updatedTempTDT > TempHighPCM ) {
				phaseChangeState = PhaseChangeStates::LIQUID;
				Tc = this->peakTempMelting;
				Tau1 = this->deltaTempMeltingLow;
				Tau2 = this->deltaTempMeltingHigh;
			}
		} else if ( phaseChangeDeltaT > 0 ) {
			if ( updatedTempTDT < TempLowPCF ) {
				phaseChangeState = PhaseChangeStates::CRYSTALLIZED;
				Tc = this->peakTempFreezing;
				Tau1 = this->deltaTempFreezingLow;
				Tau2 = this->deltaTempFreezingHigh;
			} else if ( updatedTempTDT >= TempLowPCF && updatedTempTDT <= TempHighPCF ) {
				phaseChangeState = PhaseChangeStates::FREEZING;
				Tc = this->peakTempFreezing;
				Tau1 = this->deltaTempFreezingLow;
				Tau2 = this->deltaTempFreezingHigh;
			}
			if ( ( prevPhaseChangeState == PhaseChangeStates::MELTING && phaseChangeState == PhaseChangeStates::FREEZING ) || ( prevPhaseChangeState == PhaseChangeStates::TRANSITION && phaseChangeState == PhaseChangeStates::FREEZING ) ) {
				phaseChangeState = PhaseChangeStates::TRANSITION;
			} else if ( updatedTempTDT > TempHighPCF ) {
				phaseChangeState = PhaseChangeStates::LIQUID;
				Tc = this->peakTempFreezing;
				Tau1 = this->deltaTempFreezingLow;
				Tau2 = this->deltaTempFreezingHigh;
			}
		}
		if ( prevPhaseChangeState == PhaseChangeStates::TRANSITION && phaseChangeState == PhaseChangeStates::CRYSTALLIZED ) {
			this->phaseChangeTransition = PhaseChangeStates::FREEZING;
		} else if ( prevPhaseChangeState == PhaseChangeStates::TRANSITION && phaseChangeState == PhaseChangeStates::FREEZING ) {
			this->phaseChangeTransition = PhaseChangeStates::FREEZING;
			// this->phaseChangeState = 0; ?????
		} else if ( prevPhaseChangeState == PhaseChangeStates::FREEZING && phaseChangeState == PhaseChangeStates::TRANSITION ) {
			this->phaseChangeTransition = PhaseChangeStates::FREEZING;
		} else if ( prevPhaseChangeState == PhaseChangeStates::CRYSTALLIZED && phaseChangeState == PhaseChangeStates::TRANSITION ) {
			this->phaseChangeTransition = PhaseChangeStates::FREEZING;
		} else {
			this->phaseChangeTransition = PhaseChangeStates::TRANSITION;
		}

		// if ( hysteresis flag == 1 )  -- implied by this derived class
		if ( this->phaseChangeTransition == 0 ) {
			this->enthOld = this->getEnthalpy(prevTempTD, Tc, Tau1, Tau2, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
			this->enthNew = this->getEnthalpy(updatedTempTDT, Tc, Tau1, Tau2, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
		} else if ( this->phaseChangeTransition == 1 ) {
			if ( prevPhaseChangeState == 1 && phaseChangeState == 0 ) {
				this->enthRev = this->getEnthalpy(this->TR, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthOld - ( this->specHeatTransition * prevTempTD) );
				this->enthalpyM = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthalpyF = this->getEnthalpy(updatedTempTDT, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				if ( this->enthNew < this->enthRev && this->enthNew >= this->enthalpyF && updatedTempTDT <= prevTempTD) {
					phaseChangeState = 1;
					this->enthNew = this->getEnthalpy(updatedTempTDT, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				} else if ( this->enthNew < this->enthalpyF && this->enthNew > this->enthalpyM ) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthOld - ( this->specHeatTransition * prevTempTD) );
				} else if ( this->enthNew < this->enthalpyF && updatedTempTDT > this->TR ) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				} else if ( this->enthNew <= this->enthalpyM && updatedTempTDT <= this->TR ) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				}
			} else if ( prevPhaseChangeState == 0 && phaseChangeState == 0 ) {
				if ( updatedTempTDT < this->TR ) {
					Tc = this->peakTempMelting;
					Tau1 = this->deltaTempMeltingLow;
					Tau2 = this->deltaTempMeltingHigh;
				} else if ( updatedTempTDT > this->TR ) {
					Tc = this->peakTempFreezing;
					Tau1 = this->deltaTempFreezingLow;
					Tau2 = this->deltaTempFreezingHigh;
				}
				this->enthRev = this->getEnthalpy(this->TR, Tc, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthOld - ( this->specHeatTransition * prevTempTD) );
				this->enthalpyM = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthalpyF = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				if ( updatedTempTDT < this->TR && this->enthNew > this->enthalpyF ) {
					phaseChangeState = 1;
					this->enthNew = this->getEnthalpy(updatedTempTDT, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				} else if ( this->enthNew < this->enthalpyF && this->enthNew > this->enthalpyM && ( updatedTempTDT < prevTempTD|| updatedTempTDT > prevTempTD) ) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				} else if ( this->enthNew <= this->enthalpyM && updatedTempTDT >= prevTempTD&& this->enthNew > this->enthOld ) {
					phaseChangeState = -1;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				}
			} else if ( prevPhaseChangeState == 0 && phaseChangeState == 2 ) {
				this->enthRev = this->getEnthalpy(this->TR, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				this->enthalpyM = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthalpyF = this->getEnthalpy(updatedTempTDT, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				if ( this->enthNew < this->enthalpyF && this->enthNew > this->enthalpyM ) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				} else if ( this->enthNew <= this->enthalpyM && updatedTempTDT >= prevTempTD) {
					phaseChangeState = -1;
					this->enthNew = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				}
			} else if ( prevPhaseChangeState == -1 && phaseChangeState == 0 ) {
				this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthOld - ( this->specHeatTransition * prevTempTD) );
				this->enthalpyM = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthalpyF = this->getEnthalpy(updatedTempTDT, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				if ( this->enthNew < this->enthOld && updatedTempTDT < prevTempTD) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthOld - ( this->specHeatTransition * prevTempTD) );
				} else if ( this->enthNew < this->enthalpyF && this->enthNew > this->enthalpyM && updatedTempTDT < prevTempTD) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				} else if ( this->enthNew >= this->enthalpyF && updatedTempTDT <= this->TR ) {
					phaseChangeState = 0;
					this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
				}
			} else if ( prevPhaseChangeState == 0 && phaseChangeState == 1 ) {
				this->enthalpyM = this->getEnthalpy(updatedTempTDT, this->peakTempMelting, this->deltaTempMeltingLow, this->deltaTempMeltingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthalpyF = this->getEnthalpy(updatedTempTDT, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthRev = this->getEnthalpy(this->TR, this->peakTempFreezing, this->deltaTempFreezingLow, this->deltaTempFreezingHigh, DeltaH, this->specificHeatSolid, this->specificHeatLiquid );
				this->enthNew = ( this->specHeatTransition * updatedTempTDT ) + ( this->enthRev - ( this->specHeatTransition * this->TR ) );
			}
		}
		if ( this->phaseChangeTransition == 0 ) {
			if ( this->enthNew == this->enthOld ) {
				Cp = this->CpOld;
			} else {
				Cp = this->specHeat( prevTempTD, updatedTempTDT, Tc, Tau1, Tau2, DeltaH, this->specificHeatSolid, this->specificHeatLiquid, this->enthOld, this->enthNew );
			}
		} else if ( this->phaseChangeTransition == 1 ) {
			Cp = this->specHeatTransition;
		} else {
			Cp = this->CpOld;
		}
		this->CpOld = Cp;
		return Cp;
	}

	Real64 HysteresisPhaseChange::specHeat( Real64 temperaturePrev, Real64 temperatureCurrent, Real64 criticalTemperature, Real64 tau1, Real64 tau2, Real64 deltaH, Real64 CpSolid, Real64 CpLiquid, Real64 EnthalpyOld, Real64 EnthalpyNew ) {


//		REAL(r64), INTENT(IN)   ::  Tc                  ! Critical (Melting/Freezing) Temperature of PCM
//		REAL(r64), INTENT(IN)   ::  Tau1                ! Width of Melting Zone low
//		REAL(r64), INTENT(IN)   ::  Tau2                ! Width of Melting Zone high
//		REAL(r64), INTENT(IN)   ::  DeltaH              ! Latent Heat Stored in PCM During Phase Change
//		REAL(r64), INTENT(IN)   ::  CpSolid             ! Specific Heat of PCM in Solid State
//		REAL(r64), INTENT(IN)   ::  CpLiquid            ! Specific Heat of PCM in Liquid State
//		REAL(r64), INTENT(IN)   ::  EnthalpyOld         ! Previos Timestep Nodal Enthalpy
//		REAL(r64), INTENT(IN)   ::  EnthalpyNew         ! Current Timestep Nodal Enthalpy


        Real64 Cp1 = CpSolid;
        Real64 Cp2 = CpLiquid;
        Real64 T = temperatureCurrent;

        Real64 DEta1 = - ( deltaH * ( T - criticalTemperature ) * exp( -2 * abs( T - criticalTemperature ) / tau1 ) ) / ( tau1 * abs( T - criticalTemperature ) );
        Real64 DEta2 = ( deltaH * ( T - criticalTemperature ) * exp( -2 * abs( T - criticalTemperature ) / tau2 ) ) / ( tau2 * abs( T - criticalTemperature ) );

        if ( T < criticalTemperature ) {
            return (Cp1 + DEta1);
        } else if ( T == criticalTemperature ) {
            return (EnthalpyNew - EnthalpyOld) / (temperatureCurrent - temperaturePrev);
        } else if ( T > criticalTemperature ) {
            return Cp2 + DEta2;
        } else {
            return 0;
        }

	}

	void readAllHysteresisModels() {
		DataIPShortCuts::cCurrentModuleObject = "MaterialProperty:PhaseChangeHysteresis";
		numHysteresisModels = InputProcessor::GetNumObjectsFound( DataIPShortCuts::cCurrentModuleObject );

		if ( numHysteresisModels > 0 ) {
			for ( int hmNum = 1; hmNum <= numHysteresisModels; ++ hmNum ) {
				bool errorsFound = false;
				int ioStatus;
				int numAlphas;
				int numNumbers;
				InputProcessor::GetObjectItem(
					DataIPShortCuts::cCurrentModuleObject,
					hmNum,
					DataIPShortCuts::cAlphaArgs,
					numAlphas,
					DataIPShortCuts::rNumericArgs,
					numNumbers,
					ioStatus,
					DataIPShortCuts::lNumericFieldBlanks,
					DataIPShortCuts::lAlphaFieldBlanks,
					DataIPShortCuts::cAlphaFieldNames,
					DataIPShortCuts::cNumericFieldNames
				);
				bool isNotOK = false;
				bool isBlank = false;
				// InputProcessor::VerifyName( DataIPShortCuts::cAlphaArgs( 1 ),
				if ( isNotOK ) {
					errorsFound = true;
					if ( isBlank ) {
						DataIPShortCuts::cAlphaArgs( 1 ) = "xxxxx";
					}
				}
				HysteresisPhaseChange thisHM;
				thisHM.name = DataIPShortCuts::cAlphaArgs( 1 );
				thisHM.totalLatentHeat = DataIPShortCuts::rNumericArgs( 2 );
				thisHM.specificHeatLiquid = DataIPShortCuts::rNumericArgs( 3 );
				thisHM.deltaTempMeltingHigh = DataIPShortCuts::rNumericArgs( 4 );
				thisHM.peakTempMelting = DataIPShortCuts::rNumericArgs( 5 );
				thisHM.deltaTempMeltingLow = DataIPShortCuts::rNumericArgs( 6 );
				thisHM.specificHeatSolid = DataIPShortCuts::rNumericArgs( 7 );
				thisHM.deltaTempFreezingHigh = DataIPShortCuts::rNumericArgs( 8 );
				thisHM.peakTempFreezing = DataIPShortCuts::rNumericArgs( 9 );
				thisHM.deltaTempFreezingLow = DataIPShortCuts::rNumericArgs( 10 );
				thisHM.specHeatTransition = ( thisHM.specificHeatSolid + thisHM.specificHeatLiquid ) / 2.0;
				thisHM.CpOld = thisHM.specificHeatSolid;
				if ( errorsFound ) {
					ShowFatalError( "Error processing " + DataIPShortCuts::cCurrentModuleObject + " named \"" +DataIPShortCuts::cAlphaArgs( 1 ) + "\"" );
				}
				hysteresisPhaseChangeModels.push_back( thisHM );
				// SetupOutputVariable( "Phase Change State [ ]", thisHM.phaseChangeState, "Zone", "Average", thisHM.name );
			}
		}
	}

	void clear_state() {
		numHysteresisModels = 0;
		getHysteresisModels = true;
		hysteresisPhaseChangeModels.clear();
	}

} // namespace HysteresisPhaseChange

} // namespace EnergyPlus
