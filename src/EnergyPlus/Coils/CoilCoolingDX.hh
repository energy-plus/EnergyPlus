#ifndef ENERGYPLUS_COILS_COILCOOLINGDX
#define ENERGYPLUS_COILS_COILCOOLINGDX

#include <string>
#include <vector>

#include <Coils/CoilCoolingDXCurveFitPerformance.hh>
#include <EnergyPlus.hh>

namespace EnergyPlus {

class CoilCoolingDXInputSpecification
{

public:
    std::string name;
    std::string evaporator_inlet_node_name;
    std::string evaporator_outlet_node_name;
    std::string availability_schedule_name;
    std::string condenser_zone_name;
    std::string condenser_inlet_node_name;
    std::string condenser_outlet_node_name;
    std::string performance_object_name;
    std::string condensate_collection_water_storage_tank_name;
    std::string evaporative_condenser_supply_water_storage_tank_name;
};

class CoilCoolingDX
{

    std::string const object_name = "Coil:Cooling:DX";
    CoilCoolingDXInputSpecification original_input_specs;

public:
    CoilCoolingDX(std::string name);

    void instantiateFromInputSpec(CoilCoolingDXInputSpecification input_data);
    void onetimeinit();
    void simulate(bool useAlternateMode, Real64 PLR, int speedNum, Real64 speedRatio, int fanOpMode);

    std::string name;
    bool myOneTimeInitFlag;
    int evapInletNodeIndex;
    int evapOutletNodeIndex;
    int availScheduleIndex;
    int condZoneIndex;
    int condInletNodeIndex;
    int condOutletNodeIndex;
    CoilCoolingDXCurveFitPerformance performance; // TODO: Change to unique pointer when we have base class for performance object
    int condensateTankIndex;
    int evaporativeCondSupplyTankIndex;
    int airLoopNum; // Add for AFN compatibility, revisit at a later date

    // report variables
    Real64 totalCoolingEnergyRate;
    Real64 totalCoolingEnergy;
    Real64 sensCoolingEnergyRate;
    Real64 sensCoolingEnergy;
    Real64 latCoolingEnergyRate;
    Real64 latCoolingEnergy;
    Real64 elecCoolingPower;
    Real64 elecCoolingConsumption;
    Real64 coolingCoilRuntimeFraction;

    int getDXCoilCapFTCurveIndex();

    Real64 getRatedGrossTotalCapacity();

    void inline passThroughNodeData(DataLoopNode::NodeData &in, DataLoopNode::NodeData &out);
};

extern std::vector<CoilCoolingDX> coilCoolingDXs;

} // namespace EnergyPlus

#endif // ENERGYPLUS_COILS_COILCOOLINGDX