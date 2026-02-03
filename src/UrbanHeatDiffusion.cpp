// Urban Heat Diffusion Implementation
// 1D heat diffusion dynamics for Urban Surfaces
// Based on ELM SoilTemperature.F90

#include "private/UrbanHeatDiffusionImpl.h"
#include "private/UrbanThermalFunctions.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {

// Compute ground net energy flux for a surface
KOKKOS_INLINE_FUNCTION
Real ComputeGroundNetEnergyFlux(Real netSw, Real netLw, Real eflxShGrnd,
                                Real qflxEvapSoil = 0.0,
                                Real qflxTranEvap = 0.0) {
  return netSw - netLw -
         (eflxShGrnd + qflxEvapSoil * SHR_CONST_LATVAP +
          qflxTranEvap * SHR_CONST_LATVAP);
}

// Compute derivative of ground net energy flux w.r.t. temperature
KOKKOS_INLINE_FUNCTION
Real ComputeGroundNetEnergyFluxDerivative(Real cgrnds, Real cgrndl, Real emiss,
                                          Real temp) {
  const Real cgrnd = cgrnds + cgrndl * SHR_CONST_LATVAP;
  const Real dlwrdDTemp = 4.0 * emiss * STEBOL * Kokkos::pow(temp, 3.0);
  return -cgrnd - dlwrdDTemp;
}

// Compute 1D heat diffusion for all urban surfaces
void ComputeHeatDiffusion(URBANXX::_p_UrbanType &urban) {

  const int numLandunits = urban.numLandunits;
  const int numSoilLayers = urban.numSoilLayers;
  const int numUrbanLayers = urban.numUrbanLayers;

  // Access pervious road soil property views
  auto &perv_tk_minerals = urban.perviousRoad.soil.TkMinerals;
  auto &perv_tk_dry = urban.perviousRoad.soil.TkDry;
  auto &perv_tkLayer = urban.perviousRoad.TkLayer;
  auto &perv_tkInterface = urban.perviousRoad.TkInterface;
  auto &perv_watsat = urban.perviousRoad.soil.WatSat;
  auto &perv_water_liquid = urban.perviousRoad.soil.WaterLiquid;
  auto &perv_water_ice = urban.perviousRoad.soil.WaterIce;
  auto &perv_dz = urban.perviousRoad.Dz;
  auto &perv_zc = urban.perviousRoad.Zc;
  auto &perv_zi = urban.perviousRoad.Zi;
  auto &perv_temp = urban.perviousRoad.Temperature;
  auto &perv_cv_solids = urban.perviousRoad.soil.CvSolids;
  auto &perv_cv_times_dz = urban.perviousRoad.CvTimesDz;
  auto &perv_netLw = urban.perviousRoad.NetLongRad;
  auto &perv_netSw = urban.perviousRoad.NetShortRad;
  auto &perv_EflxShGrnd = urban.perviousRoad.EflxShGrnd;
  auto &perv_QflxEvapSoil = urban.perviousRoad.QflxEvapSoil;
  auto &perv_QflxTranEvap = urban.perviousRoad.QflxTranEvap;
  auto &perv_Cgrndl = urban.perviousRoad.Cgrndl;
  auto &perv_Cgrnds = urban.perviousRoad.Cgrnds;
  auto &perv_Temp = urban.perviousRoad.EffectiveSurfTemp;
  auto &perv_emiss = urban.urbanParams.emissivity.PerviousRoad;

  // Access roof property views
  auto &roof_netLw = urban.roof.NetLongRad;
  auto &roof_netSw = urban.roof.NetShortRad;
  auto &roof_EflxShGrnd = urban.roof.EflxShGrnd;
  auto &roof_QflxEvapSoil = urban.roof.QflxEvapSoil;
  auto &roof_QflxTranEvap = urban.roof.QflxTranEvap;
  auto &roof_Cgrnds = urban.roof.Cgrnds;
  auto &roof_Cgrndl = urban.roof.Cgrndl;
  auto &roof_Temp = urban.roof.EffectiveSurfTemp;
  auto &roof_emiss = urban.urbanParams.emissivity.Roof;

  // Access impervious road property views
  auto &imperv_netLw = urban.imperviousRoad.NetLongRad;
  auto &imperv_netSw = urban.imperviousRoad.NetShortRad;
  auto &imperv_EflxShGrnd = urban.imperviousRoad.EflxShGrnd;
  auto &imperv_QflxEvapSoil = urban.imperviousRoad.QflxEvapSoil;
  auto &imperv_QflxTranEvap = urban.imperviousRoad.QflxTranEvap;
  auto &imperv_Cgrnds = urban.imperviousRoad.Cgrnds;
  auto &imperv_Cgrndl = urban.imperviousRoad.Cgrndl;
  auto &imperv_Temp = urban.imperviousRoad.EffectiveSurfTemp;
  auto &imperv_emiss = urban.urbanParams.emissivity.ImperviousRoad;

  // Access sunlit wall property views
  auto &sunwall_netLw = urban.sunlitWall.NetLongRad;
  auto &sunwall_netSw = urban.sunlitWall.NetShortRad;
  auto &sunwall_EflxShGrnd = urban.sunlitWall.EflxShGrnd;
  auto &sunwall_Cgrnds = urban.sunlitWall.Cgrnds;
  auto &sunwall_Cgrndl = urban.sunlitWall.Cgrndl;
  auto &sunwall_Temp = urban.sunlitWall.EffectiveSurfTemp;
  auto &wall_emiss = urban.urbanParams.emissivity.Wall;

  // Access shaded wall property views
  auto &shadewall_netLw = urban.shadedWall.NetLongRad;
  auto &shadewall_netSw = urban.shadedWall.NetShortRad;
  auto &shadewall_EflxShGrnd = urban.shadedWall.EflxShGrnd;
  auto &shadewall_Cgrnds = urban.shadedWall.Cgrnds;
  auto &shadewall_Cgrndl = urban.shadedWall.Cgrndl;
  auto &shadewall_Temp = urban.shadedWall.EffectiveSurfTemp;

  // Single parallel kernel over all landunits
  Kokkos::parallel_for(
      "ComputeHeatDiffusion", numLandunits, KOKKOS_LAMBDA(int l) {
        // Step 1: Compute thermal conductivity for pervious road soil layers
        ComputeSoilThermalConductivity(
            l, numSoilLayers, perv_tk_minerals, perv_tk_dry, perv_tkLayer,
            perv_watsat, perv_water_liquid, perv_water_ice, perv_dz, perv_temp);

        // Step 2: Compute thermal conductivity at layer interfaces
        // (Note: Interface thermal conductivity for impervious road, walls, and
        // roof
        //  are computed during initialization and don't need to be recomputed)

        // Pervious road: all soil layers are active
        ComputeInterfaceThermalConductivity(l, numSoilLayers, numSoilLayers,
                                            perv_tkLayer, perv_tkInterface,
                                            perv_zc, perv_zi, 0.0);

        // Step 3: Compute heat capacity times layer thickness for pervious road
        ComputeSoilHeatCapacityTimesDz(
            l, numSoilLayers, perv_cv_solids, perv_watsat, perv_water_liquid,
            perv_water_ice, perv_dz, perv_cv_times_dz);

        // Step 4: Compute ground net energy flux and its derivative for all
        // surfaces

        // Pervious road (with evaporation)
        const Real perv_EflxGnet = ComputeGroundNetEnergyFlux(
            perv_netSw(l), perv_netLw(l), perv_EflxShGrnd(l),
            perv_QflxEvapSoil(l), perv_QflxTranEvap(l));
        const Real perv_DEflxGnet_DTemp = ComputeGroundNetEnergyFluxDerivative(
            perv_Cgrnds(l), perv_Cgrndl(l), perv_emiss(l), perv_Temp(l));

        // Roof (with evaporation)
        const Real roof_EflxGnet = ComputeGroundNetEnergyFlux(
            roof_netSw(l), roof_netLw(l), roof_EflxShGrnd(l),
            roof_QflxEvapSoil(l), roof_QflxTranEvap(l));
        const Real roof_DEflxGnet_DTemp = ComputeGroundNetEnergyFluxDerivative(
            roof_Cgrnds(l), roof_Cgrndl(l), roof_emiss(l), roof_Temp(l));

        // Impervious road (with evaporation)
        const Real imperv_EflxGnet = ComputeGroundNetEnergyFlux(
            imperv_netSw(l), imperv_netLw(l), imperv_EflxShGrnd(l),
            imperv_QflxEvapSoil(l), imperv_QflxTranEvap(l));
        const Real imperv_DEflxGnet_DTemp =
            ComputeGroundNetEnergyFluxDerivative(
                imperv_Cgrnds(l), imperv_Cgrndl(l), imperv_emiss(l),
                imperv_Temp(l));

        // Sunlit wall (no evaporation - defaults to 0)
        const Real sunwall_EflxGnet = ComputeGroundNetEnergyFlux(
            sunwall_netSw(l), sunwall_netLw(l), sunwall_EflxShGrnd(l));
        const Real sunwall_DEflxGnet_DTemp =
            ComputeGroundNetEnergyFluxDerivative(
                sunwall_Cgrnds(l), sunwall_Cgrndl(l), wall_emiss(l),
                sunwall_Temp(l));

        // Shaded wall (no evaporation - defaults to 0)
        const Real shadewall_EflxGnet = ComputeGroundNetEnergyFlux(
            shadewall_netSw(l), shadewall_netLw(l), shadewall_EflxShGrnd(l));
        const Real shadewall_DEflxGnet_DTemp =
            ComputeGroundNetEnergyFluxDerivative(
                shadewall_Cgrnds(l), shadewall_Cgrndl(l), wall_emiss(l),
                shadewall_Temp(l));

        // TODO: Add remaining heat diffusion steps:
        // Step 5: Setup tridiagonal system for heat conduction equation
        // Step 6: Apply boundary conditions (surface flux, bottom temperature)
        // Step 7: Solve tridiagonal system for new temperatures
        // Step 8: Update surface temperatures
      });
  Kokkos::fence();
}

} // namespace URBANXX
