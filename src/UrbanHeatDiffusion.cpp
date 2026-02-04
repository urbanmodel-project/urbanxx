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

// Solve tridiagonal system using Thomas algorithm
// Solves: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = r[i]
// for i = 0 to n-1
// Note: a[0] and c[n-1] are not used
KOKKOS_INLINE_FUNCTION
void SolveTridiagonal(int n, const Real *a, const Real *b, const Real *c,
                      const Real *r, Real *x) {
  // Working arrays for modified coefficients
  Real cp[NUM_SOIL_LAYERS]; // Modified upper diagonal
  Real rp[NUM_SOIL_LAYERS]; // Modified right-hand side

  // Forward elimination
  cp[0] = c[0] / b[0];
  rp[0] = r[0] / b[0];

  for (int i = 1; i < n; i++) {
    const Real denom = b[i] - a[i] * cp[i - 1];
    cp[i] = c[i] / denom;
    rp[i] = (r[i] - a[i] * rp[i - 1]) / denom;
  }

  // Back substitution
  x[n - 1] = rp[n - 1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = rp[i] - cp[i] * x[i + 1];
  }
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
  auto &imperv_tkLayer = urban.imperviousRoad.TkLayer;
  auto &imperv_tkInterface = urban.imperviousRoad.TkInterface;
  auto &imperv_dz = urban.imperviousRoad.Dz;
  auto &imperv_zc = urban.imperviousRoad.Zc;
  auto &imperv_zi = urban.imperviousRoad.Zi;
  auto &imperv_temp = urban.imperviousRoad.Temperature;
  auto &imperv_cv_times_dz = urban.imperviousRoad.CvTimesDz;
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

        Real factPervRoad[NUM_SOIL_LAYERS];
        Real fnPervRoad[NUM_SOIL_LAYERS];
        Real aPervRoad[NUM_SOIL_LAYERS];
        Real bPervRoad[NUM_SOIL_LAYERS];
        Real cPervRoad[NUM_SOIL_LAYERS];
        Real rPervRoad[NUM_SOIL_LAYERS];
        Real newTempPervRoad[NUM_SOIL_LAYERS];

        Real factImpervRoad[NUM_SOIL_LAYERS];
        Real fnImpervRoad[NUM_SOIL_LAYERS];
        Real aImpervRoad[NUM_SOIL_LAYERS];
        Real bImpervRoad[NUM_SOIL_LAYERS];
        Real cImpervRoad[NUM_SOIL_LAYERS];
        Real rImpervRoad[NUM_SOIL_LAYERS];
        Real newTempImpervRoad[NUM_SOIL_LAYERS];

        const Real dtime = 30.0 * 60.0; // seconds
        int level;

        //
        // Compute factors and interface fluxes for pervious road
        //

        // Top layer
        level = 0;

        const Real capr =
            0.34; // Turing factor to turn first layer T into surface T
        const Real dz1 = perv_zc(l, level) - perv_zi(l, level);
        const Real dz2 = perv_zc(l, level + 1) - perv_zi(l, level);
        const Real dz_eff = 0.5 * (dz1 + capr * dz2);

        factPervRoad[level] =
            dtime / perv_cv_times_dz(l, level) * perv_dz(l, level) / dz_eff;

        fnPervRoad[level] = perv_tkInterface(l, level) *
                            (perv_temp(l, level + 1) - perv_temp(l, level)) /
                            (perv_zc(l, level + 1) - perv_zc(l, level));

        // Internal layers
        for (level = 1; level < numSoilLayers - 1; ++level) {
          factPervRoad[level] = dtime / perv_cv_times_dz(l, level);
          fnPervRoad[level] = perv_tkInterface(l, level) *
                              (perv_temp(l, level + 1) - perv_temp(l, level)) /
                              (perv_zc(l, level + 1) - perv_zc(l, level));
        }

        // Bottom layer
        level = numSoilLayers - 1;
        factPervRoad[level] = dtime / perv_cv_times_dz(l, level);
        fnPervRoad[level] = 0.0;

        //
        // Compute RHS vector for pervious road
        //

        // Top layer
        level = 0;
        Real dzp, dzm;

        rPervRoad[level] =
            perv_temp(l, level) +
            factPervRoad[level] *
                (perv_EflxGnet - perv_DEflxGnet_DTemp * perv_temp(l, level) +
                 CRANK_NICONSON_FACTOR * fnPervRoad[level]);

        dzp = perv_zc(l, level + 1) - perv_zc(l, level);
        aPervRoad[level] = 0.0;
        bPervRoad[level] = 1.0 +
                           (1.0 - CRANK_NICONSON_FACTOR) * factPervRoad[level] *
                               perv_tkInterface(l, level) / dzp -
                           factPervRoad[level] * perv_DEflxGnet_DTemp;
        cPervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                           factPervRoad[level] * perv_tkInterface(l, level) /
                           dzp;
        // Internal layers
        for (level = 1; level < numSoilLayers - 1; ++level) {
          rPervRoad[level] = perv_temp(l, level) +
                             factPervRoad[level] * CRANK_NICONSON_FACTOR *
                                 (fnPervRoad[level] - fnPervRoad[level - 1]);
          dzm = perv_zc(l, level) - perv_zc(l, level - 1);
          dzp = perv_zc(l, level + 1) - perv_zc(l, level);
          aPervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                             factPervRoad[level] *
                             perv_tkInterface(l, level - 1) / dzm;
          bPervRoad[level] = 1.0 + (1.0 - CRANK_NICONSON_FACTOR) *
                                       factPervRoad[level] *
                                       (perv_tkInterface(l, level) / dzp +
                                        perv_tkInterface(l, level - 1) / dzm);
          cPervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                             factPervRoad[level] * perv_tkInterface(l, level) /
                             dzp;
        }

        // Bottom layer
        level = numSoilLayers - 1;
        rPervRoad[level] = perv_temp(l, level) -
                           CRANK_NICONSON_FACTOR * factPervRoad[level - 1] *
                               fnPervRoad[level] +
                           factPervRoad[level] * fnPervRoad[level];
        dzm = perv_zc(l, level) - perv_zc(l, level - 1);
        aPervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                           factPervRoad[level] * perv_tkLayer(l, level - 1) /
                           dzm;
        bPervRoad[level] = 1.0 + (1.0 - CRANK_NICONSON_FACTOR) *
                                     factPervRoad[level] *
                                     perv_tkLayer(l, level - 1) / dzm;
        cPervRoad[level] = 0.0;

        // Step 5: Solve tridiagonal system for new temperatures
        SolveTridiagonal(numSoilLayers, aPervRoad, bPervRoad, cPervRoad,
                         rPervRoad, newTempPervRoad);

        // Step 6: Update pervious road temperature
        for (int j = 0; j < numSoilLayers; ++j) {
          perv_temp(l, j) = newTempPervRoad[j];
        }

        //
        // Solve heat diffusion for impervious road
        //

        //
        // Compute factors and interface fluxes for impervious road
        //

        // Top layer
        level = 0;

        const Real dz1_imperv = imperv_zc(l, level) - imperv_zi(l, level);
        const Real dz2_imperv = imperv_zc(l, level + 1) - imperv_zi(l, level);
        const Real dz_eff_imperv = 0.5 * (dz1_imperv + capr * dz2_imperv);

        factImpervRoad[level] = dtime / imperv_cv_times_dz(l, level) *
                                imperv_dz(l, level) / dz_eff_imperv;

        fnImpervRoad[level] =
            imperv_tkInterface(l, level) *
            (imperv_temp(l, level + 1) - imperv_temp(l, level)) /
            (imperv_zc(l, level + 1) - imperv_zc(l, level));

        // Internal layers
        for (level = 1; level < numSoilLayers - 1; ++level) {
          factImpervRoad[level] = dtime / imperv_cv_times_dz(l, level);
          fnImpervRoad[level] =
              imperv_tkInterface(l, level) *
              (imperv_temp(l, level + 1) - imperv_temp(l, level)) /
              (imperv_zc(l, level + 1) - imperv_zc(l, level));
        }

        // Bottom layer
        level = numSoilLayers - 1;
        factImpervRoad[level] = dtime / imperv_cv_times_dz(l, level);
        fnImpervRoad[level] = 0.0;

        //
        // Compute RHS vector for impervious road
        //

        // Top layer
        level = 0;

        rImpervRoad[level] =
            imperv_temp(l, level) +
            factImpervRoad[level] *
                (imperv_EflxGnet -
                 imperv_DEflxGnet_DTemp * imperv_temp(l, level) +
                 CRANK_NICONSON_FACTOR * fnImpervRoad[level]);

        dzp = imperv_zc(l, level + 1) - imperv_zc(l, level);
        aImpervRoad[level] = 0.0;
        bImpervRoad[level] = 1.0 +
                             (1.0 - CRANK_NICONSON_FACTOR) *
                                 factImpervRoad[level] *
                                 imperv_tkInterface(l, level) / dzp -
                             factImpervRoad[level] * imperv_DEflxGnet_DTemp;
        cImpervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                             factImpervRoad[level] *
                             imperv_tkInterface(l, level) / dzp;

        // Internal layers
        for (level = 1; level < numSoilLayers - 1; ++level) {
          rImpervRoad[level] =
              imperv_temp(l, level) +
              factImpervRoad[level] * CRANK_NICONSON_FACTOR *
                  (fnImpervRoad[level] - fnImpervRoad[level - 1]);
          dzm = imperv_zc(l, level) - imperv_zc(l, level - 1);
          dzp = imperv_zc(l, level + 1) - imperv_zc(l, level);
          aImpervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                               factImpervRoad[level] *
                               imperv_tkInterface(l, level - 1) / dzm;
          bImpervRoad[level] =
              1.0 + (1.0 - CRANK_NICONSON_FACTOR) * factImpervRoad[level] *
                        (imperv_tkInterface(l, level) / dzp +
                         imperv_tkInterface(l, level - 1) / dzm);
          cImpervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                               factImpervRoad[level] *
                               imperv_tkInterface(l, level) / dzp;
        }

        // Bottom layer
        level = numSoilLayers - 1;
        rImpervRoad[level] = imperv_temp(l, level) -
                             CRANK_NICONSON_FACTOR * factImpervRoad[level - 1] *
                                 fnImpervRoad[level] +
                             factImpervRoad[level] * fnImpervRoad[level];
        dzm = imperv_zc(l, level) - imperv_zc(l, level - 1);
        aImpervRoad[level] = -(1.0 - CRANK_NICONSON_FACTOR) *
                             factImpervRoad[level] *
                             imperv_tkInterface(l, level - 1) / dzm;
        bImpervRoad[level] = 1.0 + (1.0 - CRANK_NICONSON_FACTOR) *
                                       factImpervRoad[level] *
                                       imperv_tkInterface(l, level - 1) / dzm;
        cImpervRoad[level] = 0.0;

        // Solve tridiagonal system for new temperatures
        if (l == 0) {
          printf("\nImpervious Road Tridiagonal System:\n");
          for (int j = 0; j < numSoilLayers; ++j) {
            printf("Layer %2d: a=%18.16f, b=%18.16f, c=%18.16f, r=%18.16f\n", j,
                   aImpervRoad[j], bImpervRoad[j], cImpervRoad[j],
                   rImpervRoad[j]);
          }
        }

        SolveTridiagonal(numSoilLayers, aImpervRoad, bImpervRoad, cImpervRoad,
                         rImpervRoad, newTempImpervRoad);

        // Update impervious road temperature
        for (int j = 0; j < numSoilLayers; ++j) {
          imperv_temp(l, j) = newTempImpervRoad[j];
        }

        // TODO: Add remaining heat diffusion:
        // - Roof (numUrbanLayers)
        // - Sunlit wall (numUrbanLayers)
        // - Shaded wall (numUrbanLayers)
        // - Update surface temperatures for all surfaces
      });
  Kokkos::fence();
}

} // namespace URBANXX
