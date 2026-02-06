#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanThermalFunctions.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

using namespace URBANXX;

// Internal function to initialize temperatures
// Note: Parameters are pre-validated by caller (UrbanSetup)
// Exception handling is done by the caller
static void UrbanInitializeTemperature(UrbanType urban) {
  // Get references to temperature views
  auto &roofTemp = urban->roof.EffectiveSurfTemp;
  auto &imperviousRoadTemp = urban->imperviousRoad.EffectiveSurfTemp;
  auto &perviousRoadTemp = urban->perviousRoad.EffectiveSurfTemp;
  auto &sunlitWallTemp = urban->sunlitWall.EffectiveSurfTemp;
  auto &shadedWallTemp = urban->shadedWall.EffectiveSurfTemp;

  // Get references to layer temperature views
  auto &roofLayerTemp = urban->roof.Temperature;
  auto &imperviousRoadLayerTemp = urban->imperviousRoad.Temperature;
  auto &perviousRoadLayerTemp = urban->perviousRoad.Temperature;
  auto &sunlitWallLayerTemp = urban->sunlitWall.Temperature;
  auto &shadedWallLayerTemp = urban->shadedWall.Temperature;
  auto &buildingTemp = urban->building.Temperature;

  // Get references to canyon air properties
  auto &taf = urban->urbanCanyon.Taf;
  auto &qaf = urban->urbanCanyon.Qaf;

  const int numLandunits = urban->numLandunits;
  const int numUrbanLayers = urban->numUrbanLayers;
  const int numSoilLayers = urban->numSoilLayers;

  // Temperature initialization constants
  constexpr Real TEMP_ROOF_INIT = 292.0;
  constexpr Real TEMP_WALL_INIT = 292.0;
  constexpr Real TEMP_ROAD_INIT = 274.0;
  constexpr Real TEMP_CANYON_AIR_INIT = 283.0;
  constexpr Real QAF_INIT = 1.e-4; // kg/kg

  // Initialize surface temperatures and canyon air properties
  Kokkos::parallel_for(
      "InitializeSurfaceTemperatures", numLandunits, KOKKOS_LAMBDA(int l) {
        // Initialize effective surface temperatures
        roofTemp(l) = TEMP_ROOF_INIT;
        imperviousRoadTemp(l) = TEMP_ROAD_INIT;
        perviousRoadTemp(l) = TEMP_ROAD_INIT;
        sunlitWallTemp(l) = TEMP_WALL_INIT;
        shadedWallTemp(l) = TEMP_WALL_INIT;
        buildingTemp(l) = TEMP_WALL_INIT;

        // Initialize layer temperatures for urban surfaces (numUrbanLayers)
        for (int k = 0; k < numUrbanLayers; ++k) {
          roofLayerTemp(l, k) = TEMP_ROOF_INIT;
          sunlitWallLayerTemp(l, k) = TEMP_WALL_INIT;
          shadedWallLayerTemp(l, k) = TEMP_WALL_INIT;
        }

        // Initialize road layer temperatures (numSoilLayers)
        for (int k = 0; k < numSoilLayers; ++k) {
          imperviousRoadLayerTemp(l, k) = TEMP_ROAD_INIT;
          perviousRoadLayerTemp(l, k) = TEMP_ROAD_INIT;
        }

        // Initialize canyon air properties
        taf(l) = TEMP_CANYON_AIR_INIT; // Initialize to canyon air temperature
        qaf(l) = QAF_INIT; // Initialize to reasonable humidity value (kg/kg)
      });
  Kokkos::fence();
}

// Helper function to compute vertical discretization for a surface with uniform
// layers
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeVertDiscretizationForRoofOrWall(
    const Real thickness, const int numLevels, const int l, const ViewType &zc,
    const ViewType &dz, const ViewType &zi) {
  // Compute cell centers
  for (int k = 0; k < numLevels; ++k) {
    zc(l, k) = (k + 0.5) * (thickness / numLevels);
  }

  // Compute layer thickness
  dz(l, 0) = 0.5 * (zc(l, 0) + zc(l, 1));
  for (int k = 1; k < numLevels - 1; ++k) {
    dz(l, k) = 0.5 * (zc(l, k + 1) - zc(l, k - 1));
  }
  dz(l, numLevels - 1) = zc(l, numLevels - 1) - zc(l, numLevels - 2);

  // Compute interface depths
  zi(l, 0) = 0.0;
  for (int k = 1; k < numLevels; ++k) {
    zi(l, k) = 0.5 * (zc(l, k - 1) + zc(l, k));
  }
  zi(l, numLevels) = zc(l, numLevels - 1) + 0.5 * dz(l, numLevels - 1);
}

// Helper function to compute vertical discretization for roads with
// exponential layers
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void ComputeVertDiscretizationForRoad(
    const int numLevels, const int l, const ViewType &zc, const ViewType &dz,
    const ViewType &zi, const Real scalez, const Real zecoeff) {
  // Compute cell centers (node depths) with exponential spacing
  for (int k = 0; k < numLevels; ++k) {
    zc(l, k) = scalez * (Kokkos::exp(zecoeff * ((double)k + 0.5)) - 1.0);
  }

  // Compute layer thickness
  dz(l, 0) = 0.5 * (zc(l, 0) + zc(l, 1));
  for (int k = 1; k < numLevels - 1; ++k) {
    dz(l, k) = 0.5 * (zc(l, k + 1) - zc(l, k - 1));
  }
  dz(l, numLevels - 1) = zc(l, numLevels - 1) - zc(l, numLevels - 2);

  // Compute interface depths
  zi(l, 0) = 0.0;
  for (int k = 1; k < numLevels; ++k) {
    zi(l, k) = 0.5 * (zc(l, k - 1) + zc(l, k));
  }
  zi(l, numLevels) = zc(l, numLevels - 1) + 0.5 * dz(l, numLevels - 1);
}

static void UrbanInitializeVerticalDiscretization(UrbanType urban) {
  const int numLandunits = urban->numLandunits;
  const int numUrbanLayers = urban->numUrbanLayers;
  const int numSoilLayers = urban->numSoilLayers;

  // Access vertical discretization views
  auto &thick_wall = urban->urbanParams.building.WallThickness;
  auto &thick_roof = urban->urbanParams.building.RoofThickness;

  auto &zc_sunlit_wall = urban->sunlitWall.Zc;
  auto &dz_sunlit_wall = urban->sunlitWall.Dz;
  auto &zi_sunlit_wall = urban->sunlitWall.Zi;
  auto &depth_sunlit_wall = urban->sunlitWall.TotalDepth;

  auto &zc_shaded_wall = urban->shadedWall.Zc;
  auto &dz_shaded_wall = urban->shadedWall.Dz;
  auto &zi_shaded_wall = urban->shadedWall.Zi;
  auto &depth_shaded_wall = urban->shadedWall.TotalDepth;

  auto &zc_pervious_road = urban->perviousRoad.Zc;
  auto &dz_pervious_road = urban->perviousRoad.Dz;
  auto &zi_pervious_road = urban->perviousRoad.Zi;
  auto &depth_pervious_road = urban->perviousRoad.TotalDepth;

  auto &zc_impervious_road = urban->imperviousRoad.Zc;
  auto &dz_impervious_road = urban->imperviousRoad.Dz;
  auto &zi_impervious_road = urban->imperviousRoad.Zi;
  auto &depth_impervious_road = urban->imperviousRoad.TotalDepth;

  auto &zc_roof = urban->roof.Zc;
  auto &dz_roof = urban->roof.Dz;
  auto &zi_roof = urban->roof.Zi;
  auto &depth_roof = urban->roof.TotalDepth;

  // Initialize vertical discretization
  Kokkos::parallel_for(
      "UrbanInitializeVerticalDiscretization", numLandunits,
      KOKKOS_LAMBDA(int l) {
        // Road discretization parameters
        constexpr Real scalez = 0.025;
        constexpr Real zecoeff = 0.5;

        // Sunlit wall - uniform discretization
        ComputeVertDiscretizationForRoofOrWall(thick_wall(l), numUrbanLayers, l,
                                               zc_sunlit_wall, dz_sunlit_wall,
                                               zi_sunlit_wall);
        depth_sunlit_wall(l) = thick_wall(l);

        // Shaded wall - uniform discretization
        ComputeVertDiscretizationForRoofOrWall(thick_wall(l), numUrbanLayers, l,
                                               zc_shaded_wall, dz_shaded_wall,
                                               zi_shaded_wall);
        depth_shaded_wall(l) = thick_wall(l);

        // Roof - uniform discretization
        ComputeVertDiscretizationForRoofOrWall(thick_roof(l), numUrbanLayers, l,
                                               zc_roof, dz_roof, zi_roof);
        depth_roof(l) = thick_roof(l);

        // Pervious road - exponential discretization
        ComputeVertDiscretizationForRoad(numSoilLayers, l, zc_pervious_road,
                                         dz_pervious_road, zi_pervious_road,
                                         scalez, zecoeff);
        depth_pervious_road(l) = zi_pervious_road(l, numSoilLayers);

        // Impervious road - exponential discretization
        ComputeVertDiscretizationForRoad(numSoilLayers, l, zc_impervious_road,
                                         dz_impervious_road, zi_impervious_road,
                                         scalez, zecoeff);
        depth_impervious_road(l) = zi_impervious_road(l, numSoilLayers);
      });
  Kokkos::fence();
}

// Helper function to copy thermal properties from parameters to surface data
template <typename ViewType>
KOKKOS_INLINE_FUNCTION void
CopyThermalProperties(const int numLevels, const int l, const ViewType &src_tk,
                      const ViewType &src_cv, const ViewType &dst_tk,
                      const ViewType &dst_cv, const ViewType &dst_cv_times_dz,
                      const ViewType &dz) {
  for (int k = 0; k < numLevels; ++k) {
    dst_tk(l, k) = src_tk(l, k);
    dst_cv(l, k) = src_cv(l, k);
    dst_cv_times_dz(l, k) = src_cv(l, k) * dz(l, k);
  }
}

static void UrbanInitializeThermalProperties(UrbanType urban) {
  const int numLandunits = urban->numLandunits;
  const int numUrbanLayers = urban->numUrbanLayers;
  const int numSoilLayers = urban->numSoilLayers;

  // Access parameter thermal properties
  auto &tk_wall_params = urban->urbanParams.tk.Wall;
  auto &cv_wall_params = urban->urbanParams.cv.Wall;
  auto &tk_roof_params = urban->urbanParams.tk.Roof;
  auto &cv_roof_params = urban->urbanParams.cv.Roof;
  auto &tk_road_params = urban->urbanParams.tk.Road;
  auto &cv_road_params = urban->urbanParams.cv.Road;

  // Access surface thermal properties
  auto &tk_sunlit_wall = urban->sunlitWall.TkLayer;
  auto &cv_sunlit_wall = urban->sunlitWall.Cv;
  auto &cv_times_dz_sunlit_wall = urban->sunlitWall.CvTimesDz;
  auto &dz_sunlit_wall = urban->sunlitWall.Dz;
  auto &tk_shaded_wall = urban->shadedWall.TkLayer;
  auto &cv_shaded_wall = urban->shadedWall.Cv;
  auto &cv_times_dz_shaded_wall = urban->shadedWall.CvTimesDz;
  auto &dz_shaded_wall = urban->shadedWall.Dz;
  auto &tk_roof = urban->roof.TkLayer;
  auto &cv_roof = urban->roof.Cv;
  auto &cv_times_dz_roof = urban->roof.CvTimesDz;
  auto &dz_roof = urban->roof.Dz;
  auto &tk_pervious_road = urban->perviousRoad.TkLayer;
  auto &cv_pervious_road = urban->perviousRoad.Cv;
  auto &cv_times_dz_pervious_road = urban->perviousRoad.CvTimesDz;
  auto &dz_pervious_road = urban->perviousRoad.Dz;
  auto &tk_impervious_road = urban->imperviousRoad.TkLayer;
  auto &cv_impervious_road = urban->imperviousRoad.Cv;
  auto &cv_times_dz_impervious_road = urban->imperviousRoad.CvTimesDz;
  auto &dz_impervious_road = urban->imperviousRoad.Dz;

  // Copy thermal properties
  Kokkos::parallel_for(
      "UrbanInitializeThermalProperties", numLandunits, KOKKOS_LAMBDA(int l) {
        // Copy wall thermal properties to both sunlit and shaded walls
        CopyThermalProperties(numUrbanLayers, l, tk_wall_params, cv_wall_params,
                              tk_sunlit_wall, cv_sunlit_wall,
                              cv_times_dz_sunlit_wall, dz_sunlit_wall);
        CopyThermalProperties(numUrbanLayers, l, tk_wall_params, cv_wall_params,
                              tk_shaded_wall, cv_shaded_wall,
                              cv_times_dz_shaded_wall, dz_shaded_wall);

        // Copy roof thermal properties
        CopyThermalProperties(numUrbanLayers, l, tk_roof_params, cv_roof_params,
                              tk_roof, cv_roof, cv_times_dz_roof, dz_roof);

        // Copy road thermal properties to both pervious and impervious roads
        CopyThermalProperties(numSoilLayers, l, tk_road_params, cv_road_params,
                              tk_pervious_road, cv_pervious_road,
                              cv_times_dz_pervious_road, dz_pervious_road);
        CopyThermalProperties(numSoilLayers, l, tk_road_params, cv_road_params,
                              tk_impervious_road, cv_impervious_road,
                              cv_times_dz_impervious_road, dz_impervious_road);
      });
  Kokkos::fence();
}

// Helper function: Pedotransfer function based on Cosby et al. 1984, Table 5
// Converts sand and clay percentages to soil hydraulic properties
KOKKOS_INLINE_FUNCTION void PedotransferCosbyTable5(const Real sand,
                                                    const Real clay,
                                                    Real &watsat, Real &bsw,
                                                    Real &sucsat, Real &xksat) {
  // Cosby et al. 1984, Table 5
  // Input: sand, clay in percent (0-100)
  // Output: watsat (v/v), bsw (-), sucsat (mm), xksat (mm/s)
  watsat = 0.489 - 0.00126 * sand;
  bsw = 2.91 + 0.159 * clay;
  sucsat = 10.0 * Kokkos::pow(10.0, 1.88 - 0.0131 * sand);
  xksat = 0.0070556 * Kokkos::pow(10.0, -0.884 + 0.0153 * sand);
}

// Helper function: Depth-dependent organic matter properties
// Based on Letts et al. 2000 parameterization
KOKKOS_INLINE_FUNCTION void
OrganicPropertiesDepthDependent(const Real depth, const Real zsapric,
                                Real &om_watsat, Real &om_b, Real &om_sucsat,
                                Real &om_hksat) {
  // Depth-dependent properties for organic soil
  // depth: soil depth (m)
  // zsapric: depth at which organic matter becomes sapric peat (m)
  const Real depth_ratio = depth / zsapric;
  om_watsat = Kokkos::fmax(0.93 - 0.1 * depth_ratio, 0.83);
  om_b = Kokkos::fmin(2.7 + 9.3 * depth_ratio, 12.0);
  om_sucsat = Kokkos::fmin(10.3 - 0.2 * depth_ratio, 10.1);
  om_hksat = Kokkos::fmax(0.28 - 0.2799 * depth_ratio, 0.0001);
}

// Helper function: Calculate percolation fraction for organic soil
KOKKOS_INLINE_FUNCTION void CalculatePercolationFraction(const Real om_frac,
                                                         const Real pcalpha,
                                                         const Real pcbeta,
                                                         Real &perc_frac,
                                                         Real &uncon_frac) {
  // Percolation calculation for organic soil hydraulic conductivity
  if (om_frac > pcalpha) {
    const Real perc_norm = Kokkos::pow(1.0 - pcalpha, -pcbeta);
    perc_frac = perc_norm * Kokkos::pow(om_frac - pcalpha, pcbeta);
  } else {
    perc_frac = 0.0;
  }
  uncon_frac = (1.0 - om_frac) + (1.0 - perc_frac) * om_frac;
}

static void UrbanInitializePerviousRoadSoils(UrbanType urban) {
  const int numLandunits = urban->numLandunits;
  const int numSoilLayers = urban->numSoilLayers;
  const int numUrbanLayers = urban->numUrbanLayers;

  // Access soil property views for pervious road
  auto &sand = urban->perviousRoad.soil.Sand;
  auto &clay = urban->perviousRoad.soil.Clay;
  auto &organic = urban->perviousRoad.soil.Organic;
  auto &watsat = urban->perviousRoad.soil.WatSat;
  auto &tk_minerals = urban->perviousRoad.soil.TkMinerals;
  auto &tk_dry = urban->perviousRoad.soil.TkDry;
  auto &tk_saturated = urban->perviousRoad.soil.TkSaturated;
  auto &tkLayer = urban->perviousRoad.TkLayer;
  auto &cv_solids = urban->perviousRoad.soil.CvSolids;
  auto &water_liquid = urban->perviousRoad.soil.WaterLiquid;
  auto &water_ice = urban->perviousRoad.soil.WaterIce;
  auto &water_vol = urban->perviousRoad.soil.WaterVolumetric;
  auto &dz = urban->perviousRoad.Dz;
  auto &temp = urban->perviousRoad.Temperature;

  // Soil property constants from ELM (SoilStateType.F90)
  constexpr Real om_tkm = 0.25; // Thermal conductivity of organic soil [W/m-K]
  constexpr Real om_csol = 2.5e6; // Heat capacity of peat soil [J/(m³·K)]
  constexpr Real om_tkd =
      0.05; // Thermal conductivity of dry organic soil [W/m-K]
  constexpr Real zsapric = 0.5;  // Depth for sapric peat characteristics [m]
  constexpr Real pcalpha = 0.5;  // Percolation threshold
  constexpr Real pcbeta = 0.139; // Percolation exponent
  constexpr Real organic_max = 130.0;   // Maximum organic matter [kg/m³]
  constexpr Real soil_density = 2700.0; // Mineral soil density [kg/m³]

  // Initialize soil properties based on sand, clay, and organic matter inputs
  Kokkos::parallel_for(
      "UrbanInitializeSoilProperties", numLandunits, KOKKOS_LAMBDA(int l) {
        // Soil vertical discretization parameters (same as ELM)
        constexpr Real scalez = 0.025;
        constexpr Real zecoeff = 0.5;

        for (int k = 0; k < numSoilLayers; ++k) {
          // Get input soil properties (sand, clay in %, organic in kg/m³)
          const Real sand_pct = sand(l, k);
          const Real clay_pct = clay(l, k);

          // The soils within urban areas are assumed to have no organic matter.
          // Thus, organic(l, k) is not used.
          const Real organic_kgm3 = 0.0;

          // Compute soil depth using exponential discretization (same as ELM)
          // zsoi(j) = scalez * (exp(zecoeff * (j - 0.5)) - 1.0)
          const Real depth =
              scalez * (Kokkos::exp(zecoeff * ((double)k + 0.5)) - 1.0);

          // Calculate organic matter fraction (squared for non-lake soils)
          const Real om_frac_raw = organic_kgm3 / organic_max;
          const Real om_frac = Kokkos::fmin(om_frac_raw * om_frac_raw, 1.0);

          // Step 1: Pedotransfer function for mineral soil properties
          Real watsat_mineral, bsw_mineral, sucsat_mineral, xksat_mineral;
          PedotransferCosbyTable5(sand_pct, clay_pct, watsat_mineral,
                                  bsw_mineral, sucsat_mineral, xksat_mineral);

          // Step 2: Depth-dependent organic matter properties
          Real om_watsat, om_b, om_sucsat, om_hksat;
          OrganicPropertiesDepthDependent(depth, zsapric, om_watsat, om_b,
                                          om_sucsat, om_hksat);

          // Step 3: Calculate bulk density
          const Real bd = (1.0 - watsat_mineral) * soil_density;

          // Step 4: Mix mineral and organic properties
          const Real watsat_mixed =
              (1.0 - om_frac) * watsat_mineral + om_watsat * om_frac;

          // Thermal conductivity of soil minerals
          const Real tkm = (1.0 - om_frac) *
                               (8.80 * sand_pct + 2.92 * clay_pct) /
                               (sand_pct + clay_pct) +
                           om_tkm * om_frac;

          const Real bsw_mixed = (1.0 - om_frac) * bsw_mineral + om_frac * om_b;
          const Real sucsat_mixed =
              (1.0 - om_frac) * sucsat_mineral + om_sucsat * om_frac;

          // Step 5: Calculate percolation fraction
          Real perc_frac, uncon_frac;
          CalculatePercolationFraction(om_frac, pcalpha, pcbeta, perc_frac,
                                       uncon_frac);

          // Step 6: Calculate hydraulic conductivity with percolation
          Real uncon_hksat;
          if (om_frac < 1.0) {
            uncon_hksat =
                uncon_frac / ((1.0 - om_frac) / xksat_mineral +
                              ((1.0 - perc_frac) * om_frac) / om_hksat);
          } else {
            uncon_hksat = 0.0;
          }
          const Real hksat =
              uncon_frac * uncon_hksat + (perc_frac * om_frac) * om_hksat;

          // Step 7: Calculate thermal conductivities
          // Thermal conductivity of soil minerals (function of porosity)
          const Real tkmg = Kokkos::pow(tkm, 1.0 - watsat_mixed);

          // Saturated thermal conductivity
          const Real tksatu = tkmg * Kokkos::pow(0.57, watsat_mixed);

          // Dry thermal conductivity
          const Real tkdry =
              ((0.135 * bd + 64.7) / (soil_density - 0.947 * bd)) *
                  (1.0 - om_frac) +
              om_tkd * om_frac;

          // Step 8: Calculate heat capacity
          const Real csol =
              ((1.0 - om_frac) * (2.128 * sand_pct + 2.385 * clay_pct) /
                   (sand_pct + clay_pct) +
               om_csol * om_frac) *
              1.0e6; // Convert to J/(m³·K)

          // Store computed properties
          watsat(l, k) = watsat_mixed;
          tk_minerals(l, k) = tkmg;
          tk_dry(l, k) = tkdry;
          tk_saturated(l, k) = tksatu;
          tkLayer(l, k) = tkdry; // Initially set to dry thermal conductivity
          cv_solids(l, k) = csol;

          // Initialize volumetric water content
          // Layers 0-9 (first 10 layers): min(0.3, watsat)
          // Layers 10-14 (last 5 layers): 0.0 (hydrologically inactive)
          if (k < 10) {
            water_vol(l, k) = Kokkos::fmin(0.3, watsat_mixed);
          } else {
            water_vol(l, k) = 0.0;
            cv_solids(l, k) = CSOL_BEDROCK;
          }

          // Initialize liquid and ice water based on layer temperature
          constexpr Real tfrz = SHR_CONST_TKFRZ;
          constexpr Real denice = SHR_CONST_RHOICE;
          constexpr Real denh2o = SHR_CONST_RHOWATER;

          const Real layer_temp = temp(l, k);

          if (layer_temp <= tfrz) {
            water_ice(l, k) = dz(l, k) * denice * water_vol(l, k);
            water_liquid(l, k) = 0.0;
          } else {
            water_ice(l, k) = 0.0;
            water_liquid(l, k) = dz(l, k) * denh2o * water_vol(l, k);
          }
        }
      });
  Kokkos::fence();

  // Initialize hydrology boundary conditions and state variables
  auto &qflx_infl = urban->perviousRoad.QflxInfl;
  auto &qflx_tran = urban->perviousRoad.QflxTran;
  auto &zwt = urban->perviousRoad.Zwt;
  auto &h2osoi_liq = urban->perviousRoad.H2OSoiLiq;
  auto &h2osoi_ice = urban->perviousRoad.H2OSoiIce;
  auto &h2osoi_vol = urban->perviousRoad.H2OSoiVol;
  auto &qcharge = urban->perviousRoad.Qcharge;
  auto &qflx_deficit = urban->perviousRoad.QflxDeficit;

  Kokkos::parallel_for(
      "UrbanInitializeHydrologyBCs", numLandunits, KOKKOS_LAMBDA(int l) {
        // Initialize boundary conditions
        qflx_infl(l) = 0.0;          // No infiltration initially
        zwt(l) = 4.8018819123227204; // Water table depth [m]
        qcharge(l) = 0.0;            // No aquifer recharge initially
        qflx_deficit(l) = 0.0;       // No water deficit initially

        // Initialize transpiration to zero for all layers
        for (int k = 0; k < numSoilLayers; ++k) {
          qflx_tran(l, k) = 0.0;
        }

        // Copy water content from soil to hydrology state variables
        for (int k = 0; k < numSoilLayers; ++k) {
          h2osoi_liq(l, k) = water_liquid(l, k);
          h2osoi_ice(l, k) = water_ice(l, k);
          h2osoi_vol(l, k) = water_vol(l, k);
        }
      });
  Kokkos::fence();
}

static void UrbanInitializeInterfaceThermalProperties(UrbanType urban) {
  const int numLandunits = urban->numLandunits;
  const int numUrbanLayers = urban->numUrbanLayers;
  const int numSoilLayers = urban->numSoilLayers;

  // Access impervious road views
  auto &imperv_tkLayer = urban->imperviousRoad.TkLayer;
  auto &imperv_tkInterface = urban->imperviousRoad.TkInterface;
  auto &imperv_zc = urban->imperviousRoad.Zc;
  auto &imperv_zi = urban->imperviousRoad.Zi;
  auto &imperv_numActiveLayers = urban->imperviousRoad.NumberOfActiveLayers;

  // Access sunlit wall views
  auto &sunlit_tkLayer = urban->sunlitWall.TkLayer;
  auto &sunlit_tkInterface = urban->sunlitWall.TkInterface;
  auto &sunlit_zc = urban->sunlitWall.Zc;
  auto &sunlit_zi = urban->sunlitWall.Zi;

  // Access shaded wall views
  auto &shaded_tkLayer = urban->shadedWall.TkLayer;
  auto &shaded_tkInterface = urban->shadedWall.TkInterface;
  auto &shaded_zc = urban->shadedWall.Zc;
  auto &shaded_zi = urban->shadedWall.Zi;

  // Access roof views
  auto &roof_tkLayer = urban->roof.TkLayer;
  auto &roof_tkInterface = urban->roof.TkInterface;
  auto &roof_zc = urban->roof.Zc;
  auto &roof_zi = urban->roof.Zi;

  // Compute interface thermal conductivity for all surfaces
  Kokkos::parallel_for(
      "UrbanInitializeInterfaceThermalProperties", numLandunits,
      KOKKOS_LAMBDA(int l) {
        // Impervious road: variable active layers based on density class
        ComputeInterfaceThermalConductivity(
            l, numSoilLayers, imperv_numActiveLayers(l), imperv_tkLayer,
            imperv_tkInterface, imperv_zc, imperv_zi, 0.0);

        // Sunlit wall: all urban layers are active
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, numUrbanLayers, sunlit_tkLayer,
            sunlit_tkInterface, sunlit_zc, sunlit_zi,
            sunlit_tkLayer(l, numUrbanLayers - 1));

        // Shaded wall: all urban layers are active
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, numUrbanLayers, shaded_tkLayer,
            shaded_tkInterface, shaded_zc, shaded_zi,
            shaded_tkLayer(l, numUrbanLayers - 1));

        // Roof: all urban layers are active
        ComputeInterfaceThermalConductivity(
            l, numUrbanLayers, numUrbanLayers, roof_tkLayer, roof_tkInterface,
            roof_zc, roof_zi, roof_tkLayer(l, numUrbanLayers - 1));
      });
  Kokkos::fence();
}

extern "C" {
// Public setup function that performs all initialization steps
void UrbanSetup(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Initialize surface temperatures
    UrbanInitializeTemperature(urban);
    UrbanInitializeVerticalDiscretization(urban);
    UrbanInitializeThermalProperties(urban);
    UrbanInitializeInterfaceThermalProperties(urban);
    UrbanInitializePerviousRoadSoils(urban);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
