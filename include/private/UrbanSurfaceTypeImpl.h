#ifndef URBAN_SURFACE_TYPE_IMPL_H
#define URBAN_SURFACE_TYPE_IMPL_H
#include "private/DataTypesImpl.h"
#include "private/KokkosViewMacros.h"
#include <Kokkos_Core.hpp>

namespace URBANXX {
// Composite data for combined pervious and impervious road surfaces
struct CompositeRoadSurfaceData {
  DECLARE_DEVICE_VIEW(3DR8, DownwellingShortRad) // downwelling shortwave
                                                 // radiation per unit road area
  DECLARE_DEVICE_VIEW(3DR8, SnowAlbedo)          // snow albedo
  DECLARE_DEVICE_VIEW(
      3DR8, AlbedoWithSnowEffects) // albedo of road including snow effects

  CompositeRoadSurfaceData(int numLandunits, int numRadBands, int numRadTypes) {
    ALLOCATE_DEVICE_VIEW(DownwellingShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(SnowAlbedo, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(AlbedoWithSnowEffects, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

struct BuildingDataType {
  DECLARE_DEVICE_VIEW(1DR8, Temperature;)    // building temperature (K)
  DECLARE_DEVICE_VIEW(1DR8, EFlxForHeating;) // building heat flux (W/m^2)
  DECLARE_DEVICE_VIEW(1DR8, EFluxForAC;)     // building cool flux (W/m^2)

  BuildingDataType(int numLandunits) {
    ALLOCATE_DEVICE_VIEW(Temperature, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EFlxForHeating, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EFluxForAC, Array1DR8, numLandunits)
  }
};

struct SoilDataType {
  DECLARE_DEVICE_VIEW(2DR8,
                      TkMinerals)  // thermal conductivity, soil
                                   // minerals [W/m-K]
  DECLARE_DEVICE_VIEW(2DR8, TkDry) // thermal conductivity, dry soil [W/m-K]
  DECLARE_DEVICE_VIEW(2DR8,
                      TkSaturated)    // thermal conductivity,
                                      // saturated soil [W/m-K]
  DECLARE_DEVICE_VIEW(2DR8, CvSolids) // heat capacity, soil solids [J/m^3/K]
  DECLARE_DEVICE_VIEW(2DR8,
                      WatSat)                // volumetric soil water at
                                             // saturation (porosity) [-]
  DECLARE_DEVICE_VIEW(2DR8, WaterLiquid)     // liquid water [kg/m^2]
  DECLARE_DEVICE_VIEW(2DR8, WaterIce)        // ice lens [kg/m^2]
  DECLARE_DEVICE_VIEW(2DR8, WaterVolumetric) // volumetric water content [-]
  DECLARE_DEVICE_VIEW(2DR8, Sand)            // soil texture: percent sand [-]
  DECLARE_DEVICE_VIEW(2DR8, Clay)            // soil texture: percent clay [-]
  DECLARE_DEVICE_VIEW(2DR8, Organic)         // organic matter [kg/m³]

  // Hydraulic properties
  DECLARE_DEVICE_VIEW(2DR8, HkSat)  // saturated hydraulic conductivity [mm/s]
  DECLARE_DEVICE_VIEW(2DR8, Bsw)    // Clapp and Hornberger "b" parameter [-]
  DECLARE_DEVICE_VIEW(2DR8, SucSat) // saturated suction [mm]

  // Tridiagonal matrix arrays for hydrology solver (nlevbed+1 to include
  // aquifer layer)
  DECLARE_DEVICE_VIEW(2DR8, Amx)  // left off-diagonal of tridiagonal matrix
  DECLARE_DEVICE_VIEW(2DR8, Bmx)  // diagonal of tridiagonal matrix
  DECLARE_DEVICE_VIEW(2DR8, Cmx)  // right off-diagonal of tridiagonal matrix
  DECLARE_DEVICE_VIEW(2DR8, Rmx)  // right-hand side forcing vector
  DECLARE_DEVICE_VIEW(2DR8, Dwat) // solution vector (change in water content)

  SoilDataType(int numLandunits, int numSoilLayers) {
    ALLOCATE_DEVICE_VIEW(TkMinerals, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(TkDry, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(CvSolids, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(WatSat, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(TkSaturated, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(WaterLiquid, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(WaterIce, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(WaterVolumetric, Array2DR8, numLandunits,
                         numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Sand, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Clay, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Organic, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(HkSat, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Bsw, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(SucSat, Array2DR8, numLandunits, numSoilLayers)

    // Allocate tridiagonal arrays (+1 for aquifer layer)
    ALLOCATE_DEVICE_VIEW(Amx, Array2DR8, numLandunits, numSoilLayers + 1)
    ALLOCATE_DEVICE_VIEW(Bmx, Array2DR8, numLandunits, numSoilLayers + 1)
    ALLOCATE_DEVICE_VIEW(Cmx, Array2DR8, numLandunits, numSoilLayers + 1)
    ALLOCATE_DEVICE_VIEW(Rmx, Array2DR8, numLandunits, numSoilLayers + 1)
    ALLOCATE_DEVICE_VIEW(Dwat, Array2DR8, numLandunits, numSoilLayers + 1)
  }
};

//
// struct hierarchy:
//
// SurfaceDataBase (base)
// ├── SnowCoveredSurfaceData
// │   ├── ImperviousRoadDataType
// │   ├── PerviousRoadDataType
// │   └── RoofDataType
// └── WallDataType
//

// Base struct for all surface types with common radiative properties
struct SurfaceDataBase {
  DECLARE_DEVICE_VIEW(3DR8,
                      ReflectedShortRad) // reflected solar radiation per unit
                                         // surface area per unit incident flux
  DECLARE_DEVICE_VIEW(3DR8,
                      AbsorbedShortRad) // absorbed solar radiation per unit
                                        // surface area per unit incident flux
  DECLARE_DEVICE_VIEW(1DR8, Emissivity) // emissivity
  DECLARE_DEVICE_VIEW(1DR8, EffectiveSurfTemp) // effective temperature
  DECLARE_DEVICE_VIEW(1DR8, NetLongRad)        // net longwave radiation
  DECLARE_DEVICE_VIEW(1DR8, UpwardLongRad)     // upward longwave radiation
  DECLARE_DEVICE_VIEW(1DR8, NetShortRad)       // net shortwave radiation
  DECLARE_DEVICE_VIEW(1DR8,
                      Cgrnds) // d(sensible heat flux)/dT for implicit solver
  DECLARE_DEVICE_VIEW(1DR8,
                      Cgrndl) // d(latent heat flux)/dT for implicit solver
  DECLARE_DEVICE_VIEW(1DR8, EflxShGrnd) // sensible heat flux from ground
                                        // (W/m**2) [+ to atm]

  // Saturation humidity variables
  DECLARE_DEVICE_VIEW(1DR8, Es)   // saturation vapor pressure (Pa)
  DECLARE_DEVICE_VIEW(1DR8, EsdT) // d(es)/dT (Pa/K)
  DECLARE_DEVICE_VIEW(1DR8, Qs)   // saturation specific humidity (kg/kg)
  DECLARE_DEVICE_VIEW(1DR8, QsdT) // d(qs)/dT (1/K)

  // Information about layers
  DECLARE_DEVICE_VIEW(2DR8, Zc)         // depth at the center of layer (m)
  DECLARE_DEVICE_VIEW(2DR8, Zi)         // depth at the layer interface (m)
  DECLARE_DEVICE_VIEW(2DR8, Dz)         // layer thickness (m)
  DECLARE_DEVICE_VIEW(1DR8, TotalDepth) // total depth of surface layers (m)

  // Thermal properties
  DECLARE_DEVICE_VIEW(2DR8, TkLayer) // thermal conductivity (W/m/K)
  DECLARE_DEVICE_VIEW(
      2DR8, TkInterface) // thermal conductivity at layer interface (W/m/K)
  DECLARE_DEVICE_VIEW(2DR8, Cv) // volumetric heat capacity (J/m^3/K)
  DECLARE_DEVICE_VIEW(
      2DR8, CvTimesDz) // heat capacity times layer thickness (J/m^2/K)
  DECLARE_DEVICE_VIEW(2DR8, Temperature) // layer temperature (K)

  SurfaceDataBase(int numLandunits, int numRadBands, int numRadTypes,
                  int numLayers) {
    ALLOCATE_DEVICE_VIEW(ReflectedShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(AbsorbedShortRad, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(Emissivity, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EffectiveSurfTemp, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(NetLongRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(UpwardLongRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(NetShortRad, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Cgrnds, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Cgrndl, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EflxShGrnd, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Es, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(EsdT, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Qs, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(QsdT, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Zc, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(Zi, Array2DR8, numLandunits, numLayers + 1)
    ALLOCATE_DEVICE_VIEW(Dz, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(TotalDepth, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(TkLayer, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(TkInterface, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(Cv, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(CvTimesDz, Array2DR8, numLandunits, numLayers)
    ALLOCATE_DEVICE_VIEW(Temperature, Array2DR8, numLandunits, numLayers)
  }
};

// For surfaces with snow coverage (roads and roofs)
struct SnowCoveredSurfaceData : SurfaceDataBase {
  DECLARE_DEVICE_VIEW(3DR8, SnowAlbedo) // snow albedo
  DECLARE_DEVICE_VIEW(3DR8,
                      AlbedoWithSnowEffects) // albedo including snow effects
  DECLARE_DEVICE_VIEW(1DR8,
                      QflxEvapSoil) // soil evaporation (mm H2O/s) (+ = to atm)
  DECLARE_DEVICE_VIEW(1DR8, QflxTranEvap) // vegetation evaporation (mm H2O/s)
                                          // (+ = to atm)

  SnowCoveredSurfaceData(int numLandunits, int numRadBands, int numRadTypes,
                         int numLayers)
      : SurfaceDataBase(numLandunits, numRadBands, numRadTypes, numLayers) {
    ALLOCATE_DEVICE_VIEW(SnowAlbedo, Array3DR8, numLandunits, numRadBands,
                         numRadTypes)
    ALLOCATE_DEVICE_VIEW(AlbedoWithSnowEffects, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
    ALLOCATE_DEVICE_VIEW(QflxEvapSoil, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(QflxTranEvap, Array1DR8, numLandunits)
  }
};

struct ImperviousRoadDataType : SnowCoveredSurfaceData {
  DECLARE_DEVICE_VIEW(1DI4, NumberOfActiveLayers) // number of active layers
  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase
  ImperviousRoadDataType(int numLandunits, int numRadBands, int numRadTypes,
                         int numLayers)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes,
                               numLayers) {
    // Don't initialize to 0 - must be set via
    // UrbanSetNumberOfActiveLayersImperviousRoad
    ALLOCATE_VIEW_NO_INIT(NumberOfActiveLayers, Array1DI4, numLandunits)
  }
};

struct PerviousRoadDataType : SnowCoveredSurfaceData {
  SoilDataType soil;

  // Hydrology state variables (layer-level)
  DECLARE_DEVICE_VIEW(2DR8, H2OSoiLiq) // liquid water content [kg/m²]
  DECLARE_DEVICE_VIEW(2DR8, H2OSoiIce) // ice water content [kg/m²]
  DECLARE_DEVICE_VIEW(2DR8, H2OSoiVol) // volumetric water content [m³/m³]
  DECLARE_DEVICE_VIEW(2DR8,
                      EffPorosity) // effective porosity (porosity - ice) [-]
  DECLARE_DEVICE_VIEW(2DR8, Hk)    // hydraulic conductivity [mm/s]
  DECLARE_DEVICE_VIEW(2DR8, Smp)   // soil matrix potential [mm]

  // Hydrology flux variables (layer-level)
  DECLARE_DEVICE_VIEW(1DR8, QflxInfl) // infiltration flux [mm/s]
  DECLARE_DEVICE_VIEW(2DR8, QflxTran) // transpiration flux [mm/s]

  // Column-level hydrology variables
  DECLARE_DEVICE_VIEW(1DR8, Zwt)         // water table depth [m]
  DECLARE_DEVICE_VIEW(1DR8, Qcharge)     // aquifer recharge rate [mm/s]
  DECLARE_DEVICE_VIEW(1DI4, Jwt)         // layer index above water table [-]
  DECLARE_DEVICE_VIEW(1DR8, QflxDeficit) // water deficit flux [mm/s]

  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase
  // Also includes soil data for pervious road
  PerviousRoadDataType(int numLandunits, int numRadBands, int numRadTypes,
                       int numSoilLayers)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes,
                               numSoilLayers),
        soil(numLandunits, numSoilLayers) {
    ALLOCATE_DEVICE_VIEW(H2OSoiLiq, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(H2OSoiIce, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(H2OSoiVol, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(EffPorosity, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Hk, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Smp, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(QflxInfl, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(QflxTran, Array2DR8, numLandunits, numSoilLayers)
    ALLOCATE_DEVICE_VIEW(Zwt, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Qcharge, Array1DR8, numLandunits)
    ALLOCATE_DEVICE_VIEW(Jwt, Array1DI4, numLandunits)
    ALLOCATE_DEVICE_VIEW(QflxDeficit, Array1DR8, numLandunits)
  }
};

struct WallDataType : SurfaceDataBase {
  DECLARE_DEVICE_VIEW(
      3DR8,
      DownwellingShortRad) // downwelling shortwave radiation per unit wall area
  // Inherits common radiative fields from SurfaceDataBase

  WallDataType(int numLandunits, int numRadBands, int numRadTypes,
               int numLayers)
      : SurfaceDataBase(numLandunits, numRadBands, numRadTypes, numLayers) {
    ALLOCATE_DEVICE_VIEW(DownwellingShortRad, Array3DR8, numLandunits,
                         numRadBands, numRadTypes)
  }
};

struct RoofDataType : SnowCoveredSurfaceData {
  // Inherits all fields from SnowCoveredSurfaceData and SurfaceDataBase

  RoofDataType(int numLandunits, int numRadBands, int numRadTypes,
               int numLayers)
      : SnowCoveredSurfaceData(numLandunits, numRadBands, numRadTypes,
                               numLayers) {}
};
} // namespace URBANXX

#endif // URBAN_SURFACE_TYPE_IMPL_H