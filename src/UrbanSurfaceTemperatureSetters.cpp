#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanSetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

// Layer temperature setter functions (2D)
void UrbanSetLayerTempRoof(UrbanType urban, const double *values,
                           const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->roof.Temperature, values, size, status);
  if (*status == URBAN_SUCCESS) {
    auto col0 = Kokkos::subview(urban->roof.Temperature, Kokkos::ALL(), 0);
    Kokkos::deep_copy(urban->roof.EffectiveSurfTemp, col0);
    Kokkos::deep_copy(urban->roof.TGrnd0, col0);
  }
}

void UrbanSetLayerTempImperviousRoad(UrbanType urban, const double *values,
                                     const int size[2],
                                     UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->imperviousRoad.Temperature, values, size, status);
  if (*status == URBAN_SUCCESS) {
    auto col0 =
        Kokkos::subview(urban->imperviousRoad.Temperature, Kokkos::ALL(), 0);
    Kokkos::deep_copy(urban->imperviousRoad.EffectiveSurfTemp, col0);
    Kokkos::deep_copy(urban->imperviousRoad.TGrnd0, col0);
  }
}

void UrbanSetLayerTempPerviousRoad(UrbanType urban, const double *values,
                                   const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->perviousRoad.Temperature, values, size, status);
  if (*status == URBAN_SUCCESS) {
    auto col0 =
        Kokkos::subview(urban->perviousRoad.Temperature, Kokkos::ALL(), 0);
    Kokkos::deep_copy(urban->perviousRoad.EffectiveSurfTemp, col0);
    Kokkos::deep_copy(urban->perviousRoad.TGrnd0, col0);
  }
}

void UrbanSetLayerTempSunlitWall(UrbanType urban, const double *values,
                                 const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->sunlitWall.Temperature, values, size, status);
  if (*status == URBAN_SUCCESS) {
    auto col0 =
        Kokkos::subview(urban->sunlitWall.Temperature, Kokkos::ALL(), 0);
    Kokkos::deep_copy(urban->sunlitWall.EffectiveSurfTemp, col0);
    Kokkos::deep_copy(urban->sunlitWall.TGrnd0, col0);
  }
}

void UrbanSetLayerTempShadedWall(UrbanType urban, const double *values,
                                 const int size[2], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->shadedWall.Temperature, values, size, status);
  if (*status == URBAN_SUCCESS) {
    auto col0 =
        Kokkos::subview(urban->shadedWall.Temperature, Kokkos::ALL(), 0);
    Kokkos::deep_copy(urban->shadedWall.EffectiveSurfTemp, col0);
    Kokkos::deep_copy(urban->shadedWall.TGrnd0, col0);
  }
}

// Canyon air property setter functions (1D)
void UrbanSetCanyonAirTemperature(UrbanType urban, const double *values,
                                  int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanCanyon.Taf, values, length, status);
}

void UrbanSetCanyonSpecificHumidity(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanCanyon.Qaf, values, length, status);
}

// Building temperature setter function (1D)
void UrbanSetBuildingTemperature(UrbanType urban, const double *values,
                                 int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->building.Temperature, values, length, status);
}

// Surface wetness setter functions (1D)
void UrbanSetFractionWetImperviousRoad(UrbanType urban, const double *values,
                                       int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  // Validate that fraction wet values are within [0, 1]
  if (!ValidateRange(values, length, 0.0, 1.0, status))
    return;

  SetView1D(urban->imperviousRoad.FractionWet, values, length, status);
}

void UrbanSetFractionWetRoof(UrbanType urban, const double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  // Validate that fraction wet values are within [0, 1]
  if (!ValidateRange(values, length, 0.0, 1.0, status))
    return;

  SetView1D(urban->roof.FractionWet, values, length, status);
}

} // extern "C"
