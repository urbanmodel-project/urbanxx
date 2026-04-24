#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanGetterHelpers.h"
#include "private/UrbanSetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

void UrbanSetInfiltrationFluxForPerviousRoad(UrbanType urban,
                                             const double *values, int length,
                                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.QflxInfl, values, length, status);
}

void UrbanSetTranspirationFluxForPerviousRoad(UrbanType urban,
                                              const double *values,
                                              const int size[2],
                                              UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->perviousRoad.QflxTran, values, size, status);
}

void UrbanSetWaterTableDepth(UrbanType urban, const double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.Zwt, values, length, status);
}

void UrbanSetSoilLiquidWaterForPerviousRoad(UrbanType urban,
                                            const double *values,
                                            const int size[2],
                                            UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->perviousRoad.H2OSoiLiq, values, size, status);
}

void UrbanSetSoilIceContentForPerviousRoad(UrbanType urban,
                                           const double *values,
                                           const int size[2],
                                           UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView2D(urban->perviousRoad.H2OSoiIce, values, size, status);
}

// Infiltration input setters
void UrbanSetSurfaceRunoffForPerviousRoad(UrbanType urban, const double *values,
                                          int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.QflxSurf, values, length, status);
}

void UrbanSetGroundEvapFluxForPerviousRoad(UrbanType urban,
                                           const double *values, int length,
                                           UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.QflxEvapGrnd, values, length, status);
}

// Compute
void UrbanComputeInfiltration(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  const int nlandunits = urban->numLandunits;
  auto qflx_top_soil = urban->atmosphereData.ForcRain;
  auto qflx_surf = urban->perviousRoad.QflxSurf;
  auto qflx_evap_grnd = urban->perviousRoad.QflxEvapGrnd;
  auto qflx_infl = urban->perviousRoad.QflxInfl;

  Kokkos::parallel_for(
      "UrbanComputeInfiltration", nlandunits, KOKKOS_LAMBDA(const int l) {
        qflx_infl(l) = qflx_top_soil(l) - qflx_surf(l) - qflx_evap_grnd(l);
      });
  Kokkos::fence();
  *status = URBAN_SUCCESS;
}

// Getter
void UrbanGetInfiltrationFluxPerviousRoad(UrbanType urban, double *values,
                                          int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  GetView1D(urban->perviousRoad.QflxInfl, values, length, status);
}

// ============================================================================
// WaterTable setter functions
// ============================================================================

void UrbanSetAquiferWaterForPerviousRoad(UrbanType urban, const double *values,
                                         int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.Wa, values, length, status);
}

void UrbanSetFracH2osfcForPerviousRoad(UrbanType urban, const double *values,
                                       int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.FracH2osfc, values, length, status);
}

void UrbanSetDewGrndFluxForPerviousRoad(UrbanType urban, const double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.QflxDewGrnd, values, length, status);
}

void UrbanSetDewSnowFluxForPerviousRoad(UrbanType urban, const double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.QflxDewSnow, values, length, status);
}

void UrbanSetSubSnowFluxForPerviousRoad(UrbanType urban, const double *values,
                                        int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.QflxSubSnow, values, length, status);
}

void UrbanSetQchargeForPerviousRoad(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->perviousRoad.Qcharge, values, length, status);
}

} // extern "C"
