#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSetterHelpers.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"
#include "private/UrbanViewFactorImpl.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

// Helper function to compute view factors from canyon height-to-width ratio
static void ComputeViewFactors(UrbanType urban, UrbanErrorCode *status) {
  using namespace URBANXX;

  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    ViewFactorComputer computer(urban);
    computer.run();
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

extern "C" {

using namespace URBANXX;

void UrbanSetFracPervRoadOfTotalRoad(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  // Validate that the fraction values are within [0, 1]
  for (int i = 0; i < length; ++i) {
    const double v = values[i];
    if (v < 0.0 || v > 1.0) {
      if (status) {
        *status = URBAN_ERR_INVALID_ARGUMENT;
      }
      return;
    }
  }

  SetView1D(urban->urbanParams.FracPervRoadOfTotalRoad, values, length, status);
}

void UrbanSetCanyonHwr(UrbanType urban, const double *values, int length,
                       UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  // Set the canyon height-to-width ratio using the template function
  SetView1D(urban->urbanParams.CanyonHwr, values, length, status);

  // If the set operation failed, return early
  if (*status != URBAN_SUCCESS) {
    return;
  }

  // Compute the view factors based on the canyon geometry
  ComputeViewFactors(urban, status);
}

void UrbanSetWtRoof(UrbanType urban, const double *values, int length,
                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.WtRoof, values, length, status);
}

void UrbanSetAlbedoPerviousRoad(UrbanType urban, const double *values,
                                const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView3D(urban->urbanParams.albedo.PerviousRoad, values, size, status);
}

void UrbanSetAlbedoImperviousRoad(UrbanType urban, const double *values,
                                  const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView3D(urban->urbanParams.albedo.ImperviousRoad, values, size, status);
}

void UrbanSetAlbedoSunlitWall(UrbanType urban, const double *values,
                              const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView3D(urban->urbanParams.albedo.SunlitWall, values, size, status);
}

void UrbanSetAlbedoShadedWall(UrbanType urban, const double *values,
                              const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView3D(urban->urbanParams.albedo.ShadedWall, values, size, status);
}

void UrbanSetAlbedoRoof(UrbanType urban, const double *values,
                        const int size[3], UrbanErrorCode *status) {
  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  SetView3D(urban->urbanParams.albedo.Roof, values, size, status);
}

// Emissivity setter functions
void UrbanSetEmissivityPerviousRoad(UrbanType urban, const double *values,
                                    int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.emissivity.PerviousRoad, values, length, status);
}

void UrbanSetEmissivityImperviousRoad(UrbanType urban, const double *values,
                                      int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.emissivity.ImperviousRoad, values, length,
            status);
}

void UrbanSetEmissivityWall(UrbanType urban, const double *values, int length,
                            UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.emissivity.Wall, values, length, status);
}

void UrbanSetEmissivityRoof(UrbanType urban, const double *values, int length,
                            UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.emissivity.Roof, values, length, status);
}

// Thermal conductivity setter functions
void UrbanSetThermalConductivityRoad(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.tk.Road, values, length, status);
}

void UrbanSetThermalConductivityWall(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.tk.Wall, values, length, status);
}

void UrbanSetThermalConductivityRoof(UrbanType urban, const double *values,
                                     int length, UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.tk.Roof, values, length, status);
}

// Heat capacity setter functions
void UrbanSetHeatCapacityRoad(UrbanType urban, const double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.cv.Road, values, length, status);
}

void UrbanSetHeatCapacityWall(UrbanType urban, const double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.cv.Wall, values, length, status);
}

void UrbanSetHeatCapacityRoof(UrbanType urban, const double *values, int length,
                              UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.cv.Roof, values, length, status);
}

// Height parameter setter functions
void UrbanSetForcHgtT(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.heights.ForcHgtT, values, length, status);
}

void UrbanSetForcHgtU(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.heights.ForcHgtU, values, length, status);
}

void UrbanSetZDTown(UrbanType urban, const double *values, int length,
                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.heights.ZDTown, values, length, status);
}

void UrbanSetZ0Town(UrbanType urban, const double *values, int length,
                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.heights.Z0Town, values, length, status);
}

void UrbanSetHtRoof(UrbanType urban, const double *values, int length,
                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.heights.HtRoof, values, length, status);
}

void UrbanSetWindHgtCanyon(UrbanType urban, const double *values, int length,
                           UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->urbanParams.heights.WindHgtCanyon, values, length, status);
}

} // extern "C"
