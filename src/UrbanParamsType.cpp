#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

// Base template function for setting 3D views
template <typename ViewType>
static void SetView3D(ViewType &view, const double *values, const int size[3],
                      UrbanErrorCode *status) {
  using namespace URBANXX;

  if (values == nullptr || size == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Check if each dimension matches the view extent
    if (size[0] != static_cast<int>(view.extent(0)) ||
        size[1] != static_cast<int>(view.extent(1)) ||
        size[2] != static_cast<int>(view.extent(2))) {
      *status = URBAN_ERR_SIZE_MISMATCH;
      return;
    }

    // Create an unmanaged host view from the input array
    auto values_view = Kokkos::View<const double ***, Kokkos::HostSpace,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        values, size[0], size[1], size[2]);

    // Deep copy from the temporary host view to the device view
    Kokkos::deep_copy(view, values_view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

// Base template function for setting 1D views
template <typename ViewType>
static void SetView1D(ViewType &view, const double *values, int length,
                      UrbanErrorCode *status) {
  using namespace URBANXX;

  if (values == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    // Check if the length matches the view extent
    if (length != static_cast<int>(view.extent(0))) {
      *status = URBAN_ERR_SIZE_MISMATCH;
      return;
    }

    // Create an unmanaged host view from the input array
    auto values_view =
        Kokkos::View<const double *, Kokkos::HostSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>(values, length);

    // Deep copy from the temporary host view to the device view
    Kokkos::deep_copy(view, values_view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

// Helper function to compute view factors from canyon height-to-width ratio
static void ComputeViewFactors(UrbanType urban, UrbanErrorCode *status) {
  using namespace URBANXX;

  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    auto &CanyonHwr = urban->urbanParams.CanyonHwr;
    auto &sr = urban->urbanParams.viewFactor.SkyFrmRoad;
    auto &sw = urban->urbanParams.viewFactor.SkyFrmWall;
    auto &rw = urban->urbanParams.viewFactor.RoadFrmWall;
    auto &wr = urban->urbanParams.viewFactor.WallFrmRoad;
    auto &ww = urban->urbanParams.viewFactor.OtherWallFrmWall;

    Kokkos::parallel_for(
        "ComputingViewFactor", urban->numLandunits, KOKKOS_LAMBDA(int l) {
          const Real hwr = CanyonHwr(l);
          const Real sqrt_term = sqrtf(hwr * hwr + 1.0);

          sr(l) = sqrt_term - hwr;                     // eqn 2.25
          wr(l) = 0.5 * (1.0 - sr(l));                 // eqn 2.27
          sw(l) = 0.5 * (hwr + 1.0 - sqrt_term) / hwr; // eqn 2.24
          rw(l) = sw(l);                               // eqn 2.27
          ww(l) = 1.0 - sw(l) - rw(l);                 // eqn 2.28
        });
    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

extern "C" {

using namespace URBANXX;

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

} // extern "C"
