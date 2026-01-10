#include "Urban.h"
#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

// Base template function for setting 1D host views
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

    // Deep copy from the temporary host view to the host view
    Kokkos::deep_copy(view, values_view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

extern "C" {

using namespace URBANXX;

void UrbanSetAtmTemp(UrbanType urban, const double *values, int length,
                     UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcTempH, values, length, status);
}

void UrbanSetAtmPotTemp(UrbanType urban, const double *values, int length,
                        UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcPotTempH, values, length, status);
}

void UrbanSetAtmRho(UrbanType urban, const double *values, int length,
                    UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcRhoH, values, length, status);
}

void UrbanSetAtmSpcHumd(UrbanType urban, const double *values, int length,
                        UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcSpcHumdH, values, length, status);
}

void UrbanSetAtmPress(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcPressH, values, length, status);
}

void UrbanSetAtmWindU(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcWindUH, values, length, status);
}

void UrbanSetAtmWindV(UrbanType urban, const double *values, int length,
                      UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcWindVH, values, length, status);
}

void UrbanSetAtmCoszen(UrbanType urban, const double *values, int length,
                       UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.CoszenH, values, length, status);
}

void UrbanSetAtmFracSnow(UrbanType urban, const double *values, int length,
                         UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.FracSnowH, values, length, status);
}

void UrbanSetAtmLongwaveDown(UrbanType urban, const double *values, int length,
                             UrbanErrorCode *status) {
  if (!ValidateInputsWithData(urban, values, status))
    return;

  SetView1D(urban->atmosphereData.ForcLRadH, values, length, status);
}

void UrbanSetAtmShortwaveDown(UrbanType urban, const double *values,
                              const int size[3], UrbanErrorCode *status) {
  using namespace URBANXX;

  if (!ValidateInputsWithSize(urban, values, size, status))
    return;

  try {
    auto &h_ForcSRad = urban->atmosphereData.ForcSRadH;

    // Check if each dimension matches the view extent
    if (size[0] != static_cast<int>(h_ForcSRad.extent(0)) ||
        size[1] != static_cast<int>(h_ForcSRad.extent(1)) ||
        size[2] != static_cast<int>(h_ForcSRad.extent(2))) {
      *status = URBAN_ERR_SIZE_MISMATCH;
      return;
    }

    // Create an unmanaged host view from the input array
    auto values_view = Kokkos::View<const double ***, Kokkos::HostSpace,
                                    Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        values, size[0], size[1], size[2]);

    // Deep copy from the temporary host view to the host view
    Kokkos::deep_copy(h_ForcSRad, values_view);

    *status = URBAN_SUCCESS;
  } catch (...) {
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
