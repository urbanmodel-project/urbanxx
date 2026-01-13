#include "private/AtmosphereTypeImpl.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanConstants.h"
#include "private/UrbanDebugUtils.h"
#include "private/UrbanParamsTypeImpl.h"
#include "private/UrbanSurfaceFluxesImpl.h"
#include "private/UrbanSurfaceTypeImpl.h"
#include "private/UrbanTypeImpl.h"
#include <Kokkos_Core.hpp>
#include <cmath>
#include <iostream>

namespace URBANXX {

#define SHR_CONST_TKFRZ 273.15

KOKKOS_INLINE_FUNCTION
void QSat(Real T, Real p, Real &es, Real &esdT, Real &qs, Real &qsdT) {
  // For water vapor (temperature range 0C-100C)
  Real a0 = 6.11213476;
  Real a1 = 0.444007856;
  Real a2 = 0.143064234e-01;
  Real a3 = 0.264461437e-03;
  Real a4 = 0.305903558e-05;
  Real a5 = 0.196237241e-07;
  Real a6 = 0.892344772e-10;
  Real a7 = -0.373208410e-12;
  Real a8 = 0.209339997e-15;
  // For derivative:water vapor
  Real b0 = 0.444017302;
  Real b1 = 0.286064092e-01;
  Real b2 = 0.794683137e-03;
  Real b3 = 0.121211669e-04;
  Real b4 = 0.103354611e-06;
  Real b5 = 0.404125005e-09;
  Real b6 = -0.788037859e-12;
  Real b7 = -0.114596802e-13;
  Real b8 = 0.381294516e-16;
  // For ice (temperature range -75C-0C)
  Real c0 = 6.11123516;
  Real c1 = 0.503109514;
  Real c2 = 0.188369801e-01;
  Real c3 = 0.420547422e-03;
  Real c4 = 0.614396778e-05;
  Real c5 = 0.602780717e-07;
  Real c6 = 0.387940929e-09;
  Real c7 = 0.149436277e-11;
  Real c8 = 0.262655803e-14;
  // For derivative:ice
  Real d0 = 0.503277922;
  Real d1 = 0.377289173e-01;
  Real d2 = 0.126801703e-02;
  Real d3 = 0.249468427e-04;
  Real d4 = 0.313703411e-06;
  Real d5 = 0.257180651e-08;
  Real d6 = 0.133268878e-10;
  Real d7 = 0.394116744e-13;
  Real d8 = 0.498070196e-16;

  Real T_limit = T - SHR_CONST_TKFRZ;
  if (T_limit > 100.0)
    T_limit = 100.0;
  if (T_limit < -75.0)
    T_limit = -75.0;

  Real td = T_limit;

  if (td >= 0.0) {
    es =
        a0 +
        td * (a1 +
              td * (a2 +
                    td * (a3 +
                          td * (a4 +
                                td * (a5 + td * (a6 + td * (a7 + td * a8)))))));

    esdT =
        b0 +
        td * (b1 +
              td * (b2 +
                    td * (b3 +
                          td * (b4 +
                                td * (b5 + td * (b6 + td * (b7 + td * b8)))))));

  } else {
    es =
        c0 +
        td * (c1 +
              td * (c2 +
                    td * (c3 +
                          td * (c4 +
                                td * (c5 + td * (c6 + td * (c7 + td * c8)))))));

    esdT =
        d0 +
        td * (d1 +
              td * (d2 +
                    td * (d3 +
                          td * (d4 +
                                td * (d5 + td * (d6 + td * (d7 + td * d8)))))));
  }

  es = es * 100.0;     // pa
  esdT = esdT * 100.0; // pa/K

  Real vp = 1.0 / (p - 0.378 * es);
  Real vp1 = 0.622 * vp;
  Real vp2 = vp1 * vp;

  qs = es * vp1;         // kg/kg
  qsdT = esdT * vp2 * p; // 1 / K
}

// Compute surface fluxes for all urban surfaces
void ComputeSurfaceFluxes(URBANXX::_p_UrbanType &urban) {
  const int numLandunits = urban.numLandunits;

  std::cout << "ComputeSurfaceFluxes: Processing " << numLandunits
            << " landunits (placeholder)" << std::endl;
}

} // namespace URBANXX
