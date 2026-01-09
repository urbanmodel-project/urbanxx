#ifndef URBAN_CONSTANTS_H
#define URBAN_CONSTANTS_H

namespace URBANXX {

// Radiation band indices
enum RadiationBand {
  VIS = 0, // Visible band
  NIR = 1  // Near-infrared band
};

// Constants
constexpr int NUM_RAD_BANDS = 2; // vis and nir
constexpr int NUM_RAD_TYPES = 2; // direct and diffuse

} // namespace URBANXX

#endif // URBAN_CONSTANTS_H
