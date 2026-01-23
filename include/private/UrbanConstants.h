#ifndef URBAN_CONSTANTS_H
#define URBAN_CONSTANTS_H

namespace URBANXX {

// Radiation band indices
enum RadiationBand {
  VIS = 0, // Visible band
  NIR = 1  // Near-infrared band
};

// Radiation type indices
enum RadiationType {
  DIRECT = 0, // Direct radiation
  DIFFUSE = 1 // Diffuse radiation
};

// Constants
constexpr int NUM_RAD_BANDS = 2;        // vis and nir
constexpr int NUM_RAD_TYPES = 2;        // direct and diffuse
constexpr int NUM_VERTICAL_LEVELS = 15; // number of vertical levels

#define STEBOL 5.670374419e-8

} // namespace URBANXX

#endif // URBAN_CONSTANTS_H
