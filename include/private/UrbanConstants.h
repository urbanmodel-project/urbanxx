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
constexpr int NUM_RAD_BANDS = 2;    // vis and nir
constexpr int NUM_RAD_TYPES = 2;    // direct and diffuse
constexpr int NUM_URBAN_LAYERS = 5; // number of vertical levels
constexpr int NUM_SOIL_LAYERS = 15; // number of soil layers

// Physical constants
constexpr double GRAVITY = 9.80616; // acceleration due to gravity [m/s^2]
constexpr double CPAIR = 1004.64;   // specific heat of dry air [J/kg/K]
constexpr double SHR_CONST_TKFRZ = 273.15; // freezing temperature of water [K]
constexpr double STEBOL =
    5.670374419e-8; // Stefan-Boltzmann constant [W/m^2/K^4]
constexpr double SHR_CONST_RHOICE = 0.917e3; // density of ice [kg/m^3]
constexpr double SHR_CONST_RHOWATER = 1.0e3; // density of water [kg/m^3]

} // namespace URBANXX

#endif // URBAN_CONSTANTS_H
