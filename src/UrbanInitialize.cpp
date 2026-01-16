#include "Urban.h"
#include "private/DataTypesImpl.h"
#include "private/UrbanInitializeImpl.h"
#include "private/UrbanTypeImpl.h"
#include "private/UrbanValidation.h"
#include <iostream>

// Define the C struct to match the C++ class
struct _p_UrbanType : public URBANXX::_p_UrbanType {
  using URBANXX::_p_UrbanType::_p_UrbanType;
};

extern "C" {

using namespace URBANXX;

void UrbanInitializeTemperature(UrbanType urban, UrbanErrorCode *status) {
  if (urban == nullptr || status == nullptr) {
    if (status)
      *status = URBAN_ERR_INVALID_ARGUMENT;
    return;
  }

  try {
    std::cout << "DEBUG: Creating UrbanTemperatureInitializer..." << std::endl;
    // Create initializer object and run
    UrbanTemperatureInitializer initializer(urban);
    std::cout << "DEBUG: Calling initializer.run()..." << std::endl;
    initializer.run();
    std::cout << "DEBUG: initializer.run() completed" << std::endl;

    *status = URBAN_SUCCESS;
  } catch (...) {
    std::cout << "DEBUG: Exception caught in UrbanInitializeTemperature"
              << std::endl;
    *status = URBAN_ERR_INTERNAL;
  }
}

} // extern "C"
