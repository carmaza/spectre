# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY NewtonianEuler)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Characteristics.cpp
  ComputeItems.cpp
  ConservativeFromPrimitive.cpp
  Fluxes.cpp
  PrimitiveFromConservative.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Characteristics.hpp
  ComputeItems.hpp
  ConservativeFromPrimitive.hpp
  Fluxes.hpp
  PrimitiveFromConservative.hpp
  Sources.hpp
  System.hpp
  Tags.hpp
  TagsDeclarations.hpp
  TimeDerivativeTerms.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC DataStructures
  INTERFACE ErrorHandling
  INTERFACE Hydro
  )

add_subdirectory(Limiters)
add_subdirectory(NumericalFluxes)
add_subdirectory(Sources)
