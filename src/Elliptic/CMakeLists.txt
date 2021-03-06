# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY Elliptic)

add_spectre_library(${LIBRARY} INTERFACE)

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  FirstOrderComputeTags.hpp
  FirstOrderOperator.hpp
  Tags.hpp
  )

target_link_libraries(
  ${LIBRARY}
  INTERFACE
  DataStructures
  Domain
  LinearOperators
  Utilities
  )

add_subdirectory(Actions)
add_subdirectory(DiscontinuousGalerkin)
add_subdirectory(Executables)
add_subdirectory(Systems)
add_subdirectory(Triggers)
