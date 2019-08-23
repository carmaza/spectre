# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY NewtonianEulerAnalyticData)

set(LIBRARY_SOURCES
  RayleighTaylorInst.cpp
  )

add_spectre_library(${LIBRARY} ${LIBRARY_SOURCES})

target_link_libraries(
  ${LIBRARY}
  INTERFACE DataStructures
  INTERFACE ErrorHandling
  INTERFACE Hydro
  )
