# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_CurvedScalarWave")

set(LIBRARY_SOURCES
  Test_Constraints.cpp
  Test_Characteristics.cpp
  Test_DuDt.cpp
  Test_Tags.cpp
  Test_UpwindFlux.cpp
  )

add_test_library(
  ${LIBRARY}
  "Evolution/Systems/CurvedScalarWave/"
  "${LIBRARY_SOURCES}"
  "CurvedScalarWave;CurvedScalarWaveHelpers;GeneralRelativityHelpers"
  )
