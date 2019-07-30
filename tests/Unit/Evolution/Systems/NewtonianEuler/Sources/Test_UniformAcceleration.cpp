// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "Evolution/Systems/NewtonianEuler/Sources/UniformAcceleration.hpp"
#include "tests/Unit/Pypp/CheckWithRandomValues.hpp"
#include "tests/Unit/Pypp/SetupLocalPythonEnvironment.hpp"

namespace {

template <size_t Dim>
void test_sources(const std::array<double, Dim>& acceleration_field,
                  const DataVector& used_for_size) noexcept {
  NewtonianEuler::Sources::UniformAcceleration<Dim> source(acceleration_field);
  pypp::check_with_random_values<2>(
      &NewtonianEuler::Sources::UniformAcceleration<Dim>::apply, source,
      "UniformAcceleration",
      {"momentum_density_source", "energy_density_source"},
      {{{0.0, 3.0}, {-1.0, 1.0}}}, std::make_tuple(acceleration_field),
      used_for_size);
}

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.NewtonianEuler.Sources.UniformAcceleration",
    "[Unit][Evolution]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "Evolution/Systems/NewtonianEuler/Sources"};

  test_sources<1>({{-2.0}}, DataVector(5));
  test_sources<2>({{-1.0, 0.4}}, DataVector(5));
  test_sources<3>({{0.7, -1.2, 4.5}}, DataVector(5));
}
