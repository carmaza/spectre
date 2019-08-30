// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/NewtonianEuler/Sources.hpp"
#include "Evolution/Systems/NewtonianEuler/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "tests/Unit/Pypp/CheckWithRandomValues.hpp"
#include "tests/Unit/Pypp/SetupLocalPythonEnvironment.hpp"

namespace {

namespace Tags {
template <typename>
struct Source;
}  // namespace Tags

struct FirstArg : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() noexcept { return "FirstArg"; }
};

template <size_t Dim>
struct SecondArg : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim>;
  static std::string name() noexcept { return "SecondArg"; }
};

struct ThirdArg : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() noexcept { return "ThirdArg"; }
};

template <size_t Dim>
struct FourthArg : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim>;
  static std::string name() noexcept { return "FourthArg"; }
};

template <size_t Dim>
struct SomeSourceType {
  static constexpr size_t volume_dim = Dim;

  using return_tags = tmpl::list<
      Tags::Source<NewtonianEuler::Tags::MomentumDensity<DataVector, Dim>>,
      Tags::Source<NewtonianEuler::Tags::EnergyDensity<DataVector>>>;
  using argument_tags =
      tmpl::list<FirstArg, SecondArg<Dim>, ThirdArg, FourthArg<Dim>>;

  void apply(
      const gsl::not_null<tnsr::I<DataVector, Dim>*> source_momentum_density,
      const gsl::not_null<Scalar<DataVector>*> source_energy_density,
      const Scalar<DataVector>& first_arg,
      const tnsr::I<DataVector, Dim>& second_arg,
      const Scalar<DataVector>& third_arg,
      const tnsr::i<DataVector, Dim>& fourth_arg) const noexcept {
    for (size_t i = 0; i < Dim; ++i) {
      source_momentum_density->get(i) =
          (get(first_arg) - 1.5 * get(third_arg)) * second_arg.get(i);
    }
    get(*source_energy_density) =
        get(dot_product(second_arg, fourth_arg)) + 3.0 * get(third_arg);
  }
};

template <typename SourceTermType>
struct ComputeSourcesProxy {
  static constexpr size_t dim = SourceTermType::volume_dim;
  static void apply(
      const gsl::not_null<tnsr::I<DataVector, dim>*> source_momentum_density,
      const gsl::not_null<Scalar<DataVector>*> source_energy_density,
      const Scalar<DataVector>& first_arg,
      const tnsr::I<DataVector, dim>& second_arg,
      const Scalar<DataVector>& third_arg,
      const tnsr::i<DataVector, dim>& fourth_arg) noexcept {
    SourceTermType source_computer;
    NewtonianEuler::ComputeSources<SourceTermType>::apply(
        source_momentum_density, source_energy_density, source_computer,
        first_arg, second_arg, third_arg, fourth_arg);
  }
};

template <size_t Dim>
void test_sources(const DataVector& used_for_size) {
  pypp::check_with_random_values<4>(
      &ComputeSourcesProxy<SomeSourceType<Dim>>::apply, "TestFunctions",
      {"momentum_density_source", "energy_density_source"},
      {{{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}}}, used_for_size);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.NewtonianEuler.Sources",
                  "[Unit][Evolution]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "Evolution/Systems/NewtonianEuler"};

  GENERATE_UNINITIALIZED_DATAVECTOR;
  CHECK_FOR_DATAVECTORS(test_sources, (1, 2, 3))
}
