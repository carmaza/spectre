// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>
#include <tuple>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "ErrorHandling/Error.hpp"
#include "Evolution/Systems/NewtonianEuler/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/AnalyticData/NewtonianEuler/KhInstability.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {

template <size_t Dim>
struct KhInstabilityProxy : NewtonianEuler::AnalyticData::KhInstability<Dim> {
  using NewtonianEuler::AnalyticData::KhInstability<Dim>::KhInstability;

  template <typename DataType>
  using variables_tags =
      tmpl::list<NewtonianEuler::Tags::MassDensity<DataType>,
                 NewtonianEuler::Tags::Velocity<DataType, Dim, Frame::Inertial>,
                 NewtonianEuler::Tags::SpecificInternalEnergy<DataType>,
                 NewtonianEuler::Tags::Pressure<DataType>>;

  template <typename DataType>
  tuples::tagged_tuple_from_typelist<variables_tags<DataType>>
  primitive_variables(const tnsr::I<DataType, Dim, Frame::Inertial>& x) const
      noexcept {
    return this->variables(x, variables_tags<DataType>{});
  }
};

template <size_t Dim, typename DataType>
void test_analytic_data(const DataType& used_for_size) noexcept {
  const double adiabatic_index = 1.43;
  const double strip_bimedian = 0.5;
  const double strip_thickness = 0.4;
  const double strip_density = 2.1;
  const double strip_velocity = 0.3;
  const double background_density = 2.0;
  const double background_velocity = -0.2;
  const double pressure = 1.1;
  const double perturbation_amplitude = 0.1;
  const double damping_factor = 0.01;

  KhInstabilityProxy<Dim> kh_inst(
      adiabatic_index, strip_bimedian, strip_thickness, strip_density,
      strip_velocity, background_density, background_velocity, pressure,
      perturbation_amplitude, damping_factor);
  pypp::check_with_random_values<
      1, typename KhInstabilityProxy<Dim>::template variables_tags<DataType>>(
      &KhInstabilityProxy<Dim>::template primitive_variables<DataType>, kh_inst,
      "KhInstability",
      {"mass_density", "velocity", "specific_internal_energy", "pressure"},
      {{{0.0, 1.0}}},
      std::make_tuple(adiabatic_index, strip_bimedian, strip_thickness,
                      strip_density, strip_velocity, background_density,
                      background_velocity, pressure, perturbation_amplitude,
                      damping_factor),
      used_for_size);

  const auto kh_inst_from_options = TestHelpers::test_creation<
      NewtonianEuler::AnalyticData::KhInstability<Dim>>(
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  CHECK(kh_inst_from_options == kh_inst);

  KhInstabilityProxy<Dim> kh_inst_to_move(
      adiabatic_index, strip_bimedian, strip_thickness, strip_density,
      strip_velocity, background_density, background_velocity, pressure,
      perturbation_amplitude, damping_factor);
  test_move_semantics(std::move(kh_inst_to_move), kh_inst);  //  NOLINT

  test_serialization(kh_inst);
}

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability",
    "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/AnalyticData/NewtonianEuler"};

  GENERATE_UNINITIALIZED_DOUBLE_AND_DATAVECTOR;
  CHECK_FOR_DOUBLES_AND_DATAVECTORS(test_analytic_data, (2, 3));
}

template <size_t Dim>
struct Instability {
  using type = NewtonianEuler::AnalyticData::KhInstability<Dim>;
  static constexpr OptionString help = {
      "Initial data to simulate the KH instability."};
};

// [[OutputRegex, In string:.*At line 5 column 17:.Value -2.1 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.RhoOut2d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<2>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: -2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<2>>();
}

// [[OutputRegex, In string:.*At line 5 column 17:.Value -2.1 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.RhoOut3d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<3>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: -2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<3>>();
}

// [[OutputRegex, In string:.*At line 7 column 22:.Value -2 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.RhoIn2d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<2>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: -2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<2>>();
}

// [[OutputRegex, In string:.*At line 7 column 22:.Value -2 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.RhoIn3d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<3>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: -2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<3>>();
}

// [[OutputRegex, In string:.*At line 9 column 13:.Value -1.1 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.Pressure2d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<2>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: -1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<2>>();
}

// [[OutputRegex, In string:.*At line 9 column 13:.Value -1.1 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.Pressure3d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<3>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: -1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<3>>();
}

// [[OutputRegex, In string:.*At line 11 column 18:.Value -0.01 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.DampFactor2d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<2>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: -0.01");
  test_options.get<Instability<2>>();
}

// [[OutputRegex, In string:.*At line 11 column 18:.Value -0.01 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.DampFactor3d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<3>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: 0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: -0.01");
  test_options.get<Instability<3>>();
}

// [[OutputRegex, In string:.*At line 4 column 19:.Value -0.4 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.StripThick2d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<2>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: -0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<2>>();
}

// [[OutputRegex, In string:.*At line 4 column 19:.Value -0.4 is below the
// lower bound of 0]]
SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticData.NewtEuler.KhInstability.StripThick3d",
    "[Unit][PointwiseFunctions]") {
  ERROR_TEST();
  Options<tmpl::list<Instability<3>>> test_options("");
  test_options.parse(
      "Instability:\n"
      "  AdiabaticIndex: 1.43\n"
      "  StripBimedian: 0.5\n"
      "  StripThickness: -0.4\n"
      "  StripDensity: 2.1\n"
      "  StripVelocity: 0.3\n"
      "  BackgroundDensity: 2.0\n"
      "  BackgroundVelocity: -0.2\n"
      "  Pressure: 1.1\n"
      "  PerturbAmplitude: 0.1\n"
      "  DampingFactor: 0.01");
  test_options.get<Instability<3>>();
}
