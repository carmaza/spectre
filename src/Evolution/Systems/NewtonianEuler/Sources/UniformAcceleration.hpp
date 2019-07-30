// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/NewtonianEuler/TagsDeclarations.hpp"
#include "Options/Options.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;

namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl

namespace Tags {
template <typename>
struct Source;
}  // namespace Tags
/// \endcond

namespace NewtonianEuler {
namespace Sources {

/*!
 * \brief Compute the source terms for an external uniform acceleration.
 */
template <size_t Dim>
struct UniformAcceleration {
  /// The uniform external acceleration
  struct AccelerationField {
    using type = std::array<double, Dim>;
    static constexpr OptionString help = {"The uniform external acceleration."};
  };

  using options = tmpl::list<AccelerationField>;

  static constexpr size_t volume_dim = Dim;

  static constexpr OptionString help = {
      "Sources generated from a uniform external acceleration"};

  explicit UniformAcceleration(
      std::array<double, Dim> acceleration_field) noexcept;

  UniformAcceleration() noexcept = default;
  UniformAcceleration(const UniformAcceleration& /*rhs*/) = delete;
  UniformAcceleration& operator=(const UniformAcceleration& /*rhs*/) = delete;
  UniformAcceleration(UniformAcceleration&& /*rhs*/) noexcept = default;
  UniformAcceleration& operator=(UniformAcceleration&& /*rhs*/) noexcept =
      default;
  ~UniformAcceleration() = default;

  // clang-tidy: google-runtime-references
  void pup(PUP::er& /*p*/) noexcept;  // NOLINT

  using return_tags = tmpl::list<
      ::Tags::Source<NewtonianEuler::Tags::MomentumDensity<DataVector, Dim>>,
      ::Tags::Source<NewtonianEuler::Tags::EnergyDensity<DataVector>>>;

  using argument_tags =
      tmpl::list<NewtonianEuler::Tags::MassDensityCons<DataVector>,
                 NewtonianEuler::Tags::MomentumDensity<DataVector, Dim>>;

  void apply(gsl::not_null<tnsr::I<DataVector, Dim>*> source_momentum_density,
             gsl::not_null<Scalar<DataVector>*> source_energy_density,
             const Scalar<DataVector>& mass_density_cons,
             const tnsr::I<DataVector, Dim>& momentum_density) const noexcept;

 private:
  std::array<double, Dim> acceleration_field_ =
      make_array<Dim>(std::numeric_limits<double>::signaling_NaN());
};
}  // namespace Sources
}  // namespace NewtonianEuler
