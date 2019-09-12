// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/NewtonianEuler/TagsDeclarations.hpp"
#include "Time/Tags.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

// IWYU pragma: no_forward_declare Tensor

/// \cond
class DataVector;

namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl

namespace detail {
template <size_t Dim>
struct get_sourced_variables;

template <>
struct get_sourced_variables<2> {
  using type = tmpl::list<>;
};

template <>
struct get_sourced_variables<3> {
  using type = tmpl::list<NewtonianEuler::Tags::MassDensityCons<DataVector>,
                          NewtonianEuler::Tags::MomentumDensity<DataVector, 3>,
                          NewtonianEuler::Tags::EnergyDensity<DataVector>>;
};

template <size_t Dim>
struct get_argument_tags;

template <>
struct get_argument_tags<2> {
  using type = tmpl::list<>;
};

template <>
struct get_argument_tags<3> {
  using type =
      tmpl::list<::Tags::Coordinates<3, Frame::Inertial>, ::Tags::Time>;
};
}  // namespace detail
/// \endcond

namespace NewtonianEuler {
namespace Sources {

/*!
 * \brief Source generating a modified isentropic vortex.
 *
 * If Solutions::IsentropicVortex is modifed so that the flow velocity along
 * the \f$z-\f$axis is not a constant but a function of \f$z\f$, the new vortex
 * will be a solution to the 3-D Newtonian Euler equations with a source term,
 *
 * \f{align*}
 * \partial_t\rho + \partial_i F^i(\rho) &= S(\rho)\\
 * \partial_t S^i + \partial_j F^{j}(S^i) &= S(S^i)\\
 * \partial_t e + \partial_i F^i(e) &= S(e),
 * \f}
 *
 * where \f$F^i(u)\f$ is the volume flux of the conserved quantity \f$u\f$
 * (see ComputeFluxes), and
 *
 * \f{align*}
 * S(\rho) &= \rho \dfrac{dv_z}{dz}\\
 * S(S_x) &= S_x \dfrac{dv_z}{dz}\\
 * S(S_y) &= S_y \dfrac{dv_z}{dz}\\
 * S(S_z) &= 2S_z \dfrac{dv_z}{dz}\\
 * S(e) &= \left(e + p + v_z S_z\right)\dfrac{dv_z}{dz},
 * \f}
 *
 * where \f$\rho\f$ is the mass density of the vortex, \f$S_i\f$ is
 * its momentum density, \f$e\f$ is its energy density,
 * \f$v_z = v_z(z)\f$ is the \f$z-\f$component of its velocity,
 * and \f$p\f$ is its pressure. These quantities are readily obtained
 * from the primitive variables, whose expressions are those in
 * Solutions::IsentropicVortex
 */
template <size_t Dim>
struct VortexPerturbation {
  VortexPerturbation() noexcept = default;
  VortexPerturbation(const VortexPerturbation& /*rhs*/) = default;
  VortexPerturbation& operator=(const VortexPerturbation& /*rhs*/) = default;
  VortexPerturbation(VortexPerturbation&& /*rhs*/) noexcept = default;
  VortexPerturbation& operator=(VortexPerturbation&& /*rhs*/) noexcept =
      default;
  ~VortexPerturbation() = default;

  VortexPerturbation(double adiabatic_index, double perturbation_amplitude,
                     const std::array<double, Dim>& vortex_center,
                     const std::array<double, Dim>& vortex_mean_velocity,
                     double vortex_strength) noexcept;

  // clang-tidy: google-runtime-references
  void pup(PUP::er& /*p*/) noexcept;  // NOLINT

  using sourced_variables = typename detail::get_sourced_variables<Dim>::type;

  using argument_tags = typename detail::get_argument_tags<Dim>::type;

  // Overload required for 2-D simulations, where no variable is sourced.
  void apply() const noexcept;

  // Function to be used in 3-D.
  void apply(gsl::not_null<Scalar<DataVector>*> source_mass_density_cons,
             gsl::not_null<tnsr::I<DataVector, Dim>*> source_momentum_density,
             gsl::not_null<Scalar<DataVector>*> source_energy_density,
             const tnsr::I<DataVector, Dim>& x, double time) const noexcept;

 private:
  template <size_t SpatialDim>
  friend bool
  operator==(  // NOLINT (clang-tidy: readability-redundant-declaration)
      const VortexPerturbation<SpatialDim>& lhs,
      const VortexPerturbation<SpatialDim>& rhs) noexcept;

  double adiabatic_index_ = std::numeric_limits<double>::signaling_NaN();
  double perturbation_amplitude_ = std::numeric_limits<double>::signaling_NaN();
  std::array<double, Dim> vortex_center_ =
      make_array<Dim>(std::numeric_limits<double>::signaling_NaN());
  std::array<double, Dim> vortex_mean_velocity_ =
      make_array<Dim>(std::numeric_limits<double>::signaling_NaN());
  double vortex_strength_ = std::numeric_limits<double>::signaling_NaN();
};

template <size_t Dim>
bool operator!=(const VortexPerturbation<Dim>& lhs,
                const VortexPerturbation<Dim>& rhs) noexcept;
}  // namespace Sources
}  // namespace NewtonianEuler
