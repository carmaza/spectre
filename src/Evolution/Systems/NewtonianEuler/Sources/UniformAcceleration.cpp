// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/NewtonianEuler/Sources/UniformAcceleration.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/NewtonianEuler/Tags.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace NewtonianEuler {
namespace Sources {

template <size_t Dim>
UniformAcceleration<Dim>::UniformAcceleration(
    std::array<double, Dim> acceleration_field) noexcept
    : acceleration_field_(acceleration_field) {}

template <size_t Dim>
void UniformAcceleration<Dim>::pup(PUP::er& p) noexcept {
  p | acceleration_field_;
}

template <size_t Dim>
void UniformAcceleration<Dim>::apply(
    const gsl::not_null<tnsr::I<DataVector, Dim>*> source_momentum_density,
    const gsl::not_null<Scalar<DataVector>*> source_energy_density,
    const Scalar<DataVector>& mass_density_cons,
    const tnsr::I<DataVector, Dim>& momentum_density) const noexcept {
  get(*source_energy_density) = 0.0;
  for (size_t i = 0; i < Dim; ++i) {
    source_momentum_density->get(i) =
        get(mass_density_cons) * gsl::at(acceleration_field_, i);
    get(*source_energy_density) +=
        gsl::at(acceleration_field_, i) * momentum_density.get(i);
  }
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data) template struct UniformAcceleration<DIM(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

}  // namespace Sources
}  // namespace NewtonianEuler
/// \endcond
