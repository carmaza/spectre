// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/NewtonianEuler/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace NewtonianEuler {

/*!
 * \brief Compute the source terms for the NewtonianEuler evolution
 * using a problem-specific source of type `SourceTermType`
 */
template <typename SourceTermType>
struct ComputeSources {
  using return_tags = typename SourceTermType::return_tags;

  using argument_tags = tmpl::push_front<typename SourceTermType::argument_tags,
                                         Tags::SourceTerm<SourceTermType>>;

  static constexpr size_t volume_dim = SourceTermType::volume_dim;

  template <class... Args>
  static void apply(
      gsl::not_null<tnsr::I<DataVector, volume_dim>*> source_momentum_density,
      gsl::not_null<Scalar<DataVector>*> source_energy_density,
      const db::item_type<Tags::SourceTerm<SourceTermType>>& source,
      Args&... args) noexcept {
    source.apply(source_momentum_density, source_energy_density, args...);
  }
};

}  // namespace NewtonianEuler
