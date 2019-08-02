// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/NewtonianEuler/TagsDeclarations.hpp"
#include "Evolution/Tags.hpp"
#include "Options/Options.hpp"

namespace NewtonianEuler {
/// %Tags for the conservative formulation of the Newtonian Euler system
namespace Tags {

/// The characteristic speeds.
template <size_t Dim>
struct CharacteristicSpeeds : db::SimpleTag {
  using type = std::array<DataVector, Dim + 2>;
  static std::string name() noexcept { return "CharacteristicSpeeds"; }
};

/// The mass density of the fluid.
template <typename DataType>
struct MassDensity : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "MassDensity"; }
};

/// The mass density of the fluid (as a conservative variable).
template <typename DataType>
struct MassDensityCons : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "MassDensityCons"; }
};

/// The momentum density of the fluid.
template <typename DataType, size_t Dim, typename Fr>
struct MomentumDensity : db::SimpleTag {
  using type = tnsr::I<DataType, Dim, Fr>;
  static std::string name() noexcept {
    return Frame::prefix<Fr>() + "MomentumDensity";
  }
};

/// The energy density of the fluid.
template <typename DataType>
struct EnergyDensity : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "EnergyDensity"; }
};

/// The macroscopic or flow velocity of the fluid.
template <typename DataType, size_t Dim, typename Fr>
struct Velocity : db::SimpleTag {
  using type = tnsr::I<DataType, Dim, Fr>;
  static std::string name() noexcept {
    return Frame::prefix<Fr>() + "Velocity";
  }
};

/// The specific internal energy of the fluid.
template <typename DataType>
struct SpecificInternalEnergy : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "SpecificInternalEnergy"; }
};

/// The fluid pressure.
template <typename DataType>
struct Pressure : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "Pressure"; }
};

/// The sound speed.
template <typename DataType>
struct SoundSpeed : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "SoundSpeed"; }
};

/// The square of the sound speed.
template <typename DataType>
struct SoundSpeedSquared : db::SimpleTag {
  using type = Scalar<DataType>;
  static std::string name() noexcept { return "SoundSpeedSquared"; }
};

/// The source term in the evolution equations
template <typename SourceTermType>
struct SourceTerm : db::SimpleTag {
  using type = SourceTermType;
  static std::string name() noexcept { return "SourceTerm"; }
};

}  // namespace Tags

namespace OptionTags {
/// \ingroup OptionGroupsGroup
/// Groups option tags related to the NewtonianEuler evolution system.
struct NewtonianEulerGroup {
  static std::string name() noexcept { return "NewtonianEuler"; }
  static constexpr OptionString help{"Options for the evolution system"};
  using group = ::OptionTags::EvolutionSystemGroup;
};

/// \ingroup OptionGroupsGroup
/// \brief Holds the `OptionTags::SourceTerm` option in the input file
struct SourceTermGroup {
  static std::string name() noexcept { return "SourceTerm"; }
  static constexpr OptionString help{"The source term used in the evolution."};
  using group = NewtonianEulerGroup;
};

/// \brief Source terms for the evolution equations
template <typename SourceTermType>
struct SourceTerm : Tags::SourceTerm<SourceTermType> {
  static std::string name() noexcept { return option_name<SourceTermType>(); }
  static constexpr OptionString help = SourceTermType::help;
  using type = SourceTermType;
  using group = SourceTermGroup;
};

}  // namespace OptionTags
}  // namespace NewtonianEuler
