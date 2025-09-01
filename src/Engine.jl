"""
    LineCableModels.Engine

The [`Engine`](@ref) module provides the main functionalities of the [`LineCableModels.jl`](index.md) package. This module implements data structures, methods and functions for calculating frequency-dependent electrical parameters (Z/Y matrices) of line and cable systems with uncertainty quantification. 

# Overview

- Calculation of frequency-dependent series impedance (Z) and shunt admittance (Y) matrices.
- Uncertainty propagation for geometric and material parameters using `Measurements.jl`.
- Internal impedance computation for solid, tubular and multi-layered coaxial conductors.
- Earth return impedances/admittances for overhead lines and underground cables (valid up to 10 MHz).
- Support for frequency-dependent soil properties.
- Handling of arbitrary polyphase systems with multiple conductors per phase.
- Phase and sequence domain calculations with uncertainty quantification.
- Novel N-layer concentric cable formulation with semiconductor modeling.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Engine

# Export public API
export LineParametersProblem, LineParameters, SeriesImpedance, ShuntAdmittance, per_km, per_m
export CoaxialFormulation

export compute!

# Module-specific dependencies
using Reexport, ForceImport
using Measurements
using LinearAlgebra
using SpecialFunctions
using ..Commons
import ..Commons: get_description, FormulationSet

using ..Utils
using ..Materials
using ..EarthProps: EarthModel
using ..DataModel: LineCableSystem

include("engine/types.jl")

# Problem definitions
include("engine/problemdefs.jl")
include("engine/lineparams.jl")

# Submodule `InternalImpedance`
include("engine/InternalImpedance.jl")
@force using .InternalImpedance

# Submodule `InsulationImpedance`
include("engine/InsulationImpedance.jl")
@force using .InsulationImpedance

# Submodule `EarthImpedance`
include("engine/EarthImpedance.jl")
@force using .EarthImpedance

# Submodule `InsulationAdmittance`
include("engine/InsulationAdmittance.jl")
@force using .InsulationAdmittance

# Submodule `EarthAdmittance`
include("engine/EarthAdmittance.jl")
@force using .EarthAdmittance

# Submodule `EHEM`
include("engine/EHEM.jl")
@force using .EHEM

# Helpers
include("engine/helpers.jl")

# Workspace definition
include("engine/workspace.jl")

# include all .jl files from src/legacy if the folder exists
isdir(joinpath(@__DIR__, "legacy")) && map(f -> endswith(f, ".jl") && include(joinpath(@__DIR__, "legacy", f)),
    sort(readdir(joinpath(@__DIR__, "legacy"))))

# Computation methods
include("engine/solver.jl")

# Override I/O methods
include("engine/base.jl")

# Submodule `FEM`
include("engine/FEM.jl")

@reexport using .InternalImpedance, .InsulationImpedance, .EarthImpedance,
    .InsulationAdmittance, .EarthAdmittance, .EHEM

end # module Engine