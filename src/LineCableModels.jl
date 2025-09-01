module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, FormulationSet, set_logger!

# Materials:
export Material, MaterialsLibrary

# Data model (design + system):
export Thickness, Diameter, WireArray, Strip, Tubular, Semicon, Insulator
export ConductorGroup, InsulatorGroup
export CableComponent, CableDesign, NominalData
export CablesLibrary
export CablePosition, LineCableSystem
export trifoil_formation, flat_formation, preview, simplify

# Earth properties:
export EarthModel

# Engine:
export LineParametersProblem, compute!

# Import/Export:
export export_data, read_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions

# Submodule `Commons`
include("Commons.jl")
using .Commons: IMPORTS, EXPORTS, add!, FormulationSet

# Submodule `Utils`
include("Utils.jl")
using .Utils: set_logger!

# Submodule `Validation`
include("Validation.jl")

# Submodule `Materials`
include("Materials.jl")
using .Materials: Material, MaterialsLibrary

# Submodule `EarthProps`
include("EarthProps.jl")
using .EarthProps: EarthModel

# Submodule `DataModel`
include("DataModel.jl")
using .DataModel: Thickness, Diameter, WireArray, Strip, Tubular, Semicon, Insulator, ConductorGroup, InsulatorGroup, CableComponent, CableDesign, NominalData, CablesLibrary, CablePosition, LineCableSystem, trifoil_formation, flat_formation, preview, simplify

# Submodule `Engine`
include("Engine.jl")
using .Engine: LineParametersProblem, compute!

# Submodule `ImportExport`
include("ImportExport.jl")
using .ImportExport: export_data, load!, save

end