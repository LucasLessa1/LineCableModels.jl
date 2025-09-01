# -----------------------------------------------------------------------------
# SETUP: Load packages and define materials/dimensions from the tutorial.
# -----------------------------------------------------------------------------
using DataFrames
using LineCableModels


fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_logger!(0); #hide
# Initialize a materials library with default material properties
materials = MaterialsLibrary(add_defaults=true)

# Define the cable dimensions, taken directly from your tutorial code.
# This ensures all our cables are based on the same physical specs.
d_w = 4.7e-3      # nominal strand diameter of the core
t_sct = .3e-3     # nominal thickness of the semiconductive tape
t_sc_in = 0.6e-3  # nominal internal semicon thickness
t_ins = 8e-3      # nominal main insulation thickness
t_sc_out = 0.3e-3 # nominal external semicon thickness
d_ws = .95e-3     # nominal wire screen diameter
num_sc_wires = 49 # number of screen wires
t_cut = 0.1e-3    # nominal thickness of the copper tape
w_cut = 10e-3     # nominal width of copper tape
t_wbt = .3e-3     # nominal thickness of the water blocking tape
t_alt = .15e-3    # nominal thickness of the aluminum tape
t_pet = .05e-3    # nominal thickness of the pe face in the aluminum tape
t_jac = 2.4e-3    # nominal PE jacket thickness

# -----------------------------------------------------------------------------
# FUNCTION 1: Create a complete cable with Core, Sheath, and Armor.
# -----------------------------------------------------------------------------
function create_full_cable(materials::MaterialsLibrary, cable_id::String)
    # 1. CORE CONDUCTOR (Aluminum)
    # Built from 5 layers of stranded wires
    mat_aluminum = get(materials, "aluminum")
    core_conductor = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, mat_aluminum))
    add!(core_conductor, WireArray, Diameter(d_w), 6, 15.0, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 12, 13.5, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 18, 12.5, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 24, 11.0, mat_aluminum)

    # 2. CORE INSULATION (Semicon tapes, XLPE)
    # Wraps around the core conductor
    core_insulation = InsulatorGroup(Semicon(core_conductor, Thickness(t_sct), get(materials, "polyacrylate")))
    add!(core_insulation, Semicon, Thickness(t_sc_in), get(materials, "semicon1"))
    add!(core_insulation, Insulator, Thickness(t_ins), get(materials, "pe"))
    add!(core_insulation, Semicon, Thickness(t_sc_out), get(materials, "semicon2"))
    add!(core_insulation, Semicon, Thickness(t_sct), get(materials, "polyacrylate"))
    
    # 3. BUNDLE into a "core" component and initialize the CableDesign
    core_component = CableComponent("core", core_conductor, core_insulation)
    cable_design = CableDesign(cable_id, core_component)

    # 4. SHEATH CONDUCTOR (Copper Screen)
    # Built on top of the core's insulation
    mat_copper = get(materials, "copper")
    sheath_conductor = ConductorGroup(WireArray(core_insulation, Diameter(d_ws), num_sc_wires, 10.0, mat_copper))
    add!(sheath_conductor, Strip, Thickness(t_cut), w_cut, 10.0, mat_copper)
    
    # 5. SHEATH INSULATION (Water-blocking tape)
    sheath_insulation = InsulatorGroup(Semicon(sheath_conductor, Thickness(t_wbt), get(materials, "polyacrylate")))
    
    # 6. BUNDLE into a "sheath" component and add to the design
    sheath_component = CableComponent("sheath", sheath_conductor, sheath_insulation)
    add!(cable_design, sheath_component)

    # 7. ARMOR/JACKET CONDUCTOR (Aluminum Tape)
    # Built on top of the sheath's insulation
    armor_conductor = ConductorGroup(Tubular(sheath_insulation, Thickness(t_alt), mat_aluminum))
    
    # 8. ARMOR/JACKET INSULATION (PE Layers)
    armor_insulation = InsulatorGroup(Insulator(armor_conductor, Thickness(t_pet), get(materials, "pe")))
    add!(armor_insulation, Insulator, Thickness(t_jac), get(materials, "pe"))
    
    # 9. BUNDLE into an "armor" component and add to the design
    add!(cable_design, "armor", armor_conductor, armor_insulation)
    
    return cable_design
end

# -----------------------------------------------------------------------------
# FUNCTION 2: Create a cable with only Core and Sheath.
# -----------------------------------------------------------------------------
function create_core_sheath_cable(materials::MaterialsLibrary, cable_id::String)
    # Steps 1-6 are identical to the full cable
    mat_aluminum = get(materials, "aluminum")
    core_conductor = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, mat_aluminum))
    add!(core_conductor, WireArray, Diameter(d_w), 6, 15.0, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 12, 13.5, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 18, 12.5, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 24, 11.0, mat_aluminum)

    core_insulation = InsulatorGroup(Semicon(core_conductor, Thickness(t_sct), get(materials, "polyacrylate")))
    add!(core_insulation, Semicon, Thickness(t_sc_in), get(materials, "semicon1"))
    add!(core_insulation, Insulator, Thickness(t_ins), get(materials, "pe"))
    add!(core_insulation, Semicon, Thickness(t_sc_out), get(materials, "semicon2"))
    add!(core_insulation, Semicon, Thickness(t_sct), get(materials, "polyacrylate"))
    
    core_component = CableComponent("core", core_conductor, core_insulation)
    cable_design = CableDesign(cable_id, core_component)

    mat_copper = get(materials, "copper")
    sheath_conductor = ConductorGroup(WireArray(core_insulation, Diameter(d_ws), num_sc_wires, 10.0, mat_copper))
    add!(sheath_conductor, Strip, Thickness(t_cut), w_cut, 10.0, mat_copper)
    
    sheath_insulation = InsulatorGroup(Semicon(sheath_conductor, Thickness(t_wbt), get(materials, "polyacrylate")))
    
    sheath_component = CableComponent("sheath", sheath_conductor, sheath_insulation)
    add!(cable_design, sheath_component)
    
    # We stop here. No armor is added.
    return cable_design
end

# -----------------------------------------------------------------------------
# FUNCTION 3: Create the simplest cable with only a Core.
# -----------------------------------------------------------------------------
function create_core_only_cable(materials::MaterialsLibrary, cable_id::String)
    # This function only performs the first 3 steps.
    mat_aluminum = get(materials, "aluminum")
    core_conductor = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, mat_aluminum))
    add!(core_conductor, WireArray, Diameter(d_w), 6, 15.0, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 12, 13.5, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 18, 12.5, mat_aluminum)
    add!(core_conductor, WireArray, Diameter(d_w), 24, 11.0, mat_aluminum)

    core_insulation = InsulatorGroup(Semicon(core_conductor, Thickness(t_sct), get(materials, "polyacrylate")))
    add!(core_insulation, Semicon, Thickness(t_sc_in), get(materials, "semicon1"))
    add!(core_insulation, Insulator, Thickness(t_ins), get(materials, "pe"))
    add!(core_insulation, Semicon, Thickness(t_sc_out), get(materials, "semicon2"))
    add!(core_insulation, Semicon, Thickness(t_sct), get(materials, "polyacrylate"))
    
    core_component = CableComponent("core", core_conductor, core_insulation)
    cable_design = CableDesign(cable_id, core_component)
    
    # We stop here. No sheath or armor.
    return cable_design
end

# -----------------------------------------------------------------------------
# MAIN EXECUTION: Create the cables and store them in a library.
# -----------------------------------------------------------------------------

# Create the three different cable designs by calling our functions
cable1_full = create_full_cable(materials, "MV_Cable_Full")
cable2_core_sheath = create_core_sheath_cable(materials, "MV_Cable_CoreSheath")
cable3_core_only = create_core_only_cable(materials, "MV_Cable_CoreOnly")

xp, xn, y0 = -0.5, 0.5, -1.0;
cablepos = CablePosition(cable1_full, -0.5, -1.0,
    Dict("core" => 1, "sheath" => 0, "armor" => 0))
cable_system = LineCableSystem("multipleComp", 1000.0, cablepos)

add!(cable_system, cable2_core_sheath, -1.0, -1.0,
    Dict("core" => 2, "sheath" => 0))

add!(cable_system, cable3_core_only, -1.5, -1.0,
    Dict("core" => 3))