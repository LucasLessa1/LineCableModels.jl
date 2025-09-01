using DataFrames
using LineCableModels
using Printf
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_logger!(0); #hide


# Initialize library and the required materials for this design:
materials = MaterialsLibrary(add_defaults=true)
lead = Material(21.4e-8, 1.0, 0.999983, 20.0, 0.00400) # Lead or lead alloy
add!(materials, "lead", lead)
steel = Material(13.8e-8, 1.0, 300.0, 20.0, 0.00450) # Steel
add!(materials, "steel", steel)
pp = Material(1e15, 2.8, 1.0, 20.0, 0.0) # Laminated paper propylene
add!(materials, "pp", pp)
pe = Material(1e16, 2.3, 0.99, 20.0, 0.0) # high density polyethylene
add!(materials, "HDPE", pe)

# Inspect the contents of the materials library:
materials_df = DataFrame(materials)

## Cable dimensions
# Core
num_co_wires = 61 # number of core wires
num_ar_wires = 62  # number of armor wires
d_core = 57.8e-3   # nominal core overall diameter
d_w = 0.0031915382432114617# 4mm² 12AWG    # nominal strand diameter of the core (minimum value to match datasheet)
n_layers = 5 # Layers of strands

t_sc_in = 1.5e-3     # nominal internal semicon thickness 
t_ins = 21.3e-3      # nominal main insulation thickness
t_sc_out = 1.4e-3  # nominal external semicon thickness

# Sheath
t_wbt = .7e-3      # nominal thickness of the water blocking tape
t_sc = 3e-3      # nominal lead screen thickness
t_pe = 2.5e-3        # nominal PE inner sheath thickness

# Jacket
t_bed = 0.6e-3       # nominal thickness of the PP bedding
d_wa = 5e-3    # nominal armor wire diameter
t_jac = 4e-3      # nominal PP jacket thickness

d_overall = d_core #hide



layers = [] #hide
push!(layers, ("Conductor", missing, d_overall * 1000)) #hide
d_overall += 2 * t_sc_in #hide
push!(layers, ("Inner semiconductor", t_sc_in * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_ins #hide
push!(layers, ("Main insulation", t_ins * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sc_out #hide
push!(layers, ("Outer semiconductor", t_sc_out * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_wbt #hide
push!(layers, ("Swellable tape", t_wbt * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_sc #hide
push!(layers, ("Lead screen", t_sc * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_pe #hide
push!(layers, ("PE inner sheath", t_pe * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_bed #hide
push!(layers, ("PP bedding", t_bed * 1000, d_overall * 1000)) #hide
d_overall += 2 * d_wa #hide
push!(layers, ("Stranded wire armor", d_wa * 1000, d_overall * 1000)) #hide
d_overall += 2 * t_jac #hide
push!(layers, ("PP jacket", t_jac * 1000, d_overall * 1000)); #hide


# The cable structure is summarized in a table for better visualization, with dimensions in milimiters:
df = DataFrame( #hide
    layer=first.(layers), #hide
    thickness=[ #hide
        ismissing(t) ? "-" : round(t, sigdigits=2) for t in getindex.(layers, 2) #hide
    ], #hide
    diameter=[round(d, digits=2) for d in getindex.(layers, 3)], #hide
) #hide

#=
## Core and main insulation

Initialize the conductor object and assign the central wire:
=#

material = get(materials, "copper")
core = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, material))

# Add the subsequent layers of wires and inspect the object:
n_strands = 6 # Strands per layer

for i in 1:n_layers
    add!(core, WireArray, Diameter(d_w), i * n_strands, 11.0, material)
end
core

#=
### Inner semiconductor

Inner semiconductor (1000 Ω.m as per IEC 840):
=#

material = get(materials, "semicon1")
main_insu = InsulatorGroup(Semicon(core, Thickness(t_sc_in), material))

#=
### Main insulation

Add the insulation layer:
=#

material = get(materials, "pe")
add!(main_insu, Insulator, Thickness(t_ins), material)

#=
### Outer semiconductor

Outer semiconductor (500 Ω.m as per IEC 840):
=#

material = get(materials, "semicon2")
add!(main_insu, Semicon, Thickness(t_sc_out), material)

# Water blocking (swellable) tape:
material = get(materials, "polyacrylate")
add!(main_insu, Semicon, Thickness(t_wbt), material)

# Group core-related components:
core_cc = CableComponent("core", core, main_insu)

cable_id = "525kV_2500mm2"
datasheet_info = NominalData(
    designation_code="SingleCoreSubmarineLeadSheath",
    U0=500.0,                        # Phase (pole)-to-ground voltage [kV]
    U=525.0,                         # Phase (pole)-to-phase (pole) voltage [kV]
    conductor_cross_section=2500.0,  # [mm²]
    screen_cross_section=1000.0,     # [mm²]
    resistance=nothing,              # DC resistance [Ω/km]
    capacitance=nothing,             # Capacitance [μF/km]
    inductance=nothing,              # Inductance in trifoil [mH/km]
)
cable_design = CableDesign(cable_id, core_cc, nominal_data=datasheet_info)

#=
### Lead screen/sheath

Build the wire screens on top of the previous layer:
=#

material = get(materials, "lead")
screen_con = ConductorGroup(Tubular(main_insu, Thickness(t_sc), material))

# PE inner sheath:
material = get(materials, "HDPE")
screen_insu = InsulatorGroup(Insulator(screen_con, Thickness(t_pe), material))

# PP bedding:
material = get(materials, "pp")
add!(screen_insu, Insulator, Thickness(t_bed), material)

# Group sheath components and assign to design:
sheath_cc = CableComponent("sheath", screen_con, screen_insu)
add!(cable_design, sheath_cc)

#=
### Armor and outer jacket components

=#

# Add the armor wires on top of the previous layer:
lay_ratio = 10.0 # typical value for wire screens
material = get(materials, "steel")
armor_con = ConductorGroup(
    WireArray(screen_insu, Diameter(d_wa), num_ar_wires, lay_ratio, material))

# PP layer after armor:
material = get(materials, "pp")
armor_insu = InsulatorGroup(Insulator(armor_con, Thickness(t_jac), material))

# Assign the armor parts directly to the design:
add!(cable_design, "armor", armor_con, armor_insu)

# Inspect the finished cable design:
plt3 = preview(cable_design)

#=
## Examining the cable parameters (RLC)

=#

# Summarize DC lumped parameters (R, L, C):
core_df = DataFrame(cable_design, :baseparams)

# Obtain the equivalent electromagnetic properties of the cable:
components_df = DataFrame(cable_design, :components)

#=
## Saving the cable design

Load an existing [`CablesLibrary`](@ref) file or create a new one:
=#


library = CablesLibrary()
library_file = fullfile("cables_library.json")
load!(library, file_name=library_file)
add!(library, cable_design)
library_df = DataFrame(library)

# Save to file for later use:
save(library, file_name=library_file);

#=
## Defining a cable system

=#

#=
### Earth model 

Define a constant frequency earth model:
=#

f = 1e-3 # Near DC frequency for the analysis
earth_params = EarthModel([f], 100.0, 10.0, 1.0)  # 100 Ω·m resistivity, εr=10, μr=1

# Earth model base (DC) properties:
earthmodel_df = DataFrame(earth_params)

#=
### Underground bipole configuration

=#

# Define the coordinates for both cables:
xp, xn, y0 = -0.5, 0.5, -1.0;

# Initialize the `LineCableSystem` with positive pole:
cablepos = CablePosition(cable_design, xp, y0,
    Dict("core" => 1, "sheath" => 0, "armor" => 0))
cable_system = LineCableSystem("525kV_1600mm2_bipole", 1000.0, cablepos)

# Add the other pole (negative) to the system:
add!(cable_system, cable_design, xn, y0,
    Dict("core" => 2, "sheath" => 0, "armor" => 0))

#=
### Cable system preview

In this section the complete bipole cable system is examined.
=#

# Display system details:
system_df = DataFrame(cable_system)

# Visualize the cross-section of the three-phase system:
# plt4 = preview(cable_system, zoom_factor=0.15)
