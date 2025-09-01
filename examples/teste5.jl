using DataFrames
using LineCableModels
using LineCableModels.Commons: COMPLEXSCALAR, REALSCALAR
using LineCableModels.Engine.FEM
using LineCableModels.Engine: LineParameters
using EzXML
using Revise
using Printf
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_logger!(0); #hide

include("C:/Users/Amauri/OneDrive/Documentos/UnB/Mestrado/LineCableModels.jl/examples/build_atp_file.jl")


# 1. Custom Materials Library
default_lib = MaterialsLibrary(add_defaults=false)

add!(default_lib, "copper",          Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393))
add!(default_lib, "semicon1",        Material(1000.0, 1000.0, 1.0, 20.0, 0.0))
add!(default_lib, "polyethylene",    Material(1.97e14, 2.3, 1.0, 20.0, 0.0))
add!(default_lib, "semicon2",        Material(500.0, 1000.0, 1.0, 20.0, 0.0))
add!(default_lib, "polyacrylate",    Material(5.3e3, 32.3, 1.0, 20.0, 0.0))
add!(default_lib, "lead",            Material(21.4e-8, 1.0, 0.999983, 20.0, 0.00400))
add!(default_lib, "polypropylene",   Material(1e15, 2.8, 1.0, 20.0, 0.0))
add!(default_lib, "high_density_pe", Material(1e16, 2.3, 0.99, 20.0, 0.0))
add!(default_lib, "steel",           Material(13.8e-8, 1.0, 300.0, 20.0, 0.00450))


# 2. Input set for the cable design
cable_id = "HVDC_525KV_2500mm2"
input_set = (
    # --- Identificação ---
    cable_id = cable_id,
    datasheet_info = NominalData( # Dados de referência (opcional, mas útil para comparação)
    designation_code=nothing,
    U0=0, U=525.0,
    conductor_cross_section=1000.0, screen_cross_section=35.0,
        resistance=0.0291, capacitance=0.39, inductance=0.3,
        ),
        
    # --- Geometry ---
    num_co_wires = 61,     # number of core wires
    num_ar_wires = 62,      # number of armor wires
    d_core = 57.8e-3,       # nominal core overall diameter
    d_w = 0.0031915382432114617,        # nominal strand diameter of the core

    t_sc_in = 1.5e-3,       # nominal internal semicon thickness 
    t_ins = 21.3e-3,        # nominal main insulation thickness
    t_sc_out = 1.4e-3,      # nominal external semicon thickness

    # Sheath
    t_wbt = .7e-3,          # nominal thickness of the water blocking tape
    t_sc = 3e-3,            # nominal lead screen thickness
    t_pe = 2.5e-3,          # nominal PE inner sheath thickness

    # Jacket
    t_bed = 0.6e-3,         # nominal thickness of the PP bedding
    d_wa = 5e-3,            # nominal armor wire diameter
    t_jac = 4e-3,           # nominal PP jacket thickness
    n_layers = 5,           # Number of conductor wire layers
    
    # --- Materiais ---
    materials = default_lib
)
    
# Results of base parameters
R, L, C, cable_design = build_new_cable(input_set)
# plt1 = preview(cable_design)


# Results of base parameters with modification
custom_materials = rebuild_library((id="semicon1", rho=2000), base_library=default_lib)
modified_input_set = merge(input_set, (t_sc_in=2*0.6e-3, cable_id = cable_id*"_modified", materials=custom_materials))
R_mod, L_mod, C_mod, cable_design_mod = build_new_cable(modified_input_set)

cable_design_eq = simplify(cable_design_mod)

f = 1e-3 # Near DC frequency for the analysis
earth_params = EarthModel([f], 100.0, 10.0, 1.0)  # 100 Ω·m resistivity, εr=10, μr=1
xp, xn, y0 = -0.5, 0.5, -1.0;

# Initialize the `LineCableSystem` with the first cable (phase A):
cablepos = CablePosition(cable_design_eq, xp, y0,
    Dict("core" => 1, "sheath" => 0, "armor" => 0))
cable_system = LineCableSystem(cable_id*"_equivalent", 1000.0, cablepos)

# Add remaining cables (phases B and C):
add!(cable_system, cable_design_eq, xn, y0,
    Dict("core" => 2, "sheath" => 0, "armor" => 0))

# plt4 = preview(cable_system, zoom_factor=0.15)

#=
## FEM calculations
=#

# Define a LineParametersProblem with the cable system and earth model
problem = LineParametersProblem(
    cable_system,
    temperature=20.0,  # Operating temperature
    earth_props=earth_params,
    frequencies=[f],   # Frequency for the analysis
);

# Estimate domain size based on skin depth in the earth
domain_radius = calc_domain_size(earth_params, [f]);

# Define custom mesh transitions around each cable
mesh_transition1 = MeshTransition(
    cable_system,
    [1, 2],
    r_min=0.6,
    r_length=0.25,
    mesh_factor_min=0.01 / (domain_radius / 5),
    mesh_factor_max=0.25 / (domain_radius / 5),
    n_regions=2)


# Define runtime options 
opts = (
    force_remesh=true,                # Force remeshing
    force_overwrite=true,             # Overwrite existing files
    plot_field_maps=false,            # Do not compute/ plot field maps
    mesh_only=false,                  # Preview the mesh
    save_path=fullfile("fem_output"), # Results directory
    keep_run_files=true,              # Archive files after each run
    verbosity=0,                      # Verbosity
);

# Define the FEM formulation with the specified parameters
formulation = FormulationSet(:FEM,
    impedance=Darwin(),
    admittance=Electrodynamics(), 
    domain_radius=domain_radius,
    domain_radius_inf=domain_radius * 1.25,
    elements_per_length_conductor=1,
    elements_per_length_insulator=2,
    elements_per_length_semicon=1,
    elements_per_length_interfaces=5,
    points_per_circumference=16,
    mesh_size_min=1e-6,
    mesh_size_max=domain_radius / 5,
    mesh_transitions=[mesh_transition1],
    mesh_size_default=domain_radius / 10,
    mesh_algorithm=5,
    mesh_max_retries=20,
    materials=custom_materials,
    options=opts,
);
# output_file = fullfile("$(cable_system.system_id)_export.xml")
# xml_path = export_data(Val(:atp), cable_system, earth_params, file_name=output_file)

output_file = fullfile("$(cable_system.system_id).dat")
file = generate_atp_file(cable_system, earth_params, file_name=output_file)
lis_content = run_dat_file(file)




Z, Y = read_data(Val(:atp), lis_content, cable_system)
# result = read_data(Val(:atp), lis_content, cable_system)
# Y = zeros(ComplexF64, size(Z))


# Run the FEM model
# @time workspace, line_params = compute!(problem, formulation);
