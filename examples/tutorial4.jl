using DataFrames
using LineCableModels
using LineCableModels.Engine.FEM
using EzXML
using Printf
fullfile(filename) = joinpath(@__DIR__, filename); #hide
set_logger!(0); #hide



"""
    create_equivalent_cable(original_design::CableDesign; new_id::String = "")

Builds a simplified `CableDesign` from an existing one by replacing each
component with an equivalent single-layer component.

This function iterates through each `CableComponent` of the `original_design`. 
For each component, it extracts the effective properties of its conductor and
insulator groups and uses them to create a simple `Tubular` conductor and
a single `Insulator` layer. The result is a new `CableDesign` that is
geometrically and electrically equivalent but structurally much simpler.

# Arguments
- `original_design`: The detailed `CableDesign` object to be simplified.
- `new_id`: An optional string for the new cable's ID. If not provided,
  it defaults to the original ID with an "_equivalent" suffix.

# Returns
- A new, simplified `CableDesign` object.
"""
function create_equivalent_cable(original_design::CableDesign; new_id::String = "")::CableDesign
    # Determine the ID for the new equivalent cable.
    equivalent_id = isempty(new_id) ? original_design.cable_id * "_equivalent" : new_id
    
    equivalent_design = nothing

    for (i, original_component) in enumerate(original_design.components)
        # --- 1. Build the Equivalent Conductor Group ---
        original_cond_group = original_component.conductor_group
        cond_props = original_component.conductor_props
        equivalent_cond_material = Material(
            cond_props.rho, cond_props.eps_r, cond_props.mu_r,
            cond_props.T0, cond_props.alpha
        )
        equivalent_conductor_part = Tubular(
            original_cond_group.radius_in,
            original_cond_group.radius_ext,
            equivalent_cond_material
        )
        new_cond_group = ConductorGroup(equivalent_conductor_part)

        # --- 2. Build the Equivalent Insulator Group ---
        original_ins_group = original_component.insulator_group
        ins_props = original_component.insulator_props
        equivalent_ins_material = Material(
            ins_props.rho, ins_props.eps_r, ins_props.mu_r,
            ins_props.T0, ins_props.alpha
        )
        insulator_thickness = original_ins_group.radius_ext - original_cond_group.radius_ext
        equivalent_insulator_part = Insulator(
            new_cond_group,
            Thickness(insulator_thickness),
            equivalent_ins_material
        )
        new_ins_group = InsulatorGroup(equivalent_insulator_part)

        if i == 1
            new_component = CableComponent(original_component.id, new_cond_group, new_ins_group)
            equivalent_design = CableDesign(
                equivalent_id,
                new_component,
                nominal_data=original_design.nominal_data
            )
        else
            add!(equivalent_design, original_component.id, new_cond_group, new_ins_group)
        end
    end

    return equivalent_design
end


"""
    rebuild_library(modification::NamedTuple; base_library=MaterialsLibrary(add_defaults=true))

Creates a new, modified `MaterialsLibrary` by applying changes to a single material.

This optimized function performs a deep copy of the `base_library` and directly
modifies the properties of the specified material, making it faster and more concise
than iterating through the entire collection.

# Arguments
- `modification`: A `NamedTuple` containing an `id::String` key for the material to
  change, plus any properties to override (e.g., `rho`, `eps_r`).
- `base_library`: The source `MaterialsLibrary`. Defaults to the package's default library.

# Returns
- A new `MaterialsLibrary` object with the modified material.
"""
function rebuild_library(modification::NamedTuple; base_library=MaterialsLibrary(add_defaults=true))::MaterialsLibrary

    if !haskey(modification, :id)
        error("Modification NamedTuple must contain an 'id' key.")
    end

    material_id = modification.id
    
    # faster than rebuilding.
    new_library = deepcopy(base_library)

    if !haskey(new_library.data, material_id)
        @warn "Material '$material_id' not found. Returning an unmodified library copy."
        return new_library
    end

    original_material = get(new_library, material_id)
    changes = Base.structdiff(modification, (id=nothing,))

    modified_material = Material(
        get(changes, :rho,   original_material.rho),
        get(changes, :eps_r, original_material.eps_r),
        get(changes, :mu_r,  original_material.mu_r),
        get(changes, :T0,    original_material.T0),
        get(changes, :alpha, original_material.alpha)
    )

    new_library.data[material_id] = modified_material

    return new_library
end


"""
    build_new_cable(input_set::NamedTuple) -> Tuple{Float64, Float64, Float64, CableDesign}

Builds and configures a `CableDesign` object based on the parameters provided in the `input_set`.

This function programmatically constructs a detailed cable model, including its conductive
and insulating layers. It also computes the fundamental electrical parameters (R, L, C)
for the primary conductor.

# Arguments
- `input_set::NamedTuple`: A named tuple containing all necessary parameters, such as
  layer dimensions, material properties, and component definitions.

# Returns
- `Tuple{BASE_FLOAT, BASE_FLOAT, BASE_FLOAT, CableDesign}`: A tuple containing:
    - `R`: The resistance of the core conductor (Ω/km).
    - `L`: The inductance of the core conductor (mH/km).
    - `C`: The capacitance of the core conductor (μF/km).
    - `cable_design`: The fully constructed `CableDesign` object.
"""
function build_new_cable(input_set::NamedTuple)::Tuple{Float64, Float64, Float64, CableDesign}

    materials = input_set.materials

    #--- Core Conductor ---
    material_core = get(materials, "copper")
    core = ConductorGroup(WireArray(0.0, Diameter(input_set.d_w), 1, 0.0, material_core))
    n_strands = 6 # Strands per layer
    for i in 1:input_set.n_layers
        add!(core, WireArray, Diameter(input_set.d_w), i * n_strands, 11.0, material_core)
    end

    #--- Main Insulation and Semiconductors ---
    material_sc1 = get(materials, "semicon1")
    main_insu = InsulatorGroup(Semicon(core, Thickness(input_set.t_sc_in), material_sc1))

    material_pe = get(materials, "polyethylene")
    add!(main_insu, Semicon, Thickness(input_set.t_ins), material_pe)

    material_sc2 = get(materials, "semicon2")
    add!(main_insu, Insulator, Thickness(input_set.t_sc_out), material_sc2)

    # Outer semiconductor
    material = get(materials, "polyacrylate")
    add!(main_insu, Semicon, Thickness(input_set.t_wbt), material)

    # Group the core components
    core_cc = CableComponent("core", core, main_insu)

    # Instantiate the CableDesign object
    cable_design = CableDesign(input_set.cable_id, core_cc, nominal_data=input_set.datasheet_info)

    # SHEATHS

    material = get(materials, "lead")
    screen_con = ConductorGroup(Tubular(main_insu, Thickness(input_set.t_sc), material))

    material_hdpe = get(materials, "high_density_pe")
    screen_insu = InsulatorGroup(Insulator(screen_con, Thickness(input_set.t_pe), material_hdpe))

    material_polyprop = get(materials, "polypropylene")
    add!(screen_insu, Insulator, Thickness(input_set.t_bed), material_polyprop)

    sheath_cc = CableComponent("sheath", screen_con, screen_insu)
    add!(cable_design, sheath_cc)

    # ARMOR
    lay_ratio = 10.0 # typical value for wire screens
    material_steel = get(materials, "steel")
    armor_con = ConductorGroup(
        WireArray(screen_insu, Diameter(input_set.d_wa), input_set.num_ar_wires, lay_ratio, material_steel))

    # PP layer after armor:
    armor_insu = InsulatorGroup(Insulator(armor_con, Thickness(input_set.t_jac), material_polyprop))

    # Assign the armor parts directly to the design:
    add!(cable_design, "armor", armor_con, armor_insu)

    core_df = DataFrame(cable_design, :baseparams)

    R = core_df[1,:computed]      # Ω/km
    L = core_df[2,:computed]      # mH/km
    C = core_df[3,:computed]      # μF/km

    return R, L, C, cable_design
end


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
cable_id = "Single_core_submarine_lead_sheath"
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
plt1 = preview(cable_design)

println("="^50)
println("Results of $(input_set.cable_id)")
println("-"^50)
println("R: ", round(R, digits=5), " Ω/km")
println("L: ", round(L, digits=4), " mH/km")
println("C: ", round(C, digits=4), " μF/km")
println("="^50)

# Results of base parameters with modification
custom_materials = rebuild_library((id="semicon1", rho=2000), base_library=default_lib)
modified_input_set = merge(input_set, (t_sc_in=2*0.6e-3, cable_id = cable_id*"_modified", materials=custom_materials))
R_mod, L_mod, C_mod, cable_design_mod = build_new_cable(modified_input_set)
plt2 = preview(cable_design_mod)

println("\nResults of $(modified_input_set.cable_id)")
println("-"^50)
println("R: ", round(R_mod, digits=5), " Ω/km")
println("L: ", round(L_mod, digits=4), " mH/km")
println("C: ", round(C_mod, digits=4), " μF/km")
println("="^50)

cable_design_mod.nominal_data = NominalData( # Dados de referência (opcional, mas útil para comparação)
    designation_code=cable_id*"_modified",
    U0=18.0, U=30.0, conductor_cross_section=1000.0, screen_cross_section=35.0,
    resistance=R_mod, capacitance=C_mod, inductance=L_mod)


# Obtain the equivalent electromagnetic properties of the cable:
components_df = DataFrame(cable_design, :components)


cable_design_eq = create_equivalent_cable(cable_design_mod)
plt3 = preview(cable_design_eq)

f = 1e-3 # Near DC frequency for the analysis
earth_params = EarthModel([f], 100.0, 10.0, 1.0)  # 100 Ω·m resistivity, εr=10, μr=1
# Define the coordinates for both cables:
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
output_file = fullfile("$(cable_system.system_id)_export.xml")
xml_path = export_data(Val(:atp), cable_system, earth_params, file_name=output_file)


# # Run the FEM model
# @time workspace, line_params = compute!(problem, formulation);
# println()
# println("="^50)
# println("Equivalent Cable Design FEM Results:")
# # Display primary core results
# if !opts.mesh_only
#     Z = line_params.Z[1, 1, 1]
#     Y = line_params.Y[1, 1, 1]
#     R = real(Z) * 1000
#     L = imag(Z) / (2π * f) * 1e6
#     C = imag(Y) / (2π * f) * 1e9
#     G = real(Y) * 1000
#     println("R = $(@sprintf("%.6g", R)) Ω/km")
#     println("L = $(@sprintf("%.6g", L)) mH/km")
#     println("C = $(@sprintf("%.6g", C)) μF/km/n")
#     println("G = $(@sprintf("%.6g", G)) Ω/km")
# end

# # # Z_self = line_params.Z[1, 1, 1]
# # # Y_self = line_params.Y[1, 1, 1]
# # # rho_eq = calc_equivalent_rho(real(Z_self), , 0.01)  # Expected output: ~9.42e-4 [Ω·m]
# # # println("="^50)
# # # println("Equivalent Cable Design Nominal Data:")	
# # # println("R = $(@sprintf("%.6g", cable_design_eq.nominal_data.resistance)) Ω/km")
# # # println("L = $(@sprintf("%.6g", cable_design_eq.nominal_data.inductance)) mH/km")
# # # println("C = $(@sprintf("%.6g", cable_design_eq.nominal_data.capacitance)) μF/km")

# # #  ok    ok  nope, nope, nope, nope,nope,nope  
# # # Rin, Rout, rho, muC, muI, epsI, Cext, Gext




# # # Z_matrix = randn(ComplexF64, 3, 3, length([f]))
# # # Y_matrix = randn(ComplexF64, 3, 3, length([f]))
# # # line_params = LineParameters(Z_matrix, Y_matrix)
# output_file = fullfile("$(cable_system.system_id)_ZY.xml")
# export_data(Val(:atp), line_params, [f];  file_name=output_file)
# # "C:\ATP\atpdraw\ATPDraw.exe" -runLCC "C:\\teste_atp_lcc\\Single_core_submarine_lead_sheath_equivalent_export_acp.acp"
# # "C:\ATP\tools\runATP.exe" "C:\teste_atp_lcc\Single_core_submarine_lead_sheath_equivalent_export_acp.acp"