using DataFrames
using LineCableModels
using EzXML
using Printf
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
function create_equivalent_cable(original_design::CableDesign; new_id::String = "")
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
function rebuild_library(modification::NamedTuple; base_library=MaterialsLibrary(add_defaults=true))

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
    build_new_cable(input_set::NamedTuple)

Builds, configures, and analyzes a cable system based on the parameters
provided in the `input_set`.

This function encapsulates the tutorial's logic, allowing for easy variation
of geometric, material, and system configuration parameters.

Returns the base parameters of the cable core (R, L, C).
"""
function build_new_cable(input_set::NamedTuple)

    # 1. CORE AND MAIN INSULATION CONSTRUCTION

    materials = input_set.materials

    #--- Core Conductor ---
    material_core = get(materials, "aluminum")
    core = ConductorGroup(WireArray(0.0, Diameter(input_set.d_w), 1, 0.0, material_core))
    n_strands = 6 # Strands per layer
    for i in 1:input_set.n_layers
        add!(core, WireArray, Diameter(input_set.d_w), i * n_strands, 11.0, material_core)
    end

    #--- Main Insulation and Semiconductors ---
    # Inner semiconducting tape
    material_tape = get(materials, "polyacrylate")
    main_insu = InsulatorGroup(Semicon(core, Thickness(input_set.t_sct), material_tape))

    # Inner semiconductor
    material_sc1 = get(materials, "semicon1")
    add!(main_insu, Semicon, Thickness(input_set.t_sc_in), material_sc1)

    # Main insulation (XLPE)
    material_ins = get(materials, "pe")
    add!(main_insu, Insulator, Thickness(input_set.t_ins), material_ins)

    # Outer semiconductor
    material_sc2 = get(materials, "semicon2")
    add!(main_insu, Semicon, Thickness(input_set.t_sc_out), material_sc2)

    # Outer semiconducting tape
    add!(main_insu, Semicon, Thickness(input_set.t_sct), material_tape)

    # Group the core components
    core_cc = CableComponent("core", core, main_insu)

    # Instantiate the CableDesign object
    cable_design = CableDesign(input_set.cable_id, core_cc, nominal_data=input_set.datasheet_info)

    # 2. SHIELD AND SHEATHS

    #--- Copper Wire Screen ---
    material_screen = get(materials, "copper")
    screen_con = ConductorGroup(WireArray(main_insu, Diameter(input_set.d_ws), input_set.num_sc_wires, 10.0, material_screen))
    add!(screen_con, Strip, Thickness(input_set.t_cut), input_set.w_cut, 10.0, material_screen)

    # Water-blocking tape over the screen
    screen_insu = InsulatorGroup(Semicon(screen_con, Thickness(input_set.t_wbt), material_tape))

    add!(cable_design, "sheath", screen_con, screen_insu)

    #--- Outer Jacket ---
    # Aluminum tape (moisture barrier)
    jacket_con = ConductorGroup(Tubular(screen_insu, Thickness(input_set.t_alt), material_core)) # Reuse aluminum
    jacket_insu = InsulatorGroup(Insulator(jacket_con, Thickness(input_set.t_pet), material_ins)) # Reuse PE
    add!(jacket_insu, Insulator, Thickness(input_set.t_jac), material_ins) # Reuse PE

    # Add the jacket component
    add!(cable_design, "jacket", jacket_con, jacket_insu)

    core_df = DataFrame(cable_design, :baseparams)

    R = core_df[1,:computed]      # Ω/km
    L = core_df[2,:computed]      # mH/km
    C = core_df[3,:computed]      # μF/km

    return R, L, C, cable_design
end


# 1. Custom Materials Library
default_lib = MaterialsLibrary(add_defaults=true)


# 2. Input set for the cable design
input_set = (
    # --- Identificação ---
    cable_id = "18kV_1000mm2_Original",
    datasheet_info = NominalData( # Dados de referência (opcional, mas útil para comparação)
    designation_code="NA2XS(FL)2Y",
    U0=18.0, U=30.0,
    conductor_cross_section=1000.0, screen_cross_section=35.0,
        resistance=0.0291, capacitance=0.39, inductance=0.3,
        ),
        
    # --- Geometry ---
    num_co_wires = 61,      # Number of conductor strands
    num_sc_wires = 49,      # Number of shield strands
    d_w = 4.7e-3,           # Conductor wire diameter [m]
    d_ws = 0.95e-3,         # Shield wire diameter [m]
    t_sc_in = 0.6e-3,       # Inner semiconductor thickness [m]
    t_ins = 8.0e-3,         # Core insulation thickness [m]
    t_sc_out = 0.3e-3,      # Outer semiconductor thickness [m]
    t_cut = 0.1e-3,         # Copper tape thickness [m]
    w_cut = 10e-3,          # Copper tape width [m]
    t_wbt = 0.3e-3,         # Water-blocking tape thickness [m]
    t_sct = 0.3e-3,         # Semiconductor tape thickness [m]
    t_alt = 0.15e-3,        # Aluminum tape thickness [m]
    t_pet = 0.05e-3,        # PE face thickness on aluminum tape [m]
    t_jac = 2.4e-3,         # PE outer sheath thickness [m]
    n_layers = 6,           # Number of conductor wire layers
    
    # --- Materiais ---
    materials = default_lib
)
    
# Results of base parameters
R, L, C, cable_design = build_new_cable(input_set)

println("="^50)
println("Results of $(input_set.cable_id)")
println("-"^50)
println("R: ", round(R, digits=5), " Ω/km")
println("L: ", round(L, digits=4), " mH/km")
println("C: ", round(C, digits=4), " μF/km")
println("="^50)
    
# Results of base parameters with modification
custom_materials = rebuild_library((id="semicon1", rho=2000), base_library=default_lib)
modified_input_set = merge(input_set, (t_sc_in=2*0.6e-3, cable_id = "18kV_1000mm2_modified", materials=custom_materials))
R_mod, L_mod, C_mod, cable_design_mod = build_new_cable(modified_input_set)

println("/nResults of $(modified_input_set.cable_id)")
println("-"^50)
println("R: ", round(R_mod, digits=5), " Ω/km")
println("L: ", round(L_mod, digits=4), " mH/km")
println("C: ", round(C_mod, digits=4), " μF/km")
println("="^50)

cable_design_mod.nominal_data = NominalData( # Dados de referência (opcional, mas útil para comparação)
    designation_code="18kV_1000mm2_modified",
    U0=18.0, U=30.0, conductor_cross_section=1000.0, screen_cross_section=35.0,
    resistance=R_mod, capacitance=C_mod, inductance=L_mod)


# Obtain the equivalent electromagnetic properties of the cable:
components_df = DataFrame(cable_design, :components)

# plt3 = preview(cable_design)

cable_design_eq = create_equivalent_cable(cable_design_mod)
f = 1e-3
earth_params = EarthModel([f], 100.0, 10.0, 1.0)  # 100 Ω·m resistivity, εr=10, μr=1
Z = 1.343427906989891e-5 + 1.795286625550381e-8im
Y = 1.023542748680291e-12 + 3.587644237476014e-12im 
R = real(Z)
L = imag(Z) / (2π * f)
C = imag(Y) / (2π * f)
G = real(Y)
for (i, component) in enumerate(cable_design_eq.components)
    # --- 1. Build the Equivalent Conductor Group ---
    cond_group = component.conductor_group
    cond_props = component.conductor_props
    println("Rin [m] = $(cond_group.radius_in)")
    println("Rout [m] = $(cond_group.radius_ext)")	
    rho_eq = calc_equivalent_rho(R, cond_group.radius_ext, cond_group.radius_in)
    mu_r_cond = calc_equivalent_mu(cond_group.gmr, cond_group.radius_ext, cond_group.radius_in)
    println("Rho [ohm.m] = $(rho_eq)")
    println("mu = $(mu_r_cond)")
        
    # # --- 2. Build the Equivalent Insulator Group ---
    ins_group = component.insulator_group
    ins_props = component.insulator_props
    mu_r_ins = 1
    println("mu [ins] = $(mu_r_cond)")
    eps_eq = calc_equivalent_eps(C, ins_group.radius_in, ins_group.radius_ext)
    println("eps [ins] = $(eps_eq)")
    println("G [S/m] = $(G)")
    println("C [F/m] = $(C)")
    # gmr_eq = calc_equivalent_gmr(cond_group, ins_group) 
    println("-"^50)
end





