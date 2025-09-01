module CableDatabase

# Import required libraries
using Measurements
using DataFrames
using Base.MathConstants: pi

# Export database functions
export create_materials_database, add_material_to_database!, get_material_property
export create_sc_cable_model, add_part_to_cable!, create_cables_database, add_cable_to_database!, create_system_cross_section, add_cable_to_cross_section!, create_cable_system, define_cable_system!
export convert_resistance_to_rho, convert_stranded_to_tubular, convert_capacitance_to_epsilon, correct_epsilon_for_semicon, apply_temperature_correction
export calc_equivalent_sheath_radius, calc_stranded_conductor_area
export build_geometry_matrix
export trifoil_formation, flat_formation

"""
    build_geometry_matrix(materials_df, cables_df, cable_system_df)

Constructs the geometry matrix for a given cable system by assembling the key physical properties of 
cable components (core, sheath, armor) including radii, material resistivity, permeability, and permittivity. 
The function also applies corrections for temperature and incorporates RDC and capacitance data if provided.

# Arguments
- `materials_df::DataFrame`: A DataFrame containing material properties, such as resistivity, relative permeability, 
  temperature coefficient, and reference temperature.
- `cables_df::DataFrame`: A DataFrame containing the cable configurations, including information about parts 
  (core, sheath, armor), nominal resistance (RDC), capacitance (C), and other cross-sectional properties.
- `cable_system_df::DataFrame`: A DataFrame representing the cable system, including the layout of cables in 
  the system cross-section, phase numbering, grounding conditions for sheaths and armors, and the operating 
  temperature (`T`).

# Returns
- `Geom::Matrix{Measurement{Float64}}`: A matrix where each row represents a distinct part of the cable (core, sheath, 
  or armor), with the following columns:
  - `cabID`: Integer identifier for the cable.
  - `phID`: Phase ID; 0 for grounded parts, or a unique identifier for ungrounded parts.
  - `horz`: Horizontal position of the cable core in meters (NaN for sheath/armor).
  - `vert`: Vertical position of the cable core in meters (NaN for sheath/armor).
  - `r_in`: Internal radius of the conductor/part.
  - `r_ext`: External radius of the conductor/part.
  - `rho_cond`: Resistivity of the conductor (corrected for temperature if applicable).
  - `mu_cond`: Relative permeability of the conductor material.
  - `r_ins`: External radius including insulation and semiconductor layers.
  - `mu_ins`: Relative permeability of the insulation material.
  - `eps_ins`: Relative permittivity of the insulation material (corrected for capacitance if provided).

# Notes
- The function iterates over the parts of each cable, checks for NaN values, and applies corrections based 
  on material properties and operating conditions (e.g., temperature, RDC, and capacitance).
- The phase ID is incremented for each ungrounded core, sheath, or armor. Grounded parts are assigned phase ID 0.
- Temperature corrections are applied using the formula `R = R0 * (1 + alpha * (T - T0))` if the operating temperature is known.
- Capacitance corrections are applied using the function `convert_capacitance_to_epsilon`.

# Example
```julia
materials_df = create_materials_database()
cables_df = create_cables_database()
cable_system_df = create_cable_system()

# Build the geometry matrix
Geom = build_geometry_matrix(materials_df, cables_df, cable_system_df)

# Category: Cable and materials database

"""
function build_geometry_matrix(materials_df, cables_df, cable_system_df)
	# First pass: count valid parts to preallocate matrix
	total_rows = 0
	for (_, row) in enumerate(eachrow(cable_system_df))
		cross_section = row.cross_section
		for cable_row in eachrow(cross_section)
			cable_id = cable_row.cable_id
			cable_info = first(cables_df[cables_df.cable_id.==cable_id, :])
			parts = cable_info.cable_parts

			# Count parts that exist (non-NaN radius)
			core_row = findfirst(parts.part .== "core")
			total_rows += 1  # Core always exists

			sheath_row = findfirst(parts.part .== "sheath")
			if !isnothing(sheath_row) && !isnan(parts.radius_in[sheath_row])
				total_rows += 1
			end

			armor_row = findfirst(parts.part .== "armor")
			if !isnothing(armor_row) && !isnan(parts.radius_in[armor_row])
				total_rows += 1
			end
		end
	end

	# Initialize the geometry matrix
	Geom = Matrix{Any}(undef, total_rows, 11)
	current_row = 1
	phase_counter = 1

	# Second pass: fill the matrix
	for (sys_idx, system) in enumerate(eachrow(cable_system_df))

		# Define operating temperature for corrections
		T_system = cable_system_df[sys_idx, :T]

		# Define cross-section and calculate params
		cross_section = system.cross_section

		for (cable_idx, cable_pos) in enumerate(eachrow(cross_section))
			cable_id = cable_pos.cable_id
			horz = cable_pos.horz
			vert = cable_pos.vert
			ground_sheath = cable_pos.ground_sheath
			ground_armor = cable_pos.ground_armor

			# Get cable specifications
			cable_info = first(cables_df[cables_df.cable_id.==cable_id, :])
			parts = cable_info.cable_parts

			# Handle core (always present)
			core_row = findfirst(parts.part .== "core")
			rho_core = get_material_property(materials_df, parts.material_cond[core_row], :rho)

			# Check if a valid RDC is provided in cable_info
			if !isnan(cable_info.RDC)
				# Convert RDC to rho using the predefined function
				rho_core = convert_resistance_to_rho(cable_info.RDC, parts.radius_ext[core_row]; r_in = parts.radius_in[core_row])
			end

			if !isnan(T_system)
				# Get material temperature coefficient and reference temperature from materials_df
				alpha_core = get_material_property(materials_df, parts.material_cond[core_row], :alpha)
				T0_core = get_material_property(materials_df, parts.material_cond[core_row], :T0)

				# Apply temperature correction: R = R0 * (1 + alpha * (T - T0))
				rho_core = apply_temperature_correction(rho_core, alpha_core, T_system, T0_core)
			end

			# Handle main insulation
			epsr_main_ins = get_material_property(materials_df, parts.material_ins[core_row], :epsr_r)


			r_ins_core = parts.radius_ext[core_row] + parts.width_ins[core_row] + sum(x -> !isnan(x) ? x : 0.0, [
							 parts.width_semicon_in[core_row],
							 parts.width_semicon_ext[core_row],
						 ]) # width_ins will override the sum if NaN. I have spoken.

			# Check if a valid C is provided in cable_info
			if !isnan(cable_info.C)
				# Convert C to epsr using the predefined function
				epsr_main_ins = convert_capacitance_to_epsilon(cable_info.C, r_ins_core; r_in = parts.radius_ext[core_row])
			end


			Geom[current_row, :] = [
				cable_idx,                    # cabID
				phase_counter,                # phID for core
				horz, vert,                   # core geometry
				parts.radius_in[core_row],
				parts.radius_ext[core_row],
				rho_core,
				get_material_property(materials_df, parts.material_cond[core_row], :mu_r),
				r_ins_core,
				get_material_property(materials_df, parts.material_ins[core_row], :mu_r),
				epsr_main_ins,
			]
			current_row += 1
			phase_counter += 1

			# Handle sheath
			sheath_row = findfirst(parts.part .== "sheath")
			rho_sheath = get_material_property(materials_df, parts.material_cond[sheath_row], :rho)

			radius_in_screen = r_ins_core #parts.radius_in[sheath_row]

			if !isnan(T_system)
				# Get material temperature coefficient and reference temperature from materials_df
				alpha_sheath = get_material_property(materials_df, parts.material_cond[sheath_row], :alpha)
				T0_sheath = get_material_property(materials_df, parts.material_cond[sheath_row], :T0)

				# Apply temperature correction: R = R0 * (1 + alpha * (T - T0))
				rho_sheath = apply_temperature_correction(rho_sheath, alpha_sheath, T_system, T0_sheath)
			end

			# Check if a valid screen_cross_section is provided in cable_info
			radius_ext_screen = parts.radius_ext[sheath_row]

			if !isnan(cable_info.screen_cross_section)
				# Adjust external radius to match the screen_cross_section
				radius_ext_screen = calc_equivalent_sheath_radius(cable_info.screen_cross_section, radius_in_screen)
			end

			r_ins_sheath = radius_ext_screen + parts.width_ins[sheath_row] + sum(x -> !isnan(x) ? x : 0.0, [
				parts.width_semicon_in[sheath_row],
				parts.width_semicon_ext[sheath_row],
			])

			if !isnothing(sheath_row) && !isnan(parts.radius_in[sheath_row])
				Geom[current_row, :] = [
					cable_idx,
					ground_sheath ? 0 : phase_counter,
					NaN, NaN,
					radius_in_screen,
					radius_ext_screen,
					rho_sheath,
					get_material_property(materials_df, parts.material_cond[sheath_row], :mu_r),
					r_ins_sheath,
					get_material_property(materials_df, parts.material_ins[sheath_row], :mu_r),
					epsr_main_ins,
				]
				current_row += 1
				phase_counter += ground_sheath ? 0 : 1
			end

			# Handle armor
			armor_row = findfirst(parts.part .== "armor")
			rho_armor = get_material_property(materials_df, parts.material_cond[armor_row], :rho) #TESTME
			
			radius_in_armor = r_ins_sheath #parts.radius_in[armor_row]

			if !isnan(T_system)
				# Get material temperature coefficient and reference temperature from materials_df
				alpha_armor = get_material_property(materials_df, parts.material_cond[armor_row], :alpha)
				T0_armor = get_material_property(materials_df, parts.material_cond[armor_row], :T0)

				# Apply temperature correction: R = R0 * (1 + alpha * (T - T0))
				rho_armor = apply_temperature_correction(rho_armor, alpha_armor, T_system, T0_armor)
			end

			# Check if a valid screen_cross_section is provided in cable_info
			radius_ext_armor = parts.radius_ext[armor_row]

			if !isnan(cable_info.armor_cross_section)
				# Adjust external radius to match the screen_cross_section
				radius_ext_armor = calc_equivalent_sheath_radius(cable_info.armor_cross_section, radius_in_armor)
			end

			r_ins_armor = radius_ext_armor + parts.width_ins[armor_row] + sum(x -> !isnan(x) ? x : 0.0, [
				parts.width_semicon_in[armor_row],
				parts.width_semicon_ext[armor_row],
			])

			if !isnothing(armor_row) && !isnan(parts.radius_in[armor_row])
				Geom[current_row, :] = [
					cable_idx,
					ground_armor ? 0 : phase_counter,
					NaN, NaN,
					radius_in_armor,
					radius_ext_armor,
					rho_armor,
					get_material_property(materials_df, parts.material_cond[armor_row], :mu_r),
					r_ins_armor,
					get_material_property(materials_df, parts.material_ins[armor_row], :mu_r),
					epsr_main_ins,
				]
				current_row += 1
				phase_counter += ground_armor ? 0 : 1
			end
		end
	end

	return Measurement{Float64}.(Geom)

end

"""
	create_materials_database()

Creates a `DataFrame` containing material properties for various conductive and insulating
materials, with support for uncertainties in resistivity, permittivity, and permeability.
The materials included are copper, aluminum, and XLPE (cross-linked polyethylene). The
DataFrame stores the following properties for each material:

# Returns
- `materials_df::DataFrame`: A DataFrame with the following columns:
	- `name::String`: The name of the material.
	- `rho::Measurement{Float64}`: Electrical resistivity of the material in ohm⋅m, with uncertainty.
	- `epsr_r::Measurement{Float64}`: Relative permittivity (dielectric constant) of the material, with uncertainty.
	- `mu_r::Measurement{Float64}`: Relative permeability of the material, with uncertainty.
	- `T0::Float64`: Reference temperature in °C (typically used for calculating resistivity at different temperatures).
	- `alpha::Float64`: Temperature coefficient of resistivity (1/K), used for adjusting resistivity with temperature changes.

# Category: Cable and materials database

"""
function create_materials_database()
	# Create a DataFrame with columns for material properties
	materials_df = DataFrame(
		name = String[],                      # Name of the material
		rho = Measurement{Float64}[],         # Resistivity with uncertainty in ohm*m
		epsr_r = Measurement{Float64}[],      # Relative permittivity with uncertainty
		mu_r = Measurement{Float64}[],        # Relative permeability with uncertainty
		T0 = Float64[],                       # Reference temperature in Celsius
		alpha = Float64[],                    # Temperature coefficient for resistivity corrections
	)

	# Populate the DataFrame with copper and aluminum properties
	# Values based on IEC 28 and 889 at T0 = 20ºC
	push!(materials_df, ("annealed_copper", 1.7241e-8, 1.0, 0.999994, 20.0, 0.00393))     # Copper α ≈ 0.00393 1/K
	push!(materials_df, ("aluminum", 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429))   # Aluminum α ≈ 0.00429 1/K
	push!(materials_df, ("xlpe", 1.97e14, 2.3, 1.0, 30, 0))

	return materials_df
end

"""
	add_new_material!(materials_df, name, rho, epsr_r, mu_r, T0, alpha)

Adds a new material to the provided materials `DataFrame` with its key electrical
and physical properties.

# Arguments:
- `materials_df`: The DataFrame containing the materials data.
- `name`: Name of the material (String).
- `rho`: Resistivity with uncertainty in ohm*m (Measurement{Float64}).
- `epsr_r`: Relative permittivity with uncertainty (Measurement{Float64}).
- `mu_r`: Relative permeability with uncertainty (Measurement{Float64}).
- `T0`: Reference temperature in Celsius (Float64).
- `alpha`: Temperature coefficient for resistivity corrections (Float64).

# Category: Cable and materials database

"""
function add_material_to_database!(materials_df::DataFrame, name::String, rho::Measurement{Float64}, epsr_r::Measurement{Float64}, mu_r::Measurement{Float64}, T0::Float64, alpha::Float64)
	# Add the new material as a new row to the DataFrame
	push!(materials_df, (name, rho, epsr_r, mu_r, T0, alpha))
end

"""
	get_material_property(materials_df, material_name, property_symbol)

Returns the electrical property described in the provided materials `DataFrame`.

# Category: Cable and materials database

"""
function get_material_property(materials_df, material_name::String, property_symbol::Symbol)
	row = findfirst(materials_df.name .== material_name)

	if isnothing(row)
		return NaN  # Return NaN if material_name is not found
	else
		return materials_df[row, property_symbol]
	end
end

"""
	create_sc_cable_model()

Creates a `DataFrame` template that represents the components of a single-core cable, 
which includes the core, sheath, and armor. Each row corresponds to one of these
parts, and the columns represent the physical properties of each part such as radii,
materials, and insulation thicknesses.

# Returns
- `DataFrame`: A template DataFrame for the cable components with columns: `part`, `radius_in`, 
  `radius_ext`, `material_cond`, `width_ins`, `material_ins`, `width_semicon_in`, `width_semicon_ext`.

# Category: Cable and materials database

"""
function create_sc_cable_model()
	# Inner DataFrame template for layers (core, sheath, armor)
	parts_df_template = DataFrame(
		part = ["core", "sheath", "armor"],
		radius_in = Any[NaN, NaN, NaN],     # Internal radius of conductor
		radius_ext = Any[NaN, NaN, NaN],    # External radius of conductor
		material_cond = String["none", "none", "none"],  # Conductor material (reference to materials_df)
		width_ins = Any[NaN, NaN, NaN],     # Thickness of insulation, leave NaN for bare
		material_ins = String["none", "none", "none"],  # Insulation material, leave NaN for bare
		width_semicon_in = Any[NaN, NaN, NaN],  # Inner semiconductor thickness
		width_semicon_ext = Any[NaN, NaN, NaN],  # Outer semiconductor thickness
	)
	return parts_df_template
end

"""
	add_part_to_cable!(cable_parts, part, radius_in, radius_ext, material_cond, width_ins, material_ins, width_semicon_in, width_semicon_ext)

Updates the cable parts `DataFrame` by adding or modifying one of the parts (core, sheath,
armor). It assigns values to the specified part for parameters like internal and external
radius, conductor and insulation materials, and insulation thickness.

# Arguments
- `cable_parts::DataFrame`: The DataFrame representing the cable's components.
- `part::String`: The part to update (e.g., "core", "sheath", "armor").
- `radius_in`: Internal radius of the part.
- `radius_ext`: External radius of the part.
- `material_cond::String`: Material of the conductor.
- `width_ins`: Thickness of the insulation layer.
- `material_ins::String`: Insulation material for the part.
- `width_semicon_in`: Thickness of the inner semiconductor (optional).
- `width_semicon_ext`: Thickness of the outer semiconductor (optional).

# Returns
- `DataFrame`: The updated DataFrame with the new part values.

# Category: Cable and materials database

"""
function add_part_to_cable!(
	cable_parts::DataFrame, part::String, radius_in, radius_ext, material_cond::String,
	width_ins, material_ins::String, width_semicon_in = NaN, width_semicon_ext = NaN,
)
	# Convert floats to Measurement types with zero uncertainty if needed
	radius_in = radius_in isa Measurement ? radius_in : measurement(radius_in, 0.0)
	radius_ext = radius_ext isa Measurement ? radius_ext : measurement(radius_ext, 0.0)
	width_ins = width_ins isa Measurement ? width_ins : measurement(width_ins, 0.0)
	width_semicon_in = width_semicon_in isa Measurement ? width_semicon_in : measurement(width_semicon_in, 0.0)
	width_semicon_ext = width_semicon_ext isa Measurement ? width_semicon_ext : measurement(width_semicon_ext, 0.0)

	# Find the row to update based on the 'part' column
	row_idx = findfirst(cable_parts.part .== part)

	if isnothing(row_idx)
		println("Part $part not found.")
	else
		# Assign values to specific columns for the matching row
		cable_parts[row_idx, :radius_in] = radius_in
		cable_parts[row_idx, :radius_ext] = radius_ext
		cable_parts[row_idx, :material_cond] = material_cond
		cable_parts[row_idx, :width_ins] = width_ins
		cable_parts[row_idx, :material_ins] = material_ins
		cable_parts[row_idx, :width_semicon_in] = width_semicon_in
		cable_parts[row_idx, :width_semicon_ext] = width_semicon_ext
	end

	return cable_parts
end

"""
	create_cables_database()

Creates a template `DataFrame` for storing information about different cables. This `DataFrame`
holds the cable identifier, cross-section areas, resistance, capacitance, and an associated
`DataFrame` for the cable parts (core, sheath, armor).

# Returns
- `DataFrame`: A DataFrame template for storing cable information, with columns: `cable_id`, 
  `cable_cross_section`, `screen_cross_section`, `RDC`, `C`, and `cable_parts`.

# Category: Cable and materials database

"""
function create_cables_database()
	# Create an empty DataFrame to store the cable information
	cables_df_template = DataFrame(
		cable_id = String[],               # Cable ID
		cable_cross_section = Measurement{Float64}[], # Nominal cross-section area of core, in mm²
		screen_cross_section = Measurement{Float64}[], # Nominal cross-section area of screen conductor, in mm²
		armor_cross_section = Measurement{Float64}[], # Nominal cross-section area of armor conductor, in mm²
		RDC = Measurement{Float64}[],  # Measured/nominal DC resistance of core @ 20ºC
		C = Measurement{Float64}[],    # Measured/nominal capacitance of main insulation
		cable_parts = DataFrame[],         # A DataFrame for each part (core, sheath, armor)
	)

	return cables_df_template

end

"""
	add_cable_to_database!(cables_df, cable_id, cable_cross_section, screen_cross_section, cable_parts_df, cable_rdc, cable_c)

Adds a new cable to the cables database, updating it with the parameters of cables such as
cross-section areas, resistance, and capacitance. It also stores the detailed parts
information in a nested DataFrame.

# Arguments
- `cables_df::DataFrame`: The DataFrame holding all cables' data.
- `cable_id::String`: Unique identifier for the cable.
- `cable_cross_section`: Cross-section area of the core conductor.
- `screen_cross_section`: Cross-section area of the screen conductor.
- `cable_parts_df::DataFrame`: DataFrame for the cable parts (core, sheath, armor).
- `cable_rdc::Any`: DC resistance of the cable (optional).
- `cable_c::Any`: Capacitance of the cable (optional).

# Returns
- `DataFrame`: The updated DataFrame containing the new cable information.

# Category: Cable and materials database

"""
function add_cable_to_database!(
	cables_df::DataFrame, cable_id::String, cable_cross_section, screen_cross_section, armor_cross_section, cable_parts_df::DataFrame, cable_rdc = NaN,
	cable_c = NaN,
)

	# Convert floats to Measurement types with zero uncertainty if needed
	cable_cross_section = cable_cross_section isa Measurement ? cable_cross_section : measurement(cable_cross_section, 0.0)
	screen_cross_section = screen_cross_section isa Measurement ? screen_cross_section : measurement(screen_cross_section, 0.0)
	armor_cross_section = armor_cross_section isa Measurement ? armor_cross_section : measurement(armor_cross_section, 0.0)
	cable_rdc = cable_rdc isa Measurement ? cable_rdc : measurement(cable_rdc, 0.0)
	cable_c = cable_c isa Measurement ? cable_c : measurement(cable_c, 0.0)

	push!(cables_df, (cable_id, cable_cross_section, screen_cross_section, armor_cross_section, cable_rdc, cable_c, cable_parts_df))

	return cables_df

end

"""
	create_system_cross_section()

Creates a `DataFrame` template for representing the geometric layout and grounding information
of a multi-conductor cable system. Each row represents a different cable, and the columns
store its position and grounding details.

# Returns
- `DataFrame`: A DataFrame template with columns: `cable_id`, `horz`, `vert`, `ground_sheath`, 
  `ground_armor`.

# Category: Cable and materials database

"""
function create_system_cross_section()

	# Cross-section DataFrame template for each cable geometric and grounding info
	cross_section_df_template = DataFrame(
		cable_id = String[],         # Identifier string of the cable
		horz = Measurement{Float64}[],      # Horizontal coordinate of the cable with uncertainty
		vert = Measurement{Float64}[],      # Vertical coordinate of the cable with uncertainty
		ground_sheath = Bool[],     # Boolean to indicate whether the sheath is grounded
		ground_armor = Bool[],       # Boolean to indicate whether the armor is grounded
	)

	return cross_section_df_template

end

"""
	add_cable_to_cross_section!(cross_section_df, cable_id, horz, vert, ground_sheath, ground_armor)

Adds a cable to the cross-section `DataFrame`, including its position (horizontal and vertical
coordinates) and the grounding statuses for the sheath and armor.

# Arguments
- `cross_section_df::DataFrame`: DataFrame representing the cable system's cross-section.
- `cable_id::String`: The identifier for the cable.
- `horz`: Horizontal position of the cable.
- `vert`: Vertical position of the cable.
- `ground_sheath::Bool`: Whether the sheath is grounded.
- `ground_armor::Bool`: Whether the armor is grounded.

# Returns
- `DataFrame`: The updated cross-section DataFrame with the new cable.

# Category: Cable and materials database

"""
function add_cable_to_cross_section!(
	cross_section_df::DataFrame, cable_id::String, horz, vert, ground_sheath::Bool = true, ground_armor::Bool = true,
)

	# Convert floats to Measurement types with zero uncertainty if needed
	horz = horz isa Measurement ? horz : measurement(horz, 0.0)
	vert = vert isa Measurement ? vert : measurement(vert, 0.0)


	push!(cross_section_df, (cable_id, horz, vert, ground_sheath, ground_armor))

	return cross_section_df

end

"""
	create_cable_system()

Creates a `DataFrame` template to hold information about a cable system, including operation
temperature, soil properties, and cross-sectional layout of the cables.

# Returns
- `DataFrame`: A DataFrame template with columns: `case_id`, `T`, `rho_g`, `epsr_g`, `line_length`, `cross_section`.

# Category: Cable and materials database

"""
function create_cable_system()

	cable_system_df_template = DataFrame(
		case_id = String[],
		T = Measurement{Float64}[],                       # Operation temperature with uncertainty
		rho_g = Measurement{Float64}[],                # Soil resistivity with uncertainty
		epsr_g = Measurement{Float64}[],              # Soil permittivity with uncertainty
		line_length = Measurement{Float64}[], 				# Transmission line length in meters
		cross_section = DataFrame[],   			# The nested cross-section DataFrame
	)

	return cable_system_df_template

end

"""
	define_cable_system!(cable_system_df, case_id, T, rho_g, epsr_g, line_length, cross_section_df)

Defines a cable system by adding a new entry with temperature, soil resistivity, soil
permittivity, and the cross-sectional layout.

# Arguments
- `cable_system_df::DataFrame`: The DataFrame holding all cable system data.
- `case_id::String`: Unique identifier for the cable system.
- `T`: Operating temperature of the cable system.
- `rho_g`: Soil resistivity.
- `epsr_g`: Soil relative permittivity.
- `line_length`: Transmission line length in meters.
- `cross_section_df::DataFrame`: DataFrame representing the cross-sectional layout of the cables.

# Returns
- `DataFrame`: The updated cable system DataFrame.

# Category: Cable and materials database

"""
function define_cable_system!(
	cable_system_df::DataFrame, case_id::String, T, rho_g, epsr_g, line_length, cross_section_df::DataFrame,
)

	# Convert floats to Measurement types with zero uncertainty if needed
	T = T isa Measurement ? T : measurement(T, 0.0)
	rho_g = rho_g isa Measurement ? rho_g : measurement(rho_g, 0.0)
	epsr_g = epsr_g isa Measurement ? epsr_g : measurement(epsr_g, 0.0)
	line_length = line_length isa Measurement ? line_length : measurement(line_length, 0.0)


	push!(cable_system_df, (case_id, T, rho_g, epsr_g, line_length, cross_section_df))

	return cable_system_df

end

"""
	convert_resistance_to_rho(R_DC, r_ext; r_in=0, L=1)

Calculates the resistivity `rho` from the measured DC resistance `R_DC`, conductor
dimensions, and length. This function supports both uncertainty-aware values (`Measurement`)
and regular floats.

# Arguments
- `R_DC::T`: The measured DC resistance (supports both `Measurement` and regular `Float64`).
- `r_ext::T`: The outer radius of the conductor (in meters, supports both `Measurement` and regular `Float64`).
- `r_in::T`: (Optional) The inner radius of the conductor, defaults to `0` for solid conductors.
- `L::T`: (Optional) The length of the conductor, defaults to `1` meter if not provided.

# Returns
- `rho::Measurement`: The calculated resistivity (will always return uncertainty-aware `Measurement`).

# Category: Cable and materials database

"""
function convert_resistance_to_rho(R_DC::Union{Measurement{T}, T}, r_ext::Union{Measurement{T}, T}; r_in::Union{Measurement{T}, T} = 0.0, L::Union{Measurement{T}, T} = 1.0) where {T <: Real}

	# Convert regular floats to Measurement type with zero uncertainty
	R_DC = R_DC isa Measurement ? R_DC : measurement(R_DC, 0.0)
	r_ext = r_ext isa Measurement ? r_ext : measurement(r_ext, 0.0)
	r_in = r_in isa Measurement ? r_in : measurement(r_in, 0.0)
	L = L isa Measurement ? L : measurement(L, 0.0)

	# Calculate the resistivity
	rho = R_DC * pi * (r_ext^2 - r_in^2) / L

	return rho
end

"""
	convert_stranded_to_tubular(rho_c, r_ext, A_c; r_in=0)

Converts the resistivity of a stranded conductor to an equivalent tubular resistivity. This
function supports both uncertainty-aware values (`Measurement`) and regular floats.

# Arguments
- `rho_c::T`: The resistivity of the conductor material (supports both `Measurement` and regular `Float64`).
- `r_ext::T`: The outer radius of the conductor (in meters, supports both `Measurement` and regular `Float64`).
- `A_c::T`: The cross-sectional area of the stranded conductor (in square meters).
- `r_in::T`: (Optional) The inner radius of the conductor, defaults to `0` for solid conductors.

# Returns
- `rho_tubular::Measurement`: The calculated equivalent tubular resistivity (will always return a `Measurement`).

# Category: Cable and materials database

"""
function convert_stranded_to_tubular(rho_c::Union{Measurement{T}, T}, r_ext::Union{Measurement{T}, T}, A_c::Union{Measurement{T}, T}; r_in::Union{Measurement{T}, T} = 0.0) where {T <: Real}

	# Convert regular floats to Measurement type with zero uncertainty
	rho_c = rho_c isa Measurement ? rho_c : measurement(rho_c, 0.0)
	r_ext = r_ext isa Measurement ? r_ext : measurement(r_ext, 0.0)
	A_c = A_c isa Measurement ? A_c : measurement(A_c, 0.0)
	r_in = r_in isa Measurement ? r_in : measurement(r_in, 0.0)

	# Calculate the equivalent tubular resistivity
	rho_tubular = rho_c * pi * (r_ext^2 - r_in^2) / A_c

	return rho_tubular
end

const eps0 = 8.854187817e-12  # Permittivity of free space (F/m)

"""
	convert_capacitance_to_epsilon(C, r_ext; r_in=0, L=1)

Calculates the relative permittivity of insulation based on the capacitance per unit length
of the conductor. This function supports both uncertainty-aware values (`Measurement`) and
regular floats.

# Arguments
- `C::T`: The capacitance per unit length (supports both `Measurement` and regular `Float64`).
- `r_ext::T`: The outer radius of the conductor (in meters, supports both `Measurement` and regular `Float64`).
- `r_in::T`: (Optional) The inner radius of the conductor, defaults to `0` for solid conductors.
- `L::T`: (Optional) The length of the conductor, defaults to `1` meter if not provided.

# Returns
- `eps_ins::Measurement`: The calculated relative permittivity of the insulation (will always return a `Measurement`).

# Category: Cable and materials database

"""
function convert_capacitance_to_epsilon(C::Union{Measurement{T}, T}, r_ext::Union{Measurement{T}, T}; r_in::Union{Measurement{T}, T} = 0.0, L::Union{Measurement{T}, T} = 1.0) where {T <: Real}

	# Convert regular floats to Measurement type with zero uncertainty
	C = C isa Measurement ? C : measurement(C, 0.0)
	r_ext = r_ext isa Measurement ? r_ext : measurement(r_ext, 0.0)
	r_in = r_in isa Measurement ? r_in : measurement(r_in, 0.0)
	L = L isa Measurement ? L : measurement(L, 0.0)


	if r_in == 0.0
		throw(ArgumentError("Inner radius cannot be zero for capacitance-based calculation"))
	end

	# Calculate the relative permittivity of the insulation
	eps_ins = C / (2 * pi * eps0) * log(r_ext / r_in) / L

	return eps_ins
end

"""
	correct_epsilon_for_semicon(eps_ins, radius_inner_conductor, insulation_thickness, radius_over_insulation, semicon_inner, semicon_outer)

Corrects the relative permittivity of insulation by considering the effects of
semiconducting layers in the cable. This function supports both uncertainty-aware inputs
(`Measurement`) and regular floats.

# Arguments
- `eps_ins::T`: The original relative permittivity of the insulation (supports `Measurement` and regular `Float64`).
- `radius_inner_conductor::T`: Outer radius of the inner conductor (supports `Measurement` and regular `Float64`).
- `insulation_thickness::T`: Thickness of the insulation layer (supports `Measurement` and regular `Float64`).
- `radius_over_insulation::T`: Outer radius of the insulation layer (supports `Measurement` and regular `Float64`).
- `semicon_inner::T`: Thickness of the inner semiconducting screen (supports `Measurement` and regular `Float64`).
- `semicon_outer::T`: Thickness of the outer semiconducting screen (supports `Measurement` and regular `Float64`).

# Returns
- `eps_corrected::Measurement`: The corrected relative permittivity of the insulation.

# Category: Cable and materials database

"""
function correct_epsilon_for_semicon(eps_ins::Union{Measurement{T}, T},
	radius_inner_conductor::Union{Measurement{T}, T},
	insulation_thickness::Union{Measurement{T}, T},
	radius_over_insulation::Union{Measurement{T}, T};
	semicon_inner::Union{Measurement{T}, T, Nothing} = nothing,
	semicon_outer::Union{Measurement{T}, T, Nothing} = nothing) where {T <: Real}

	# Ensure all input values are Measurements
	eps_ins = eps_ins isa Measurement ? eps_ins : measurement(eps_ins, 0.0)
	radius_inner_conductor = radius_inner_conductor isa Measurement ? radius_inner_conductor : measurement(radius_inner_conductor, 0.0)
	insulation_thickness = insulation_thickness isa Measurement ? insulation_thickness : measurement(insulation_thickness, 0.0)
	radius_over_insulation = radius_over_insulation isa Measurement ? radius_over_insulation : measurement(radius_over_insulation, 0.0)

	# If semicon thicknesses are not provided, default them to 5% of the insulation thickness
	semicon_inner = isnothing(semicon_inner) ? insulation_thickness * 0.05 : measurement(semicon_inner, 0.0)
	semicon_outer = isnothing(semicon_outer) ? insulation_thickness * 0.05 : measurement(semicon_outer, 0.0)

	# Calculate radii
	a = radius_inner_conductor + semicon_inner  # Inner radius with semicon
	b = a + insulation_thickness  # Outer radius up to insulation

	# Apply the epsilon correction formula
	eps_corrected = eps_ins * (log(radius_over_insulation / radius_inner_conductor) / log(b / a))

	return eps_corrected
end

"""
	calc_equivalent_sheath_radius(A_s, r_in)

Calculates the equivalent sheath radius for a wire screen. The wire screen is replaced by a
fictitious tubular conductor, where the radius is set to match the total wire area `A_s`.

# Arguments
- `A_s::Union{Measurement{T}, T}`: Total wire area of the screen wires, assumed in mm². Can be a `Measurement` 
  (to account for uncertainty) or a regular float.
- `r_in::Union{Measurement{T}, T}`: Inner radius of the sheath (radius of the screen's core). 
  Can be a `Measurement` (to account for uncertainty) or a regular float.

# Returns
- `r_ext::Measurement{T}`: The equivalent sheath radius (with uncertainty), in meters.

# Category: Cable and materials database

"""
function calc_equivalent_sheath_radius(A_s::Union{Measurement{T}, T}, r_in::Union{Measurement{T}, T}) where {T <: Real}

	# Ensure inputs are Measurements
	A_s = A_s isa Measurement ? A_s : measurement(A_s, 0.0)
	r_in = r_in isa Measurement ? r_in : measurement(r_in, 0.0)

	# Calculate equivalent sheath radius
	r_ext = sqrt((A_s / pi / 1e6) + r_in^2)

	return r_ext
end

"""
	calc_stranded_conductor_area(N, D)

Calculates the total cross-sectional area of a stranded conductor, where each wire has a
diameter `D` (with possible uncertainty), and there are `N` wires in total.

# Arguments
- `N::Int`: The number of individual wires in the stranded conductor (certain value).
- `D::Union{Measurement{T}, T}`: The diameter of each individual wire. Can be a `Measurement`
  (with uncertainty) or a regular float.

# Returns
- `A_total::Measurement{T}`: The total cross-sectional area of the stranded conductor.

# Category: Cable and materials database

"""
function calc_stranded_conductor_area(N::Int, D::Union{Measurement{T}, T}) where {T <: Real}

	# Ensure inputs are Measurements
	D = D isa Measurement ? D : measurement(D, 0.0)

	# Calculate total conductor area
	A_total = N * (pi * D^2 / 4)

	return A_total
end

"""
	apply_temperature_correction(rho_base, alpha, T_system, T0)

Applies temperature correction to the base resistivity value `rho_base` using the formula:

	rho = rho_0 * (1 + alpha * (T_system - T0))

# Arguments:
- `rho_base`: The base resistivity of the material at reference temperature `T0`.
- `alpha`: Temperature coefficient of the material.
- `T_system`: Actual operating temperature.
- `T0`: Reference temperature (default to 20°C if not specified in the material properties).

# Returns:
- `rho_corrected`: The corrected resistivity at temperature `T_system`.

# Category: Cable and materials database

"""
function apply_temperature_correction(rho_base, alpha::Float64, T_system, T0::Float64)
	# Ensure inputs are Measurements
	rho_base = rho_base isa Measurement ? rho_base : measurement(rho_base, 0.0)
	T_system = T_system isa Measurement ? T_system : measurement(T_system, 0.0)

	return rho_base * (1 + alpha * (T_system - T0))
end

"""
    trifoil_formation(xc, yc, r_ext)

Calculate the coordinates of three cables laid out in a trifoil pattern, given the center `(xc, yc)` and an external radius `r_ext`.

# Arguments:
- `xc`, `yc`: Coordinates of the trifoil center.
- `r_ext`: External radii of conductors.

# Returns:
    (xa, ya, xb, yb, xc, yc): Coordinates of the three circle centers.

# Category: Cable and materials database

"""
function trifoil_formation(xc, yc, r_ext)
    # Horizontal distance between centers of adjacent circles (equal to twice the radius of each circle)
    d = 2 * r_ext
    # Vertical distance from top circle center to the line between bottom two circles
    h = sqrt(3) * r_ext

    # Calculate the top circle coordinates (centered directly above the midpoint of the bottom two circles)
    xa = xc
    ya = yc + h / 2

    # Calculate the coordinates of the bottom two circles
    xb = xc - d / 2
    yb = yc - h / 2
    xc = xc + d / 2
    yc = yc - h / 2

    return xa, ya, xb, yb, xc, yc
end

"""
    flat_formation(xc, yc, s; vertical=false)

Generates coordinates for three conductors in a flat formation (horizontal or vertical) with a given spacing `s`.

# Arguments:
- `xc`, `yc`: Coordinates of the first conductor.
- `s`: Spacing between conductors.
- `vertical`: Boolean flag indicating whether the layout should be vertical. 
                 If `false`, the layout will be horizontal.

# Returns:
- Tuple of coordinates `(xa, ya, xb, yb, xc, yc)` for the three conductors.

# Category: Cable and materials database

"""
function flat_formation(xc, yc, s; vertical=false)
    if vertical
        # Layout is vertical; adjust only y-coordinates
        xa, ya = xc, yc
        xb, yb = xc, yc - s
        xc, yc = xc, yc - 2s
    else
        # Layout is horizontal; adjust only x-coordinates
        xa, ya = xc, yc
        xb, yb = xc + s, yc
        xc, yc = xc + 2s, yc
    end

    return xa, ya, xb, yb, xc, yc
end

end
