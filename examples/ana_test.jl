using LineCableModels
using LineCableModels.Engine
using LineCableModels.Engine.FEM
# using Printf
# using Measurements

fullfile(filename) = joinpath(@__DIR__, filename); #hide

materials_db = MaterialsLibrary(add_defaults=false)
load!(materials_db, file_name=fullfile("materials_library.json"))

cables_library = CablesLibrary()
load!(cables_library, file_name=fullfile("cables_library.json"))

# Retrieve the reloaded design
cable_design = get(cables_library, "18kV_1000mm2")
cable_design_s = simplify(cable_design)

# plt1 = preview(cable_design, display_id=true)
# plt2 = preview(cable_design_s, display_id=true)

display(cable_design)
display(cable_design_s)

cable_design_use = cable_design_s

# f = 10.0 .^ range(0, stop=6, length=10)  # Frequency range
f = [50.0]  # Frequency for the analysis
earth_params = EarthModel(f, 100.0, 10.0, 1.0)  # 100 Ω·m resistivity, εr=10, μr=1

# Define system center point (underground at 1 m depth) and the trifoil positions
x0 = 0.0
y0 = -1.0
xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, 0.035)

# Initialize the `LineCableSystem` with the first cable (phase A):
cablepos = CablePosition(cable_design_use, xa, ya, Dict("core" => 1, "sheath" => 0, "jacket" => 0))
cable_system = LineCableSystem("tutorial", 1000.0, cablepos)

# Add remaining cables (phases B and C):
add!(cable_system, cable_design_use, xb, yb,
    Dict("core" => 2, "sheath" => 0, "jacket" => 0),
)
add!(
    cable_system, cable_design_use, xc, yc,
    Dict("core" => 3, "sheath" => 0, "jacket" => 0),
)


# Define a LineParametersProblem with the cable system and earth model
problem = LineParametersProblem(
    cable_system,
    temperature=measurement(20.0),  # Operating temperature
    earth_props=earth_params,
    frequencies=f,  # Frequency for the analysis
);

# Define runtime options 
opts = (
    force_overwrite=true,                    # Overwrite existing files
    save_path=fullfile("lineparams_output"), # Results directory
    verbosity=1,                             # Verbosity
);

# Define the Coaxial model with the specified formulations
formulation = FormulationSet(:Coaxial,
    internal_impedance=InternalImpedance.ScaledBessel(),
    insulation_impedance=InsulationImpedance.Standard(),
    earth_impedance=EarthImpedance.Papadopoulos(),
    insulation_admittance=InsulationAdmittance.Lossless(),
    earth_admittance=EarthAdmittance.Papadopoulos(),
    equivalent_earth=EHEM.EnforceLayer(layer=-1),  # Use the last layer as effective earth
    options=opts
);

# Run solver
@time workspace, ZY, ZY012 = compute!(problem, formulation);

# display(per_km(ZY, 1; mode=:RLCG, freq=f, tol=1e-9))
# display(per_km(ZY012, 1; mode=:RLCG, freq=f, tol=1e-9))
1+1