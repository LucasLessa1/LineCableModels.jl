#=
BaraccusMonster - Solução para Análise Parametrizada de Cabos
Baseado no Tutorial 2 do LineCableModels.jl
=#

using DataFrames
using LineCableModels

# Configuração inicial (pode ser ajustada conforme necessário)
setup_logging!(0)

"""
    build_and_analyze_cable(input_set::NamedTuple)

Constrói, configura e analisa um sistema de cabos com base nos parâmetros
fornecidos no `input_set`.

Esta função encapsula a lógica do tutorial, permitindo a fácil variação
de parâmetros geométricos, de materiais e de configuração do sistema.

Retorna os parâmetros base do núcleo do cabo (R, L, C).
"""
function build_and_analyze_cable(input_set::NamedTuple)

    # 1. CONSTRUÇÃO DO NÚCLEO E ISOLAMENTO PRINCIPAL

    # Obtém a biblioteca de materiais customizada do input
    materials = input_set.materials

    #--- Condutor do núcleo (Core Conductor) ---
    material_core = get(materials, "aluminum")
    core = ConductorGroup(WireArray(0.0, Diameter(input_set.d_w), 1, 0.0, material_core))
    n_strands = 6 # Strands per layer
    for i in 1:input_set.n_layers
        add!(core, WireArray, Diameter(d_w), i * n_strands, 11.0, material)
    end

    #--- Isolamento principal e Semicondutores (Main Insulation) ---
    # Fita semicondutora interna
    material_tape = get(materials, "polyacrylate")
    main_insu = InsulatorGroup(Semicon(core, Thickness(input_set.t_sct), material_tape))

    # Semicondutor interno
    material_sc1 = get(materials, "semicon1")
    add!(main_insu, Semicon, Thickness(input_set.t_sc_in), material_sc1)

    # Isolamento principal (XLPE)
    material_ins = get(materials, "pe")
    add!(main_insu, Insulator, Thickness(input_set.t_ins), material_ins)

    # Semicondutor externo
    material_sc2 = get(materials, "semicon2")
    add!(main_insu, Semicon, Thickness(input_set.t_sc_out), material_sc2)

    # Fita semicondutora externa
    add!(main_insu, Semicon, Thickness(input_set.t_sct), material_tape)

    # Agrupa os componentes do núcleo
    core_cc = CableComponent("core", core, main_insu)

    # Instancia o objeto CableDesign
    cable_design = CableDesign(input_set.cable_id, core_cc, nominal_data=input_set.datasheet_info)

    # 2. BLINDAGEM E COBERTURAS

    #--- Blindagem de fios de cobre (Wire Screen) ---
    material_screen = get(materials, "copper")
    screen_con = ConductorGroup(WireArray(main_insu, Diameter(input_set.d_ws), input_set.num_sc_wires, 10.0, material_screen))
    add!(screen_con, Strip, Thickness(input_set.t_cut), input_set.w_cut, 10.0, material_screen)

    # Fita bloqueadora de água sobre a blindagem
    screen_insu = InsulatorGroup(Semicon(screen_con, Thickness(input_set.t_wbt), material_tape))

    # Adiciona o componente da blindagem (sheath)
    add!(cable_design, "sheath", screen_con, screen_insu)


    #--- Cobertura externa (Jacket) ---
    # Fita de alumínio (barreira de umidade)
    jacket_con = ConductorGroup(Tubular(screen_insu, Thickness(input_set.t_alt), material_core)) # Reutiliza alumínio
    jacket_insu = InsulatorGroup(Insulator(jacket_con, Thickness(input_set.t_pet), material_ins)) # Reutiliza PE
    add!(jacket_insu, Insulator, Thickness(input_set.t_jac), material_ins) # Reutiliza PE

    # Adiciona o componente da cobertura (jacket)
    add!(cable_design, "jacket", jacket_con, jacket_insu)

    # 3. ANÁLISE DO SISTEMA TRIFÁSICO

    #--- Modelo de Terra ---
    earth_params = EarthModel(input_set.f_range, input_set.rho_g, input_set.eps_g, input_set.mu_g)

    #--- Sistema Trifásico em Trifoil ---
    xa, ya, xb, yb, xc, yc = trifoil_formation(input_set.x0, input_set.y0, input_set.trifoil_gap)

    # Inicializa o sistema com a fase A
    cablepos = CablePosition(cable_design, xa, ya, Dict("core" => 1, "sheath" => 0, "jacket" => 0))
    cable_system = LineCableSystem("$(input_set.cable_id)_trifoil_system", 1000.0, cablepos)

    # Adiciona as fases B e C
    add!(cable_system, cable_design, xb, yb, Dict("core" => 2, "sheath" => 0, "jacket" => 0))
    add!(cable_system, cable_design, xc, yc, Dict("core" => 3, "sheath" => 0, "jacket" => 0))

    # 4. EXTRAÇÃO DE RESULTADOS

    # Extrai o DataFrame com os parâmetros base (R, L, C)
    core_df = DataFrame(cable_design, :baseparams)

    # Extrai os valores numéricos
    # Acessa a primeira (e única) linha e as colunas por nome
    R = core_df[1,:computed]  # Ω/km
    L = core_df[2,:computed]     # mH/km
    C = core_df[3,:computed]    # μF/km

    # Opcional: Visualizar o cabo finalizado
    # display(preview(cable_design))
    # display(preview(cable_system, zoom_factor=0.15))

    return R, L, C
end


# ÁREA DE CONFIGURAÇÃO E EXECUÇÃO

# 1. Crie uma biblioteca de materiais customizada
custom_materials = MaterialsLibrary(add_defaults=false)
base_aluminum = get(MaterialsLibrary(add_defaults=true), "aluminum")
modified_aluminum = Material(
    base_aluminum.rho * 1.1, # Aumento de 10%
    base_aluminum.eps_r,
    base_aluminum.mu_r,
    base_aluminum.T0,
    base_aluminum.alpha
)

add!(custom_materials, "aluminum",     get(MaterialsLibrary(add_defaults=true), "aluminum"))
add!(custom_materials, "copper",       get(MaterialsLibrary(add_defaults=true), "copper"))
add!(custom_materials, "pe",           get(MaterialsLibrary(add_defaults=true), "pe"))
add!(custom_materials, "semicon1",     get(MaterialsLibrary(add_defaults=true), "semicon1"))
add!(custom_materials, "semicon2",     get(MaterialsLibrary(add_defaults=true), "semicon2"))
add!(custom_materials, "polyacrylate", get(MaterialsLibrary(add_defaults=true), "polyacrylate"))


# 2. Defina o conjunto de parâmetros de entrada (input_set)
input_set = (
    # --- Identificação ---
    cable_id = "18kV_1000mm2_modificado",
    datasheet_info = NominalData( # Dados de referência (opcional, mas útil para comparação)
        designation_code="NA2XS(FL)2Y",
        U0=18.0, U=30.0,
        conductor_cross_section=1000.0, screen_cross_section=35.0,
        resistance=0.0291, capacitance=0.39, inductance=0.3,
    ),

    # --- Geometria ---
    num_co_wires = 61,      # Número de fios do condutor
    num_sc_wires = 49,      # Número de fios da blindagem
    d_w = 4.7e-3,           # Diâmetro do fio do condutor [m]
    d_ws = 0.95e-3,         # Diâmetro do fio da blindagem [m]
    t_sc_in = 0.6e-3,       # Espessura do semicondutor interno [m]
    t_ins = 8.0e-3,         # Espessura do isolamento principal [m]
    t_sc_out = 0.3e-3,      # Espessura do semicondutor externo [m]
    t_cut = 0.1e-3,         # Espessura da fita de cobre [m]
    w_cut = 10e-3,          # Largura da fita de cobre [m]
    t_wbt = 0.3e-3,         # Espessura da fita bloqueadora de água [m]
    t_sct = 0.3e-3,         # Espessura da fita semicondutora [m]
    t_alt = 0.15e-3,        # Espessura da fita de alumínio [m]
    t_pet = 0.05e-3,        # Espessura da face PE na fita de alumínio [m]
    t_jac = 2.4e-3,         # Espessura da cobertura externa de PE [m]
    n_layers = 6,           # Número de camadas de fios do condutor

    # --- Materiais ---
    materials = custom_materials,

    # --- Configuração do Sistema ---
    f_range = 1e-3, # Faixa de frequência [Hz]
    rho_g = 100.0,          # Resistividade do solo [Ω·m]
    eps_g = 10.0,           # Permissividade relativa do solo
    mu_g = 1.0,             # Permeabilidade relativa do solo
    x0 = 0.0,               # Posição central x do sistema [m]
    y0 = -1.0,              # Posição central y (profundidade) do sistema [m]
    trifoil_gap = 0.035,    # Espaçamento entre cabos na formação trifoil [m]
)

# 3. Chame a função e exiba os resultados
R_dc, L, C = build_and_analyze_cable(input_set)

println("="^50)
println("Resultados da Análise para: $(input_set.cable_id)")
println("-"^50)
println("Resistência DC (R_dc): ", round(R_dc, digits=5), " Ω/km")
println("Indutância (L):       ", round(L, digits=4), " mH/km")
println("Capacitância (C):     ", round(C, digits=4), " μF/km")
println("="^50)

# Para um segundo teste, você pode simplesmente modificar o input_set
# Exemplo: Aumentando a espessura do isolamento
modified_input_set = merge(input_set, (t_ins = 9.0e-3, cable_id = "18kV_1000mm2_isol_9mm"))
R_dc_mod, L_mod, C_mod = build_and_analyze_cable(modified_input_set)

println("\nResultados da Análise para: $(modified_input_set.cable_id)")
println("-"^50)
println("Resistência DC (R_dc): ", round(R_dc_mod, digits=5), " Ω/km")
println("Indutância (L):       ", round(L_mod, digits=4), " mH/km")
println("Capacitância (C):     ", round(C_mod, digits=4), " μF/km")
println("="^50)