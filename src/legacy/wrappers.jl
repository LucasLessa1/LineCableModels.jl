
import .CableToolbox: compute_impedance_matrix
using Measurements

# WRAPPERS TO LEGACY CODES
# ---------------------------------------------------------------------------------
"""
Internal: extract fields from a workspace-like object and delegate to the
legacy implementation. Accepts any object with compatible properties.
"""
function _compute_impedance_matrix_from_ws(ws)
    # Strict mapping to CoaxialWorkspace fields
    freq = ws.freq
    horz = ws.horz
    vert = ws.vert
    r_in = ws.r_in
    r_ext = ws.r_ext
    r_ins_ext = ws.r_ins_ext
    rho_cond = ws.rho_cond
    mu_cond = ws.mu_cond
    mu_ins = ws.mu_ins
    eps_ins = ws.eps_ins
    phase_map = ws.phase_map
    cable_map = ws.cable_map

    # Build Geom matrix with columns as specified in LineCableData.jl
    Geom = hcat(
        cable_map,        # column 1: cabID
        phase_map,        # column 2: phID
        horz,             # column 3: horz
        vert,             # column 4: vert
        r_in,             # column 5: r_in
        r_ext,            # column 6: r_ext
        rho_cond,         # column 7: rho_cond
        mu_cond,          # column 8: mu_cond
        r_ins_ext,        # column 9: r_ins (external radius of insulation)
        mu_ins,           # column 10: mu_ins
        eps_ins,          # column 11: eps_ins
    )

    # Ancillary parameters
    Ncables = ws.n_cables
    Nph = count(!=(0), phase_map)
    ph_order = ws.phase_map
    hvec = vert
    dmat = abs.(horz .- horz')

    # ws.earth is a vector (1 per layer); pick last layer
    last_layer = ws.earth[end]

    sigma_g_total = 1.0 ./ last_layer.rho_g    # conductivity vector
    e_g_total = last_layer.eps_g               # absolute permittivity vector
    m_g_total = last_layer.mu_g               # permeability vector (assumed per freq)

    is_meas_input = eltype(ws) <: Measurements.Measurement
    freq_float = Float64.(to_nominal.(freq))
    if !is_meas_input
        # Promote inputs to Measurement with zero uncertainty
        promoteM = x -> Measurements.measurement.(x, 0.0)
        sigma_g_total = promoteM(sigma_g_total)
        e_g_total = promoteM(e_g_total)
        m_g_total = promoteM(m_g_total)
        Geom = promoteM(Geom)
        hvec = promoteM(hvec)
        dmat = promoteM(dmat)
    end
    Zphase = CableToolbox.compute_impedance_matrix(
        freq_float,
        sigma_g_total,
        e_g_total,
        m_g_total,
        Geom,
        Ncables,
        Nph,
        ph_order,
        hvec,
        dmat,
    )
    if !is_meas_input
        # Convert Complex{Measurement} → Complex{Float64} by taking nominal parts
        to_nom = LineCableModels.Utils.to_nominal
        Zphase = map(z -> complex(to_nom(real(z)), to_nom(imag(z))), Zphase)
    end
    return Zphase
end


function _compute_admittance_matrix_from_ws(ws)
    # Strict mapping to CoaxialWorkspace fields
    freq = ws.freq
    horz = ws.horz
    vert = ws.vert
    r_in = ws.r_in
    r_ext = ws.r_ext
    r_ins_ext = ws.r_ins_ext
    rho_cond = ws.rho_cond
    mu_cond = ws.mu_cond
    mu_ins = ws.mu_ins
    eps_ins = ws.eps_ins
    phase_map = ws.phase_map
    cable_map = ws.cable_map

    # Build Geom matrix with columns as specified in LineCableData.jl
    Geom = hcat(
        cable_map,        # column 1: cabID
        phase_map,        # column 2: phID
        horz,             # column 3: horz
        vert,             # column 4: vert
        r_in,             # column 5: r_in
        r_ext,            # column 6: r_ext
        rho_cond,         # column 7: rho_cond
        mu_cond,          # column 8: mu_cond
        r_ins_ext,        # column 9: r_ins (external radius of insulation)
        mu_ins,           # column 10: mu_ins
        eps_ins,          # column 11: eps_ins
    )

    # Ancillary parameters
    Ncables = ws.n_cables
    Nph = count(!=(0), phase_map)
    ph_order = ws.phase_map
    hvec = vert
    dmat = abs.(horz .- horz')

    # ws.earth is a vector (1 per layer); pick last layer
    last_layer = ws.earth[end]

    sigma_g_total = 1.0 ./ last_layer.rho_g    # conductivity vector
    e_g_total = last_layer.eps_g               # absolute permittivity vector
    m_g_total = last_layer.mu_g               # permeability vector (assumed per freq)

    is_meas_input = eltype(ws) <: Measurements.Measurement
    freq_float = Float64.(to_nominal.(freq))
    if !is_meas_input
        # Promote inputs to Measurement with zero uncertainty
        promoteM = x -> Measurements.measurement.(x, 0.0)
        sigma_g_total = promoteM(sigma_g_total)
        e_g_total = promoteM(e_g_total)
        m_g_total = promoteM(m_g_total)
        Geom = promoteM(Geom)
        hvec = promoteM(hvec)
        dmat = promoteM(dmat)
    end
    Yphase = CableToolbox.compute_admittance_matrix(
        freq_float,
        sigma_g_total,
        e_g_total,
        m_g_total,
        Geom,
        Ncables,
        Nph,
        ph_order,
        hvec,
        dmat,
    )
    if !is_meas_input
        # Convert Complex{Measurement} → Complex{Float64} by taking nominal parts
        to_nom = LineCableModels.Utils.to_nominal
        Yphase = map(z -> complex(to_nom(real(z)), to_nom(imag(z))), Yphase)
    end
    return Yphase
end