function compute!(problem::LineParametersProblem{T},
    formulation::CoaxialFormulation) where {T}

    # Flatten the system into 1D arrays for computation
    workspace = init_workspace(problem, formulation)

    ComplexT = Complex{T}

    @info "Preallocating arrays"
    n_frequencies = workspace.n_frequencies
    n_phases = workspace.n_phases
    Z = zeros(ComplexT, n_phases, n_phases, n_frequencies)
    Y = zeros(ComplexT, n_phases, n_phases, n_frequencies)

    @info "Starting line parameters computation"

    Z = _compute_impedance_matrix_from_ws(workspace)
    Y = _compute_admittance_matrix_from_ws(workspace)

    @info "Computing symmetrical components"

    Z012 = CableToolbox.apply_fortescue_transform(Z)
    Y012 = CableToolbox.apply_fortescue_transform(Y)

    ZY = LineParameters(Z, Y)
    ZY012 = LineParameters(Z012, Y012)

    @info "Line parameters computation completed successfully"
    return workspace, ZY, ZY012
end

