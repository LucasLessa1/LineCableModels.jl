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

    # @info "Computing symmetrical components"

    # Z012 = CableToolbox.apply_fortescue_transform(Z)
    # Y012 = CableToolbox.apply_fortescue_transform(Y)

    ZY = LineParameters(Z, Y)

    @info "Line parameters computation completed successfully"
    return workspace, ZY
end

# Helper for Kron reduction
function _reorder_by_phase_indices(
    Z::Matrix{Complex{T}},
    ph_order::Vector{Int},
) where {T<:REALSCALAR}
    # Initialize an empty array to store the reordered indices
    reordered_indices = Int[]

    # Get the unique phases from ph_order, ignoring phase 0
    unique_phases = unique(ph_order[ph_order.>0])

    # First, process one row for each unique phase
    for phase in unique_phases
        # Get the first row index corresponding to the current phase
        idx = findfirst(x -> x == phase, ph_order)
        if idx !== nothing
            push!(reordered_indices, idx)
        end
    end

    # Now, append all remaining rows of each phase
    for phase in unique_phases
        # Get all row indices corresponding to the current phase
        idxs = findall(x -> x == phase, ph_order)

        # Remove the first row (it was already added above)
        if length(idxs) > 1
            append!(reordered_indices, idxs[2:end])
        end
    end

    # Finally, append all rows/columns corresponding to phase 0 at the end
    zero_phase_indices = findall(x -> x == 0, ph_order)
    append!(reordered_indices, zero_phase_indices)

    # Reorder both rows and columns of Z based on the reordered indices
    Z_reordered = Z[reordered_indices, reordered_indices]

    # Also reorder the phase order based on the same indices
    ph_reordered = ph_order[reordered_indices]

    return Z_reordered, ph_reordered
end

"""
	do_kron(Z_in, ph_order_in)
    Bundle reduction + Kron elimination
"""
function _do_kron(
    Z_in::Matrix{Complex{T}},
    ph_order_in::Vector{Int},
) where {T<:REALSCALAR}
    # Reorder the matrix Z_in and the phase order based on ph_order_in
    Z, ph_order = _reorder_by_phase_indices(Z_in, ph_order_in)

    num_ph = maximum(ph_order)  # Maximum phase number

    # First Matrix Operation (Z1)
    Z1 = copy(Z)

    for i ∈ 0:num_ph
        ph_pos = findall(x -> x == i, ph_order)
        cond_per_ph = length(ph_pos)
        if !isempty(ph_pos)
            if cond_per_ph > 1
                cond_col = ph_pos[1]
                for j ∈ 2:cond_per_ph
                    subcond_col = ph_pos[j]
                    Z1[:, subcond_col] -= Z[:, cond_col]
                end
            end
        end
    end

    # Second Matrix Operation (Z2)
    Z2 = copy(Z1)

    for i ∈ 0:num_ph
        ph_pos = findall(x -> x == i, ph_order)
        cond_per_ph = length(ph_pos)
        if !isempty(ph_pos)
            if cond_per_ph > 1
                cond_row = ph_pos[1]
                for j ∈ 2:cond_per_ph
                    subcond_row = ph_pos[j]
                    Z2[subcond_row, :] -= Z1[cond_row, :]
                end
            end
        end
    end

    # Apply Kron reduction (ZC is the final matrix before reduction)
    nf = num_ph
    ng = size(Z, 2) - num_ph
    ZC = Z2

    # Kron reduction formula: ZR = ZC11 - ZC12 * inv(ZC22) * ZC21
    ZR =
        ZC[1:nf, 1:nf] -
        ZC[1:nf, (nf+1):(nf+ng)] *
        (ZC[(nf+1):(nf+ng), (nf+1):(nf+ng)] \ ZC[(nf+1):(nf+ng), 1:nf])

    return ZR
end

function kron_reduce(lineparams::LineParameters{T}, ws::CoaxialWorkspace{U}) where {T<:COMPLEXSCALAR,U<:REALSCALAR}
    nfreq = ws.n_frequencies
    phase_map = ws.phase_map
    nph = count(!=(0), phase_map)

    Z_reduced = zeros(T, nph, nph, nfreq)
    Y_reduced = zeros(T, nph, nph, nfreq)

    for f in 1:nfreq
        Zf = lineparams.Z.values[:, :, f]
        Yf = lineparams.Y.values[:, :, f]
        w = 2 * pi * ws.freq[f]
        Zr = _do_kron(Zf, phase_map)
        Pr = _do_kron(inv(Yf / (1im * w)), phase_map)
        Yr = (1im * w) * inv(Pr)
        Z_reduced[:, :, f] = Zr
        Y_reduced[:, :, f] = Yr
    end

    return LineParameters(Z_reduced, Y_reduced)

end