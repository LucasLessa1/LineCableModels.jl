module CableToolbox

# Import required libraries
using Measurements
using Calculus
using SpecialFunctions  # For Bessel functions
using LinearAlgebra     # For sqrt(-1) as `im`
using Base.MathConstants: pi
using QuadGK
using Printf
using Statistics
using DataFrames
using Plots
using Plots.PlotMeasures

# Export toolbox functions
export compute_impedance_matrix, compute_admittance_matrix, compute_ZY_summary
export calc_outer_skin_effect_impedance,
    calc_inner_skin_effect_impedance,
    calc_mutual_skin_effect_impedance,
    calc_outer_insulation_impedance,
    calc_outer_potential_coefficient
export transform_loop_impedance, perform_bundle_reduction, apply_fortescue_transform
export calc_self_impedance_papadopoulos,
    calc_mutual_impedance_papadopoulos,
    calc_self_potential_coeff_papadopoulos,
    calc_mutual_potential_coeff_papadopoulos,
    compute_frequency_dependent_soil_props
export remove_small_values,
    print_matrices, print_sym_components, uncertain_from_interval, uncertain_from_percent

"""
	@uncertain_bessel(function_name(order, x))

A macro to wrap Bessel functions (or similar) with uncertainty propagation. This computes
the nominal value of the Bessel function and attaches its gradient to account for
uncertainties in the complex input.

# Arguments
- `function_name`: The Bessel function to be evaluated (`besselix`, `besselkx`, etc.).
- `order`: The order of the Bessel function.
- `x`: The input value, which should be a `Complex{Measurements.Measurement{T}}` type to include uncertainties.

# Returns
- `Measurement{T}`: The result of the Bessel function evaluated with uncertainty propagation.

# Category: Misc

"""
macro uncertain_bessel(expr::Expr)
    f = esc(expr.args[1]) # Function name
    order = expr.args[2]  # First argument (order), no need to escape this
    a = esc(expr.args[3]) # Second argument (complex number with uncertainties)

    return :(Measurements.result(
        $f($order, Measurements.value($a)),
        vcat(
            Calculus.gradient(
                x -> real($f($order, complex(x[1], x[2]))),
                [reim(Measurements.value($a))...],
            ),
            Calculus.gradient(
                x -> imag($f($order, complex(x[1], x[2]))),
                [reim(Measurements.value($a))...],
            ),
        ),
        $a,
    ))
end

"""
	calc_outer_skin_effect_impedance(radius_ex, radius_in, sigma_c, mur_c, f; SimplifiedFormula = false)

Calculates the outer surface impedance of a tubular conductor, considering skin effect and
uncertainties.

# Arguments
- `radius_ex::Measurement{T}`: External radius of the conductor.
- `radius_in::Measurement{T}`: Internal radius of the conductor.
- `sigma_c::Measurement{T}`: Electrical conductivity of the conductor.
- `mur_c::Measurement{T}`: Relative permeability of the conductor.
- `f::T`: Frequency of operation (Hz).
- `SimplifiedFormula::Bool`: Use the simplified formula if `true`, full solution otherwise.

# Returns
- `zin::Complex{Measurement{T}}`: The calculated outer surface impedance.

# Category: Cable internal impedance calculations

"""
function calc_outer_skin_effect_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    sigma_c::Measurement{T},
    mur_c::Measurement{T},
    f::T;
    SimplifiedFormula=false,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_c = m0 * mur_c
    TOL = 1e-6
    omega = 2 * pi * f

    # Calculate the reciprocal of the skin depth
    m = sqrt(im * omega * mu_c * sigma_c)

    # Approximated skin effect
    if radius_in == 0
        radius_in = eps()  # Avoid division by zero
    end

    if SimplifiedFormula
        if radius_in < TOL
            cothTerm = coth(m * radius_ex * 0.733)
        else
            cothTerm = coth(m * (radius_ex - radius_in))
        end
        Z1 = (m / sigma_c) / (2 * pi * radius_ex) * cothTerm

        if radius_in < TOL
            Z2 = 0.3179 / (sigma_c * pi * radius_ex^2)
        else
            Z2 = 1 / (sigma_c * 2 * pi * radius_ex * (radius_in + radius_ex))
        end
        zin = Z1 + Z2
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ex
        w_in = m * radius_in

        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besseli() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        N =
            @uncertain_bessel(besselix(0, w_out)) * @uncertain_bessel(besselkx(1, w_in)) +
            sc *
            @uncertain_bessel(besselkx(0, w_out)) *
            @uncertain_bessel(besselix(1, w_in))

        D =
            @uncertain_bessel(besselix(1, w_out)) * @uncertain_bessel(besselkx(1, w_in)) -
            sc *
            @uncertain_bessel(besselkx(1, w_out)) *
            @uncertain_bessel(besselix(1, w_in))

        # Final impedance calculation
        zin = (im * omega * mu_c / (2 * pi)) * (1 / w_out) * (N / D)
    end

    return zin
end

"""
	calc_inner_skin_effect_impedance(radius_ex, radius_in, sigma_c, mur_c, f; SimplifiedFormula = false)

Calculates the inner surface impedance of a tubular conductor, considering skin effect and
uncertainties.

# Arguments
- `radius_ex::Measurement{T}`: External radius of the conductor.
- `radius_in::Measurement{T}`: Internal radius of the conductor.
- `sigma_c::Measurement{T}`: Electrical conductivity of the conductor.
- `mur_c::Measurement{T}`: Relative permeability of the conductor.
- `f::T`: Frequency of operation (Hz).
- `SimplifiedFormula::Bool`: Use the simplified formula if `true`, full solution otherwise.

# Returns
- `zin::Complex{Measurement{T}}`: The calculated inner surface impedance.

# Category: Cable internal impedance calculations

"""
function calc_inner_skin_effect_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    sigma_c::Measurement{T},
    mur_c::Measurement{T},
    f::T;
    SimplifiedFormula=false,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_c = m0 * mur_c
    omega = 2 * pi * f

    # Calculate the reciprocal of the skin depth
    m = sqrt(im * omega * mu_c * sigma_c)

    if SimplifiedFormula

        if isaprox(radius_in, 0)
            zin = (m / (2 * pi * radius_ex * sigma_c)) * coth(0.733 * m * radius_ex) + 0.3179 / (pi * radius_ex^2 * sigma_c)
        else
            cothTerm = coth(m * (radius_ex - radius_in))
            Z1 = (m / sigma_c) / (2 * pi * radius_in) * cothTerm
            Z2 = 1 / (2 * pi * radius_in * (radius_in + radius_ex) * sigma_c)
            zin = Z1 + Z2
        end

    else
        w2 = radius_ex * m

        if radius_in > 0
            # Tubular conductor
            w1 = radius_in * m
            sc1 = exp(abs(real(w1)) - abs(real(w2)) + w1 - w2)  # aligns I/K scaling

            num = @uncertain_bessel(besselix(0, w2)) * @uncertain_bessel(besselkx(1, w1)) +
                  sc1 *
                  @uncertain_bessel(besselix(1, w1)) * @uncertain_bessel(besselkx(0, w2))
            den = @uncertain_bessel(besselix(1, w2)) * @uncertain_bessel(besselkx(1, w1)) -
                  sc1 *
                  @uncertain_bessel(besselix(1, w1)) * @uncertain_bessel(besselkx(1, w2))

            zin = (im * omega * mu_c) / (2 * pi * w2) * (num / den)
        else
            # Solid conductor
            num = @uncertain_bessel(besselix(0, w2))
            den = @uncertain_bessel(besselix(1, w2))
            zin = (im * omega * mu_c) / (2 * pi * w2) * (num / den)
        end
    end

    return zin
end

"""
	calc_mutual_skin_effect_impedance(radius_ex, radius_in, sigma_c, mur_c, f; SimplifiedFormula = false)

Calculates the mutual impedance between outer and inner surfaces of a tubular conductor,
considering skin effect and uncertainties.

# Arguments
- `radius_ex::Measurement{T}`: External radius of the conductor.
- `radius_in::Measurement{T}`: Internal radius of the conductor.
- `sigma_c::Measurement{T}`: Electrical conductivity of the conductor.
- `mur_c::Measurement{T}`: Relative permeability of the conductor.
- `f::T`: Frequency of operation (Hz).
- `SimplifiedFormula::Bool`: Use the simplified formula if `true`, full solution otherwise.

# Returns
- `zm::Complex{Measurement{T}}`: The calculated mutual impedance.

# Category: Cable internal impedance calculations

"""
function calc_mutual_skin_effect_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    sigma_c::Measurement{T},
    mur_c::Measurement{T},
    f::T;
    SimplifiedFormula=false,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_c = m0 * mur_c
    omega = 2 * pi * f

    # Calculate the reciprocal of the skin depth
    m = sqrt(im * omega * mu_c * sigma_c)

    # Approximated skin effect
    if radius_in == 0
        radius_in = eps()  # Avoid division by zero
    end

    if SimplifiedFormula
        cschTerm = csch(m * (radius_ex - radius_in))
        zm = m / (sigma_c * pi * (radius_in + radius_ex)) * cschTerm
    else
        # More detailed solution with Bessel functions and uncertainty
        w_out = m * radius_ex
        w_in = m * radius_in

        s_in = exp(abs(real(w_in)) - w_out)
        s_out = exp(abs(real(w_out)) - w_in)
        sc = s_in / s_out  # Should be applied to all besselix() involving w_in

        # Bessel function terms with uncertainty handling using the macro
        D =
            (@uncertain_bessel besselix(1, w_out)) * (@uncertain_bessel besselkx(1, w_in)) -
            sc *
            (@uncertain_bessel besselkx(1, w_out)) *
            (@uncertain_bessel besselix(1, w_in))

        # Final mutual impedance calculation
        zm = 1 / (2 * pi * radius_ex * radius_in * sigma_c * D * s_out)
    end

    return zm
end

"""
	calc_outer_insulation_impedance(radius_ex, radius_in, mur_ins, f)

Calculates the impedance of the insulation layer in a tubular conductor.

# Arguments
- `radius_ex::Measurement{T}`: External radius of the insulation.
- `radius_in::Measurement{T}`: Internal radius of the insulation.
- `mur_ins::Measurement{T}`: Relative permeability of the insulation.
- `f::T`: Frequency of operation (Hz).

# Returns
- `zinsu::Complex{Measurement{T}}`: The calculated insulation impedance.

# Category: Cable internal impedance calculations

"""
function calc_outer_insulation_impedance(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    mur_ins::Measurement{T},
    f::T,
) where {T<:Real}

    # Constants
    m0 = 4 * pi * 1e-7
    mu_ins = m0 * mur_ins
    omega = 2 * pi * f

    # Avoid division by zero if radius_in is 0
    if radius_in == 0
        radius_in = eps()
    end

    # Impedance of the insulation layer
    zinsu = im * omega * mu_ins * log(radius_ex / radius_in) / (2 * pi)

    return zinsu
end

"""
	transform_loop_impedance(Z)

Transforms the loop-based impedance matrix into phase-domain impedance quantities (for
core/sheath/armor).

# Arguments
- `Z::Matrix{Complex{Measurements.Measurement{T}}}`: The loop-based impedance matrix.

# Returns
- `Z_prime::Matrix{Complex{Measurements.Measurement{T}}}`: The transformed phase-domain impedance matrix.

# Category: Matrix transformations, bundle and Kron reduction

"""
function transform_loop_impedance(
    Z::Matrix{Complex{Measurements.Measurement{T}}},
) where {T<:Real}
    # Check the size of the Z matrix (assuming Z is NxN)
    N = size(Z, 1)

    # Build the voltage transformation matrix T_V
    T_V = Matrix{T}(I, N, N + 1)  # Start with an identity matrix
    for i ∈ 1:N
        T_V[i, i+1] = -1  # Set the -1 in the next column
    end
    T_V = T_V[:, 1:N]  # Remove the last column

    # Build the current transformation matrix T_I
    T_I = tril(ones(T, N, N))  # Lower triangular matrix of ones

    # Compute the new impedance matrix Z_prime
    Z_prime = T_V \ Z * T_I

    return Z_prime
end

"""
	reorder_by_phase_indices(Z, ph_order)

Reorders a given impedance matrix and the corresponding phase order based on phase indices.

# Arguments
- `Z::Matrix{Complex{Measurements.Measurement{T}}}`: The impedance matrix to reorder.
- `ph_order::Vector{Int}`: Vector indicating the phase order of the matrix.

# Returns
- `Z_reordered::Matrix{Complex{Measurements.Measurement{T}}}`: The reordered impedance matrix.
- `ph_reordered::Vector{Int}`: The reordered phase order.

# Category: Matrix transformations, bundle and Kron reduction

"""
function reorder_by_phase_indices(
    Z::Matrix{Complex{Measurements.Measurement{T}}},
    ph_order::Vector{Int},
) where {T<:Real}
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
	perform_bundle_reduction(Z_in, ph_order_in)

Performs bundle and Kron reduction on an impedance matrix to reduce it to the primary
phases.

# Arguments
- `Z_in::Matrix{Complex{Measurements.Measurement{T}}}`: The full impedance matrix.
- `ph_order_in::Vector{Int}`: Initial phase order of the impedance matrix.

# Returns
- `ZR::Matrix{Complex{Measurements.Measurement{T}}}`: The reduced impedance matrix after bundle and Kron reduction.

# Category: Matrix transformations, bundle and Kron reduction

"""
function perform_bundle_reduction(
    Z_in::Matrix{Complex{Measurements.Measurement{T}}},
    ph_order_in::Vector{Int},
) where {T<:Real}
    # Reorder the matrix Z_in and the phase order based on ph_order_in
    Z, ph_order = reorder_by_phase_indices(Z_in, ph_order_in)

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

"""
	calc_outer_potential_coefficient(radius_ex, radius_in, epsr_ins, f, rho_ins = Inf)

Calculates the potential coefficient of an insulation layer in a tubular conductor.

# Arguments
- `radius_ex::Measurement{T}`: External radius of the insulation layer.
- `radius_in::Measurement{T}`: Internal radius of the insulation layer.
- `epsr_ins::Measurement{T}`: Relative permittivity of the insulation layer.
- `f::T`: Frequency of operation (Hz).
- `rho_ins::Union{Measurement{T},T}`: Optional resistivity of the insulation layer (defaults to infinite).

# Returns
- `pm::Complex{Measurement{T}}`: The calculated potential coefficient.

# Category: Cable internal impedance calculations

"""
function calc_outer_potential_coefficient(
    radius_ex::Measurement{T},
    radius_in::Measurement{T},
    epsr_ins::Measurement{T},
    f::T,
    rho_ins::Union{Measurement{T},T}=Inf,
) where {T<:Real}

    # Constant: permittivity of free space
    e0 = 8.854187817e-12
    omega = 2 * pi * f

    # Calculate the permittivity of the insulation layer
    eps_ins = e0 * epsr_ins + 1 / (im * omega * rho_ins)

    # Calculate the potential coefficient
    pm = log(radius_ex / radius_in) / (2 * pi * eps_ins)

    return pm
end

"""
	calc_self_impedance_papadopoulos(h, r, eps_g, mu_g, sigma_g, f, con; kx=0)

Calculates the self-earth impedance of conductors using the Papadopoulos formula,
considering the properties of the ground and the frequency of the system.

# Arguments
- `h`: Vector of vertical distances to ground for conductors (must be negative).
- `r`: Vector of outermost radii of conductors.
- `eps_g`: Relative permittivity of the earth.
- `mu_g`: Permeability of the earth.
- `sigma_g`: Conductivity of the earth.
- `f`: Frequency of the system (in Hz).
- `con`: Number of conductors in the system.
- `kx`: An optional flag for selecting propagation constant:
  - `0`: Default, no propagation constant.
  - `1`: Propagation constant for air.
  - `2`: Propagation constant for earth.

# Returns
- `Matrix{Complex{Measurement{Float64}}}`: The self-earth impedance matrix for the given conductors.

# Category: Earth-return impedances and admittances

"""
function calc_self_impedance_papadopoulos(
    h::Vector{Measurement{Float64}},
    r::Vector{Measurement{Float64}},
    eps_g::Measurement{Float64},
    mu_g::Measurement{Float64},
    sigma_g::Measurement{Float64},
    f::Float64,
    con::Int,
    kx::Int=0,
)

    # Constants
    sig0 = 0.0
    eps0 = 8.8541878128e-12
    mu0 = 4 * pi * 1e-7
    w = 2 * pi * f  # Angular frequency

    # Define k_x based on input kx type
    # 0 = neglect propagation constant
    # 1 = use value of layer 1 (air)
    # 2 = use value of layer 2 (earth)
    k_x = if kx == 2
        ω -> ω * sqrt(mu_g * (eps_g - im * (sigma_g / ω)))
    elseif kx == 1
        ω -> ω * sqrt(mu0 * eps0)
    else
        ω -> ω * 0.0  # Default to zero
    end

    # Define gamma_0 and gamma_1
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))

    # Define a_0 and a_1
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    # Initialize Zg_self matrix (complex numbers)
    Zg_self = zeros(Complex{Measurement{Float64}}, con, con)

    for k ∈ 1:con
        if h[k] < 0  # Only process if h(k) < 0 (as per the original MATLAB logic)

            # Define the function zz
            zz =
                (a0, a1, hi, hj, λ, mu0, mu_g, ω, y) -> (
                    (mu_g * ω * exp(-a1 * abs(hi - hj + 1e-3)) * cos(λ * y) * 0.5im) /
                    (a1 * pi) -
                    (
                        mu_g *
                        ω *
                        exp(-a1 * (hi - hj)) *
                        cos(λ * y) *
                        (a0 * mu_g + a1 * mu0 * sign(hi)) *
                        0.5im
                    ) / (a1 * pi * (a0 * mu_g + a1 * mu0))
                )

            # Define zfun based on lambda and omega (as in the MATLAB code)
            zfun = λ -> begin
                a0 = a_0(λ, w)
                a1 = a_1(λ, w)
                zz(a0, a1, h[k], h[k], λ, mu0, mu_g, w, r[k])
            end

            # Perform the numerical integration (over complex numbers)
            Js, _ = quadgk(zfun, 0.0, Inf; rtol=1e-6)

            # Store the result (which is complex) in Zg_self
            Zg_self[k, k] = Js
        end
    end

    return Zg_self
end

"""
	calc_mutual_impedance_papadopoulos(h, d, eps_g, mu_g, sigma_g, f, con; kx=0)

Calculates the mutual earth impedance between conductors using the Papadopoulos formula,
considering the properties of the ground and the frequency of the system.

# Arguments
- `h`: Vector of vertical distances to ground for conductors.
- `d`: Matrix of distances between conductors.
- `eps_g`: Relative permittivity of the earth.
- `mu_g`: Permeability of the earth.
- `sigma_g`: Conductivity of the earth.
- `f`: Frequency of the system (in Hz).
- `con`: Number of conductors in the system.
- `kx`: An optional flag for selecting propagation constant:
  - `0`: Default, no propagation constant.
  - `1`: Propagation constant for air.
  - `2`: Propagation constant for earth.

# Returns
- `Matrix{Complex{Measurement{Float64}}}`: The mutual earth impedance matrix for the given conductors.

# Category: Earth-return impedances and admittances

"""
function calc_mutual_impedance_papadopoulos(
    h::Vector{Measurement{Float64}},
    d::Matrix{Measurement{Float64}},
    eps_g::Measurement{Float64},
    mu_g::Measurement{Float64},
    sigma_g::Measurement{Float64},
    f::Float64,
    con::Int,
    kx::Int=0,
)

    # Constants
    sig0 = 0.0
    eps0 = 8.8541878128e-12
    mu0 = 4 * pi * 1e-7
    w = 2 * pi * f  # Angular frequency

    # Define k_x based on input kx type
    # 0 = neglect propagation constant
    # 1 = use value of layer 1 (air)
    # 2 = use value of layer 2 (earth)
    k_x = if kx == 2
        ω -> ω * sqrt(mu_g * (eps_g - im * (sigma_g / ω)))
    elseif kx == 1
        ω -> ω * sqrt(mu0 * eps0)
    else
        ω -> ω * 0.0  # Default to zero
    end

    # Define gamma_0 and gamma_1
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))

    # Define a_0 and a_1
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    # Initialize Zg_mutual matrix (complex numbers)
    Zg_mutual = zeros(Complex{Measurement{Float64}}, con, con)

    # Mutual Impedance
    for x ∈ 1:con
        for y ∈ (x+1):con
            if x != y
                h1 = h[x]
                h2 = h[y]

                if h1 < 0 && h2 < 0
                    # Define the function zz
                    zz =
                        (a0, a1, hi, hj, λ, mu0, mu_g, ω, y) -> (
                            (
                                mu_g *
                                ω *
                                exp(-a1 * abs(hi - hj + 1e-3)) *
                                cos(λ * y) *
                                0.5im
                            ) / (a1 * pi) -
                            (
                                mu_g *
                                ω *
                                exp(a1 * (hi + hj)) *
                                cos(λ * y) *
                                (a0 * mu_g + a1 * mu0 * sign(hi)) *
                                0.5im
                            ) / (a1 * pi * (a0 * mu_g + a1 * mu0))
                        )

                    # Define zfun based on lambda and omega
                    zfun = λ -> begin
                        a0 = a_0(λ, w)
                        a1 = a_1(λ, w)
                        zz(a0, a1, h1, h2, λ, mu0, mu_g, w, d[x, y])
                    end

                    # Perform the numerical integration (over complex numbers)
                    Jm, _ = quadgk(zfun, 0.0, Inf; rtol=1e-6)

                    # Store the result (which is complex) in Zg_mutual
                    Zg_mutual[x, y] = Jm
                    Zg_mutual[y, x] = Zg_mutual[x, y]
                end
            end
        end
    end

    return Zg_mutual
end

"""
	calc_self_potential_coeff_papadopoulos(h, r, eps_g, mu_g, sigma_g, f, con; kx=0)

Calculates the self-earth potential coefficient of conductors using the Papadopoulos
formula, considering the ground properties and system frequency.

# Arguments
- `h`: Vector of vertical distances to ground for conductors (must be negative).
- `r`: Vector of outermost radii of conductors.
- `eps_g`: Relative permittivity of the earth.
- `mu_g`: Permeability of the earth.
- `sigma_g`: Conductivity of the earth.
- `f`: Frequency of the system (in Hz).
- `con`: Number of conductors in the system.
- `kx`: An optional flag for selecting propagation constant:
  - `0`: Default, no propagation constant.
  - `1`: Propagation constant for air.
  - `2`: Propagation constant for earth.

# Returns
- `Matrix{Complex{Measurement{Float64}}}`: The self-earth potential coefficient matrix for the given conductors.

# Category: Earth-return impedances and admittances

"""
function calc_self_potential_coeff_papadopoulos(
    h::Vector{Measurement{Float64}},
    r::Vector{Measurement{Float64}},
    eps_g::Measurement{Float64},
    mu_g::Measurement{Float64},
    sigma_g::Measurement{Float64},
    f::Float64,
    con::Int,
    kx::Int=0,
)

    # Constants
    sig0 = 0.0
    eps0 = 8.8541878128e-12
    mu0 = 4 * pi * 1e-7
    w = 2 * pi * f

    # Define k_x based on input kx type
    # 0 = neglect propagation constant
    # 1 = use value of layer 1 (air)
    # 2 = use value of layer 2 (earth)
    k_x = if kx == 2
        ω -> ω * sqrt(mu_g * (eps_g - im * (sigma_g / ω)))
    elseif kx == 1
        ω -> ω * sqrt(mu0 * eps0)
    else
        ω -> ω * 0.0  # Default to zero
    end

    # Define gamma and a functions
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    Pg_self = zeros(Complex{Measurement{Float64}}, con, con)
    TOL = 1e-3

    for k ∈ 1:con
        if h[k] < 0
            yy =
                (a0, a1, gamma0, gamma1, hi, hj, λ, mu0, mu_g, ω, y) ->
                    (
                        1.0 / gamma1^2 *
                        mu_g *
                        ω *
                        exp(-a1 * abs(hi - hj + TOL)) *
                        cos(λ * y) *
                        0.5im
                    ) / (a1 * pi) -
                    (
                        1.0 / gamma1^2 *
                        mu_g *
                        ω *
                        exp(a1 * (hi + hj)) *
                        cos(λ * y) *
                        (a0 * mu_g + a1 * mu0 * sign(hi)) *
                        0.5im
                    ) / (a1 * pi * (a0 * mu_g + a1 * mu0)) +
                    (
                        a1 * 1.0 / gamma1^2 *
                        mu0 *
                        mu_g^2 *
                        ω *
                        exp(a1 * (hi + hj)) *
                        cos(λ * y) *
                        (sign(hi) - 1.0) *
                        (gamma0^2 - gamma1^2) *
                        0.5im
                    ) / (
                        pi *
                        (a0 * gamma1^2 * mu0 + a1 * gamma0^2 * mu_g) *
                        (a0 * mu_g + a1 * mu0)
                    )

            yfun =
                λ -> yy(
                    a_0(λ, w),
                    a_1(λ, w),
                    gamma_0(w),
                    gamma_1(w),
                    h[k],
                    h[k],
                    λ,
                    mu0,
                    mu_g,
                    w,
                    r[k],
                )
            Qs, _ = quadgk(yfun, 0, Inf, rtol=1e-6)
            Pg_self[k, k] = (im * w * Qs)
        end
    end

    return Pg_self
end

"""
	calc_mutual_potential_coeff_papadopoulos(h, d, eps_g, mu_g, sigma_g, f, con; kx=0)

Calculates the mutual earth potential coefficient between conductors using the Papadopoulos
formula, considering the properties of the ground and the frequency of the system.

# Arguments
- `h`: Vector of vertical distances to ground for conductors.
- `d`: Matrix of distances between conductors.
- `eps_g`: Relative permittivity of the earth.
- `mu_g`: Permeability of the earth.
- `sigma_g`: Conductivity of the earth.
- `f`: Frequency of the system (in Hz).
- `con`: Number of conductors in the system.
- `kx`: An optional flag for selecting propagation constant:
  - `0`: Default, no propagation constant.
  - `1`: Propagation constant for air.
  - `2`: Propagation constant for earth.

# Returns
- `Matrix{Complex{Measurement{Float64}}}`: The mutual earth potential coefficient matrix for the given conductors.

# Category: Earth-return impedances and admittances

"""
function calc_mutual_potential_coeff_papadopoulos(
    h::Vector{Measurement{Float64}},
    d::Matrix{Measurement{Float64}},
    eps_g::Measurement{Float64},
    mu_g::Measurement{Float64},
    sigma_g::Measurement{Float64},
    f::Float64,
    con::Int,
    kx::Int=0,
)

    # Constants
    sig0 = 0.0
    eps0 = 8.8541878128e-12
    mu0 = 4 * pi * 1e-7
    w = 2 * pi * f

    # Define k_x based on input kx type
    # 0 = neglect propagation constant
    # 1 = use value of layer 1 (air)
    # 2 = use value of layer 2 (earth)
    k_x = if kx == 2
        ω -> ω * sqrt(mu_g * (eps_g - im * (sigma_g / ω)))
    elseif kx == 1
        ω -> ω * sqrt(mu0 * eps0)
    else
        ω -> ω * 0.0  # Default to zero
    end

    # Define gamma and a functions
    gamma_0 = ω -> sqrt(im * ω * mu0 * (sig0 + im * ω * eps0))
    gamma_1 = ω -> sqrt(im * ω * mu_g * (sigma_g + im * ω * eps_g))
    a_0 = (λ, ω) -> sqrt(λ^2 + gamma_0(ω)^2 + k_x(ω)^2)
    a_1 = (λ, ω) -> sqrt(λ^2 + gamma_1(ω)^2 + k_x(ω)^2)

    Pg_mutual = zeros(Complex{Measurement{Float64}}, con, con)
    TOL = 1e-3

    # Mutual potential coefficient
    for x ∈ 1:con
        for y ∈ (x+1):con
            if x != y
                h1 = h[x]
                h2 = h[y]
                if abs(h2 - h1) < TOL
                    h2 += TOL
                end

                if h1 < 0 && h2 < 0
                    yy =
                        (a0, a1, gamma0, gamma1, hi, hj, λ, mu0, mu_g, ω, y) ->
                            (
                                1.0 / gamma1^2 *
                                mu_g *
                                ω *
                                exp(-a1 * abs(hi - hj)) *
                                cos(λ * y) *
                                0.5im
                            ) / (a1 * pi) -
                            (
                                1.0 / gamma1^2 *
                                mu_g *
                                ω *
                                exp(a1 * (hi + hj)) *
                                cos(λ * y) *
                                (a0 * mu_g + a1 * mu0 * sign(hi)) *
                                0.5im
                            ) / (a1 * pi * (a0 * mu_g + a1 * mu0)) +
                            (
                                a1 * 1.0 / gamma1^2 *
                                mu0 *
                                mu_g^2 *
                                ω *
                                exp(a1 * (hi + hj)) *
                                cos(λ * y) *
                                (sign(hi) - 1.0) *
                                (gamma0^2 - gamma1^2) *
                                0.5im
                            ) / (
                                pi *
                                (a0 * gamma1^2 * mu0 + a1 * gamma0^2 * mu_g) *
                                (a0 * mu_g + a1 * mu0)
                            )

                    yfun =
                        λ -> yy(
                            a_0(λ, w),
                            a_1(λ, w),
                            gamma_0(w),
                            gamma_1(w),
                            h1,
                            h2,
                            λ,
                            mu0,
                            mu_g,
                            w,
                            d[x, y],
                        )
                    Qm, _ = quadgk(yfun, 0, Inf, rtol=1e-3)
                    Pg_mutual[x, y] = (im * w * Qm)
                    Pg_mutual[y, x] = Pg_mutual[x, y]
                end
            end
        end
    end

    return Pg_mutual
end

"""
	compute_frequency_dependent_soil_props(soil, FD_flag, f_total)

Computes the frequency-dependent soil properties for multiple soil layers. The returned
dictionary contains the relative permittivity, conductivity, and absolute permittivity for
each soil layer across all frequency samples.

# Arguments
- `soil`: A dictionary containing the base properties of the soil (including optional layers).
- `FD_flag`: A flag to select the type of frequency-dependent model:
  - `0`: Constant soil properties (default).
- `f_total`: Vector of frequencies for which the soil properties will be computed.

# Returns
- `Dict{Symbol,Any}`: A dictionary with frequency-dependent soil properties:
  - `:erg_total`: Relative permittivity for each layer and frequency.
  - `:sigma_g_total`: Conductivity for each layer and frequency.
  - `:e_g_total`: Absolute permittivity for each layer and frequency.
  - `:m_g_total`: Permeability (constant for all frequencies).

# Category: Frequency-dependent (FD) soil properties

"""
function compute_frequency_dependent_soil_props(
    soil,
    FD_flag::Int,
    f_total::Vector{Float64},
)

    # Initialize the output soil structure
    soilFD = Dict{Symbol,Any}()

    # Constants
    e0 = 8.854187817e-12  # Farads/meter
    siz = length(f_total)

    # Initialize number of layers
    num_layers = 1
    if haskey(soil, :layer)
        num_layers += length(soil[:layer])
    end

    # Initialize matrices for soil properties across frequencies
    erg_total = fill(measurement(0.0, 0.0), siz, num_layers)  # Relative permittivity of soil layers
    sigma_g_total = fill(measurement(0.0, 0.0), siz, num_layers)  # Conductivity of soil layers
    e_g_total = fill(measurement(0.0, 0.0), siz, num_layers)  # Absolute permittivity of soil layers

    for l ∈ 1:num_layers
        # Extract soil properties for each layer (or final layer if no more layers)
        if l == num_layers
            epsr = soil[:erg]      # Relative permittivity of earth
            m_g = soil[:m_g]       # Permeability of earth
            sigma_g = soil[:sigma_g]  # Conductivity of earth
        else
            epsr = soil[:layer][l][:erg]  # Relative permittivity of earth for this layer
            m_g = soil[:layer][l][:m_g]   # Permeability of earth for this layer
            sigma_g = soil[:layer][l][:sigma_g]  # Conductivity of earth for this layer
        end

        # For FD_flag == 0, the parameters are constant over frequencies
        if FD_flag == 0
            erg_total[:, l] .= epsr .* ones(siz)
            sigma_g_total[:, l] .= sigma_g .* ones(siz)
        else
            Base.error("Currently only constant properties (CP) soil models are available.") #TODO: port CIGRE and Longmire-Smith models
            #Issue URL: https://github.com/Electa-Git/LineCableModels.jl/issues/15
        end

        # Calculate absolute permittivity (ε = ε₀ * εr)
        e_g_total[:, l] .= e0 .* erg_total[:, l]
    end

    # Add frequency-dependent soil properties to the output dictionary
    soilFD[:erg_total] = erg_total
    soilFD[:sigma_g_total] = sigma_g_total
    soilFD[:e_g_total] = e_g_total
    soilFD[:m_g_total] = fill(soil[:m_g], siz, num_layers)  # Permeability is constant

    return soilFD
end


"""
	remove_small_values(data)

Cleans up the impedance/admittance values by replacing values smaller than machine epsilon
with zero.

# Arguments
- `data`: Input scalar, 2D matrix, or nD array to clean.

# Returns
- `Complex{Measurements.Measurement{Float64}}`: Cleaned data with small values replaced by zero.

# Category: Misc

"""
function remove_small_values(data)
    eps_value = eps(Float64)

    if isa(data, Complex{Measurements.Measurement{Float64}})
        # Handle scalar case
        return Complex(
            abs(Measurements.value(real(data))) < eps_value ?
            Measurements.measurement(0.0, 0.0) : real(data),
            abs(Measurements.value(imag(data))) < eps_value ?
            Measurements.measurement(0.0, 0.0) : imag(data),
        )
    elseif isa(data, AbstractArray)
        # Handle any-dimensional array case
        return Complex{Measurements.Measurement{Float64}}[
            Complex(
                abs(Measurements.value(real(x))) < eps_value ?
                Measurements.measurement(0.0, 0.0) : real(x),
                abs(Measurements.value(imag(x))) < eps_value ?
                Measurements.measurement(0.0, 0.0) : imag(x),
            ) for x in data
        ]
    else
        Base.error("Unsupported data type.")
    end
end


"""
	compute_impedance_matrix(freq, soilFD, Geom, Ncables, Nph, ph_order, h, d)

Computes the full impedance matrix for a multi-conductor system over a range of frequencies.

# Arguments
- `freq::Vector{Float64}`: Vector of frequencies to evaluate.
- `soilFD::Dict{Symbol,Any}`: Soil parameters including frequency-dependent properties.
- `Geom::Matrix{Float64}`: Geometry of the cable system.
- `Ncables::Int`: Number of cables in the system.
- `Nph::Int`: Number of phases in the system.
- `ph_order::Vector{Int}`: Phase order of the system.
- `h::Vector{Measurement{Float64}}`: Heights of the conductors above the ground.
- `d::Matrix{Measurement{Float64}}`: Distance matrix between conductors.

# Returns
- `Zphase::Array{Complex{Measurements.Measurement{Float64}}, 3}`: The computed impedance matrix for each frequency.

# Category: Cable constants matrices
"""
function compute_impedance_matrix(freq, sigma_g_total, e_g_total, m_g_total, Geom, Ncables, Nph, ph_order, h, d)

    freq_siz = length(freq)
    # Preallocate the 3D array for impedance with size Nph x Nph x freq_siz
    Zphase = Array{Complex{Measurements.Measurement{Float64}},3}(undef, Nph, Nph, freq_siz)

    for k ∈ 1:freq_siz

        # Initialize Z as a matrix of Any type (to hold matrices of varying sizes)
        Z = Matrix{Any}(undef, Ncables, Ncables)

        f = freq[k]

        # Get soil parameters for this frequency
        sigma_g = sigma_g_total[k]
        m_g = m_g_total[k]
        e_g = e_g_total[k]

        # Calculate outermost radii for each cable
        outermost_radii = map(
            cable_num -> maximum(
                [Geom[Geom[:, 1].==cable_num, 6]; Geom[Geom[:, 1].==cable_num, 9]],
            ),
            1:Ncables,
        )

        # Calculate external (earth return) impedance matrices
        Zg_self = calc_self_impedance_papadopoulos(
            h[ph_order.!=0],
            outermost_radii,
            0 * e_g,
            m_g,
            sigma_g,
            f,
            Ncables,
            0,
        )
        Zg_mutual = calc_mutual_impedance_papadopoulos(
            h[ph_order.!=0],
            d[ph_order.!=0, ph_order.!=0],
            0 * e_g,
            m_g,
            sigma_g,
            f,
            Ncables,
            0,
        )
        Zext = Zg_self + Zg_mutual

        # Compute internal impedance matrices
        for i ∈ 1:Ncables
            cabledata = Geom[Geom[:, 1].==i, 2:end]
            ncond_cable = size(cabledata, 1)  # Number of conductor layers in cable i

            for j ∈ 1:Ncables
                # Initialize loop-based impedance matrix for cable i
                Zloop = Matrix{Complex{Measurements.Measurement{Float64}}}(
                    undef,
                    ncond_cable,
                    ncond_cable,
                )
                Zloop .= 0.0 + 0.0im  # Fill with zeros

                if i == j
                    # Self-impedance
                    for k ∈ 1:ncond_cable
                        # Outer surface impedance of this layer
                        rin_thislayer = cabledata[k, 4]
                        rext_thislayer = cabledata[k, 5]
                        sig_c_thislayer = 1 / cabledata[k, 6]
                        mu_c_thislayer = cabledata[k, 7]
                        Zouter_thislayer = calc_outer_skin_effect_impedance(
                            rext_thislayer,
                            rin_thislayer,
                            sig_c_thislayer,
                            mu_c_thislayer,
                            f,
                            SimplifiedFormula=false,
                        )

                        # Insulation impedance (if present)
                        if !any(x -> isnan(Measurements.value(x)), cabledata[k, 8:end])
                            rext_thisinsu = cabledata[k, 8]
                            mu_thisinsu = cabledata[k, 9]
                            Zinsu_thislayer = calc_outer_insulation_impedance(
                                rext_thisinsu,
                                rext_thislayer,
                                mu_thisinsu,
                                f,
                            )
                        else
                            Zinsu_thislayer = 0.0 + 0.0im
                        end

                        # Inner surface impedance of the next layer
                        if k < ncond_cable
                            rin_nextlayer = cabledata[k+1, 4]
                            rext_nextlayer = cabledata[k+1, 5]
                            sig_c_nextlayer = 1 / cabledata[k+1, 6]
                            mu_c_nextlayer = cabledata[k+1, 7]
                            Zinner_nextlayer = calc_inner_skin_effect_impedance(
                                rext_nextlayer,
                                rin_nextlayer,
                                sig_c_nextlayer,
                                mu_c_nextlayer,
                                f,
                                SimplifiedFormula=false,
                            )
                            Zmutual = calc_mutual_skin_effect_impedance(
                                rext_nextlayer,
                                rin_nextlayer,
                                sig_c_nextlayer,
                                mu_c_nextlayer,
                                f,
                                SimplifiedFormula=false,
                            )
                        else
                            Zinner_nextlayer = Zext[i, j]  # Self-earth return impedance (underground conductor)
                            Zmutual = 0.0 + 0.0im
                        end

                        # Assign to Zloop
                        Zloop[k, k] = Zouter_thislayer + Zinsu_thislayer + Zinner_nextlayer

                        if k < ncond_cable
                            Zloop[k, k+1] = -Zmutual
                            Zloop[k+1, k] = -Zmutual
                        end
                    end
                else
                    # Mutual impedance
                    Zloop[end, end] = Zext[i, j]
                end

                # Transform loop impedance to phase-domain impedance
                Z[i, j] = transform_loop_impedance(Zloop)
            end
        end

        # Combine all impedance matrices into a full system impedance matrix
        Zfull = reduce(vcat, [reduce(hcat, Z[i, :]) for i ∈ 1:Ncables])

        # Perform bundle and Kron reduction
        Zph = perform_bundle_reduction(Zfull, ph_order)
        Zphase[:, :, k] = (Zph + transpose(Zph)) / 2.0 # enforce symmetry

    end

    return Zphase
end

"""
	compute_admittance_matrix(freq, soilFD, Geom, Ncables, Nph, ph_order, h, d)

Computes the full admittance matrix for a multi-conductor system over a range of
frequencies.

# Arguments
- `freq::Vector{Float64}`: Vector of frequencies to evaluate.
- `soilFD::Dict{Symbol,Any}`: Soil parameters including frequency-dependent properties.
- `Geom::Matrix{Float64}`: Geometry of the cable system.
- `Ncables::Int`: Number of cables in the system.
- `Nph::Int`: Number of phases in the system.
- `ph_order::Vector{Int}`: Phase order of the system.
- `h::Vector{Measurement{Float64}}`: Heights of the conductors above the ground.
- `d::Matrix{Measurement{Float64}}`: Distance matrix between conductors.

# Returns
- `Yphase::Array{Complex{Measurements.Measurement{Float64}}, 3}`: The computed admittance matrix for each frequency.

# Category: Cable constants matrices

"""
function compute_admittance_matrix(freq, sigma_g_total, e_g_total, m_g_total, Geom, Ncables, Nph, ph_order, h, d)

    freq_siz = length(freq)

    # Preallocate the 3D array for admittance with size Nph x Nph x freq_siz
    Yphase = Array{Complex{Measurements.Measurement{Float64}},3}(undef, Nph, Nph, freq_siz)

    for k ∈ 1:freq_siz

        # Initialize P as a matrix of Any type (to hold matrices of varying sizes)
        P = Matrix{Any}(undef, Ncables, Ncables)
        f = freq[k]

        # Get soil parameters for this frequency
        sigma_g = sigma_g_total[k]
        m_g = m_g_total[k]
        e_g = e_g_total[k]

        # Calculate outermost radii for each cable
        outermost_radii = map(
            cable_num -> maximum(
                [Geom[Geom[:, 1].==cable_num, 6]; Geom[Geom[:, 1].==cable_num, 9]],
            ),
            1:Ncables,
        )

        # Calculate external (earth return) admittance matrices
        Pg_self = calc_self_potential_coeff_papadopoulos(
            h[ph_order.!=0],
            outermost_radii,
            0 * e_g,
            m_g,
            sigma_g,
            f,
            Ncables,
            0,
        )
        Pg_mutual = calc_mutual_potential_coeff_papadopoulos(
            h[ph_order.!=0],
            d[ph_order.!=0, ph_order.!=0],
            0 * e_g,
            m_g,
            sigma_g,
            f,
            Ncables,
            0,
        )
        Pext = Pg_self + Pg_mutual

        # Compute internal admittance matrices
        for i ∈ 1:Ncables
            cabledata = Geom[Geom[:, 1].==i, 2:end]
            ncond_cable = size(cabledata, 1)  # Number of conductor layers in cable i

            for j ∈ 1:Ncables
                # Initialize loop-based admittance matrix for cable i
                Pcable = Matrix{Complex{Measurements.Measurement{Float64}}}(
                    undef,
                    ncond_cable,
                    ncond_cable,
                )
                Pcable .= 0.0 + 0.0im  # Fill with zeros

                if i == j
                    # Self-admittance
                    for k ∈ 1:ncond_cable
                        for l ∈ k:ncond_cable
                            Pkl = 0.0 + 0.0im  # Initialize as complex value
                            for m ∈ l:ncond_cable
                                if !any(
                                    x -> isnan(Measurements.value(x)),
                                    cabledata[k, 8:end],
                                )
                                    radius_in = cabledata[m, 5]    # Inner radius of insulation layer
                                    radius_ex = cabledata[m, 8]    # Outer radius of insulation layer
                                    eps_r = cabledata[m, 10]       # Relative permittivity of insulation layer

                                    # Call the calc_outer_potential_coefficient function and accumulate along rows
                                    Pkl += calc_outer_potential_coefficient(
                                        radius_ex,
                                        radius_in,
                                        eps_r,
                                        f,
                                    )
                                end
                            end
                            Pcable[k, l] = Pkl
                        end
                    end
                    # Mirror the upper triangle onto the lower triangle
                    Pcable = Pcable + transpose(triu(Pcable, 1))
                end
                P[i, j] = Pcable .+ Pext[i, j]
            end
        end

        # Combine all admittance matrices into a full system admittance matrix
        Pfull = reduce(vcat, [reduce(hcat, P[i, :]) for i ∈ 1:Ncables])

        # Perform bundle and Kron reduction
        Pphase = perform_bundle_reduction(Pfull, ph_order)

        # Compute the reduced admittance matrix
        w = 2 * pi * f
        Yph = 1im * w * inv(Pphase)
        Yphase[:, :, k] = (Yph + transpose(Yph)) / 2.0 # enforce symmetry

    end

    return Yphase
end

"""
	apply_fortescue_transform2(Z)

Applies the Fortescue transformation to the impedance matrix to convert it into modal
components (zero, positive, and negative sequence).

# Arguments
- `Z::Array{Complex{Measurement{Float64}}}`: The 3D impedance matrix with conductors and frequencies.

# Returns
- `Z012::Array{Complex{Measurement{Float64}}}`: The transformed modal impedance matrix for each frequency.

# Category: Matrix transformations, bundle and Kron reduction

"""
function apply_fortescue_transform(
    Z::Union{Array{Complex{Measurement{Float64}}},Array{Complex{Float64}}};
    transpose_line::Bool=true,
)
    num_phases = size(Z, 1)  # Fetch the number of conductors from Z dimensions
    num_freq_samples = size(Z, 3)  # Number of frequency samples (3rd dimension of Z)

    Ti = nothing

    if num_phases == 3
        # Fortescue transformation matrix for 3-phase systems
        A = exp(1im * 2 * pi / 3)  # Rotational constant
        Ti = [1 1 1; 1 A^2 A; 1 A A^2]

    else
        throw(ArgumentError("Unsupported number of phases: $num_phases"))
    end

    # Preallocate Z012 as a 3D array with full transformed 3x3 matrices for each frequency sample
    Z012 = Array{Complex{Measurement{Float64}},3}(
        undef,
        num_phases,
        num_phases,
        num_freq_samples,
    )


    # Loop over each frequency sample and apply the modal transformation
    for k ∈ 1:num_freq_samples
        # Extract the 2D matrix for the current frequency sample
        Z_slice = Z[:, :, k]
        if transpose_line
            # Enforce transposition (snaking)
            diag_mean = mean(diag(Z_slice))
            off_diag_values = Z_slice[triu(trues(size(Z_slice)), 1)]  # Upper triangular off-diagonal
            append!(off_diag_values, Z_slice[tril(trues(size(Z_slice)), -1)])  # Lower triangular off-diagonal
            off_diag_mean = mean(off_diag_values)  # Compute the mean
            Z_slice =
                diagm(fill(diag_mean, 3)) .+
                (ones(3, 3) .- Matrix(I, 3, 3)) .* off_diag_mean
        end
        # Apply transformation to the 2D matrix
        Z_transformed = inv(Ti) * Z_slice * Ti
        Z012[:, :, k] = Z_transformed

    end

    return Z012
end

"""
	print_matrices(Z, freq_range)

Formats and prints a 3D matrix `Z`, which represents frequency-dependent impedance/admittance 
values, in a clean and human-readable way. Each frequency's matrix is printed with rows and columns 
formatted neatly for easy reading, and the real and imaginary parts are shown with 6 decimal places.

# Arguments
- `Z::Array{Complex{Measurement{T}}, 3}`: A 3D array where the first two dimensions represent the rows and columns 
  of the matrix, and the third dimension corresponds to different frequency samples.
- `frequencies::Vector{T}`: A vector of frequencies corresponding to the third dimension of `Z`.

# Category: Misc

"""
function print_matrices(Z, freq_range)
    num_freqs = size(Z, 3)  # 3rd dimension is the number of frequencies

    for k in 1:num_freqs
        Zprint = remove_small_values(Z[:, :, k])
        println("  - Phase-based quantities for frequency f = $(freq_range[k]) Hz:")
        println("[")  # Open matrix bracket

        for i in 1:size(Z, 1)

            row_values = [
                imag(Zprint[i, j]) >= 0 ?
                @sprintf(
                    "%.6E ± %.6E + %.6E ± %.6Eim",
                    real(Measurements.value(Zprint[i, j])),
                    Measurements.uncertainty(real(Zprint[i, j])),
                    imag(Measurements.value(Zprint[i, j])),
                    Measurements.uncertainty(imag(Zprint[i, j]))
                ) :
                @sprintf(
                    "%.6E ± %.6E - %.6E ± %.6Eim",
                    real(Measurements.value(Zprint[i, j])),
                    Measurements.uncertainty(real(Zprint[i, j])),
                    abs(imag(Measurements.value(Zprint[i, j]))),
                    Measurements.uncertainty(imag(Zprint[i, j]))
                )
                for j in 1:size(Z, 2)
            ]

            if i < size(Z, 1)
                println("   ", join(row_values, ", "), ";")  # Add semicolon at the end of each row except the last
            else
                println("   ", join(row_values, ", "))  # No semicolon for the last row
            end
        end

        println("]")  # Close matrix bracket
        println("\n")  # Linebreak between frequency matrices
    end
end

"""
	print_sym_components(Z, freq_range)

Formats and prints the diagonal elements of a 3D matrix `Z`, which represents the sequence impedance/admittance 
values, in a clean and human-readable way.

# Arguments
- `Z::Array{Complex{Measurement{T}}, 3}`: A 3D array where the first two dimensions represent the rows and columns 
  of the matrix, and the third dimension corresponds to different frequency samples.
- `frequencies::Vector{T}`: A vector of frequencies corresponding to the third dimension of `Z`.

# Category: Misc

"""
function print_sym_components(Z, freq_range)
    num_freqs = size(Z, 3)  # 3rd dimension is the number of frequencies

    for k in 1:num_freqs
        Zprint = diag(remove_small_values(Z[:, :, k]))
        println(
            "  - Symmetrical components ordered as 012 for frequency f = $(freq_range[k]) Hz:",
        )
        println("[")  # Open matrix bracket

        for i in axes(Zprint, 1)

            row_values = [
                imag(Zprint[i, j]) >= 0 ?
                @sprintf(
                    "%.6E ± %.6E + %.6E ± %.6Eim",
                    real(Measurements.value(Zprint[i, j])),
                    Measurements.uncertainty(real(Zprint[i, j])),
                    imag(Measurements.value(Zprint[i, j])),
                    Measurements.uncertainty(imag(Zprint[i, j]))
                ) :
                @sprintf(
                    "%.6E ± %.6E - %.6E ± %.6Eim",
                    real(Measurements.value(Zprint[i, j])),
                    Measurements.uncertainty(real(Zprint[i, j])),
                    abs(imag(Measurements.value(Zprint[i, j]))),
                    Measurements.uncertainty(imag(Zprint[i, j]))
                )
                for j in axes(Zprint, 2)
            ]

            if i < size(Z, 1)
                println("   ", join(row_values, ", "), ";")  # Add semicolon at the end of each row except the last
            else
                println("   ", join(row_values, ", "))  # No semicolon for the last row
            end
        end

        println("]")  # Close matrix bracket
        println("\n")  # Linebreak between frequency matrices
    end
end

"""
	uncertain_from_interval(max, min)

Creates a Measurement object from max/min interval. 

# Arguments
- `max`: Maximum value in range.
- `min`: Minimum value in range.

# Returns
- `Measurement{Float64}`: Mean value with deviation.

# Category: Misc

"""
function uncertain_from_interval(max, min)
    measurement(mean([max, min]), abs(max - min) / 2)
end

"""
	uncertain_from_percent(val, perc)

Creates a Measurement object from value and percentage uncertainty. 

# Arguments
- `val`: Nominal value.
- `perc`: Percent deviation (from 0 to 100).

# Returns
- `Measurement{Float64}`: Nominal value with deviation.

# Category: Misc

"""
function uncertain_from_percent(val, perc)
    measurement(val, (perc * val) / 100)
end

"""
	compute_ZY_summary(Z_in, Y_in, freq_range, cable_number)

Computes a detailed summary of electrical parameters (impedance, admittance, and their derived quantities) for a given cable system, extracting values from the main diagonal of the impedance (Z) and admittance (Y) matrices.

# Arguments
- `Z_in::Array{Complex{Measurement{Float64}}, 3}`: 3D array representing the impedance matrix over frequency, where the first two dimensions represent the cables and the third dimension represents different frequencies. Each element contains a complex number with uncertainties.
- `Y_in::Array{Complex{Measurement{Float64}}, 3}`: 3D array representing the admittance matrix over frequency, structured similarly to `Z_in`.
- `freq_range::Vector{Float64}`: A vector containing the frequency values corresponding to the third dimension of `Z_in` and `Y_in`.
- `cable_number::Int`: The index of the cable (diagonal element) for which the summary will be computed.

# Returns
- `Z_summary_df::DataFrame`: A DataFrame containing the following columns for each frequency:
	- `frequency`: The frequency value.
	- `Z`: The diagonal element of the impedance matrix for the given cable.
	- `R`: The real part (resistance) of the impedance.
	- `delta_R`: The uncertainty of the resistance.
	- `delta_R_percent`: The percentage uncertainty of the resistance.
	- `L`: The inductance extracted from the imaginary part of the impedance.
	- `delta_L`: The uncertainty of the inductance.
	- `delta_L_percent`: The percentage uncertainty of the inductance.

- `Y_summary_df::DataFrame`: A DataFrame containing the following columns for each frequency:
	- `frequency`: The frequency value.
	- `Y`: The diagonal element of the admittance matrix for the given cable.
	- `C`: The capacitance extracted from the imaginary part of the admittance.
	- `delta_C`: The uncertainty of the capacitance.
	- `delta_C_percent`: The percentage uncertainty of the capacitance.
	- `G`: The real part (conductance) of the admittance.
	- `delta_G`: The uncertainty of the conductance.
	- `delta_G_percent`: The percentage uncertainty of the conductance.

# Notes
- The resistance (R) and inductance (L) are extracted from the real and imaginary parts of `Z`, respectively. The capacitance (C) and conductance (G) are extracted from the imaginary and real parts of `Y`, respectively.
- Uncertainties are calculated for each parameter, and percentage uncertainties are provided to give a relative sense of precision.
- The function uses the `remove_small_values` function to clean up small floating-point values below machine precision (`eps`).

# Category: Cable constants matrices

"""
function compute_ZY_summary(
    Z_in::Array{Complex{Measurement{Float64}},3},
    Y_in::Array{Complex{Measurement{Float64}},3},
    freq_range::Vector{Float64},
    cable_number::Int,
)

    # Create empty DataFrame
    Z_summary_df = DataFrame(
        frequency=Float64[],        # Frequency column
        Z=Complex{Measurement{Float64}}[],   # Z diagonal element for the cable
        R=Float64[],          # Real part (Resistance)
        delta_R=Float64[],    # Uncertainty of Resistance
        delta_R_percent=Float64[],  # % Uncertainty of R
        L=Float64[],          # Inductance
        delta_L=Float64[],    # Uncertainty of L
        delta_L_percent=Float64[],  # % Uncertainty of L
    )

    Y_summary_df = DataFrame(
        frequency=Float64[],        # Frequency column
        Y=Complex{Measurement{Float64}}[],   # Y diagonal element for the cable
        C=Float64[],          # Capacitance
        delta_C=Float64[],    # Uncertainty of C
        delta_C_percent=Float64[],  # % Uncertainty of C
        G=Float64[],          # Conductance
        delta_G=Float64[],    # Uncertainty of G
        delta_G_percent=Float64[],  # % Uncertainty of G
    )

    for k in 1:length(freq_range)
        f = freq_range[k]  # Frequency

        # Get Z_cable and Y_cable for the current frequency
        Z_cable = remove_small_values(Z_in[cable_number, cable_number, k])
        Y_cable = remove_small_values(Y_in[cable_number, cable_number, k])

        # Extract Resistance (R) from Z
        R_value = Measurements.value(real(Z_cable))
        R_uncertainty = Measurements.uncertainty(real(Z_cable))
        R_uncertainty_percent = (R_uncertainty / R_value) * 100

        # Extract Inductance (L) from imag(Z) using the formula: L = imag(Z) / (2 * pi * f)
        L_value = Measurements.value(imag(Z_cable)) / (2 * pi * f)
        L_uncertainty = Measurements.uncertainty(imag(Z_cable)) / (2 * pi * f)
        L_uncertainty_percent = (L_uncertainty / L_value) * 100

        # Extract Conductance (G) from Y
        G_value = Measurements.value(real(Y_cable))
        G_uncertainty = Measurements.uncertainty(real(Y_cable))
        G_uncertainty_percent = (G_uncertainty / G_value) * 100

        # Extract Capacitance (C) from imag(Y) using the formula: C = imag(Y) / (2 * pi * f)
        C_value = Measurements.value(imag(Y_cable)) / (2 * pi * f)
        C_uncertainty = Measurements.uncertainty(imag(Y_cable)) / (2 * pi * f)
        C_uncertainty_percent = (C_uncertainty / C_value) * 100

        # Add row to DataFrame
        push!(
            Z_summary_df,
            (
                f,                # Frequency
                Z_cable,          # Z diagonal element
                R_value,          # Resistance value
                R_uncertainty,    # Resistance uncertainty
                R_uncertainty_percent,  # Resistance % uncertainty
                L_value,          # Inductance value
                L_uncertainty,    # Inductance uncertainty
                L_uncertainty_percent,  # Inductance % uncertainty
            ),
        )

        push!(
            Y_summary_df,
            (
                f,                # Frequency
                Y_cable,          # Y diagonal element
                C_value,          # Capacitance value
                C_uncertainty,    # Capacitance uncertainty
                C_uncertainty_percent,  # Capacitance % uncertainty
                G_value,          # Conductance value
                G_uncertainty,    # Conductance uncertainty
                G_uncertainty_percent,  # Conductance % uncertainty
            ),
        )

    end

    return Z_summary_df, Y_summary_df
end

end  # End of module
