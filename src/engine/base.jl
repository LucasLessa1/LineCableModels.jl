
Base.eltype(::LineParametersProblem{T}) where {T} = T
Base.eltype(::Type{LineParametersProblem{T}}) where {T} = T

Base.eltype(::LineParameters{T}) where {T} = T
Base.eltype(::Type{LineParameters{T}}) where {T} = T

Base.eltype(::CoaxialWorkspace{T}) where {T} = T
Base.eltype(::Type{CoaxialWorkspace{T}}) where {T} = T


abstract type UnitLen end
struct PerMeter <: UnitLen end
struct PerKilometer <: UnitLen end
_len_scale(::PerMeter) = 1.0
_len_scale(::PerKilometer) = 1_000.0
_len_label(::PerMeter) = "m"
_len_label(::PerKilometer) = "km"

abstract type DisplayMode end
struct AsZY <: DisplayMode end
struct AsRLCG <: DisplayMode
    ω::Float64  # angular frequency for L,C extraction
end

using Printf
using Measurements: value, uncertainty

"""
ResultsView: pretty, *non-mutating* renderer with zero-clipping and units.

- `mode = AsZY()` prints Z [Ω/len], Y [S/len].
- `mode = AsRLCG(ω)` prints R [Ω/len], L [mH/len], G [S/len], C [µF/len].
- `tol` clips tiny magnitudes to 0.0 in display only (value & uncertainty).
"""
struct ResultsView{LP<:LineParameters,U<:UnitLen,M<:DisplayMode}
    lp::LP
    unit::U
    mode::M
    tol::Float64
end

# Builders
resultsview(lp::LineParameters; per::Symbol=:km, mode::Symbol=:ZY, ω::Union{Nothing,Real}=nothing, tol::Real=sqrt(eps(Float64))) =
    ResultsView(
        lp,
        per === :km ? PerKilometer() : PerMeter(),
        mode === :ZY ? AsZY() : AsRLCG(ω === nothing ? throw(ArgumentError("ω required for mode=:RLCG")) : float(ω)),
        float(tol)
    )

per_km(lp::LineParameters; mode::Symbol=:ZY, ω::Union{Nothing,Real}=nothing, tol::Real=sqrt(eps(Float64))) =
    resultsview(lp; per=:km, mode=mode, ω=ω, tol=tol)

per_m(lp::LineParameters; mode::Symbol=:ZY, ω::Union{Nothing,Real}=nothing, tol::Real=sqrt(eps(Float64))) =
    resultsview(lp; per=:m, mode=mode, ω=ω, tol=tol)

# --- Scalar formatting with zero-clipping -------------------------------------

# Zero-clip helpers (display only)
_clip(x::Real, tol) = (abs(x) < tol ? 0.0 : x)

_format_real(io, x::Real, tol) = @printf(io, "%.6g", _clip(x, tol))

_format_meas(io, m, tol) = begin
    v = _clip(value(m), tol)
    u = _clip(uncertainty(m), tol)
    @printf(io, "%.6g±%.6g", v, u)
end

_format_complex(io, z, tol) = begin
    # z may be Complex{<:Real} or Complex{<:Measurement}
    print(io, "")
    if z.re isa Real
        _format_real(io, real(z), tol)
    else
        _format_meas(io, real(z), tol)
    end
    print(io, "+")
    if z.im isa Real
        _format_real(io, imag(z), tol)
    else
        _format_meas(io, imag(z), tol)
    end
    print(io, "im")
end

_format_any(io, x, tol) = x isa Complex ? _format_complex(io, x, tol) :
                          x isa Measurements.Measurement ? _format_meas(io, x, tol) :
                          _format_real(io, x, tol)



# --- Show methods --------------------------------------------------------------

function _show_matrix(io::IO, A::AbstractArray; tol::Float64, map::Function=identity)
    n1, n2 = size(A, 1), size(A, 2)
    for i in 1:n1
        for j in 1:n2
            j > 1 && print(io, "  ")
            _format_any(io, map(A[i, j]), tol)
        end
        i < n1 && print(io, '\n')
    end
end

function Base.show(io::IO, ::MIME"text/plain", rv::ResultsView)
    lp = rv.lp
    unit = rv.unit
    tol = rv.tol
    scale = _len_scale(unit)
    ulabel = _len_label(unit)
    _, _, nf = size(lp.Z)

    print(io, "LineParameters ResultsView  |  mode = ")
    print(io, rv.mode isa AsZY ? "ZY" : @sprintf("RLCG @ ω=%.6g", rv.mode.ω))
    print(io, "  |  units per ", ulabel, "  |  tol = ", tol, "\n")

    @views for k in 1:nf
        print(io, "\n[:, :, ", k, "]\n")
        Zk = lp.Z.values[:, :, k]
        Yk = lp.Y.values[:, :, k]

        if rv.mode isa AsZY
            println(io, "Z [Ω/", ulabel, "] =")
            _show_matrix(io, Zk; tol=tol, map=x -> scale * x)

            print(io, "\n\nY [S/", ulabel, "] =\n")
            _show_matrix(io, Yk; tol=tol, map=x -> scale * x)
        else
            ω = rv.mode.ω

            println(io, "R [Ω/", ulabel, "] =")
            _show_matrix(io, Zk; tol=tol, map=x -> scale * real(x))

            print(io, "\n\nL [mH/", ulabel, "] =\n")
            _show_matrix(io, Zk; tol=tol, map=x -> (scale * 1e3 / ω) * imag(x))

            print(io, "\n\nG [S/", ulabel, "] =\n")
            _show_matrix(io, Yk; tol=tol, map=x -> scale * real(x))

            print(io, "\n\nC [µF/", ulabel, "] =\n")
            _show_matrix(io, Yk; tol=tol, map=x -> (scale * 1e6 / ω) * imag(x))
        end

        k < nf && print(io, "\n", "---"^10, "\n")
    end
end

function Base.show(io::IO, ::MIME"text/plain", Z::SeriesImpedance)
    println(io, "SeriesImpedance [Ω/m], n×n×nfreq = ", size(Z.values))
    @views _show_matrix(io, Z.values[:, :, 1]; tol=sqrt(eps(Float64)))
    size(Z, 3) > 1 && print(io, "\n… (", size(Z, 3) - 1, " more slice", size(Z, 3) - 1 == 1 ? "" : "s", ")")
end

function Base.show(io::IO, ::MIME"text/plain", Y::ShuntAdmittance)
    println(io, "ShuntAdmittance [S/m], n×n×nfreq = ", size(Y.values))
    @views _show_matrix(io, Y.values[:, :, 1]; tol=sqrt(eps(Float64)))
    size(Y, 3) > 1 && print(io, "\n… (", size(Y, 3) - 1, " more slice", size(Y, 3) - 1 == 1 ? "" : "s", ")")
end



# ---- SeriesImpedance array-ish interface ----
Base.size(Z::SeriesImpedance) = size(Z.values)
Base.size(Z::SeriesImpedance, d::Int) = size(Z.values, d)
Base.axes(Z::SeriesImpedance) = axes(Z.values)
Base.ndims(::Type{SeriesImpedance{T}}) where {T} = 3
Base.eltype(::Type{SeriesImpedance{T}}) where {T} = T
Base.getindex(Z::SeriesImpedance, I...) = @inbounds Z.values[I...]

# ---- ShuntAdmittance array-ish interface ----
Base.size(Y::ShuntAdmittance) = size(Y.values)
Base.size(Y::ShuntAdmittance, d::Int) = size(Y.values, d)
Base.axes(Y::ShuntAdmittance) = axes(Y.values)
Base.ndims(::Type{ShuntAdmittance{T}}) where {T} = 3
Base.eltype(::Type{ShuntAdmittance{T}}) where {T} = T
Base.getindex(Y::ShuntAdmittance, I...) = @inbounds Y.values[I...]

# --- Frequency-slice sugar ----------------------------------------------------
@inline Base.getindex(lp::LineParameters, k::Integer) = LineParameters(
    SeriesImpedance(@view lp.Z.values[:, :, k:k]),
    ShuntAdmittance(@view lp.Y.values[:, :, k:k])
)

# --- One-argument k, derive ω from freq (or accept ω directly) ---------------
function per_km(lp::LineParameters, k::Integer;
    mode::Symbol=:ZY,
    freq::AbstractVector{<:Real}=Float64[],
    ω=nothing,
    tol::Real=sqrt(eps(Float64)))
    lpk = lp[k]
    if mode === :ZY
        return resultsview(lpk; per=:km, mode=:ZY, tol=tol)
    elseif mode === :RLCG
        ωv = ω === nothing ? (isempty(freq) ? throw(ArgumentError("Provide ω or freq")) : 2π * float(freq[k])) : float(ω)
        return resultsview(lpk; per=:km, mode=:RLCG, ω=ωv, tol=tol)
    else
        throw(ArgumentError("mode must be :ZY or :RLCG"))
    end
end

function per_m(lp::LineParameters, k::Integer;
    mode::Symbol=:ZY,
    freq::AbstractVector{<:Real}=Float64[],
    ω=nothing,
    tol::Real=sqrt(eps(Float64)))
    lpk = lp[k]
    if mode === :ZY
        return resultsview(lpk; per=:m, mode=:ZY, tol=tol)
    elseif mode === :RLCG
        ωv = ω === nothing ? (isempty(freq) ? throw(ArgumentError("Provide ω or freq")) : 2π * float(freq[k])) : float(ω)
        return resultsview(lpk; per=:m, mode=:RLCG, ω=ωv, tol=tol)
    else
        throw(ArgumentError("mode must be :ZY or :RLCG"))
    end
end

# Helper: detect uncertainties in element type
_has_uncertainty_type(::Type{Complex{S}}) where {S} = S <: Measurement
_has_uncertainty_type(::Type) = false

# Terse summary (used inside collections)
function Base.show(io::IO, lp::LineParameters)
    n, _, nf = size(lp.Z)
    T = eltype(lp.Z)
    print(io, "LineParameters{$(T)} ", n, "×", n, "×", nf, "  [Z:Ω/m, Y:S/m]")
    _has_uncertainty_type(T) && print(io, " (±)")
end

# Pretty REPL display
function Base.show(io::IO, ::MIME"text/plain", lp::LineParameters)
    n, _, nf = size(lp.Z)
    T = eltype(lp.Z)
    tol = sqrt(eps(Float64))
    scale = 1_000.0   # per km preview
    ulabel = "km"

    print(io, "LineParameters  |  n = ", n, "  |  nf = ", nf,
        "  |  eltype = ", T)
    _has_uncertainty_type(T) && print(io, "  |  uncertainties: yes")
    print(io, "\n")

    # Preview: slice 1, per km, Z then Y
    @views begin
        Z1 = view(lp.Z.values, :, :, 1)
        Y1 = view(lp.Y.values, :, :, 1)

        println(io, "\nPreview (slice 1/", nf, ")  per ", ulabel)
        println(io, "Z [Ω/", ulabel, "] =")
        _show_matrix(io, Z1; tol=tol, map=x -> scale * x)

        print(io, "\n\nY [S/", ulabel, "] =\n")
        _show_matrix(io, Y1; tol=tol, map=x -> scale * x)
    end

    if nf > 1
        print(io, "\n\n… (", nf - 1, " more frequency slice", nf - 1 == 1 ? "" : "s", ")")
    end
end