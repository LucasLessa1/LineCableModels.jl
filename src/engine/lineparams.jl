struct SeriesImpedance{T}
    values::Array{T,3}   # n×n×nfreq, units: Ω/m
end

struct ShuntAdmittance{T}
    values::Array{T,3}   # n×n×nfreq, units: S/m
end

"""
$(TYPEDEF)

Represents the frequency-dependent line parameters (series impedance and shunt admittance matrices) for a cable or line system.

$(TYPEDFIELDS)
"""
struct LineParameters{T<:COMPLEXSCALAR}
    "Series impedance matrices \\[Ω/m\\]."
    Z::SeriesImpedance{T}
    "Shunt admittance matrices \\[S/m\\]."
    Y::ShuntAdmittance{T}

    @doc """
    $(TYPEDSIGNATURES)

    Constructs a [`LineParameters`](@ref) instance.

    # Arguments

    - `Z`: Series impedance matrices \\[Ω/m\\].
    - `Y`: Shunt admittance matrices \\[S/m\\].

    # Returns

    - A [`LineParameters`](@ref) object with prelocated impedance and admittance matrices.

    # Examples

    ```julia
    params = $(FUNCTIONNAME)(Z, Y)
    ```
    """
    function LineParameters(Z::SeriesImpedance{T}, Y::ShuntAdmittance{T}) where {T<:COMPLEXSCALAR}
        size(Z, 1) == size(Z, 2) || throw(DimensionMismatch("Z must be square"))
        size(Y, 1) == size(Y, 2) || throw(DimensionMismatch("Y must be square"))
        size(Z) == size(Y) || throw(DimensionMismatch("Z and Y must have same dimensions (n×n×nfreq)"))
        new{T}(Z, Y)
    end
end

SeriesImpedance(A::AbstractArray{T,3}) where {T} = SeriesImpedance{T}(Array(A))
ShuntAdmittance(A::AbstractArray{T,3}) where {T} = ShuntAdmittance{T}(Array(A))

# Convenience ctor for LineParameters from raw arrays
function LineParameters(Z::Array{T,3}, Y::Array{T,3}) where {T<:COMPLEXSCALAR}
    size(Z) == size(Y) || throw(DimensionMismatch("Z and Y must have same dims"))
    return LineParameters(SeriesImpedance{T}(Z), ShuntAdmittance{T}(Y))
end