module ParametricBuilder

# Export public API
export make_stranded, make_screened
export Conductor, Insulator, Material

# Module-specific dependencies
using ..Commons
import ..Commons: add!
using ..Materials: Materials
using ..DataModel: DataModel
using Measurements
using Base.Iterators: product

# normalize input to (spec, pct)
_spec(x) = (x isa Tuple && length(x)==2) ? x : (x, nothing)

# Use lib/material nominal; kw is either percent-only or (value,pct)
_pair_from_nominal(nom, x) =
	x === nothing                 ? (nom, nothing) :
	(x isa Tuple && length(x)==2) ? x :
	(nom, x)

# -------------------- material spec --------------------

"""
MaterialSpec: pass specs for fields (value spec + optional %unc)

Example:
  MaterialSpec(; rho=(2.826e-8, nothing),
				 eps_r=(1.0, nothing),
				 mu_r=(1.0, nothing),
				 T0=(20.0, nothing),
				 alpha=(4.0e-3, nothing))
"""
struct MaterialSpec
	rho::Any;
	eps_r::Any;
	mu_r::Any;
	T0::Any;
	alpha::Any
end
MaterialSpec(; rho, eps_r, mu_r, T0, alpha) = MaterialSpec(rho, eps_r, mu_r, T0, alpha)

# --- 1) Ad-hoc numeric: values (or (value,pct)) ---
Material(; rho, eps = 1.0, mu = 1.0, T = 20.0, alpha = 0.0) =
	MaterialSpec(
		rho = _spec(rho),
		eps_r = _spec(eps),
		mu_r = _spec(mu),
		T0 = _spec(T),
		alpha = _spec(alpha),
	)

# --- 2) From an existing Material: append %unc by default, or override with (value,pct) ---
function Material(
	m::Materials.Material;
	rho = nothing,
	eps = nothing,
	mu = nothing,
	T = nothing,
	alpha = nothing,
)
	MaterialSpec(
		rho   = _pair_from_nominal(m.rho, rho),
		eps_r = _pair_from_nominal(m.eps_r, eps),
		mu_r  = _pair_from_nominal(m.mu_r, mu),
		T0    = _pair_from_nominal(m.T0, T),
		alpha = _pair_from_nominal(m.alpha, alpha),
	)
end

# --- 3) From a MaterialsLibrary + name ---
Material(lib::Materials.MaterialsLibrary, name::AbstractString; kwargs...) =
	Material(get(lib, name); kwargs...)
Material(lib::Materials.MaterialsLibrary, name::Symbol; kwargs...) =
	Material(lib, String(name); kwargs...)


function _make_range(ms::MaterialSpec)
	ρs = _make_range(ms.rho[1]; pct = ms.rho[2])
	εs = _make_range(ms.eps_r[1]; pct = ms.eps_r[2])
	μs = _make_range(ms.mu_r[1]; pct = ms.mu_r[2])
	Ts  = _make_range(ms.T0[1]; pct = ms.T0[2])
	αs = _make_range(ms.alpha[1]; pct = ms.alpha[2])
	[Materials.Material(ρ, ε, μ, T, α) for (ρ, ε, μ, T, α) in product(ρs, εs, μs, Ts, αs)]
end


# spec: (value_spec, pct_spec) — pct_spec can be `nothing | number | vector | (lo,hi,n)`
"""
PartSpec:
- component::Symbol           # e.g. :core, :sheath, :jacket
- part_type::Type             # WireArray, Tubular, Strip, Insulator, Semicon, …
- n_layers::Int               # how many stacked layers of this part_type
- dim::Tuple                  # diameter OR thickness OR radius (spec, pct)
- args::Tuple                 # specialized ctor positional args
- material::MaterialSpec
"""
struct PartSpec
	component::Symbol
	part_type::Type
	n_layers::Int
	dim::Tuple # (spec, pct)
	args::Tuple # positional args; each entry is either a number or (spec, pct)
	material::MaterialSpec
end

PartSpec(component::Symbol, part_type::Type, n_layers::Int;
	dim, args = (), material::MaterialSpec) =
	PartSpec(component, part_type, n_layers, dim, args, material)

"""
CableBuilderSpec:
- cable_id::String
- parts::Vector{PartSpec}  # may interleave conductor/insulator arbitrarily
- nominal::Union{Nothing,DataModel.NominalData}
"""
struct CableBuilderSpec
	cable_id::String
	parts::Vector{PartSpec}
	nominal::Union{Nothing, DataModel.NominalData}
end
CableBuilderSpec(id::AbstractString, parts::Vector{PartSpec}; nominal = nothing) =
	CableBuilderSpec(String(id), parts, nominal)

# --- minimal flattening helpers (accept PartSpec or collections of them) -----
function _collect_parts!(acc::Vector{PartSpec}, x)
	if x isa PartSpec
		push!(acc, x)
	elseif x isa AbstractVector
		@inbounds for y in x
			_collect_parts!(acc, y)
		end
	else
		error("Expected PartSpec or a collection of PartSpec; got $(typeof(x))")
	end
	return acc
end

# ctor that accepts a vector with possible nested vectors (no splat needed)
function CableBuilderSpec(id::AbstractString, parts_any::AbstractVector; nominal = nothing)
	acc = PartSpec[]
	_collect_parts!(acc, parts_any)
	return CableBuilderSpec(String(id), acc, nominal)  # calls your primary ctor
end

# ctor that accepts varargs (mixed PartSpec and vectors), plus nominal kw
function CableBuilderSpec(id::AbstractString, parts...; nominal = nothing)
	acc = PartSpec[]
	@inbounds for p in parts
		_collect_parts!(acc, p)
	end
	return CableBuilderSpec(String(id), acc, nominal)
end

# ----- realizers -----
_values(x::Number) = (x,)
_values(v::AbstractVector) = collect(v)
_values(t::Tuple{<:Number, <:Number, <:Integer}) = range(t[1], t[2]; length = t[3])

_pcts(::Nothing) = (0.0,)
_pcts(p::Number) = (float(p),)
_pcts(v::AbstractVector) = map(float, collect(v))
_pcts(t::Tuple{<:Number, <:Number, <:Integer}) =
	range(float(t[1]), float(t[2]); length = t[3])

function _make_range(spec; pct = nothing)
	vs, ps = collect(_values(spec)), collect(_pcts(pct))
	if all(p->p==0.0, ps)
		;
		return vs;
	end
	out = Any[]
	for v in vs, p in ps
		push!(out, measurement(v, abs(v)*(p/100)))
	end
	out
end

# expand positional args tuple → iterator of resolved tuples
function _expand_args(args::Tuple)
	spaces =
		map(a -> (a isa Tuple && length(a)==2 ? _make_range(a[1]; pct = a[2]) : (a,)), args)
	return (tuple(vals...) for vals in Iterators.product(spaces...))
end

# ----- anchor: last physical layer, not the container -----
@inline _anchor(x::Real) = x
@inline _anchor(x::DataModel.AbstractConductorPart) = x
@inline _anchor(x::DataModel.AbstractInsulatorPart) = x

@inline function _anchor(g::DataModel.ConductorGroup)
	L = getfield(g, :layers)
	@assert !isempty(L) "ConductorGroup has no layers to anchor on."
	return L[end]
end
@inline function _anchor(g::DataModel.InsulatorGroup)
	L = getfield(g, :layers)
	@assert !isempty(L) "InsulatorGroup has no layers to anchor on."
	return L[end]
end

# ----- proxy for radius_ext by CONTRACT -----
@inline function _resolve_dim(T::Type, is_abs_first::Bool)
	return T <: DataModel.AbstractWireArray ? :diameter :
		   (is_abs_first && T === DataModel.Tubular ? :diameter : :thickness)
end

@inline _make_dim(::Val{:diameter}, d)  = DataModel.Diameter(d)
@inline _make_dim(::Val{:thickness}, d) = DataModel.Thickness(d)
@inline _make_dim(::Val{:radius}, r)    = r  # if direct radius
@inline _make_dim(sym::Symbol, d)       = _make_dim(Val(sym), d)


function _init_cg(T::Type, base, dim_val, args_pos::Tuple, mat; abs_first::Bool)

	r_in = _anchor(base)
	sym  = _resolve_dim(T, abs_first)
	_    = _make_dim(sym, dim_val)   # keeps intent (WireArray ignores this)

	if T <: DataModel.AbstractWireArray
		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
		n   = args_pos[1]
		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
		return DataModel.ConductorGroup(
			DataModel.WireArray(r_in, DataModel.Diameter(dim_val), n, lay, mat),
		)
	else
		return DataModel.ConductorGroup(T(r_in, _make_dim(sym, dim_val), args_pos..., mat))
	end
end


function _add_conductor!(
	cg::DataModel.ConductorGroup,
	T::Type,
	dim_val,
	args_pos::Tuple,
	mat;
	layer::Int,
)
	if T <: DataModel.AbstractWireArray
		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
		n   = args_pos[1]
		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
		add!(cg, DataModel.WireArray, DataModel.Diameter(dim_val), layer*n, lay, mat)
	else
		# thickness by contract for all non-wire additions
		add!(cg, T, _make_dim(:thickness, dim_val), args_pos..., mat)
	end
end


function _init_ig(T::Type, cg::DataModel.ConductorGroup, dim_val, args_pos::Tuple, mat)
	c_last = _anchor(cg)
	obj = T(c_last, DataModel.Thickness(dim_val), args_pos..., mat)  # insulators use THICKNESS
	return DataModel.InsulatorGroup(obj)
end


function _add_insulator!(
	ig::DataModel.InsulatorGroup,
	T::Type,
	dim_val,
	args_pos::Tuple,
	mat,
)
	add!(ig, T, DataModel.Thickness(dim_val), args_pos..., mat)
end


# Build all variants of ONE component, anchored at `base` (0.0 for the very first)
function _make_variants(ps::Vector{PartSpec}, base)
	cond = [p for p in ps if p.part_type <: DataModel.AbstractConductorPart]
	insu = [p for p in ps if p.part_type <: DataModel.AbstractInsulatorPart]
	isempty(cond) && error("component has no conductors")
	isempty(insu) && error("component has no insulators")

	variants = Tuple{DataModel.CableComponent, DataModel.InsulatorGroup}[]

	# --- conductors ---
	p1c   = cond[1]
	mats1 = _make_range(p1c.material)
	dims1 = _make_range(p1c.dim[1]; pct = p1c.dim[2])

	for mat1 in mats1, d1 in dims1, a1 in _expand_args(p1c.args)
		cg = _init_cg(p1c.part_type, base, d1, a1, mat1; abs_first = base == 0.0)
		for k ∈ 2:p1c.n_layers, a in _expand_args(p1c.args)
			_add_conductor!(cg, p1c.part_type, d1, a, mat1; layer = k)
		end
		# remaining conductor specs
		for pc in cond[2:end]
			matsc = _make_range(pc.material);
			dimsc = _make_range(pc.dim[1]; pct = pc.dim[2])
			for mc in matsc,
				dc in dimsc,
				a in _expand_args(pc.args),
				k ∈ 1:pc.n_layers

				_add_conductor!(cg, pc.part_type, dc, a, mc; layer = k)
			end
		end

		# --- insulators (wrap around the conductor group) ---
		p1i   = insu[1]
		matsi = _make_range(p1i.material);
		dimsi = _make_range(p1i.dim[1]; pct = p1i.dim[2])
		for mi in matsi, di in dimsi, ai in _expand_args(p1i.args)
			ig = _init_ig(p1i.part_type, cg, di, ai, mi)
			for pi in insu[2:end]
				m2 = _make_range(pi.material);
				d2 = _make_range(pi.dim[1]; pct = pi.dim[2])
				for m2i in m2,
					d2i in d2,
					a2 in _expand_args(pi.args),
					k ∈ 1:pi.n_layers

					_add_insulator!(ig, pi.part_type, d2i, a2, m2i)
				end
			end
			push!(variants, (DataModel.CableComponent(String(ps[1].component), cg, ig), ig))
		end
	end

	return variants
end

function build_designs(cbs::CableBuilderSpec)
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	# partials: (built_components, last_ig_or_nothing)
	partials =
		Tuple{Vector{DataModel.CableComponent}, Union{Nothing, DataModel.InsulatorGroup}}[
			(Vector{DataModel.CableComponent}(), nothing)
		]

	for cname in comp_names
		ps = by_comp[cname]
		new_partials = Tuple{
			Vector{DataModel.CableComponent},
			Union{Nothing, DataModel.InsulatorGroup},
		}[]
		for (built, last_ig) in partials
			base = last_ig === nothing ? 0.0 : last_ig  # will be resolved to last layer where needed
			for (comp, ig) in _make_variants(ps, base)
				push!(new_partials, (vcat(built, comp), ig))
			end
		end
		partials = new_partials
	end

	designs = DataModel.CableDesign[]
	for (comps, _) in partials
		des = DataModel.CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
		for k in Iterators.drop(eachindex(comps), 1)
			add!(des, comps[k])
		end
		push!(designs, des)
	end
	return designs
end

module Conductor

	using ..ParametricBuilder: PartSpec, _spec
	using ...DataModel: DataModel

	# wire: args are (n, lay)
	Wires(component::Symbol; layers::Int, d, n::Int, lay = 11.0, mat) =
		PartSpec(component, DataModel.WireArray, layers;
			dim = _spec(d), args = (n, _spec(lay)), material = mat)

	# tube: no extra args
	Tubular(component::Symbol; layers::Int, t, mat) =
		PartSpec(component, DataModel.Tubular, layers;
			dim = _spec(t), args = (), material = mat)

	# strip: args are (width, lay)
	Strip(component::Symbol; layers::Int, t, w, lay = 0.0, mat) =
		PartSpec(component, DataModel.Strip, layers;
			dim = _spec(t), args = (_spec(w), _spec(lay)), material = mat)

	# central + hex rings sugar
	function Stranded(component::Symbol; layers::Int, d, n::Int, lay = 11.0, mat)
		@assert layers >= 1 "stranded: layers must be ≥ 1 (includes the central wire)."
		specs = PartSpec[]

		# 1) central wire: 1 layer, n=1, lay=0.0
		push!(
			specs,
			PartSpec(component, DataModel.WireArray, 1;
				dim = _spec(d), args = (1, (0.0, nothing)), material = mat),
		)

		# 2) rings: (layers-1) layers, base n, common lay
		if layers > 1
			push!(
				specs,
				PartSpec(component, DataModel.WireArray, layers - 1;
					dim = _spec(d), args = (n, _spec(lay)), material = mat),
			)
		end

		return specs
	end

end

module Insulator

	using ..ParametricBuilder: PartSpec, _spec
	using ...DataModel: DataModel

	Tubular(component::Symbol; layers::Int, t, mat) =
		PartSpec(component, DataModel.Insulator, layers;
			dim = _spec(t), args = (), material = mat)

	Semicon(component::Symbol; layers::Int, t, mat) =
		PartSpec(component, DataModel.Semicon, layers;
			dim = _spec(t), args = (), material = mat)
end


# Submodule `WirePatterns`
include("wirepatterns/WirePatterns.jl")
using .WirePatterns

end # module ParametricBuilder
