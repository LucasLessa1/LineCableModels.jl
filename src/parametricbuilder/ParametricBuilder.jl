module ParametricBuilder

# Export public API
export make_stranded, make_screened
export Conductor, Insulator, Material
export cardinality, show_trace

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
Material(; rho, eps_r = 1.0, mu_r = 1.0, T0 = 20.0, alpha = 0.0) =
	MaterialSpec(
		rho = _spec(rho),
		eps_r = _spec(eps_r),
		mu_r = _spec(mu_r),
		T0 = _spec(T0),
		alpha = _spec(alpha),
	)

# --- 2) From an existing Material: append %unc by default, or override with (value,pct) ---
function Material(
	m::Materials.Material;
	rho = nothing,
	eps_r = nothing,
	mu_r = nothing,
	T0 = nothing,
	alpha = nothing,
)
	MaterialSpec(
		rho   = _pair_from_nominal(m.rho, rho),
		eps_r = _pair_from_nominal(m.eps_r, eps_r),
		mu_r  = _pair_from_nominal(m.mu_r, mu_r),
		T0    = _pair_from_nominal(m.T0, T0),
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

struct PartChoice
	idx::Int                  # index in ps vector (1-based)
	role::Symbol              # :cond or :insu
	T::Type
	dim::Any                  # chosen scalar (Diameter/Thickness proxy input)
	args::Tuple               # chosen positional args (scalars)
	mat::Materials.Material   # concrete material used
	layers::Int               # n_layers replicated with that choice
end

struct ComponentTrace
	name::String
	choices::Vector{PartChoice}
end

struct DesignTrace
	cable_id::String
	components::Vector{ComponentTrace}
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

	variants = Tuple{DataModel.CableComponent, DataModel.InsulatorGroup, ComponentTrace}[]

	# ---------------- first conductor choice spaces ----------------
	p1c    = cond[1]
	mats1  = _make_range(p1c.material)
	dims1  = _make_range(p1c.dim[1]; pct = p1c.dim[2])
	args1s = collect(_expand_args(p1c.args))  # Vector{<:Tuple}

	# remaining conductors — spaces, with COUPLING flags to p1c
	# Tuple layout: (pc, mcs_or_nothing, dcs_or_nothing, acs_or_nothing)
	rest_cond_spaces = Tuple{PartSpec, Any, Union{Nothing, Any}, Union{Nothing, Any}}[]
	for pc in cond[2:end]
		same_mat  = (pc.material == p1c.material)
		same_dim  = (pc.dim == p1c.dim)
		same_args = (pc.args == p1c.args)

		mcs = same_mat ? nothing : _make_range(pc.material)
		dcs = same_dim ? nothing : _make_range(pc.dim[1]; pct = pc.dim[2])
		acs = same_args ? nothing : collect(_expand_args(pc.args))

		push!(rest_cond_spaces, (pc, mcs, dcs, acs))
	end

	# ---------------- first insulator choice spaces ----------------
	p1i    = insu[1]
	matsi  = _make_range(p1i.material)
	dimsi  = _make_range(p1i.dim[1]; pct = p1i.dim[2])
	args1i = collect(_expand_args(p1i.args))

	# remaining insulators — spaces, with COUPLING flags to p1i
	rest_ins_spaces = Tuple{PartSpec, Any, Union{Nothing, Any}, Union{Nothing, Any}}[]
	for pi in insu[2:end]
		same_mat  = (pi.material == p1i.material)
		same_dim  = (pi.dim == p1i.dim)
		same_args = (pi.args == p1i.args)

		m2 = same_mat ? nothing : _make_range(pi.material)
		d2 = same_dim ? nothing : _make_range(pi.dim[1]; pct = pi.dim[2])
		a2 = same_args ? nothing : collect(_expand_args(pi.args))

		push!(rest_ins_spaces, (pi, m2, d2, a2))
	end

	# ---------------- selection stacks (resolved tuples) -------------
	# chosen_c stores (pc, mc, dc, ac)
	# chosen_i stores (pi, m2i, d2i, a2i)
	chosen_c = Vector{NTuple{4, Any}}()
	chosen_i = Vector{NTuple{4, Any}}()

	# ---------------- build with current resolved choices ------------
	function build_with_current_selection(mat1, d1, a1, mi, di, ai)
		# 1) conductors
		cg = _init_cg(p1c.part_type, base, d1, a1, mat1; abs_first = base == 0.0)
		for k in 2:p1c.n_layers
			_add_conductor!(cg, p1c.part_type, d1, a1, mat1; layer = k)
		end
		for (pc, mc, dc, ac) in chosen_c
			for k in 1:pc.n_layers
				_add_conductor!(cg, pc.part_type, dc, ac, mc; layer = k)
			end
		end

		# 2) insulators
		ig = _init_ig(p1i.part_type, cg, di, ai, mi)
		for (pi, m2i, d2i, a2i) in chosen_i
			for k in 1:pi.n_layers
				_add_insulator!(ig, pi.part_type, d2i, a2i, m2i)
			end
		end

		# assemble trace
		choices = PartChoice[]
		# first conductor spec
		push!(choices, PartChoice(1, :cond, p1c.part_type, d1, a1, mat1, p1c.n_layers))
		# remaining conductors
		for (j, (pc, mc, dc, ac)) in enumerate(chosen_c)
			push!(choices, PartChoice(1 + j, :cond, pc.part_type, dc, ac, mc, pc.n_layers))
		end
		# first insulator spec
		push!(
			choices,
			PartChoice(length(choices)+1, :insu, p1i.part_type, di, ai, mi, p1i.n_layers),
		)
		# remaining insulators
		for (pi, m2i, d2i, a2i) in chosen_i
			push!(
				choices,
				PartChoice(
					length(choices)+1,
					:insu,
					pi.part_type,
					d2i,
					a2i,
					m2i,
					pi.n_layers,
				),
			)
		end
		ctrace = ComponentTrace(String(ps[1].component), choices)

		push!(
			variants,
			(DataModel.CableComponent(String(ps[1].component), cg, ig), ig, ctrace),
		)

		# push!(variants, (DataModel.CableComponent(String(ps[1].component), cg, ig), ig))
	end

	# ---------------- enumerate insulators with coupling -------------
	function choose_ins(idx::Int, mi, di, ai, mat1, d1, a1)
		if idx > length(rest_ins_spaces)
			build_with_current_selection(mat1, d1, a1, mi, di, ai)
			return
		end
		pi, m2, d2, a2 = rest_ins_spaces[idx]

		Ms = (m2 === nothing) ? (mi,) : m2
		Ds = (d2 === nothing) ? (di,) : d2
		As = (a2 === nothing) ? (ai,) : a2

		for m2i in Ms, d2i in Ds, a2i in As
			push!(chosen_i, (pi, m2i, d2i, a2i))
			choose_ins(idx + 1, mi, di, ai, mat1, d1, a1)
			pop!(chosen_i)
		end
	end

	# ---------------- enumerate conductors with coupling -------------
	function choose_cond(idx::Int, mat1, d1, a1, mi, di, ai)
		if idx > length(rest_cond_spaces)
			empty!(chosen_i)
			choose_ins(1, mi, di, ai, mat1, d1, a1)
			return
		end
		pc, mcs, dcs, acs = rest_cond_spaces[idx]

		Ms = (mcs === nothing) ? (mat1,) : mcs
		Ds = (dcs === nothing) ? (d1,) : dcs
		As = (acs === nothing) ? (a1,) : acs

		for mc in Ms, dc in Ds, ac in As
			push!(chosen_c, (pc, mc, dc, ac))
			choose_cond(idx + 1, mat1, d1, a1, mi, di, ai)
			pop!(chosen_c)
		end
	end

	# ---------------- top-level selection loops ----------------------
	for mat1 in mats1, d1 in dims1, a1 in args1s
		for mi in matsi, di in dimsi, ai in args1i
			empty!(chosen_c)
			empty!(chosen_i)
			choose_cond(1, mat1, d1, a1, mi, di, ai)
		end
	end

	return variants
end


function build_designs(cbs::CableBuilderSpec; trace::Bool = false)
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	# partials: (built_components, last_ig_or_nothing)
	partials = Tuple{
		Vector{DataModel.CableComponent},
		Union{Nothing, DataModel.InsulatorGroup},
		Vector{ComponentTrace},
	}[(DataModel.CableComponent[], nothing, ComponentTrace[])]

	for cname in comp_names
		ps = by_comp[cname]
		new_partials = Tuple{
			Vector{DataModel.CableComponent},
			Union{Nothing, DataModel.InsulatorGroup},
			Vector{ComponentTrace},
		}[]
		for (built, last_ig, tr) in partials
			base = last_ig === nothing ? 0.0 : last_ig
			for (comp, ig, ctrace) in _make_variants(ps, base)
				push!(new_partials, (vcat(built, comp), ig, [tr...; ctrace]))
			end
		end
		partials = new_partials
	end



	if !trace
		designs = DataModel.CableDesign[]
		for (comps, _) in ((x[1], x[2]) for x in partials)
			des = DataModel.CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
			for k in Iterators.drop(eachindex(comps), 1)
				add!(des, comps[k])
			end
			push!(designs, des)
		end
		return designs
	else
		out = Tuple{DataModel.CableDesign, DesignTrace}[]
		for (comps, _, ctraces) in partials
			des = DataModel.CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
			for k in Iterators.drop(eachindex(comps), 1)
				add!(des, comps[k])
			end
			push!(out, (des, DesignTrace(cbs.cable_id, ctraces)))
		end
		return out
	end
end

"""
	iterate_designs(cbs) -> Channel{DataModel.CableDesign}

Lazy stream of `CableDesign`s built from `CableBuilderSpec` without allocating all of them.
Works with `for d in iterate_designs(cbs)`.
"""
function iterate_designs(cbs::CableBuilderSpec; trace::Bool = false)
	# group parts by first occurrence of component symbol
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	return Channel(32) do ch
		built = DataModel.CableComponent[]
		lastig = Ref{Union{Nothing, DataModel.InsulatorGroup}}(nothing)
		tr_stack = ComponentTrace[]

		function dfs(idx::Int)
			if idx > length(comp_names)
				des = DataModel.CableDesign(
					cbs.cable_id,
					built[1];
					nominal_data = cbs.nominal,
				)
				for k in Iterators.drop(eachindex(built), 1)
					add!(des, built[k])
				end
				if trace
					put!(ch, (des, DesignTrace(cbs.cable_id, copy(tr_stack))))
				else
					put!(ch, des)
				end
				return
			end
			cname = comp_names[idx];
			ps    = by_comp[cname]
			base  = (lastig[] === nothing) ? 0.0 : lastig[]
			for (comp, ig, ctrace) in _make_variants(ps, base)
				push!(built, comp);
				prev = lastig[];
				lastig[] = ig
				push!(tr_stack, ctrace)
				dfs(idx + 1)
				pop!(tr_stack);
				lastig[] = prev;
				pop!(built)
			end
		end
		dfs(1)
	end
end

function show_trace(tr::DesignTrace)
	println("Design: ", tr.cable_id)
	for comp in tr.components
		println("  Component: ", comp.name)
		for c in comp.choices
			mat = c.mat
			println("    [", c.role, "] ", c.T,
				" layers=", c.layers,
				" dim=", c.dim,
				" args=", c.args,
				" ρ=", mat.rho, " εr=", mat.eps_r, " μr=", mat.mu_r)
		end
	end
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
		dspec = _spec(d)

		# 1) central wire: 1 layer, n=1, lay=0.0
		push!(
			specs,
			PartSpec(component, DataModel.WireArray, 1;
				dim = dspec, args = (1, (0.0, nothing)), material = mat),
		)

		# 2) rings: (layers-1) layers, base n, common lay
		if layers > 1
			push!(
				specs,
				PartSpec(component, DataModel.WireArray, layers - 1;
					dim = dspec, args = (n, _spec(lay)), material = mat),
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

# how many choices are in a "range-like" thing
_choice_count(x) =
	x === nothing                   ? 1 :
	(x isa Tuple && length(x) == 2) ? _choice_count(x[1]) * _choice_count(x[2]) :
	(x isa AbstractVector)          ? length(x) :
	(x isa Tuple && length(x) == 3) ? last(x) : 1

# count choices for a MaterialSpec (rho/eps/mu/T/α product)
_choice_count(ms::MaterialSpec) = length(_make_range(ms))

# args: each entry can be scalar | vector | (lo,hi,n) | (value_spec, pct_spec)
_arg_choice_count(a) =
	(a isa Tuple && length(a) == 2) ? (_choice_count(a[1]) * _choice_count(a[2])) :
	_choice_count(a)

_args_choice_count(args::Tuple) =
	isempty(args) ? 1 : prod(_arg_choice_count(a) for a in args)

function cardinality(cbs::CableBuilderSpec)
	comp_names = unique(p.component for p in cbs.parts)
	by_comp = Dict{Symbol, Vector{PartSpec}}()
	for p in cbs.parts
		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
	end

	total = 1
	for cname in comp_names
		ps = by_comp[cname]
		cond = [p for p in ps if p.part_type <: DataModel.AbstractConductorPart]
		insu = [p for p in ps if p.part_type <: DataModel.AbstractInsulatorPart]
		isempty(cond) && error("component '$cname' has no conductors")
		isempty(insu) && error("component '$cname' has no insulators")

		# first conductor axes
		p1c    = cond[1]
		c_dim  = _choice_count(p1c.dim[1]) * _choice_count(p1c.dim[2])
		c_args = _args_choice_count(p1c.args)
		c_mat  = _choice_count(p1c.material)

		# uncoupled extras from later conductors (couple when tuples compare equal)
		for pc in cond[2:end]
			pc_dim_same  = (pc.dim == p1c.dim)
			pc_args_same = (pc.args == p1c.args)
			pc_mat_same  = (pc.material == p1c.material)

			c_dim  *= pc_dim_same ? 1 : (_choice_count(pc.dim[1]) * _choice_count(pc.dim[2]))
			c_args *= pc_args_same ? 1 : _args_choice_count(pc.args)
			c_mat  *= pc_mat_same ? 1 : _choice_count(pc.material)
		end
		cond_factor = c_dim * c_args * c_mat

		# first insulator axes
		p1i    = insu[1]
		i_dim  = _choice_count(p1i.dim[1]) * _choice_count(p1i.dim[2])
		i_args = _args_choice_count(p1i.args)
		i_mat  = _choice_count(p1i.material)

		for pi in insu[2:end]
			pi_dim_same  = (pi.dim == p1i.dim)
			pi_args_same = (pi.args == p1i.args)
			pi_mat_same  = (pi.material == p1i.material)

			i_dim  *= pi_dim_same ? 1 : (_choice_count(pi.dim[1]) * _choice_count(pi.dim[2]))
			i_args *= pi_args_same ? 1 : _args_choice_count(pi.args)
			i_mat  *= pi_mat_same ? 1 : _choice_count(pi.material)
		end
		insu_factor = i_dim * i_args * i_mat

		total *= cond_factor * insu_factor
	end
	return total
end


end # module ParametricBuilder
