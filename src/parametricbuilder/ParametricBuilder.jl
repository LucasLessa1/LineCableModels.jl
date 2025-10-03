module ParametricBuilder

# Export public API
export make_stranded, make_screened
export conductor, insulator, material

# Module-specific dependencies
using ..Commons
import ..Commons: add!
using ..Materials: Material
using ..DataModel: NominalData, AbstractConductorPart, AbstractInsulatorPart,
	AbstractWireArray, ConductorGroup, InsulatorGroup, CableComponent, CableDesign,
	Thickness, Diameter, WireArray, Tubular, Strip, Insulator, Semicon
using Measurements
using Base.Iterators: product

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
material(; rho, eps = 1.0, mu = 1.0, T = 20.0, α = 0.0) =
	MaterialSpec(rho = rho, eps_r = (eps, nothing), mu_r = (mu, nothing),
		T0 = (T, nothing), alpha = (α, nothing))

function _make_range(ms::MaterialSpec)
	ρs = _make_range(ms.rho[1]; pct = ms.rho[2])
	εs = _make_range(ms.eps_r[1]; pct = ms.eps_r[2])
	μs = _make_range(ms.mu_r[1]; pct = ms.mu_r[2])
	Ts  = _make_range(ms.T0[1]; pct = ms.T0[2])
	αs = _make_range(ms.alpha[1]; pct = ms.alpha[2])
	[Material(ρ, ε, μ, T, α) for (ρ, ε, μ, T, α) in product(ρs, εs, μs, Ts, αs)]
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
	dim::Tuple                      # (spec, pct)
	args::Tuple                     # positional args; each entry is either a number or (spec, pct)
	material::MaterialSpec
end

PartSpec(component::Symbol, part_type::Type, n_layers::Int;
	dim, args = (), material::MaterialSpec) =
	PartSpec(component, part_type, n_layers, dim, args, material)

"""
CableBuilderSpec:
- cable_id::String
- parts::Vector{PartSpec}  # may interleave conductor/insulator arbitrarily
- nominal::Union{Nothing,NominalData}
"""
struct CableBuilderSpec
	cable_id::String
	parts::Vector{PartSpec}
	nominal::Union{Nothing, NominalData}
end
CableBuilderSpec(id::AbstractString, parts::Vector{PartSpec}; nominal = nothing) =
	CableBuilderSpec(String(id), parts, nominal)

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
@inline _anchor(x::AbstractConductorPart) = x
@inline _anchor(x::AbstractInsulatorPart) = x

@inline function _anchor(g::ConductorGroup)
	L = getfield(g, :layers)
	@assert !isempty(L) "ConductorGroup has no layers to anchor on."
	return L[end]
end
@inline function _anchor(g::InsulatorGroup)
	L = getfield(g, :layers)
	@assert !isempty(L) "InsulatorGroup has no layers to anchor on."
	return L[end]
end

# ----- proxy for radius_ext by CONTRACT -----
@inline function _which_dim(T::Type, is_abs_first::Bool)
	return T <: AbstractWireArray ? :diameter :
		   (is_abs_first && T === Tubular ? :diameter : :thickness)
end

@inline _make_dim(::Val{:diameter}, d)  = Diameter(d)
@inline _make_dim(::Val{:thickness}, d) = Thickness(d)
@inline _make_dim(::Val{:radius}, r)    = r  # if you ever feed a direct radius
@inline _make_dim(sym::Symbol, d)       = _make_dim(Val(sym), d)


function _init_cg(T::Type, base, dim_val, args_pos::Tuple, mat; abs_first::Bool)

	r_in = _anchor(base)
	sym  = _which_dim(T, abs_first)
	_    = _make_dim(sym, dim_val)   # keeps intent (WireArray ignores this)

	if T <: AbstractWireArray
		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
		n   = args_pos[1]
		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
		return ConductorGroup(WireArray(r_in, Diameter(dim_val), n, lay, mat))
	else
		return ConductorGroup(T(r_in, _make_dim(sym, dim_val), args_pos..., mat))
	end
end


function _add_conductor!(
	cg::ConductorGroup,
	T::Type,
	dim_val,
	args_pos::Tuple,
	mat;
	layer::Int,
)
	if T <: AbstractWireArray
		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
		n   = args_pos[1]
		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
		add!(cg, WireArray, Diameter(dim_val), layer*n, lay, mat)
	else
		# thickness by contract for all non-wire additions
		add!(cg, T, _make_dim(:thickness, dim_val), args_pos..., mat)
	end
end


function _init_ig(T::Type, cg::ConductorGroup, dim_val, args_pos::Tuple, mat)
	c_last = _anchor(cg)
	obj = T(c_last, Thickness(dim_val), args_pos..., mat)  # insulators use THICKNESS
	return InsulatorGroup(obj)
end


function _add_insulator!(
	ig::InsulatorGroup,
	T::Type,
	dim_val,
	args_pos::Tuple,
	mat,
)
	add!(ig, T, Thickness(dim_val), args_pos..., mat)
end


# # seed a ConductorGroup for the FIRST conductor spec of a component
# function _init_cg(T::Type, base, dim_val, args_pos::NTuple, mat; abs_first::Bool)
# 	r_in = _anchor(base)
# 	sym = _which_dim(T, abs_first)
# 	r_ex = _make_dim(sym, dim_val)
# 	if T <: AbstractWireArray
# 		@assert length(args_pos) ≥ 1 "WireArray needs (n, [lay])."
# 		n   = args_pos[1]
# 		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
# 		return ConductorGroup(WireArray(r_in, Diameter(dim_val), n, lay, mat))
# 	else
# 		return ConductorGroup(T(r_in, r_ex, args_pos..., mat))
# 	end
# end

# # stack one conductor layer into an existing cg
# function _add_conductor!(
# 	cg::ConductorGroup,
# 	T::Type,
# 	dim_val,
# 	args_pos::NTuple,
# 	mat;
# 	layer::Int,
# )
# 	if T <: AbstractWireArray
# 		@assert length(args_pos) ≥ 1
# 		n   = args_pos[1]
# 		lay = length(args_pos) ≥ 2 ? args_pos[2] : 0.0
# 		add!(cg, WireArray, Diameter(dim_val), layer*n, lay, mat)
# 	else
# 		# thickness by contract for all non-wire additions
# 		add!(cg, T, _make_dim(:thickness, dim_val), args_pos..., mat)
# 	end
# end

# # seed an InsulatorGroup wrapping a conductor group `cg`
# function _init_ig(T::Type, cg::ConductorGroup, dim_val, args_pos::NTuple, mat)
# 	c_last = _anchor(cg)
# 	# insulators use THICKNESS by contract
# 	obj = T(c_last, Thickness(dim_val), args_pos..., mat)
# 	return InsulatorGroup(obj)
# end

# # stack one insulator layer
# @inline function _add_insulator!(
# 	ig::InsulatorGroup,
# 	T::Type,
# 	dim_val,
# 	args_pos::NTuple,
# 	mat,
# )
# 	add!(ig, T, Thickness(dim_val), args_pos..., mat)
# end

# Build all variants of ONE component, anchored at `base` (0.0 for the very first)
function _make_variants(ps::Vector{PartSpec}, base)
	cond = [p for p in ps if p.part_type <: AbstractConductorPart]
	insu = [p for p in ps if p.part_type <: AbstractInsulatorPart]
	isempty(cond) && error("component has no conductors")
	isempty(insu) && error("component has no insulators")

	variants = Tuple{CableComponent, InsulatorGroup}[]

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
			push!(variants, (CableComponent(String(ps[1].component), cg, ig), ig))
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
	partials = Tuple{Vector{CableComponent}, Union{Nothing, InsulatorGroup}}[
		(Vector{CableComponent}(), nothing)
	]

	for cname in comp_names
		ps = by_comp[cname]
		new_partials = Tuple{Vector{CableComponent}, Union{Nothing, InsulatorGroup}}[]
		for (built, last_ig) in partials
			base = last_ig === nothing ? 0.0 : last_ig  # will be resolved to last layer where needed
			for (comp, ig) in _make_variants(ps, base)
				push!(new_partials, (vcat(built, comp), ig))
			end
		end
		partials = new_partials
	end

	designs = CableDesign[]
	for (comps, _) in partials
		des = CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
		for k in Iterators.drop(eachindex(comps), 1)
			add!(des, comps[k])
		end
		push!(designs, des)
	end
	return designs
end

# normalize input to (spec, pct)
_spec(x) = (x isa Tuple && length(x)==2) ? x : (x, nothing)

module conductor

	using ..ParametricBuilder: PartSpec, _spec
	using ...DataModel: WireArray, Tubular, Strip

	# wire: args are (n, lay)
	wire(component::Symbol; layers::Int, d, n::Int, lay = 11.0, mat) =
		PartSpec(component, WireArray, layers;
			dim = _spec(d), args = (n, _spec(lay)), material = mat)

	# tube: no extra args
	tube(component::Symbol; layers::Int, t, mat) =
		PartSpec(component, Tubular, layers;
			dim = _spec(t), args = (), material = mat)

	# strip: args are (width, lay)
	strip(component::Symbol; layers::Int, t, w, lay = 0.0, mat) =
		PartSpec(component, Strip, layers;
			dim = _spec(t), args = (_spec(w), _spec(lay)), material = mat)

	# central + hex rings sugar
	function stranded(component::Symbol; layers::Int, d, n::Int, lay = 11.0, mat)
		@assert layers >= 1 "stranded: layers must be ≥ 1 (includes the central wire)."
		specs = PartSpec[]

		# 1) central wire: 1 layer, n=1, lay=0.0
		push!(
			specs,
			PartSpec(component, WireArray, 1;
				dim = _spec(d), args = (1, (0.0, nothing)), material = mat),
		)

		# 2) rings: (layers-1) layers, base n, common lay
		if layers > 1
			push!(
				specs,
				PartSpec(component, WireArray, layers - 1;
					dim = _spec(d), args = (n, _spec(lay)), material = mat),
			)
		end

		return specs
	end

end

module insulator

	using ..ParametricBuilder: PartSpec, _spec
	using ...DataModel: Insulator, Semicon

	tube(component::Symbol; layers::Int, t, mat) =
		PartSpec(component, Insulator, layers;
			dim = _spec(t), args = (), material = mat)

	semi(component::Symbol; layers::Int, t, mat) =
		PartSpec(component, Semicon, layers;
			dim = _spec(t), args = (), material = mat)
end


# # -------------------- tiny value DSL --------------------

# # numeric “specs”: number, (lo,hi,n), or Vector
# _values(x::Number) = (x,)
# _values(v::AbstractVector) = collect(v)
# _values(t::Tuple{<:Number, <:Number, <:Integer}) = range(t[1], t[2]; length = t[3])

# # percent-unc “specs”: nothing (0%), number, (lo,hi,n), or Vector
# _pcts(::Nothing) = (0.0,)
# _pcts(p::Number) = (float(p),)
# _pcts(v::AbstractVector) = map(float, collect(v))
# _pcts(t::Tuple{<:Number, <:Number, <:Integer}) =
# 	range(float(t[1]), float(t[2]); length = t[3])

# # _make_range numeric + pct → Float64 or Measurement
# function _make_range(spec; pct = nothing)
# 	vs  = collect(_values(spec))
# 	pcs = collect(_pcts(pct))
# 	if all(p -> p == 0.0, pcs)
# 		return vs  # deterministic
# 	end
# 	out = Any[]
# 	for v in vs, p in pcs
# 		σ = abs(v) * (p/100)
# 		push!(out, measurement(v, σ))
# 	end
# 	out
# end



# # -------------------- part / component / cable specs --------------------

# """
# PartSpec:
# - component::Symbol           # e.g. :core, :sheath, :jacket
# - part_type::Type             # WireArray, Tubular, Strip, Insulator, Semicon, …
# - n_layers::Int               # how many stacked layers of this part_type
# - dim::(spec, pct)            # diameter OR thickness OR radius (see geom)
# - args::NamedTuple            # other ctor args as specs, e.g. (nstrands_base=6, lay_ratio=(11.0,nothing), width=(10e-3,nothing))
# - material::MaterialSpec
# - geom::Union{Nothing,Symbol} # :diameter | :thickness | :radius | nothing(default)
# """
# struct PartSpec
# 	component::Symbol
# 	part_type::Type
# 	n_layers::Int
# 	dim::Any
# 	args::NamedTuple
# 	material::MaterialSpec
# 	geom::Union{Nothing, Symbol}
# end

# PartSpec(component::Symbol, part_type::Type, n_layers::Int;
# 	dim, args = (;), material::MaterialSpec, geom::Union{Nothing, Symbol} = nothing) =
# 	PartSpec(component, part_type, n_layers, dim, args, material, geom)

# """
# CableBuilderSpec:
# - cable_id::String
# - parts::Vector{PartSpec}  # may interleave conductor/insulator arbitrarily
# - nominal::Union{Nothing,NominalData}
# """
# struct CableBuilderSpec
# 	cable_id::String
# 	parts::Vector{PartSpec}
# 	nominal::Union{Nothing, NominalData}
# end
# CableBuilderSpec(id::AbstractString, parts::Vector{PartSpec}; nominal = nothing) =
# 	CableBuilderSpec(String(id), parts, nominal)

# # -------------------- geometry proxy defaults --------------------

# # decide Diameter vs Thickness vs Radius for a given part in its component sequence
# function _geom_choice(
# 	part_type::Type,
# 	is_first_conductor::Bool,
# 	user::Union{Nothing, Symbol},
# )
# 	if user !== nothing
# 		return user  # honor explicit override (:diameter | :thickness | :radius)
# 	end
# 	if part_type <: AbstractWireArray
# 		return :diameter
# 	elseif part_type <: AbstractConductorPart
# 		return is_first_conductor ? :diameter : :thickness
# 	else
# 		# insulators are always thickness in your framework
# 		return :thickness
# 	end
# end

# # -------------------- build helpers --------------------

# _is_conductor(T::Type) = T <: AbstractConductorPart
# _is_insulator(T::Type) = T <: AbstractInsulatorPart

# # expand a NamedTuple of specs into a NamedTuple iterator (cartesian product)
# function _expand_args(nt::NamedTuple)
# 	keys_ = collect(keys(nt))
# 	vals_ = map(k -> nt[k], keys_)
# 	# each entry can be:
# 	#   - plain Int/Number  -> treat as singleton (no pct)
# 	#   - (spec, pct)       -> _make_range(spec; pct)
# 	function expand_one(v)
# 		if v isa Tuple && length(v) == 2 &&
# 		   (v[1] isa Number || v[1] isa AbstractVector || v[1] isa Tuple)
# 			_make_range(v[1]; pct = v[2])
# 		else
# 			(v,)  # singleton
# 		end
# 	end
# 	spaces = map(expand_one, vals_)
# 	(NamedTuple{Tuple(keys_)}(tpl) for tpl in product(spaces...))
# end

# # -------------------- main: build_designs --------------------
# # Build all variants of a single component, anchored at `base` (0.0 or an InsulatorGroup)
# function _make_variants(ps::Vector{PartSpec}, base)
# 	conductors = [p for p in ps if p.part_type <: AbstractConductorPart]
# 	insulators = [p for p in ps if p.part_type <: AbstractInsulatorPart]
# 	isempty(conductors) && error("component has no conductors")
# 	isempty(insulators) && error("component has no insulators")

# 	# --- conductor group ---
# 	p1c   = conductors[1]
# 	mats1 = _make_range(p1c.material)
# 	dims1 = _make_range(p1c.dim[1]; pct = p1c.dim[2])

# 	comps_with_ig = Tuple{CableComponent, InsulatorGroup}[]

# 	for mat1 in mats1, d1 in dims1
# 		for a1 in _expand_args(p1c.args)
# 			geom = _geom_choice(p1c.part_type, true, p1c.geom)
# 			cg   = _make_first_conductor_group(p1c.part_type, geom, d1, a1, mat1, base)

# 			# stack remaining layers of first conductor spec
# 			for layer ∈ 2:p1c.n_layers
# 				for a in _expand_args(p1c.args)
# 					_add_conductor!(
# 						cg,
# 						p1c.part_type,
# 						_geom_choice(p1c.part_type, false, p1c.geom),
# 						d1,
# 						a,
# 						mat1;
# 						layer,
# 					)
# 				end
# 			end

# 			# other conductor specs
# 			for pc in conductors[2:end]
# 				matsc = _make_range(pc.material)
# 				dimsc = _make_range(pc.dim[1]; pct = pc.dim[2])
# 				for mc in matsc, dc in dimsc
# 					for layer ∈ 1:pc.n_layers
# 						for a in _expand_args(pc.args)
# 							_add_conductor!(
# 								cg,
# 								pc.part_type,
# 								_geom_choice(pc.part_type, false, pc.geom),
# 								dc,
# 								a,
# 								mc;
# 								layer,
# 							)
# 						end
# 					end
# 				end
# 			end

# 			# --- insulator group anchored on cg ---
# 			p1i   = insulators[1]
# 			matsi = _make_range(p1i.material)
# 			dimsi = _make_range(p1i.dim[1]; pct = p1i.dim[2])
# 			for mi in matsi, di in dimsi
# 				ig = _make_first_insulator_group(p1i.part_type, di, p1i.args, mi, cg)
# 				# stack remaining insulators
# 				for pi in insulators[2:end]
# 					m2 = _make_range(pi.material)
# 					d2 = _make_range(pi.dim[1]; pct = pi.dim[2])
# 					for mi2 in m2, di2 in d2
# 						for layer ∈ 1:pi.n_layers
# 							_add_insulator!(ig, pi.part_type, di2, pi.args, mi2)
# 						end
# 					end
# 				end
# 				push!(comps_with_ig, (CableComponent(String(ps[1].component), cg, ig), ig))
# 			end
# 		end
# 	end

# 	return comps_with_ig
# end

# function build_designs(cbs::CableBuilderSpec)
# 	# group parts by component in declared order
# 	comp_names = unique(p.component for p in cbs.parts)
# 	parts_by   = Dict{Symbol, Vector{PartSpec}}()
# 	for p in cbs.parts
# 		get!(parts_by, p.component, PartSpec[]) |> v -> push!(v, p)
# 	end

# 	# each partial: (built_components::Vector{CableComponent}, last_ig::Union{Nothing,InsulatorGroup})
# 	partials = Tuple{Vector{CableComponent}, Union{Nothing, InsulatorGroup}}[(
# 		Vector{CableComponent}(),
# 		nothing,
# 	)]

# 	for cname in comp_names
# 		ps = parts_by[cname]
# 		new_partials = Tuple{Vector{CableComponent}, Union{Nothing, InsulatorGroup}}[]
# 		for (built, last_ig) in partials
# 			base = last_ig === nothing ? 0.0 : last_ig
# 			for (comp, ig) in _make_variants(ps, base)
# 				push!(new_partials, (vcat(built, comp), ig))
# 			end
# 		end
# 		partials = new_partials
# 	end

# 	# finalize designs
# 	designs = CableDesign[]
# 	for (comps, _ig) in partials
# 		@assert !isempty(comps)
# 		des = CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
# 		for k in Iterators.drop(eachindex(comps), 1)
# 			add!(des, comps[k])
# 		end

# 		push!(designs, des)
# 	end
# 	return designs
# end


# # function build_designs(cbs::CableBuilderSpec)::Vector{CableDesign}
# # 	# group by component, preserve original order within each group
# # 	by_comp = Dict{Symbol, Vector{PartSpec}}()
# # 	for p in cbs.parts
# # 		get!(by_comp, p.component, PartSpec[]) |> v -> push!(v, p)
# # 	end
# # 	comp_names = unique(p.component for p in cbs.parts)
# # 	comp_variants = Vector{Vector{CableComponent}}(undef, length(comp_names))

# # 	for (ci, cname) in enumerate(comp_names)
# # 		ps = by_comp[cname]
# # 		# split by role, preserve order within class
# # 		conductors = [p for p in ps if _is_conductor(p.part_type)]
# # 		insulators = [p for p in ps if _is_insulator(p.part_type)]
# # 		isempty(conductors) && error("Component $cname has no conductors.")
# # 		isempty(insulators) && error("Component $cname has no insulators.")

# # 		# --- build conductor group variants ---
# # 		cg_variants = ConductorGroup[]
# # 		begin
# # 			# seed with first conductor
# # 			p1 = conductors[1]
# # 			mats = _make_range(p1.material)
# # 			dims = _make_range(p1.dim[1]; pct = p1.dim[2])
# # 			# additional args (e.g. nstrands_base, lay_ratio, width)
# # 			for mat in mats, d in dims
# # 				# choose geom (Diameter vs Thickness) by contract/defaults
# # 				geom = _geom_choice(p1.part_type, true, p1.geom)
# # 				# layer 1 extras
# # 				for a1 in _expand_args(p1.args)
# # 					cg = _make_first_conductor_group(p1.part_type, geom, d, a1, mat)
# # 					# stack remaining layers of p1 (if n_layers>1)
# # 					for layer ∈ 2:p1.n_layers
# # 						for a in _expand_args(p1.args)
# # 							_add_conductor!(
# # 								cg,
# # 								p1.part_type,
# # 								_geom_choice(p1.part_type, false, p1.geom),
# # 								d,
# # 								a,
# # 								mat;
# # 								layer = layer,
# # 							)
# # 						end
# # 					end
# # 					push!(cg_variants, cg)
# # 				end
# # 			end
# # 		end

# # 		# append remaining conductor specs on top
# # 		if length(conductors) > 1
# # 			for p in conductors[2:end]
# # 				mats = _make_range(p.material)
# # 				dims = _make_range(p.dim[1]; pct = p.dim[2])
# # 				for cg in cg_variants, mat in mats, d in dims
# # 					for layer ∈ 1:p.n_layers
# # 						for a in _expand_args(p.args)
# # 							_add_conductor!(
# # 								cg,
# # 								p.part_type,
# # 								_geom_choice(p.part_type, false, p.geom),
# # 								d,
# # 								a,
# # 								mat;
# # 								layer = layer,
# # 							)
# # 						end
# # 					end
# # 				end
# # 			end
# # 		end

# # 		# --- build insulator group variants anchored to each cg ---
# # 		ig_variants = InsulatorGroup[]
# # 		begin
# # 			p1 = insulators[1]
# # 			mats = _make_range(p1.material)
# # 			dims = _make_range(p1.dim[1]; pct = p1.dim[2])
# # 			for cg in cg_variants, mat in mats, d in dims
# # 				ig = _make_first_insulator_group(p1.part_type, d, p1.args, mat, cg)
# # 				# stack remaining layers of p1 if any
# # 				for layer ∈ 2:p1.n_layers
# # 					_add_insulator!(ig, p1.part_type, d, p1.args, mat)
# # 				end
# # 				push!(ig_variants, ig)
# # 			end
# # 		end

# # 		if length(insulators) > 1
# # 			for p in insulators[2:end]
# # 				mats = _make_range(p.material)
# # 				dims = _make_range(p.dim[1]; pct = p.dim[2])
# # 				for ig in ig_variants, mat in mats, d in dims
# # 					for layer ∈ 1:p.n_layers
# # 						_add_insulator!(ig, p.part_type, d, p.args, mat)
# # 					end
# # 				end
# # 			end
# # 		end

# # 		# zip cg/ig (parallel variants)
# # 		comps = CableComponent[]
# # 		for (cg, ig) in zip(cg_variants, ig_variants)
# # 			push!(comps, CableComponent(String(cname), cg, ig))
# # 		end
# # 		comp_variants[ci] = comps
# # 	end

# # 	# product over components → CableDesign
# # 	out = CableDesign[]
# # 	for combo in product(comp_variants...)
# # 		comps = collect(combo)
# # 		design = CableDesign(cbs.cable_id, comps[1]; nominal_data = cbs.nominal)
# # 		for k in Iterators.drop(eachindex(comps), 1)
# # 			add!(design, comps[k])
# # 		end
# # 		push!(out, design)
# # 	end
# # 	out
# # end

# """
# 	iterate_designs(cbs) -> generator of CableDesign

# Similar behavior to `build_designs` but yields one `CableDesign` at a time.
# """
# function iterate_designs(cbs::CableBuilderSpec)
# 	comp_names = unique(p.component for p in cbs.parts)
# 	comp_sets  = [_component_variants(cbs, name) for name in comp_names]
# 	prod       = Iterators.product(comp_sets...)             # lazy Cartesian product
# 	return (_design_from_combo(cbs.cable_id, cbs.nominal, combo) for combo in prod)
# end


# # first conductor group
# @inline function _make_first_conductor_group(
# 	T::Type,
# 	geom::Symbol,
# 	dim,
# 	args::NamedTuple,
# 	mat,
# 	base,
# )
# 	if T <: AbstractWireArray
# 		nbase = args[:nstrands_base]
# 		lay   = get(args, :lay_ratio, 0.0)
# 		return ConductorGroup(WireArray(base, _to_geom(geom, dim), nbase, lay, mat))
# 	elseif T === Tubular
# 		return ConductorGroup(Tubular(base, _to_geom(geom, dim), mat))
# 	elseif T === Strip
# 		width = args[:width]
# 		lay   = get(args, :lay_ratio, 0.0)
# 		return ConductorGroup(Strip(base, _to_geom(:thickness, dim), width, lay, mat))
# 	else
# 		return ConductorGroup(T(base, _to_geom(geom, dim), mat))
# 	end
# end

# # stacking conductors
# function _add_conductor!(
# 	cg::ConductorGroup,
# 	T::Type,
# 	geom::Symbol,
# 	dim,
# 	args::NamedTuple,
# 	mat;
# 	layer::Int,
# )
# 	if T <: AbstractWireArray
# 		nbase = args[:nstrands_base]
# 		lay   = get(args, :lay_ratio, 0.0)
# 		add!(cg, WireArray, _to_geom(geom, dim), layer*nbase, lay, mat)
# 	elseif T === Tubular
# 		add!(cg, Tubular, _to_geom(geom, dim), mat)
# 	elseif T === Strip
# 		width = args[:width]
# 		lay   = get(args, :lay_ratio, 0.0)
# 		add!(cg, Strip, _to_geom(:thickness, dim), width, lay, mat)
# 	else
# 		add!(cg, T, _to_geom(geom, dim), mat)
# 	end
# end

# # first insulator group around cg
# function _make_first_insulator_group(
# 	T::Type,
# 	dim,
# 	args::NamedTuple,
# 	mat,
# 	cg::ConductorGroup,
# )
# 	if T === Insulator
# 		return InsulatorGroup(Insulator(cg, Thickness(dim), mat))
# 	elseif T === Semicon
# 		return InsulatorGroup(Semicon(cg, Thickness(dim), mat))
# 	else
# 		return InsulatorGroup(T(cg, Thickness(dim), mat))
# 	end
# end

# # stacking insulators
# function _add_insulator!(ig::InsulatorGroup, T::Type, dim, args::NamedTuple, mat)
# 	if T === Insulator
# 		add!(ig, Insulator, Thickness(dim), mat)
# 	elseif T === Semicon
# 		add!(ig, Semicon, Thickness(dim), mat)
# 	else
# 		add!(ig, T, Thickness(dim), mat)
# 	end
# end

# # proxy conversion
# @inline function _to_geom(s::Symbol, dim)
# 	s === :diameter && return Diameter(dim)
# 	s === :thickness && return Thickness(dim)
# 	s === :radius && return dim  # already a radius
# 	error("unsupported geometry symbol '$s'")
# end

# """
# 	count_combinations(cbs) -> Int
# Rough upper bound of the cartesian product size.
# """
# function count_combinations(cbs::CableBuilderSpec)
# 	total = 1
# 	for p in cbs.parts
# 		# dimension values × pct choices
# 		dim_vals =
# 			length(collect((_ -> _make_range(p.dim[1]; pct = p.dim[2])))(nothing))
# 		mats = _make_range(p.material)
# 		# args space
# 		nargs = 1
# 		for (_k, v) in pairs(p.args)
# 			if v isa Tuple && length(v)==2
# 				nargs *= length(_make_range(v[1]; pct = v[2]))
# 			else
# 				nargs *= 1
# 			end
# 		end
# 		# for each layer
# 		total *= dim_vals * length(mats) * nargs ^ p.n_layers
# 	end
# 	total
# end

# """
# 	map_baseparams(designs) -> Vector{NamedTuple}

# Accepts either `Vector{CableDesign}` or the lazy iterator.
# Returns a vector of (id, R, L, C, design) for quick filtering/comparison.
# R in Ω/km, L in mH/km, C in μF/km as per your DataFrame(..., :baseparams) order.
# """
# function map_baseparams(designs)
# 	out = NamedTuple[]
# 	for des in designs
# 		df = DataFrame(des, :baseparams)
# 		R = df[1, :computed]      # Ω/km
# 		L = df[2, :computed]      # mH/km
# 		C = df[3, :computed]      # μF/km
# 		push!(out, (id = des.cable_id, R = R, L = L, C = C, design = des))
# 	end
# 	out
# end

# # ===== Namespaced DSL =========================================================
# # shared tiny guards: accept number OR (spec, pct)
# _dim(x) = (x isa Tuple && length(x) == 2) ? x : (x, nothing)
# _arg(x) = (x isa Tuple && length(x) == 2) ? x : (x, nothing)

# # -- conductor namespace
# module conductor
# 	using ..ParametricBuilder: PartSpec, _dim, _arg
# 	using ...DataModel: WireArray, Tubular, Strip

# 	# WireArray (first conductor defaults to Diameter; stacked becomes Thickness automatically by builder)
# 	wire(
# 		component::Symbol;
# 		layers::Int,
# 		d,
# 		n::Int,
# 		lay = 11.0,
# 		mat,
# 		geom::Union{Nothing, Symbol} = nothing,
# 	) =
# 		PartSpec(component, WireArray, layers;
# 			dim = _dim(d),
# 			args = (nstrands_base = n, lay_ratio = _arg(lay)),
# 			material = mat, geom = geom)

# 	# Tubular conductor
# 	tube(component::Symbol; layers::Int, t, mat, geom::Union{Nothing, Symbol} = nothing) =
# 		PartSpec(component, Tubular, layers;
# 			dim = _dim(t),
# 			material = mat, geom = geom)

# 	# Strip conductor (namespaced; no clash with Base.strip)
# 	strip(component::Symbol; layers::Int, t, w, lay = 0.0, mat) =
# 		PartSpec(component, Strip, layers;
# 			dim = _dim(t),
# 			args = (width = _arg(w), lay_ratio = _arg(lay)),
# 			material = mat)

# end # submodule conductor

# # -- insulator namespace
# module insulator
# 	using ..ParametricBuilder: PartSpec, _dim
# 	using ...DataModel: Insulator, Semicon

# 	# “tubular” insulator (current coaxial form)
# 	tube(component::Symbol; layers::Int, t, mat) =
# 		PartSpec(component, Insulator, layers; dim = _dim(t), material = mat)

# 	# semiconductor (also tubular in current implementation)
# 	semi(component::Symbol; layers::Int, t, mat) =
# 		PartSpec(component, Semicon, layers; dim = _dim(t), material = mat)

# end # submodule insulator
# # ============================================================================== 

# Submodule `WirePatterns`
include("wirepatterns/WirePatterns.jl")
using .WirePatterns

end # module ParametricBuilder
