const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(arg)
        isa(arg, Symbol) && return arg
        @assert isa(arg, GlobalRef)
        return arg.name
    end

    f = get(__bodyfunction__, mnokw, nothing)
    if f === nothing
        fmod = mnokw.module
        # The lowered code for `mnokw` should look like
        #   %1 = mkw(kwvalues..., #self#, args...)
        #        return %1
        # where `mkw` is the name of the "active" keyword body-function.
        ast = Base.uncompressed_ast(mnokw)
        if isa(ast, Core.CodeInfo) && length(ast.code) >= 2
            callexpr = ast.code[end-1]
            if isa(callexpr, Expr) && callexpr.head == :call
                fsym = callexpr.args[1]
                if isa(fsym, Symbol)
                    f = getfield(fmod, fsym)
                elseif isa(fsym, GlobalRef)
                    if fsym.mod === Core && fsym.name === :_apply
                        f = getfield(mnokw.module, getsym(callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(callexpr.args[3]))
                    else
                        f = getfield(fsym.mod, fsym.name)
                    end
                else
                    f = missing
                end
            else
                f = missing
            end
        else
            f = missing
        end
        __bodyfunction__[mnokw] = f
    end
    return f
end

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.4968179
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{_A, Matrix{ComplexF64}} where _A<:Lattice,SquareLattice{2},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.4957111
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64,Float64})   # time: 0.459173
    Base.precompile(Tuple{typeof(site_coords),HoneycombLattice,LatticeSite{2}})   # time: 0.4100649
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{_A, Matrix{ComplexF64}} where _A<:Lattice,SquareLattice{2},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.3777554
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Float64})   # time: 0.315348
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.2522875
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.2191144
    Base.precompile(Tuple{typeof(!),BondSet{SquareLattice{2}}})   # time: 0.1928145
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.1538048
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.137659
    Base.precompile(Tuple{typeof(^),BondSet{SquareLattice{2}},Int64})   # time: 0.1367733
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.1069798
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0963936
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},Function,Hopping{Matrix{Int64}},LandauField})   # time: 0.0921112
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0844346
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0843398
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0818063
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0704339
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0645587
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Adjoint{Float64, Vector{Float64}}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeArray{SquareLattice{2}, Vector{Float64}}})   # time: 0.0590546
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0492399
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64,Int64})   # time: 0.0471388
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0307529
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0298411
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}}})   # time: 0.0214686
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0193318
    Base.precompile(Tuple{Type{HoneycombLattice},Function,Int64,Int64})   # time: 0.0179485
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0170474
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.015958
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.01417
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0141179
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{SquareLattice{2}, Vector{Float64}}})   # time: 0.0132733
    Base.precompile(Tuple{typeof(|),BondSet{SquareLattice{2}},BondSet{SquareLattice{2}}})   # time: 0.0103077
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.010027
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0098479
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0085282
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.0076961
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0075966
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0066418
    Base.precompile(Tuple{typeof(^),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.00506
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}},Type{Float64}})   # time: 0.0042761
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0040795
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}, LatticeArray{SquareLattice{2}, Vector{Float64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}, LatticeArray{SquareLattice{2}, Vector{Float64}}},))
        end
    end   # time: 0.0035878
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0033719
    isdefined(LatticeModels, Symbol("#18#19")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#18#19")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0030225
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{LatticeArray{SquareLattice{2}, Vector{Float64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{LatticeArray{SquareLattice{2}, Vector{Float64}}},))
        end
    end   # time: 0.0029516
    Base.precompile(Tuple{typeof(_lazy_tp),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0029336
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0028374
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0026476
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Tuple{}})   # time: 0.0024595
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.002374
    Base.precompile(Tuple{typeof(randn),SquareLattice{2}})   # time: 0.002291
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Int64,Tuple{Int64},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Int64,Tuple{Int64},))
        end
    end   # time: 0.0020438
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0019431
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0019339
    Base.precompile(Tuple{typeof(coord_operators),SquareLattice{2},Int64})   # time: 0.0018323
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Float64}})   # time: 0.0016946
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0016573
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0014062
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{Int64}})   # time: 0.0011363
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Matrix{ComplexF64},Int64,Int64})   # time: 0.0010322
end
