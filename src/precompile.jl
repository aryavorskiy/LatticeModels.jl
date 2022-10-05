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
    Base.precompile(Tuple{typeof(exp),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.7767232
    Base.precompile(Tuple{typeof(^),LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}},Int64})   # time: 0.5976823
    Base.precompile(Tuple{typeof(^),BondSet{SquareLattice},Int64})   # time: 0.4571465
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},SquareLattice,Function,Hopping{Matrix{Float64}},Landau})   # time: 0.4400253
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},typeof(_always_true_on_lattice),Hopping{Matrix{Float64}},Landau})   # time: 0.4128046
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SubLattice{SquareLattice},Matrix{ComplexF64}}})   # time: 0.4127565
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},SquareLattice,Function,Hopping{Matrix{ComplexF64}},Landau})   # time: 0.3488409
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},typeof(_always_true_on_lattice),Hopping{Matrix{ComplexF64}},Landau})   # time: 0.3379313
    Base.precompile(Tuple{typeof(hopping_operator),SubLattice{SquareLattice},Hopping{Matrix{Int64}}})   # time: 0.293779
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}}})   # time: 0.2479526
    Base.precompile(Tuple{typeof(spectrum),LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}})   # time: 0.1582789
    Base.precompile(Tuple{typeof(coords),SquareLattice,LatticeIndex})   # time: 0.1418465
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice}})   # time: 0.1270032
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.1144147
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat,LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}}})   # time: 0.1097949
    Base.precompile(Tuple{typeof(bonds),SquareLattice,Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.1092177
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.1076165
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Vararg{Int64}})   # time: 0.0798423
    Base.precompile(Tuple{typeof(+),LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}},LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}},LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}})   # time: 0.0791366
    Base.precompile(Tuple{typeof(diag_operator),Function,Basis{SubLattice{SquareLattice}}})   # time: 0.0745952
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}},Int64})   # time: 0.0732307
    Base.precompile(Tuple{typeof(+),LatticeVecOrMat{SquareLattice,Matrix{Float64}},LatticeVecOrMat{SquareLattice,Matrix{Float64}}})   # time: 0.0633114
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}},Any})   # time: 0.0605246
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0526841
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:axis,),Tuple{Int64}},Type{Hopping}})   # time: 0.0500769
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}}})   # time: 0.045432
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0432408
    Base.precompile(Tuple{typeof(*),Int64,LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}}})   # time: 0.0389378
    Base.precompile(Tuple{typeof(_diag_from_macro),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},SquareLattice,Function})   # time: 0.0329321
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0321736
    Base.precompile(Tuple{typeof(_diag_from_macro),SquareLattice,Function})   # time: 0.0297966
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:axis,),Tuple{Int64}},Type{Hopping},Matrix{ComplexF64}})   # time: 0.0219545
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVStyle})   # time: 0.0195809
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64},Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(*),Tuple{LatticeValue{Float64},Int64}}})   # time: 0.0176338
    Base.precompile(Tuple{typeof(diag_operator),Function,Basis{SquareLattice}})   # time: 0.0175454
    Base.precompile(Tuple{typeof(evolution_operator),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},Float64})   # time: 0.0136178
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}})   # time: 0.0114729
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0106231
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64},Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(|>),Tuple{LatticeValue{ComplexF64},Base.RefValue{typeof(real)}}}})   # time: 0.0105725
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat,LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0099635
    Base.precompile(Tuple{typeof(==),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0099262
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64},LatticeValue{Float64}})   # time: 0.0098581
    Base.precompile(Tuple{typeof(==),LatticeVecOrMat{SquareLattice,Matrix{Float64}},LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0097662
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},LatticeVecOrMat{SquareLattice,Adjoint{ComplexF64,Matrix{ComplexF64}}}})   # time: 0.00972
    Base.precompile(Tuple{typeof(-),LatticeVecOrMat{SquareLattice,Matrix{Float64}},UniformScaling{Bool}})   # time: 0.0091452
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64},Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(*),Tuple{LatticeValue{Float64},LatticeValue{Float64}}}})   # time: 0.0080298
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64},Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(*),Tuple{Int64,LatticeValue{Float64}}}})   # time: 0.0077927
    Base.precompile(Tuple{typeof(+),LatticeVecOrMat{SquareLattice,Matrix{Float64}},LatticeVecOrMat{SquareLattice,Matrix{Float64}},LatticeVecOrMat{SquareLattice,Matrix{Float64}}})   # time: 0.0065302
    Base.precompile(Tuple{typeof(spectrum),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0064976
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0059513
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}},LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}}})   # time: 0.0058373
    Base.precompile(Tuple{typeof(_diag_from_macro),SquareLattice,Matrix{Int64}})   # time: 0.0057676
    Base.precompile(Tuple{typeof(_wrap_smart!),Expr})   # time: 0.0056397
    Base.precompile(Tuple{typeof(*),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0056357
    Base.precompile(Tuple{typeof(is_adjacent),BondSet{SquareLattice},LatticeIndex,LatticeIndex})   # time: 0.0046569
    Base.precompile(Tuple{typeof(setindex!),LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}},Matrix{ComplexF64},Int64,Int64})   # time: 0.0044994
    Base.precompile(Tuple{Core.kwftype(typeof(Type)),NamedTuple{(:axis,),Tuple{Int64}},Type{Hopping},Matrix{Float64}})   # time: 0.0038045
    isdefined(LatticeModels, Symbol("#15#16")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#15#16")),SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}})   # time: 0.0036949
    Base.precompile(Tuple{typeof(adjoint),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0035785
    Base.precompile(Tuple{Type{LatticeValue},Function,SquareLattice})   # time: 0.0035679
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice},Tuple{Matrix{Float64}},Tuple{LatticeVecOrMat{SquareLattice,Matrix{Float64}},LatticeVecOrMat{SquareLattice,Matrix{Float64}}}})   # time: 0.003272
    Base.precompile(Tuple{typeof(/),LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}},Int64})   # time: 0.003174
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeVecOrMat{SquareLattice,SubArray{Float64,2,Array{Float64,3},Tuple{Base.Slice{Base.OneTo{Int64}},Base.Slice{Base.OneTo{Int64}},Int64},true}},Vararg{Any}})   # time: 0.0025584
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeVecOrMat{SquareLattice,Matrix{Float64}},Vararg{LatticeVecOrMat{SquareLattice,Matrix{Float64}}}})   # time: 0.0022784
    let fbody = try
            __lookup_kwbody__(which(_unwrap_wlattice, (Function, Basis{SubLattice{SquareLattice}}, Tuple{Matrix{ComplexF64}}, LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}, Tuple{LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}},)))
        catch missing
        end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol,Union{},Tuple{},NamedTuple{(),Tuple{}}}, typeof(_unwrap_wlattice), Function, Basis{SubLattice{SquareLattice}}, Tuple{Matrix{ComplexF64}}, LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}, Tuple{LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}},))
        end
    end   # time: 0.0019719
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(*),Tuple{LatticeValue{Float64},LatticeValue{Float64}}},Type{Float64}})   # time: 0.001955
    Base.precompile(Tuple{typeof(_match),Hopping{Matrix{Int64}},SquareLattice,LatticeIndex,LatticeIndex})   # time: 0.0019057
    let fbody = try
            __lookup_kwbody__(which(_unwrap_wlattice, (Function, Basis{SquareLattice}, Tuple{Matrix{ComplexF64}}, LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}, Tuple{LatticeVecOrMat{SquareLattice,Adjoint{ComplexF64,Matrix{ComplexF64}}}},)))
        catch missing
        end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol,Union{},Tuple{},NamedTuple{(),Tuple{}}}, typeof(_unwrap_wlattice), Function, Basis{SquareLattice}, Tuple{Matrix{ComplexF64}}, LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}, Tuple{LatticeVecOrMat{SquareLattice,Adjoint{ComplexF64,Matrix{ComplexF64}}}},))
        end
    end   # time: 0.001757
    Base.precompile(Tuple{typeof(getindex),LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}},Int64,Int64})   # time: 0.001692
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeVecOrMat{SquareLattice,Matrix{ComplexF64}}})   # time: 0.0015691
    Base.precompile(Tuple{typeof(getindex),LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}},Int64,Int64})   # time: 0.0014973
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(|>),Tuple{LatticeValue{ComplexF64},Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0013647
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,Int64,Vararg{Int64}})   # time: 0.0011742
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(*),Tuple{LatticeValue{Float64},Int64}},Type{Float64}})   # time: 0.0010771
    let fbody = try
            __lookup_kwbody__(which(_unwrap_wlattice, (Function, Basis{SquareLattice}, Tuple{Matrix{Float64}}, LatticeVecOrMat{SquareLattice,Matrix{Float64}}, Tuple{LatticeVecOrMat{SquareLattice,Matrix{Float64}}},)))
        catch missing
        end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol,Union{},Tuple{},NamedTuple{(),Tuple{}}}, typeof(_unwrap_wlattice), Function, Basis{SquareLattice}, Tuple{Matrix{Float64}}, LatticeVecOrMat{SquareLattice,Matrix{Float64}}, Tuple{LatticeVecOrMat{SquareLattice,Matrix{Float64}}},))
        end
    end   # time: 0.001051
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle,Tuple{Base.OneTo{Int64}},typeof(*),Tuple{Int64,LatticeValue{Float64}}},Type{Float64}})   # time: 0.0010443
    let fbody = try
            __lookup_kwbody__(which(_unwrap_wlattice, (Function, Basis{SubLattice{SquareLattice}}, Tuple{Matrix{ComplexF64}}, Tuple{LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}},LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}},)))
        catch missing
        end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol,Union{},Tuple{},NamedTuple{(),Tuple{}}}, typeof(_unwrap_wlattice), Function, Basis{SubLattice{SquareLattice}}, Tuple{Matrix{ComplexF64}}, Tuple{LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}},LatticeVecOrMat{SubLattice{SquareLattice},Matrix{ComplexF64}}},))
        end
    end   # time: 0.0010064
end
