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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:reduce_fn, :sort), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0620564
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{1}})   # time: 0.8366427
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Int64})   # time: 0.5234054
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SingleParticleBasis{SquareLattice{2}}, Matrix{ComplexF64}}})   # time: 0.4791743
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{1}}},Nothing,Hopping{1},NoField})   # time: 0.4354512
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.3918187
    Base.precompile(Tuple{typeof(Haldane),HoneycombLattice,Int64,Int64,Int64})   # time: 0.3603056
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:j1, :j2, :index), Tuple{Int64, Int64, Int64}},typeof(getindex),HoneycombLattice})   # time: 0.3309681
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.2984177
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.2662462
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),Function,LatticeValue{Float64, :square}})   # time: 0.2587308
    Base.precompile(Tuple{typeof(pade_exp),Matrix{ComplexF64},Int64})   # time: 0.2294125
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeSite{2}})   # time: 0.2277251
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.2246569
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.214818
    Base.precompile(Tuple{typeof(ldos),Spectrum{SingleParticleBasis{SquareLattice{2}}, Matrix{ComplexF64}},Int64,Float64})   # time: 0.2000251
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.1965443
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1859396
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Symbol})   # time: 0.1852091
    Base.precompile(Tuple{typeof(coords),SingleParticleBasis{HoneycombLattice}})   # time: 0.1779326
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:pbc,), Tuple{Bool}},typeof(TightBinding),SquareLattice{1}})   # time: 0.1748882
    Base.precompile(Tuple{typeof(ldos),Spectrum{SingleParticleBasis{SquareLattice{2}}, Matrix{ComplexF64}},Float64})   # time: 0.1695119
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.1682023
    Base.precompile(Tuple{typeof(diff),LatticeValueRecord})   # time: 0.1551369
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1471355
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.1300735
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.1278806
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.1268993
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},SMatrix{2, 2, Float64, 4},Int64,Int64})   # time: 0.122939
    Base.precompile(Tuple{typeof(coords),SingleParticleBasis{SquareLattice{2}}})   # time: 0.1205561
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LandauField})   # time: 0.1153717
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.1102889
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1077307
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.1035592
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.1010283
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{2},LandauField})   # time: 0.0997515
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.0966628
    Base.precompile(Tuple{Type{SquareLattice},Int64})   # time: 0.0962338
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.0942269
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Bool}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0934221
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0897517
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.0897024
    Base.precompile(Tuple{typeof(integrate),LatticeValueRecord})   # time: 0.0894642
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Nothing,Hopping{2},LandauField})   # time: 0.0871479
    Base.precompile(Tuple{typeof(line_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0854924
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.0818278
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Any})   # time: 0.0809664
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.0806011
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Float64})   # time: 0.0797621
    isdefined(LatticeModels, Symbol("#88#89")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#88#89")),Int64})   # time: 0.0797146
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SingleParticleBasis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0791378
    Base.precompile(Tuple{Type{LatticeTimeSequence},Vector{LatticeValue{Float64, :square}},Vector{Float64}})   # time: 0.076106
    Base.precompile(Tuple{typeof(line_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0742832
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeValue{Float64, :square}})   # time: 0.0736744
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},SMatrix{2, 2, Int64, 4},Int64,Int64})   # time: 0.0727854
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:x, :x2), Tuple{Int64, Int64}},typeof(getindex),LatticeValue{Float64, :square}})   # time: 0.0694017
    isdefined(LatticeModels, Symbol("#84#85")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#84#85")),Int64})   # time: 0.0626882
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0626159
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0623239
    Base.precompile(Tuple{typeof(==),LatticeValueRecord,LatticeValueRecord})   # time: 0.062036
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),Function,SquareLattice{2}})   # time: 0.061175
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SingleParticleBasis{SquareLattice{2}}, Matrix{ComplexF64}},Float64})   # time: 0.0609035
    Base.precompile(Tuple{typeof(line_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0595919
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Nothing,Hopping,LandauField})   # time: 0.0594814
    let fbody = try __lookup_kwbody__(which(map_currents, (typeof(site_distance),DensityCurrents,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (typeof(sum),Bool,typeof(map_currents),typeof(site_distance),DensityCurrents,))
    end
end   # time: 0.0586801
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},Nothing,Hopping{1},NoField})   # time: 0.0570815
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{SingleParticleBasis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0554889
    Base.precompile(Tuple{typeof(line_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.055067
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0534229
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeValue{Bool, :square}})   # time: 0.052836
    Base.precompile(Tuple{LatticeValueRecord,Float64,Float64})   # time: 0.0525196
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0515842
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(TightBinding),LatticeValue{Float64, :honeycomb}})   # time: 0.0500472
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0492647
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Float64,Nothing,Bool})   # time: 0.0469417
    Base.precompile(Tuple{typeof(site_density),LatticeArray{Vector{ComplexF64}, SingleParticleBasis{HoneycombLattice}, 1}})   # time: 0.0468291
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{1},Vararg{Hopping{1}}})   # time: 0.0430492
    Base.precompile(Tuple{typeof(taylor_exp),Matrix{ComplexF64},Int64})   # time: 0.0425134
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Float64,Int64,Bool})   # time: 0.0418916
    Base.precompile(Tuple{typeof(+),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0416425
    Base.precompile(Tuple{typeof(*),LatticeOperator{Adjoint{Float64, Vector{Float64}}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}})   # time: 0.0415583
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0409443
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0406948
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Int64})   # time: 0.0377933
    Base.precompile(Tuple{typeof(materialize),SparseMatrixBuilder{ComplexF64}})   # time: 0.037503
    Base.precompile(Tuple{typeof(dot),LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1},LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}})   # time: 0.0365709
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.035384
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.035135
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0340383
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.0335298
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0318886
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0318433
    Base.precompile(Tuple{typeof(site_index),Lattice,LatticeSite{2}})   # time: 0.0305633
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.0305158
    Base.precompile(Tuple{typeof(==),LatticeOperator{SparseMatrixCSC{ComplexF64, Int64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0299555
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2}})   # time: 0.0290581
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0287301
    Base.precompile(Tuple{typeof(diag_operator),LatticeValue{Float64, :square},Int64})   # time: 0.0281376
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{SparseMatrixBuilder{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Function,Hopping{2},LandauField})   # time: 0.0279213
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{HoneycombLattice}},Matrix{Float64},Int64,Int64})   # time: 0.0275876
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0272062
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVWStyle})   # time: 0.0263291
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0258713
    Base.precompile(Tuple{typeof(==),LatticeSite,LatticeSite})   # time: 0.0256958
    Base.precompile(Tuple{typeof(zero_on_basis),SingleParticleBasis{SquareLattice{2}},Type{SparseMatrixBuilder{ComplexF64}}})   # time: 0.0230609
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(exp), Tuple{LatticeValue{Float64, :honeycomb}}}}}})   # time: 0.0229143
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0225641
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},Nothing,Hopping{1},LandauField})   # time: 0.0217988
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{1},LandauField})   # time: 0.0213338
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0206368
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0201774
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Function,Hopping{2},LandauField})   # time: 0.0201281
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:pade, :k), Tuple{Bool, Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0195245
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0190175
    let fbody = try __lookup_kwbody__(which(SpinTightBinding, (LatticeValue{Float64, :square},Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(SpinTightBinding),LatticeValue{Float64, :square},Vararg{Any},))
    end
end   # time: 0.0187215
    isdefined(LatticeModels, Symbol("#102#103")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#102#103")),Float64})   # time: 0.0178242
    Base.precompile(Tuple{LatticeValueRecord,Float64})   # time: 0.0173827
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0172975
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.017262
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0169836
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0156987
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.014591
    Base.precompile(Tuple{typeof(getindex),Spectrum{SingleParticleBasis{HoneycombLattice}, Matrix{ComplexF64}},Int64})   # time: 0.0145357
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0139526
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0138717
    Base.precompile(Tuple{typeof(line_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64}})   # time: 0.013861
    Base.precompile(Tuple{typeof(dot),LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}})   # time: 0.0135801
    Base.precompile(Tuple{typeof(_kws_to_mask),SquareLattice{2},Any})   # time: 0.0128924
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0127437
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0126052
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}})   # time: 0.0125302
    Base.precompile(Tuple{typeof(-),MaterializedCurrents,MaterializedCurrents})   # time: 0.0124721
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0124325
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0122385
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{1}})   # time: 0.0120914
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0120006
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Vector{Bool}}},typeof(hopping)})   # time: 0.0119517
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0117366
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0116902
    Base.precompile(Tuple{typeof(==),LatticeArray{Vector{ComplexF64}, SingleParticleBasis{HoneycombLattice}, 1},LatticeArray{Vector{ComplexF64}, SingleParticleBasis{HoneycombLattice}, 1}})   # time: 0.011684
    Base.precompile(Tuple{typeof(zero_on_basis),SingleParticleBasis{SquareLattice{1}},Type{Matrix{ComplexF64}}})   # time: 0.0110465
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}}})   # time: 0.0110303
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.0107115
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.0105927
    Base.precompile(Tuple{Type{LatticeOperator},SingleParticleBasis{SquareLattice{2}},UniformScaling{Bool}})   # time: 0.0105129
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0102778
    isdefined(LatticeModels, Symbol("#108#110")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#108#110")),Vector{Float64}})   # time: 0.0100866
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0100554
    Base.precompile(Tuple{typeof(line_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.009822
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0098096
    Base.precompile(Tuple{typeof(getindex),Spectrum{SingleParticleBasis{HoneycombLattice}, Matrix{ComplexF64}},BitVector})   # time: 0.0097496
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}}})   # time: 0.0096538
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0094955
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0093453
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.0089195
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0088847
    Base.precompile(Tuple{typeof(==),Hopping{1},Hopping{1}})   # time: 0.0088474
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0084789
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0084075
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0083269
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.007907
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0077823
    Base.precompile(Tuple{Core.kwftype(typeof(Base.Broadcast.dotview)),NamedTuple{(:x1,), Tuple{Int64}},typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square}})   # time: 0.0076002
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.0073738
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0073004
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SparseMatrixBuilder{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Matrix{Complex{Int64}}})   # time: 0.0072157
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0072117
    Base.precompile(Tuple{typeof(-),LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},UniformScaling{Bool}})   # time: 0.0068011
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0065288
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Matrix{ComplexF64},Int64,Int64})   # time: 0.0062569
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Nothing,Hopping{2},NoField})   # time: 0.006218
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0062055
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0059602
    Base.precompile(Tuple{typeof(==),LatticeValue,LatticeValue{Float64, :square}})   # time: 0.0058057
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{SparseMatrixBuilder{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.005801
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0057283
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0055673
    Base.precompile(Tuple{typeof(insert!),LatticeValueRecord,Int64,LatticeValue{Float64, :square}})   # time: 0.0054498
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0054389
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.0051966
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Float64, :square}})   # time: 0.0051898
    Base.precompile(Tuple{typeof(_expr_depends_on),Expr,Symbol})   # time: 0.0050973
    Base.precompile(Tuple{typeof(materialize),TensorProduct{LatticeValue{Float64, :square}, 2, Bool}})   # time: 0.0050313
    Base.precompile(Tuple{typeof(materialize),TensorProduct{LatticeValue{Float64, :square}, 2, Int64}})   # time: 0.00487
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{1}}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.0048557
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{1}})   # time: 0.0048456
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},TensorProduct{LatticeValue{Float64, :honeycomb}, 2, Int64}})   # time: 0.004824
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0046812
    Base.precompile(Tuple{typeof(materialize),TensorProduct{LatticeValue{Float64, :honeycomb}, 2, Int64}})   # time: 0.0046479
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0046357
    Base.precompile(Tuple{typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0045436
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0045226
    Base.precompile(Tuple{typeof(zero_on_basis),SingleParticleBasis{SquareLattice{2}},Type{Matrix{ComplexF64}}})   # time: 0.004434
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0043585
    Base.precompile(Tuple{typeof(getindex),Spectrum{SingleParticleBasis{HoneycombLattice}, Matrix{ComplexF64}},UnitRange{Int64}})   # time: 0.0043388
    Base.precompile(Tuple{typeof(_make_wrapper),SparseMatrixCSC{ComplexF64, Int64},SingleParticleBasis{SquareLattice{2}}})   # time: 0.0042716
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0042389
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64},Bravais{1, 1}})   # time: 0.0042238
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},TensorProduct{_A, _B, Int64} where {_A<:(LatticeValue{<:Number}), _B}})   # time: 0.0041351
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0041155
    Base.precompile(Tuple{Type{LatticeArray},SingleParticleBasis{SquareLattice{2}},Vector{Float64}})   # time: 0.0040917
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Matrix{Int64}})   # time: 0.0040324
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0039264
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0037426
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0037143
    isdefined(LatticeModels, Symbol("#55#56")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#55#56")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0036115
    isdefined(LatticeModels, Symbol("#108#110")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#108#110")),Array})   # time: 0.0035658
    isdefined(LatticeModels, Symbol("#55#56")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#55#56")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0034975
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0034763
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.0034742
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.0034069
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0033883
    isdefined(LatticeModels, Symbol("#109#111")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#109#111")),Tuple{Float64, Vector{Float64}}})   # time: 0.0033857
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0033201
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},SMatrix{2, 2, Float64, 4},Int64,Int64})   # time: 0.003317
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0033092
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.0032455
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0032323
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(exp), Tuple{LatticeValue{Float64, :honeycomb}}}}},Type{Float64}})   # time: 0.0031473
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}})   # time: 0.0031205
    Base.precompile(Tuple{typeof(line_integral),LandauField,SVector{2, Float64},SVector{2, Float64}})   # time: 0.0031049
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0030444
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},Matrix{Float64},Int64,Int64})   # time: 0.0030374
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0030219
    Base.precompile(Tuple{typeof(SpinTightBinding),LatticeValue{Float64, :square},Int64})   # time: 0.0029957
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{1}})   # time: 0.0029925
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :axis), Tuple{Tuple{Int64, Int64}, Int64}},typeof(hopping),Int64})   # time: 0.0029879
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0029697
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0029074
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0028471
    Base.precompile(Tuple{typeof(_wrap_eye),String,Matrix{Bool}})   # time: 0.0028426
    Base.precompile(Tuple{typeof(zero_on_basis),SingleParticleBasis{HoneycombLattice}})   # time: 0.0028314
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.002818
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},TensorProduct{LatticeValue{Float64, :square}, 2, Bool}})   # time: 0.0027943
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0027454
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}, LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.0027308
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},TensorProduct{LatticeValue{Float64, :square}, 2, Int64}})   # time: 0.0026825
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.002678
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0026484
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0026261
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0026229
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0025985
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.0024902
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},TensorProduct{LatticeValue{Int64, :square}, 2, Int64}})   # time: 0.0024515
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.0024447
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0024316
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Vararg{LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.0023751
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0023539
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},Vararg{LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.0023009
    Base.precompile(Tuple{typeof(init_record),LatticeValue{Float64, :square}})   # time: 0.0022687
    Base.precompile(Tuple{Type{Basis},SquareLattice{2},Int64})   # time: 0.0022519
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,SingleParticleBasis{SquareLattice{2}},Tuple{Vector{Float64}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Vector{Float64}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0022402
    Base.precompile(Tuple{Type{LatticeTimeSequence},Vector{LatticeValue{Float64, :square}},Vector{Int64}})   # time: 0.0022257
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Vararg{Any}})   # time: 0.0022228
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{LatticeOperator{SparseMatrixBuilder{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.0021868
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0021309
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,SingleParticleBasis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0021157
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}},))
    end
end   # time: 0.002083
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.002045
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, Int64}})   # time: 0.0020332
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.001983
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0019738
    Base.precompile(Tuple{typeof(zero_on_basis),SingleParticleBasis{SquareLattice{2}}})   # time: 0.0019691
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0019116
    Base.precompile(Tuple{typeof(iterate),LatticeValueRecord})   # time: 0.0018835
    Base.precompile(Tuple{typeof(copy),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0018275
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},Tuple{Int64}})   # time: 0.0018235
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}},))
    end
end   # time: 0.0017903
    Base.precompile(Tuple{typeof(projector),Spectrum{SingleParticleBasis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0017727
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0017449
    Base.precompile(Tuple{typeof(line_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0017181
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{HoneycombLattice}}},))
    end
end   # time: 0.0017105
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0016552
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,SingleParticleBasis{SquareLattice{2}},Tuple{Matrix{Bool}},Tuple{LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Matrix{Bool}},Tuple{LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}},))
    end
end   # time: 0.0016175
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Bool}})   # time: 0.0016114
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0015984
    Base.precompile(Tuple{typeof(line_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0015966
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0015952
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0015835
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.00158
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0015796
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}},Type{Float64}})   # time: 0.0015712
    Base.precompile(Tuple{Type{LatticeTimeSequence{LatticeValue{Float64, :square}}},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.0015691
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Complex{Int64}, Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Complex{Int64}, Float64},))
    end
end   # time: 0.0015643
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{1}}},LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{1}}}})   # time: 0.0015584
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, :square}})   # time: 0.0015296
    Base.precompile(Tuple{PairLhsGraph,HoneycombLattice,LatticeSite{2},LatticeSite{2}})   # time: 0.0015222
    Base.precompile(Tuple{typeof(check_lattice_fits),PairLhsGraph,SquareLattice{2}})   # time: 0.0015083
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector})   # time: 0.0015043
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0015041
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}},))
    end
end   # time: 0.0014861
    Base.precompile(Tuple{typeof(diag_operator),Function,SingleParticleBasis{SquareLattice{2}}})   # time: 0.0014849
    Base.precompile(Tuple{Type{Basis},SquareLattice{1},Int64})   # time: 0.001363
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}})   # time: 0.0013285
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0013153
    Base.precompile(Tuple{Type{LatticeValue},HoneycombLattice,Vector{Bool}})   # time: 0.0012605
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0012384
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0012257
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0011998
    Base.precompile(Tuple{typeof(*),Int64,MaterializedCurrents})   # time: 0.0011953
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Bool}})   # time: 0.0011842
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,SingleParticleBasis{SquareLattice{2}},Tuple{Vector{Float64}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Vector{Float64}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0011623
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,SingleParticleBasis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, SingleParticleBasis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0011614
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0011562
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Tuple{}})   # time: 0.0011543
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.0011499
    Base.precompile(Tuple{typeof(isless),LatticeSite{2},LatticeSite{2}})   # time: 0.0011305
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},ComplexF64,Tuple{LatticeOperator{Matrix{ComplexF64}, SingleParticleBasis{SquareLattice{2}}}}})   # time: 0.0011302
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,SingleParticleBasis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},Tuple{LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}},Tuple{LatticeOperator{Matrix{Float64}, SingleParticleBasis{SquareLattice{2}}}},))
    end
end   # time: 0.0011187
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2},Vector{Bool}})   # time: 0.0011087
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,SingleParticleBasis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{Int64}})   # time: 0.0010513
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0010319
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0010056
end
