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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:reduce_fn, :sort), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0305424
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.4200073
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4136177
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.3439597
    Base.precompile(Tuple{typeof(site_coords),SquareLattice{2},LatticeSite{2}})   # time: 0.3081142
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.298647
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.2290371
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.2275542
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.1904641
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.183018
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.1508282
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1492565
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.1409026
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1408945
    Base.precompile(Tuple{typeof(pade_exp),Matrix{ComplexF64},Int64})   # time: 0.139345
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeSite{2}})   # time: 0.1392032
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.137785
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64,Float64})   # time: 0.1371236
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.1298081
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Symbol})   # time: 0.1270841
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.1134716
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.1127129
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.1112905
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1091272
    Base.precompile(Tuple{Type{LatticeOperator},UniformScaling{Bool},Basis{SquareLattice{2}}})   # time: 0.1052683
    Base.precompile(Tuple{typeof(diff),LatticeValueRecord})   # time: 0.0994838
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.0977625
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.0968306
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.0952186
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Float64})   # time: 0.0904905
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0863897
    Base.precompile(Tuple{typeof(==),LatticeValueRecord,LatticeValueRecord})   # time: 0.0845814
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.0825248
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0806036
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.076743
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0764915
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0728126
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0719597
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LandauField})   # time: 0.0708746
    Base.precompile(Tuple{typeof(integrate),LatticeValueRecord})   # time: 0.0696686
    isdefined(LatticeModels, Symbol("#78#79")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#78#79")),Int64})   # time: 0.0684332
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Float64})   # time: 0.0682128
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0652193
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Any})   # time: 0.0643046
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{Matrix{Int64}},LandauField})   # time: 0.0625055
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.0604758
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2}})   # time: 0.0604303
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Float64}})   # time: 0.0603286
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0600143
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0598343
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0587946
    let fbody = try __lookup_kwbody__(which(map_currents, (typeof(site_distance),DensityCurrents,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (typeof(sum),Bool,typeof(map_currents),typeof(site_distance),DensityCurrents,))
    end
end   # time: 0.0526763
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.0523907
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0522133
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0492345
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0460139
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.045719
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0452323
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0426816
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.0419356
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.0415687
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Nothing,Bool})   # time: 0.0409249
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Matrix{Int64}})   # time: 0.0405652
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeValue{Bool, :square}})   # time: 0.0400744
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0382832
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0365286
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0343501
    Base.precompile(Tuple{typeof(+),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0341939
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0341545
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0337949
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0335027
    Base.precompile(Tuple{typeof(diag_operator),Basis{SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0329977
    Base.precompile(Tuple{typeof(*),LatticeOperator{Adjoint{Float64, Vector{Float64}}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0327886
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0323538
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{Matrix{Float64}}})   # time: 0.0317191
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.0311939
    Base.precompile(Tuple{typeof(taylor_exp),Matrix{ComplexF64},Int64})   # time: 0.0308781
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0305511
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0305472
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0272015
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.0267144
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Int64,Bool})   # time: 0.0265174
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0263469
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0263273
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Matrix{Float64},Int64,Int64})   # time: 0.0259609
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0255934
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVWStyle})   # time: 0.0236212
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0235625
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0222218
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.0221396
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0218056
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0206679
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0195602
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.0192415
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite{2}})   # time: 0.0180436
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0160806
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0154891
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0152837
    isdefined(LatticeModels, Symbol("#94#95")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#94#95")),Float64})   # time: 0.01513
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0151025
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:pade, :k), Tuple{Bool, Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0132633
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.0131356
    isdefined(LatticeModels, Symbol("#54#55")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#54#55")),FluxField})   # time: 0.0131313
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0130919
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0116185
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0114797
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0113908
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0112721
    Base.precompile(Tuple{LatticeValueRecord,Float64})   # time: 0.011249
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0105496
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0104322
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0099677
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0099591
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}}})   # time: 0.0099139
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0096119
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0095714
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0091281
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.0091116
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0089504
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.008915
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.008658
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0083743
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0082844
    Base.precompile(Tuple{typeof(==),LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1},LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1}})   # time: 0.0080611
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0080298
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0080118
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0079825
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0079446
    isdefined(LatticeModels, Symbol("#100#102")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#100#102")),Vector{Float64}})   # time: 0.0077905
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0075326
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0074838
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.007359
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0073564
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0073518
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0071135
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{HoneycombLattice}})   # time: 0.006665
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0059966
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.0059375
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0057613
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0057536
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0056478
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0053752
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0052044
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0051188
    Base.precompile(Tuple{typeof(-),LatticeOperator{Matrix{Float64}, SquareLattice{2}},UniformScaling{Bool}})   # time: 0.0049062
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0048479
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0048185
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0047895
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0046262
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Float64, :square}})   # time: 0.0046115
    isdefined(LatticeModels, Symbol("#100#102")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#100#102")),Array})   # time: 0.0045864
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.00446
    Base.precompile(Tuple{typeof(diag_operator),Function,SquareLattice{2}})   # time: 0.0044316
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0044246
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0044224
    Base.precompile(Tuple{typeof(insert!),LatticeValueRecord,Int64,LatticeValue{Float64, :square}})   # time: 0.0042325
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0042078
    Base.precompile(Tuple{typeof(Base.maybeview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.004125
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.0039157
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0039067
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0038691
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.003846
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0038372
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0037923
    Base.precompile(Tuple{typeof(_expr_depends_on),Expr,Symbol})   # time: 0.0037635
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0035036
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.0034724
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Function})   # time: 0.003246
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0031799
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0030798
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0030643
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0030533
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0030511
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0029239
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.0028671
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0027499
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0027451
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Function})   # time: 0.0026542
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0025826
    isdefined(LatticeModels, Symbol("#101#103")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#101#103")),Tuple{Float64, Vector{Float64}}})   # time: 0.0024524
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0023984
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.002381
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}, LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0023788
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0023551
    Base.precompile(Tuple{typeof(copy),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0023481
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0023305
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.002255
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0022383
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64,Int64})   # time: 0.0022311
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}},Type{Float64}})   # time: 0.0022289
    isdefined(LatticeModels, Symbol("#22#23")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#22#23")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0021566
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0020914
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{Float64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0020716
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{Vector{Bool}, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, CartesianIndex{3}, Int64}})   # time: 0.0020533
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0020178
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}}})   # time: 0.0020152
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0020093
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0019954
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0019846
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0019825
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0019799
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0019707
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, :square}})   # time: 0.0019179
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Int64}})   # time: 0.0018489
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0018374
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Matrix{Float64}})   # time: 0.0017798
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0017671
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0017603
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Union{Nothing, Int64}}})   # time: 0.0017531
    Base.precompile(Tuple{typeof(randn),SquareLattice{2}})   # time: 0.0017528
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}}})   # time: 0.0017517
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))
    end
end   # time: 0.0017396
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0016879
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0016692
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Vararg{Any}})   # time: 0.0016384
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0016079
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0015891
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0015685
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0015509
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0015181
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0014823
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.001466
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0014091
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.001406
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.001362
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0013437
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))
    end
end   # time: 0.0013321
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0013284
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Complex{Int64}, Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Complex{Int64}, Float64},))
    end
end   # time: 0.0013027
    Base.precompile(Tuple{typeof(init_record),LatticeValue{Float64, :square}})   # time: 0.0012895
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Nothing}})   # time: 0.0012599
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0012328
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0012297
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Tuple{Int64}})   # time: 0.0012001
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.001168
    Base.precompile(Tuple{typeof(_lazy_tp),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0011421
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))
    end
end   # time: 0.0011222
    Base.precompile(Tuple{typeof(check_lattice_fits),PairLhsSelector,SquareLattice{2}})   # time: 0.0010496
    Base.precompile(Tuple{PairLhsSelector,SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0010318
    Base.precompile(Tuple{typeof(isless),LatticeSite{2},LatticeSite{2}})   # time: 0.0010302
end
