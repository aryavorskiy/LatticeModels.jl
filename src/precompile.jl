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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:aggr_fn, :sorted), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0633302
    Base.precompile(Tuple{typeof(exp),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.5750092
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4521567
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.4513247
    Base.precompile(Tuple{typeof(site_coords),SquareLattice{2},LatticeSite{2}})   # time: 0.3359604
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.3222846
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.307499
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.2278221
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.2195619
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.2041106
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.1682073
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1548712
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.140857
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Symbol})   # time: 0.136401
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1338562
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64,Float64})   # time: 0.1273172
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.1192321
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1186837
    isdefined(LatticeModels, Symbol("#76#77")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#76#77")),Int64})   # time: 0.1158286
    Base.precompile(Tuple{Type{LatticeOperator},UniformScaling{Bool},Basis{SquareLattice{2}}})   # time: 0.1119948
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.1035008
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{Matrix{Int64}},LandauField})   # time: 0.1019112
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.099893
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Any})   # time: 0.0967942
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.0957148
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.0935687
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.0921084
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.0913804
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0846888
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0808721
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0796038
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0782803
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0768116
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Float64})   # time: 0.0766956
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.0763451
    Base.precompile(Tuple{typeof(+),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0750568
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0742071
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0711431
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0706146
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.070158
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LandauField})   # time: 0.0681859
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0669078
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0603449
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0591372
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.0576358
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0507199
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0501967
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0481567
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0462475
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0453685
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.0446499
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0445114
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0390745
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.038606
    Base.precompile(Tuple{typeof(diag_operator),Basis{SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.03816
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0377132
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Matrix{Int64}})   # time: 0.0374982
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0373836
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0370007
    Base.precompile(Tuple{typeof(taylor_exp),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64})   # time: 0.0366132
    Base.precompile(Tuple{typeof(*),LatticeOperator{Adjoint{Float64, Vector{Float64}}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0349825
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0344955
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.0343447
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0330376
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0322957
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0321685
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0293778
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0293653
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0290281
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.028985
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{Matrix{Float64}}})   # time: 0.028263
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVWStyle})   # time: 0.0273783
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Matrix{Float64},Int64,Int64})   # time: 0.0273112
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0267254
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.02668
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.02656
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0242642
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0238334
    isdefined(LatticeModels, Symbol("#14#15")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#14#15")),Int64})   # time: 0.0237461
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0232993
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0219952
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0194728
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0191543
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.01914
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite{2}})   # time: 0.0182753
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0178993
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.0169304
    isdefined(LatticeModels, Symbol("#92#93")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#92#93")),Float64})   # time: 0.0166516
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0164778
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0162907
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0144665
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0139085
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0129469
    isdefined(LatticeModels, Symbol("#52#53")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#52#53")),FluxField})   # time: 0.0127045
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0123414
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}}})   # time: 0.0121262
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Nothing})   # time: 0.0118967
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0117804
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0117411
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0114341
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Int64})   # time: 0.0111649
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.0109522
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0108405
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{HoneycombLattice}})   # time: 0.0106688
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0106654
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0100703
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0097436
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0096934
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0095826
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0094838
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0093391
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.009286
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0092579
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0089676
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0084178
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0083011
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.0082691
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0081425
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0081313
    Base.precompile(Tuple{typeof(==),LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1},LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1}})   # time: 0.0079459
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0072302
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0071792
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.00696
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0065917
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0065274
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0059432
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0059253
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0058914
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0057847
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.005539
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0052576
    Base.precompile(Tuple{typeof(-),LatticeOperator{Matrix{Float64}, SquareLattice{2}},UniformScaling{Bool}})   # time: 0.004973
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0048025
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0047576
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0047293
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0046631
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0046514
    Base.precompile(Tuple{typeof(Base.maybeview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0045754
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0045657
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0045613
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0045203
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0044929
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.0043935
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.004327
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0042085
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0039235
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.0038527
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.003808
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.0037024
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0036688
    Base.precompile(Tuple{PairLhsSelector,SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.00362
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0035531
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Function})   # time: 0.0035487
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0035311
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0032869
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0032425
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0031462
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0030633
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0030447
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{Vector{Bool}, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, CartesianIndex{3}, Int64}})   # time: 0.0029745
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0029341
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.002845
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}, LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0028419
    isdefined(LatticeModels, Symbol("#20#21")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#20#21")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0027937
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.0027258
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0027238
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0026367
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0025405
    Base.precompile(Tuple{typeof(/),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.0025308
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Function})   # time: 0.0024907
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0024845
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0024697
    Base.precompile(Tuple{typeof(view),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64,Int64})   # time: 0.0024141
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0023993
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Vararg{Any}})   # time: 0.0023327
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{Float64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0023018
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0022791
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0022253
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0021269
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0020584
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0020484
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0020353
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0020294
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}},Type{Float64}})   # time: 0.0020117
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0019941
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, :square}})   # time: 0.0019235
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))
    end
end   # time: 0.0018326
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0018324
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0018223
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0017873
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0017846
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.001777
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0017483
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0017386
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0017221
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.0016885
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0016868
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}}})   # time: 0.001638
    Base.precompile(Tuple{typeof(copy),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0016344
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}}})   # time: 0.0015632
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Union{Nothing, Int64}}})   # time: 0.0015548
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0015224
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0014632
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0014478
    Base.precompile(Tuple{typeof(randn),SquareLattice{2}})   # time: 0.0014239
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0013999
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0013812
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))
    end
end   # time: 0.0013587
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0013256
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0013079
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Nothing}})   # time: 0.0012984
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64,Int64})   # time: 0.0012938
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0012502
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Float64},Tuple{Int64, Float64}})   # time: 0.0011985
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0011984
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :square}})   # time: 0.0011657
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0011474
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Adjoint{ComplexF64, SubArray{ComplexF64, 2, Matrix{ComplexF64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}},Int64,Int64})   # time: 0.0011
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0010861
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))
    end
end   # time: 0.0010854
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Float64}})   # time: 0.001085
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0010368
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.001034
    Base.precompile(Tuple{typeof(isless),LatticeSite{2},LatticeSite{2}})   # time: 0.0010332
end
