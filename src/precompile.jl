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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:aggr_fn, :sorted), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0983264
    Base.precompile(Tuple{typeof(exp),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.8518205
    Base.precompile(Tuple{typeof(^),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.4695441
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4639641
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.3723907
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.3721642
    Base.precompile(Tuple{typeof(site_coords),SquareLattice{2},LatticeSite{2}})   # time: 0.3526428
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.2745803
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.2548061
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.2245743
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.1891651
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1874825
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.186154
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64,Float64})   # time: 0.1710599
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.1636907
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1531998
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1499515
    Base.precompile(Tuple{Type{LatticeOperator},UniformScaling{Bool},Basis{SquareLattice{2}}})   # time: 0.1303069
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.1249843
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.1220025
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.1104754
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.1078105
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.1043542
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.087994
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0869468
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0852528
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{Matrix{Int64}},LandauField})   # time: 0.0842324
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Float64})   # time: 0.0818659
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0800794
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LandauField})   # time: 0.0756355
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0730712
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0722998
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.072118
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Any})   # time: 0.0699906
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0683587
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0641868
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.063594
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeValue{Float64, :square}})   # time: 0.0627227
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0598705
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0595419
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0563742
    isdefined(LatticeModels, Symbol("#66#67")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#66#67")),Int64})   # time: 0.0537684
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0525943
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0522731
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0513941
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0466321
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Adjoint{Float64, Vector{Float64}}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeArray{SquareLattice{2}, Vector{Float64}}})   # time: 0.0433748
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Matrix{Int64}})   # time: 0.0427449
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0423419
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0418451
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0415956
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0408739
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0386114
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0369609
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0364851
    Base.precompile(Tuple{typeof(taylor_exp),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64})   # time: 0.0363762
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.0354586
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.033647
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0333626
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.031753
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0314589
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0312025
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0310258
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0305951
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0296091
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0293288
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{Matrix{Float64}}})   # time: 0.0291221
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.029018
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.028973
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVStyle})   # time: 0.0283534
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0246921
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0214591
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0202496
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.0202418
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.0198423
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0197978
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0173334
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0172365
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0171909
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Float64,Nothing})   # time: 0.0170846
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.0161207
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0153449
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0146911
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0146117
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.014434
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0143589
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0141667
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Float64},SVector{2, Float64}})   # time: 0.0140514
    isdefined(LatticeModels, Symbol("#82#83")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#82#83")),Float64})   # time: 0.0139961
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0135531
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0134627
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.0134349
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Float64,Int64})   # time: 0.0127026
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0125759
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.012108
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0116882
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.0114893
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0113728
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0110888
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0104181
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0103808
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0100046
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{HoneycombLattice}})   # time: 0.0098253
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0094244
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0092094
    isdefined(LatticeModels, Symbol("#46#47")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#46#47")),FluxField})   # time: 0.009036
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Int64,Int64})   # time: 0.0090146
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.0089892
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0089389
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0087875
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0085895
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}}})   # time: 0.0085768
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0083684
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{SquareLattice{2}, Vector{Float64}}})   # time: 0.0081932
    Base.precompile(Tuple{typeof(+),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0080873
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0077468
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.00748
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0074721
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0074371
    Base.precompile(Tuple{typeof(==),LatticeArray{HoneycombLattice, Vector{ComplexF64}},LatticeArray{HoneycombLattice, Vector{ComplexF64}}})   # time: 0.0068516
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0066889
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.006593
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0062228
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0061463
    Base.precompile(Tuple{typeof(diag_operator),Basis{SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.006086
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0059866
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0057531
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},Function})   # time: 0.0055563
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0053384
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, Matrix{Float64}},UniformScaling{Bool}})   # time: 0.0052903
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0052127
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0050371
    Base.precompile(Tuple{PairLhsSelector,SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0050281
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.005009
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0049572
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0047654
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.0046183
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0045061
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0043777
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0043583
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0042866
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0042466
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0041376
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0040898
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0039784
    Base.precompile(Tuple{typeof(copy),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0037171
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.0036089
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0036026
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0035272
    Base.precompile(Tuple{typeof(/),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0034165
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.003404
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0032817
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0029727
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0029477
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0029368
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0028501
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.0028385
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0027187
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0026655
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}}})   # time: 0.0026064
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.002436
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{Vector{Bool}, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, CartesianIndex{3}, Int64}})   # time: 0.0023752
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.0023597
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.002358
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{Float64}},Vararg{LatticeOperator{SquareLattice{2}, Matrix{Float64}}}})   # time: 0.002301
    Base.precompile(Tuple{typeof(coord_operators),SquareLattice{2},Int64})   # time: 0.0022851
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0022511
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.0022172
    isdefined(LatticeModels, Symbol("#20#21")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#20#21")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0022082
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0021559
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},))
    end
end   # time: 0.0021015
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Adjoint{ComplexF64, SubArray{ComplexF64, 2, Matrix{ComplexF64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}},Int64,Int64})   # time: 0.002004
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}},Type{Float64}})   # time: 0.0020023
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0019756
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Vararg{Any}})   # time: 0.0019398
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0019298
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.001928
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.001884
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Function})   # time: 0.0018465
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0018129
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0018127
    Base.precompile(Tuple{typeof(_extract_lattice_s),SquareLattice{2},SquareLattice{2},Tuple{}})   # time: 0.001797
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0017602
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0017238
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}, LatticeArray{SquareLattice{2}, Vector{Float64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}, LatticeArray{SquareLattice{2}, Vector{Float64}}},))
    end
end   # time: 0.001649
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))
    end
end   # time: 0.0016048
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0015715
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0015657
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0015461
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0015122
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0015094
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Vararg{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}}})   # time: 0.0014736
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0014733
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}})   # time: 0.0014101
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0013577
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Int64,Tuple{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Int64,Tuple{Int64},))
    end
end   # time: 0.0013285
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64,Int64})   # time: 0.0013137
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Tuple{}})   # time: 0.0012932
    Base.precompile(Tuple{typeof(randn),SquareLattice{2}})   # time: 0.0012928
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.001255
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0012442
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.00124
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0012145
    isdefined(LatticeModels, Symbol("#51#52")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#51#52")),Symbol})   # time: 0.0011937
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))
    end
end   # time: 0.0011631
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0011441
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.001138
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0011278
    Base.precompile(Tuple{typeof(view),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64,Int64})   # time: 0.0011042
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0011032
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0010907
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0010389
end
