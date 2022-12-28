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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:aggr_fn, :sorted), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0353428
    Base.precompile(Tuple{typeof(exp),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.5261819
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.4594854
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4326023
    Base.precompile(Tuple{typeof(site_coords),SquareLattice{2},LatticeSite{2}})   # time: 0.3150098
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.307991
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.2978594
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.2168073
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.2065066
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.1604156
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.158693
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.1511231
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1321561
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Symbol})   # time: 0.1302176
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1251539
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64,Float64})   # time: 0.1251182
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.124457
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.1230046
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1146487
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.1129842
    Base.precompile(Tuple{typeof(diff),LatticeValueRecord})   # time: 0.1089987
    Base.precompile(Tuple{Type{LatticeOperator},UniformScaling{Bool},Basis{SquareLattice{2}}})   # time: 0.106829
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.097934
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.091587
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.089772
    Base.precompile(Tuple{typeof(ldos),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Float64})   # time: 0.0852989
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0812867
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0784026
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.0780314
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0779161
    Base.precompile(Tuple{typeof(integrate),LatticeValueRecord})   # time: 0.0759631
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0750847
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.0748542
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0742839
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Float64}})   # time: 0.073299
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.0669591
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0660402
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0649245
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0644026
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{Matrix{Int64}},LandauField})   # time: 0.0633882
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0617555
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeSite{2}})   # time: 0.0610472
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Any})   # time: 0.0605902
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LandauField})   # time: 0.0594187
    Base.precompile(Tuple{typeof(==),LatticeValueRecord,LatticeValueRecord})   # time: 0.0570546
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0558829
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.053947
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0535916
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.0534228
    isdefined(LatticeModels, Symbol("#76#77")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#76#77")),Int64})   # time: 0.0527535
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0482815
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0468298
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0468057
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.043617
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0411642
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0391781
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0384359
    Base.precompile(Tuple{typeof(*),LatticeOperator{Adjoint{Float64, Vector{Float64}}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0381709
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0380677
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.0377113
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0373308
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0372126
    Base.precompile(Tuple{typeof(+),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0367157
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0363964
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Matrix{Int64}})   # time: 0.0360237
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.0350687
    Base.precompile(Tuple{typeof(taylor_exp),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64})   # time: 0.0344572
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0335451
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0332106
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.032848
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0323973
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{Matrix{Float64}}})   # time: 0.0318608
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0317084
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.0316575
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0310538
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0300498
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0294716
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0293229
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.029262
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0280632
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2}})   # time: 0.0247273
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVWStyle})   # time: 0.023998
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0214524
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0212855
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0194381
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.0188963
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0184997
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite{2}})   # time: 0.018461
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0168866
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0166774
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0161927
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0158987
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0157756
    isdefined(LatticeModels, Symbol("#92#93")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#92#93")),Float64})   # time: 0.0152589
    Base.precompile(Tuple{LatticeValueRecord,Float64})   # time: 0.0150906
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.0150636
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeValue{Bool, :square}})   # time: 0.0149777
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.01458
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0130236
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0127726
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Int64})   # time: 0.0127296
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0123388
    isdefined(LatticeModels, Symbol("#52#53")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#52#53")),FluxField})   # time: 0.0122604
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0120865
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0117921
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0115011
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0114162
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0112803
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.0110982
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0109401
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.0108297
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0104619
    isdefined(LatticeModels, Symbol("#96#98")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#96#98")),Vector{Float64}})   # time: 0.0104567
    isdefined(LatticeModels, Symbol("#97#99")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#97#99")),Array})   # time: 0.0103205
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}}})   # time: 0.009416
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0093598
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Nothing})   # time: 0.0092647
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0092084
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0092013
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0088452
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0088351
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0085185
    Base.precompile(Tuple{typeof(==),LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1},LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1}})   # time: 0.0084206
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0083648
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0082898
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0080382
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0077797
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0077789
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Matrix{Float64},Int64,Int64})   # time: 0.0076875
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{HoneycombLattice}})   # time: 0.0076824
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0074011
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0071717
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0070345
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0067561
    Base.precompile(Tuple{typeof(diag_operator),Basis{SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0067202
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0058142
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0057495
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0057272
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0056738
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.005369
    Base.precompile(Tuple{typeof(-),LatticeOperator{Matrix{Float64}, SquareLattice{2}},UniformScaling{Bool}})   # time: 0.0052007
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0051498
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0049295
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0048112
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0047913
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0046886
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0045757
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Function})   # time: 0.0045581
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.004495
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0042278
    Base.precompile(Tuple{typeof(Base.maybeview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0041942
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.0041139
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0040473
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.003895
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.0038749
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0037147
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :honeycomb}}},Type{Float64}})   # time: 0.0036601
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0035921
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}, LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0035297
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0034663
    Base.precompile(Tuple{typeof(insert!),LatticeValueRecord,Int64,LatticeValue{Float64, :square}})   # time: 0.0033674
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.003291
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.003246
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0032266
    Base.precompile(Tuple{typeof(/),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.0032209
    isdefined(LatticeModels, Symbol("#97#99")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#97#99")),Vector{Float64}})   # time: 0.0031019
    isdefined(LatticeModels, Symbol("#96#98")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#96#98")),Array})   # time: 0.0030686
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.0030207
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0029784
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0029683
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0029526
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.0029
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.0027523
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0026875
    isdefined(LatticeModels, Symbol("#20#21")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#20#21")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0026548
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0026473
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.002633
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0025641
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0025442
    Base.precompile(Tuple{typeof(init_record),LatticeValue{Float64, :square}})   # time: 0.0025259
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0023051
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0022983
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0022473
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0022398
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0021605
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0021218
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Int64}})   # time: 0.0021218
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0020649
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0020523
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{Float64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0020244
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0019813
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0019285
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Function})   # time: 0.0019277
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, :square}})   # time: 0.0019157
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{Vector{Bool}, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, CartesianIndex{3}, Int64}})   # time: 0.0019001
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0018845
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0018658
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0018108
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0017986
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0017891
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0017234
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))
    end
end   # time: 0.0017066
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.001697
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0016943
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Vararg{Any}})   # time: 0.0016942
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0016845
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0016605
    Base.precompile(Tuple{typeof(copy),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0015586
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))
    end
end   # time: 0.001548
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0015069
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}}})   # time: 0.0014788
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0014759
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{Int64}})   # time: 0.001455
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0014511
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0014465
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0013889
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}}})   # time: 0.0013872
    Base.precompile(Tuple{typeof(randn),SquareLattice{2}})   # time: 0.0013862
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64,Int64})   # time: 0.0013093
    Base.precompile(Tuple{Type{LatticeRecord{LatticeValue{Float64, :square}}},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.0012665
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0012551
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0012512
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Union{Nothing, Int64}}})   # time: 0.0012335
    Base.precompile(Tuple{PairLhsSelector,SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.001231
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Nothing}})   # time: 0.0012304
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0012209
    isdefined(LatticeModels, Symbol("#57#58")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#57#58")),Symbol})   # time: 0.0012181
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0011928
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0011736
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Float64}})   # time: 0.0011611
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Adjoint{ComplexF64, SubArray{ComplexF64, 2, Matrix{ComplexF64}, Tuple{UnitRange{Int64}, UnitRange{Int64}}, false}},Int64,Int64})   # time: 0.0011298
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0011283
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))
    end
end   # time: 0.0011213
    Base.precompile(Tuple{typeof(check_lattice_fits),PairLhsSelector,SquareLattice{2}})   # time: 0.0011192
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0011138
    Base.precompile(Tuple{typeof(view),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64,Int64})   # time: 0.0011136
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Tuple{}})   # time: 0.0010904
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0010803
    Base.precompile(Tuple{typeof(isless),LatticeSite{2},LatticeSite{2}})   # time: 0.0010429
end
