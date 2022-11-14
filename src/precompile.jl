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
    Base.precompile(Tuple{typeof(exp),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.5584582
    Base.precompile(Tuple{typeof(^),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.470652
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4244771
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.3543766
    Base.precompile(Tuple{typeof(^),BondSet{SquareLattice{2}},Int64})   # time: 0.3517945
    Base.precompile(Tuple{typeof(site_coords),SquareLattice{2},LatticeSite{2}})   # time: 0.3422401
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.2532276
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.2419901
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.2411334
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.2258302
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.1775341
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1735842
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.1697463
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},SquareLattice{2},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.1666825
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1583045
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{_A, Matrix{ComplexF64}} where _A<:Lattice,SquareLattice{2},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.1506742
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.1446634
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{_A, Matrix{ComplexF64}} where _A<:Lattice,SquareLattice{2},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.1377605
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},SquareLattice{2},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.1295377
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1287129
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.1117757
    Base.precompile(Tuple{Type{LatticeOperator},UniformScaling{Bool},Basis{SquareLattice{2}}})   # time: 0.1003606
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LandauField})   # time: 0.0985253
    Base.precompile(Tuple{typeof(_diag_from_macro),SquareLattice{2},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.0929339
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0891972
    Base.precompile(Tuple{typeof(trip_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.088906
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0873567
    Base.precompile(Tuple{typeof(|),BondSet{SquareLattice{2}},BondSet{SquareLattice{2}}})   # time: 0.0821412
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0792278
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0730827
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Any})   # time: 0.0691874
    Base.precompile(Tuple{typeof(site_index),LatticeSite{2},SquareLattice{2}})   # time: 0.0663705
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0618636
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0615309
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0602833
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0566039
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.055669
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.051319
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0509149
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.0502475
    Base.precompile(Tuple{typeof(_diag_from_macro),SquareLattice{2},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.0500541
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0475025
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0413073
    Base.precompile(Tuple{Core.kwftype(typeof(trip_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(trip_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0385175
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0381245
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0376275
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0375882
    Base.precompile(Tuple{typeof(_hops_from_macro),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},HoneycombLattice,Function,Hopping{Matrix{Int64}},LandauField})   # time: 0.0366356
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0362062
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0360518
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0354093
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0348526
    Base.precompile(Tuple{typeof(taylor_exp),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64})   # time: 0.0338588
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.0334828
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0332396
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0312264
    Base.precompile(Tuple{Core.kwftype(typeof(trip_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(trip_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0309801
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0307969
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVStyle})   # time: 0.0306759
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0306572
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0303516
    Base.precompile(Tuple{typeof(_diag_from_macro),SquareLattice{2},Matrix{Int64}})   # time: 0.0303407
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},BondSet})   # time: 0.027145
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0263852
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0234399
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0234259
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0206915
    Base.precompile(Tuple{Core.kwftype(typeof(trip_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(trip_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0179207
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Float64},SVector{2, Float64}})   # time: 0.0177258
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0174036
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.016435
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0141907
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0138324
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.0128915
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.012266
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Float64,Int64})   # time: 0.0122081
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0115875
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0104925
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0102759
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0101158
    Base.precompile(Tuple{Core.kwftype(typeof(trip_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(trip_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0098445
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0097051
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0096645
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0094728
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0094237
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0093867
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,Int64})   # time: 0.0091911
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.009144
    Base.precompile(Tuple{typeof(_diag_from_macro),HoneycombLattice,TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0090036
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0089203
    isdefined(LatticeModels, Symbol("#46#47")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#46#47")),FluxField})   # time: 0.0088167
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0087767
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0080891
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0080799
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Int64,Int64})   # time: 0.0076148
    Base.precompile(Tuple{typeof(==),LatticeArray{HoneycombLattice, Vector{ComplexF64}},LatticeArray{HoneycombLattice, Vector{ComplexF64}}})   # time: 0.0075334
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0075202
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0073723
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0071156
    Base.precompile(Tuple{typeof(_diag_from_macro),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},HoneycombLattice,TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0070535
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0068536
    Base.precompile(Tuple{typeof(diag_operator),Basis{SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0065793
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, Matrix{Float64}},UniformScaling{Bool}})   # time: 0.0063163
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0062944
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.006035
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Vector{Bool}}},typeof(hopping)})   # time: 0.0056529
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.005248
    Base.precompile(Tuple{typeof(trip_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0051802
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0051528
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0046206
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0045697
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.004444
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0044367
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0043262
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0041738
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0041684
    Base.precompile(Tuple{typeof(_zero_on_basis),HoneycombLattice,Matrix{Int64}})   # time: 0.004166
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping),Matrix{Int64}})   # time: 0.004109
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0039419
    Base.precompile(Tuple{typeof(_extract_lattice_s),HoneycombLattice,HoneycombLattice,Tuple{}})   # time: 0.0038484
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0038012
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Function})   # time: 0.0036282
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0033693
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.00336
    Base.precompile(Tuple{typeof(_diag_from_macro),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},SquareLattice{2},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.0032765
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0032067
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0031609
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0030868
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0029815
    Base.precompile(Tuple{typeof(/),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0027967
    Base.precompile(Tuple{typeof(_diag_from_macro),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},HoneycombLattice,Function})   # time: 0.0027603
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{Float64}},Vararg{LatticeOperator{SquareLattice{2}, Matrix{Float64}}}})   # time: 0.0027336
    Base.precompile(Tuple{typeof(_extract_lattice_s),SquareLattice{2},SquareLattice{2},Tuple{}})   # time: 0.0027005
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0026532
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{Vector{Bool}, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, CartesianIndex{3}, Int64}})   # time: 0.0025933
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}}})   # time: 0.002587
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0025256
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0021533
    Base.precompile(Tuple{typeof(is_adjacent),BondSet{SquareLattice{2}},LatticeSite{2},LatticeSite{2}})   # time: 0.002015
    Base.precompile(Tuple{typeof(copy),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0019711
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64,Int64})   # time: 0.0019011
    isdefined(LatticeModels, Symbol("#18#19")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#18#19")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0017217
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0017004
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0016949
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0016922
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.0016876
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}})   # time: 0.0016839
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0016786
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Vararg{Any}})   # time: 0.0016765
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0016636
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0016419
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.0016173
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0015555
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0015204
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))
    end
end   # time: 0.0014874
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0014503
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0014482
    Base.precompile(Tuple{typeof(trip_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0014397
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Vararg{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}}})   # time: 0.0014345
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0014217
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0014137
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},))
    end
end   # time: 0.0013993
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}})   # time: 0.0013789
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0013428
    Base.precompile(Tuple{typeof(trip_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0013421
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0013371
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0013096
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0013057
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0012134
    isdefined(LatticeModels, Symbol("#51#52")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#51#52")),Symbol})   # time: 0.001212
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))
    end
end   # time: 0.0011901
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{SquareLattice{2}}})   # time: 0.0011817
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0011388
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0011349
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
    end
end   # time: 0.0010763
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0010294
end
