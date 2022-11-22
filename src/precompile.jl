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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:aggr_fn, :sorted), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 0.9653658
    Base.precompile(Tuple{typeof(exp),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4704647
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.4217755
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.3158322
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.2597242
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.253144
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.1955212
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.1696183
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.1488565
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.1419429
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.1362987
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1361772
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.0960598
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.0889045
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LandauField})   # time: 0.0878909
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0760401
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0719688
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0699118
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.0659278
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.064317
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0599529
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Any})   # time: 0.0593863
    Base.precompile(Tuple{typeof(site_index),LatticeSite{2},SquareLattice{2}})   # time: 0.0581251
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0562383
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0519945
    Base.precompile(Tuple{typeof(^),BondSet{SquareLattice{2}},Int64})   # time: 0.051596
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0510117
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeValue{Float64, :square}})   # time: 0.0506949
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0505968
    Base.precompile(Tuple{typeof(_hamiltonian_block),Expr})   # time: 0.0492173
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.0473351
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.046888
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0445491
    Base.precompile(Tuple{Type{LatticeOperator},UniformScaling{Bool},Basis{SquareLattice{2}}})   # time: 0.0434101
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0384439
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0360132
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Float64,Int64})   # time: 0.0358036
    Base.precompile(Tuple{typeof(!),BondSet{SquareLattice{2}}})   # time: 0.0355855
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0353595
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{Matrix{Int64}},LandauField})   # time: 0.0352357
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0335429
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0333541
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0330742
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0328509
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0322773
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0321802
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Matrix{Int64}})   # time: 0.0316987
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0316721
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0306554
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0306408
    Base.precompile(Tuple{typeof(+),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0285418
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0272195
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0257477
    Base.precompile(Tuple{typeof(taylor_exp),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64})   # time: 0.0253573
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},BondSet})   # time: 0.0245803
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0245611
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0240718
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVStyle})   # time: 0.02342
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0221222
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0189989
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0184493
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1},Vector{Bool}})   # time: 0.0171102
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0166423
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0165699
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0156014
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.0141427
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type})   # time: 0.0135693
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{Matrix{Float64}}})   # time: 0.0133555
    Base.precompile(Tuple{typeof(diag_aggregate),typeof(tr),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0127499
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.012614
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.011542
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0109437
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Float64},SVector{2, Float64}})   # time: 0.0109
    Base.precompile(Tuple{typeof(_zero_on_basis),HoneycombLattice,Int64,Type})   # time: 0.0107505
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.0099512
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.0097696
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0092302
    Base.precompile(Tuple{Core.kwftype(typeof(path_integral)),NamedTuple{(:n_steps,), Tuple{Int64}},typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.009133
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0090704
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,Int64})   # time: 0.0089666
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0086644
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0083132
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0082287
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0081208
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0078878
    Base.precompile(Tuple{typeof(==),LatticeArray{HoneycombLattice, Vector{ComplexF64}},LatticeArray{HoneycombLattice, Vector{ComplexF64}}})   # time: 0.0078369
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0077733
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0077567
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0075606
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0075248
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0073776
    Base.precompile(Tuple{typeof(==),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0073705
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.006905
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0065401
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0058159
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{HoneycombLattice}})   # time: 0.0057425
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Matrix{Float64},Int64,Int64})   # time: 0.0054593
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0053184
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0051039
    Base.precompile(Tuple{typeof(diag_operator),Basis{SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0050944
    Base.precompile(Tuple{typeof(-),LatticeOperator{SquareLattice{2}, Matrix{Float64}},UniformScaling{Bool}})   # time: 0.0049664
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0048434
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0047748
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0045598
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0045293
    Base.precompile(Tuple{typeof(path_integral),FieldSum{2},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0043498
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0042493
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0042147
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.004108
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0040825
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0040434
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.0039916
    Base.precompile(Tuple{typeof(|),BondSet{SquareLattice{2}},BondSet{SquareLattice{2}}})   # time: 0.0037295
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}})   # time: 0.0036305
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{HoneycombLattice, Matrix{ComplexF64}},Function})   # time: 0.0034643
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}}})   # time: 0.0029469
    Base.precompile(Tuple{typeof(*),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0028081
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}}})   # time: 0.0027893
    Base.precompile(Tuple{typeof(_extract_lattice_s),SquareLattice{2},SquareLattice{2},Tuple{}})   # time: 0.0027665
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0027302
    Base.precompile(Tuple{typeof(+),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0026148
    Base.precompile(Tuple{typeof(_extract_lattice_s),HoneycombLattice,HoneycombLattice,Tuple{}})   # time: 0.0025948
    Base.precompile(Tuple{typeof(copy),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0025861
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0025624
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0025093
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0023563
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0023536
    Base.precompile(Tuple{typeof(/),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0022978
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{Float64}},Vararg{LatticeOperator{SquareLattice{2}, Matrix{Float64}}}})   # time: 0.0022781
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Union{Nothing, Int64}}})   # time: 0.0021685
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.002098
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0020716
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0020706
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Vararg{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}}})   # time: 0.0020475
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0020051
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Function})   # time: 0.0018865
    Base.precompile(Tuple{typeof(^),LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Int64})   # time: 0.0018468
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}})   # time: 0.0018366
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{SquareLattice{2}}})   # time: 0.001784
    Base.precompile(Tuple{typeof(is_adjacent),BondSet{SquareLattice{2}},LatticeSite{2},LatticeSite{2}})   # time: 0.0017781
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.0017777
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Vararg{Any}})   # time: 0.0017533
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0017332
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0017198
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0017081
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0016826
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
        end
    end   # time: 0.0016472
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0016121
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}},))
        end
    end   # time: 0.0014993
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0014572
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0014448
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0014231
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0013843
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0013711
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
        end
    end   # time: 0.0013596
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Int64,Int64})   # time: 0.0013299
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0012979
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
        end
    end   # time: 0.0012937
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0012891
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.001205
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0012002
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0011817
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Float64}})   # time: 0.001106
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{SquareLattice{2}, Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))
        end
    end   # time: 0.0010975
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.0010968
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0010869
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0010575
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{LatticeOperator{SquareLattice{2}, SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}})   # time: 0.0010517
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{ComplexF64}}},))
        end
    end   # time: 0.001021
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0010189
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))) catch missing end
        if !ismissing(fbody)
            precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{SquareLattice{2}, Matrix{Float64}}, LatticeOperator{SquareLattice{2}, Matrix{Float64}}},))
        end
    end   # time: 0.0010066
end
