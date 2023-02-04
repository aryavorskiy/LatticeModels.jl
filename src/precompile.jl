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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:reduce_fn, :sort), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0197494
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.4488046
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.4142342
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}}})   # time: 0.3938517
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(TightBinding),LatticeValue{Float64, :honeycomb}})   # time: 0.3833641
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),LatticeValue{Float64, :square}})   # time: 0.2999187
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.2829747
    Base.precompile(Tuple{typeof(site_coords),SquareLattice{2},LatticeSite{2}})   # time: 0.2797755
    Base.precompile(Tuple{typeof(Haldane),HoneycombLattice,Int64,Int64,Int64})   # time: 0.2401493
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, HoneycombLattice}})   # time: 0.2244773
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{1}},Nothing,Hopping{Matrix{Int64}},NoField})   # time: 0.180226
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Symbol})   # time: 0.1505805
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.146619
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.1446031
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1297785
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeSite{2}})   # time: 0.1285647
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1269412
    Base.precompile(Tuple{typeof(pade_exp),Matrix{ComplexF64},Int64})   # time: 0.1253843
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.1242001
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1239449
    Base.precompile(Tuple{typeof(coord_operators),Basis{HoneycombLattice}})   # time: 0.1202688
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.1142478
    Base.precompile(Tuple{typeof(ldos),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Int64,Float64})   # time: 0.1089329
    Base.precompile(Tuple{Type{LatticeOperator},Basis{SquareLattice{2}},UniformScaling{Bool}})   # time: 0.1022987
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.1008829
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0990202
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:pbc,), Tuple{Bool}},typeof(TightBinding),SquareLattice{1}})   # time: 0.0973984
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.0970121
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0925084
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{Float64}},LandauField})   # time: 0.0891315
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Float64}})   # time: 0.0839465
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{ComplexF64}},LandauField})   # time: 0.0838854
    Base.precompile(Tuple{typeof(diff),LatticeValueRecord})   # time: 0.0826484
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.0758616
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.0742703
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0721807
    Base.precompile(Tuple{typeof(ldos),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0703236
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0650201
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LandauField})   # time: 0.0618731
    Base.precompile(Tuple{typeof(path_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0615591
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0598587
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0578784
    let fbody = try __lookup_kwbody__(which(map_currents, (typeof(site_distance),DensityCurrents,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (typeof(sum),Bool,typeof(map_currents),typeof(site_distance),DensityCurrents,))
    end
end   # time: 0.0575611
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Any})   # time: 0.0573208
    Base.precompile(Tuple{typeof(integrate),LatticeValueRecord})   # time: 0.0572199
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.0560245
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.053984
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0507485
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0495589
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0492777
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.0491171
    Base.precompile(Tuple{Type{SquareLattice},Int64})   # time: 0.0480395
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0477394
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Nothing,Hopping{Matrix{Complex{Int64}}},NoField})   # time: 0.0474395
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.0471994
    Base.precompile(Tuple{typeof(==),LatticeValueRecord,LatticeValueRecord})   # time: 0.0462345
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0446059
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0441365
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0434325
    Base.precompile(Tuple{LatticeValueRecord,Float64})   # time: 0.0409649
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0389533
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeValue{Bool, :square}})   # time: 0.0383388
    Base.precompile(Tuple{LatticeValueRecord,Float64,Float64})   # time: 0.0374203
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Nothing,Bool})   # time: 0.0356954
    Base.precompile(Tuple{typeof(+),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0344809
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeValue{Float64, :square}})   # time: 0.0336469
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.0336469
    Base.precompile(Tuple{typeof(*),LatticeOperator{Adjoint{Float64, Vector{Float64}}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0332629
    Base.precompile(Tuple{typeof(taylor_exp),Matrix{ComplexF64},Int64})   # time: 0.0330796
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0327116
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0316967
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{SquareLattice{2}, Matrix{ComplexF64}},Float64})   # time: 0.0315468
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0305626
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Float64})   # time: 0.0298375
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0298017
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVWStyle})   # time: 0.0296063
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0295763
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0292743
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.0286487
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0272212
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0261477
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0258092
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2}})   # time: 0.025093
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Float64,Int64,Bool})   # time: 0.0243343
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0241656
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{Matrix{Int64}},LandauField})   # time: 0.0241432
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0236859
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Int64, :square}, Matrix{Int64}}})   # time: 0.0236694
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.023334
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0220593
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0205821
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Nothing,Hopping{Matrix{Int64}},NoField})   # time: 0.0197484
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{Matrix{Float64}}})   # time: 0.0195914
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0181672
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0172899
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{Matrix{Int64}},LandauField})   # time: 0.0171368
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0171332
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite{2}})   # time: 0.0167781
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Nothing,Hopping{Matrix{Int64}},LandauField})   # time: 0.016404
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.01532
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0150566
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0129863
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:pade, :k), Tuple{Bool, Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0129794
    Base.precompile(Tuple{typeof(diag_operator),LatticeValue{Float64, :square},Int64})   # time: 0.0125559
    Base.precompile(Tuple{typeof(SpinTightBinding),LatticeValue{Float64, :square}})   # time: 0.0114197
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.01123
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0109892
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0103345
    Base.precompile(Tuple{typeof(path_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0100598
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.0100186
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0098252
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0094847
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(exp), Tuple{LatticeValue{Float64, :honeycomb}}}}}})   # time: 0.0094365
    Base.precompile(Tuple{typeof(-),MaterializedCurrents,MaterializedCurrents})   # time: 0.0091529
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0090763
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},Int64})   # time: 0.008754
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0086346
    Base.precompile(Tuple{typeof(==),Hopping{Matrix{Int64}},Hopping{Matrix{Int64}}})   # time: 0.0084334
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0083736
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0083508
    Base.precompile(Tuple{typeof(==),LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1},LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1}})   # time: 0.0081176
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0081143
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0080571
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0078703
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0078041
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{_A, Matrix{Int64}} where _A<:(LatticeValue{<:Number})})   # time: 0.007779
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0076032
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}}})   # time: 0.0075466
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0074624
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0074303
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0073478
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, HoneycombLattice},Matrix{Float64},Int64,Int64})   # time: 0.0071847
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0071137
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0068679
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}}})   # time: 0.0067342
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{1},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0066753
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0066336
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},BitVector})   # time: 0.0062876
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.0061334
    Base.precompile(Tuple{typeof(Base.maybeview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.00595
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0058405
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.005607
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0053555
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0052692
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.004951
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0049252
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Float64, :square}})   # time: 0.0049143
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}},Vararg{Hopping{Matrix{Int64}}}})   # time: 0.0048034
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{Float64}},NoField})   # time: 0.0046767
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0045301
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0042732
    Base.precompile(Tuple{typeof(-),LatticeOperator{Matrix{Float64}, SquareLattice{2}},UniformScaling{Bool}})   # time: 0.0042308
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0042161
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{Matrix{Int64}}})   # time: 0.0041227
    Base.precompile(Tuple{typeof(insert!),LatticeValueRecord,Int64,LatticeValue{Float64, :square}})   # time: 0.0039819
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0039797
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0039669
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0039537
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{Matrix{ComplexF64}},NoField})   # time: 0.003843
    Base.precompile(Tuple{typeof(_expr_depends_on),Expr,Symbol})   # time: 0.0037399
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0036777
    isdefined(LatticeModels, Symbol("#112#114")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#112#114")),Array})   # time: 0.0035948
    Base.precompile(Tuple{typeof(diag_operator),Function,SquareLattice{2},Int64})   # time: 0.00342
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},TensorProduct{LatticeValue{Float64, :honeycomb}, Matrix{Int64}}})   # time: 0.0034173
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0032793
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0032527
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.0031423
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0030297
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0030204
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.0029592
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.0029461
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.00291
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0027911
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Matrix{Int64}})   # time: 0.0027263
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}, LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0027065
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},TensorProduct{LatticeValue{Float64, :square}, Matrix{Int64}}})   # time: 0.0026989
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0026172
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0026045
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, :square}})   # time: 0.0025747
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))
    end
end   # time: 0.0025361
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.002506
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0024989
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.0023834
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0023667
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.002339
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0022858
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0022828
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64},Bravais{1, 1}})   # time: 0.0022041
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{Float64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{Float64}, SquareLattice{2}}}})   # time: 0.0021922
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0021055
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0020976
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :axis), Tuple{Tuple{Int64, Int64}, Int64}},typeof(hopping),Int64})   # time: 0.0020754
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0020676
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0020592
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0020583
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0020534
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{Vector{Bool}, CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, CartesianIndex{3}, Int64}})   # time: 0.0020469
    Base.precompile(Tuple{typeof(_wrap_eye),String,Matrix{Bool}})   # time: 0.0020362
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.002018
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0019627
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0019541
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0018786
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0018738
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, SquareLattice{2}},LatticeOperator{Matrix{Float64}, SquareLattice{2}}})   # time: 0.0018705
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Int64}})   # time: 0.0018578
    Base.precompile(Tuple{typeof(iterate),LatticeValueRecord})   # time: 0.0018328
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0017994
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.001797
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.001792
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Float64},SVector{2, Float64}})   # time: 0.0017615
    Base.precompile(Tuple{typeof(pairs),LatticeValueRecord})   # time: 0.0016593
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Vararg{Any}})   # time: 0.0016584
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0016341
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.0016059
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.001521
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.001516
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(exp), Tuple{LatticeValue{Float64, :honeycomb}}}}},Type{Float64}})   # time: 0.0015042
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0015036
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}}})   # time: 0.001496
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.001469
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Vararg{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}}})   # time: 0.001444
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0014376
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0013746
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.001342
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0013143
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},ComplexF64,Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}}})   # time: 0.0013125
    Base.precompile(Tuple{typeof(copy),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0013087
    Base.precompile(Tuple{typeof(getindex),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Int64,Int64})   # time: 0.0012929
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0012922
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0012822
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}},))
    end
end   # time: 0.0012814
    Base.precompile(Tuple{typeof(init_record),LatticeValue{Float64, :square}})   # time: 0.0012808
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0012805
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}},))
    end
end   # time: 0.0012513
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Complex{Int64}, Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Complex{Int64}, Float64},))
    end
end   # time: 0.0012464
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}},))
    end
end   # time: 0.0012454
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Complex{Int64},Tuple{Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Complex{Int64},Tuple{Float64},))
    end
end   # time: 0.0011938
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.001187
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0011867
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.001182
    Base.precompile(Tuple{Type{LatticeValue},HoneycombLattice,Vector{Bool}})   # time: 0.0011613
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0011503
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0011496
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.001138
    Base.precompile(Tuple{typeof(isless),LatticeSite{2},LatticeSite{2}})   # time: 0.0011216
    Base.precompile(Tuple{PairLhsSelector,SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0011176
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0011092
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64,Type{Matrix{ComplexF64}}})   # time: 0.0011087
    Base.precompile(Tuple{typeof(projector),Spectrum{HoneycombLattice, Matrix{ComplexF64}}})   # time: 0.0010979
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Tuple{Int64}})   # time: 0.0010912
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Bool}})   # time: 0.0010864
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0010649
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Tuple{LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Tuple{LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))
    end
end   # time: 0.0010546
end
