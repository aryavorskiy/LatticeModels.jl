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
    Base.precompile(Tuple{Core.kwftype(typeof(map_currents)),NamedTuple{(:reduce_fn, :sort), Tuple{typeof(sum), Bool}},typeof(map_currents),Function,DensityCurrents})   # time: 1.0616351
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{1}})   # time: 0.8064909
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Int64})   # time: 0.4704298
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{Basis{SquareLattice{2}}, Matrix{ComplexF64}}})   # time: 0.3798119
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.3791504
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{1}}},Nothing,Hopping{1},NoField})   # time: 0.3112514
    Base.precompile(Tuple{typeof(Haldane),HoneycombLattice,Int64,Int64,Int64})   # time: 0.2971772
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{HoneycombLattice}}})   # time: 0.2300329
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),Function,LatticeValue{Float64, :square}})   # time: 0.226362
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:j1, :j2, :index), Tuple{Int64, Int64, Int64}},typeof(getindex),HoneycombLattice})   # time: 0.2126204
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.185525
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1733509
    Base.precompile(Tuple{typeof(_angle),SVector{2, Int64},SVector{2, Int64}})   # time: 0.1639154
    Base.precompile(Tuple{typeof(ldos),Spectrum{Basis{SquareLattice{2}}, Matrix{ComplexF64}},Int64,Float64})   # time: 0.1606383
    Base.precompile(Tuple{Type{LatticeOperator},Basis{SquareLattice{2}},UniformScaling{Bool}})   # time: 0.148086
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeSite{2}})   # time: 0.1448105
    Base.precompile(Tuple{typeof(coord_operators),Basis{HoneycombLattice}})   # time: 0.1434976
    Base.precompile(Tuple{typeof(pade_exp),Matrix{ComplexF64},Int64})   # time: 0.1434541
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.1390004
    Base.precompile(Tuple{typeof(ptrace),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Symbol})   # time: 0.1355539
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.1278094
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.1264102
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1230059
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:pbc,), Tuple{Bool}},typeof(TightBinding),SquareLattice{1}})   # time: 0.1209705
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.1195617
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.1107136
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.1070198
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},SMatrix{2, 2, Float64, 4},Int64,Int64})   # time: 0.0997829
    Base.precompile(Tuple{typeof(apply_field!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LandauField})   # time: 0.0954588
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.0950082
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.0946857
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.094471
    Base.precompile(Tuple{typeof(diff),LatticeValueRecord})   # time: 0.0940187
    Base.precompile(Tuple{typeof(getindex),DensityCurrents,LatticeValue{Bool, :square}})   # time: 0.0926651
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.0901182
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Nothing,Hopping{2},LandauField})   # time: 0.0880999
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Nothing,Hopping,LandauField})   # time: 0.0828373
    Base.precompile(Tuple{typeof(==),LatticeValueRecord,LatticeValueRecord})   # time: 0.0818084
    Base.precompile(Tuple{typeof(ldos),Spectrum{Basis{SquareLattice{2}}, Matrix{ComplexF64}},Float64})   # time: 0.0802158
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.0724432
    Base.precompile(Tuple{typeof(|),PairSet{SquareLattice{2}},PairSet{SquareLattice{2}}})   # time: 0.0686861
    Base.precompile(Tuple{typeof(integrate),LatticeValueRecord})   # time: 0.0680098
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0668233
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:x, :x2), Tuple{Int64, Int64}},typeof(getindex),LatticeValue{Float64, :square}})   # time: 0.0656036
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},SMatrix{2, 2, Int64, 4},Int64,Int64})   # time: 0.0650302
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Float64}})   # time: 0.0647749
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{2},LandauField})   # time: 0.0631442
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},PairSet})   # time: 0.063047
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Any})   # time: 0.062839
    Base.precompile(Tuple{typeof(*),LatticeOperator{Adjoint{Float64, Vector{Float64}}, Basis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}})   # time: 0.0625743
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeValue{Float64, :square}})   # time: 0.0615102
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Float64,Nothing,Bool})   # time: 0.0606819
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0606556
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{Basis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0603859
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0601829
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.059619
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2}})   # time: 0.0579321
    isdefined(LatticeModels, Symbol("#86#87")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#86#87")),Int64})   # time: 0.0573458
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0569097
    let fbody = try __lookup_kwbody__(which(map_currents, (typeof(site_distance),DensityCurrents,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (typeof(sum),Bool,typeof(map_currents),typeof(site_distance),DensityCurrents,))
    end
end   # time: 0.0566414
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}})   # time: 0.0506145
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0474937
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0468285
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0454402
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0445395
    Base.precompile(Tuple{LatticeValueRecord,Float64,Float64})   # time: 0.040653
    Base.precompile(Tuple{typeof(getindex),Spectrum{Basis{HoneycombLattice}, Matrix{ComplexF64}},Int64})   # time: 0.0404714
    Base.precompile(Tuple{typeof(path_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.039813
    isdefined(LatticeModels, Symbol("#82#83")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#82#83")),Int64})   # time: 0.0397923
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),Function,SquareLattice{2}})   # time: 0.0391828
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeValue{Bool, :square}})   # time: 0.0391058
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},Nothing,Hopping{1},NoField})   # time: 0.038124
    Base.precompile(Tuple{typeof(+),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}})   # time: 0.0378402
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.0374888
    Base.precompile(Tuple{typeof(*),Int64,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0361414
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0353866
    Base.precompile(Tuple{typeof(site_density),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0344904
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:E,), Tuple{Int64}},typeof(getindex),Spectrum{Basis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0343712
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0342565
    Base.precompile(Tuple{typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0340447
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Int64})   # time: 0.0339121
    Base.precompile(Tuple{typeof(filled_projector),Spectrum{Basis{SquareLattice{2}}, Matrix{ComplexF64}},Float64})   # time: 0.0338848
    Base.precompile(Tuple{Type{SquareLattice},Int64})   # time: 0.0337739
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Float64})   # time: 0.033357
    Base.precompile(Tuple{typeof(taylor_exp),Matrix{ComplexF64},Int64})   # time: 0.0319879
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping)})   # time: 0.0308918
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0308044
    Base.precompile(Tuple{typeof(site_index),Lattice,LatticeSite{2}})   # time: 0.0306267
    Base.precompile(Tuple{typeof(site_density),LatticeArray{Vector{ComplexF64}, Basis{HoneycombLattice}, 1}})   # time: 0.0302111
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0295941
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(TightBinding),LatticeValue{Float64, :honeycomb}})   # time: 0.0294462
    Base.precompile(Tuple{typeof(-),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0288124
    Base.precompile(Tuple{typeof(evolution_operator),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Float64,Int64,Bool})   # time: 0.0285874
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LVWStyle})   # time: 0.0262221
    Base.precompile(Tuple{typeof(==),LatticeSite,LatticeSite})   # time: 0.0261692
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}})   # time: 0.0250988
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}})   # time: 0.0246539
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0243793
    Base.precompile(Tuple{typeof(materialize),SparseMatrixBuilder{ComplexF64}})   # time: 0.023859
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0231054
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0206027
    Base.precompile(Tuple{typeof(eltype),LatticeValue{ComplexF64, :square}})   # time: 0.0186639
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{1},LandauField})   # time: 0.018243
    Base.precompile(Tuple{typeof(materialize),DensityCurrents})   # time: 0.0176935
    Base.precompile(Tuple{typeof(==),LatticeOperator{SparseMatrixCSC{ComplexF64, Int64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0172895
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0172703
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{SparseMatrixBuilder{ComplexF64}, Basis{SquareLattice{2}}},Function,Hopping{2},LandauField})   # time: 0.016564
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0160233
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0159108
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{SquareLattice{2}},Type{SparseMatrixBuilder{ComplexF64}}})   # time: 0.0158809
    isdefined(LatticeModels, Symbol("#100#101")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#100#101")),Float64})   # time: 0.0154534
    let fbody = try __lookup_kwbody__(which(SpinTightBinding, (LatticeValue{Float64, :square},Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(SpinTightBinding),LatticeValue{Float64, :square},Vararg{Any},))
    end
end   # time: 0.0138286
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0137431
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.0135053
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Function,Hopping{2},LandauField})   # time: 0.0134879
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:pade, :k), Tuple{Bool, Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0134828
    Base.precompile(Tuple{typeof(_kws_to_mask),SquareLattice{2},Any})   # time: 0.0134701
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},Nothing,Hopping{1},LandauField})   # time: 0.0130574
    Base.precompile(Tuple{typeof(coord_values),HoneycombLattice})   # time: 0.0123858
    Base.precompile(Tuple{LatticeValueRecord,Float64})   # time: 0.0123271
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0122657
    Base.precompile(Tuple{Type{LatticeSite},Vector{Int64},Int64})   # time: 0.0122079
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0114389
    Base.precompile(Tuple{typeof(dot),LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}})   # time: 0.0111756
    Base.precompile(Tuple{typeof(path_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0111224
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(exp), Tuple{LatticeValue{Float64, :honeycomb}}}}}})   # time: 0.0109307
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0108944
    Base.precompile(Tuple{typeof(diag_operator),LatticeValue{Float64, :square},Int64})   # time: 0.0108596
    Base.precompile(Tuple{typeof(adjoint),LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}})   # time: 0.0107788
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0105783
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Float64}})   # time: 0.0102167
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0100894
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{1}})   # time: 0.0099678
    Base.precompile(Tuple{typeof(-),MaterializedCurrents,MaterializedCurrents})   # time: 0.0098105
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0097862
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0094128
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0092169
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0090597
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.009036
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0090151
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.008888
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc,), Tuple{Vector{Int64}}},typeof(hopping)})   # time: 0.0087734
    Base.precompile(Tuple{typeof(==),LatticeArray{Vector{ComplexF64}, Basis{HoneycombLattice}, 1},LatticeArray{Vector{ComplexF64}, Basis{HoneycombLattice}, 1}})   # time: 0.0087522
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0087005
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.00864
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0084168
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0081237
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Vector{Bool}}},typeof(hopping)})   # time: 0.0080401
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.007958
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}}})   # time: 0.0079562
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0079257
    isdefined(LatticeModels, Symbol("#106#108")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#106#108")),Vector{Float64}})   # time: 0.0079098
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}}})   # time: 0.0078828
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{SquareLattice{1}},Type{Matrix{ComplexF64}}})   # time: 0.0078391
    Base.precompile(Tuple{Core.kwftype(typeof(dotview)),NamedTuple{(:x1,), Tuple{Int64}},typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square}})   # time: 0.0077695
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :honeycomb}})   # time: 0.0076783
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0076341
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0076258
    Base.precompile(Tuple{typeof(==),Hopping{1},Hopping{1}})   # time: 0.0074234
    Base.precompile(Tuple{typeof(getindex),Spectrum{Basis{HoneycombLattice}, Matrix{ComplexF64}},BitVector})   # time: 0.0071594
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{HoneycombLattice}},Matrix{Float64},Int64,Int64})   # time: 0.0071274
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0069643
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0069561
    Base.precompile(Tuple{typeof(⊗),LatticeValue{Float64, :square},Matrix{Int64}})   # time: 0.0069149
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}})   # time: 0.0066583
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},DensityCurrents})   # time: 0.006168
    Base.precompile(Tuple{Type{LatticeValueRecord},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.0061063
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0060523
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Float64})   # time: 0.005838
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64},Matrix{Int64}})   # time: 0.0058107
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{1},Vararg{Hopping{1}}})   # time: 0.0056357
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Int64, :square}})   # time: 0.0053501
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{SparseMatrixBuilder{ComplexF64}, Basis{SquareLattice{2}}},Matrix{Complex{Int64}}})   # time: 0.0049631
    Base.precompile(Tuple{typeof(-),LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},UniformScaling{Bool}})   # time: 0.0049409
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Nothing,Hopping{2},NoField})   # time: 0.0048424
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Matrix{Int64},Int64,Int64})   # time: 0.0047514
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64, Int64}}},typeof(hopping)})   # time: 0.0046152
    Base.precompile(Tuple{typeof(==),LatticeValue,LatticeValue{Float64, :square}})   # time: 0.0046059
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{Int64}})   # time: 0.0045973
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{1}})   # time: 0.0045685
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0044898
    Base.precompile(Tuple{Type{LatticeArray},Basis{SquareLattice{2}},Vector{Float64}})   # time: 0.0044031
    Base.precompile(Tuple{typeof(_expr_depends_on),Expr,Symbol})   # time: 0.0042893
    Base.precompile(Tuple{typeof(copy),TensorProduct{LatticeValue{Float64, :square}, 2, Int64}})   # time: 0.0042851
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Matrix{ComplexF64},Int64,Int64})   # time: 0.0042184
    Base.precompile(Tuple{typeof(insert!),LatticeValueRecord,Int64,LatticeValue{Float64, :square}})   # time: 0.0042003
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}})   # time: 0.004118
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0040878
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :pbc), Tuple{Vector{Int64}, Vector{Bool}}},typeof(hopping),Matrix{Int64}})   # time: 0.0040841
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0040663
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{SparseMatrixBuilder{ComplexF64}, Basis{SquareLattice{2}}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.0040015
    Base.precompile(Tuple{typeof(evolved),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0038474
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Bool, :honeycomb},LatticeSite{2}})   # time: 0.0037975
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0037356
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0037067
    Base.precompile(Tuple{typeof(_make_wrapper),SparseMatrixCSC{ComplexF64, Int64},Basis{SquareLattice{2}}})   # time: 0.0035758
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :pbc), Tuple{Tuple{Int64, Int64}, Vector{Bool}}},typeof(hopping)})   # time: 0.0034732
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0034579
    Base.precompile(Tuple{typeof(copy),TensorProduct{LatticeValue{Float64, :honeycomb}, 2, Int64}})   # time: 0.0034095
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0033708
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Matrix{Int64}})   # time: 0.0032035
    isdefined(LatticeModels, Symbol("#28#29")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#28#29")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.003135
    isdefined(LatticeModels, Symbol("#107#109")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#107#109")),Tuple{Float64, Vector{Float64}}})   # time: 0.0030459
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{ComplexF64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0029374
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.002909
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.0028845
    Base.precompile(Tuple{typeof(*),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}})   # time: 0.0028744
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:translate_uc, :site_indices), Tuple{Vector{Int64}, Int64}},typeof(hopping)})   # time: 0.0028339
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, Int64}})   # time: 0.0028274
    isdefined(LatticeModels, Symbol("#106#108")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#106#108")),Array})   # time: 0.0028078
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},TensorProduct{_A, _B, Int64} where {_A<:(LatticeValue{<:Number}), _B}})   # time: 0.0026652
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0026098
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}})   # time: 0.0026006
    Base.precompile(Tuple{typeof(*),ComplexF64,LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0025766
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}, LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}}})   # time: 0.0025683
    Base.precompile(Tuple{typeof(dot),LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1},LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}})   # time: 0.0025594
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0025561
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices,), Tuple{Tuple{Int64, Int64}}},typeof(hopping)})   # time: 0.002546
    Base.precompile(Tuple{typeof(diag_operator),Function,SquareLattice{2},Int64})   # time: 0.0025391
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},TensorProduct{LatticeValue{Float64, :honeycomb}, 2, Int64}})   # time: 0.0025131
    Base.precompile(Tuple{typeof(-),UniformScaling{Bool},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0024941
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},TensorProduct{LatticeValue{Int64, :square}, 2, Int64}})   # time: 0.0024908
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:site_indices, :axis), Tuple{Tuple{Int64, Int64}, Int64}},typeof(hopping),Int64})   # time: 0.002464
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{SquareLattice{2}},Type{Matrix{ComplexF64}}})   # time: 0.0024625
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64},Bravais{1, 1}})   # time: 0.0024504
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},TensorProduct{LatticeValue{Float64, :square}, 2, Int64}})   # time: 0.0024476
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}})   # time: 0.0024254
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Bool}},typeof(hopping)})   # time: 0.0024111
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{1}})   # time: 0.0024041
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0023902
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},Int64})   # time: 0.0023572
    Base.precompile(Tuple{PairSet{SquareLattice{2}},SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0023508
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0023456
    Base.precompile(Tuple{typeof(getindex),Spectrum{Basis{HoneycombLattice}, Matrix{ComplexF64}},UnitRange{Int64}})   # time: 0.0023445
    isdefined(LatticeModels, Symbol("#28#29")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#28#29")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0023326
    Base.precompile(Tuple{typeof(_wrap_eye),String,Matrix{Bool}})   # time: 0.0023315
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0023085
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{1}}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.0022947
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.0022452
    Base.precompile(Tuple{Type{LatticeRecord},Vector{LatticeValue{Float64, :square}},Vector{Int64}})   # time: 0.0022447
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},Vararg{LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}}})   # time: 0.0022339
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},Matrix{Float64},Int64,Int64})   # time: 0.0022229
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},SMatrix{2, 2, Float64, 4},Int64,Int64})   # time: 0.0022136
    Base.precompile(Tuple{typeof(materialize),PairSet{SquareLattice{2}},SubCurrents{DensityCurrents}})   # time: 0.0021706
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0021152
    Base.precompile(Tuple{typeof(site_distance),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0020962
    Base.precompile(Tuple{typeof(materialize),SubCurrents{DensityCurrents}})   # time: 0.002055
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0020521
    Base.precompile(Tuple{typeof(+),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}})   # time: 0.0020384
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0019943
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Float64}})   # time: 0.0019938
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0018977
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Vector{Float64}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Vector{Float64}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}, LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0018945
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Vararg{Any}})   # time: 0.0018935
    Base.precompile(Tuple{typeof(increment!),LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.0018866
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, :square}})   # time: 0.0018281
    Base.precompile(Tuple{Core.kwftype(typeof(_evolution_block)),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0018241
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0018127
    Base.precompile(Tuple{typeof(*),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0018117
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}})   # time: 0.0017898
    Base.precompile(Tuple{typeof(iterate),LatticeValueRecord})   # time: 0.0017524
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Vararg{LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}}})   # time: 0.0017438
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.0017181
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Int64}}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0017054
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},Matrix{Float64},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}}})   # time: 0.0016679
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},ComplexF64,Tuple{LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}}})   # time: 0.0016666
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :honeycomb}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(exp), Tuple{LatticeValue{Float64, :honeycomb}}}}},Type{Float64}})   # time: 0.0016658
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}},))
    end
end   # time: 0.0016489
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Int64, Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Int64, Int64},))
    end
end   # time: 0.0016339
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0015937
    Base.precompile(Tuple{typeof(copy),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.001577
    Base.precompile(Tuple{typeof(copy),SquareLattice{3}})   # time: 0.001566
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0015506
    Base.precompile(Tuple{typeof(_lazy_tp),LatticeValue{Int64, :square},Matrix{Int64}})   # time: 0.0015315
    Base.precompile(Tuple{typeof(_zero_on_basis),Basis{HoneycombLattice}})   # time: 0.0015282
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0015147
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.001506
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector})   # time: 0.001494
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}, LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},Tuple{LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}, LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}},))
    end
end   # time: 0.0014925
    Base.precompile(Tuple{typeof(SpinTightBinding),LatticeValue{Float64, :square},Int64})   # time: 0.0014912
    Base.precompile(Tuple{typeof(copy),SquareLattice{2}})   # time: 0.0014693
    Base.precompile(Tuple{typeof(init_record),LatticeValue{Float64, :square}})   # time: 0.001468
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}},Type{Float64}})   # time: 0.0014658
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{LatticeOperator{SparseMatrixBuilder{ComplexF64}, Basis{SquareLattice{2}}}}})   # time: 0.0014609
    Base.precompile(Tuple{typeof(projector),Spectrum{Basis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0013555
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0013402
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}},))
    end
end   # time: 0.0013295
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0013184
    Base.precompile(Tuple{typeof(path_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0013158
    Base.precompile(Tuple{Type{LatticeValue},HoneycombLattice,Vector{Bool}})   # time: 0.0013037
    Base.precompile(Tuple{typeof(_unwrap),Function,Tuple{},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},Tuple{Int64}})   # time: 0.0012906
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{ComplexF64, LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}},))
    end
end   # time: 0.0012796
    Base.precompile(Tuple{Type{DensityCurrents},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0012712
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{Complex{Int64}, Float64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{Complex{Int64}, Float64},))
    end
end   # time: 0.0012627
    Base.precompile(Tuple{typeof(_unwrap_from_macro),Function,LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{2}}}})   # time: 0.0012585
    let fbody = try __lookup_kwbody__(which(_unwrap, (Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap),Function,Tuple{},Tuple{UniformScaling{Bool}, LatticeOperator{Matrix{ComplexF64}, Basis{HoneycombLattice}}},))
    end
end   # time: 0.0012561
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{1}}},LatticeOperator{Matrix{ComplexF64}, Basis{SquareLattice{1}}}})   # time: 0.0012359
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Bool}})   # time: 0.0012214
    Base.precompile(Tuple{PairLhsSelector,HoneycombLattice,LatticeSite{2},LatticeSite{2}})   # time: 0.001218
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0012029
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Float64},SVector{2, Float64}})   # time: 0.0012028
    Base.precompile(Tuple{typeof(isless),LatticeSite{2},LatticeSite{2}})   # time: 0.0011976
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},Tuple{LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}},Tuple{LatticeOperator{Matrix{Float64}, Basis{SquareLattice{2}}}},))
    end
end   # time: 0.0011972
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0011272
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector{Int64}})   # time: 0.0010836
    Base.precompile(Tuple{typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}}}})   # time: 0.0010727
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}},Type{Float64}})   # time: 0.0010532
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :square}})   # time: 0.0010516
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Adjoint{Float64, Vector{Float64}}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0010282
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Vector{Float64}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Vector{Float64}},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, Basis{SquareLattice{2}}},Tuple{LatticeArray{Vector{Float64}, Basis{SquareLattice{2}}, 1}},))
    end
end   # time: 0.0010228
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0010181
end
