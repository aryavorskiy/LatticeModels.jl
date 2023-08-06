const __bodyfunction__ = Dict{Method,Any}()

# Find keyword "body functions" (the function that contains the body
# as written by the developer, called after all missing keyword-arguments
# have been assigned values), in a manner that doesn't depend on
# gensymmed names.
# `mnokw` is the method that gets called when you invoke it without
# supplying any keywords.
function __lookup_kwbody__(mnokw::Method)
    function getsym(ast, arg)
        isa(arg, Symbol) && return arg
        isa(arg, GlobalRef) && return arg.name
        if isa(arg, Core.SSAValue)
            arg = ast.code[arg.id]
            return getsym(ast, arg)
        end
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
                        f = getfield(mnokw.module, getsym(ast, callexpr.args[2]))
                    elseif fsym.mod === Core && fsym.name === :_apply_iterate
                        f = getfield(mnokw.module, getsym(ast, callexpr.args[3]))
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
    Base.precompile(Tuple{typeof(site_distance),Lattice,LatticeSite,LatticeSite})   # time: 0.6202028
    Base.precompile(Tuple{typeof(diagonalize),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{ComplexF64, Int64}}})   # time: 0.4914219
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:statistics,), Tuple{ParticleStatistics}},typeof(densitymatrix),Eigensystem{LatticeBasis{SquareLattice{2}}, Matrix{ComplexF64}}})   # time: 0.4558818
    Base.precompile(Tuple{typeof(hoppings),SquareLattice{2},Bonds{Nothing, 1}})   # time: 0.3522327
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.282558
    Base.precompile(Tuple{typeof(^),AdjacencyMatrix{SquareLattice{2}},Int64})   # time: 0.2786342
    Base.precompile(Tuple{typeof(haldane),HoneycombLattice,Int64,Int64,Int64})   # time: 0.2783684
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.2716752
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1}})   # time: 0.2513833
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.2165287
    Base.precompile(Tuple{typeof(line_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.2115584
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1713283
    Base.precompile(Tuple{typeof(shift_site),PeriodicBoundaryConditions,SquareLattice{2},LatticePointer{2}})   # time: 0.1625207
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1516442
    Base.precompile(Tuple{typeof(evolution_operator!),SparseMatrixCSC{ComplexF64, Int64},SparseMatrixCSC{ComplexF64, Int64},Float64})   # time: 0.1135576
    Base.precompile(Tuple{typeof(coord_operator),SquareLattice{2},Symbol})   # time: 0.0839081
    Base.precompile(Tuple{typeof(coord_operators),LatticeBasis{HoneycombLattice}})   # time: 0.146227
    let fbody = try __lookup_kwbody__(which(haldane, (Nothing,HoneycombLattice,Int64,Vararg{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(haldane),Nothing,HoneycombLattice,Int64,Vararg{Int64},))
    end
end   # time: 0.1038413
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0969747
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:statistics,), Tuple{ParticleStatistics}},typeof(densitymatrix),HamiltonianEigensystem{LatticeBasis{HoneycombLattice}, Matrix{ComplexF64}, FilledZones{Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}}}}})   # time: 0.0921558
    Base.precompile(Tuple{typeof(|),AdjacencyMatrix{SquareLattice{2}},AdjacencyMatrix{SquareLattice{2}}})   # time: 0.086397
    Base.precompile(Tuple{Type{AdjacencyMatrix},SquareLattice{2},Matrix{Bool}})   # time: 0.0820857
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:j1, :j2, :index), Tuple{Int64, Int64, Int64}},typeof(getindex),HoneycombLattice})   # time: 0.0763927
    Base.precompile(Tuple{Type{AdjacencyMatrix},SquareLattice{2},BitMatrix})   # time: 0.075774
    Base.precompile(Tuple{typeof(!),AdjacencyMatrix{SquareLattice{2}}})   # time: 0.0707524
    Base.precompile(Tuple{typeof(qwz),SquareLattice{2}})   # time: 0.0667708
    Base.precompile(Tuple{typeof(adjacency_matrix),SquareLattice{2},Bonds{Nothing, 2}})   # time: 0.065828
    let fbody = try __lookup_kwbody__(which(qwz, (Nothing,SquareLattice{2},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(qwz),Nothing,SquareLattice{2},))
    end
end   # time: 0.0594488
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{NoField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Int64, Bonds{Pair{Int64, Int64}, 0}}})   # time: 0.0556176
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.0536449
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0523084
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{LandauField}},typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0495927
    Base.precompile(Tuple{typeof(line_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0455025
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.045227
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0443286
    Base.precompile(Tuple{typeof(==),Bonds{Nothing, 1},Bonds{Nothing, 1}})   # time: 0.0441518
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:x, :x2), Tuple{Int64, Int64}},typeof(getindex),LatticeValue{Float64, :square}})   # time: 0.0430959
    Base.precompile(Tuple{typeof(qwz),LatticeValue{Float64, :square}})   # time: 0.0415561
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0415535
    Base.precompile(Tuple{typeof(hoppings),SquareLattice{2},Bonds{Nothing, 2}})   # time: 0.0415334
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{Float64, Int64}}})   # time: 0.0408384
    Base.precompile(Tuple{typeof(apply_field!),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{ComplexF64, Int64}},LandauField})   # time: 0.0387842
    Base.precompile(Tuple{typeof(one_hot),Int64,Val{1}})   # time: 0.0380971
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{NoField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Int64, Bonds{Pair{Int64, Int64}, 2}}})   # time: 0.0352311
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0347003
    isdefined(LatticeModels, Symbol("#110#111")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#110#111")),Int64})   # time: 0.0342264
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Operator{LatticeBasis{HoneycombLattice}, LatticeBasis{HoneycombLattice}, SparseMatrixCSC{ComplexF64, Int64}},Pair{Int64, Bonds{Pair{Int64, Int64}, 0}},Pair{Int64, Bonds{Pair{Int64, Int64}, 1}},Vararg{Any}})   # time: 0.0339475
    isdefined(LatticeModels, Symbol("#108#109")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#108#109")),Int64})   # time: 0.0336241
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{PairLhsGraph, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1}})   # time: 0.0332918
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},SparseMatrixCSC{ComplexF64, Int64},Bonds{Nothing, 1},LandauField,BoundaryConditions{Tuple{}}})   # time: 0.0325176
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:axis,), Tuple{Int64}},Type{SiteOffset}})   # time: 0.0321993
    Base.precompile(Tuple{TimeSequence{LatticeValue{Float64, :square}},Float64,Float64})   # time: 0.0271582
    Base.precompile(Tuple{Type{SiteOffset},Nothing,SVector{2, Int64}})   # time: 0.0263974
    Base.precompile(Tuple{typeof(diagonalize),Hamiltonian{FilledZones{Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}}}, LatticeBasis{HoneycombLattice}, SparseMatrixCSC{ComplexF64, Int64}}})   # time: 0.0263201
    Base.precompile(Tuple{Type{SiteOffset},Nothing,SVector{1, Int64}})   # time: 0.0259551
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0258823
    Base.precompile(Tuple{typeof(iterate),Lattice{:square},Tuple{Any, Int64}})   # time: 0.0248116
    let fbody = try __lookup_kwbody__(which(qwz, (Nothing,SquareLattice{2},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Pairs{Symbol, LandauField, Tuple{Symbol}, NamedTuple{(:field,), Tuple{LandauField}}},typeof(qwz),Nothing,SquareLattice{2},))
    end
end   # time: 0.0247594
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{NoField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Complex{Int64}, Bonds{Pair{Int64, Int64}, 1}}})   # time: 0.0223335
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0221769
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0210203
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0204903
    Base.precompile(Tuple{typeof(site_density),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{Float64, Int64}}})   # time: 0.0190599
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 2}})   # time: 0.0188378
    Base.precompile(Tuple{typeof(-),MaterializedCurrents,MaterializedCurrents})   # time: 0.0184858
    Base.precompile(Tuple{typeof(coord_operators),SquareLattice{2}})   # time: 0.0173641
    Base.precompile(Tuple{Type{BoundaryConditions},Pair{Int64, Bool}})   # time: 0.0169853
    Base.precompile(Tuple{typeof(*),Int64,MaterializedCurrents})   # time: 0.0167367
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{LandauField}},typeof(qwz),SquareLattice{2}})   # time: 0.0166815
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{LandauField}},typeof(hoppings),SquareLattice{2},Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0165251
    let fbody = try __lookup_kwbody__(which(qwz, (Nothing,SquareLattice{2},LatticeValue{Float64, :square},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(qwz),Nothing,SquareLattice{2},LatticeValue{Float64, :square},))
    end
end   # time: 0.0162602
    Base.precompile(Tuple{typeof(getindex),TimeSequence{LatticeValue},LatticeValue{Bool, :square}})   # time: 0.0155244
    Base.precompile(Tuple{typeof(==),LatticePointer{2},LatticePointer{2}})   # time: 0.0153249
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.014736
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},Any,SiteOffset,LandauField,BoundaryConditions{Tuple{}}})   # time: 0.0144237
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},SparseMatrixCSC{Float64, Int64},Bonds{Nothing, 2},LandauField,BoundaryConditions{Tuple{}}})   # time: 0.0142497
    Base.precompile(Tuple{typeof(iterate),LatticeSite{2}})   # time: 0.0133907
    Base.precompile(Tuple{Type{TimeSequence{LatticeValue}}})   # time: 0.0126296
    Base.precompile(Tuple{typeof(_kws_to_mask),SquareLattice{2},Any})   # time: 0.0126052
    Base.precompile(Tuple{typeof(site_index),SquareLattice{2},LatticeSite{2}})   # time: 0.0123649
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:pade, :k), Tuple{Bool, Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0121062
    Base.precompile(Tuple{typeof(line_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0119509
    Base.precompile(Tuple{typeof(which(add_assuming_zeros,(Union{Tuple{Vararg{Any, M}}, SArray{Tuple{M}, T, 1, M} where T} where M,Union{Tuple{Vararg{Any, N}}, SArray{Tuple{N}, T, 1, N} where T} where N,)).generator.gen),Any,Any,Any,Any,Any})   # time: 0.0116406
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{NoField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Int64, Bonds{Pair{Int64, Int64}, 1}}})   # time: 0.0114963
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0114819
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractGraph})   # time: 0.0112451
    Base.precompile(Tuple{TimeSequence{LatticeValue{Float64, :square}},Float64})   # time: 0.0107779
    Base.precompile(Tuple{typeof(get_coord),LatticeSite{2},Tuple{Symbol, Nothing}})   # time: 0.009724
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.0096676
    Base.precompile(Tuple{typeof(_extract_lattice),SquareLattice{2},Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}})   # time: 0.0096562
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite})   # time: 0.0093994
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0091569
    Base.precompile(Tuple{typeof(differentiate),TimeSequence{LatticeValue}})   # time: 0.0090699
    Base.precompile(Tuple{typeof(qwz),SquareLattice{2},LatticeValue{Float64, :square}})   # time: 0.0090503
    Base.precompile(Tuple{typeof(add_diagonal!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{Int64, Int64},Vector{Float64}})   # time: 0.0089359
    Base.precompile(Tuple{typeof(hoppings),SquareLattice{2},Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0087628
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite{2}})   # time: 0.008572
    isdefined(LatticeModels, Symbol("#95#96")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#95#96")),Int64})   # time: 0.008525
    Base.precompile(Tuple{typeof(mm_assuming_zeros),SMatrix{2, 2, Float64, 4},SVector{1, Int64}})   # time: 0.0083718
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :honeycomb},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}}})   # time: 0.0082964
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0081744
    Base.precompile(Tuple{typeof(getindex),TimeSequence{LatticeValue},LatticeSite{2}})   # time: 0.0081623
    Base.precompile(Tuple{typeof(line_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0080935
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0080914
    Base.precompile(Tuple{typeof(add_diagonal!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64},Vector{Float64}})   # time: 0.0079846
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}}})   # time: 0.0075587
    Base.precompile(Tuple{typeof(hoppings),Function,SquareLattice{2},Bonds{Nothing, 1}})   # time: 0.0074993
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :honeycomb},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}}})   # time: 0.0073543
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{LandauField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample,Pair{SparseMatrixCSC{ComplexF64, Int64}, LatticeValue{Int64, :square}}})   # time: 0.0072326
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{NoField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Complex{Int64}, Bonds{Pair{Int64, Int64}, 2}}})   # time: 0.0071359
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0071326
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{LandauField}},typeof(build_operator!),SparseMatrixBuilder{ComplexF64},Sample,Pair{<:Any, <:Union{Pair{LatticeSite{N}, LatticeSite{N}} where N, SiteOffset}}})   # time: 0.0070027
    Base.precompile(Tuple{typeof(==),TimeSequence{LatticeValue{Float64, :square}},TimeSequence{LatticeValue{Float64, :square}}})   # time: 0.0069814
    Base.precompile(Tuple{typeof(cartesian_index),LatticePointer{2}})   # time: 0.0068769
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0068625
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}}})   # time: 0.0068306
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0067118
    Base.precompile(Tuple{typeof(differentiate),TimeSequence{Float64}})   # time: 0.0066706
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}}})   # time: 0.0064249
    Base.precompile(Tuple{Type{SiteOffset},Vector{Int64}})   # time: 0.0063515
    Base.precompile(Tuple{typeof(hoppings),PairLhsGraph,SquareLattice{2},Bonds{Nothing, 1}})   # time: 0.0062008
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :honeycomb}, LatticeValue{Float64, :honeycomb}}},Type{Bool}})   # time: 0.0061716
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0059259
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:x1,), Tuple{Int64}},typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square}})   # time: 0.0059145
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, :square}}}})   # time: 0.0058754
    Base.precompile(Tuple{typeof(shift_site),BoundaryConditions{Tuple{TwistedBoundary}},SquareLattice{2},LatticePointer{2}})   # time: 0.0057739
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0057392
    Base.precompile(Tuple{typeof(check_issublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0056351
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticePointer{2}})   # time: 0.0053724
    Base.precompile(Tuple{typeof(get_site_periodic),SquareLattice{2},LatticePointer{2}})   # time: 0.0052379
    Base.precompile(Tuple{typeof(==),LatticeValue,LatticeValue{Float64, :square}})   # time: 0.0051372
    Base.precompile(Tuple{typeof(getproperty),LatticeSite,Symbol})   # time: 0.0051139
    Base.precompile(Tuple{typeof(==),Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0051108
    Base.precompile(Tuple{typeof(check_issublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0047814
    Base.precompile(Tuple{typeof(getproperty),LatticeSite{2},Symbol})   # time: 0.0046708
    Base.precompile(Tuple{typeof(site_distance),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0046205
    Base.precompile(Tuple{typeof(differentiate!),TimeSequence{LatticeValue{Float64, :square}}})   # time: 0.0045409
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},Any,Pair{LatticeSite{N}, LatticeSite{N}} where N,LandauField,BoundaryConditions{Tuple{}}})   # time: 0.0045129
    Base.precompile(Tuple{typeof(_process_increment_recursive!),Expr})   # time: 0.0044009
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0043963
    Base.precompile(Tuple{typeof(mm_assuming_zeros),SMatrix{2, 2, Float64, 4},SVector{2, Int64}})   # time: 0.0041563
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Function,SquareLattice{2},Int64,Pair{LatticeSite{2}, LatticeSite{2}},NoField,BoundaryConditions{Tuple{}}})   # time: 0.0041427
    Base.precompile(Tuple{typeof(add_assuming_zeros),SVector{2, Int64},SVector{1, Int64}})   # time: 0.0041216
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{Float64, Int64},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{Float64, Int64},Int64,Int64,))
    end
end   # time: 0.0038539
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.003845
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64},Int64,Int64,))
    end
end   # time: 0.0036623
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Bool}})   # time: 0.0036276
    Base.precompile(Tuple{typeof(_expr_depends_on),Expr,Symbol})   # time: 0.0034854
    Base.precompile(Tuple{typeof(integrate!),TimeSequence{LatticeValue{Float64, :square}}})   # time: 0.0032508
    Base.precompile(Tuple{typeof(insert!),TimeSequence{LatticeValue},Int64,LatticeValue{Float64, :square}})   # time: 0.0031664
    Base.precompile(Tuple{typeof(integrate),TimeSequence{LatticeValue}})   # time: 0.0031059
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Float64, :square}})   # time: 0.0028922
    Base.precompile(Tuple{typeof(preprocess_argument),Sample,LatticeValue{Float64, :square}})   # time: 0.002881
    Base.precompile(Tuple{typeof(iterate),Lattice{:square}})   # time: 0.0028093
    Base.precompile(Tuple{typeof(adjacency_matrix),SquareLattice{2},Bonds{Nothing, 1},Vararg{SiteOffset}})   # time: 0.0027292
    let fbody = try __lookup_kwbody__(which(build_hamiltonian, (Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Operator{LatticeBasis{HoneycombLattice}, LatticeBasis{HoneycombLattice}, SparseMatrixCSC{ComplexF64, Int64}},Vararg{Any},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(build_hamiltonian),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Operator{LatticeBasis{HoneycombLattice}, LatticeBasis{HoneycombLattice}, SparseMatrixCSC{ComplexF64, Int64}},Vararg{Any},))
    end
end   # time: 0.0026945
    Base.precompile(Tuple{typeof(line_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0026941
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},Adjoint{Float64, SparseMatrixCSC{Float64, Int64}},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},Adjoint{Float64, SparseMatrixCSC{Float64, Int64}},Int64,Int64,))
    end
end   # time: 0.0025415
    Base.precompile(Tuple{typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0025322
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:μ, :statistics), Tuple{Int64, ParticleStatistics}},Type{System},Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}}})   # time: 0.0025304
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Int64},Tuple{Int64, Any}})   # time: 0.0024129
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0023363
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},Adjoint{ComplexF64, SparseMatrixCSC{ComplexF64, Int64}},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},Adjoint{ComplexF64, SparseMatrixCSC{ComplexF64, Int64}},Int64,Int64,))
    end
end   # time: 0.0023198
    Base.precompile(Tuple{typeof(_extract_lattice_s),HoneycombLattice,HoneycombLattice,Tuple{}})   # time: 0.0022625
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LatticeStyle})   # time: 0.0022219
    Base.precompile(Tuple{typeof(_extract_lattice_s),SquareLattice{2},SquareLattice{2},Tuple{}})   # time: 0.0020796
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0020776
    Base.precompile(Tuple{typeof(site_density),Operator{LatticeBasis{HoneycombLattice}, LatticeBasis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0020083
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:field,), Tuple{LandauField}},typeof(build_hamiltonian),FilledZones,Vararg{Any}})   # time: 0.0019026
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:k,), Tuple{Int64}},typeof(_evolution_block),Expr,Expr})   # time: 0.0018879
    Base.precompile(Tuple{typeof(match),AdjacencyMatrix{SquareLattice{2}},LatticeSite{2},LatticeSite{2}})   # time: 0.0018631
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :honeycomb}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0018427
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0018143
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Float64},Tuple{Int64, Any}})   # time: 0.0017843
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0017748
    isdefined(LatticeModels, Symbol("#126#127")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#126#127")),Float64})   # time: 0.0017709
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64})   # time: 0.0017373
    Base.precompile(Tuple{typeof(preprocess_argument),Sample,Pair})   # time: 0.0016903
    Base.precompile(Tuple{typeof(line_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0016519
    Base.precompile(Tuple{Type{TimeSequence},LatticeValue{Float64, :square}})   # time: 0.0016014
    Base.precompile(Tuple{Type{BoundaryConditions},Tuple{}})   # time: 0.0015725
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, Int64}})   # time: 0.0015139
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, :square}})   # time: 0.001513
    Base.precompile(Tuple{typeof(+),LatticePointer{2},Bonds{Pair{Int64, Int64}, 0}})   # time: 0.0014979
    Base.precompile(Tuple{Type{SiteOffset},Pair{Int64, Int64},Vector{Int64}})   # time: 0.0014978
    let fbody = try __lookup_kwbody__(which(Sample, (Function,SquareLattice{2},Nothing,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Type{Sample},Function,SquareLattice{2},Nothing,))
    end
end   # time: 0.0014928
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, :square}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0014496
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0014438
    Base.precompile(Tuple{typeof(_extract_lattice),HoneycombLattice,Tuple{LatticeValue{Float64, :honeycomb}}})   # time: 0.0014372
    Base.precompile(Tuple{typeof(preprocess_argument),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Int64, Bonds{Pair{Int64, Int64}, 0}}})   # time: 0.0014339
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Int64, :square},LatticeSite{2}})   # time: 0.0013308
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, :square},LatticeValue{Float64, :square}})   # time: 0.001252
    Base.precompile(Tuple{Type{LatticeValue},Lattice{:plot_fallback, 2, 1},Vector{Float64}})   # time: 0.0011967
    Base.precompile(Tuple{Type{BoundaryConditions},Tuple{TwistedBoundary}})   # time: 0.0011481
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, :square}})   # time: 0.0010876
    Base.precompile(Tuple{typeof(preprocess_argument),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Int64, Bonds{Pair{Int64, Int64}, 2}}})   # time: 0.0010876
    Base.precompile(Tuple{typeof(add_assuming_zeros),SVector{2, Int64},SVector{2, Int64}})   # time: 0.0010754
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0010434
    Base.precompile(Tuple{typeof(insert!),TimeSequence{LatticeValue{Float64, :square}},Int64,LatticeValue{Float64, :square}})   # time: 0.0010135
    Base.precompile(Tuple{typeof(zeros),SquareLattice{2}})   # time: 0.0010128
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Float64, :square}}},Type{Float64}})   # time: 0.0010086
end
