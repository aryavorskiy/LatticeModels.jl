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
    Base.precompile(Tuple{typeof(diagonalize),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{ComplexF64, Int64}}})   # time: 0.6480974
    Base.precompile(Tuple{typeof(site_distance),Lattice,LatticeSite,LatticeSite})   # time: 0.6472142
    Base.precompile(Tuple{typeof(hoppings),SquareLattice{2},Bonds{Nothing, 1}})   # time: 0.3846873
    Base.precompile(Tuple{typeof(evolution_operator!),Matrix{ComplexF64},Matrix{ComplexF64},Int64})   # time: 0.3700524
    Base.precompile(Tuple{Type{SquareLattice},Function,Int64,Int64})   # time: 0.3029684
    Base.precompile(Tuple{typeof(^),AdjacencyMatrix{SquareLattice{2}},Int64})   # time: 0.295184
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1}})   # time: 0.2886488
    Base.precompile(Tuple{typeof(line_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.2212095
    Base.precompile(Tuple{typeof(evolution_operator!),SparseMatrixCSC{ComplexF64, Int64},SparseMatrixCSC{ComplexF64, Int64},Float64})   # time: 0.2056302
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2}})   # time: 0.2021825
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.1854825
    Base.precompile(Tuple{typeof(shift_site),PeriodicBoundaryConditions,SquareLattice{2},LatticePointer{2}})   # time: 0.1829159
    Base.precompile(Tuple{typeof(coord_operators),LatticeBasis{HoneycombLattice}})   # time: 0.1827454
    Base.precompile(Tuple{typeof(pop!),SquareLattice{2}})   # time: 0.1790193
    Base.precompile(Tuple{typeof(haldane),HoneycombLattice,Int64,Int64,Int64})   # time: 0.1377649
    Base.precompile(Tuple{typeof(|),AdjacencyMatrix{SquareLattice{2}},AdjacencyMatrix{SquareLattice{2}}})   # time: 0.1369989
    Base.precompile(Tuple{typeof(qwz),SquareLattice{2}})   # time: 0.0983236
    Base.precompile(Tuple{typeof(coord_operator),SquareLattice{2},Symbol})   # time: 0.0940104
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},AbstractCurrents})   # time: 0.0891704
    let fbody = try __lookup_kwbody__(which(qwz, (Nothing,SquareLattice{2},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(qwz),Nothing,SquareLattice{2},))
    end
end   # time: 0.08321
    Base.precompile(Tuple{Type{AdjacencyMatrix},SquareLattice{2},Matrix{Bool}})   # time: 0.0760446
    Base.precompile(Tuple{typeof(collect_coords),SquareLattice{3}})   # time: 0.0744948
    let fbody = try __lookup_kwbody__(which(haldane, (Nothing,HoneycombLattice,Int64,Vararg{Int64},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(haldane),Nothing,HoneycombLattice,Int64,Vararg{Int64},))
    end
end   # time: 0.0703386
    Base.precompile(Tuple{typeof(!),AdjacencyMatrix{SquareLattice{2}}})   # time: 0.0652936
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64,Int64})   # time: 0.0555559
    Base.precompile(Tuple{typeof(qwz),LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0545447
    Base.precompile(Tuple{typeof(one_hot),Int64,Val{1}})   # time: 0.0463551
    Base.precompile(Tuple{typeof(line_integral),LandauField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0462783
    Base.precompile(Tuple{typeof(diag_reduce),typeof(tr),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{Float64, Int64}}})   # time: 0.0440815
    Base.precompile(Tuple{typeof(lattice_density),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{Float64, Int64}}})   # time: 0.0429766
    Base.precompile(Tuple{Type{SiteOffset},Nothing,SVector{1, Int64}})   # time: 0.0421869
    Base.precompile(Tuple{typeof(-),MaterializedCurrents,MaterializedCurrents})   # time: 0.0405042
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},SparseMatrixCSC{ComplexF64, Int64},Bonds{Nothing, 1},LandauField,BoundaryConditions{Tuple{}}})   # time: 0.0390868
    Base.precompile(Tuple{typeof(coord_operators),SquareLattice{2}})   # time: 0.0384461
    Base.precompile(Tuple{typeof(apply_field!),Operator{LatticeBasis{SquareLattice{2}}, LatticeBasis{SquareLattice{2}}, SparseMatrixCSC{ComplexF64, Int64}},LandauField})   # time: 0.037036
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, <:Lattice{:square}}})   # time: 0.0364986
    Base.precompile(Tuple{typeof(iterate),LatticeSite{2}})   # time: 0.0362646
    isdefined(LatticeModels, Symbol("#119#120")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#119#120")),Int64})   # time: 0.0355846
    isdefined(LatticeModels, Symbol("#117#118")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#117#118")),Int64})   # time: 0.0354442
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{PairLhsGraph, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1}})   # time: 0.0346874
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, SquareLattice{2}},LatticeValue{Float64, SquareLattice{2}},LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0344494
    Base.precompile(Tuple{Type{SiteOffset},Nothing,SVector{2, Int64}})   # time: 0.0311749
    Base.precompile(Tuple{typeof(site_index),SquareLattice{2},LatticeSite{2}})   # time: 0.0305971
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0304483
    Base.precompile(Tuple{typeof(diagonalize),Hamiltonian{FilledZones{Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}}}, LatticeBasis{HoneycombLattice}, SparseMatrixCSC{ComplexF64, Int64}}})   # time: 0.0274774
    let fbody = try __lookup_kwbody__(which(qwz, (Nothing,SquareLattice{2},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, LandauField, Tuple{Symbol}, NamedTuple{(:field,), Tuple{LandauField}}},typeof(qwz),Nothing,SquareLattice{2},))
    end
end   # time: 0.0262093
    Base.precompile(Tuple{typeof(==),Bonds{Nothing, 1},Bonds{Nothing, 1}})   # time: 0.0260114
    Base.precompile(Tuple{typeof(getindex),TimeSequence{LatticeValue{Float64, SquareLattice{2}}},ClosedInterval{Float64}})   # time: 0.0258955
    Base.precompile(Tuple{typeof(_evolution_block),Expr,Expr})   # time: 0.0258532
    Base.precompile(Tuple{Type{BoundaryConditions},Tuple{TwistedBoundary}})   # time: 0.0250294
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},Int64})   # time: 0.0241701
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, SquareLattice{2}},Symbol})   # time: 0.0234393
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, _A} where _A<:Lattice,Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0224039
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 2}})   # time: 0.0222521
    Base.precompile(Tuple{typeof(coord_values),SquareLattice{2}})   # time: 0.0193175
    Base.precompile(Tuple{typeof(increment!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64}})   # time: 0.0190963
    Base.precompile(Tuple{typeof(hoppings),SquareLattice{2},Bonds{Nothing, 2}})   # time: 0.0190509
    Base.precompile(Tuple{typeof(adjacency_matrix),SquareLattice{2},Bonds{Nothing, 2}})   # time: 0.0181176
    Base.precompile(Tuple{typeof(==),LatticePointer{2},LatticePointer{2}})   # time: 0.0175775
    Base.precompile(Tuple{typeof(*),Int64,MaterializedCurrents})   # time: 0.0175421
    let fbody = try __lookup_kwbody__(which(qwz, (Nothing,SquareLattice{2},LatticeValue{Float64, SquareLattice{2}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(qwz),Nothing,SquareLattice{2},LatticeValue{Float64, SquareLattice{2}},))
    end
end   # time: 0.0174606
    Base.precompile(Tuple{Type{BoundaryConditions},Pair{Int64, Bool}})   # time: 0.0165401
    Base.precompile(Tuple{typeof(hoppings),SquareLattice{2},Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.016344
    Base.precompile(Tuple{typeof(collect_coords),Lattice{:plot_fallback, 2, 1}})   # time: 0.0150389
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},SparseMatrixCSC{Float64, Int64},Bonds{Nothing, 2},LandauField,BoundaryConditions{Tuple{}}})   # time: 0.0145549
    Base.precompile(Tuple{typeof(basis),Sample})   # time: 0.0137924
    Base.precompile(Tuple{typeof(line_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0126127
    Base.precompile(Tuple{Type{HoneycombLattice},Int64,Int64})   # time: 0.0125443
    Base.precompile(Tuple{typeof(getproperty),LatticeSite{2},Symbol})   # time: 0.0110878
    Base.precompile(Tuple{typeof(mm_assuming_zeros),SMatrix{2, 2, Float64, 4},SVector{1, Int64}})   # time: 0.0108318
    Base.precompile(Tuple{typeof(getindex),TimeSequence{LatticeValue{Float64, SquareLattice{2}}},Float64})   # time: 0.0106733
    Base.precompile(Tuple{typeof(setindex!),TimeSequence{LatticeValue},LatticeValue{Float64, SquareLattice{2}},Int64})   # time: 0.0104865
    Base.precompile(Tuple{typeof(differentiate),TimeSequence{LatticeValue}})   # time: 0.0101799
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0100369
    Base.precompile(Tuple{typeof(getproperty),LatticeSite{3},Symbol})   # time: 0.0099567
    Base.precompile(Tuple{typeof(evolution_operator!),SparseMatrixCSC{ComplexF64, Int64},SparseMatrixCSC{ComplexF64, Int64},Int64})   # time: 0.0098654
    Base.precompile(Tuple{typeof(adjacency_matrix),SquareLattice{2},Bonds{Nothing, 1},Vararg{SiteOffset}})   # time: 0.0098116
    Base.precompile(Tuple{typeof(build_hamiltonian),Sample{Nothing, SquareLattice{2}, Nothing, BoundaryConditions{Tuple{}}},Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0094902
    Base.precompile(Tuple{Type{MaterializedCurrents},SquareLattice{2}})   # time: 0.009332
    Base.precompile(Tuple{typeof(cartesian_index),LatticePointer{2}})   # time: 0.0091336
    Base.precompile(Tuple{typeof(hoppings),PairLhsGraph,SquareLattice{2},Bonds{Nothing, 1}})   # time: 0.0086111
    isdefined(LatticeModels, Symbol("#104#105")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#104#105")),Int64})   # time: 0.0083957
    Base.precompile(Tuple{typeof(line_integral),FieldSum{Tuple{FluxField, SymmetricField}},SVector{2, Int64},SVector{2, Int64}})   # time: 0.0082945
    Base.precompile(Tuple{typeof(hoppings),Function,SquareLattice{2},Bonds{Nothing, 1}})   # time: 0.0081856
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticeSite{2}})   # time: 0.0080277
    Base.precompile(Tuple{typeof(add_diagonal!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64},Vector{Float64}})   # time: 0.0080007
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, HoneycombLattice},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, HoneycombLattice}, LatticeValue{Float64, HoneycombLattice}}}})   # time: 0.0079453
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SquareLattice{2}, SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{LatticeValue{Float64, SquareLattice{2}}, Int64}}})   # time: 0.0079097
    Base.precompile(Tuple{typeof(add_diagonal!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{Int64, Int64},Vector{Float64}})   # time: 0.0078389
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, SquareLattice{2}}}}})   # time: 0.0076357
    Base.precompile(Tuple{typeof(coord_operator),SquareLattice{2},Int64})   # time: 0.0076221
    Base.precompile(Tuple{typeof(==),TimeSequence{LatticeValue{Float64, SquareLattice{2}}},TimeSequence{LatticeValue{Float64, SquareLattice{2}}}})   # time: 0.0076047
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, SquareLattice{2}},LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0072366
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, SquareLattice{2}}, LatticeValue{Float64, SquareLattice{2}}}}})   # time: 0.0070367
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, SquareLattice{2}}, LatticeValue{Float64, SquareLattice{2}}}}})   # time: 0.0066987
    Base.precompile(Tuple{Type{SiteOffset},Vector{Int64}})   # time: 0.0066776
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, HoneycombLattice},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, HoneycombLattice}, Base.RefValue{typeof(real)}}}})   # time: 0.0066735
    Base.precompile(Tuple{typeof(differentiate),TimeSequence{Float64}})   # time: 0.0065741
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0064549
    Base.precompile(Tuple{typeof(shift_site),BoundaryConditions{Tuple{TwistedBoundary}},SquareLattice{2},LatticePointer{2}})   # time: 0.0064271
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{typeof(real)}}}})   # time: 0.0063886
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, SquareLattice{2}}, LatticeValue{Float64, SquareLattice{2}}}}})   # time: 0.0062885
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, SquareLattice{2}}, Int64}}})   # time: 0.006188
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, SquareLattice{2}},Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, SquareLattice{2}}}}})   # time: 0.0060837
    Base.precompile(Tuple{typeof(==),LatticeSite{2},LatticePointer{2}})   # time: 0.0057562
    Base.precompile(Tuple{typeof(add_assuming_zeros),SVector{2, Int64},SVector{1, Int64}})   # time: 0.0055314
    Base.precompile(Tuple{typeof(get_site_periodic),SquareLattice{2},LatticePointer{2}})   # time: 0.0054188
    Base.precompile(Tuple{typeof(getproperty),LatticeSite,Symbol})   # time: 0.0053374
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Nothing,SquareLattice{2},SparseMatrixCSC{ComplexF64, Int64},Bonds{Nothing, 2},NoField,BoundaryConditions{Tuple{}}})   # time: 0.005279
    Base.precompile(Tuple{typeof(_extract_lattice),Int64,Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LatticeStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}})   # time: 0.0052385
    Base.precompile(Tuple{typeof(mm_assuming_zeros),SMatrix{2, 2, Float64, 4},SVector{2, Int64}})   # time: 0.0051898
    Base.precompile(Tuple{typeof(==),Bonds{Nothing, 1},Bonds{Nothing, 2}})   # time: 0.0051664
    Base.precompile(Tuple{typeof(Base.Broadcast.dotview),LatticeValue{Float64, SquareLattice{2}},LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0050896
    Base.precompile(Tuple{typeof(coord_value),SquareLattice{2},Symbol})   # time: 0.004931
    Base.precompile(Tuple{typeof(site_distance),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0045203
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{Float64, Int64},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{Float64, Int64},Int64,Int64,))
    end
end   # time: 0.0043952
    Base.precompile(Tuple{typeof(differentiate!),TimeSequence{LatticeValue{Float64, SquareLattice{2}}}})   # time: 0.0043902
    Base.precompile(Tuple{typeof(add_hoppings!),SparseMatrixBuilder{ComplexF64},Function,SquareLattice{2},Int64,Pair{LatticeSite{2}, LatticeSite{2}},NoField,BoundaryConditions{Tuple{}}})   # time: 0.00433
    Base.precompile(Tuple{typeof(integrate),TimeSequence{LatticeValue}})   # time: 0.0042309
    Base.precompile(Tuple{typeof(parse_timesequences!),Expr,Symbol})   # time: 0.0042196
    Base.precompile(Tuple{typeof(==),LatticeValue{T, LT} where {T, LT<:Lattice},LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0037912
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},SparseMatrixCSC{ComplexF64, Int64},Int64,Int64,))
    end
end   # time: 0.0037334
    Base.precompile(Tuple{typeof(_process_increment_recursive!),Expr})   # time: 0.0036023
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, HoneycombLattice}})   # time: 0.0035827
    Base.precompile(Tuple{typeof(integrate!),TimeSequence{LatticeValue{Float64, SquareLattice{2}}}})   # time: 0.0035476
    Base.precompile(Tuple{typeof(getindex),LatticeValueWrapper{SquareLattice{2}},Int64})   # time: 0.0035063
    Base.precompile(Tuple{typeof(_expr_depends_on),Expr,Symbol})   # time: 0.0034813
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Int64},Tuple{Int64, Any}})   # time: 0.0032869
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0030956
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 2}})   # time: 0.0028248
    Base.precompile(Tuple{typeof(line_integral),FluxField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0027725
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, HoneycombLattice}, LatticeValue{Float64, HoneycombLattice}}},Type{Bool}})   # time: 0.0026748
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},Adjoint{ComplexF64, SparseMatrixCSC{ComplexF64, Int64}},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},Adjoint{ComplexF64, SparseMatrixCSC{ComplexF64, Int64}},Int64,Int64,))
    end
end   # time: 0.0026153
    Base.precompile(Tuple{Type{Base.Broadcast.BroadcastStyle},Base.Broadcast.DefaultArrayStyle{1},LatticeStyle})   # time: 0.0025926
    Base.precompile(Tuple{typeof(preprocess_argument),Sample,LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0024939
    let fbody = try __lookup_kwbody__(which(increment!, (SparseMatrixBuilder{ComplexF64},Adjoint{Float64, SparseMatrixCSC{Float64, Int64}},Int64,Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (ComplexF64,typeof(increment!),SparseMatrixBuilder{ComplexF64},Adjoint{Float64, SparseMatrixCSC{Float64, Int64}},Int64,Int64,))
    end
end   # time: 0.0024822
    let fbody = try __lookup_kwbody__(which(coord_operator, (Nothing,SquareLattice{2},Int64,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(coord_operator),Nothing,SquareLattice{2},Int64,))
    end
end   # time: 0.0024086
    Base.precompile(Tuple{typeof(preprocess_argument),Sample,SparseMatrixCSC{ComplexF64, Int64}})   # time: 0.0023461
    Base.precompile(Tuple{typeof(get_inner),LatticeValue{Float64, SquareLattice{2}},LatticeValue{Bool, SquareLattice{2}}})   # time: 0.0022355
    let fbody = try __lookup_kwbody__(which(coord_operators, (Nothing,SquareLattice{2},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(coord_operators),Nothing,SquareLattice{2},))
    end
end   # time: 0.0021877
    Base.precompile(Tuple{Type{SiteOffset},Pair{Int64, Int64},Vector{Int64}})   # time: 0.0021677
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64, Int64},Bravais{3, 1}})   # time: 0.0021363
    Base.precompile(Tuple{typeof(lattice_density),Operator{LatticeBasis{HoneycombLattice}, LatticeBasis{HoneycombLattice}, Matrix{ComplexF64}}})   # time: 0.0021046
    Base.precompile(Tuple{Type{Lattice},Symbol,Tuple{Int64, Int64},Bravais{2, 1}})   # time: 0.0020588
    Base.precompile(Tuple{Type{SquareLattice},Int64,Int64})   # time: 0.002001
    Base.precompile(Tuple{typeof(match),AdjacencyMatrix{SquareLattice{2}},LatticeSite{2},LatticeSite{2}})   # time: 0.0019427
    Base.precompile(Tuple{typeof(dot_assuming_zeros),SVector{2, Float64},Tuple{Int64, Any}})   # time: 0.0019307
    Base.precompile(Tuple{typeof(build_hamiltonian),FilledZones,Vararg{Any}})   # time: 0.001928
    Base.precompile(Tuple{typeof(setindex!),LatticeValueWrapper{SquareLattice{2}},Any,Int64})   # time: 0.0019259
    Base.precompile(Tuple{typeof(plot_fallback),LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0018477
    Base.precompile(Tuple{typeof(line_integral),SymmetricField,SVector{2, Int64},SVector{2, Int64},Int64})   # time: 0.0018468
    let fbody = try __lookup_kwbody__(which(coord_operator, (Nothing,SquareLattice{2},Symbol,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(coord_operator),Nothing,SquareLattice{2},Symbol,))
    end
end   # time: 0.0018299
    Base.precompile(Tuple{typeof(+),LatticePointer{2},Bonds{Pair{Int64, Int64}, 0}})   # time: 0.0018145
    Base.precompile(Tuple{typeof(setindex!),TimeSequence{LatticeValue{Float64, SquareLattice{2}}},LatticeValue{Float64, SquareLattice{2}},Int64})   # time: 0.0018098
    Base.precompile(Tuple{typeof(eltype),LatticeValueWrapper{SquareLattice{2}}})   # time: 0.0018071
    Base.precompile(Tuple{typeof(add_assuming_zeros),SVector{2, Int64},SVector{2, Int64}})   # time: 0.001803
    Base.precompile(Tuple{Type{TimeSequence{LatticeValue}}})   # time: 0.0017146
    isdefined(LatticeModels, Symbol("#130#131")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#130#131")),Float64})   # time: 0.0016448
    Base.precompile(Tuple{Type{TimeSequence},Int64,LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0016245
    Base.precompile(Tuple{typeof(macro_cell_values),LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0016069
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, SquareLattice{2}}, LatticeValue{Float64, SquareLattice{2}}}},Type{Float64}})   # time: 0.0015952
    Base.precompile(Tuple{typeof(iterate),SquareLattice{2},Tuple{CartesianIndices{3, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}}, Int64}})   # time: 0.0015297
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.0015174
    let fbody = try __lookup_kwbody__(which(Sample, (Function,SquareLattice{2},Nothing,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (BoundaryConditions{Tuple{}},Type{Sample},Function,SquareLattice{2},Nothing,))
    end
end   # time: 0.0014981
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0014877
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, SquareLattice{2}}, LatticeValue{Float64, SquareLattice{2}}}},Type{Bool}})   # time: 0.0014269
    Base.precompile(Tuple{typeof(==),LatticeValue{Float64, SquareLattice{2}},LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0013976
    Base.precompile(Tuple{typeof(setindex!),LatticeValueWrapper{SquareLattice{2}},Float64,Int64})   # time: 0.0013184
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, HoneycombLattice}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0013037
    Base.precompile(Tuple{Type{BoundaryConditions},Tuple{}})   # time: 0.0012392
    Base.precompile(Tuple{typeof(basis),Sample{AdjMatT, LT, Nothing} where {AdjMatT, LT}})   # time: 0.0012066
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SquareLattice{2}, SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0011794
    Base.precompile(Tuple{Type{FluxField},Float64})   # time: 0.0011527
    Base.precompile(Tuple{typeof(iterate),TimeSequence{LatticeValue}})   # time: 0.001098
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(|>), Tuple{LatticeValue{Float64, SquareLattice{2}}, Base.RefValue{typeof(real)}}},Type{Float64}})   # time: 0.0010964
    Base.precompile(Tuple{typeof(preprocess_argument),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Complex{Int64}, Bonds{Pair{Int64, Int64}, 1}}})   # time: 0.0010855
    Base.precompile(Tuple{typeof(preprocess_argument),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Complex{Int64}, Bonds{Pair{Int64, Int64}, 2}}})   # time: 0.0010585
    Base.precompile(Tuple{typeof(zero),LatticeValue{Float64, SquareLattice{2}}})   # time: 0.0010487
    Base.precompile(Tuple{typeof(rand),SquareLattice{2}})   # time: 0.0010376
    Base.precompile(Tuple{typeof(preprocess_argument),Sample{Nothing, HoneycombLattice, Nothing, BoundaryConditions{Tuple{}}},Pair{Int64, Bonds{Pair{Int64, Int64}, 0}}})   # time: 0.0010304
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LatticeStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{Int64, LatticeValue{Float64, SquareLattice{2}}}},Type{Float64}})   # time: 0.0010087
    Base.precompile(Tuple{typeof(iterate),LatticeSite{2},Int64})   # time: 0.0010019
end
