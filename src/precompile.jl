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
    let fbody = try __lookup_kwbody__(which(map_currents, (typeof(site_distance),DensityCurrents,))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (typeof(sum),Bool,typeof(map_currents),typeof(site_distance),DensityCurrents,))
    end
end   # time: 0.9735952
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{1}})   # time: 0.5386716
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:j1, :j2, :index), Tuple{Int64, Int64, Int64}},typeof(getindex),HoneycombLattice})   # time: 0.3183274
    Base.precompile(Tuple{typeof(project),LatticeValue{Float64, :square},Symbol})   # time: 0.2301702
    Base.precompile(Tuple{Core.kwftype(typeof(TightBinding)),NamedTuple{(:pbc,), Tuple{Bool}},typeof(TightBinding),SquareLattice{1}})   # time: 0.1546252
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice})   # time: 0.1395991
    Base.precompile(Tuple{typeof(coord_operators),Basis{HoneycombLattice}})   # time: 0.1307498
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{1}},Nothing,Hopping{1},NoField})   # time: 0.1207519
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{2},LandauField})   # time: 0.0814459
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),Function,LatticeValue{Float64, :square}})   # time: 0.0797958
    isdefined(LatticeModels, Symbol("#102#103")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#102#103")),Int64})   # time: 0.0796542
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0743964
    Base.precompile(Tuple{typeof(coord_operators),Basis{SquareLattice{2}}})   # time: 0.07072
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},Lattice,Any})   # time: 0.069989
    Base.precompile(Tuple{typeof(getindex),SquareLattice{2},LatticeValue{Bool, :square}})   # time: 0.0616314
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Nothing,Hopping{1},NoField})   # time: 0.0607948
    Base.precompile(Tuple{typeof(hopping_operator),Function,HoneycombLattice,Hopping{2},LandauField})   # time: 0.057081
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Nothing,Hopping{1},LandauField})   # time: 0.0562994
    Base.precompile(Tuple{Type{PairSet},SquareLattice{2},BitMatrix})   # time: 0.0551008
    Base.precompile(Tuple{Core.kwftype(typeof(getindex)),NamedTuple{(:x, :x2), Tuple{Int64, Int64}},typeof(getindex),LatticeValue{Float64, :square}})   # time: 0.0547369
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping,LandauField})   # time: 0.0517717
    Base.precompile(Tuple{typeof(bonds),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}})   # time: 0.0501984
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Matrix{Complex{Int64}}})   # time: 0.0501556
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}}})   # time: 0.0452136
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.0399979
    Base.precompile(Tuple{typeof(!),PairSet{SquareLattice{2}}})   # time: 0.0388518
    Base.precompile(Tuple{Core.kwftype(typeof(SpinTightBinding)),NamedTuple{(:field,), Tuple{LandauField}},typeof(SpinTightBinding),Function,SquareLattice{2}})   # time: 0.0379513
    isdefined(LatticeModels, Symbol("#98#99")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#98#99")),Int64})   # time: 0.0378991
    Base.precompile(Tuple{typeof(^),PairSet{SquareLattice{2}},Int64})   # time: 0.034993
    Base.precompile(Tuple{typeof(site_density),LatticeArray{Vector{ComplexF64}, HoneycombLattice, 1}})   # time: 0.0330594
    Base.precompile(Tuple{Type{LatticeSite},Vector{Int64},Int64})   # time: 0.0318464
    Base.precompile(Tuple{typeof(==),LatticeSite,LatticeSite})   # time: 0.0286208
    Base.precompile(Tuple{LatticeValueRecord,Float64,Float64})   # time: 0.0265472
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(*),Any,LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}})   # time: 0.022242
    Base.precompile(Tuple{typeof(getindex),LatticeValueRecord,LatticeValue{Bool, :square}})   # time: 0.0187908
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, :square}})   # time: 0.0186307
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},SMatrix{2, 2, ComplexF64, 4},Int64,Int64})   # time: 0.0170394
    Base.precompile(Tuple{Type{HoneycombLattice},Function,Int64,Int64})   # time: 0.0169604
    Base.precompile(Tuple{typeof(site_coords),HoneycombLattice,LatticeSite{2}})   # time: 0.0169382
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Function,Hopping{2},LandauField})   # time: 0.0164494
    Base.precompile(Tuple{typeof(hopping_operator),SquareLattice{2},Hopping{1},LandauField})   # time: 0.0160619
    isdefined(LatticeModels, Symbol("#118#119")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#118#119")),Float64})   # time: 0.0159276
    Base.precompile(Tuple{typeof(path_integral),FluxField,SVector{2, Int64},SVector{2, Int64}})   # time: 0.0150953
    Base.precompile(Tuple{typeof(_kws_to_mask),SquareLattice{2},Any})   # time: 0.0145448
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, Int64}}})   # time: 0.0117615
    Base.precompile(Tuple{typeof(dot),LatticeArray{Vector{Float64}, SquareLattice{2}, 1},LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0117156
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Float64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0110625
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis,), Tuple{Int64}},typeof(hopping),Matrix{ComplexF64}})   # time: 0.011016
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(>=), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0107059
    Base.precompile(Tuple{Core.kwftype(typeof(hopping)),NamedTuple{(:axis, :pbc), Tuple{Int64, Vector{Bool}}},typeof(hopping)})   # time: 0.0103254
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Bool, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(<), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0099521
    Base.precompile(Tuple{typeof(hopping_operator),Function,SquareLattice{2},Hopping{1}})   # time: 0.009941
    isdefined(LatticeModels, Symbol("#124#126")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#124#126")),Vector{Float64}})   # time: 0.0090709
    Base.precompile(Tuple{typeof(spectrum),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice}})   # time: 0.0080255
    Base.precompile(Tuple{typeof(getindex),HoneycombLattice,LatticeValue{Bool, :honeycomb}})   # time: 0.0076141
    Base.precompile(Tuple{typeof(copyto!),LatticeValue{Int64, :square},Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}}})   # time: 0.0068864
    Base.precompile(Tuple{Core.kwftype(typeof(Base.Broadcast.dotview)),NamedTuple{(:x1,), Tuple{Int64}},typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square}})   # time: 0.0066486
    Base.precompile(Tuple{typeof(getindex),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0065708
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(+), Tuple{Int64, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(*), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(<=), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(sqrt), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(+), Tuple{Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}, Base.Broadcast.Broadcasted{LVWStyle, Nothing, typeof(Base.literal_pow), Tuple{Base.RefValue{typeof(^)}, LatticeValue{Float64, :square}, Base.RefValue{Val{2}}}}}}}}, Int64}}, Int64}}}},Type{Int64}})   # time: 0.0063367
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{1},Vararg{Hopping{1}}})   # time: 0.0062053
    Base.precompile(Tuple{typeof(radius_vector),SquareLattice{2},LatticeSite{2},LatticeSite{2}})   # time: 0.0058683
    Base.precompile(Tuple{typeof(==),Hopping{1},Hopping{1}})   # time: 0.0056631
    Base.precompile(Tuple{typeof(pairs),LatticeValue{Int64, :square}})   # time: 0.0056102
    Base.precompile(Tuple{typeof(_hopping_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Nothing,Hopping{2},NoField})   # time: 0.0055397
    Base.precompile(Tuple{typeof(Base.Broadcast.dotview),LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0053981
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Matrix{Int64},Int64,Int64})   # time: 0.0046977
    Base.precompile(Tuple{typeof(bonds),SquareLattice{2},Hopping{1}})   # time: 0.0041099
    isdefined(LatticeModels, Symbol("#125#127")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#125#127")),Tuple{Float64, Vector{Float64}}})   # time: 0.0035945
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Function})   # time: 0.0031453
    Base.precompile(Tuple{typeof(check_is_sublattice),HoneycombLattice,HoneycombLattice})   # time: 0.0029196
    Base.precompile(Tuple{typeof(check_is_sublattice),SquareLattice{2},SquareLattice{2}})   # time: 0.0028477
    Base.precompile(Tuple{typeof(getindex),Spectrum{HoneycombLattice, Matrix{ComplexF64}},UnitRange{Int64}})   # time: 0.0026996
    isdefined(LatticeModels, Symbol("#124#126")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#124#126")),Array})   # time: 0.002699
    Base.precompile(Tuple{typeof(_diag_operator!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},Function})   # time: 0.00247
    Base.precompile(Tuple{typeof(_zero_on_basis),SquareLattice{2},SMatrix{1, 1, ComplexF64, 1}})   # time: 0.0024009
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Vector{Float64}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Vector{Float64}},Tuple{LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}}, LatticeArray{Vector{Float64}, SquareLattice{2}, 1}},))
    end
end   # time: 0.0022994
    Base.precompile(Tuple{typeof(_zero_on_basis),HoneycombLattice,SMatrix{2, 2, ComplexF64, 4}})   # time: 0.002273
    Base.precompile(Tuple{typeof(dot),LatticeArray{Vector{Float64}, SquareLattice{2}, 1},LatticeArray{Vector{Float64}, SquareLattice{2}, 1}})   # time: 0.0022513
    Base.precompile(Tuple{typeof(copyto!),LatticeValueWrapper{SubArray{Float64, 1, Vector{Float64}, Tuple{Vector{Int64}}, false}, :square},Base.Broadcast.Broadcasted{Base.Broadcast.DefaultArrayStyle{0}, Tuple{Base.OneTo{Int64}}, typeof(identity), Tuple{Int64}}})   # time: 0.0021889
    Base.precompile(Tuple{typeof(⊗),Matrix{Int64},LatticeValue{Float64, :square}})   # time: 0.0019365
    Base.precompile(Tuple{typeof(^),LatticeOperator{SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}, SquareLattice{2}},Int64})   # time: 0.0019088
    Base.precompile(Tuple{typeof(path_integral),LandauField,SVector{2, Int64},SVector{2, Float64}})   # time: 0.0017033
    isdefined(LatticeModels, Symbol("#30#31")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#30#31")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.00166
    Base.precompile(Tuple{typeof(hopping_operator),HoneycombLattice,Hopping{1}})   # time: 0.0016362
    Base.precompile(Tuple{typeof(*),Int64,MaterializedCurrents})   # time: 0.0016237
    Base.precompile(Tuple{typeof(SpinTightBinding),LatticeValue{Float64, :square},Int64})   # time: 0.0015976
    Base.precompile(Tuple{typeof(setindex!),LatticeValue{Float64, :square},LatticeValue{Float64, :square},LatticeValue{Bool, :square}})   # time: 0.0015815
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, Int64}},Type{Float64}})   # time: 0.0015804
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(*), Tuple{LatticeValue{Float64, :square}, LatticeValue{Int64, :square}}},Type{Float64}})   # time: 0.0015668
    Base.precompile(Tuple{typeof(_zero_on_basis),HoneycombLattice,SMatrix{1, 1, ComplexF64, 1}})   # time: 0.0015558
    isdefined(LatticeModels, Symbol("#30#31")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#30#31")),SubArray{Float64, 2, Array{Float64, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Base.Slice{Base.OneTo{Int64}}, Int64}, true}})   # time: 0.0015476
    Base.precompile(Tuple{typeof(RecipesBase.apply_recipe),AbstractDict{Symbol, Any},LatticeValue{<:Number, LT} where LT})   # time: 0.0014139
    Base.precompile(Tuple{Type{LatticeValue},SquareLattice{2},Vector})   # time: 0.0013263
    isdefined(LatticeModels, Symbol("#128#129")) && Base.precompile(Tuple{getfield(LatticeModels, Symbol("#128#129")),Float64,Vector{Float64}})   # time: 0.0012599
    Base.precompile(Tuple{Type{LatticeRecord{LatticeValue{Float64, :square}}},SquareLattice{2},Vector{Vector{Float64}},Vector{Float64}})   # time: 0.0011731
    Base.precompile(Tuple{typeof(getindex),MaterializedCurrents,LatticeValue{Bool, :square}})   # time: 0.0011649
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{ComplexF64}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}},Tuple{LatticeOperator{Matrix{ComplexF64}, SquareLattice{2}}},))
    end
end   # time: 0.0011599
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, HoneycombLattice},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.001149
    let fbody = try __lookup_kwbody__(which(_unwrap_wlattice, (Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}, LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))) catch missing end
    if !ismissing(fbody)
        precompile(fbody, (Base.Iterators.Pairs{Symbol, Union{}, Tuple{}, NamedTuple{(), Tuple{}}},typeof(_unwrap_wlattice),Function,Basis{SquareLattice{2}},Tuple{Matrix{Float64}},Tuple{LatticeOperator{Matrix{Float64}, SquareLattice{2}}, LatticeOperator{Matrix{Float64}, SquareLattice{2}}},))
    end
end   # time: 0.0011077
    Base.precompile(Tuple{typeof(==),LatticeOperator{Matrix{ComplexF64}, SquareLattice{1}},LatticeOperator{Matrix{ComplexF64}, SquareLattice{1}}})   # time: 0.0010856
    Base.precompile(Tuple{typeof(similar),Base.Broadcast.Broadcasted{LVWStyle, Tuple{Base.OneTo{Int64}}, typeof(imag), Tuple{LatticeValue{ComplexF64, :square}}},Type{Float64}})   # time: 0.0010831
    Base.precompile(Tuple{typeof(setindex!),LatticeOperator{Matrix{ComplexF64}, SquareLattice{1}},SMatrix{1, 1, ComplexF64, 1},Int64,Int64})   # time: 0.0010767
    Base.precompile(Tuple{typeof(_unwrap_from_macro),typeof(ones),Int64,Int64})   # time: 0.001014
end
