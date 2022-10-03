using Logging
import LinearAlgebra: eigen, Hermitian
import Base: length, getindex

_diag_from_macro(lop::LatticeOperator, ::AbstractLattice, m::AbstractMatrix) =
    _diag_operator!(lop, m)
_diag_from_macro(lop::LatticeOperator, l::AbstractLattice, f::Function) =
    _diag_operator!(lop, _propagate_lattice_args(f, l))
_diag_from_macro(::LatticeOperator, ::AbstractLattice, ::T) where T =
    error("unextected argument type $T in @diag")
_diag_from_macro(l::AbstractLattice, arg::Any) =
    _diag_from_macro(_zero_on_basis(l, arg), l, arg)

_hops_from_macro(lop::LatticeOperator, l::AbstractLattice, pr_fun::Function, hop::Hopping, field::AbstractField=NoField()) =
    _hopping_operator!(lop, pr_fun, l, hop, field)
_hops_from_macro(l::AbstractLattice, pr_fun::Function, hop::Hopping, field::AbstractField=NoField()) =
    _hops_from_macro(_zero_on_basis(l, hop.hop_operator), l, pr_fun, hop, field)

macro hamiltonian(block)
    lattice_sym = nothing
    field_sym = nothing
    ham_block = Expr(:block)
    assign_flag = true
    if Meta.isexpr(block, :block)
        for line in block.args
            if Meta.isexpr(line, :(:=))
                key, value = line.args
                if key === :lattice
                    lattice_sym !== nothing &&
                        error("cannot overwrite key $key")
                    lattice_sym = :($(esc(value)))
                elseif key === :field
                    field_sym !== nothing &&
                        error("cannot overwrite key $key")
                    field_sym = :($(esc(value)))
                else
                    @warn "skipping unused key $key"
                end
            elseif Meta.isexpr(line, :macrocall)
                lattice_sym === nothing && error("define lattice first")
                macro_name, macro_args... = line.args
                macro_args = [a for a in macro_args if (a isa Expr || a isa Symbol)]
                if macro_name === Symbol("@diag")
                    length(macro_args) != 1 &&
                        error("@diag accepts only one argument")
                    macro_arg = only(macro_args)
                    push!(ham_block.args, :(
                        _diag_from_macro($lattice_sym, $(esc(macro_arg)))
                    ))
                elseif macro_name === Symbol("@hop") || macro_name === Symbol("@hopping")
                    lambda_i = findfirst(a -> Meta.isexpr(a, (:function, :->)), macro_args)
                    if lambda_i !== nothing
                        pr_lambda = :(_propagate_lattice_args($(esc(macro_args[lambda_i])), $lattice_sym))
                        popat!(macro_args, lambda_i)
                    else
                        pr_lambda = :_always_true_on_lattice
                    end
                    for arg in macro_args
                        if Meta.isexpr(arg, :(=))
                            arg.head = :kw
                        end
                    end

                    hopcall = :(Hopping($(esc.(macro_args)...)))
                    push!(ham_block.args, :(
                            _hops_from_macro($lattice_sym, $pr_lambda, $hopcall)
                    ))
                else error("unexpected macro call $macro_name in @hamiltonian")
                end
                assign_flag = false
            end

        end
        isfirst_statement = true
        for statement in ham_block.args
            !Meta.isexpr(statement, :call) && continue
            if statement.args[1] === :_hops_from_macro && field_sym !== nothing
                push!(statement.args, field_sym)
            end
            if isfirst_statement
                ham_block.args[1] = :(H = $statement)
                isfirst_statement = false
            else
                insert!(statement.args, 2, :H)
            end
        end
        return ham_block
    end
end

struct Spectrum
    basis::Basis
    states::AbstractMatrix
    energies::AbstractVector
    function Spectrum(basis::Basis, states::AbstractMatrix, energies::AbstractVector)
        @assert length(energies) == size(states)[1] "inconsistent energies list length"
        @assert length(basis) == size(states)[2] "inconsistent basis dimensionality"
        new(basis, states, energies)
    end
end

function spectrum(lop::LatticeOperator{Matrix{T}} where T)
    @assert all(isfinite.(lop.operator)) "operator matrix has non-finite values (NaN of Inf)"
    vals, vecs = eigen(Hermitian(lop.operator))
    return Spectrum(lop.basis, vecs, vals)
end

length(sp::Spectrum) = length(sp.energies)
getindex(sp::Spectrum, i::Int) = LatticeVecOrMat(sp.basis, sp.states[:, i])
getindex(sp::Spectrum, mask::AbstractVector{Bool}) =
    Spectrum(sp.basis, sp.states[:, mask], sp.energies[mask])
function getindex(sp::Spectrum; E::Number)
    min_e_dst = abs(E - sp.energies[1])
    i = 1
    for j in 2:length(sp)
        if abs(E - sp.energies[j]) < min_e_dst
            i = j
            min_e_dst = abs(E - sp.energies[j])
        end
    end
    LatticeVecOrMat(sp.basis, sp.states[:, i])
end

function show(io::IO, ::MIME"text/plain", sp::Spectrum)
    println(io, "Spectrum with $(length(sp)) eigenstates")
    println(io, "Eigenvalues in range $(minimum(sp.energies)) .. $(maximum(sp.energies))")
end

function projector(sp::Spectrum)
    return LatticeVecOrMat(sp.basis, sp.states * sp.states')
end

function projector(f::Function, sp::Spectrum)
    return LatticeVecOrMat(sp.basis, sp.states * (f.(sp.energies) .* sp.states'))
end

function filled_projector(sp::Spectrum, fermi_level=0)
    return projector(E -> E < fermi_level, sp)
end
