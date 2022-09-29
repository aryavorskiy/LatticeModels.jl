using Logging
import LinearAlgebra: eigen, Hermitian
import Base: length, getindex

macro hamiltonian(block)
    lattice_sym = nothing
    field_sym = nothing
    ham_block = :(+())
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
                        error("wrong arguments for diag part; expected 1, got $(length(macro_args))")
                    macro_arg = only(macro_args)
                    if Meta.isexpr(macro_arg, (:->, :function))
                        push!(ham_block.args, :(
                            diag_operator($(esc(macro_arg)), $lattice_sym)
                        ))
                    else
                        push!(ham_block.args, :(
                            diag_operator($lattice_sym, $(esc(macro_arg)))
                        ))
                    end
                elseif macro_name === Symbol("@hop") || macro_name === Symbol("@hopping")
                    hopp_kwargs = Dict{Symbol, Any}(a.args for a in macro_args if Meta.isexpr(a, :(=)))
                    lambda_i = findfirst(a -> Meta.isexpr(a, (:function, :->)), macro_args)
                    op_i = findfirst(a -> !Meta.isexpr(a, (:function, :->, :(=))), macro_args)
                    if op_i !== nothing && !(:hop_operator in keys(hopp_kwargs))
                        hopp_kwargs[:hop_operator] = macro_args[op_i]
                    end
                    hopcall = :(Hopping())
                    append!(hopcall.args, [Expr(:kw, k, :($(esc(v)))) for (k, v) in hopp_kwargs])
                    if lambda_i === nothing
                        push!(ham_block.args, :(
                            hopping_operator($lattice_sym, $hopcall)
                        ))
                    else
                        push!(ham_block.args, :(
                            hopping_operator($(esc(macro_args[lambda_i])), $lattice_sym, $hopcall)
                        ))
                    end
                elseif macro_name === Symbol("@hopping_operator")
                    error("usage of @hopping_operator macro is forbidden in hamiltonian generator; please consult the manual")
                else error("unexpected macro call $macro_name")
                end
            end

        end
        if field_sym !== nothing
            for statement in ham_block.args[2:end]
                if Meta.isexpr(statement, :call) && statement.args[1] === :hopping_operator
                    push!(
                        statement.args,
                        Expr(:kw, :field, field_sym)
                    )
                end
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

function Spectrum(lop::LatticeOperator{Matrix{T}} where T)
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
