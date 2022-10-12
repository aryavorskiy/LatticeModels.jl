using Logging
import LinearAlgebra: eigen, Hermitian, eigvals, eigvecs
import Base: length, getindex

_diag_from_macro(lop::LatticeOperator, ::Lattice, selector) =
    _diag_operator!(lop, selector)
_diag_from_macro(l::Lattice, selector) =
    _diag_from_macro(_zero_on_basis(l, selector), l, selector)

_hops_from_macro(lop::LatticeOperator, ::Lattice, selector, hop::Hopping, field::AbstractField=NoField()) =
    _hopping_operator!(lop, selector, hop, field)
_hops_from_macro(l::Lattice, selector, hop::Hopping, field::AbstractField=NoField()) =
    _hops_from_macro(_zero_on_basis(l, hop.hop_operator), l, selector, hop, field)

_lazy_tp(m::AbstractMatrix, lv::LatticeValue) = TensorProduct(lv, m)
_lazy_tp(lv::LatticeValue, m::AbstractMatrix) = TensorProduct(lv, m)

function _hamiltonian_block(block::Expr)
    lattice_sym = nothing
    field_sym = nothing
    ham_block = Expr(:block)
    assign_flag = true
    for line in block.args
        if Meta.isexpr(line, :(:=))
            key, value = line.args
            if key === :lattice
                lattice_sym !== nothing &&
                    error("cannot overwrite key :$key")
                lattice_sym = :($(esc(value)))
            elseif key === :field
                field_sym !== nothing &&
                    error("cannot overwrite key :$key")
                field_sym = :($(esc(value)))
            else
                @warn "skipping unused key :$key"
            end
        elseif Meta.isexpr(line, :macrocall)
            lattice_sym === nothing && error("define lattice first")
            macro_name, macro_args... = line.args
            macro_args = [a for a in macro_args if (a isa Expr || a isa Symbol)]
            if macro_name === Symbol("@diag")
                length(macro_args) != 1 &&
                    error("@diag accepts only one argument")
                macro_arg = only(macro_args)
                if Meta.isexpr(macro_arg, :call, 3) && macro_arg.args[1] == :⊗
                    macro_arg.args[1] = :(LatticeModels._lazy_tp)
                end
                push!(ham_block.args, :(
                    _diag_from_macro($lattice_sym, $(esc(macro_arg)))
                ))
            elseif macro_name === Symbol("@hop") || macro_name === Symbol("@hopping")
                args_is = findall(a -> !Meta.isexpr(a, :(=), 2), macro_args)
                hop_operator = 1
                local selector_sym = :nothing
                if length(args_is) ≥ 1
                    hop_operator = esc(macro_args[args_is[1]])
                end
                if length(args_is) ≥ 2
                    selector_sym = esc(macro_args[args_is[2]])
                end
                deleteat!(macro_args, args_is)
                for arg in macro_args
                    arg.head = :kw
                end
                hopcall = :(hopping($hop_operator, $(esc.(macro_args)...)))
                push!(ham_block.args, :(
                    _hops_from_macro($lattice_sym, $selector_sym, $hopcall)
                ))
            else
                error("unexpected macro call $macro_name in @hamiltonian")
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
    ham_block
end

macro hamiltonian(expr)
    _hamiltonian_block(expr)
end

struct Spectrum{LT<:Lattice,MT<:AbstractMatrix}
    basis::Basis{LT}
    states::MT
    energies::Vector{Float64}
    function Spectrum(basis::Basis{LT}, states::MT, energies::AbstractVector) where {LT,MT}
        length(basis) != size(states)[1] && error("inconsistent basis dimensionality")
        length(energies) != size(states)[2] && error("inconsistent energies list length")
        new{LT,MT}(basis, states, energies)
    end
end

function spectrum(lop::LatticeOperator{LT, Matrix{T}} where {LT,T})
    !all(isfinite.(lop.operator)) && error("NaN of Inf in operator matrix")
    vals, vecs = eigen(Hermitian(lop.operator))
    Spectrum(lop.basis, vecs, vals)
end

eigvals(sp::Spectrum) = sp.energies
eigvecs(sp::Spectrum) = sp.states

length(sp::Spectrum) = length(sp.energies)
getindex(sp::Spectrum, i::Int) = LatticeArray(sp.basis, sp.states[:, i])
function getindex(sp::Spectrum; E::Number)
    min_e_dst = abs(E - sp.energies[1])
    i = 1
    for j in 2:length(sp)
        if abs(E - sp.energies[j]) < min_e_dst
            i = j
            min_e_dst = abs(E - sp.energies[j])
        end
    end
    LatticeArray(sp.basis, sp.states[:, i])
end
getindex(sp::Spectrum, mask) =
    Spectrum(sp.basis, sp.states[:, mask], sp.energies[mask])

function show(io::IO, ::MIME"text/plain", sp::Spectrum)
    println(io, "Spectrum with $(length(sp)) eigenstates")
    println(io, "Eigenvalues in range $(minimum(sp.energies)) .. $(maximum(sp.energies))")
end

projector(sp::Spectrum) = LatticeArray(sp.basis, sp.states * sp.states')

projector(f::Function, sp::Spectrum) =
    LatticeArray(sp.basis, sp.states * (f.(sp.energies) .* sp.states'))

filled_projector(sp::Spectrum, fermi_level=0) = projector(E -> E < fermi_level, sp)
