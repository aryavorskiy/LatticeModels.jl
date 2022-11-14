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

"""
    @hamiltonian block

Creates a hamiltonian according to the rules defined in the `block`.

Each line in the block must be a `:=` assignment or a macro-like diagonal/hopping operator description.

The lattice on and the magnetic field for the hamiltonian can be set by assigning `lattice` and `field`.

`@diag` stands for a diagonal part of the hamiltonian - after this you can use a matrix
(representing the operator affecting the internal state),
a function or a `⊗` tensor product notation.

`@hop` stands for a hopping part of the hamiltonian - list arguments to pass to the `hopping` function after this macrocall.

## Examples

This is how a Chern insulator hamiltonian is generated.
```julia
l = SquareLattice(10, 10)
x, y = coord_values(l)
H = @hamiltonian begin
    lattice := l
    field := LandauField(0.5)   # Landau-calibrated uniform magnetic field, 0.5 flux quanta per 1×1
    @diag (@. abs(x) < 2) ⊗ [1 0; 0 -1]
    @hop [1 im; im -1] / 2 axis = 1
    @hop [1 1; -1 -1] / 2 axis = 2
end
```
"""
macro hamiltonian(expr)
    _hamiltonian_block(expr)
end

"""
    Spectrum{LT, MT} where {LT<:Lattice, MT<:AbstractMatrix}

Eigenvalues and eigenvectors for some operator.
"""
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

const LatticeOperatorMT{MT} = LatticeOperator{LT,<:MT} where {LT}

"""
    spectrum(operator)

Finds eigenvalues and eigenvectors for a `LatticeOperator` and stores in in a Spectrum.

!!! note
    This method finds eigenvalues and eigenvectors using `LinearAlgebra.eigen`, which can be not defined for some array types.
    Consider redefining it for your array type or constructing the Spectrum object explicitly.
"""
function spectrum(lop::LatticeOperator)
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
    i = argmin(@. abs(E - sp.energies))
    LatticeArray(sp.basis, sp.states[:, i])
end
getindex(sp::Spectrum, mask) =
    Spectrum(sp.basis, sp.states[:, mask], sp.energies[mask])

function show(io::IO, ::MIME"text/plain", sp::Spectrum)
    println(io, "Spectrum with $(length(sp)) eigenstates")
    println(io, "Eigenvalues in range $(minimum(sp.energies)) .. $(maximum(sp.energies))")
end

"""
    projector(spectrum)

Creates a `LatticeOperator` that projects onto the eigenvectors of the spectrum
"""
projector(sp::Spectrum) = LatticeArray(sp.basis, sp.states * sp.states')

"""
    projector(fun, spectrum)

Creates a `LatticeOperator` that projects onto the eigenvectors of the spectrum
with amplitude defined by the `fun` functions, which takes the eigenvalue and returns a number (or a boolean).
"""
projector(f::Function, sp::Spectrum) =
    LatticeArray(sp.basis, sp.states * (f.(sp.energies) .* sp.states'))

"""
    filled_projector(spectrum[, fermi_level])

Creates a `LatticeOperator` that projects onto the eigenvectors which have eigenvalues less than `fermi_level` (0 by default).
"""
filled_projector(sp::Spectrum, fermi_level=0) = projector(E -> E < fermi_level, sp)
