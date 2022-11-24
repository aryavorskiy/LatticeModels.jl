using Logging
import LinearAlgebra: eigen, Hermitian, eigvals, eigvecs
import Base: length, getindex

_lazy_tp(m::AbstractMatrix, lv::LatticeValue) = TensorProduct(lv, m)
_lazy_tp(lv::LatticeValue, m::AbstractMatrix) = TensorProduct(lv, m)

function _hamiltonian_block(block::Expr)
    assignments = Dict{Symbol, Any}()
    for line in block.args
        if Meta.isexpr(line, :(:=), 2)
            k, v = line.args
            if k in keys(assignments)
                error("cannot overwrite key :$k with value '$v' ('$(assignments[k])' assigned before)")
            end
            assignments[k] = esc(v)
        end
    end
    !(:lattice in keys(assignments)) &&
        error("lattice not defined")
    lattice_sym = assignments[:lattice]
    field_sym = get(assignments, :field, :(NoField()))
    arrtype_sym = get(assignments, :arrtype, :(Matrix{ComplexF64}))
    dim_int_sym = get(assignments, :dims_internal, 1)
    ham_block = quote
        H = _zero_on_basis($lattice_sym, $dim_int_sym, $arrtype_sym)
    end
    for line in block.args
        if Meta.isexpr(line, :macrocall)
            macro_name, macro_args... = line.args
            macro_args = [a for a in macro_args if !(a isa QuoteNode || a isa LineNumberNode)]
            if macro_name === Symbol("@diag")
                length(macro_args) != 1 &&
                    error("@diag accepts only one argument")
                macro_arg = only(macro_args)
                if Meta.isexpr(macro_arg, :call, 3) && macro_arg.args[1] == :⊗
                    macro_arg.args[1] = :(LatticeModels._lazy_tp)
                end
                push!(ham_block.args, :(
                    _diag_operator!(H, $(esc(macro_arg)))
                ))
            elseif macro_name === Symbol("@hop")
                args_is = findall(a -> !Meta.isexpr(a, :(=), 2), macro_args)
                hop_operator = 1
                selector_sym = :nothing
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
                    _hopping_operator!(H, $selector_sym, $hopcall, $field_sym)
                ))
            else
                error("unexpected macro call $macro_name in @hamiltonian")
            end
        end
    end
    ham_block
end

"""
    @hamiltonian block

Creates a hamiltonian according to the rules defined in the `block`.

Each line in the block must be a `:=` assignment or a macro-like diagonal/hopping operator description.

The lattice on and the magnetic field for the hamiltonian can be set by assigning `lattice` and `field`.
You may also need to set the internal dimension count if it is not equal to 1 - use the `dims_internal` keyword.
The `arrtype` key sets the type of the returned array.

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
    dims_internal := 2
    @diag (@. abs(x) < 2) ⊗ [1 0; 0 -1]
    @diag randn(l) .* 0.1       # Add Gaussian noise
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
basis(sp::Spectrum) = sp.basis

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

@doc raw"""
    dos(spectrum, δ)

Generates a function to calculate density of states, which is defined as
$\text{tr}\left(\frac{1}{\hat{H} - E - i\delta}\right)$ and can be understood as a sum
of Lorenz distributions with width equal to $\delta$.
"""
dos(sp::Spectrum, δ::Real) = (E -> imag(sum(1 ./ (eigvals(sp) .- (E + im * δ)))))

@doc raw"""
    ldos(spectrum, E, δ)

Calculates local density of states, which is defined as the imaginary part of partial trace of
$\frac{1}{\hat{H} - E - i\delta}$ operator.
"""
function ldos(sp::Spectrum, E::Real, δ::Real)
    Es = eigvals(sp)
    Vs = eigvecs(sp)
    l = lattice(sp)
    N = dims_internal(sp)
    inves = imag.(1 ./ (Es .- (E + im * δ)))'
    LatticeValue(l, [sum(abs2.(Vs[(i-1)*N+1:i*N, :]) .* inves) for i in 1:length(l)])
end

"""
    ldos(spectrum, δ)

Generates a function that accepts the energy `E` and returns `ldos(spectrum, E, δ)`.
Use this if you want to find the LDOS for many different values of `E` -
the produced function is optimized and reduces overall computation time dramatically.
"""
function ldos(sp::Spectrum, δ::Real)
    Es = eigvals(sp)
    l = lattice(sp)
    N = dims_internal(sp)
    density_sums = reshape(
        sum(reshape(abs2.(eigvecs(sp)), (N, :, length(sp))), dims=1), (:, length(sp)))
    E -> LatticeValue(l, vec(sum(density_sums .* imag.(1 ./ (Es .- (E + im * δ)))', dims=2)))
end