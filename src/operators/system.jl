import QuantumOpticsBase: basis, check_samebases, OccupationNumbers, FermionStatistics, BosonStatistics

struct Sample{LT, BasisT}
    lat::LT
    internal::BasisT
end
function Sample(l::AbstractLattice, internal::IT=nothing) where {IT<:Nullable{Basis}}
    Sample{typeof(l), IT}(l, internal)
end
const SampleWithoutInternal{LT} = Sample{LT, Nothing}
const SampleWithInternal{LT} = Sample{LT, <:Basis}

Base.:(==)(s1::Sample, s2::Sample) = s1.lat == s2.lat && s1.internal == s2.internal
Base.length(sample::Sample) = length(sample.lat) * length(sample.internal)
Base.length(sample::SampleWithoutInternal) = length(sample.lat)
QuantumOpticsBase.basis(sample::SampleWithInternal) = sample.internal ⊗ LatticeBasis(sample.lat)
QuantumOpticsBase.basis(sample::SampleWithoutInternal) = LatticeBasis(sample.lat)

Base.getindex(sample::Sample, args...; kw...) =
    Sample(getindex(sample.lat, args...; kw...), sample.internal)

function Base.show(io::IO, mime::MIME"text/plain", sample::Sample)
    hasinternal(sample) && print(io, "(")
    summary(io, lattice(sample))
    if hasinternal(sample)
        print(io, ") ⊗ ")
        show(io, mime, internal_basis(sample))
    end
end

"""
Returns the `Sample` of the object.

Define this function for your type to implement `Sample` API.

!!! info
    This function can be considered stable internal API. Feel free to use it in your packages.
"""
sample(lb::LatticeBasis) = Sample(lb.lat)
sample(b::CompositeLatticeBasis) = Sample(b.bases[2].lat, b.bases[1])
sample(mb::ManyBodyBasis) = sample(mb.onebodybasis)
sample(b::Basis) = throw(ArgumentError("Basis has no lattice: $b"))
sample(state::StateType) = sample(basis(state))
sample(op::AbstractLatticeOperator) = sample(basis(op))
lattice(sample::Sample) = sample.lat
lattice(any) = lattice(sample(any))
internal_basis(sample::SampleWithInternal) = sample.internal
internal_basis(::SampleWithoutInternal) = throw(ArgumentError("Sample has no internal basis"))
internal_basis(any) = internal_basis(sample(any))
internal_length(sample::SampleWithInternal) = length(sample.internal)
internal_length(sample::SampleWithoutInternal) = 1
internal_length(any) = internal_length(sample(any))
internal_one(sample::SampleWithInternal) = one(sample.internal)
internal_one(sample::SampleWithoutInternal) = 1
internal_one(any) = internal_one(sample(any))
hasinternal(::SampleWithInternal) = true
hasinternal(::SampleWithoutInternal) = false
hasinternal(any) = hasinternal(sample(any))

@enum ParticleStatistics begin
    FermiDirac = 1
    BoseEinstein = -1
end

abstract type System{SampleT} end
abstract type OneParticleBasisSystem{SampleT} <: System{SampleT} end
sample(sys::System) = sys.sample
struct OneParticleSystem{SampleT} <: OneParticleBasisSystem{SampleT}
    sample::SampleT
    T::Float64
    statistics::ParticleStatistics
    OneParticleSystem(sample::SampleT, T::Real=0, statistics=FermiDirac) where SampleT =
        new{SampleT}(sample, Float64(T), statistics)
end
OneParticleSystem(l::AbstractLattice, b::Nullable{Basis}=nothing; T=0, statistics=FermiDirac) =
    OneParticleSystem(Sample(l, b), T, statistics)

QuantumOpticsBase.tensor(l::AbstractLattice, b::Basis) = OneParticleSystem(l, b)
QuantumOpticsBase.tensor(b::Basis, l::AbstractLattice) = OneParticleSystem(l, b)

struct FixedMu{SampleT} <: OneParticleBasisSystem{SampleT}
    sample::SampleT
    chempotential::Float64
    statistics::ParticleStatistics
    T::Float64
end
function FixedMu(sample::SampleT, mu; statistics=FermiDirac, T = 0) where SampleT
    FixedMu{SampleT}(sample, mu, statistics, T)
end

struct FixedN{SampleT} <: OneParticleBasisSystem{SampleT}
    sample::SampleT
    nparticles::Int
    statistics::ParticleStatistics
    T::Float64
end
function FixedN(sample::SampleT, N; statistics=FermiDirac, T = 0) where SampleT
    FixedN{SampleT}(sample, N, statistics, T)
end

function Base.:(==)(sys1::OneParticleBasisSystem, sys2::OneParticleBasisSystem)
    (sys1.sample == sys2.sample && sys1.T == sys2.T) || return false
    if sys1 isa OneParticleSystem && sys2 isa OneParticleSystem
        return true
    elseif sys1 isa FixedMu && sys2 isa FixedMu
        return sys1.chempotential == sys2.chempotential && sys1.statistics == sys2.statistics
    elseif sys1 isa FixedN && sys2 isa FixedN
        return sys1.nparticles == sys2.nparticles && sys1.statistics == sys2.statistics
    else
        return false
    end
end
function Base.show(io::IO, mime::MIME"text/plain", sys::OneParticleBasisSystem)
    if sys isa OneParticleSystem
        print(io, "One particle on ")
    elseif sys isa FixedN
        noun = sys.statistics == FermiDirac ? "fermion" : "boson"
        print(io, fmtnum(sys.nparticles, "non-interacting " * noun), " on ")
    elseif sys isa FixedMu
        noun = sys.statistics == FermiDirac ? "fermion" : "boson"
        print(io, "Non-interactng $(noun)s with fixed μ=",
            trunc(sys.chempotential, digits=2), " on ")
    end
    show(io, mime, sys.sample)
end

abstract type ManyBodySystem{SampleT} <: System{SampleT} end

struct NParticles{OccT,SampleT,NPT} <: ManyBodySystem{SampleT}
    sample::SampleT
    nparticles::NPT
    statistics::ParticleStatistics
    T::Float64
end

"""
    NParticles(lat[, internal], N[; T=0, statistics=FermiDirac, occupations_type])
    NParticles(sys, N[; T=0, statistics=FermiDirac, occupations_type])

Create a manybody system with a given lattice and a given number of particles.

## Arguments
- `lat`: the lattice of the system.
- `internal`: The basis for the internal degrees of freedom.
- `sys`: a one-particle system.
- `N`: the number of particles in the system.

## Keyword Arguments
- `T`: the temperature of the system. Default is `0`.
- `statistics`: the statistics of the particles. Default is `FermiDirac`.
- `occupations_type`: The occupations type for the many-body operator. By default, the
    occupation numbers are stored in vectors, but you can use, for example, set it to
    `FermionBitstring`s for better performance on fermion systems.

## Example
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(3, 3);

julia> NParticles(lat, 4, statistics=BoseEinstein)
NParticles(4 bosons) on 9-site SquareLattice in 2D space
```
"""
NParticles(sample::SampleT, nparticles; statistics = FermiDirac, T = 0,
    occupations_type=nothing) where {SampleT<:Sample} =
    NParticles{occupations_type, SampleT, typeof(nparticles)}(sample, nparticles, statistics, T)
NParticles(onep::OneParticleSystem, n; kw...) = NParticles(onep.sample, n; T = onep.T, kw...)
NParticles(l::AbstractLattice, n; kw...) = NParticles(Sample(l), n; kw...)
NParticles(l::AbstractLattice, b::Basis, n; kw...) = NParticles(Sample(l, b), n; kw...)
Base.:(==)(sys1::NParticles, sys2::NParticles) =
    sys1.sample == sys2.sample && sys1.nparticles == sys2.nparticles &&
    sys1.statistics == sys2.statistics && sys1.T == sys2.T
function Base.show(io::IO, mime::MIME"text/plain", sys::NParticles{OccT}) where OccT
    noun = sys.statistics == FermiDirac ? "fermion" : "boson"
    n = sys.nparticles
    print(io, "NParticles(", n isa Int ? fmtnum(n, noun) : string(n) * " $(noun)(s)",
        OccT === nothing ? "" : ", occupations_type=$OccT", ") on ")
    show(io, mime, sys.sample)
end

struct ManyBodyBasisSystem{SampleT, MBT<:ManyBodyBasis} <: ManyBodySystem{SampleT}
    sample::SampleT
    mb::MBT
    T::Float64
    function ManyBodyBasisSystem(mb::MBT; T=0) where {MBT<:ManyBodyBasis}
        s = sample(mb)
        new{typeof(s), MBT}(s, mb, T)
    end
end
function Base.show(io::IO, mime::MIME"text/plain", sys::ManyBodyBasisSystem)
    print(io, "Many-body system on ")
    show(io, mime, sys.sample)
    print(io, " ($(length(sys.mb)) states, T=$(sys.T))")
end
function Base.union(mbs1::ManyBodyBasisSystem, mbs2::ManyBodyBasisSystem)
    sample(mbs1) == sample(mbs2) || throw(ArgumentError("Incompatible systems"))
    mbs1.T == mbs2.T || throw(ArgumentError("Incompatible temperatures"))
    new_basis = ManyBodyBasis(onebodybasis(mbs1), union(mbs1.mb.occupations, mbs2.mb.occupations))
    return ManyBodyBasisSystem(new_basis, T=mbs1.T)
end

"""
    System(lat[, internal; T, μ, N, statistics])

Create a system with a given lattice and optionally internal degrees of freedom.


## Arguments
- `lat`: the lattice of the system.
- `internal`: The basis for the internal degrees of freedom.

## Keyword Arguments
- `T`: the temperature of the system. Default is `0`.
- `μ`: the chemical potential of the system. Use `mu` synonym if Unicode input is not available.
- `N`: the number of particles in the system.
- `statistics`: the statistics of the particles. Default is `FermiDirac`.

## Example
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(3, 3);

julia> System(lat)
One particle on 9-site SquareLattice in 2D space

julia> System(lat, N=4, statistics=BoseEinstein)
4 non-interacting bosons on 9-site SquareLattice in 2D space

julia> System(lat, mu=0, statistics=BoseEinstein)
Non-interactng bosons with fixed μ=0.0 on 9-site SquareLattice in 2D space
```
"""
function System(sample::Sample; μ = nothing, mu = μ, N = nothing, T = 0, statistics=FermiDirac)
    if mu !== nothing && N === nothing
        return FixedMu(sample, mu, statistics=statistics, T=T)
    elseif N !== nothing && mu === nothing
        return FixedN(sample, N, statistics=statistics, T=T)
    elseif N === mu === nothing
        return OneParticleSystem(sample, T)
    else
        throw(ArgumentError("cannot specify both N and μ"))
    end
end
System(onep::OneParticleSystem; kw...) = System(onep.sample; T=onep.T, kw...)
function System(args...; μ = nothing, mu = μ, N = nothing, statistics=FermiDirac, T=0, kw...)
    isempty(kw) || throw(ArgumentError("Unsupported keyword arguments " * join(keys(kw), ", ")))
    System(Sample(args...), mu=mu, N=N, T=T, statistics=statistics)
end

"""
    System(mb[; T])

Create a system with a given many-body basis `mb`.

This function is used to create a many-body system from an arbitrary many-body basis with a
lattice.

## Example
```jldoctest
julia> using LatticeModels

julia> lat = SquareLattice(3, 3);

julia> bas = SpinBasis(1//2) ⊗ LatticeBasis(lat);

julia> mbas = ManyBodyBasis(bas, fermionstates(bas, 2));

julia> System(mbas, T=2)
Many-body system on (9-site SquareLattice in 2D space) ⊗ Spin(1/2) (153 states, T=2.0)
```
"""
System(mb::ManyBodyBasis; T=0) = ManyBodyBasisSystem(mb; T=T)
_occtype(::NParticles{OccT}) where {OccT} = OccT
function occupations(np::NParticles, occupations_type::Union{Type,Nothing}=_occtype(np))
    n = np.nparticles
    if np.statistics == FermiDirac
        new_occ = occupations_type !== nothing ? occupations_type : OccupationNumbers{FermionStatistics, Int}
        fermionstates(new_occ, length(np.sample), n isa Int ? n : collect(n))
    elseif np.statistics == BoseEinstein
        new_occ = occupations_type !== nothing ? occupations_type : OccupationNumbers{BosonStatistics, Int}
        bosonstates(new_occ, length(np.sample), n isa Int ? n : collect(n))
    else
        throw(ArgumentError("Unsupported statistics: $(np.statistics)"))
    end
end
occupations(sys::ManyBodyBasisSystem) = sys.mb.occupations

onebodybasis(sys::System) = basis(sys.sample)
QuantumOpticsBase.basis(sys::OneParticleBasisSystem) = onebodybasis(sys)
QuantumOpticsBase.basis(sys::NParticles) = ManyBodyBasis(basis(sys.sample), occupations(sys))
QuantumOpticsBase.basis(sys::ManyBodyBasisSystem) = sys.mb

Base.zero(sys::System) = zero(basis(sys))

macro accepts_system(fname, default_basis=nothing)
    esc(quote
        $fname(sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(l::AbstractLattice, bas::Basis, args...; kw...) =
            $fname(Sample(l, bas), args...; kw...)
        $fname(l::AbstractLattice, args...; kw...) =
            $fname(Sample(l, $default_basis), args...; kw...)
    end)
end

macro accepts_t(fname)
    esc(quote
    $fname(args...; kw...) = $fname(ComplexF64, args...; kw...)
    $fname(::Type, ::Type, args...; kw...) =
        throw(MethodError($fname, args))
    end)
end

macro accepts_system_t(fname, default_basis=nothing)
    esc(quote
        $fname(type::Type, sample::Sample, args...; T = 0, μ = nothing, mu = μ, N = nothing, statistics = nothing, kw...) =
            $fname(type, System(sample, T=T, mu=mu, N=N, statistics=statistics), args...; kw...)
        $fname(type::Type, l::AbstractLattice, bas::Basis, args...; kw...) =
            $fname(type, Sample(l, bas), args...; kw...)
        $fname(type::Type, l::AbstractLattice, args...; kw...) =
            $fname(type, Sample(l, $default_basis), args...; kw...)
        $fname(args...; kw...) = $fname(ComplexF64, args...; kw...)
        $fname(::Type, ::Type, args...; kw...) =
            throw(MethodError($fname, args))
    end)
end

"""
    Hamiltonian <: QuantumOpticsBase.DataOperator

A wrapper for a Hamiltonian operator. Contains the operator matrix and the system it acts on.

---
    Hamiltonian(sys, op)

Create a Hamiltonian operator for a given system and a given operator.

## Arguments
- `sys`: the system the Hamiltonian acts on.
- `op`: the operator matrix.

## Example
```jldoctest
julia> using LatticeModels

julia> l = SquareLattice(4, 4);

julia> H = tightbinding_hamiltonian(l)
Hamiltonian(dim=16x16)
System: One particle on 16-site SquareLattice in 2D space
16×16 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 48 stored entries:
⎡⠪⡢⠑⢄⠀⠀⠀⠀⎤
⎢⠑⢄⠪⡢⠑⢄⠀⠀⎥
⎢⠀⠀⠑⢄⠪⡢⠑⢄⎥
⎣⠀⠀⠀⠀⠑⢄⠪⡢⎦

julia> l2 = SquareLattice(5, 5);

julia> H2 = tightbinding_hamiltonian(l2)
Hamiltonian(dim=25x25)
System: One particle on 25-site SquareLattice in 2D space
25×25 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 80 stored entries:
⎡⠪⡢⡈⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠢⡈⠠⡢⡈⠢⡀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠈⠢⡈⠊⡠⡈⠢⡀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠢⡈⠪⠂⡈⠢⡀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠢⡈⠪⡢⠈⠢⡀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠢⡀⠪⡢⡀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠀⠈⠀⎦

julia> H + H2
ERROR: Incompatible Hamiltonians:
  #1: One particle on 16-site SquareLattice in 2D space
  #2: One particle on 25-site SquareLattice in 2D space
[...]
```
"""
struct Hamiltonian{SystemT, BasisT, T} <: DataOperator{BasisT, BasisT}
    sys::SystemT
    basis_l::BasisT
    basis_r::BasisT
    data::T
end
Hamiltonian(sys::System, data::AbstractMatrix) = Hamiltonian(sys, basis(sys), basis(sys), data)
Hamiltonian(sys::System, op::Operator) = Hamiltonian(sys, op.data)

QuantumOpticsBase.Operator(ham::Hamiltonian) = Operator(ham.basis_l, ham.data)
sample(ham::Hamiltonian) = sample(ham.sys)
lattice(ham::Hamiltonian) = lattice(sample(ham))

Base.:(*)(op::Operator{B1, B2}, ham::Hamiltonian{Sys, B2}) where {Sys, B1, B2} = op * Operator(ham)
Base.:(*)(ham::Hamiltonian{Sys, B2}, op::Operator{B1, B2}) where {Sys, B1, B2} = Operator(ham) * op
function Base.:(+)(op::Operator{B, B}, ham::Hamiltonian{Sys, B}) where {Sys, B}
    QuantumOpticsBase.check_samebases(basis(op), ham.basis_l)
    return Hamiltonian(ham.sys, op.data + ham.data)
end
Base.:(+)(ham::Hamiltonian{Sys, B}, op::Operator{B, B}) where {Sys, B} = op + ham
Base.:(-)(ham::Hamiltonian) = Hamiltonian(ham.sys, -ham.data)
function Base.:(-)(op::Operator{B, B}, ham::Hamiltonian{Sys, B}) where {Sys, B}
    QuantumOpticsBase.check_samebases(basis(op), ham.basis_l)
    return Hamiltonian(ham.sys, op.data - ham.data)
end
function Base.:(-)(ham::Hamiltonian{Sys, B}, op::Operator{B, B}) where {Sys, B}
    QuantumOpticsBase.check_samebases(basis(op), ham.basis_l)
    return Hamiltonian(ham.sys, ham.data - op.data)
end

struct IncompatibleSystems <: Exception
    sys1::System
    sys2::System
end
function Base.showerror(io::IO, e::IncompatibleSystems)
    println(io, "Incompatible Hamiltonians:")
    io = IOContext(io, :compact => true)
    print(io, "  #1: ")
    show(io, "text/plain", e.sys1)
    print(io, "\n  #2: ")
    show(io, "text/plain", e.sys2)
end
function Base.:(+)(ham::Hamiltonian{Sys, B}, ham2::Hamiltonian{Sys, B}) where {Sys, B}
    ham.sys == ham2.sys || throw(IncompatibleSystems(ham.sys, ham2.sys))
    return Hamiltonian(ham.sys, Operator(ham) + Operator(ham2))
end
function Base.:(-)(ham::Hamiltonian{Sys, B}, ham2::Hamiltonian{Sys, B}) where {Sys, B}
    ham.sys == ham2.sys || throw(IncompatibleSystems(ham.sys, ham2.sys))
    return Hamiltonian(ham.sys, Operator(ham) - Operator(ham2))
end

QuantumOpticsBase.dense(ham::Hamiltonian) = Hamiltonian(ham.sys, dense(Operator(ham)))
QuantumOpticsBase.sparse(ham::Hamiltonian) = Hamiltonian(ham.sys, sparse(Operator(ham)))

Base.summary(io::IO, ham::Hamiltonian) =
    print(io, "Hamiltonian(dim=", length(ham.basis_l), "x", length(ham.basis_r), ")")
function Base.show(io::IO, mime::MIME"text/plain", ham::Hamiltonian)
    summary(io, ham)
    print(io, "\nSystem: ")
    show(io, mime, ham.sys)
    print(io, "\n")
    show(io, mime, ham.data)
end
