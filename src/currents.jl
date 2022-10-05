import Base: getindex
using LinearAlgebra

abstract type AbstractCurrents end
current_lambda(::AbstractCurrents) = error("not implemented")
lattice(::AbstractCurrents) = error("not implemented")

struct ChargeCurrents <: AbstractCurrents
    hamiltonian::LatticeOperator
    density::LatticeOperator
    function ChargeCurrents(ham::LatticeOperator, dens::LatticeOperator)
        @assert ham.basis == dens.basis "basis mismatch"
        new(ham, dens)
    end
end

function current_lambda(curr::ChargeCurrents)
    return (i::Int, j::Int) -> 2imag(tr(curr.density[i, j] * curr.hamiltonian[j, i]))
end
lattice(curr::ChargeCurrents) = curr.hamiltonian.basis.lattice

struct MaterializedCurrents <: AbstractCurrents
    lattice::AbstractLattice
    currents::Matrix{Float64}
    function MaterializedCurrents(l::AbstractLattice, curs::Matrix{Float64})
        @assert all(length(l) .== size(curs)) "dimension mismatch"
        new(l, curs)
    end
end

MaterializedCurrents(l::AbstractLattice) =
    MaterializedCurrents(l, zeros(length(l), length(l)))

lattice(mcurr::MaterializedCurrents) = mcurr.lattice
current_lambda(mcurr::MaterializedCurrents) = (i::Int, j::Int) -> mcurr.currents[i, j]

function materialize(curr::AbstractCurrents)
    l = lattice(curr)
    m = MaterializedCurrents(l)
    curr_fn = current_lambda(curr)
    for i in 1:length(l), j in 1:i-1
        ij_curr = curr_fn(i, j)
        m.currents[i, j] = ij_curr
        m.currents[j, i] = -ij_curr
    end
    m
end

function materialize(f::Function, curr::AbstractCurrents)
    l = lattice(curr)
    m = MaterializedCurrents(l)
    curr_fn = current_lambda(curr)
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            !(f(site1, site2)) && continue
            ij_curr = curr_fn(i, j)
            m.currents[i, j] = ij_curr
            m.currents[j, i] = -ij_curr
            j += 1
        end
        i += 1
    end
    m
end

@recipe function f(curr::AbstractCurrents)
    l = lattice(curr)
    @assert dims(l) == 2 "2D lattice expected"
    Xs = Float64[]
    Ys = Float64[]
    Qs = NTuple{2,Float64}[]
    curr_fn = current_lambda(curr)
    crd = zeros(2)
    bvs = bravais(l)
    buf = zeros(2)
    arrows_scale --> 1
    arrows_rtol --> 1e-2
    seriestype := :quiver
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            j ≥ i && continue
            ij_curr = curr_fn(i, j)::Real
            coords!(crd, l, (ij_curr > 0 ? site1 : site2), bvs, buf)
            vc = radius_vector(l, site2, site1)
            vc_n = norm(vc)
            if vc_n < abs(ij_curr * plotattributes[:arrows_scale] / plotattributes[:arrows_rtol])
                push!(Xs, crd[1])
                push!(Ys, crd[2])
                push!(Qs, Tuple(vc * (ij_curr * plotattributes[:arrows_scale] / vc_n)))
            end
            j += 1
        end
        i += 1
    end
    quiver := Qs
    Xs, Ys
end
