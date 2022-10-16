import Base: getindex
using LinearAlgebra

abstract type AbstractCurrents end
current_lambda(::AbstractCurrents) = error("not implemented")
lattice(::AbstractCurrents) = error("not implemented")

struct DensityCurrents <: AbstractCurrents
    hamiltonian::LatticeOperator
    density::LatticeOperator
    function DensityCurrents(ham::LatticeOperator, dens::LatticeOperator)
        ham.basis != dens.basis && error("basis mismatch")
        new(ham, dens)
    end
end

current_lambda(curr::DensityCurrents) =
    (i::Int, j::Int) -> 2imag(tr(curr.density[i, j] * curr.hamiltonian[j, i]))
lattice(curr::DensityCurrents) = curr.hamiltonian.basis.lattice

struct SubCurrents{CT<:AbstractCurrents} <: AbstractCurrents
    parent_currents::CT
    lattice::Lattice
    indices::Vector{Int}
end
function SubCurrents(parent_currents::CT, indices::Vector{Int}) where {CT}
    l = lattice(parent_currents)
    m = zeros(Bool, length(l))
    m[indices] .= true
    new_mask = zero(l.mask)
    new_mask[l.mask] = m
    SubCurrents{CT}(parent_currents,
        Lattice(lattice_type(l), size(l), bravais(l), new_mask), indices)
end

function current_lambda(scurr::SubCurrents)
    in_fn = current_lambda(scurr.parent_currents)
    (i::Int, j::Int) -> in_fn(scurr.indices[i], scurr.indices[j])
end
lattice(scurr::SubCurrents) = scurr.lattice

struct MaterializedCurrents <: AbstractCurrents
    lattice::Lattice
    currents::Matrix{Float64}
    function MaterializedCurrents(l::Lattice, curs::Matrix{Float64})
        !all(length(l) .== size(curs)) && error("dimension mismatch")
        new(l, curs)
    end
end

MaterializedCurrents(l::Lattice) =
    MaterializedCurrents(l, zeros(length(l), length(l)))

lattice(mcurr::MaterializedCurrents) = mcurr.lattice
current_lambda(mcurr::MaterializedCurrents) = (i::Int, j::Int) -> mcurr.currents[i, j]

function getindex(curr::AbstractCurrents, lvm::LatticeValue{Bool})
    lattice(curr) != lattice(lvm) && error("lattice mismatch")
    indices = findall(lvm.values)
    SubCurrents(curr, indices)
end

function getindex(curr::MaterializedCurrents, lvm::LatticeValue{Bool})
    lattice(curr) != lattice(lvm) && error("lattice mismatch")
    indices = findall(lvm.values)
    MaterializedCurrents(curr.lattice[lvm], curr.currents[indices, indices])
end

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
            if i == j || !(f(site1, site2))
                j += 1
                continue
            end
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
    dims(l) != 2 && error("2D lattice expected")
    Xs = Float64[]
    Ys = Float64[]
    Qs = NTuple{2,Float64}[]
    curr_fn = current_lambda(curr)
    arrows_scale --> 1
    arrows_rtol --> 1e-2
    seriestype := :quiver
    i = 1
    for site1 in l
        j = 1
        for site2 in l
            j ≥ i && continue
            ij_curr = curr_fn(i, j)::Real
            crd = coords(l, (ij_curr > 0 ? site1 : site2))
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
