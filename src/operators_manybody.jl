const AbstractLatticeBasis = Union{LatticeBasis, CompositeLatticeBasis}
# Override this function with more performant one
function QuantumOpticsBase.manybodyoperator_1(basis::ManyBodyBasis,
        op::Operator{<:AbstractLatticeBasis, <:AbstractLatticeBasis, <:SparseMatrixCSC})
    N = length(basis)
    builder = SparseMatrixBuilder{ComplexF64}(N, N)
    M = op.data
    @inbounds for colindex = 1:M.n
        for i=M.colptr[colindex]:M.colptr[colindex+1]-1
            row = M.rowval[i]
            value = M.nzval[i]
            for m=1:N, n=1:N
                C = my_coefficient(basis.occupations[m], basis.occupations[n], row, colindex)
                if C != 0.
                    increment!(builder, C * value, m, n)
                end
            end
        end
    end
    return SparseOperator(basis, to_matrix(builder))
end

# And also this
Base.@propagate_inbounds function my_coefficient(occ_m, occ_n, at_indices, a_indices)
    any(==(0), (occ_m[m] for m in at_indices)) && return 0.
    any(==(0), (occ_n[n] for n in a_indices)) && return 0.
    C = prod(√, @view occ_m[at_indices]) * prod(√, @view occ_n[a_indices])
    for i in 1:length(occ_m)
        vm = occ_m[i]
        vn = occ_n[i]
        i in at_indices && (vm -= 1)
        i in a_indices && (vn -= 1)
        vm != vn && return zero(C)
    end
    return C
end

function fills_iterator end

function midsite_interaction(f::Function, sample::Sample)
    for site1 in lattice(sample)
        for site2 in lattice(sample)

        end
    end
end
