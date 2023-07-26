using LatticeModels, Plots
l = SquareLattice(16, 16)
spin = SpinBasis(1//2)
H = tightbinding_hamiltonian(l, tn = 1, tnn = 0.25, tnnn = -1)
eig = diagonalize(H)
P = densitymatrix(eig, μ = 0, T = 1)
dens = site_density(P)
plot(dens)

# m = ones(l)
# m[x = 1, y = 2] = -1
# m[x = 2, y = 1] = -1
# d = diagonaloperator(m)

# spin = SpinBasis(1//2)
# H = sigmaz(spin) ⊗ d
# H = tightbinding_hamiltonian(l, sigmaz(spin) => m)

# b = LatticeBasis(l)
# op = zero(b)
# for (site1, site2) in bond_pairs(l, Bonds(axis = 1))
#     i = site_index(l, site1)
#     j = site_index(l, site2)
#     op.data[i, j] = t(i, j)
# end
# opb = OperatorBuilder(l, spin)
# op = zero(l ⊗ spin)
# for site1 in l, site2 in l
#     # e1 = Bonds(axis = 1)
#     if site1 == site2 + basis_vector(l, 1)
#         op[site1, site2] = f(site1, site2) * sigmaz(spin)
#     end
# end
# op = to_operator(opb)
# for (site1, site2) in shifted(l, )
#     op[site1, site2]
# end

# amplitudes = ones(l)
# site1 = l[2]
# site2 = l[4]
# amplitudes[site] = 0

#==
Define hopping coefficients
==#
# l = SquareLattice(16, 16)
# sample = Sample(l, N=3, statistics=BoseEinstein)
# pbc = BoundaryConditions(:x => true, :y => false)

# site1 = l[x = 2, y = 3]
# site2 = l[23]
# hop = tightbinding_hamiltonian(l,
#     field = LandauField(0.5),
#     0.5 => Bonds(axis = 1),
#     site1 => site2,
#     boundaries = pbc)
# add_tb!(hop, 0.5 => Bonds([1, 1]))
# U = 10
# inter = interaction(sample) do site1, site2
#     1 / norm(site1.coords - site2.coords)
# end
# H = hop + inter

# ...
# field = MagneticField() do x, y
#     [B * y, 0]
# end
