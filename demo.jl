using LatticeModels
l = SquareLattice(6, 6)
spin = SpinBasis(1//2)
sample = Sample(l, N=2, statistics=BoseEinstein)

ms = 1 .+ rand(l) .* 0.1
# H = qwz(sample, ms, field = LandauField)

H = tight_binding(sample,
    [1;;] => Bonds(axis = 1),
    [1;;] => Bonds(axis = 2),
    field = LandauField(0.5))
U = 10
# inter = interaction(sample) do site1, site2
#     site1 == site2 ? U : zero(U)
# end
# H = hop + inter

# ...
