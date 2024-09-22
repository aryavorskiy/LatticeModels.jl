using PrecompileTools

if VERSION >= v"1.8"
    @compile_workload begin
        l1 = SquareLattice(3, 3)
        l2 = GrapheneRibbon(3, 3)
        for l in (l1, l2)
            m = LatticeValue(l, :x)
            H1 = tightbinding_hamiltonian(l, m, field = LandauGauge(0.1))
            d = diagonalize(dense(H1))
            dens = localdensity(d[1])
            rho = densitymatrix(d, N = 4, info=false)
            rho2 = densitymatrix(d, mu = 3, T = 1, info=false)
            curr = Currents(DensityCurrents(H1, rho))
            for state in Evolution(H1, rho, showprogress=false, timedomain=0:0.1:0.2)
                P, H, t = state
            end
        end
    end
end
