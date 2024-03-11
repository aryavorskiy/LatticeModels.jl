@testset "Workflows" begin
    @test begin
        l = SquareLattice(10, 10) do (x, y)
            √(x^2 + y^2) < 5
        end
        b = LatticeBasis(l)
        d = identityoperator(b)
        hx = construct_operator(l, BravaisTranslation(axis=1))
        hy = construct_operator(l, BravaisTranslation(axis=2))
        H = d + hx + hy
        eig = diagonalize(H)
        P = densitymatrix(eig, statistics=FermiDirac, N = 3)
        all(isfinite, P.data)
    end

    @test begin
        l = HoneycombLattice(10, 10, boundaries = (:axis1 => true))
        H = haldane(l, 1, 1, 1)
        P = densitymatrix(diagonalize(H), statistics=BoseEinstein)
        X, Y = coordoperators(basis(H))
        d = localdensity(4π * im * P * X * (one(P) - P) * Y * P)
        rd = d .|> real
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        spin = SpinBasis(1//2)
        phasespace = l ⊗ spin
        H0 = construct_hamiltonian(phasespace,
            [1 0; 0 -1] => rand(l),
            [1 im; im -1] / 2 => BravaisTranslation(axis = 1),
            [1 1; -1 -1] / 2 => BravaisTranslation(axis = 2),
            field = LandauGauge(0.5)
        )
        function h(t)
            x, y = coordvalues(l)
            ms = @. 3 + (√(x^2 + y^2) ≤ 2) * -2
            construct_hamiltonian(l, spin,
                sigmaz(spin) => ms,
                randn(l),
                [1 im; im -1] / 2 => BravaisTranslation(axis = 1),
                [1 1; -1 -1] / 2 => BravaisTranslation(axis = 2),
                field = LandauGauge(t)
            )
        end
        P0 = densitymatrix(diagonalize(H0), μ = 3)
        X, Y = coordoperators(basis(H0))
        evol = Evolution(h, P = P0)
        densities = TimeSequence{LatticeValue}()
        for state in evol(0:0.1:10)
            P, H, t = state
            d = localdensity(4π * im * P * X * (one(P) - P) * Y * P)
            ch = Currents(DensityCurrents(H, P))
            densities[t] = d .|> real
        end
        true
    end

    @test begin
        b = FunctionBoundary([0, 10]) do site
            exp(im * site.x)
        end
        l = SquareLattice(10, 10, boundaries=(:axis1 => true, b))
        spin = SpinBasis(1//2)
        X, Y = coordoperators(l)
        x, y = coordvalues(l)
        xy = x .* y
        p = plot(layout=4)
        plot!(p[1], xy)
        H = construct_hamiltonian(l, spin,
            sigmaz(spin),
            [1 im; im -1] / 2 => BravaisTranslation(axis = 1),
            [1 1; -1 -1] / 2 => BravaisTranslation(axis = 2),
            field = LandauGauge(0.5)
        )
        P = densitymatrix(diagonalize(H), N = 3, T = 1, statistics=FermiDirac)

        dc = DensityCurrents(H, P)
        quiver!(p[1], dc[x.<y])
        scatter!(p[1], l)
        scatter!(p[1], l[x.<y], high_contrast=true)
        scatter!(p[1], xy[x.≥y])
        plot!(p[1], adjacencymatrix(H))
        plot!(p[1], adjacencymatrix(BravaisTranslation(l, [1, 1])))
        surface!(p[2], xy)
        scatter!(p[3], SquareLattice(3, 4, 5))
        plot!(p[4], project(xy, :x))
        plot!(p[4], project(xy, :j1))
        mpcs = mapgroup_currents(sitedistance, sum, dc, sortresults=true)
        plot!(p[4], mpcs)
        true
    end
end