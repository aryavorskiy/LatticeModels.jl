@testset "Workflows" begin
    begin
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
        @test all(isfinite, P.data)
    end

    @test begin
        l = GrapheneRibbon(10, 5, rotate=pi/2)
        @test_throws ArgumentError haldane(SquareLattice(10, 10), 1, 1, 1)
        H = haldane(l, 1, 1, 1)
        P = densitymatrix(diagonalize(H), statistics=BoseEinstein)
        X, Y = coordoperators(basis(H))
        d = localdensity(4π * im * P * X * (one(P) - P) * Y * P)
        rd = d .|> real
        true
    end

    @test begin
        l = SquareLattice(10, 10, boundaries = (:axis1 => true))
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
        p = plot(layout=5)
        plot!(p[1], xy, st=:shape, shape=:circle, markerscale=true)
        plot!(p[1], xy, st=:shape, shape=:polygon, markerscale=true)
        H = construct_hamiltonian(l, spin,
            sigmaz(spin),
            [1 im; im -1] / 2 => BravaisTranslation(axis = 1),
            [1 1; -1 -1] / 2 => BravaisTranslation(axis = 2),
            field = LandauGauge(0.5)
        )
        P = densitymatrix(diagonalize(H), N = 3, T = 1, statistics=FermiDirac)

        dc = DensityCurrents(H, P)
        plot!(p[1], dc[x.<y])
        plot!(p[1], Currents(dc)[x.<y])
        scatter!(p[1], l, shownumbers=true)
        scatter!(p[1], l[2])
        scatter!(p[1], l[x.<y], :high_contrast)
        scatter!(p[1], xy[x.≥y])
        plot!(p[1], AdjacencyMatrix(H))
        plot!(p[1], AdjacencyMatrix(BravaisTranslation(l, [1, 1])))
        @test_throws ErrorException surface!(p[2], xy)
        plot!(p[2], l)
        plot!(p[2], Circle(10), Rectangle(10 .. 20, -5 .. 5), Hexagon(10, [-10, 0]),
            SiteAt([0, 0]), Path([-10, -15], [10, 15]))
        plot!(p[3], SquareLattice(3, 4, 5))
        plot!(p[4], axis=:x)
        plot!(p[5], UnitCell([[1, 0] [0.1, 1]], [[0, 0.1] [-0.2, 0]]))
        true
    end
end
