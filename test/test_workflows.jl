@testset "Workflows" begin
    @test begin
        l = SquareLattice(10, 10) do (x, y)
            √(x^2 + y^2) < 5
        end
        b = LatticeBasis(l)
        d = identityoperator(b)
        hx = build_operator(l, Bonds(axis=1))
        hy = build_operator(l, Bonds(axis=2))
        H = d + hx + hy
        eig = diagonalize(H)
        P = densitymatrix(eig, statistics=FermiDirac)
        all(isfinite, P.data)
    end

    @test begin
        l = HoneycombLattice(10, 10)
        H = haldane(l, 1, 1, 1)
        P = densitymatrix(diagonalize(H), statistics=FermiDirac)
        X, Y = coord_operators(basis(H))
        d = lattice_density(4π * im * P * X * (one(P) - P) * Y * P)
        rd = d .|> real
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        spin = SpinBasis(1//2)
        phasespace = l ⊗ spin
        H0 = build_hamiltonian(phasespace,
            [1 0; 0 -1] => rand(l),
            [1 im; im -1] / 2 => Bonds(axis = 1),
            [1 1; -1 -1] / 2 => Bonds(axis = 2),
            field = LandauField(0.5)
        )
        function h(t)
            x, y = coord_values(l)
            ms = @. 3 + (√(x^2 + y^2) ≤ 2) * -2
            build_hamiltonian(l, spin,
                sigmaz(spin) => ms,
                randn(l),
                [1 im; im -1] / 2 => Bonds(axis = 1),
                [1 1; -1 -1] / 2 => Bonds(axis = 2),
                field = LandauField(t)
            )
        end
        P0 = densitymatrix(diagonalize(H0), statistics=FermiDirac)
        X, Y = coord_operators(basis(H0))
        @evolution {H := h(t), P0 --> H --> P} for t in 0:0.1:10
            d = lattice_density(4π * im * P * X * (one(P) - P) * Y * P)
            ch = Currents(DensityCurrents(H, P))
            rd = d .|> real
        end
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        spin = SpinBasis(1//2)
        X, Y = coord_operators(l)
        x, y = coord_values(l)
        xy = x .* y
        p = plot(layout=4)
        plot!(p[1], xy)
        H = build_hamiltonian(l, spin,
            sigmaz(spin),
            [1 im; im -1] / 2 => Bonds(axis = 1),
            [1 1; -1 -1] / 2 => Bonds(axis = 2),
            field = LandauField(0.5)
        )
        P = densitymatrix(diagonalize(H), μ = 0.1, statistics=FermiDirac)
        @evolution k = 2 pade = true {P --> H --> PP} for t in 0:0.1:1
        end
        @evolution k = 2 {P --> H --> PP} for t in 0:0.1:1
        end

        dc = DensityCurrents(H, P)
        quiver!(p[1], dc[x.<y])
        scatter!(p[1], l)
        scatter!(p[1], l[x.<y], high_contrast=true)
        scatter!(p[1], xy[x.≥y])
        plot!(p[1], adjacency_matrix(H))
        plot!(p[1], adjacency_matrix(l, Bonds([1, 1])))
        surface!(p[2], xy)
        scatter!(p[3], SquareLattice(3, 4, 5))
        plot!(p[4], project(xy, :x))
        plot!(p[4], project(xy, p"j1"))
        mpcs = map_currents(site_distance, dc, reduce_fn=sum, sort=true)
        plot!(p[4], mpcs)
        true
    end
end
