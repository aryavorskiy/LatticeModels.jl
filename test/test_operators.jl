@testset "Operators" begin
    @testset "Basics" begin
        l = SquareLattice(10, 10)
        H_0 = qwz(l)
        H_1 = qwz(l, field=LandauGauge(0.1))
        P = densitymatrix(H_0, statistics=FermiDirac)

        # Check Heisenberg and von Neumann equation
        site = l[11]
        state = basisstate(l, site)
        dens_op = one(SpinBasis(1//2)) ⊗ (state ⊗ state')
        dens_dt = tr(im * (H_1 * dens_op - dens_op * H_1) * P)

        dts = localdensity(-im * (H_1 * P - P * H_1))
        @test dens_dt ≈ dts[site]
    end

    @testset "Operator builder" begin
        l = SquareLattice(10, 10)
        spin = SpinBasis(1//2)
        builder = OperatorBuilder(l, spin, auto_hermitian=true, field=LandauGauge(0.1))
        builder2 = OperatorBuilder(l ⊗ spin, field=LandauGauge(0.1))
        builder3 = OperatorBuilder(ComplexF16, l, spin, auto_hermitian=true, field=LandauGauge(0.1))

        for site in l
            site_hx = site + BravaisTranslation(l, axis = 1)
            site_hy = site + BravaisTranslation(l, axis = 2)

            builder[site, site] = sigmaz(spin)
            builder[site, site_hx] = (sigmaz(spin) - im * sigmax(spin)) / 2
            builder[site, site_hy] += (sigmaz(spin) - im * sigmay(spin)) / 2

            builder2[site, site] += sigmaz(spin)
            builder2[site, site_hx] += (sigmaz(spin) - im * sigmax(spin)) / 2
            builder2[site_hx, site] += (sigmaz(spin) + im * sigmax(spin)) / 2
            builder2[site, site_hy] += (sigmaz(spin) - im * sigmay(spin)) / 2
            builder2[site_hy, site] += (sigmaz(spin) + im * sigmay(spin)) / 2

            builder3[site, site] += sigmaz(spin)
            builder3[site, site_hx] += (sigmaz(spin) - im * sigmax(spin)) / 2
            builder3[site, site_hy] += (sigmaz(spin) - im * sigmay(spin)) / 2
        end
        H = qwz(l, field=LandauGauge(0.1))
        H1 = Hamiltonian(builder)
        H2 = Hamiltonian(builder2)
        H3 = Hamiltonian(builder3)
        H4 = construct_hamiltonian(l, spin, field=LandauGauge(0.1), sigmaz(spin) => 1,
            (sigmaz(spin) - im * sigmax(spin)) / 2 => BravaisTranslation(axis = 1),
            (sigmaz(spin) - im * sigmay(spin)) / 2 => BravaisTranslation(axis = 2))
        @test H ≈ H1
        @test H ≈ H2
        @test H ≈ H3
        @test H ≈ H4
    end

    @testset "Operator builtins" begin
        l = SquareLattice(10, 10)
        X, Y = coordoperators(l)
        X2 = coordoperator(l, :x)
        X3 = diagonaloperator(l, LatticeModels.Coord(1))
        X4 = diagonaloperator(coordvalue(l, :x))
        @test X == X2
        @test X == X3
        @test X == X4

        spin = SpinBasis(1//2)
        Xs = one(spin) ⊗ X
        Xs1 = diagonaloperator(l, spin, :x)
        Xs2 = diagonaloperator(l ⊗ spin, LatticeModels.Coord(1))
        @test Xs == Xs1
        @test Xs == Xs2
    end

    @testset "Interaction" begin
        l = SquareLattice(2, 2)
        # Bose-Hubbard interaction
        sys = NParticles(l, 3, statistics=BoseEinstein)
        I1 = interaction(sys) do site1, site2
            site1 == site2 ? 1.0 : 0.0
        end
        I2 = interaction(sys, 2) do (site1, site2), (site3, site4)
            site1 == site2 == site3 == site4 ? 0.5 : 0.0
        end
        @test I1 == I2

        # Fermi-Hubbard interaction
        sys2 = NParticles(l ⊗ SpinBasis(1//2), 3, statistics=FermiDirac)
        IS1 = interaction(sys2) do site1, site2
            site1 == site2 ? 1.0 : 0.0
        end
        IS2 = interaction(sys2, 2, affect_internal = false) do (site1, site2), (site3, site4)
            site1 == site2 == site3 == site4 ? 1.0 : 0.0
        end
        IS3 = interaction(sys2, 2) do (site1, site2), is1, (site3, site4), is2
            (site1 == site2 == site3 == site4 && is1 == is2) ? 1.0 : 0.0
        end
        @test IS1 == IS2
        @test IS2 == IS3
    end

    @testset "Diagonalize" begin
        l = SquareLattice(5, 5)
        eig = diagonalize(qwz(l))
        exact_eig = eigen(Array(qwz(l).data))
        @test eig.values ≈ exact_eig.values
        @test eig.states ≈ exact_eig.vectors

        e1 = eig[1:4]
        @test e1.values ≈ eig.values[1:4]
        @test e1.states ≈ eig.states[:, 1:4]

        e2 = eig[3:6]
        e3 = union(e1, e2)
        @test e3.values ≈ eig.values[1:6]
        @test e3.states ≈ eig.states[:, 1:6]
    end

    @testset "Green's function" begin
        l = SquareLattice(10, 10)
        eig = diagonalize(qwz(ones(l)))
        G = greenfunction(eig)
        x, y = coordvalues(l)
        G2 = G[x .< 3 .&& y .< 3]
        E = 2
        site1 = l[1]
        site2 = l[2]
        slice = G[site1, site2]
        point = slice[1, 2]
        val = point(E)
        @test slice(E)[1, 2] == val
        @test G[1, 4](E) == val
        @test G2[1, 4](E) == val
        @test G(E)[1, 4] == val
        @test G(E)[site1, site2] == slice(E)
        @test G2(E)[site1, site2] == slice(E)
        δ = 0.2
        Es = eig.values
        Vs = eig.states
        ld1 = imag.(diag_reduce(tr, Operator(basis(eig), Vs * (@.(1 / (Es - E - im * δ)) .* Vs'))))
        ld2 = ldos(G, E, broaden=δ)
        @test ld2.values ≈ ld1.values
        @test dos(eig, E, broaden=δ) ≈ sum(ld1)
        @test dos(G, E, broaden=δ) ≈ sum(ld1)
    end
end