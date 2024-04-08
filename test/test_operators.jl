import LatticeModels: Sample

@testset "Operators" begin
    @testset "Basics" begin
        l = SquareLattice(10, 10)
        x = coordvalue(l, :x)
        psi = ketstate(@. exp(-x^2/2 + im * x))
        psi_t = brastate(@. exp(-x^2/2 - im * x))
        @test LatticeValue(psi) == @. exp(-x^2/2 + im * x)
        @test LatticeValue(psi_t) == @. exp(-x^2/2 - im * x)
        @test_throws ArgumentError LatticeValue(basisstate(SpinBasis(1//2), 1) ⊗ psi)

        spin = SpinBasis(1//2)
        s = basisstate(spin, 1)
        psic = s ⊗ psi
        psic2 = @.(exp(-x^2/2 + im * x)) ⊗ s
        @test psic.data ≈ psic2.data
        psic_t = @.(exp(-x^2/2 - im * x)) ⊗ s'
        @test psic_t * psic ≈ psi_t * psi

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
        @test_throws ArgumentError localdensity(basisstate(spin, 1))

        # ptrace
        @test localdensity(ptrace(P, :internal)) ≈ localdensity(P)
        @test tr(P) ≈ tr(ptrace(P, :lattice))
        @test_throws ArgumentError ptrace(P, :invalid_subsp)
    end

    @testset "Operator builder" begin
        l = SquareLattice(10, 10)
        spin = SpinBasis(1//2)
        builder = OperatorBuilder(l, spin, auto_hermitian=true, field=LandauGauge(0.1))
        builder2 = FastOperatorBuilder(l ⊗ spin, field=LandauGauge(0.1))
        builder3 = OperatorBuilder(ComplexF16, l, spin, auto_hermitian=true, field=LandauGauge(0.1))

        site = l[10]
        builder[site, site] = 1234  # This should be overwritten
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
        @test H - Operator(H1) ≈ Operator(H2) - H4 atol=1e-10
        @test H - H1 ≈ -H2 + H4 atol=1e-10

        site1, site2 = l[10], l[25]
        sigmaydi10 = sigmay(spin) ⊗ diagonaloperator(l .== Ref(site1))
        @test 2sigmaydi10 == construct_operator(l, spin, sigmaydi10, sigmay(spin) => site1)
        onedi25 = diagonaloperator(l .== Ref(site2))
        @test 3 * one(spin) ⊗ onedi25 == construct_operator(l, spin, onedi25, site2, [1 0; 0 1] => site2)
        @test one(spin) ⊗ one(LatticeBasis(l)) == construct_operator(l, spin, [1 0; 0 1])
        @test transition(l, site1, site2) == construct_operator(l, site1 => site2)
        @test_throws ArgumentError construct_operator(l, spin, one(spin) ⊗ one(spin))
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

        @test bosehubbard(l, 2, U = 10) ≈
            bosehubbard(l, 2, t1=0.1, U = 0) + bosehubbard(l, 2, t1=0.9, U = 10)

        H = bosehubbard(l, 1)
        d1 = localdensity(groundstate(H))
        H2 = bosehubbard(l, 2)
        d2 = localdensity(groundstate(H2))
        @test d1.values * 2 ≈ d2.values
    end

    @testset "Diagonalize" begin
        l = SquareLattice(5, 5)
        eig = diagonalize(qwz(l))
        exact_eig = eigen(Array(qwz(l).data))
        @test eig.values ≈ exact_eig.values
        @test eig.states ≈ exact_eig.vectors

        val = eig.values[1]
        @test eig[value=val] == eig[1]

        e1 = eig[1:4]
        @test e1.values ≈ eig.values[1:4]
        @test e1.states ≈ eig.states[:, 1:4]

        e2 = eig[3:6]
        e3 = union(e1, e2)
        @test e3.values ≈ eig.values[1:6]
        @test e3.states ≈ eig.states[:, 1:6]

        d1 = diagonalize(qwz(l, field=PointFlux(0.1, (3.5, 3.5))))
        d2 = diagonalize(qwz(l, field=PointFlux(0.1, (3.5, 3.5), gauge=:singular)))
        @test d1.values ≈ d2.values

        @test_throws ErrorException diagonalize(qwz(l), :invalid_routine)
        eigk = diagonalize(qwz(l), :krylovkit)
        @test abs(eigk[1]' * eig[1]) ≈ 1

        spin = SpinBasis(1//2)
        @test System(l ⊗ spin, N = 3) == LatticeModels.FixedN(Sample(l, spin), 3)
        @test System(l, mu = 0) == LatticeModels.FixedMu(Sample(l), 0)
        @test System(l, N = 3) != System(l, mu = 0)
        @test_throws ArgumentError System(l, N = 3, mu = 0)

        sys = System(l, N = 3, T = 1, statistics=FermiDirac)
        H = tightbinding_hamiltonian(sys)
        dm = densitymatrix(H)
        @test tr(dm) ≈ 3

        H2 = tightbinding_hamiltonian(l)
        P = densitymatrix(H2, T = 0)
        psi = groundstate(H2)
        @test P ≈ psi ⊗ psi'
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
        @test ld2 ≈ ld1
        @test dos(eig, E, broaden=δ) ≈ sum(ld1)
        @test dos(G, E, broaden=δ) ≈ sum(ld1)

        l2 = SquareLattice(4, 4)
        U = 0
        G1 = greenfunction(diagonalize(tightbinding_hamiltonian(l2)))
        E0, gs = findgroundstate(bosehubbard(l2, 2, U = U))
        Hp = bosehubbard(l2, 3, U = U)
        Hm = bosehubbard(l2, 1, U = U)
        G2 = greenfunction(gs, Hp, Hm, E₀ = E0)
        @test G1(E + 0.1im).values ≈ G2(E + 0.1im).values
    end
end
