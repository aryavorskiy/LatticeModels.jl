@testset "Currents" begin
    l = SquareLattice(4, 4)
    x, y = coordvalues(l)
    H_0 = qwz(l)
    H_1 = qwz(l, field=LandauGauge(0.1))
    dg = diagonalize(H_0)
    gs = dg[1]
    P = densitymatrix(dg, statistics=FermiDirac)
    dc = DensityCurrents(H_1, P)

    @testset "Basics" begin
        s1 = l[6]
        s2 = l[11]
        @test dc[s1, s2] == -dc[s2, s1]
        @test dc[s1, s1] ≈ 0 atol=eps()

        # Check Heisenberg equation
        t1 = s1 + BravaisTranslation(l, axis=1)
        t2 = s1 + BravaisTranslation(l, axis=2)
        t3 = s1 + BravaisTranslation(l, axis=1, dist=-1)
        t4 = s1 + BravaisTranslation(l, axis=2, dist=-1)
        state = basisstate(l, s1)
        dens_op = one(SpinBasis(1//2)) ⊗ (state ⊗ state')
        dens_dt = tr(im * (H_1 * dens_op - dens_op * H_1) * P)
        @test dc[s1, t1] + dc[s1, t2] + dc[s1, t3] + dc[s1, t4] ≈ dens_dt
        @test currentsfromto(dc, s1) ≈ dens_dt

        # Test state representations
        ground_state = dg[1]
        c1 = DensityCurrents(H_1, ground_state) |> Currents
        c2 = DensityCurrents(H_1, ground_state ⊗ ground_state') |> Currents
        @test c1 ≈ c2

        # 'materialized' currents
        l2 = SquareLattice(1, 2)
        site1, site2 = l2
        m = Currents(l2)
        m[site1, site2] = 2
        m[site1, site1] = 1
        @test m.currents == [0 2; -2 0]
        mc = Currents(dc)
        @test mc + mc == 2mc + zero(mc)
    end

    @testset "Subcurrents" begin
        bs = AdjacencyMatrix(H_1)
        m1 = Currents(dc)[x.<y]
        m2 = Currents(dc[x.<y])
        m3 = Currents(dc, bs)[x.<y]
        m4 = Currents(dc[x.<y], bs)
        @test m1 ≈ m2
        @test m1 ≈ m3
        @test m1 ≈ m4
    end

    sys = NParticles(l ⊗ SpinBasis(1//2), 2, statistics=BoseEinstein)
    H2_0 = qwz(sys)
    H2_1 = qwz(sys, field=LandauGauge(0.1))
    gs2 = diagonalize(H2_0)[1]

    @testset "Manybody" begin
        gs_curr = DensityCurrents(H_1, gs) |> Currents
        mb_gs_curr = DensityCurrents(H2_1, gs2) |> Currents
        @test 2 * gs_curr ≈ mb_gs_curr
    end

    @testset "Operator currents" begin
        gs_curr = DensityCurrents(H_1, gs) |> Currents
        new_gs_curr = LocalOperatorCurrents(H_1, gs, [1 0; 0 1]) |> Currents
        @test new_gs_curr ≈ gs_curr

        spinup_curr = LocalOperatorCurrents(H_1, gs, [1 0; 0 0]) |> Currents
        spindown_curr = LocalOperatorCurrents(H_1, gs, [0 0; 0 1]) |> Currents
        spin_curr = LocalOperatorCurrents(H_1, gs, [1 0; 0 -1]) |> Currents

        @test spinup_curr + spindown_curr ≈ gs_curr
        @test spinup_curr - spindown_curr ≈ spin_curr

        # Check Heisenberg equation
        site = l[6]
        state = basisstate(l, site)
        spin_op = Operator(SpinBasis(1//2), [1 0; 0 -1]) ⊗ (state ⊗ state')
        spin_dt = gs' * im * (H_1 * spin_op - spin_op * H_1) * gs
        @test currentsfromto(spin_curr, site) ≈ spin_dt

        mb_spin_curr = LocalOperatorCurrents(H2_1, gs2, [1 0; 0 -1]) |> Currents
        mb_spin_op = construct_operator(sys, spin_op)
        mb_spin_dt = gs2' * im * (H2_1 * mb_spin_op - mb_spin_op * H2_1) * gs2
        @test currentsfromto(mb_spin_curr, site) ≈ mb_spin_dt
        @test mb_spin_curr ≈ 2 * spin_curr
    end
end
