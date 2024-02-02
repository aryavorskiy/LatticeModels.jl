@testset "Currents" begin
    l = SquareLattice(4, 4)
    x, y = coord_values(l)
    H(B) = qwz(l, field=LandauField(B))
    dg = diagonalize(H(0))
    P = densitymatrix(dg, statistics=FermiDirac)
    dc = DensityCurrents(H(0.1), P)

    @testset "Basics" begin
        s1 = l[6]
        s2 = l[11]
        @test dc[s1, s2] == -dc[s2, s1]
        @test dc[s1, s1] ≈ 0 atol=eps()

        # Check Heisenberg equation
        t1 = s1 + SiteOffset(axis=1)
        t2 = s1 + SiteOffset(axis=2)
        t3 = s1 + SiteOffset(axis=1, dist=-1)
        t4 = s1 + SiteOffset(axis=2, dist=-1)
        dts = lattice_density(im * (H(0.1) * P - P * H(0.1)))
        @test dc[s1, t1] + dc[s1, t2] + dc[s1, t3] + dc[s1, t4] ≈ dts[s1]

        # Test state representations
        ground_state = dg[1]
        c1 = DensityCurrents(H(0.1), ground_state) |> Currents
        c2 = DensityCurrents(H(0.1), ground_state ⊗ ground_state') |> Currents
        @test c1 ≈ c2

        # 'materialized' currents
        l2 = SquareLattice(1, 2)
        m = Currents(l2)
        m[l2[1], l2[2]] = 2
        @test m.currents == [0 2; -2 0]
        mc = Currents(dc)
        @test mc + mc == 2mc + zero(mc)
    end

    @testset "Subcurrents" begin
        bs = adjacency_matrix(H(0.1))
        m1 = Currents(dc)[x.<y]
        m2 = Currents(dc[x.<y])
        m3 = Currents(dc, bs)[x.<y]
        m4 = Currents(dc[x.<y], bs)
        @test m1.currents == m2.currents
        @test m1.currents == m3.currents
        @test m1.currents == m4.currents
    end

    @testset "Manybody" begin
        sys = NParticles(l ⊗ SpinBasis(1//2), 2, statistics=BoseEinstein)
    end
end
