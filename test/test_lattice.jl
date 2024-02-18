@testset "Lattice" begin
    @testset "Interface" begin
        sql = SquareLattice(10, 20)
        @test [site_index(sql, s) for s in sql] == 1:length(sql)
        x, y = coord_values(sql)
        s_0 = SquareLattice(10, 20) do (x, y)
            x < y
        end
        s_1 = filter(sql) do (x, y)
            x < y
        end
        s_2 = sql[x.<y]
        @test s_0 == s_1
        @test s_1 == s_2

        hl = HoneycombLattice(10, 10)
        xh, yh = coord_values(hl)
        s_3 = HoneycombLattice(10, 10) do (x, y)
            x < y
        end
        s_4 = filter(hl) do (x, y)
            x < y
        end
        s_5 = hl[xh.<yh]
        @test s_3 == s_4
        @test s_4 == s_5
        @test_throws LatticeModels.IncompatibleLattices hl[x.<y]
        @test hl[!, LatticeCoord(1) => 3, j2 = 2, index = 1] ==
            LatticeModels.get_site(LatticeModels.stripparams(hl),
            LatticeModels.BravaisPointer(SA[3, 2], 1))
        xb, yb = coord_values(SquareLattice(5, 40))
        @test_throws LatticeModels.IncompatibleLattices sql[xb.<yb]
        sql2 = filter(sql) do site
            site ∉ (sql[1], sql[end])
        end
        pop!(sql)
        popfirst!(sql)
        @test sql == sql2
    end

    l = SquareLattice(10, 10)
    x, y = coord_values(l)
    xm2 = LatticeValue(l) do (x, y)
        2x
    end
    y = LatticeValue(l) do (x, y)
        y
    end
    xy = LatticeValue(l) do site
        site.coords[1] * site.coords[2]
    end
    idxs = LatticeValue(l) do site
        site_index(l, site)
    end

    @testset "LatticeValue" begin
        small_l = SquareLattice(2, 2)
        small_x, small_y = coord_values(small_l)
        @test small_x.values == [1, 1, 2, 2]
        @test small_y.values == [1, 2, 1, 2]

        x2 = coord_value(l, :x)
        x3 = siteproperty_value(l, LatticeModels.Coord(1))
        @test [idxs[s] for s in l] == 1:length(l)
        @test x == x2
        @test x == x3

        (mn, i) = findmin(xy)
        (mx, j) = findmax(xy)
        @test xy[i] == mn
        @test xy[j] == mx
        @test all(mn .≤ xy.values .≤ mx)
        imx = findall(==(mx), xy)
        @test all((site in imx || xy[site] < mx) for site in l)
    end

    @testset "Broadcast" begin
        x4 = l .|> (site -> site.coords[1])
        @test x == x4
        @test x .* y == xy
        @test x .* 2 == xm2
        @test 2 .* x == xm2
        @test x .|> (x -> 2x) == xm2
        y .= x .* y
        @test y == xy
        @test_throws ArgumentError x .* ones(100)
        l′ = SquareLattice(5, 20)
        x′ = LatticeValue(l′) do (x, y)
            x
        end
        @test_throws LatticeModels.IncompatibleLattices x .* x′
    end

    @testset "Indexing" begin
        z = zeros(l)
        z2 = zero(z)
        z .= 1
        z2[x.≥y] .+= 1
        z2[x.<y] = ones(l)
        @test z == ones(l)
        @test z2 == ones(l)
        @test z[!, x=1, y=1] == 1
        z3 = ones(l)
        for i in 2:10
            z3[x=i] .= i
        end
        @test z3 == x
    end
end

@testset "Bonds" begin
    @testset "Constructor" begin
        @test BravaisTranslation(axis=1) == BravaisTranslation([1])
        @test BravaisTranslation(axis=1) != BravaisTranslation([1, 0])
        @test_throws ArgumentError BravaisTranslation(2 => 2, [0, 0])
    end
    @testset "Hopping matching" begin
        l1 = SquareLattice(6, 5)
        l2 = SquareLattice(6, 5, boundaries = ([6, 0] => true))
        hx = BravaisTranslation(l1, axis=1)
        hxmy = BravaisTranslation(l2, [1, -1])
        for site in l1
            ucx, ucy = site.latcoords
            dst_dx = LatticeModels.resolve_site(l1, site + hx)
            dst_dxmy = LatticeModels.resolve_site(l2, site + hxmy)
            @test (dst_dx === nothing) == (ucx == 6)
            @test (dst_dxmy === nothing) == (ucy == 1)
        end
    end
    @testset "Bonds" begin
        l = SquareLattice(2, 2)
        ls1, _, ls3, ls4 = l
        bs = adjacency_matrix(l, BravaisTranslation(axis=1), BravaisTranslation(axis=2))
        bs1 = union(adjacency_matrix(l, BravaisTranslation(axis=1)), adjacency_matrix(l, BravaisTranslation(axis=2)))
        @test bs.mat == bs1.mat
        @test bs[ls1, ls3]
        @test !bs[ls1, ls4]
    end
end
