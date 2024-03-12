@testset "Lattice" begin
    @testset "Interface" begin
        sql = SquareLattice(10, 20)
        @test [site_index(sql, s) for s in sql] == 1:length(sql)
        x, y = coordvalues(sql)
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
        xh, yh = coordvalues(hl)
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
        xb, yb = coordvalues(SquareLattice(5, 40))
        @test_throws LatticeModels.IncompatibleLattices sql[xb.<yb]
        s = sql[!, x = 1, y = 2]
        sql2 = filter(sql) do site
            site ∉ (sql[1], sql[end], s)
        end
        pop!(sql)
        popfirst!(sql)
        delete!(sql, s)
        @test sql == sql2

        sql3 = sql[x = -Inf..5.2]
        filter!(site -> site.coords[1] < 5.2, sql)
        @test sql == sql3

        small_l = SquareLattice(2, 2)
        ls1, ls2, ls3, ls4 = small_l
        gl1 = GenericLattice(small_l)
        gl2 = GenericLattice{typeof(ls1)}()
        push!(gl2, ls3, ls4, ls2, ls1)
        @test gl1 == gl2
    end

    @testset "Creation" begin
        sl = SquareLattice(-1:1, -1:1)
        sl2 = SquareLattice{2}(Square(h = 1))
        @test LatticeModels.stripparams(sl) == LatticeModels.stripparams(sl2)

        circle_l = TriangularLattice(-10:10, -10:10) do site
            return site.coords[1]^2 + site.coords[2]^2 ≤ 50 &&
                (abs(site.coords[1]) > 3 || abs(site.coords[2]) > 3)
        end
        circle_l2 = TriangularLattice(Circle(√50), !Square(h = 3))
        @test LatticeModels.stripparams(circle_l) == LatticeModels.stripparams(circle_l2)

        complexsample = SquareLattice{2}(
            Circle(10), Circle(10, [20, 0]), Circle(10, [10, 5√3]),
            !Circle(5), !Circle(5, [20, 0]), !Circle(5, [10, 5√3]),
            Path([-12, -12], [-12, 32]), sites=10000
        )
        removedangling!(complexsample, maxdepth=2)
        @test length(complexsample) ≈ 10000 rtol=0.03

        @test_throws ArgumentError SquareLattice(Circle(10))
    end

    l = SquareLattice(10, 10)
    x, y = coordvalues(l)
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
        small_x, small_y = coordvalues(small_l)
        @test small_x.values == [1, 1, 2, 2]
        @test small_y.values == [1, 2, 1, 2]

        x2 = coordvalue(l, :x)
        x3 = LatticeValue(l, LatticeModels.Coord(1))
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

        gl = GenericLattice{2}()
        push!(gl, (1, 2), (2, 1), (3, 4))
        gx = coordvalue(gl, :y)
        @test gx.values == [2, 1, 4]
        @test_throws ArgumentError coordvalue(gl, :j1)
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
    @testset "AdjacencyMatrix" begin
        l = SquareLattice(2, 2)
        ls1, ls2, ls3, ls4 = l
        am = adjacencymatrix(l, BravaisTranslation(axis=1), BravaisTranslation(axis=2))
        am1 = union(adjacencymatrix(l, Bravais[1, 0]), adjacencymatrix(l, Bravais[0, 1]))
        @test am.mat == am1.mat
        @test am[ls1, ls3]
        @test !am[ls1, ls4]
        @test Set(adjacentsites(am, ls1)) == Set([ls2, ls3])

        am2 = AdjacencyMatrix(l)
        am2[ls1, ls2] = true
        am2[ls1, ls3] = true
        am2[ls2, ls4] = true
        am2[ls3, ls4] = true
        @test am.mat == am2.mat

        ps = Tuple{Int, Int}[]
        for (s1, s2) in am2
            push!(ps, (s1.index, s2.index))
        end
        @test Set(ps) == Set([(1, 2), (1, 3), (2, 4), (3, 4)])
    end

    @testset "Abstract bonds" begin
        l = HoneycombLattice(-2:2, -2:2)
        tr = Translation(l, [0, 2√3/3])
        site1 = l[!, j1 = 1, j2 = 1, index = 1]
        site2 = site1 + tr
        @test site2 == l[!, j1 = 0, j2 = 2, index = 2]
        @test only(tr.translations) == BravaisTranslation(1 => 2, [-1, 1])

        bsm = LatticeModels.BravaisSiteMapping(l,
            BravaisTranslation(axis=1),
            BravaisTranslation(axis=2),
            BravaisTranslation([1, -1]))
        @test bsm.translations == union(Bravais[1], Bravais[0, 1], Bravais[1, -1]).translations
    end
end
