import LatticeModels: stripmeta, bravaispointer_to_site, BravaisPointer, IncompatibleLattices

@testset "Lattice" begin
    @testset "Creation" begin
        sl = SquareLattice(3, 3, offset=[0.234, 1], rotate=[0 1; 1 0])
        @test sl[!, j1=1, j2=2].coords == [2.234, 2]
        su = SquareLattice(3, 3, offset=[-100, 0], postoffset=:centeralign) do (x, y)
            (x + 100) ^ 2 + y ^ 2 < 10
        end
        su2 = SquareLattice(3, 3, offset=:centeralign) do (x, y)
            (x - 0.5) ^ 2 + (y - 0.5) ^ 2 < 10
        end
        @test su == su2
        uc = UnitCell([1 0; 0 1], [[0.2, 0] [-0.3, 0]], offset=:center)
        @test uc.basissites == [0.25 -0.25; 0 0]
        @test_throws ArgumentError UnitCell([1 0; 0 1], [[0.2, 0] [-0.3, 0]], offset=:AAA)

        sl = SquareLattice(-1:1, -1:1)
        sl2 = SquareLattice{2}(Square(h = 1))
        @test stripmeta(sl) == stripmeta(sl2)

        circle_l = TriangularLattice(-10:10, -10:10) do (x, y)
            x^2 + y^2 < 50 && (abs(x) > 3 || abs(y) > 3)
        end
        circle_l2 = TriangularLattice(Circle(√50), !Square(h = 3))
        @test stripmeta(circle_l) == stripmeta(circle_l2)

        complexsample = SquareLattice{2}(
            Circle(10), Circle(10, [20, 0]), Circle(10, [10, 10√3]),
            !Circle(5), !Circle(5, [20, 0]), !Circle(5, [10, 10√3]),
            Rectangle(-5..5, -14..(-12)), Rectangle(15..25, -14..(-12)),
            Path([-12, 32], [32, 32]), sites=10000
        )
        removedangling!(complexsample, maxdepth=2)
        addshapes!(complexsample, SiteAt([0, 0]))
        @test length(complexsample) ≈ 10000 rtol=0.03

        @test_throws ArgumentError SquareLattice(Circle(10))
    end
    @testset "Indexing" begin
        sl = SquareLattice(10, 20)
        x, y = coordvalues(sl)
        s_0 = SquareLattice(10, 20) do (x, y)
            x < y
        end
        s_1 = filter(sl) do (x, y)
            x < y
        end
        s_2 = sl[x.<y]
        filter!(site -> site.coords[1] < site.coords[2], sl)
        @test sl == s_0
        @test sl == s_1
        @test sl == s_2

        x, y = coordvalues(SquareLattice(10, 20))
        hl = HoneycombLattice(10, 10)
        xh = LatticeValue(hl, :x)
        j1h = LatticeValue(hl, :j1)
        h_1 = hl[x = 4..Inf, j1 = 1..3]
        h_2 = hl[@. xh ≥ 4 && 1 ≤ j1h ≤ 3]
        @test h_1 == h_2
        @test_throws IncompatibleLattices hl[x.<y]

        @test hl[!, LatticeCoord(1) => 3, j2 = 2, index = 1] ==
            bravaispointer_to_site(stripmeta(hl), BravaisPointer(SA[3, 2], 1))
        @test_throws BoundsError hl[!, LatticeCoord(1) => 100]

        sl = SquareLattice(10, 20)
        xb, yb = coordvalues(SquareLattice(5, 40))
        @test_throws IncompatibleLattices sl[xb.<yb]
        s = sl[!, x = 1, y = 2]
        sl2 = filter(sl) do site
            site ∉ (sl[1], sl[end], s)
        end
        pop!(sl)
        popfirst!(sl)
        delete!(sl, s)
        @test sl == sl2

        lwh = SquareLattice(10, 10) do (x, y)
            x^2 + y^2 > 25
        end
        lwh2 = addlookuptable(lwh)
        gen_lwh = GenericLattice(lwh)
        for (i, site) in enumerate(lwh)
            @test site_index(lwh, site) == i
            @test site_index(lwh2, site) == i
            @test site_index(gen_lwh, site) == i
        end
    end
    @testset "GenericLattice" begin
        small_l = SquareLattice(2, 2)
        ls1, ls2, ls3, ls4 = small_l
        gl1 = GenericLattice(small_l)
        gl2 = GenericLattice{typeof(ls1)}()
        push!(gl2, ls3, ls4, ls2, ls1)
        @test gl1 == gl2

        gl3 = GenericLattice{2}()
        union!(gl3, small_l)
        gl4 = GenericLattice([GenericSite(1, 1), GenericSite(1, 2), GenericSite(2, 1), GenericSite(2, 2)])
        gl5 = GenericLattice([(1, 1), (1, 2), (2, 1), [2, 2]])
        @test gl3 == gl4
        @test gl3 == gl5
        @test_throws ArgumentError GenericLattice([(1, 2), (3, 4, 5)])
        @test_throws ArgumentError GenericLattice([(1, 2), (1, 1)])
    end
end

@testset "LatticeValue" begin
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

    @testset "Builtins" begin
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
        @test x * 2 == xm2
        @test 2 .* x == xm2
        @test 2x == xm2
        @test x - x == zero(x)
        @test -x + x == zero(x)
        @test x .|> (x -> 2x) == xm2
        y .= x .* y
        @test y == xy
        @test_throws ArgumentError x .* ones(100)
        x′ = LatticeValue(SquareLattice(5, 20), :x)
        @test_throws IncompatibleLattices x .* x′
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

        sites = [l[!, x = 2, y = i] for i in 1:10]
        @test z2[x = 2] == z2[sites]

        zz = zeros(l)
        for i in 1:10
            zz[!, x=i, y=i] =1
        end
        zzvals = zeros(100)
        zzvals[1:11:end] .= 1
        @test zz.values == zzvals
        @test_throws ArgumentError zz[x=1, y=1] = 1
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
        am = AdjacencyMatrix(l, BravaisTranslation(axis=1), BravaisTranslation(axis=2))
        am1 = union(AdjacencyMatrix(l, Bravais[1, 0]), AdjacencyMatrix(l, Bravais[0, 1]))
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
    @testset "Bravais translations" begin
        l = HoneycombLattice(-2:2, -2:2)
        tr = Translation(l, [0, 2√3/3])
        site1 = l[!, j1 = 1, j2 = 1, index = 1]
        site2 = site1 + tr
        @test site2 == l[!, j1 = 0, j2 = 2, index = 2]
        @test only(tr.translations) == BravaisTranslation(1 => 2, [-1, 1])

        bsm = LatticeModels.BravaisSiteMapping(
            BravaisTranslation(axis=1),
            BravaisTranslation(axis=2),
            BravaisTranslation([1, -1]))
        bsm2 = LatticeModels.adapt_bonds(bsm, l)
        @test bsm == union(Bravais[1], Bravais[0, 1], Bravais[1, -1])
        @test bsm == union(Bravais[0, 1], bsm)
        @test bsm == union(bsm, bsm)
        @test bsm2 == union(bsm, bsm2)

        t1 = BravaisTranslation(l, axis=1)
        t2 = BravaisTranslation(l, axis=2)
        t3 = BravaisTranslation(HoneycombLattice(2, 2), [1, -1])
        @test union(t1, t2, bsm) == bsm2
        @test_throws IncompatibleLattices union(t1, t2, t3)
    end
end
