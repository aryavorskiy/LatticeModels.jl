using Test, LinearAlgebra, StaticArrays, Plots
using LatticeModels

@testset "Workflows" begin
    @test begin
        l = SquareLattice(10, 10) do site, (x, y)
            √(x^2 + y^2) < 5
        end
        b = Basis(l, 1)
        d = LatticeOperator(I, b)
        hx = hopping_operator(l, hopping(axis=1))
        hy = hopping_operator(l, hopping(axis=2))
        H = d + hx + hy
        sp = spectrum(H)
        P = filled_projector(sp)
        d = diag_aggregate(tr, P)
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        H = @hamiltonian begin
            lattice := l
            field := LandauField(3)
            @diag [1 0; 0 -1]
            @hop axis = 1 [1 1; 1 -1] / √2
            @hop axis = 2 [1 -im; im -1] / √2
        end
        P = filled_projector(spectrum(H))
        X, Y = coord_operators(Basis(l, 2))
        d = diag_aggregate(tr, 4π * im * P * X * (I - P) * Y * P)
        rd = d .|> real
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        H0 = @hamiltonian begin
            lattice := l
            field := LandauField(0.5)
            @diag [1 0; 0 -1]
            @hop axis = 1 [1 im; im -1] / 2
            @hop axis = 2 [1 1; -1 -1] / 2
        end
        function h(t)
            x, y = coord_values(l)
            ms = @. 3 + (√(x^2 + y^2) ≤ 2) * -2
            @hamiltonian begin
                lattice := l
                field := LandauField(t)
                @diag ms ⊗ [1 0; 0 -1]
                @hop axis = 1 [1 im; im -1] / 2
                @hop axis = 2 [1 1; -1 -1] / 2
            end
        end
        P0 = filled_projector(spectrum(H0))
        X, Y = coord_operators(Basis(l, 2))
        @evolution {H := h(t), P0 --> H --> P} for t in 0:0.1:10
            d = diag_aggregate(tr, 4π * im * P * X * (I - P) * Y * P)
            ch = materialize(DensityCurrents(H, P))
            rd = d .|> real
        end
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        X, Y = coord_operators(l, 2)
        x, y = coord_values(l)
        xy = x .* y
        p = plot(layout=3)
        heatmap!(p[1], xy)
        H = @hamiltonian begin
            lattice := l
            field := LandauField(0.5)
            @diag [1 0; 0 -1]
            @hop axis = 1 [1 im; im -1] / 2
            @hop axis = 2 [1 1; -1 -1] / 2
        end
        P = filled_projector(spectrum(H), 0.1)
        @evolution k = 2 {P --> H --> PP} for t in 0:0.1:1
        end

        quiver!(p[1], DensityCurrents(H, P)[x.<y])
        scatter!(p[1], l)
        scatter!(p[1], l[x.<y], high_contrast=true)
        scatter!(p[1], xy[x.≥y])
        plot!(p[1], bonds(H))
        plot!(p[1], bonds(l, hopping(translate_uc=[1, 1])))
        surface!(p[2], xy)
        scatter!(p[3], SquareLattice(3, 4, 5))
        true
    end
end

@testset "Lattice tests" begin
    sql = SquareLattice(10, 20)
    @test [site_index(s, sql) for s in sql] == 1:length(sql)
    x, y = coord_values(sql)
    s_0 = SquareLattice(10, 20) do site, (x, y)
        x < y
    end
    s_1 = sublattice(sql) do site, (x, y)
        x < y
    end
    s_2 = sql[x.<y]
    @test s_0 == s_1
    @test s_1 == s_2

    hl = HoneycombLattice(10, 10)
    xh, yh = coord_values(hl)
    s_3 = HoneycombLattice(10, 10) do site, (x, y)
        x < y
    end
    s_4 = sublattice(hl) do site, (x, y)
        x < y
    end
    s_5 = hl[xh.<yh]
    @test s_3 == s_4
    @test s_4 == s_5
    @test_throws MethodError hl[x.<y]
    xb, yb = coord_values(SquareLattice(5, 40))
    @test_throws "lattice mismatch" sql[xb.<yb]
end

@testset "LatticeValue tests" begin
    l = SquareLattice(10, 10)
    bas = Basis(l, 1)
    X, Y = coord_operators(bas)
    x, y = coord_values(l)
    xtr = diag_aggregate(tr, X)
    xtr2 = ptrace(X)
    xm2 = LatticeValue(l) do site, (x, y)
        2x
    end
    y = LatticeValue(l) do site, (x, y)
        y
    end
    xy = LatticeValue(l) do site, (x, y)
        x * y
    end
    @test x == xtr
    @test x == xtr2
    @test x .* y == xy
    @test x .* 2 == xm2
    @test 2 .* x == xm2
    @test x .|> (x -> 2x) == xm2
    y .= x .* y
    @test y == xy
    @test_throws "cannot broadcast" x .* ones(100)
    l2 = SquareLattice(5, 20)
    x2 = LatticeValue(l2) do site, (x, y)
        x
    end
    @test_throws "lattice mismatch" x .* x2
end

@testset "LatticeOperator tests" begin
    l = SquareLattice(10, 10)
    bas = Basis(l, 2)
    X, Y = coord_operators(bas)
    x, y = coord_values(l)
    xsq = diag_operator(bas) do site, (x, y)
        x^2
    end
    xy = diag_operator(bas) do site, (x, y)
        x * y
    end
    x2p2ym1 = diag_operator(bas) do site, (x, y)
        x * 2 + 2y - 1
    end
    xd2pxypexpy = diag_operator(bas) do site, (x, y)
        x / 2 + x * y + exp(y)
    end

    l2 = SquareLattice(5, 20)
    bas2 = Basis(l2, 2)
    X2, Y2 = coord_operators(bas2)

    @testset "In-place arithmetics" begin
        @test X^2 == xsq
        @test X * Y == xy
        @test X * 2 + 2Y - I == x2p2ym1
        @test_throws "basis mismatch" X * X2
        @test_throws MethodError X * ones(200, 200)
    end
    @testset "Lattice-value constructor" begin
        @test diag_operator(bas, x .* y) == xy
        @test (x .* y) ⊗ [1 0; 0 1] == xy
        @test [1 0; 0 1] ⊗ (x .* y) == xy
    end
    @testset "Wrapper macro" begin
        @test @on_lattice(X / 2 + X * Y + exp(Y)) == xd2pxypexpy
        @test X / 2 + X * Y + @on_lattice(exp(Y)) == xd2pxypexpy
        @test_logs (
            :warn,
            "avoid using lattice operators and matrices \
in one function call"
        ) @on_lattice X * ones(200, 200)
        @test_logs (
            :warn,
            "avoid using lattice operators and matrices \
in one function call"
        ) @on_lattice ones(200, 200) * Y
    end
end

@testset "Hopping tests" begin
    l = SquareLattice(2, 2)
    ls1, _, ls2, ls3 = l
    @testset "Constructor" begin
        @test hopping(axis=1) == hopping(translate_uc=[1])
        @test hopping(axis=1) != hopping(translate_uc=[1, 0])
        @test hopping(axis=1) == LatticeModels.promote_dims!(hopping(translate_uc=[1, 0]), 1)
        @test hopping(axis=1, pbc=[true, false]) == hopping(translate_uc=[1, 0], pbc=[true, false])
        @test hopping([-1;;], axis=1) == hopping(-1, axis=1)
        @test_throws "cannot shrink" LatticeModels.promote_dims!(hopping(translate_uc=[0, 1]), 1)
    end
    @testset "Hopping matching" begin
        @test !LatticeModels._match(hopping(axis=1), l, ls1, ls2)
        @test LatticeModels._match(hopping(axis=2), l, ls1, ls2)
        @test LatticeModels._match(hopping(axis=1), l, ls2, ls3)
        @test !LatticeModels._match(hopping(axis=2), l, ls2, ls3)
    end
    @testset "Adjacency" begin
        bs = bonds(l, hopping(axis=1), hopping(axis=2))
        bs1 = bonds(l, hopping(axis=1)) | bonds(l, hopping(axis=2))
        bs2 = bs^2
        f = pairs_by_adjacent(bs)
        f2 = pairs_by_adjacent(bs2)
        @test bs.bmat == bs1.bmat
        @test is_adjacent(bs, ls1, ls2)
        @test !is_adjacent(bs, ls1, ls3)
        @test is_adjacent(bs2, ls1, ls2)
        @test is_adjacent(bs2, ls1, ls3)
        @test f(l, site_index(ls1, l), site_index(ls2, l))
        @test !f(l, site_index(ls1, l), site_index(ls3, l))
        @test f2(l, site_index(ls1, l), site_index(ls2, l))
        @test f2(l, site_index(ls1, l), site_index(ls3, l))
    end
    x, y = coord_values(l)
    op1 = hopping_operator(pairs_by_lhs(x .< y), l, hopping(axis=1))
    op2 = hopping_operator(l, hopping(axis=1)) do l, i, j
        x, y = site_coords(l, l[i])
        x < y
    end
    @test op1 == op2
end

@testset "Field tests" begin
    @field_def struct LazyLandauField(B::Number)
        vector_potential(x) = (0, x * B)
        n_steps := 100
    end
    @field_def struct StrangeLandauField(B::Number)
        vector_potential(point...) = (0, point[1] * B)
        trip_integral(p1, p2) = 123
    end
    l = SquareLattice(10, 10)
    la = LandauField(0.1)
    lla = LazyLandauField(0.1)
    sla = StrangeLandauField(0.1)
    sym = SymmetricField(0.1)
    flx = FluxField(0.1, (0, 0))
    @testset "Trip integral" begin
        p1 = SA[1, 2]
        p2 = SA[3, 4]
        @test LatticeModels.trip_integral(la, p1, p2) ≈
              LatticeModels.trip_integral(la, p1, p2; n_steps=100)
        @test LatticeModels.trip_integral(lla, p1, p2; n_steps=100) ==
              LatticeModels.trip_integral(lla, p1, p2)
        @test LatticeModels.trip_integral(la, p1, p2) ≈
              LatticeModels.trip_integral(lla, p1, p2)

        @test LatticeModels.trip_integral(sla, p1, p2; n_steps=100) ≈
              LatticeModels.trip_integral(la, p1, p2)
        @test LatticeModels.trip_integral(sla, p1, p2) == 123

        @test LatticeModels.trip_integral(sym, p1, p2; n_steps=1) ≈
              LatticeModels.trip_integral(sym, p1, p2)
        @test LatticeModels.trip_integral(sym, p1, p2; n_steps=100) ≈
              LatticeModels.trip_integral(sym, p1, p2)

        @test LatticeModels.trip_integral(flx, p1, p2; n_steps=1000) ≈
              LatticeModels.trip_integral(flx, p1, p2) atol = 1e-8
        @test LatticeModels.trip_integral(flx + sym, p1, p2; n_steps=1000) ≈
              LatticeModels.trip_integral(flx + sym, p1, p2) atol = 1e-8
    end

    @testset "Field application" begin
        H1 = hopping_operator(l, hopping(axis=1), la) +
             hopping_operator(l, hopping(axis=2), la)
        H2 = hopping_operator(l, hopping(axis=1), lla) +
             hopping_operator(l, hopping(axis=2), lla)
        H3 = hopping_operator(l, hopping(axis=1)) +
             hopping_operator(l, hopping(axis=2))
        H4 = copy(H3)
        apply_field!(H3, la)
        apply_field!(H4, lla)
        @test H1.operator ≈ H2.operator
        @test H2.operator ≈ H3.operator
        @test H3.operator ≈ H4.operator
    end
end

@testset "Hamiltonian tests" begin
    l = HoneycombLattice(10, 10) do site, (x, y)
        x < y
    end
    x, y = coord_values(l)
    H1 = @hamiltonian begin
        lattice := l
        field := LandauField(0.5)
        @diag [1 0; 0 -1] ⊗ (@. x + y)
        @hop [1 0; 0 -1] axis = 2 pairs_by_lhs(@. x + 1 < y)
    end
    H2 = @hamiltonian begin
        lattice := l
        @hop translate_uc = [0, 1] [1 0; 0 -1] (l, i, j) -> ((x, y) = site_coords(l, l[i]); x + 1 < y)
        @diag (site, (x, y)) -> [x+y 0; 0 -x-y]
        field := LandauField(0.5)
    end
    @test H1 == H2
    sp = spectrum(H1)
    Es = eigvals(sp)
    states = eigvecs(sp)
    @test length(sp) == size(states)[2]
    @test sp[1] == sp[E=-100]
    @test filled_projector(sp).operator ≈ projector(sp[Es.<0]).operator
end

@testset "Currents tests" begin
    l = SquareLattice(10, 10)
    x, y = coord_values(l)
    H(B) = @hamiltonian begin
        lattice := l
        field := LandauField(B)
        @diag [1 0; 0 -1]
        @hop axis = 1 [1 im; im -1] / 2
        @hop axis = 2 [1 1; -1 -1] / 2
    end
    P = filled_projector(spectrum(H(0)))
    dc = DensityCurrents(H(0.1), P)
    bs = bonds(H(0.1))
    m1 = materialize(dc)[x.<y]
    m2 = materialize(dc[x.<y])
    m3 = materialize(pairs_by_adjacent(bs), dc)[x.<y]
    m4 = materialize(pairs_by_adjacent(bs), dc[x.<y])
    md = materialize(pairs_by_distance(≤(2)), dc[x.<y])
    @test m1.currents == m2.currents
    @test m1.currents == m3.currents
    @test m1.currents == m4.currents
    @test m1.currents == md.currents
end