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
        heatmap(xy)
        H = @hamiltonian begin
            lattice := l
            field := LandauField(0.5)
            @diag [1 0; 0 -1]
            @hop axis = 1 [1 im; im -1] / 2
            @hop axis = 2 [1 1; -1 -1] / 2
        end
        P = filled_projector(spectrum(H), 0.1)
        quiver!(DensityCurrents(H, P)[x.<y])
        sl = l[x.<y]
        scatter!(sl)
        plot!(bonds(H))
        plot!(bonds(l, hopping(tr_vector=[1, 1])))
        true
    end
end

@testset "Lattice tests" begin
    sql = SquareLattice(10, 20)
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
    ls1 = LatticeSite(SA[1, 1], 1)
    ls2 = LatticeSite(SA[1, 2], 1)
    ls3 = LatticeSite(SA[2, 2], 1)
    @testset "Constructor" begin
        @test hopping(axis=1) == hopping(tr_vector=[1])
        @test hopping(axis=1) != hopping(tr_vector=[1, 0])
        @test hopping(axis=1) == LatticeModels.promote_dims!(hopping(tr_vector=[1, 0]), 1)
        @test hopping(axis=1, pbc=[true, false]) == hopping(tr_vector=[1, 0], pbc=[true, false])
        @test hopping([-1;;], axis=1) == hopping(-1, axis=1)
        @test_throws "cannot shrink" LatticeModels.promote_dims!(hopping(tr_vector=[0, 1]), 1)
    end
    @testset "Hopping matching" begin
        @test !LatticeModels._match(hopping(axis=1), l, ls1, ls2)
        @test LatticeModels._match(hopping(axis=2), l, ls1, ls2)
        @test LatticeModels._match(hopping(axis=1), l, ls2, ls3)
        @test !LatticeModels._match(hopping(axis=2), l, ls2, ls3)
    end
    @testset "Adjacency" begin
        bs = bonds(l, hopping(axis=1), hopping(axis=2))
        bs2 = bs^2
        f = is_adjacent(bs)
        f2 = is_adjacent(bs2)
        @test is_adjacent(bs, ls1, ls2)
        @test !is_adjacent(bs, ls1, ls3)
        @test is_adjacent(bs2, ls1, ls2)
        @test is_adjacent(bs2, ls1, ls3)
        @test f(ls1, ls2)
        @test !f(ls1, ls3)
        @test f2(ls1, ls2)
        @test f2(ls1, ls3)
    end
end

@testset "Field tests" begin
    @field_def struct LazyLandauField(B::Number)
        vector_potential(x) = (0, x * B)
        n_integrate := 100
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
    @testset "Trip integral" begin
        p1 = SA[1, 2]
        p2 = SA[3, 4]
        @test LatticeModels.trip_integral(la, p1, p2) ≈
              LatticeModels.trip_integral(lla, p1, p2)
        @test LatticeModels.trip_integral(la, p1, p2; n_integrate=100) ≈
              LatticeModels.trip_integral(lla, p1, p2)
        @test LatticeModels.trip_integral(lla, p1, p2; n_integrate=100) ==
              LatticeModels.trip_integral(lla, p1, p2)

        @test LatticeModels.trip_integral(sla, p1, p2; n_integrate=100) ≈
              LatticeModels.trip_integral(la, p1, p2)
        @test LatticeModels.trip_integral(sla, p1, p2) == 123

        @test LatticeModels.trip_integral(sym, p1, p2; n_integrate=1) ≈
              LatticeModels.trip_integral(sym, p1, p2)
        @test LatticeModels.trip_integral(sym, p1, p2; n_integrate=100) ≈
              LatticeModels.trip_integral(sym, p1, p2)
    end

    @testset "Field application" begin
        H1 = hopping_operator(l, hopping(axis=1), field=la) +
             hopping_operator(l, hopping(axis=2), field=la)
        H2 = hopping_operator(l, hopping(axis=1), field=lla) +
             hopping_operator(l, hopping(axis=2), field=lla)
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
        @hop [1 0; 0 -1] axis = 2 (@. x + 1 < y)
    end
    H2 = @hamiltonian begin
        lattice := l
        @hop tr_vector = [0, 1] [1 0; 0 -1] (site, (x, y)) -> (x + 1 < y)
        @diag (site, (x, y)) -> [x+y 0; 0 -x-y]
        field := LandauField(0.5)
    end
    @test H1 == H2
    sp = spectrum(H1)
    Es = eigvals(sp)
    states = eigvecs(sp)
    @test length(sp) == size(states)[2]
    @test sp[1] == sp[E=-100]
    @test filled_projector(sp).operator ≈ projector(sp[Es .< 0]).operator
end

@testset "Currents tests" begin
    l = SquareLattice(10, 10)
    x, y = coord_values(l)
    H = @hamiltonian begin
        lattice := l
        field := LandauField(0.5)
        @diag [1 0; 0 -1]
        @hop axis = 1 [1 im; im -1] / 2
        @hop axis = 2 [1 1; -1 -1] / 2
    end
    P = filled_projector(spectrum(H))
    dc = DensityCurrents(H, P)
    bs = bonds(H)
    m1 = materialize(dc)[x.<y]
    m2 = materialize(dc[x.<y])
    m3 = materialize(is_adjacent(bs), dc[x.<y])
    @test m1.currents == m2.currents
    @test m2.currents == m3.currents
end
