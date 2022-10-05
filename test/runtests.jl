using Test, LinearAlgebra
using LatticeModels

@testset "Workflows" begin
    @test begin
        l = SquareLattice(10, 10) do x, y; √(x^2 + y^2) < 5; end
        b = Basis(l, 1)
        d = diag_operator(b) do site; 1 end
        hx = hopping_operator(l, Hopping(axis=1))
        hy = hopping_operator(l, Hopping(axis=2))
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
            field := Landau(3)
            @diag [1 0; 0 -1]
            @hop axis=1 [1 1; 1 -1] / √2
            @hop axis=2 [1 -im; im -1] / √2
        end
        P = filled_projector(spectrum(H))
        X, Y = coord_operators(Basis(l, 2))
        d = diag_aggregate(tr, 4π * im * P * X * (I - P) * Y * P)
        rd = d .|> real
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        H = @hamiltonian begin
            lattice := l
            field := Landau(0.5)
            @diag [1 0; 0 -1]
            @hop axis=1 [1 im; im -1] / 2
            @hop axis=2 [1 1; -1 -1] / 2
        end
        h(t) = @hamiltonian begin
            lattice := l
            field := Landau(t)
            @diag (x, y) -> if √(x^2 + y^2) ≤ 2
                [1 0; 0 -1]
            else
                [3 0; 0 -3]
            end
            @hop axis=1 [1 im; im -1] / 2
            @hop axis=2 [1 1; -1 -1] / 2
        end
        P0 = filled_projector(spectrum(H))
        X, Y = coord_operators(Basis(l, 2))
        @evolution {P0 => h(t) => P} for t in 0:0.1:10
            d = diag_aggregate(tr, 4π * im * P * X * (I - P) * Y * P)
            rd = d .|> real
        end
        true
    end
end

@testset "LatticeValue tests" begin
    l = SquareLattice(10, 10)
    bas = Basis(l, 1)
    X, Y = coord_operators(bas)
    xtr = diag_aggregate(tr, X)
    x = LatticeValue(l) do x, y; x; end
    xm2 = LatticeValue(l) do x, y; 2x; end
    y = LatticeValue(l) do x, y; y; end
    xy = LatticeValue(l) do x, y; x * y; end
    @test x == xtr
    @test x .* y == xy
    @test x .* 2 == xm2
    @test 2 .* x == xm2
    @test x .|> (x -> 2x) == xm2
    y .= x .* y
    @test y == xy
    @test_throws "cannot broadcast" x .* ones(100)
    l2 = SquareLattice(5, 20)
    x2 = LatticeValue(l2) do x, y; x; end
    @test_throws "lattice mismatch" x .* x2
end

@testset "LatticeOperator tests" begin
    l = SquareLattice(10, 10)
    bas = Basis(l, 2)
    X, Y = coord_operators(bas)
    xsq = diag_operator(bas) do x, y; x^2; end
    xy = diag_operator(bas) do x, y; x * y; end
    x2p2ym1 = diag_operator(bas) do x, y; x * 2 + 2y - 1; end
    xd2pxypexpy = diag_operator(bas) do x, y; x / 2 + x * y + exp(y); end

    l2 = SquareLattice(5, 20)
    bas2 = Basis(l2, 2)
    X2, Y2 = coord_operators(bas2)

    @testset "In-place arithmetics" begin
        @test X ^ 2 == xsq
        @test X * Y == xy
        @test X * 2 + 2Y - I == x2p2ym1
        @test_throws "basis mismatch" X * X2
        @test_throws MethodError X * ones(200, 200)
    end
    @testset "Wrapper macro" begin
        @test @on_lattice(X / 2 + X * Y + exp(Y)) == xd2pxypexpy
        @test X / 2 + X * Y + @on_lattice(exp(Y)) == xd2pxypexpy
        @test_logs (:warn, "avoid using lattice operators and matrices \
            in one function call") @on_lattice X * ones(200, 200)
    end
end

@testset "Hoppings" begin
    l = SquareLattice(2, 2)
    ls1 = LatticeIndex([1,1],1)
    ls2 = LatticeIndex([1,2],1)
    ls3 = LatticeIndex([2,2],1)
    @testset "Hopping matching" begin
        @test !LatticeModels._match(Hopping(axis=1), l, ls1, ls2)
        @test LatticeModels._match(Hopping(axis=2), l, ls1, ls2)
    end
    @testset "Adjacency" begin
        bs = bonds(l, Hopping(axis=1), Hopping(axis=2))
        bs2 = bs ^ 2
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
