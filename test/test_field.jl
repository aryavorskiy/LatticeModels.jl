@testset "Field" begin
    struct LazyLandauGauge <: LatticeModels.AbstractField
        B::Number
    end
    LatticeModels.vector_potential(field::LazyLandauGauge, p1) = (0, p1[1] * field.B)
    struct StrangeLandauGauge <: LatticeModels.AbstractField end
    LatticeModels.vector_potential(::StrangeLandauGauge, point) = (0, point[1] * 0.1)
    LatticeModels.line_integral(::StrangeLandauGauge, p1, p2) = 123
    struct EmptyField <: LatticeModels.AbstractField end
    l = SquareLattice(10, 10)
    la = LandauGauge(0.1)
    lla = LazyLandauGauge(0.1)
    sla = StrangeLandauGauge()
    sym = SymmetricGauge(0.1)
    flx = PointFlux(0.1)
    @test flx.point == (0, 0)
    emf = EmptyField()
    fla = GaugeField(n = 10) do (x,)
        (0, x * 0.1)
    end
    @testset "Line integral" begin
        p1 = SA[1, 2]
        p2 = SA[3, 4]
        @test LatticeModels.line_integral(la, p1, p2) ≈
              LatticeModels.line_integral(la, p1, p2, 100)
        @test LatticeModels.line_integral(lla, p1, p2) ≈
              LatticeModels.line_integral(lla, p1, p2, 100)
        @test LatticeModels.line_integral(la, p1, p2) ≈
              LatticeModels.line_integral(lla, p1, p2)
        @test LatticeModels.line_integral(la, p1, p2) ≈
              LatticeModels.line_integral(fla, p1, p2)

        @test LatticeModels.line_integral(sla, p1, p2, 100) ≈
              LatticeModels.line_integral(la, p1, p2)
        @test LatticeModels.line_integral(sla, p1, p2) == 123

        @test LatticeModels.line_integral(sym, p1, p2, 1) ≈
              LatticeModels.line_integral(sym, p1, p2)
        @test LatticeModels.line_integral(sym, p1, p2, 100) ≈
              LatticeModels.line_integral(sym, p1, p2)

        @test LatticeModels.line_integral(flx, p1, p2, 1000) ≈
              LatticeModels.line_integral(flx, p1, p2) atol = 1e-8
        @test LatticeModels.line_integral(flx + sym, p1, p2, 1000) ≈
              LatticeModels.line_integral(flx + sym, p1, p2) atol = 1e-8

        @test_throws ArgumentError LatticeModels.vector_potential(emf, SA[1, 2, 3])
    end

    @testset "Field application" begin
        H1 = construct_operator(l, BravaisTranslation(axis=1), BravaisTranslation(axis=2), field = la)
        H2 = construct_operator(l, BravaisTranslation(axis=1), BravaisTranslation(axis=2), field = lla)
        H3 = construct_operator(l, BravaisTranslation(axis=1), BravaisTranslation(axis=2))
        H4 = copy(H3)
        apply_field!(H3, la)
        apply_field!(H4, lla)
        @test H1 ≈ H2
        @test H2 ≈ H3
        @test H3 ≈ H4
    end
end
