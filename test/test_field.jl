import LatticeModels: line_integral

@testset "Field" begin
    struct LazyLandauGauge <: LatticeModels.AbstractField
        B::Number
    end
    LatticeModels.vector_potential(field::LazyLandauGauge, p) = (0, p[1] * field.B)
    struct StrangeLandauGauge <: LatticeModels.AbstractField end
    LatticeModels.vector_potential(::StrangeLandauGauge, p) = (0, p[1] * 0.1)
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
        @test line_integral(la, p1, p2) ≈ line_integral(la, p1, p2, 100)
        @test line_integral(lla, p1, p2) ≈ line_integral(lla, p1, p2, 100)
        @test line_integral(la, p1, p2) ≈ line_integral(lla, p1, p2)
        @test line_integral(la, p1, p2) ≈ line_integral(fla, p1, p2)

        @test line_integral(sla, p1, p2, 100) ≈ line_integral(la, p1, p2)
        @test line_integral(sla, p1, p2) == 123

        @test line_integral(sym, p1, p2, 1) ≈ line_integral(sym, p1, p2)
        @test line_integral(sym, p1, p2, 100) ≈ line_integral(sym, p1, p2)

        @test line_integral(flx, p1, p2, 1000) ≈ line_integral(flx, p1, p2) atol = 1e-8
        @test line_integral(flx + sym, p1, p2, 1000) ≈ line_integral(flx + sym, p1, p2) atol = 1e-8

        @test_throws ArgumentError LatticeModels.vector_potential(emf, SA[1, 2, 3])
    end

    @testset "Point fluxes" begin
        p1 = SA[0, 0]
        p2 = SA[4, 0]

        pf1 = PointFlux(0.1, (1, 2))
        pf2 = PointFlux(0.1, (3, 4))
        pfs2 = PointFlux(0.1, (3, 4), gauge=:singular)
        ps1 = PointFluxes(0.1, [(1, 2), (3, 4)], gauge=:singular)
        ps2 = PointFluxes([pf1, pfs2], gauge=:singular)
        ps3 = PointFluxes([pf1, pf2])

        @test_throws TypeError PointFluxes([pf1, 3])
        @test_throws ArgumentError PointFluxes([pf1, pfs2])

        @test line_integral(pf1 + pf2, p1, p2) ≈ line_integral(ps3, p1, p2) atol = 1e-8
        @test line_integral(ps1, p1, p2) ≈ 0.2 atol = 1e-8
        @test line_integral(ps2, p1, p2) ≈ 0.2 atol = 1e-8

        bf1 = PointFluxes()
        push!(bf1, PointFlux(0.1, (1, 2)), PointFlux(0.1, (3, 4)))
        bf2 = PointFluxes()
        append!(bf2, ps3)
        @test line_integral(bf1, p1, p2) ≈ line_integral(ps3, p1, p2) atol = 1e-8
        @test line_integral(bf2, p1, p2) ≈ line_integral(ps3, p1, p2) atol = 1e-8

        ppfs1 = periodic_fluxes(l, PointFlux(0.1, (0.4, 0.4)))
        ppfs2 = periodic_fluxes(l, PointFlux(0.1, (1.4, 0.4)))
        ppfs3 = PointFluxes(0.1, l, offset=(0.4, 0.4))
        # Sort the points, because the order of the points in `periodic_fluxes` is not guaranteed
        @test all(zip(sort(ppfs1.points), sort(ppfs2.points))) do (p1, p2)
            sqrt(sum(abs2, p1 .- p2)) < 1e-10
        end
        @test all(zip(sort(ppfs1.points), sort(ppfs3.points))) do (p1, p2)
            sqrt(sum(abs2, p1 .- p2)) < 1e-10
        end
    end

    @testset "Fluxes and boundaries" begin
        lnb = SquareLattice(5, 5)
        lwb = setboundaries(lnb, :axis1 => true, :axis2 => true)

        flx = PointFlux(0.1, (0.5, 0.5), gauge=:singular)
        flxs = periodic_fluxes(lnb, flx)
        flxs2 = periodic_fluxes(lwb, flx)
        @test flxs.points == flxs2.points

        flx_adn = LatticeModels.adapt_field(flx, lnb)
        flxs_adn = LatticeModels.adapt_field(flxs, lnb)
        flx_adw = LatticeModels.adapt_field(flx, lwb)
        flxs_adw = LatticeModels.adapt_field(flxs, lwb)
        @test flx_adn isa PointFlux
        @test length(flxs_adn.points) == 25
        @test length(flx_adw.points) == 9
        @test length(flxs_adw.points) == 9 * 25
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
