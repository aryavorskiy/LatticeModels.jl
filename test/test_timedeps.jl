@testset "Timedeps" begin
    @testset "TimeSequence" begin
        l = SquareLattice(3, 3)
        site = l[5]
        x, y = coord_values(l)
        xy = x .* y
        xly = x .< y
        @test_throws ArgumentError TimeSequence([0.5], [xy, xy])
        rec = TimeSequence{LatticeValue}()
        rec[0] = xy
        rec[1] = xy
        rec[2] = xy
        @test rec[site, t = 0..Inf] == TimeSequence(timestamps(rec), fill(xy[site], 3))
        @test rec[xly, t = 0..Inf] == TimeSequence(0:2, fill(xy[xly], 3))
        @test collect(rec) == [0 => xy, 1 => xy, 2 => xy]
        @test differentiate(rec) == TimeSequence([0.5, 1.5], [zeros(l), zeros(l)])
        @test differentiate(rec[site, t = 0..Inf]) == TimeSequence([0.5, 1.5], [0, 0])
        rec2 = TimeSequence(0, xy .* 0)
        rec2[1] = xy .* 1
        rec2[2] = xy .* 2
        @test rec2[1e-9] == zeros(l)
        @test rec2[1 + 1e-9] == xy
        @test rec2[t = 0.9..2.1] == TimeSequence([1, 2], [xy, xy .* 2])
        @test integrate(rec) == rec2
        @test_throws KeyError rec2[0.5]
    end

    @testset "Evolution" begin
        l = SquareLattice(10, 10)
        Hs = qwz(l)     # sparse
        Hd = dense(Hs)  # dense

        τ = 10
        correct_ev = exp(-im * τ * Hd.data)
        Evs = LatticeModels.evolution_operator!(Operator(Hs), Hs, τ)
        Evd = Operator(basis(Hd), one(Hd.data))
        LatticeModels.evolution_operator!(Evd, Hd, τ)
        @test correct_ev ≈ Evs.data atol=1e-6
        @test correct_ev ≈ Evd.data atol=1e-6
    end
end
