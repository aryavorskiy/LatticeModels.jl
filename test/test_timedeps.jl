@testset "Timedeps" begin
    @testset "TimeSequence" begin
        l = SquareLattice(3, 3)
        site = l[5]
        x, y = coordvalues(l)
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
        rec2 = TimeSequence{LatticeValue}()
        rec2[0] = xy * 0
        rec2[1] = xy * 1
        rec2[2] = xy * 2
        @test rec2[1e-9] == zeros(l)
        @test rec2[1 + 1e-9] == xy
        @test rec2[t = 0.9 .. 2.1] == TimeSequence([1, 2], [xy, xy * 2])
        @test integrate(rec) == rec2
        @test_throws KeyError rec2[0.5]

        ts = TimeSequence(0:0.1:10, (0:0.1:10).^ 2)
        ts2 = empty(ts)
        @test typeof(ts2) == TimeSequence{Float64}
        @test length(ts2) == 0
        for t in 7:0.1:10
            delete!(ts, t)
        end
        delete!(ts, t=5 .. 7)
        for t in 0:0.1:4.9
            ts2[t] = t^2
        end
        @test ts == ts2
    end

    @testset "Evolution" begin
        l = SquareLattice(10, 10)
        Hs = qwz(l)     # sparse
        Hd = dense(Hs)  # dense
        psi = groundstate(Hd)

        ts = 0:0.1:10
        val1 = ComplexF64[]
        for state in Evolution(KrylovKitExp, Hs, psi)(ts)
            ψ = state[1]
            push!(val1, ψ.data[2])
        end
        val2 = ComplexF64[]
        for st in Evolution(CachedExp(threshold=1e-12), Hs, psi, timedomain=ts)
            ψ = st.state
            push!(val2, ψ.data[2])
        end
        correct_val = ComplexF64[]
        correct_ev = exp(-im * step(ts) * Hd.data)
        psidata = psi.data
        for _ in ts
            push!(correct_val, psidata[2])
            psidata = correct_ev * psidata
        end
        @test val1 ≈ correct_val atol=1e-10
        @test val2 ≈ correct_val atol=1e-10
    end
end
