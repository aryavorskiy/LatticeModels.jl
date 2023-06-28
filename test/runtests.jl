using Test, LinearAlgebra, StaticArrays, Plots
using LatticeModels

@static if VERSION < v"1.8"
    # Ignore test_throws
    macro test_throws(args...)
        :()
    end
end

@testset "Workflows" begin
    @test begin
        l = SquareLattice(10, 10) do (x, y)
            √(x^2 + y^2) < 5
        end
        b = Basis(l, 1)
        d = LatticeOperator(b, I)
        hx = hopping_operator(l, hopping(axis=1))
        hy = hopping_operator(l, hopping(axis=2))
        H = d + hx + hy
        sp = spectrum(H)
        P = filled_projector(sp)
        d = diag_reduce(tr, P)
        true
    end

    @test begin
        l = HoneycombLattice(10, 10)
        H = Haldane(l, 1, 1, 1)
        P = filled_projector(spectrum(H))
        X, Y = coord_operators(basis(H))
        d = site_density(4π * im * P * X * (I - P) * Y * P)
        rd = d .|> real
        true
    end

    @test begin
        l = SquareLattice(10, 10)
        H0 = @hamiltonian begin
            lattice := l
            field := LandauField(0.5)
            dims_internal := 2
            @diag [1 0; 0 -1]
            @diag rand(l)
            @hop axis = 1 [1 im; im -1] / 2
            @hop axis = 2 [1 1; -1 -1] / 2
        end
        function h(t)
            x, y = coord_values(l)
            ms = @. 3 + (√(x^2 + y^2) ≤ 2) * -2
            @hamiltonian begin
                lattice := l
                field := LandauField(t)
                dims_internal := 2
                @diag ms ⊗ [1 0; 0 -1]
                @diag randn(l)
                @hop axis = 1 [1 im; im -1] / 2
                @hop axis = 2 [1 1; -1 -1] / 2
            end
        end
        P0 = filled_projector(spectrum(H0))
        X, Y = coord_operators(Basis(l, 2))
        @evolution {H := h(t), P0 --> H --> P} for t in 0:0.1:10
            d = diag_reduce(tr, 4π * im * P * X * (I - P) * Y * P)
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
        p = plot(layout=4)
        plot!(p[1], xy)
        H = @hamiltonian begin
            lattice := l
            field := LandauField(0.5)
            dims_internal := 2
            @diag [1 0; 0 -1]
            @hop axis = 1 [1 im; im -1] / 2
            @hop axis = 2 [1 1; -1 -1] / 2
        end
        P = filled_projector(spectrum(H), 0.1)
        @evolution k = 2 pade = true {P --> H --> PP} for t in 0:0.1:1
        end
        @evolution k = 2 {P --> H --> PP} for t in 0:0.1:1
        end

        dc = DensityCurrents(H, P)
        quiver!(p[1], dc[x.<y])
        scatter!(p[1], l)
        scatter!(p[1], l[x.<y], high_contrast=true)
        scatter!(p[1], xy[x.≥y])
        plot!(p[1], bonds(H))
        plot!(p[1], bonds(l, hopping(translate_uc=[1, 1])))
        surface!(p[2], xy)
        scatter!(p[3], SquareLattice(3, 4, 5))
        plot!(p[4], project(xy, :x))
        plot!(p[4], project(xy, :j1))
        mpcs = map_currents(site_distance, dc, reduce_fn=sum, sort=true)
        plot!(p[4], mpcs)
        true
    end
end

@testset "Lattice" begin
    sql = SquareLattice(10, 20)
    @test [site_index(sql, s) for s in sql] == 1:length(sql)
    x, y = coord_values(sql)
    s_0 = SquareLattice(10, 20) do (x, y)
        x < y
    end
    s_1 = sublattice(sql) do (x, y)
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
    s_4 = sublattice(hl) do (x, y)
        x < y
    end
    s_5 = hl[xh.<yh]
    @test s_3 == s_4
    @test s_4 == s_5
    @test_throws MethodError hl[x.<y]
    @test hl[j1=3, j2=2, index=1] == LatticeSite(SA[3, 2], 1, SA[4.0, 1.7320508075688772])
    xb, yb = coord_values(SquareLattice(5, 40))
    @test_throws "macrocell mismatch" sql[xb.<yb]
    sql2 = sublattice(sql) do site
        site ∉ (sql[1], sql[end])
    end
    pop!(sql)
    popfirst!(sql)
    @test sql == sql2

end

@testset "LatticeValue" begin
    l = SquareLattice(2, 2)
    x, y = coord_values(l)
    @test x.values == [1, 2, 1, 2]
    @test y.values == [1, 1, 2, 2]
    l = SquareLattice(10, 10)
    bas = Basis(l, 1)
    X, Y = coord_operators(bas)
    x, y = coord_values(l)
    xtr = diag_reduce(tr, X)
    xtr2 = site_density(X)
    xm2 = LatticeValue(l) do (x, y)
        2x
    end
    y = LatticeValue(l) do (x, y)
        y
    end
    xy = LatticeValue(l) do site
        site.x * site.y
    end
    idxs = LatticeValue(l) do site
        site_index(l, site)
    end
    @test [idxs[s] for s in l] == 1:length(l)
    @test x == xtr
    @test x == xtr2
    @testset "Broadcast" begin
        @test x .* y == xy
        @test x .* 2 == xm2
        @test 2 .* x == xm2
        @test x .|> (x -> 2x) == xm2
        y .= x .* y
        @test y == xy
        @test_throws "cannot broadcast" x .* ones(100)
        l2 = SquareLattice(5, 20)
        x2 = LatticeValue(l2) do (x, y)
            x
        end
        @test_throws "lattice mismatch" x .* x2
    end
    @testset "Indexing" begin
        z = zeros(l)
        z2 = zero(z)
        z .= 1
        z2[x.≥y] .+= 1
        z2[x.<y] = ones(l)
        @test z == ones(l)
        @test z2 == ones(l)
        @test z[x=1, x2=1] == 1
        z3 = ones(l)
        for i in 2:10
            z3[x1=i] .= i
        end
        @test z3 == x
    end
    @testset "Interface" begin
        (mn, i) = findmin(xy)
        (mx, j) = findmax(xy)
        @test xy[i] == mn
        @test xy[j] == mx
        @test all(mn .≤ xy.values .≤ mx)
        imx = findall(==(mx), xy)
        @test all((site in imx || xy[site] < mx) for site in l)
    end
end

@testset "LatticeArray" begin
    l = SquareLattice(10, 10)
    bas = Basis(l, 2)
    X, Y = coord_operators(bas)
    x, y = coord_values(l)
    xsq = diag_operator(bas) do site
        site.x^2
    end
    xy = diag_operator(bas) do site
        site.x * site.y * [1 0; 0 1]
    end
    x2p2ym1 = diag_operator(bas) do (x, y)
        (x * 2 + 2y - 1) * [1 0; 0 1]
    end
    xd2pxypexpy = diag_operator(bas) do (x, y)
        (x / 2 + x * y + exp(y)) * [1 0; 0 1]
    end

    l2 = SquareLattice(5, 20)
    bas2 = Basis(l2, 2)
    X2, Y2 = coord_operators(bas2)
    vcv = zeros(length(bas))
    vcv[1] = 1
    vc = LatticeModels.LatticeArray(bas, vcv)
    xf = l[1].x

    @testset "In-place arithmetics" begin
        @test X - Y == X + -Y
        @test X^2 == xsq
        @test X * Y == xy
        @test X * 2 + 2Y - I == x2p2ym1
        @test vc' * X * vc == xf
        @test dot(vc, vc) == 1
        @test dot(vc, X, vc) == xf
        @test ptrace(X, :internal) == diagm(x.values) * 2
        @test ptrace(X, :lattice) == sum(x.values) * [1 0; 0 1]
        @test_throws "basis mismatch" X * X2
        @test_throws MethodError X * ones(200, 200)
    end
    @testset "Lattice-value constructor" begin
        @test diag_operator(x .* y, 2) == xy
        @test (x .* y) ⊗ [1 0; 0 1] == xy
        @test [1 0; 0 1] ⊗ (x .* y) == xy
        @test_throws "Lambda returned a String" diag_operator(bas2) do _
            "kek"
        end
    end
    @testset "Wrapper macro" begin
        @test @on_lattice(X / 2 + X * Y + exp(Y)) == xd2pxypexpy
        @test X * 0.5 + X * Y + @on_lattice(exp(Y)) == xd2pxypexpy
        @test_logs (
            :warn,
            "avoid using lattice arrays and unwrapped arrays in one function call"
        ) @on_lattice X * ones(200, 200)
        @test_logs (
            :warn,
            "avoid using lattice arrays and unwrapped arrays in one function call"
        ) @on_lattice ones(200, 200) * Y
    end
end

@testset "TimeSequence" begin
    l = SquareLattice(3, 3)
    site = l[5]
    x, y = coord_values(l)
    xy = x .* y
    xly = x .< y
    @test_throws "length mismatch" TimeSequence([0.5], [xy, xy])
    rec = LatticeValueSequence()
    insert!(rec, 0, xy)
    insert!(rec, 1, xy)
    insert!(rec, 2, xy)
    @test rec[site] == TimeSequence(time_domain(rec), fill(xy[site], 3))
    @test rec[xly] == TimeSequence(0:2, fill(xy[xly], 3))
    @test collect(rec) == [0 => xy, 1 => xy, 2 => xy]
    @test differentiate(rec) == TimeSequence([0.5, 1.5], [zeros(l), zeros(l)])
    @test differentiate(rec[site]) == TimeSequence([0.5, 1.5], [0, 0])
    rec2 = TimeSequence(xy .* 0)
    insert!(rec2, 1, xy .* 1)
    insert!(rec2, 2, xy .* 2)
    @test rec2(0.9) == xy
    @test rec2(0.9, 2.1) == TimeSequence([1, 2], [xy, xy .* 2])
    @test integrate(rec) == rec2
end

@testset "Hopping" begin
    @testset "Constructor" begin
        @test hopping(axis=1) == hopping(translate_uc=[1])
        @test hopping(axis=1) != hopping(translate_uc=[1, 0])
        @test hopping(axis=1) == LatticeModels.promote_dims!(hopping(translate_uc=[1, 0], pbc=[false, true]), 1)
        @test hopping(site_indices=(1, 2)) == LatticeModels.promote_dims!(hopping(site_indices=(1, 2), pbc=[false, true]), 1)
        @test hopping(axis=2, pbc=true) == hopping(translate_uc=[0, 1], pbc=[true, true])
        @test hopping(axis=2, pbc=[true, true]) == hopping(translate_uc=[0, 1], pbc=[true, true])
        @test hopping(fill(-1, 1, 1), axis=1) == hopping(-1, axis=1)
        @test_throws "to hopping indices" hopping(site_indices=(1, 2, 3))
        @test_throws "connects site to itself" hopping(translate_uc=[0, 0], site_indices=2)
        @test_throws "cannot shrink" LatticeModels.promote_dims!(hopping(translate_uc=[0, 1]), 1)
    end
    @testset "Hopping matching" begin
        l = SquareLattice(6, 5)
        hx = hopping(axis=1)
        LatticeModels.promote_dims!(hx, dims(l))
        hxmy = hopping(translate_uc=[1, -1], pbc=[true, false])
        for site in l
            ucx, ucy = site.unit_cell
            dst_dx = LatticeModels.hopping_dest(l, hx, site)
            dst_dxmy = LatticeModels.hopping_dest(l, hxmy, site)
            @test (dst_dx === nothing) == (ucx == 6)
            @test (dst_dxmy === nothing) == (ucy == 1)
        end
    end
    @testset "Bonds" begin
        l = SquareLattice(2, 2)
        ls1, _, ls3, ls4 = l
        bs = bonds(l, hopping(axis=1), hopping(axis=2))
        bs1 = bonds(l, hopping(axis=1)) | bonds(l, hopping(axis=2))
        bs2 = bs^2
        @test bs.bmat == (!!bs).bmat
        @test bs.bmat == bs1.bmat
        @test bs(l, ls1, ls3)
        @test !bs(l, ls1, ls4)
        @test bs2(l, ls1, ls3)
        @test bs2(l, ls1, ls4)
    end
    l = SquareLattice(5, 5) do site
        local (x, y) = site.coords
        abs(x + y) < 2.5
    end
    x, y = coord_values(l)
    op1 = hopping_operator(PairLhsSelector(x .< y), l, hopping(axis=1))
    op2 = hopping_operator(l, hopping(axis=1)) do l, site1, site2
        local (x, y) = site1.coords
        x < y
    end
    @test op1 == op2
end

@testset "Field" begin
    @field_def struct LazyLandauField(B::Number)
        vector_potential(x) = (0, x * B)
        n_steps := 100
        show(io::IO, ::MIME"text/plain") = print(io, "Lazy Landau calibration field; B = $B flux quanta per 1×1 plaquette")
    end
    @field_def struct StrangeLandauField
        vector_potential(point...) = (0, point[1] * 0.1)
        path_integral(p1, p2) = 123
    end
    @field_def struct EmptyField end
    l = SquareLattice(10, 10)
    la = LandauField(0.1)
    lla = LazyLandauField(0.1)
    sla = StrangeLandauField()
    sym = SymmetricField(0.1)
    flx = FluxField(0.1)
    flx2 = FluxField(0.1, (0, 0))
    @test flx.P == flx2.P
    emf = EmptyField()
    @testset "Path integral" begin
        p1 = SA[1, 2]
        p2 = SA[3, 4]
        @test LatticeModels.path_integral(la, p1, p2) ≈
              LatticeModels.path_integral(la, p1, p2, 100)
        @test LatticeModels.path_integral(lla, p1, p2, 100) ==
              LatticeModels.path_integral(lla, p1, p2)
        @test LatticeModels.path_integral(la, p1, p2) ≈
              LatticeModels.path_integral(lla, p1, p2)

        @test LatticeModels.path_integral(sla, p1, p2, 100) ≈
              LatticeModels.path_integral(la, p1, p2)
        @test LatticeModels.path_integral(sla, p1, p2) == 123

        @test LatticeModels.path_integral(sym, p1, p2, 1) ≈
              LatticeModels.path_integral(sym, p1, p2)
        @test LatticeModels.path_integral(sym, p1, p2, 100) ≈
              LatticeModels.path_integral(sym, p1, p2)

        @test LatticeModels.path_integral(flx, p1, p2, 1000) ≈
              LatticeModels.path_integral(flx, p1, p2) atol = 1e-8
        @test LatticeModels.path_integral(flx + sym, p1, p2, 1000) ≈
              LatticeModels.path_integral(flx + sym, p1, p2) atol = 1e-8

        @test_throws "no vector potential function" LatticeModels.vector_potential(emf, SA[1, 2, 3])
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
        @test H1.array ≈ H2.array
        @test H2.array ≈ H3.array
        @test H3.array ≈ H4.array
    end
end

@testset "Hamiltonian" begin
    @testset "Hamiltonian macro" begin
        l = HoneycombLattice(10, 10) do site
            local (x, y) = site.coords
            x < y
        end
        x, y = coord_values(l)
        H0 = [1 0; 0 -1] ⊗ (@. x + y) + hopping_operator(l, hopping([1 0; 0 -1], axis=2), LandauField(0.5)) do l, site1, site2
            local (x, y) = site1.coords
            x + 1 < y
        end
        H1 = @hamiltonian begin
            lattice := l
            field := LandauField(0.5)
            dims_internal := 2
            @diag [1 0; 0 -1] ⊗ (@. x + y)
            @hop [1 0; 0 -1] axis = 2 PairLhsSelector(@. x + 1 < y)
        end
        H2 = @hamiltonian begin
            lattice := l
            dims_internal := 2
            @hop translate_uc = [0, 1] pbc = [true, false] [1 0; 0 -1] (l, s1, s2) ->
                ((x, y) = s1.coords; x + 1 < y)
            @diag site -> [site.x+site.y 0; 0 -site.x-site.y]
            field := LandauField(0.5)
        end
        @test H1 == H0
        @test H2 == H0

        H3 = hopping_operator(l, hopping(-0.25, axis=1)) + hopping_operator(l, hopping(0.25, axis=2))
        H4 = @hamiltonian begin
            lattice := l
            @hop -0.25 axis = 1
            @hop fill(0.25, 1, 1) axis = 2
        end
        @test H3 == H4
        sp = spectrum(H1)
        Es = eigvals(sp)
        states = eigvecs(sp)
        @test length(sp) == size(states)[2]
        @test sp[1] == sp[E=-100]
        @test site_density(sp[1]).values ≈ site_density(projector(sp[1:1])).values
        @test filled_projector(sp).array ≈ projector(sp[Es.<0]).array
    end

    @testset "Hamiltonian builtins" begin
        l = SquareLattice(15)
        fld = LandauField(0.2)
        PBC = false
        H1 = @hamiltonian begin
            lattice := l
            @hop axis = 1 pbc = PBC
        end
        @test H1 == TightBinding(l, pbc=PBC)

        l = HoneycombLattice(5, 5)
        x, y = coord_values(l)
        H2 = @hamiltonian begin
            lattice := l
            field := fld
            @diag site -> site.x + exp(site.y)
            @hop site_indices = (2, 1)
            @hop site_indices = (2, 1) axis = 1
            @hop site_indices = (2, 1) axis = 2
        end
        lv = @. x + exp(y)
        @test H2 == TightBinding(lv, field=fld)
        @test H2 == TightBinding(nothing, lv, field=fld)

        l = SquareLattice(10, 10)
        x, y = coord_values(l)
        sel = DomainsSelector(x .< 0)
        σ = [[0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
        H3 = @hamiltonian begin
            lattice := l
            field := fld
            dims_internal := 2
            sparse := true
            @diag σ[3]
            @hop (σ[3] - im * σ[1]) / 2 axis = 1 sel
            @hop (σ[3] - im * σ[2]) / 2 axis = 2 sel
        end
        @test H3 == SpinTightBinding(sel, l, field=fld)
        @test H3 == SpinTightBinding(sel, ones(l), field=fld)
        @test_throws "no method" SpinTightBinding(ones(l), 2)
    end

    @testset "DOS & LDOS" begin
        sp = spectrum(SpinTightBinding(ones(SquareLattice(10, 10))))
        E = 2
        δ = 0.2
        Es = eigvals(sp)
        Vs = eigvecs(sp)
        ld1 = imag.(diag_reduce(tr, LatticeModels.LatticeArray(basis(sp), Vs * (@.(1 / (Es - E - im * δ)) .* Vs'))))
        ldosf = ldos(sp, δ)
        @test ldos(sp, E, δ).values ≈ ld1.values
        @test ldosf(E).values ≈ ld1.values
        @test dos(sp, δ)(E) ≈ sum(ld1)
    end
end

@testset "Currents" begin
    l = SquareLattice(10, 10)
    x, y = coord_values(l)
    H(B) = SpinTightBinding(l, field=LandauField(B))
    P = filled_projector(spectrum(H(0)))
    dc = DensityCurrents(H(0.1), P)
    bs = bonds(H(0.1))
    m1 = materialize(dc)[x.<y]
    m2 = materialize(dc[x.<y])
    m3 = materialize(bs, dc)[x.<y]
    m4 = materialize(bs, dc[x.<y])
    md = materialize(pairs_by_distance(≤(2)), dc[x.<y])
    @test m1.currents == m2.currents
    @test m1.currents == m3.currents
    @test m1.currents == m4.currents
    @test m1.currents == md.currents
    @test (m1 + m2 - m3).currents ≈ (2 * m4 * 1 - md / 1).currents
    s1 = l[23]
    s2 = l[45]
    @test dc[s1, s2] == -dc[s2, s1]
    @test abs(dc[s1, s1]) < eps()
end
