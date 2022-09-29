using Test, LinearAlgebra
using LatticeModels

@testset "Workflows" begin
    @test begin
        l = SquareLattice(10, 10)
        d = diag_operator(l) do site; [1;;] end
        hx = hopping_operator(l, Hopping(axis=1))
        hy = hopping_operator(l, Hopping(axis=2))
        H = @on_lattice d + hx + hy
        sp = Spectrum(H)
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
        P = filled_projector(Spectrum(H))
        X, Y = coord_operators(Basis(l, 2))
        d = diag_aggregate(tr, @on_lattice 4π * im * P * X * (I - P) * Y * P)
        da = @on_lattice d .|> real
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
        P0 = filled_projector(Spectrum(H))
        X, Y = coord_operators(Basis(l, 2))
        @evolution {P0 => h(t) => P} for t in 0:0.1:10
            d = diag_aggregate(tr, 4π * im * P * X * (I - P) * Y * P)
            da = @on_lattice d .|> real
        end
        true
    end
end
