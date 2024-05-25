include("common.jl")
using LatticeModels

function measure_lm(n, magnetic_field=3.0)
    @elapsed begin
        shape = LatticeModels.Circle()
        lat = HoneycombLattice(shape, sites=n)
        r = shaperadius(lat, shape)
        pn_junction_potential = LatticeValue(lat) do (x, y)
            x^2 + y^2 < r^2/4 ? 0.2 : 0.0
        end
        field = SymmetricGauge(magnetic_field * 1.519e-3)
        tightbinding_hamiltonian(addlookuptable(lat), t1=-2.8, pn_junction_potential, field=field)
    end
end

pyexec("""
import numpy as np
import kwant, cmath, math
import pybinding as pb
from pybinding.repository import graphene
def calc_radius(num_sites, lattice=graphene.monolayer()):
    unit_area = np.linalg.norm(np.cross(*lattice.vectors)) / len(lattice.sublattices)
    return math.sqrt(num_sites * unit_area / math.pi)
""", Main)

@pyexec ```
def measure_kwant(num_sites):
    def make_system(radius, potential=0.2, magnetic_field=3):
        def circle(pos):
            x, y = pos
            return x**2 + y**2 < radius**2

        def onsite(site):
            x, y = site.pos
            if x**2 + y**2 < (radius / 2)**2:
                return potential
            else:
                return 0

        def hopping(site_i, site_j):
            xi, yi = site_i.pos
            xj, yj = site_j.pos
            phi_ij = 0.5 * magnetic_field * (xi + xj) * (yi - yj)
            const = 1.519e-3  # 2*pi*e/h*[nm^2]
            return -2.8 * cmath.exp(1j * phi_ij * const)

        a_cc = 0.142
        a = a_cc * math.sqrt(3)
        graphene_lattice = kwant.lattice.general([(a, 0), (a/2, a/2 * math.sqrt(3))],
                                                 [(0, -a_cc/2), (0, a_cc/2)])
        builder = kwant.Builder()
        builder[graphene_lattice.shape(circle, (0, 0))] = onsite
        builder[graphene_lattice.neighbors()] = hopping
        builder.eradicate_dangling()
        return builder.finalized()

    with pb.utils.timed() as time:
        system = make_system(radius=calc_radius(num_sites))
        h = system.hamiltonian_submatrix(sparse=True)

    return time.elapsed
``` => measure_kwant

@pyexec ```
def measure_pybinding(num_sites, warmup=False, plot=False):
    def circular_pn_junction(v0, radius):
        @pb.onsite_energy_modifier
        def function(energy, x, y):
            energy[(x**2 + y**2) < (radius / 2)**2] = v0
            return energy
        return function

    def make_model(radius, potential=0.2, magnetic_field=3):
        return pb.Model(
            graphene.monolayer().with_min_neighbors(2),
            pb.circle(radius),
            circular_pn_junction(potential, radius),
            graphene.constant_magnetic_field(magnetic_field)
        )

    with pb.utils.timed() as time:
        model = make_model(radius=calc_radius(num_sites))
        h = model.hamiltonian

    return time.elapsed``` => measure_pybinding

Ns = 10 .^ (2:7)
benchmark_results = Dict{String, Any}()
benchmark_results["LatticeModels"] = run_benchmarks(Ns, measure_lm)
benchmark_results["Kwant"] = run_benchmarks(Ns, measure_kwant)
benchmark_results["Pybinding"] = run_benchmarks(Ns, measure_pybinding)

f, ax = plot_benchmarks(Ns, benchmark_results)
ax.title = "Hamiltonian creation"
save("benchmark_hamiltonian.svg", f)
