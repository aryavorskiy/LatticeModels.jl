include("common.jl")

function measure_lm_ev_const(n)
    H = make_hamiltonian_lm(n)
    psi = groundstate(H)
    @elapsed begin
        ev = Evolution(H, psi, timedomain=0:0.1:10, showprogress=false)
        dens = TimeSequence(ev) do moment
            localdensity(moment.state)
        end
    end
end

@pyexec ```
def measure_kwant_ev_const(num_sites):
    syst = make_model_kwant(num_sites)
    h_mat = syst.hamiltonian_submatrix(sparse=True)
    psi = sla.eigsh(h_mat, k=1, which='SA')[1][:, 0]
    dens = []
    density_op = kwant.operator.Density(syst)
    with pb.utils.timed() as time:
        state = tkwant.onebody.WaveFunction(h_mat, None, psi)
        for t in np.arange(0, 10.1, 0.1):
            state.evolve(t)
            dens.append(state.evaluate(density_op))
    return time.elapsed
``` => measure_kwant_ev_const

Ns = round.(Int, 10 .^ (2:0.5:5))
@isdefined(ev_benchmark_results) || (ev_benchmark_results = Dict{String, Any}())
ev_benchmark_results["LatticeModels"] = run_benchmarks(Ns, measure_lm_ev_const)
ev_benchmark_results["Kwant"] = run_benchmarks(Ns, measure_kwant_ev_const)

f, ax = plot_benchmarks(Ns, ev_benchmark_results)
ax.title = "Evolution, constant Hamiltonian"
save("benchmark_evolution_const.svg", f)
