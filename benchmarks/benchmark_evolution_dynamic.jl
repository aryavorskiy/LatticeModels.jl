println("\nBenchmark: Evolution, time-dependent Hamiltonian")

function measure_lm_evd_const(n)
    H = make_hamiltonian_lm(n)
    psi = groundstate(H)
    @elapsed begin
        ev = Evolution(KrylovKitExp(), t -> make_hamiltonian_lm(n, 3.0 * min(t / 5, 1)), psi,
            timedomain=0:0.1:10, showprogress=false)
        dens = TimeSequence(ev) do moment
            localdensity(moment.state)
        end
    end
end

@pyexec ```
def measure_kwant_evd_const(num_sites):
    syst = make_model_tkwant(num_sites)
    h_mat = syst.hamiltonian_submatrix(sparse=True, params={'time': 0})
    eigval, eigvec, *info = sla.eigsh(h_mat, k=1, which='SA')
    dens = []
    density_op = kwant.operator.Density(syst)
    with pb.utils.timed() as time:
        state = tkwant.onebody.WaveFunction.from_kwant(syst=syst,
            psi_init=eigvec[:, 0], energy=eigval[0])
        for t in np.arange(0, 10.1, 0.1):
            state.evolve(t)
            dens.append(state.evaluate(density_op))
    return time.elapsed
``` => measure_kwant_evd_const

Ns = round.(Int, 10 .^ (2:0.5:5))
evd_benchmark_results = Dict{String, Any}()
evd_benchmark_results["LatticeModels"] = run_benchmarks(Ns, measure_lm_evd_const)
evd_benchmark_results["Kwant"] = run_benchmarks(Ns, measure_kwant_evd_const)
add_benchmark!("Evolution, time-dependent Hamiltonian", Ns, evd_benchmark_results)
