println("\nBenchmark: Ground state density")

function measure_lm_groundstate(n)
    H = make_hamiltonian_lm(n)
    @elapsed begin
        psi = groundstate(H)
        d = localdensity(psi)
    end
end

@pyexec """
def measure_kwant_groundstate(num_sites):
    h_mat = make_model_kwant(num_sites).hamiltonian_submatrix(sparse=True)
    with pb.utils.timed() as time:
        psi = sla.eigsh(h_mat, k=1, which='SA')[1]
        d = np.abs(psi)**2
    return time.elapsed
""" => measure_kwant_groundstate

@pyexec """
def measure_pb_groundstate(num_sites):
    h = make_model_pb(num_sites)
    with pb.utils.timed() as time:
        solv = pb.solver.arpack(h, k=1, which='SA', sigma=None)
        d = solv.calc_probability(0)
    return time.elapsed
""" => measure_pb_groundstate

Ns = round.(Int, 10 .^ (2:0.5:5))
gs_benchmark_results = Dict{String, Any}()
gs_benchmark_results["LatticeModels"] = run_benchmarks(Ns, measure_lm_groundstate)
gs_benchmark_results["Kwant"] = run_benchmarks(Ns, measure_kwant_groundstate)
gs_benchmark_results["Pybinding"] = run_benchmarks(Ns, measure_pb_groundstate)
add_benchmark!("Ground state density computation", Ns, gs_benchmark_results)
