println("\nBenchmark: Hamiltonian creation")

function measure_lm(n)
    @elapsed begin
        make_hamiltonian_lm(n)
    end
end

@pyexec ```
def measure_kwant(num_sites):
    with pb.utils.timed() as time:
        m = make_model_kwant(num_sites)
        h = m.hamiltonian_submatrix(sparse=True)
    return time.elapsed
``` => measure_kwant

@pyexec ```
def measure_pb(num_sites):
    with pb.utils.timed() as time:
        m = make_model_pb(num_sites)
        h = m.hamiltonian
    return time.elapsed``` => measure_pb

Ns = 10 .^ (2:7)
benchmark_results = Dict{String, Any}()
benchmark_results["LatticeModels"] = run_benchmarks(Ns, measure_lm)
benchmark_results["Kwant"] = run_benchmarks(Ns, measure_kwant)
benchmark_results["Pybinding"] = run_benchmarks(Ns, measure_pb)
add_benchmark!("Hamiltonian creation", Ns, benchmark_results)
