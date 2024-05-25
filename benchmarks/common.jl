using PythonCall, CondaPkg, CairoMakie

# CondaPkg.update()
CondaPkg.add("kwant")
CondaPkg.add("tkwant")
CondaPkg.add("pybinding", channel="conda-forge")

_pyconvert(T, py::Py) = pyconvert(T, py)
_pyconvert(T, not_py) = convert(T, not_py)
function run_benchmarks(Ns, f; nrepeat=2)
    println("Benchmarking $(f)...")
    f(first(Ns)) # warmup
    Ts = [minimum(_pyconvert(Number, f(n)) for _ in 1:nrepeat) for n in Ns]
    println(round.(Ts, sigdigits=5))
    return Ts
end

function plot_benchmarks(Ns, benchmark_results)
    f = Figure()
    ax = Axis(f[1,1],
        xscale = log10,
        xlabel = "Number of sites",
        yscale = log10,
        ylabel = "Time [s]",
        yminorticks = IntervalsBetween(5),
        yminorticksvisible = true,
        yminorgridvisible = true,
        yminorgridstyle = :dash
    )
    for key in ["LatticeModels", "Kwant", "Pybinding"]
        key in keys(benchmark_results) || continue
        value = benchmark_results[key]
        lines!(ax, Ns, value, label=key)
        scatter!(ax, Ns, value, label=key)
    end
    axislegend(ax, position=:rb, merge=true)

    return f, ax
end
