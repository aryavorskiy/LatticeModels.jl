using CondaPkg
CondaPkg.add(["numpy", "kwant", "tkwant", "pybinding"], channel="conda-forge")

using PythonCall, CairoMakie
include("models.jl")

pyexec("""
import threading
import ctypes

def _async_raise(tid, exctype):
    '''Abort the thread'''
    if not isinstance(tid, int):
        tid = tid.ident
    res = ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), ctypes.py_object(exctype))
    if res == 0:
        raise ValueError("Invalid thread id")
    elif res > 1:
        # If it returns a number greater than one, we're in trouble, and we try to revert the effect
        ctypes.pythonapi.PyThreadState_SetAsyncExc(ctypes.c_long(tid), None)
        raise SystemError("PyThreadState_SetAsyncExc failed")""", Main)

MAXTIME_DEFAULT = 360

function time_function(f::Function, n; maxtime=MAXTIME_DEFAULT)
    c = Channel{Float64}() do ch
        sleep(maxtime)
        put!(ch, NaN)
    end
    @async begin    # Worker
        isopen(c) && put!(c, f(n))
    end
    return take!(c)
end

@pyexec ```
result = float('NaN')
def py_measure_time(func, maxtime, *args):
    def wrapper():
        global result
        result = func(*args)

    thread = threading.Thread(target=wrapper)
    thread.start()
    thread.join(timeout=maxtime)
    if thread.is_alive():
        _async_raise(thread, SystemExit)

    return result``` => py_measure_time
time_function(f::Py, n; maxtime=MAXTIME_DEFAULT) =
    pyconvert(Number, py_measure_time(f, maxtime, n))

function run_benchmarks(Ns, f; nrepeat=5, kw...)
    println("Benchmarking $(f)...")
    f(first(Ns)) # warmup
    Ts = Float64[]
    print("[")
    for n in Ns
        t = time_function(f, n; kw...)
        for _ in 2:nrepeat
            (t < 1) && (t = time_function(f, n; kw...))
        end
        print(round(t, sigdigits=5))
        n == Ns[end] || print(", ")
        push!(Ts, t)
    end
    println("]")
    return Ts
end

@isdefined(benchmark_names) || (benchmark_names = String[])
@isdefined(benchmark_results_array) || (benchmark_results_array = [])
function add_benchmark!(name, ns, results)
    i = findfirst(==(name), benchmark_names)
    if i === nothing
        push!(benchmark_names, name)
        push!(benchmark_results_array, (ns, results))
    else
        benchmark_results_array[i] = (ns, results)
    end
end
function plot_benchmarks(;kw...)
    fig = Figure(;kw...)
    n = length(benchmark_names)
    for i in eachindex(benchmark_names)
        name = benchmark_names[i]
        ns, results = benchmark_results_array[i]
        rowi = (i - 1) ÷ 2 + 1
        coli = (i - 1) % 2 + 1

        ax = Axis(fig,
            xscale = log10,
            xlabel = "Number of sites",
            yscale = log10,
            ylabel = "Time [s]",
            yminorticks = IntervalsBetween(5),
            yminorticksvisible = true,
            yminorgridvisible = true,
            yminorgridstyle = :dash,
            alignmode = Outside(),
            title = name
        )
        for key in ["LatticeModels", "Kwant", "Pybinding"]
            key in keys(results) || continue
            value = results[key]
            lines!(ax, ns, value, label=key)
            scatter!(ax, ns, value, label=key)
        end
        axislegend(ax, position=:rb, merge=true)

        if i == n && isodd(n) && n != 1
            fig[rowi, 2:3] = ax
        elseif coli == 1
            fig[rowi, 1:2] = ax
        else
            fig[rowi, 3:4] = ax
        end
    end
    fig
end
