using PythonCall, CondaPkg, CairoMakie
include("models.jl")

CondaPkg.add("kwant")
CondaPkg.add("tkwant")
CondaPkg.add("pybinding", channel="conda-forge")

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
