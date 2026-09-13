include(joinpath(@__DIR__, "..", "..", "src", "DQCode.jl"))
using .DQCode
using BenchmarkTools
using Statistics

mqt_path = normpath(joinpath(@__DIR__, "..", "..", "..", "qecc")) # assumes that PyCall was build correspondingly (see setup instructions in README.md)
isdir(mqt_path) || error("Cannot find qecc at $mqt_path. Clone qecc next to DQCode, or set mqt_path manually in this script.")

exp_label = "trivariate_6_6" # available configurations are stored and can be adapted in src/experiment/config.jl

circuit_path = "warmstart_ga/GA_circuit.jls" # circuit path (within {code}/{architecture}/ folder) pointing to the circit to simulate
method = "optimal" 

ps = [0.000175]
p_bells = [0.0009]
telegate_idle_depth = 12                                       
p_single_ratio = 1/100                                          
p_idle_ratio = 1/10     
samples =  Int.(floor.(10 .^ range(log10(1e4),log10(5e6),length=20)))
run_times = Dict{Int, Tuple{Float64, Float64}}()
for num_samples in samples
    benchmark_samples = @benchmarkable DQCode.dqc_simulation($exp_label, $mqt_path, $circuit_path, Int($num_samples), $ps, $p_bells,
                      $telegate_idle_depth, $p_single_ratio, $p_idle_ratio, $method) samples = 3 seconds=1000 evals=1 
    b = run(benchmark_samples)
    run_times[num_samples] = (mean(b.times./1e9), std(b.times./1e9)) 
    display(b)
    println(run_times)
end
df = DataFrame(
    num_samples = collect(keys(run_times)),
    mean_time_s = [run_times[k][1] for k in keys(run_times)],
    std_time_s  = [run_times[k][2] for k in keys(run_times)]
)
CSV.write("data/dqc_sim_scaling.csv", df)
