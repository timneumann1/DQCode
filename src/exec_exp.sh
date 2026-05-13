# bash script for orchestrated execution of experiments

THREADS=4
PROJECT="envs/opt"
LOG="src/results/output.log"

julia --project="$PROJECT" -t "$THREADS" -e '
include("src/DQCircuitSearch.jl")
using .DQCircuitSearch

# --- one-time data generation ---
for exp_label in [
    #"steane_4_3",
    #"shor_3_3_3",
    #"trivariate_3_3_3_3",
    # "trivariate_4_4_4",
    # "trivariate_6_6",
    # "bivariate_3_3_3_3_3_3",
    # "bivariate_6_6_6",
    # "bivariate_9_9",
]
    DQCircuitSearch.create_code_network_data(exp_label)
end

# --- Qiskit baseline ---
for exp_label in [
    # "steane_4_3",
    # "shor_3_3_3",
    # "trivariate_3_3_3_3",
    # "trivariate_4_4_4",
    # "trivariate_6_6",
    # "bivariate_3_3_3_3_3_3",
    # "bivariate_6_6_6",
    # "bivariate_9_9",
]
    DQCircuitSearch.baseline_encoding_qiskit(exp_label)
end

# --- MQT Encoding ---
for exp_label in [
    #"steane_4_3",
    #"shor_3_3_3",
    # "trivariate_3_3_3_3",
    # "trivariate_4_4_4",
    # "trivariate_6_6",
    # "bivariate_3_3_3_3_3_3",
    # "bivariate_6_6_6",
    # "bivariate_9_9",
]
    DQCircuitSearch.baseline_encoding_mqt(exp_label)
end

# --- Gottesman baseline (+ GA inside circuit_search_gott) ---
for exp_label in [
    #"steane_4_3",
    #"shor_3_3_3",
    "trivariate_3_3_3_3",
    "trivariate_4_4_4",
    "trivariate_6_6",
    "bivariate_3_3_3_3_3_3",
    "bivariate_6_6_6",
    "bivariate_9_9",
]
    DQCircuitSearch.circuit_search_gott(exp_label)
end

# --- Monte Carlo Tree Search ---
for exp_label in [
    "steane_4_3",
    # "shor_3_3_3",
    # "trivariate_3_3_3_3",
    # "trivariate_4_4_4",
    # "trivariate_6_6",
    # "bivariate_3_3_3_3_3_3",
    # "bivariate_6_6_6",
    # "bivariate_9_9",
]
    DQCircuitSearch.circuit_search_MCTS(exp_label)
end
' 2>&1 | tee "$LOG"