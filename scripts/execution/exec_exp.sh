# bash script for orchestrated execution of experiments

THREADS=4
LOG="data/output.log"

julia --project=. -t "$THREADS" -e '
include("src/DQCode.jl")
using .DQCode

# --- one-time data generation ---
for exp_label in [
    "steane_4_3",
    "shor_3_3_3",
    "trivariate_3_3_3_3",
    "trivariate_4_4_4",
    "trivariate_6_6",
    "color_8_9",
    "bivariate_3_3_3_3_3_3",
    "bivariate_6_6_6",
    "bivariate_9_9",
]
    DQCode.create_code_network_data(exp_label)
end

# --- Qiskit baseline ---
for exp_label in [
    "steane_4_3",
    "shor_3_3_3",
    "trivariate_3_3_3_3",
    "trivariate_4_4_4",
    "trivariate_6_6",
    "color_8_9",
    "bivariate_3_3_3_3_3_3",
    "bivariate_6_6_6",
    "bivariate_9_9",
]
    DQCode.baseline_encoding_qiskit(exp_label)
end

# --- MQT Encoding ---
for exp_label in [
    "steane_4_3",
    "shor_3_3_3",
]
    DQCode.baseline_encoding_mqt(exp_label, DQCode.MQT_PATH, "optimal") 
end

# --- MQT Encoding ---
for exp_label in [
    "trivariate_3_3_3_3",
    "trivariate_4_4_4",
    "trivariate_6_6",
    "color_8_9",
    "bivariate_3_3_3_3_3_3",
    "bivariate_6_6_6",
    "bivariate_9_9",
]
    DQCode.baseline_encoding_mqt(exp_label, DQCode.MQT_PATH, "heuristic") 
end


# --- Gottesman warm_start Genetic Algorithm ---
for exp_label in [
    "steane_4_3",
    "shor_3_3_3",
    "trivariate_3_3_3_3",
    "trivariate_4_4_4",
    "trivariate_6_6",
    "color_8_9",
    "bivariate_3_3_3_3_3_3",
    "bivariate_6_6_6",
    "bivariate_9_9",
]
    DQCode.circuit_search_gott_ga(exp_label)
end

# --- Monte Carlo Tree Search ---
for exp_label in [
    "steane_4_3",
    "shor_3_3_3",
    "trivariate_3_3_3_3",
    "trivariate_4_4_4",
    "trivariate_6_6",
    "color_8_9",
    "bivariate_3_3_3_3_3_3",
    "bivariate_6_6_6",
    "bivariate_9_9",
]
    DQCode.circuit_search_MCTS(exp_label)
end
' #2>&1 | tee "$LOG"