#!/bin/bash

# Check if the number of processes is provided as a command-line argument
if [ -z "$1" ]; then
  echo "Usage: $0 <number_of_processes>"
  exit 1
fi

# Get the number of processes from the command-line argument
num_procs=$1

# Define the values for a, x, and N
a_values=(0.1 0.2 0.3 0.4 0.5)   # Example values for a
x_values=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0) # Example values for x
N_values=(100 200 300 400 500 600 700 800 900 1000) # Example values for N
s_values=(1 2 3 4 5)  # Example values for seed

# Define variables
array_range=(1 2)  # Array of indices for the task
base_dir="/home/rnp/Software/MOOSE/projects/panther/problems"  # Set your base directory path
template_dir="$base_dir/template"
output_dir="$base_dir/output_dump"

# Automatically generate the array_range based on the number of values for a, x, N, and s
num_a=${#a_values[@]}  # Number of values for a
num_x=${#x_values[@]}  # Number of values for x
num_N=${#N_values[@]}  # Number of values for N
num_s=${#s_values[@]}  # Number of values for s

array_range=$((num_a * num_x * num_N * num_s))

# Loop over all combinations
for idx in $(seq 0 $((array_range - 1))); do
  # Calculate the index for a, x, N, and s
  a_idx=$((idx / (num_x * num_N * num_s)))       # Index for a
  x_idx=$((idx / (num_N * num_s) % num_x))       # Index for x
  N_idx=$(((idx / num_s) % num_N))               # Index for N
  s_idx=$((idx % num_s))                         # Index for s

  # Get the actual values for a, x, N, and s
  a=${a_values[$a_idx]}
  x=${x_values[$x_idx]}
  N=${N_values[$N_idx]}
  s=${s_values[$s_idx]}

  echo "Running simulation for: a=$a, x=$x, N=$N, s=$s"

  # Create a unique directory for this combination
  task_dir="$output_dir/DBC_$((idx+1))_a${a}_x${x}_N${N}_s${s}"
  mkdir -p "$task_dir"
  cd "$task_dir" || exit 1

  echo "Deleting all previous files in the folder..."
  rm -rf ./*  # Remove all files

  echo "Copying starter files from the Template directory..."
  cp "$template_dir"/* .

  echo "Updating parameters..."
  sed -i "s/a_val/$a/" DBC_split_nd.i
  sed -i "s/x_val/$x/" DBC_split_nd.i
  sed -i "s/N_val/$N/" DBC_split_nd.i
  sed -i "s/N_val/$s/" DBC_split_nd.i

  # Log start time
  start_time=$(date)
  echo "Simulation started at: $start_time"
  
  # Define output log files
  log_file="$task_dir/output_${idx}_a${a}_x${x}_N${N}_s${s}.log"
  err_file="$task_dir/error_${idx}_a${a}_x${x}_N${N}_s${s}.log"

  echo "Running the simulation with $num_procs processes..."
  time mpiexec -np "$num_procs" ../../../panther-opt -i DBC_split_nd.i > "$log_file" 2> "$err_file"

  echo "Simulation task number $((idx+1)) complete. Output saved to $log_file and errors to $err_file."
done

echo "All tasks completed."
