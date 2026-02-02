#!/bin/bash

# Number of iterations
ITERATIONS=100

# SLURM submission script
SLURM_SCRIPT="submit.slurm"

# Python script to call
PYTHON_SCRIPT="monte_carlo_chain_insertion_v2.py"

for ((i=1; i<=ITERATIONS; i++))
do
    echo "Starting iteration $i..."

    # Submit the SLURM job and capture the job ID
    JOB_ID=$(sbatch $SLURM_SCRIPT | awk '{print $4}')
    echo "Submitted SLURM job with ID $JOB_ID"

    # Wait for the SLURM job to finish
    echo "Waiting for SLURM job $JOB_ID to complete..."
    while [ $(squeue -j $JOB_ID | wc -l) -gt 1 ]; do
        sleep 300
    done
    
    echo "SLURM job $JOB_ID completed."

    # Make a copy of the configuration after the LAMMPS job
    cp final_config.dat final_config_${i}.dat
    echo "Configuration file copied to final_config_${i}.dat"

    # Run the Python script and wait for it to finish
    echo "Running MC insertion/deletion script..."
    python3 $PYTHON_SCRIPT \
    --trajectory all.lammpstrj \
    --configuration final_config.dat \
    --logfile MC_insertion_logs.txt \
    --cavity_bounds 950.0 1050.0

    # Make a copy of the updated config file after the Python script
    cp updated_config.dat updated_config_${i}.dat

    # Make copy of the density plot
    cp interface_location.pdf interface_loc_${i}.pdf

    echo "Iteration $i completed."
done

echo "All iterations completed successfully."