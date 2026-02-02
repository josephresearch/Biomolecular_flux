#!/bin/bash

templatesPath="/scratch/gpfs/JERELLE/yashraj/RNA-flux-project/systematic-study/synthetic-sequences/flux-simulations/Templates"


dir="production"
itr=25

# Define the host systems
hostSystems=("Y34S116" "Y35S115" "Y36S114" "Y37S113" "Y38S112" "Y1S3-36_S6" "Y3S9-12_S6" "Y4S12-9_S6")
nhostChains=(388 398 424 440 454 424 400 388)

# guestChains=("G10" "G30" "G50")
guestChains=("G20")


temp=-1

# Iterate over each host system
for host in "${hostSystems[@]}"; do
    echo "Generating initial configurations for production simulations for host: "${host}

    # Temperature for the current host
    temperature=310

    temp=$((temp + 1))
    nhost=${nhostChains[$temp]}

    # Change to the host directory
    cd host_${host}

    temp1=0
    # Iterate over each guest chain
    for guest in "${guestChains[@]}"; do
        echo "Processing guest chain: ${guest}"
    
        # Create a directory for the guest chain
        cd guest_${guest}

        nguest=${nguestChains[$temp1]}

        echo "Temperature: " ${temperature}
        echo "# host chains: " ${nhost}
        echo "# guest chains: "${nguest}

        # Create directory for production simulations
        mkdir -p ${dir}
        cd ${dir}

        # Copy and make appropriate changes to the simulation files
        cp ${templatesPath}/simulate.sh .
        # Change the number of iterations in the simulate.sh script
        sed -i "s/ITERATIONS=100/ITERATIONS=${itr}/" simulate.sh

        # Copy LAMMPS input file
        cp ${templatesPath}/flux_run.lmp .
        # Correct the LAMMPS input file
        sed -i "1s/.*/variable temperature equal ${temperature}/" flux_run.lmp
        sed -i "2s/.*/variable randomSeed equal ${RANDOM}/" flux_run.lmp
        sed -i "3s/.*/variable numHostChains equal ${nhost}/" flux_run.lmp
        
        # Copy other necessary files
        cp ${templatesPath}/monte_carlo_chain_insertion_v2.py .
        cp ${templatesPath}/submit.slurm .
        sed -i "2s/.*/#SBATCH --job-name=${host}_${guest}/" submit.slurm
        cp ${templatesPath}/repopulate_cavity.py .

        # Repopulate the cavity
        python3 repopulate_cavity.py \
        --trajectory ../equilibrate/equilibrate.lammpstrj \
        --configuration ../equilibrate/initial_config.dat \
        --target_num_guest_chains 25 \
        --logfile repopulate_cavity_logs.txt \
        --cavity_bounds 950.0 1050.0

        # Copy the updated configuration file
        cp equilibrated_config_after_repopulation.dat updated_config.dat

        cd .. # exit production directory

        cd .. # exit guest directory
        
        temp1=$((temp1 + 1))

    done

    cd .. # exit host directory

done
