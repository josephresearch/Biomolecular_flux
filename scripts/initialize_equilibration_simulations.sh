#!/bin/bash

templatesPath="/scratch/gpfs/JERELLE/yashraj/RNA-flux-project/systematic-study/synthetic-sequences/flux-simulations/Templates"

dir="equilibrate"

# Define the host systems
hostSystems=("Y34S116" "Y35S115" "Y36S114" "Y37S113" "Y38S112" "Y1S3-36_S6" "Y3S9-12_S6" "Y4S12-9_S6")
nhostChains=(388 398 424 440 454 424 400 388)

# Define the guest chains and their respective counts
# guestChains=("G10" "G30" "G50")
# seqs=(GGGGGGGGGG GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG)
# nguestChains=(100 65 40)
guestChains=("G20")
seqs=(GGGGGGGGGGGGGGGGGGGG)
nguestChains=(65)


temp=0

# Iterate over each host system
for host in "${hostSystems[@]}"; do
    echo "Generating initial configurations for host: "${host}

    # Temperature for the current host
    temperature=310
    
    nhost=${nhostChains[$temp]}

    # Change to the host directory
    mkdir -p host_${host}
    cd host_${host}

    temp1=0
    # Iterate over each guest chain
    for guest in "${guestChains[@]}"; do
        echo "Processing guest chain: ${guest}"
    
        # Create a directory for the guest chain
        mkdir -p guest_${guest}
        cd guest_${guest}

        nguest=${nguestChains[$temp1]}
        guest_seq=${seqs[$temp1]}

        echo "Temperature: " ${temperature}
        echo "# host chains: " ${nhost}
        echo "# guest chains: "${nguest}

        # Create directory for equilibration simulations
        mkdir -p ${dir}
        cd ${dir}

        # Copy and make appropriate changes to the simulation files
        cp ${templatesPath}/equilibrate.lmp .
        sed -i "2s/.*/variable randomSeed equal $RANDOM/" equilibrate.lmp
        sed -i "3s/.*/variable numHostChains equal ${nhost}/" equilibrate.lmp

        cp ${templatesPath}/run_equilibration.slurm .
        
        # Create the initial configuration
        python3 ${templatesPath}/generate_slabs.py \
        --option cavity+protein1+protein2+peptides \
        --simBox 0.0 200.0 0.0 200.0 0.0 2000.0 \
        --cavity_prot1_and_prot2_config ${templatesPath}/empty_cavity_HOST_${host}_NHOST_${nhost}/equilbrated_initial_config.dat \
        --cavity_prot1_and_prot2_traj ${templatesPath}/empty_cavity_HOST_${host}_NHOST_${nhost}/equilibrate.lammpstrj \
        --cavity_peptides_seq ${guest_seq} \
        --cavity_peptides_nchains ${nguest}

        mv config_with_peptides_inside_cavity.dat initial_config.dat

        # Submit the equilibration job
        sbatch run_equilibration.slurm

        cd .. # exit equilibration directory

        cd .. # exit guest directory
        
        temp1=$((temp1 + 1))

    done

    temp=$((temp + 1))
    cd .. # exit host directory

done
