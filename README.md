# Biomolecular flux
Data repository for "Non-equilibrium modeling of directed flux through biomolecular condensates" manuscript.

## Organization
* **Figures**<br>
  Directory containing all the data files required for reproducing the main manuscript figures -- Flux times, mean squared displacements, pore size distributions, radius of gyration, contact maps.

* **force-field**<br>
  Force field parameters compatibale with LAMMPS for the flux simulations with native mpipi parameters `potential_60_particle_types.dat`, and reduced interactions `potential_60_particle_types_handoff_mech_weaker_G.dat`. Vanilla mpipi simulations (20 particle types) can be performed with the `potential.dat` force field file. 
  
* **scripts**<br>
  Necessary script for performing flux simulations.<br>
  `equilibrate.lmp` -- LAMMPS script to equilibrate host-guest system where all the guest molecules are inside the active site. `run_equilibration.sh` is the slurm submission script for the LAMMPS equilibration script.<br>
  `flux_run/lmp` -- LAMMPS script to perform flux simulation starting from a equilibrated configuration. `submit.slurm` is the submission script for the LAMMPS flux_run script.<br>
  `generate_slabs.py` -- Python script to create an initial config for any combination of host and guest molecules.<br>
  `chain_insertion-deletion_v2.py` -- Python script to execute the insertion/deletion step of the flux simulation protocol. `chain_insertion-deletion_v2_with_angles.py` does the same but for guest molecules with nonzero stiffness.<br>
  `repopulate_cavity.py` -- Python script to repopulate the active site with additional guest molecules after equilibration IF the number of guest molecules in the active site is smaller than the target.<br>
  `simulate.sh` -- Bash script to automate the flux simulations. The script recursively executes the flux and insertion/deletion steps (in a sequence).

* **template-configs**<br>
  Template configurations for all the host proteins with empty active sites. The active sites can be populated with any required guest molecules and used to perform the flux simulations.
  
  