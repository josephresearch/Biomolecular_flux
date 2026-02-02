import numpy as np
import sys
sys.path.insert(1, '/home/yw9071/scripts')
sys.path.insert(1, '/home/yw9071/scripts/RNA_flux_project')
from RNA_flux_project import utils
import matplotlib.pyplot as plt
plt.style.use("/home/yw9071/scripts/plotting/joseph_group.mplstyle")
import argparse


"""
Description: 
This script reads the trajectory and configuration files from the equilibration simulations
and then
(1) deletes any guest chains in the dilute regime and adds them inside the cavity.
(2) Also, if the total number of guest chians inside the cavity is less than the 
    target number, then we add additional guest chains to the cavity.

    
Algorithm:
    1. Read the trajectory and configuration files from the equilibraiton simulations.
    2. Determine the interface location (z-axis) on either side of the cavity.
    3. Iterate over every guest chain and check if it is in the dilute phase.
    4. If it is in the dilute phase, then delete it and add it inside the cavity.
    5. If the total number of guest chains inside the cavity is less than the target number,
       then add additional guest chains to the cavity.
    6. Write the updated configuration file and save the insertion/deletion information 
       in a text file.

       
Input parameters:
    1. LAMMPS trajectory file
    2. LAMMPS config file
    3. Target number of guest chains inside the cavity
    4. TEXT file to save insertion/deletion information
    5. Cavity bounds (z-axis) for the insertion of guest chains
"""

# Read input parameters using argparse
parser = argparse.ArgumentParser()
parser.add_argument("--trajectory", type=str, required=True)
parser.add_argument("--configuration", type=str, required=True)
parser.add_argument("--target_num_guest_chains", type=int, required=True)
parser.add_argument("--logfile", type=str, required=True)
parser.add_argument("--cavity_bounds", type=float, nargs=2, required=False)
args = parser.parse_args()

# .................................................................................
# Read the trajectory file to get the positions of the guest chains
print(f"Reading the trajectory file: {args.trajectory}")
# NOTE: UNWRAPPED COORDINATES EXPECTED
simulation_box, coords, timestep = utils.read_trajectory_last_frame(args.trajectory)
# coords = (atom_id, mol_id, atom_type, charge, xu, yu, zu)

print(f"Reading the configuration file: {args.configuration}\n")
num_atoms, num_bonds, bond_data = utils.read_config(args.configuration)
bonds = np.array(bond_data) # (bond_id, bond_type, atom1, atom2)

# Filter the positions of guest and host chains
""" Monomers  |  Atom types
------------------------------
    Host      |     1-20
    Cavity    |     21-40
    Guest     |     41-60       """

mask_host = (coords[:,2] >= 1) & (coords[:,2] <= 20) # column 3 (index 2) contains atom type info.
mask_cavity = (coords[:,2] >= 21) & (coords[:,2] <= 40)
mask_guest = (coords[:,2] >= 41) & (coords[:,2] <= 60)


# .................................................................................
# Determine the interface location (z-axis) on either side of the cavity
print("Determining the interface location on either side of the cavity.")

# To determine the interface we need the positions of the host chains and cavity monomer
atom_positions_host = coords[:,4:][mask_host] # unwrapped coordinates
atom_positions_host_wrapped = utils.wrap_positions(atom_positions_host, simulation_box) # wrapped coordinates


# 1. Compute density profile of the host monomers.
bin_width = 10.0 # Angstrom
density_profile_host, bin_centers = utils.compute_number_density_profile(atom_positions_host_wrapped.reshape((1, -1, 3)),
                                                                         simulation_box,
                                                                         bin_width)
# Function expects positions at diff. timesteps, so we reshape the positions array.
# Reshaped to (1, num_atoms, 3) to match the expected input shape of the function.

# 2. Find the interface location by fitting a super gaussian to the density profile.
mask = (bin_centers >= args.cavity_bounds[0]) & (bin_centers <= args.cavity_bounds[1])
bin_centers_new = bin_centers[~mask]
density_profile_host_new = density_profile_host[~mask]
left_interface_coord, right_interface_coord, fine_coords, fitted_profile = utils.find_interfaces(bin_centers_new,
                                                                                                 density_profile_host_new)

print(f"Left interface coordinate: z = {left_interface_coord}")
print(f"Right interface coordinate: z = {right_interface_coord}\n")

# 4. Plot fitted profile and interface location as a sanity check.
# print("Plotting the fitted profile and interface location.\n")
# fig, ax = plt.subplots()
# ax.plot(bin_centers, density_profile_host, marker="o", linestyle="")
# ax.plot(fine_coords, fitted_profile, linestyle="-")
# ax.axvline(left_interface_coord, color="black", linestyle="--")
# ax.axvline(right_interface_coord, color="black", linestyle="--")
# ax.set_xlabel("z (Angstrom)")
# ax.set_ylabel("Number density (1/Angstrom^3)")
# fig.savefig(f"interface_location_t{timestep}.pdf", dpi=300)
# plt.close(fig)

# .................................................................................
"""
1. Iterate over every molecule in the "coords" array
2. If given molecule is a guest chain, then first determine its center of mass.
3. If COM is in the dilute regime, then we try to insert a chain in the cavity and if 
    we are able to insert it, then delete this chain.
"""
unique_mol_ids = np.unique(coords[:,1]) # Unique molecule ids in the coords array

print(f"Checking if guest chains are in the dilute phase and inserting them in the cavity if they are.")

# Create list to store mol ids of chains that were inserted (and deleted).
inserted_mol_ids = []
nguest_chains_dilute = 0
final_nguest_cavity = 0


atom_positions_cavity = coords[:,4:][mask_cavity] # unwrapped coordinates
atom_positions_cavity_wrapped = utils.wrap_positions(atom_positions_cavity, simulation_box) # wrapped coordinates


# Iterate through the mol_ids
for counter, umolid in enumerate(unique_mol_ids):
    # create mask for the coords
    mask = (coords[:,1] == umolid)
    
    # Get the atom types of the current molecule
    current_mol_atom_types = coords[mask, 2]
    
    # If the molecule is a guest chain
    if np.all((current_mol_atom_types >= 41) & (current_mol_atom_types <= 60)):
        # Get the atom positions of the current molecule
        current_mol_positions = coords[mask, 4:] # unwrapped coordinates
        # NOTE: The atom sequence can be jumbled, so we need to sort it out later.

        # Get the mass of the atoms
        current_mol_atom_masses = np.array(utils.extract_particle_masses(current_mol_atom_types))
        
        # Compute the COM of the current molecule (in unwrapped coordinates)
        # NOTE: Even if the atom sequence is jumbled, the function will compute the correct COM as the masses and positions
        # would correspond to the correct atoms.
        current_mol_com_pos = utils.compute_COM_position(current_mol_positions, current_mol_atom_masses) # unwrapped coordinates
        current_mol_com_pos = np.reshape(current_mol_com_pos, (1,3))

        # Wrap the COM position
        current_mol_com_pos_wrapped = utils.wrap_positions(current_mol_com_pos, simulation_box)[0] # wrapped coordinates
        
        # Check if the COM is in the dilute phase (with some tolerance)
        buffer = 10.0 # Angstrom
        if (current_mol_com_pos_wrapped[2] < (left_interface_coord - buffer) or current_mol_com_pos_wrapped[2] > (right_interface_coord + buffer)):
            # print(f"Chain with mol_id:{umolid} is in the dilute phase.")
            # Add to the count of guest chains in the dilute phase
            nguest_chains_dilute += 1

            ## Try to insert a new chain in the cavity
            # bounds for insertion
            box_bounds = np.array([[simulation_box[0,0],simulation_box[0,1]],
                                   [simulation_box[1,0],simulation_box[1,1]],
                                   [args.cavity_bounds[0], args.cavity_bounds[1]]])
            
            # Existing positions of guest and cavity chains
            atom_positions_guest = coords[:,4:][mask_guest] # unwrapped coordinates
            atom_positions_guest_wrapped = utils.wrap_positions(atom_positions_guest, simulation_box) # wrapped coordinates
            # NOTE: Cavity monomers included so that chains do not overlap with host monomers.
            existing_positions = np.vstack((atom_positions_guest_wrapped, atom_positions_cavity_wrapped))


            # Other parameters for chain insertion
            chain_length = current_mol_positions.shape[0] # We want to insert the same molecule but at a different location
            lattice_spacing = 3.81 # Angstrom (Based on the equilibrium bond distance in MPIPI)
            overlap_cutoff = 7.31 # Angstrom (For WF potential r_min = 1.19*sigma; Avg. sigma of the monomers = 6.14 Angstrom)
            max_attempts = 50000 # Max number of Monte Carlo attempts to insert the chain

            new_chain_position = utils.monte_carlo_insert_single_chain(box_bounds,
                                                                       existing_positions,
                                                                       chain_length,
                                                                       lattice_spacing,
                                                                       overlap_cutoff,
                                                                       simulation_box,
                                                                       max_attempts)
            
            # If the chain is successfully inserted, then delete the current chain
            if len(new_chain_position) > 0:
                # Replace the position of the current molecule with the new chain position
                
                # First, we need to determine the correct order of the monomers of the guest chain based on the bond connectivity.
                # (This is important because we need to insert the chain in the cavity in the same order.)
                
                # Get the atom ids for the current molecule
                current_mol_atom_ids = coords[mask, 0]
                
                # Get the bond data for the current molecule
                mask_bonds = np.isin(bonds[:,2], current_mol_atom_ids) | np.isin(bonds[:,3], current_mol_atom_ids)
                current_mol_bonds = bonds[mask_bonds]
                
                # Determine the correct sequence of atom ids for the current molecule
                ordered_atom_ids = utils.determine_correct_atom_id_sequence(current_mol_atom_ids, current_mol_bonds[:,2:])
                
                ## Replace the atom positions in coords array according to the correct atom id sequence
                # Iterate through the ordered atom ids
                for i, atom_id in enumerate(ordered_atom_ids):
                    # Find the index of the atom in the coords array
                    index = np.where(coords[:,0] == atom_id)[0][0]
                    # Replace the position of the atom with the new chain position
                    coords[index, 4:] = new_chain_position[i]
                
                # Store the mol id of the inserted chain
                inserted_mol_ids.append(umolid)

                # Add to the counter tracking number of guest chains inside cavity
                final_nguest_cavity += 1

            else:
                continue
        elif (current_mol_com_pos_wrapped[2] >= args.cavity_bounds[0] and current_mol_com_pos_wrapped[2] <= args.cavity_bounds[1]):
            # Add to the count of guest chains inside the cavity
            final_nguest_cavity += 1

print(f"Monte Carlo chain insertion completed.")
print(f"Number of guest chains in the dilute phase: {nguest_chains_dilute}")
print(f"Number of chains inserted: {len(inserted_mol_ids)}")
print(f"Number of unsuccessful chain insertion attempts: {nguest_chains_dilute - len(inserted_mol_ids)}\n")


# Now we need to check if the total number of guest chains inside the cavity is less than the target number
# If it is, then we need to add additional guest chains to the cavity
if final_nguest_cavity < args.target_num_guest_chains:
    print(f"Total number of guest chains ({final_nguest_cavity}) inside the cavity is less than the target number ({args.target_num_guest_chains}).")
    print(f"Adding additional guest chains to the cavity.\n")
    
    # Calculate the number of additional guest chains to add
    num_additional_chains = args.target_num_guest_chains - final_nguest_cavity
    print(f"Number of additional guest chains to add: {num_additional_chains}\n")

    # bounds for insertion
    box_bounds = np.array([[simulation_box[0,0],simulation_box[0,1]],
                           [simulation_box[1,0],simulation_box[1,1]],
                           [args.cavity_bounds[0], args.cavity_bounds[1]]])
    
    # Get the atom data for the last molecule in the coords array
    # NOTE: The logic below assumes that the last molecule in the coords array is a guest chain and 
    # the atoms are ordered in the trajecotry file.
    last_mol_id = np.max(coords[:,1])
    last_mol_atom_types = coords[coords[:,1] == last_mol_id, 2]
    last_mol_charges = coords[coords[:,1] == last_mol_id, 3]
    if np.all((last_mol_atom_types >= 41) & (last_mol_atom_types <= 60)): # Confirm if the last molecule is a guest chain
        print(f"Last molecule is a guest chain.")
    else:
        raise ValueError("Last molecule is not a guest chain (expected). Please check the input files.")


    # Add additional guest chains to the cavity
    for i in range(num_additional_chains):
        mask_guest = (coords[:,2] >= 41) & (coords[:,2] <= 60) # column 3 (index 2) contains atom type info.
        # Existing positions of guest and cavity chains
        atom_positions_guest = coords[:,4:][mask_guest] # unwrapped coordinates
        atom_positions_guest_wrapped = utils.wrap_positions(atom_positions_guest, simulation_box) # wrapped coordinates
        # NOTE: Cavity monomers included so that chains do not overlap with host monomers.
        existing_positions = np.vstack((atom_positions_guest_wrapped, atom_positions_cavity_wrapped))


        # Other parameters for chain insertion
        chain_length = current_mol_positions.shape[0] # We want to insert the same molecule but at a different location
        lattice_spacing = 3.81 # Angstrom (Based on the equilibrium bond distance in MPIPI)
        overlap_cutoff = 7.31 # Angstrom (For WF potential r_min = 1.19*sigma; Avg. sigma of the monomers = 6.14 Angstrom)
        max_attempts = 50000 # Max number of Monte Carlo attempts to insert the chain
        
        # Try to insert a new chain in the cavity
        new_chain_position = utils.monte_carlo_insert_single_chain(box_bounds,
                                                                   existing_positions,
                                                                   chain_length,
                                                                   lattice_spacing,
                                                                   overlap_cutoff,
                                                                   simulation_box,
                                                                   max_attempts)
        
        # If the chain is successfully inserted, then add it to the coords array
        if len(new_chain_position) > 0:
            # Create a new mol_id for the new chain
            new_mol_id = np.max(coords[:,1]) + 1
            
            # Create a new atom id for the new chain
            new_atom_ids = np.arange(np.max(coords[:,0]) + 1, np.max(coords[:,0]) + 1 + chain_length)
            
            # Create a new atom type for the new chain
            new_atom_types = last_mol_atom_types
            
            # Create a new charge for the new chain
            new_charges = np.zeros(chain_length) # It is okay to have zero here since the potential file resets the charges correctly
            
            # Create a new bond data for the new chain
            new_bond_data = np.zeros((chain_length - 1, 4), dtype=int)
            new_bond_data[:,0] = np.arange(np.max(bonds[:,0]) + 1, np.max(bonds[:,0]) + 1 + chain_length - 1)
            new_bond_data[:,1] = 1 # bond type
            new_bond_data[:,2] = new_atom_ids[:-1]
            new_bond_data[:,3] = new_atom_ids[1:]
            
            # Create a new coords array for the new chain
            new_coords = np.zeros((chain_length, 7))
            new_coords[:,0] = new_atom_ids
            new_coords[:,1] = new_mol_id
            new_coords[:,2] = new_atom_types
            new_coords[:,3] = new_charges
            new_coords[:,4:] = new_chain_position

            # Append the new coords to the existing coords array
            coords = np.vstack((coords, new_coords))

            # Append the new bond data to the existing bond data
            bonds = np.vstack((bonds, new_bond_data))


# .................................................................................
# Write the new configuration file.
print(f"Writing the updated configuration file: equilibrated_config_after_repopulation.dat \n")
utils.write_config("equilibrated_config_after_repopulation.dat",
                   simulation_box,
                   coords,
                   coords.shape[0],
                   bonds.shape[0],
                   bonds)


# .................................................................................
# Save the insertion/deletion information in a text file.
print(f"Saving the insertion/deletion information in the log file: {args.logfile} \n")
with open(args.logfile, "a+") as f:
    f.write(f"{timestep}\t{inserted_mol_ids}\n")

