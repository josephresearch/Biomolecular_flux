import sys
sys.path.insert(1, '/home/yw9071/scripts')
sys.path.insert(2, '/home/yw9071/scripts/RNA_flux_project')
from RNA_flux_project import utils
import numpy as np
import argparse

# Input data
parser = argparse.ArgumentParser()
# Necessary arguments
parser.add_argument("--option", type=str, required=True, choices=["cavity+protein1", "cavity+protein1+protein2", "cavity+protein1+protein2+peptides", "cavity+peptides", "cavity+protein1+protein2+peptides_with_angle"])
parser.add_argument("--simBox", nargs=6, type=float, required=True, help='Slab size: xlo xhi ylo yhi zlo zhi')
# Optional arguments
parser.add_argument("--cavity_spacing", type=float, default=5.5, required=False, help='Spacing between cavity monomers in Angstroms (default: 5.5)')
parser.add_argument("--prot1_seq", type=str, required=False)
parser.add_argument("--prot1_config", type=str, required=False)
parser.add_argument("--prot1_traj", type=str, required=False)

parser.add_argument("--cavity_and_prot1_config", type=str, required=False)
parser.add_argument("--cavity_and_prot1_traj", type=str, required=False)
parser.add_argument("--prot2_config", type=str, required=False)
parser.add_argument("--prot2_traj", type=str, required=False)

parser.add_argument("--cavity_prot1_and_prot2_config", type=str, required=False)
parser.add_argument("--cavity_prot1_and_prot2_traj", type=str, required=False)
parser.add_argument("--cavity_peptides_seq", type=str, required=False)
parser.add_argument("--cavity_peptides_nchains", type=str, required=False) # reading as string since input "0" creates problem reading as int with argparse
parser.add_argument("--kangle", type=float, required=False)
args = parser.parse_args()

# Check what option was chosen
if args.option == "cavity+protein1":
    if args.prot1_config and args.prot1_traj and args.prot1_seq and args.cavity_spacing:
        # Cavity related variables
        x_length = args.simBox[1] - args.simBox[0] # same as simulation box
        y_length = args.simBox[3] - args.simBox[2] # same as simulation box
        z_length = 100 # 10 nm width of the cavity

        # Get cavity positions and types
        cavity_positions, cavity_types = utils.generate_cuboid_cavity_with_exact_AA_composition_2(args.prot1_seq,
                                                                                                x_length,
                                                                                                y_length,
                                                                                                z_length,
                                                                                                spacing=args.cavity_spacing)
        print(f"# cavity monomers = {cavity_positions.shape[0]}")
        # Get protein coordinates from trajectory file
        prot1_box_bounds, prot1_coords, prot1_highest_mol_id, _ = utils.read_trajectory(args.prot1_traj)
        # Get protein bond date from config file
        prot1_num_atoms, prot1_num_bonds, prot1_bonds = utils.read_config(args.prot1_config)
        # Generate config file with a cuboidal cavity surrounded by proteins on both sides
        utils.write_config_cuboid_cavity_with_inner_prot([[args.simBox[0], args.simBox[1]], [args.simBox[2], args.simBox[3]], [args.simBox[4], args.simBox[5]]],
                                                         cavity_positions, cavity_types,
                                                         prot1_coords, prot1_bonds)
    else:
        raise ValueError("You didn't enter config or trajectory files or the sequence for inner protein!!")


elif args.option == "cavity+protein1+protein2":
    if args.cavity_and_prot1_config and args.cavity_and_prot1_traj and args.prot2_config and args.prot2_traj:
        # Get cavity and inner protein coordinates from trajectory file
        box_bounds1, cavity_with_inner_protein_coords, cavity_with_inner_protein_highest_molid, _ = utils.read_trajectory(args.cavity_and_prot1_traj)
        # Get inner protein bond data from config file
        cavity_and_inner_protein_num_atoms, inner_protein_num_bonds, inner_protein_bonds = utils.read_config(args.cavity_and_prot1_config)
        # Get outer protein coordinates from trajectory file
        box_bounds2, outer_protein_coords, outer_protein_highest_molid, _ = utils.read_trajectory(args.prot2_traj)
        # Get outer protein bond data from config file
        outer_protein_num_atoms, outer_protein_num_bonds, outer_protein_bonds = utils.read_config(args.prot2_config)
        # Generate config file with outer protein added on either sides of the cavity and inner proteins
        utils.write_config_cuboid_cavity_with_inner_and_outer_prot([[args.simBox[0], args.simBox[1]], [args.simBox[2], args.simBox[3]], [args.simBox[4], args.simBox[5]]],
                                                    cavity_with_inner_protein_coords, inner_protein_bonds,
                                                    outer_protein_coords, outer_protein_bonds)
    else:
        raise ValueError("Insufficient input given!!")

elif args.option == "cavity+protein1+protein2+peptides": # Here we can also have the case where there is no protein2 and it will still work
    if args.cavity_prot1_and_prot2_config and args.cavity_prot1_and_prot2_traj and args.cavity_peptides_seq and args.cavity_peptides_nchains:
        # Convert input for # peptides from str to int
        nchains = int(args.cavity_peptides_nchains)
        # Get system coordinates from trajectory file
        box_bounds1, system_coords, system_highest_molid, _ = utils.read_trajectory(args.cavity_prot1_and_prot2_traj)
        # Get system bond data from config file
        system_num_atoms, system_num_bonds, system_bonds = utils.read_config(args.cavity_prot1_and_prot2_config)
        # Generate config file with peptides inside cavity with proteins on its either sides
        previous_cavity_centre_wrapped_in_new_simBox = args.simBox[4] + (args.simBox[5] - args.simBox[4])/2
        cavity_dimensions = [[0., 200.], 
                             [0., 200.], 
                             [-50.,50.]]
        utils.add_unequilibrated_chains_to_cuboidal_cavity([[args.simBox[0], args.simBox[1]], [args.simBox[2], args.simBox[3]], [args.simBox[4], args.simBox[5]]],
                                                            system_coords, system_bonds,
                                                            previous_cavity_centre_wrapped_in_new_simBox, cavity_dimensions,
                                                            args.cavity_peptides_seq, nchains)
    else:
        raise ValueError("Insufficient input given!!")

elif args.option == "cavity+protein1+protein2+peptides_with_angle": # Here we can also have the case where there is no protein2 and it will still work
    if args.cavity_prot1_and_prot2_config and args.cavity_prot1_and_prot2_traj and args.cavity_peptides_seq and args.cavity_peptides_nchains and args.kangle:
        # Convert input for # peptides from str to int
        nchains = int(args.cavity_peptides_nchains)
        # Get system coordinates from trajectory file
        box_bounds1, system_coords, system_highest_molid, _ = utils.read_trajectory(args.cavity_prot1_and_prot2_traj)
        # Get system bond data from config file
        system_num_atoms, system_num_bonds, system_bonds = utils.read_config(args.cavity_prot1_and_prot2_config)
        # Generate config file with peptides inside cavity with proteins on its either sides
        previous_cavity_centre_wrapped_in_new_simBox = args.simBox[4] + (args.simBox[5] - args.simBox[4])/2
        cavity_dimensions = [[0., 200.], 
                             [0., 200.], 
                             [-50.,50.]]
        utils.add_unequilibrated_chains_to_cuboidal_cavity_with_angle_potential([[args.simBox[0], args.simBox[1]], [args.simBox[2], args.simBox[3]], [args.simBox[4], args.simBox[5]]],
                                                            system_coords, system_bonds,
                                                            previous_cavity_centre_wrapped_in_new_simBox, cavity_dimensions,
                                                            args.cavity_peptides_seq, nchains,
                                                            kangle = args.kangle)
    else:
        raise ValueError("Insufficient input given!!")

elif args.option == "cavity+peptides":
    if args.prot1_seq and args.cavity_peptides_seq and args.cavity_peptides_nchains:
        # Convert input for # peptides from str to int
        nchains = int(args.cavity_peptides_nchains)
        
        simBox = [[args.simBox[0], args.simBox[1]], [args.simBox[2], args.simBox[3]], [args.simBox[4], args.simBox[5]]]
        utils.write_config_cuboidal_cavity_w_peptides_wo_outer_proteins(simBox, args.prot1_seq, args.cavity_peptides_seq, nchains)
    else:
        raise ValueError("Insufficient input given!!")