import numpy as np
from ase.io import read, write

# load traj file
# repleace path with YOUR PATH. the command "pwd" tells you your path.
# index=":" argument reads all frames from the file
trajectory = read('/projectnb/me508/students/lclemens/ASE-MACE/OLDproduction_nvt_ljc.traj', index=":") 

# check how many frames were loaded
num_frames = len(trajectory)
# just in case your traj file is bad...
if num_frames == 0:
    print("Error: Trajectory file is empty or could not be read.")
    exit()
print(f"Loaded {num_frames} frames from the trajectory.")

# extract all atomic positions into a single numpy array
# shape of the array will be (num_frames, num_atoms, 3)
all_positions = np.array([atoms.get_positions() for atoms in trajectory])

# calculate average positions along the frames dimension (axis=0)
# result should be an array of shape (num_atoms, 3), representing the average x, y, z for each atom
average_positions = all_positions.mean(axis=0)

# create a new Atoms object with our average positions
# use the first frame as a template to get the cell, symbols, and pbc info
average_atoms = trajectory[0].copy()
average_atoms.set_positions(average_positions)

# write the average structure to a new file
write('298K_avg_structure_ljc.extxyz', average_atoms) # replace with YOUR initials and temp.

print("\nAverage coordinates calculated and saved to xyz file")
print("Average positions (first 5 atoms):")
print(average_positions[:5])
