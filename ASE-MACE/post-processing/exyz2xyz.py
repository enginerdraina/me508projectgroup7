# code to convert extended (extxyz) to xyz
# vmd has a hard time reading extxyz
# makes comment lines readable instead of Lattice=blah.....

from ase.io import read, write

# read all frames from the extended xyz
atoms_list = read("equilibration-empty_ljc.extxyz", ":")

with open("equilibration-empty_ljc.xyz", "w") as f:
    for i, atoms in enumerate(atoms_list, start=1):
        n_atoms = len(atoms) # finds total number of atoms
        f.write(f"{n_atoms}\n") # write total atoms as first line
        f.write(f"Frame {i}\n") # writes frame number
        for sym, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
            f.write(f"{sym:2s} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n")
            # writes as:
            # total atoms
            # Frame {i}
            # N 1 2 3 (atom type, then x y z coord)
