from ase.io import read, write

traj = read('production_nvt_ljc.traj', ':')

wrapped_traj = []
for atoms in traj:
    atoms.wrap()
    wrapped_traj.append(atoms)

write('new150Kprod_nvt_wrapped_ljc.traj', wrapped_traj)
