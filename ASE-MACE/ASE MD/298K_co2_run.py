#!/projectnb/me508/lclemens/.conda/envs/ASEMaceEnvironmentGPU/bin/python <-- !!! CHANGE TO YOURS !!!
# "conda env list" will give you this info
# last modified oct. 27, 2025

#import time # these you only need if dumping into csv format or code timing
#import csv  # these you only need if dumping into csv format or code timing 
import os
import numpy as np
from ase import units, Atoms
from ase.io import read, write
from ase.build import molecule, rotate # setup
from ase.constraints import FixAtoms # setup
from ase.md.verlet import VelocityVerlet # integrator
from ase.md.bussi import Bussi # thermostat
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary #initalize vel
from ase.md import MDLogger
from mace.calculators import mace_mp
from ase.optimize.sciopt import SciPyFminCG
from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator

# load in cif file and CO2 molecule
cif = read('900050_no_CO2.cif')
supercellNdims = (2,1,1) # trying to simulate "bulk" system
                         # i want co2 fully in a pore
system = cif.repeat(supercellNdims)
co2_2 = molecule('CO2') # load in a CO2
co2_2.center(vacuum=0)

# unit cell parameters (change these for your cell, re: VESTA)
a, b, c = 11.88310, 5.54180, 22.65150 # angstrom
beta = 101.6070  # degrees

# rotate CO2 to align with b-axis (so it points into pore)
co2_2.rotate(90, 'z', center='COP')

# position near upper-right pore
x_upperright = a * 0.75
z_upperright = c * 0.75
co2_2.translate([x_upperright, 0.0, z_upperright])

# move CO2 slightly into the pore (along b-axis)
cell = system.get_cell()
direction_b = cell[1] / np.linalg.norm(cell[1])
buffer_distance = -1.0  # angstroms
co2_2.translate(-buffer_distance * direction_b)

# add CO2 to the MOF supercell
system += co2_2
system.set_cell([2*a,b,c,90,beta,90])
system.pbc = (True, True, True)


#################
# !!! if multiple pores in your unit cell !!!
# further exploration: add second CO2 at bottom-left pore

## place it just outside pore entrance along -b direction (my pore dir)
#x_bottomleft = a * 0.25  # approximate bottom-left corner
#z_bottomleft = c * 0.25  # approximate bottom-left corner
#y_outside = -buffer_distance  # outside along -b axis
#co2_position = [x_bottomleft, y_outside, z_bottomleft]
#co2.translate(co2_position)

## align CO2 along b-axis (so it points into the pore)
## again you should align your CO2 wherever your pore sits
## click on a b c in VESTA to see that
## CO2 is linear, so rotate it to align with +y direction
#co2.rotate(90, 'z')  # rotate to align with y-axis

## combine mof and co2 into simulation
#system += co2
#################

# D3 cutoff is already ~halved from default value 
# move these around if calculation is too slow 
# they are set in TorchDFTD3Calculator
D3Cutoff   = 40*units.Bohr # was set to 50. too slow
D3CNthresh = 30*units.Bohr # was set to 40. too slow
# device shouldn't change 
device = 'cuda'

# keep precision to float32 UNLESS you are optimizing 
# optimization example below
precision64 = 'float64'
# load path to model (download manually) REPLACE W YOURS
model_path = os.path.expanduser('/projectnb/me508/students/lclemens/ASE-MACE/mace-mpa-0-medium.model')
MACEPOTENTIAL64 = mace_mp(model=model_path, device=device, default_dtype=precision64)
FULLPOTENTIAL64 = TorchDFTD3Calculator(dft=MACEPOTENTIAL64, device=device, every=1, delay=1, check=True, skin=2.0, cutoff=D3Cutoff, cnthr=D3CNthresh)

###################
# OPTIMIZATION
###################
# to keep the number of files generated to a minimum, lets not use log files here

# fix mof, only relax CO2
# let us assume mof is rigid :-)
n_mof_atoms = len(system) - 3  # 3 atoms in CO2
constraint = FixAtoms(indices=range(n_mof_atoms))
system.set_constraint(constraint)

# make sure you attach the correct 64 precision calculator for optimization
system.calc = FULLPOTENTIAL64
optcutoff=0.05 # default cutoff for forces to be under 0.05 ev per atom 
optimiser = SciPyFminCG(system)
optimiser.run(optcutoff)

# remove this constraint before simulation
system.set_constraint()

# reload calculator, but with precision32 to run dynamics
# 64 would work, but it would be 2x as slow

###################
# EQUILIBRATION (NVT)
###################
name_equil = 'equilibration' # checkpoint file 
logfile_equil = f"{name_equil}_ljc.log"
trajfile_equil = f"{name_equil}_ljc.traj"

# float32 for dynamics
precision32 = 'float32'
# load mace potential
MACEPOTENTIAL32 = mace_mp(model=model_path, device=device, default_dtype=precision32)
FULLPOTENTIAL32 = TorchDFTD3Calculator(dft=MACEPOTENTIAL32, device=device, every=1, delay=1, check=True, skin=2.0, cutoff=D3Cutoff, cnthr=D3CNthresh)

# ATTACH CALCULATOR TO SYSTEM AGAIN!
system.calc = FULLPOTENTIAL32

# initialize velocities at some temperature
temperature_K = 298 # kelvin
MaxwellBoltzmannDistribution(system, temperature_K=temperature_K)
Stationary(system)

# run NVT equilibration with thermostat (bussi)
timestep_fs = 1.0 * units.fs # 1 fs timestep
taut = 100 * units.fs
equil_steps = 50_000  # 50 ps equilibration (should be enough but can go to 100 ps)
# equilibrate with bussi thermostat, is velocity verlet + thermostat
dyn_equil = Bussi(system, timestep_fs, temperature_K=temperature_K, taut=taut, trajectory=trajfile_equil)

# set up logger
logger_equil = MDLogger(dyn_equil, system, logfile_equil, header=True, stress=False)

# attach logger to dynamics with appropriate interval
dyn_equil.attach(logger_equil, interval=100)  # log every 100 steps (0.1 ps)

# run equilibration
dyn_equil.run(equil_steps)

########################
# NVT PRODUCTION RUN
########################

## leave this commented out unless you are RERUNNING for some reason
## for example if you wanted to confirm dynamics with NVE then run NVT, or vice versa
# system_nve = read('equilibrated_run_ljc.traj') # note nvE ending, not nvT :-)
# system_nve.calc = FULLPOTENTIAL32

# keep track of files
name_nvt = 'production_nvt' # checkpoint
logfile_nvt = f"{name_nvt}_ljc.log"
trajfile_nvt = f"{name_nvt}_ljc.traj"
my_steps_nvt = 1_000_000 # 1 ns dynamics
# TODO: adjust length as needed

# set up nvt run
dyn_nvt = Bussi(system, timestep_fs, temperature_K=temperature_K, taut=taut, trajectory=trajfile_nvt)
logger_nvt = MDLogger(dyn_nvt, system, logfile_nvt, header=True, stress=False)
# attach logger to dynamics with appropriate interval !
# emphasis on appropriate interval
dyn_nvt.attach(logger_nvt, interval=1000) # log every ps (TODO: adjust as needed?)
# run nvt
dyn_nvt.run(my_steps_nvt)
