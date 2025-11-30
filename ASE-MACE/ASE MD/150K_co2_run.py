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

# D3 cutoff is already ~halved from default value 
# move these around if calculation is too slow 
# they are set in TorchDFTD3Calculator
D3Cutoff   = 40*units.Bohr # was set to 50. too slow
D3CNthresh = 30*units.Bohr # was set to 40. too slow
# device shouldn't change 
device = 'cuda'

# load path to model (download manually) REPLACE W YOURS
model_path = os.path.expanduser('/projectnb/me508/students/lclemens/ASE-MACE/mace-mpa-0-medium.model')
precision32 = 'float32'
# load mace potential
MACEPOTENTIAL32 = mace_mp(model=model_path, device=device, default_dtype=precision32)
FULLPOTENTIAL32 = TorchDFTD3Calculator(dft=MACEPOTENTIAL32, device=device, every=1, delay=1, check=True, skin=2.0, cutoff=D3Cutoff, cnthr=D3CNthresh)

# ATTACH CALCULATOR TO SYSTEM AGAIN!

########################
# NVT PRODUCTION RUN
########################

# leave this commented out unless you are RERUNNING for some reason
# for example if you wanted to confirm dynamics with NVE then run NVT, or vice versa
system = read('equilibration-empty_ljc.traj') 
co2_2 = molecule('CO2') # load in a CO2
co2_2.center(vacuum=0)

# unit cell parameters (change these for your cell, re: VESTA)
a, b, c = 11.88310, 5.54180, 22.65150 # angstrom
beta = 101.6070  # degrees

# rotate CO2 to align with b-axis (so it points into pore)
co2_2.rotate(90, 'z', center='COP')

# position near upper-right pore
x_upperright = a * 0.5
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

system.calc = FULLPOTENTIAL32

# initialize velocities at some temperature
temperature_K = 150 # kelvin
MaxwellBoltzmannDistribution(system, temperature_K=temperature_K)
Stationary(system)

# run NVT equilibration with thermostat (bussi)
timestep_fs = 1.0 * units.fs # 1 fs timestep
taut = 100 * units.fs

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
