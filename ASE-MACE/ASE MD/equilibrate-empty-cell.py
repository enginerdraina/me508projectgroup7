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
system.pbc = (True, True, True)

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
# make sure you attach the correct 64 precision calculator for optimization
system.calc = FULLPOTENTIAL64
optcutoff=0.05 # default cutoff for forces to be under 0.05 ev per atom 
optimiser = SciPyFminCG(system)
optimiser.run(optcutoff)

# reload calculator, but with precision32 to run dynamics
# 64 would work, but it would be 2x as slow

###################
# EQUILIBRATION (NVT)
###################
name_equil = 'equilibration-empty' # checkpoint file 
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
