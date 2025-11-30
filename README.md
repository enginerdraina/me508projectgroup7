# me508projectgroup7

Run MD in this order:
  1. 298K_co2_run.py
  2. equilibrate-empty-cell.py

  Using empty cell .traj, run
  
  3. 150K_co2_run.py

After this, run post-processing in this order:
  1. average-coordinates.py
  2. exyz2xyz.py
  3. avg-coord-error.ipynb
  4. oscillation_calcs.py
   
  To vizualize:
 
  5. trajectory_wrapping.py

Then, run Gaussian in any order after truncating and capping.
