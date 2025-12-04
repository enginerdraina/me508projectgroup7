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

-----------
Author Contributions:

Lucas: Completed methodology section in paper and presentation, wrote all of the MD/DFT code and most of the post-processing code. Ran all simulations as described in methodology. Updated GitHub repository.

Daniel: Ran simulations at two temperatures, wrote portions of the discussion and conclusion, edited report.

Riley: Ran simulations at one temperature, wrote portions of the discussion, conclusion and future work.  

Monika: Ran simulations at two temperatures, wrote portions of the introduction, results and conclusions.

Raina: Ran simulations at two temperatures, wrote temperature-dependent oscillation analysis code, wrote portion of results, wrote portion of discussion section, outlined portion of future work section.
