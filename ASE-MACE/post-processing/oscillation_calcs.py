import numpy as np
from ase.io import read, write

# load traj file
trajectory = read('/projectnb/me508/students/lclemens/ASE-MACE/OLDproduction_nvt_ljc.traj', index=":")

num_frames = len(trajectory)
# just in case your traj file is bad
if num_frames == 0:
    print("Error: Trajectory file is empty or could not be read.")
    exit()
print(f"Loaded {num_frames} frames from the trajectory.")

reference = trajectory[0].get_positions()   # use first frame as reference
disp_per_frame = [] # to store displacements

max_frame = 1000001  # limit to first 20,000 frames; this is arbitrary and can change
print(f"Calculating displacements for the first {min(max_frame, num_frames)} frames")
for frame in trajectory[1:max_frame]:    # loop over each frame (only some subset of the beginning), start from 1 to not include the first
    disp = np.linalg.norm(frame.get_positions() - reference, axis=1)    # calculate displacement from reference for each atom; this is magnitude of displacement so already positive
    disp_per_frame.append(disp.mean()) # list of average displacements per frame (so each element is from one frame and it's the average of all atoms)

avg_disp_per_frame = np.array(disp_per_frame)
overall_avg_disp = avg_disp_per_frame.mean()  # overall average displacement across all frames
print(f"\nOverall average displacement from reference structure: {overall_avg_disp:.4f} A")
