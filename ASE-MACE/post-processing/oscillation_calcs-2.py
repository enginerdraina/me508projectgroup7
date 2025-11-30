import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write

# load traj file
trajectory = read('/projectnb/me508/students/krtung/project/ASE/production_nvt_150_2_RT.traj', index=":")

num_frames = len(trajectory)
# just in case your traj file is bad
if num_frames == 0:
    print("Error: Trajectory file is empty or could not be read.")
    exit()
print(f"Loaded {num_frames} frames from the trajectory.")

prev_pos = trajectory[0].get_positions()   # use first frame for first previous position
disp_per_frame = [] # to store displacements

for frame in trajectory[1:]:    # loop over each frame (only some subset of the beginning), start from 1 to not include the first
    current_pos = frame.get_positions()
    disp = np.linalg.norm(current_pos - prev_pos, axis=1)    # calculate displacement from reference for each atom; this is magnitude of displacement so already positive
    disp_per_frame.append(disp.mean()) # list of average displacements per frame (so each element is from one frame and it's the average of all atoms)
    prev_pos = current_pos

avg_disp_per_frame = np.array(disp_per_frame)
overall_avg_disp = avg_disp_per_frame.mean()  # overall average displacement across all frames
print(f"\nOverall average displacement from reference structure: {overall_avg_disp:.4f} A")

# Plot average displacement vs frame index
plt.figure(figsize=(8,5))
plt.plot(avg_disp_per_frame, 'o', markersize=.5, label='Average displacement')
plt.xlabel('Frame')
plt.ylabel('Average Euclidean Distance (Å)')
plt.ylim(ymin=0) 
plt.title('150K (RT)')   # Change for each temp
plt.grid(True)

# Save figure
plt.savefig("avg_distances_vs_frame.png", dpi=300, bbox_inches='tight')
plt.close()

print("Saved plot: avg_distances_vs_frame.png")