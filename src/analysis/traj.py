import matplotlib.pyplot as plt
from src.analysis.colvar import plot_colvar_trajectories

STRIDE = 500
TIMESTEP = 0.002 * 10**-3  # in ns    

def plot_biases(colvar_df, ax, timestep, stride):
    # get all bias from columns and plot them seperately, as well as the sum
    column_biases = [col for col in colvar_df.columns if 'bias' in col]
    
    # Calculate moving average window size (100 points)
    window = 100
    
    time = colvar_df["time"] * timestep * stride
    for bias in column_biases:
        # Calculate moving average of bias
        rolling_bias = colvar_df[bias].rolling(window=window, center=True).mean()
        ax.plot(time, rolling_bias, label=bias, alpha=0.5)
        
    # Calculate moving average of total bias
    total_bias = colvar_df[column_biases].sum(axis=1)
    rolling_total = total_bias.rolling(window=window, center=True).mean()
    ax.plot(time, rolling_total, 'k--', label='total', alpha=0.3)

    ax.set_title("Bias (moving average, window = 100 points)")
    
    ax.legend()
    return ax

def plot_trajectory(colvar_df, directory, system, cvs):
    # HACK
    fig, axs = plt.subplots(3, 1, figsize=(12, 8))
    axs[:2] = plot_colvar_trajectories(colvar_df, cvs, axs[:2], timestep=TIMESTEP, stride=STRIDE)
    axs[2] = plot_biases(colvar_df, axs[2], timestep=TIMESTEP, stride=STRIDE)
    
    # Remove x-labels from upper plots
    axs[0].set_xlabel('')
    axs[1].set_xlabel('')
    axs[2].set_xlabel('Time [ns]')
    
    # Adjust spacing to accommodate labels
    plt.tight_layout()
    
    plt.savefig(
        f"{directory}/{system}_colvar_trajectories.png", 
        dpi=300, 
        bbox_inches='tight', 
        facecolor='white', 
        edgecolor='none'
    )
    plt.close()

from src.utils import get_checkpoint_interval
LOGGING_INTERVAL = 100 # ps per frame
CHECKPOINT_INTERVAL = 1000 # ps per checkpoint, 1 ns per checkpoint

from src.utils import get_restarted_files_by_extension, get_file_by_extension
import MDAnalysis as mda
import os
import numpy as np

def stitch_trajectories(directory, system, date):
    """
    Stitch multiple DCD trajectories into a list of AtomGroups.
    This final trajetory does not have the first equilibrated frame.
    TODO: One could extend this function to get that from equilibrated.cif.
    
    Parameters
    ----------
    directory : str
        Path to directory containing DCD files
    system : str
        System name
    date : str
        Date string
    
    Returns
    -------
    list
        List of AtomGroups, one for each frame
    """
    # Delete previous stitched trajectory
    try:
        full_dcd_file = get_file_by_extension(directory, '_full.dcd')
        os.remove(full_dcd_file)
    except FileNotFoundError:
        pass

    # Get topology and trajectory files
    dcd_files = get_restarted_files_by_extension(directory, '.dcd')
    pdb_file = get_file_by_extension(directory, '_fixed.pdb')
    
    if len(dcd_files) == 0:
        raise ValueError("No DCD files found in directory")
    
    assert dcd_files == sorted(dcd_files), "DCD files are not sorted"
    
    # Initialize universe from PDB for topology
    base_universe = mda.Universe(pdb_file)
    frames = []
    
    # Process each DCD file
    for i, dcd_file in enumerate(dcd_files):
        traj_universe = mda.Universe(pdb_file, dcd_file)
        
        # Verify topology matches
        assert base_universe.atoms.n_atoms == traj_universe.atoms.n_atoms, \
            f"Number of atoms mismatch in {dcd_file}"
        
        # Calculate frames to keep based on checkpoint interval
        frames_per_checkpoint = int(CHECKPOINT_INTERVAL / LOGGING_INTERVAL)
        
        # Determine how many frames to keep from this trajectory
        if i < len(dcd_files) - 1:
            last_frame_index = int(traj_universe.trajectory.n_frames - 
                                 (traj_universe.trajectory.n_frames % frames_per_checkpoint))
            print(f"Processing {dcd_file} with {last_frame_index} frames (trimmed)")
        else:
            last_frame_index = traj_universe.trajectory.n_frames
            print(f"Processing {dcd_file} with {last_frame_index} frames (untrimmed)")
        
        # Extract frames
        for ts in traj_universe.trajectory[:last_frame_index]:
            frames.append(ts.copy())
    
    print(f"Total frames collected: {len(frames)}")
    
    # Verify no consecutive frames are identical
    for i in range(len(frames)-1):
        assert not np.array_equal(frames[i].positions, frames[i+1].positions), \
            f"Consecutive frames {i} and {i+1} are identical"


    # Create combined universe with all frames
    combined_universe = mda.Universe(pdb_file)
    combined_universe.trajectory = mda.coordinates.memory.MemoryReader(np.array(frames))

    # Write the stitched trajectory
    output_path = f"{directory}/{system}_full.dcd"
    with mda.Writer(output_path, combined_universe.atoms.n_atoms) as w:
        for ts in combined_universe.trajectory:
            w.write(combined_universe)
    
    return combined_universe.trajectory

# def stitch_trajectories(directory, system, date):

#     try:
#         full_dcd_file = get_file_by_extension(directory, '_full.dcd')
#         os.remove(full_dcd_file)
#     except FileNotFoundError:
#         pass
    
#     dcd_files = get_restarted_files_by_extension(directory, '.dcd')
#     pdb_file = get_file_by_extension(directory, '_fixed.pdb')


    
    
#     atom_groups = []
#     for i, dcd_file in enumerate(dcd_files):
#         print(f"Processing {dcd_file}")
#         trajectory = mda.Universe(pdb_file, dcd_file)

#         import pdb; pdb.set_trace()
#         for atom_group in trajectory.atoms:
#             print(atom_group)
#             atom_groups.append(atom_group)


        
#         # HACK: doing checkpint trimming manually, but could do it via the actual .chk file
#         # but that requires separating the simulation intiailization into a function
#         last_frame_index = int(trajectory.n_frames - (trajectory.n_frames % (CHECKPOINT_INTERVAL / LOGGING_INTERVAL)))
#         # print(last_frame_index)
#         trajectory = trajectory[:last_frame_index]
#         _frames = [f for f in trajectory]

#         assert len(_frames) == last_frame_index, f"Expected {last_frame_index} frames, got {len(_frames)}"

#         for frame in trajectory:
#             # print(frames[-1].frame)
#             if len(frames) > 0:
#                 frame.time = frames[-1].time + 1
#             else:
#                 frame.time = 0
#             frames.append(frame)
    
#     pdb_file = get_file_by_extension(directory, '_fixed.pdb')
#     universe = mda.Universe(pdb_file)
#     print(universe.atoms)

#     # import pdb; pdb.set_trace()
#     # assert 0 == 1

#     with mda.Writer(f"{directory}/{system}_full.dcd", universe.atoms.n_atoms) as w:
#         for frame in frames:
#             w.write(frame)
        
        

# OBSOLETE
# def plot_colvar_traj_in_fes(directory, target, binder, num_runs):
#     # Create a wide figure to accommodate all runs
#     fig_width = 6 * num_runs  # 6 inches per subplot
#     fig, axes = plt.subplots(1, num_runs, figsize=(fig_width, 6))
#     if num_runs == 1:
#         axes = [axes]
    
#     # Process all runs in parallel
#     with ThreadPoolExecutor() as executor:
#         process_run = partial(process_single_run, 
#                             directory=directory, 
#                             target=target, 
#                             binder=binder, 
#                             shared_fig=fig, 
#                             shared_axes=axes)
        
#         run_data = list(executor.map(process_run, range(1, num_runs + 1)))
    
#     # Find the maximum trajectory length
#     max_frames = max(len(data['cv1_traj']) for data in run_data)
    
#     def init():
#         elements = []
#         for data in run_data:
#             data['line'].set_data([], [])
#             data['point'].set_data([], [])
#             elements.extend([data['line'], data['point']])
#         return elements

#     def animate(frame):
#         elements = []
#         for data in run_data:
#             # Handle different trajectory lengths
#             curr_frame = min(frame, len(data['cv1_traj'])-1)
            
#             # Update trajectory line
#             data['line'].set_data(data['cv1_traj'][:curr_frame], 
#                                 data['cv2_traj'][:curr_frame])
#             # Update current point
#             data['point'].set_data([data['cv1_traj'][curr_frame]], 
#                                  [data['cv2_traj'][curr_frame]])
#             elements.extend([data['line'], data['point']])
#         return elements

#     # Create animation with progress bar
#     frames = tqdm(range(max_frames), desc="Generating animation")
#     anim = FuncAnimation(
#         fig, 
#         animate, 
#         init_func=init,
#         frames=frames, 
#         interval=20,
#         blit=True
#     )
    
#     # Adjust layout and save
#     plt.tight_layout()
#     writer = PillowWriter(fps=30)
#     anim.save(f"{directory}/{target}_{binder}_all_trajectories.gif", writer=writer)
#     logger.info("Saved combined trajectory animation")
#     plt.close()