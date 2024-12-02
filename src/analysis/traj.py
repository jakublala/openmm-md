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

def plot_trajectory(colvar_df, directory, system):
    # HACK
    fig, axs = plt.subplots(3, 1, figsize=(12, 8))
    axs[:2] = plot_colvar_trajectories(colvar_df, axs[:2], timestep=TIMESTEP, stride=STRIDE)
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