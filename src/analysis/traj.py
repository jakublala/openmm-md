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