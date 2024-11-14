import pandas as pd
import matplotlib.pyplot as plt
from src.constants import kB, mol

def plot_system_data(filepath, window=100, figsize=(10, 8)):
    """
    Plot and save time series data with moving averages for a given system.
    
    :param filepath: filepath
    :param window: Size of the moving average window
    :param figsize: Size of the figure (width, height)
    """
    def plot_with_moving_average(ax, df, x_column, y_column, moving_average=True):
        """
        Plot original data and its moving average on a given axis.
        """
        ax.plot(df[x_column], df[y_column], label='Original', alpha=0.5)
        if moving_average:
            rolling_mean = df[y_column].rolling(window=window, center=True, min_periods=1).mean()
            ax.plot(df[x_column], rolling_mean, label=f'Moving Average (window={window})', color='red')
        ax.set_xlabel(x_column)
        ax.set_ylabel(y_column)
        ax.legend()

    # Read the CSV file
    file_path = f"{filepath}.out"
    filename = file_path.split("/")[-1]
    try:
        df = pd.read_csv(file_path, sep=',')
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return

    # Get column names for indices 3, 4, 5, 6
    columns_to_plot = [
        "Potential Energy / $10^3 k_B T$",
        "Kinetic Energy /  $10^3 k_B T$",
        "Total Energy /  $10^3 k_B T$",
        "Temperature (K)",
        "Box Volume (nm^3)",
        "Density (g/mL)"
        ]

    # convert to kBT
    T = df['Temperature (K)']
    df["Potential Energy / $10^3 k_B T$"] = df["Potential Energy (kJ/mole)"] / (mol * kB * T)
    df["Kinetic Energy /  $10^3 k_B T$"] = df["Kinetic Energy (kJ/mole)"] / (mol * kB * T)
    df["Total Energy /  $10^3 k_B T$"] = df["Total Energy (kJ/mole)"] / (mol * kB * T)

    # Create a single figure with four subfigures
    fig, axs = plt.subplots(3, 2, figsize=figsize)
    fig.suptitle(f'System: {filename}', fontsize=16)

    # Flatten the 2x2 array of axes for easier iteration
    axs = axs.flatten()

    # Plot for each specified column
    for i, column in enumerate(columns_to_plot):
        plot_with_moving_average(
            axs[i], 
            df, 
            'Time (ps)', 
            column,
            moving_average=True
            )

    plt.tight_layout()
    plt.subplots_adjust(top=0.95)  # Adjust to prevent title overlap

    # TODO: might need to change naming for consitency
    output_pathfile = filepath.replace(".out", ".png")
    plt.savefig(output_pathfile, dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure to free up memory
