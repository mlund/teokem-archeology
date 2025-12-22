#!/usr/bin/env python3
"""
Plot radial distribution functions (RDF) from bulk Monte Carlo simulation.

Reads rdf.csv and creates plots of g(r) for all species pairs.
Can be used as standalone script or interactively in Jupyter with ipywidgets.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path

# Try to import ipywidgets for interactive plotting
try:
    import ipywidgets as widgets
    from IPython.display import display
    IPYWIDGETS_AVAILABLE = True
except ImportError:
    IPYWIDGETS_AVAILABLE = False


def plot_rdf(csv_file='rdf.csv', output_file='rdf_plot.png'):
    """
    Plot radial distribution functions from CSV file.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing RDF data
    output_file : str
        Path to save the output plot
    """
    # Check if file exists
    if not Path(csv_file).exists():
        print(f"Error: File '{csv_file}' not found.")
        sys.exit(1)

    # Read CSV file
    print(f"Reading {csv_file}...")
    df = pd.read_csv(csv_file)

    # Extract distance column
    r = df['r'].values

    # Get all g(r) columns
    g_columns = [col for col in df.columns if col.startswith('g_')]

    if not g_columns:
        print("Error: No g(r) columns found in CSV file.")
        sys.exit(1)

    # Determine number of species from column names
    # Parse column names like 'g_1_1(r)' to extract species indices
    species_pairs = []
    for col in g_columns:
        # Extract numbers from column name like 'g_1_2(r)'
        parts = col.replace('g_', '').replace('(r)', '').split('_')
        if len(parts) == 2:
            species_pairs.append((int(parts[0]), int(parts[1])))

    num_species = max(max(pair) for pair in species_pairs)

    print(f"Found {num_species} species with {len(g_columns)} g(r) pairs")

    # Create figure
    if len(g_columns) == 1:
        # Single plot for one species pair
        fig, ax = plt.subplots(figsize=(8, 6))
        axes = [ax]
    elif num_species == 2:
        # 2x2 grid for 2 species (4 pairs)
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
    else:
        # Grid layout for multiple species
        nrows = num_species
        ncols = num_species
        fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.5*nrows))
        if num_species == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

    # Plot each g(r)
    for idx, col in enumerate(g_columns):
        if idx >= len(axes):
            break

        ax = axes[idx]
        g_r = df[col].values

        # Extract species indices from column name
        i, j = species_pairs[idx]

        # Plot g(r)
        ax.plot(r, g_r, linewidth=2, color='C0')
        ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='g(r) = 1')
        ax.set_xlabel('r (Å)', fontsize=11)
        ax.set_ylabel(f'$g_{{{i}{j}}}(r)$', fontsize=11)
        ax.set_title(f'Species {i}-{j}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(left=0)

        # Set y-axis to start at 0 if all values are positive
        if g_r.min() >= 0:
            ax.set_ylim(bottom=0)

    # Remove any unused subplots
    for idx in range(len(g_columns), len(axes)):
        fig.delaxes(axes[idx])

    # Overall title
    fig.suptitle('Radial Distribution Functions g(r)', fontsize=16, fontweight='bold', y=0.995)

    # Adjust layout
    plt.tight_layout()

    # Save figure
    print(f"Saving plot to {output_file}...")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved successfully!")

    # Also create a single plot with all g(r) on same axes
    if len(g_columns) > 1:
        fig2, ax2 = plt.subplots(figsize=(10, 7))

        # Color map for different pairs
        colors = plt.cm.tab10(np.linspace(0, 1, len(g_columns)))

        for idx, col in enumerate(g_columns):
            g_r = df[col].values
            i, j = species_pairs[idx]
            ax2.plot(r, g_r, linewidth=2, color=colors[idx],
                    label=f'$g_{{{i}{j}}}(r)$', alpha=0.8)

        ax2.axhline(y=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        ax2.set_xlabel('r (Å)', fontsize=12)
        ax2.set_ylabel('g(r)', fontsize=12)
        ax2.set_title('All Radial Distribution Functions', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        ax2.legend(fontsize=10, framealpha=0.9)
        ax2.set_xlim(left=0)

        if df[g_columns].values.min() >= 0:
            ax2.set_ylim(bottom=0)

        plt.tight_layout()

        combined_output = output_file.replace('.png', '_combined.png')
        print(f"Saving combined plot to {combined_output}...")
        plt.savefig(combined_output, dpi=300, bbox_inches='tight')
        print(f"Combined plot saved successfully!")

    # Show plots
    plt.show()


def load_rdf_data(csv_file='rdf.csv'):
    """
    Load RDF data from CSV file and extract species information.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing RDF data

    Returns
    -------
    df : pd.DataFrame
        RDF data
    species_pairs : list of tuples
        List of (i, j) species pairs
    g_columns : list of str
        Column names for g(r) data
    """
    if not Path(csv_file).exists():
        raise FileNotFoundError(f"File '{csv_file}' not found.")

    df = pd.read_csv(csv_file)

    # Get all g(r) columns
    g_columns = [col for col in df.columns if col.startswith('g_')]

    if not g_columns:
        raise ValueError("No g(r) columns found in CSV file.")

    # Parse species pairs from column names
    species_pairs = []
    for col in g_columns:
        parts = col.replace('g_', '').replace('(r)', '').split('_')
        if len(parts) == 2:
            species_pairs.append((int(parts[0]), int(parts[1])))

    return df, species_pairs, g_columns


def plot_interactive(csv_file='rdf.csv'):
    """
    Create interactive RDF plot with ipywidgets.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file containing RDF data

    Notes
    -----
    This function requires ipywidgets and is designed for use in Jupyter notebooks.
    """
    if not IPYWIDGETS_AVAILABLE:
        print("Error: ipywidgets is not available. Install with: pip install ipywidgets")
        print("Falling back to static plots...")
        plot_rdf(csv_file)
        return

    # Load data
    try:
        df, species_pairs, g_columns = load_rdf_data(csv_file)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        return

    r = df['r'].values

    # Create dropdown options
    pair_options = [('All pairs (combined)', 'all')]
    pair_options.extend([
        (f'Species {i}-{j}', f'g_{i}_{j}(r)')
        for i, j in species_pairs
    ])

    # Create widgets
    pair_dropdown = widgets.Dropdown(
        options=pair_options,
        value='all',
        description='Select pair:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='300px')
    )

    show_reference = widgets.Checkbox(
        value=True,
        description='Show g(r) = 1 reference line',
        style={'description_width': 'initial'}
    )

    linewidth_slider = widgets.FloatSlider(
        value=2.0,
        min=0.5,
        max=5.0,
        step=0.5,
        description='Line width:',
        style={'description_width': 'initial'},
        layout=widgets.Layout(width='400px')
    )

    output = widgets.Output()

    def update_plot(pair_selection, show_ref, lw):
        """Update plot based on widget values."""
        with output:
            output.clear_output(wait=True)

            fig, ax = plt.subplots(figsize=(10, 7))

            if pair_selection == 'all':
                # Plot all pairs
                colors = plt.cm.tab10(np.linspace(0, 1, len(g_columns)))

                for idx, col in enumerate(g_columns):
                    g_r = df[col].values
                    i, j = species_pairs[idx]
                    ax.plot(r, g_r, linewidth=lw, color=colors[idx],
                           label=f'$g_{{{i}{j}}}(r)$', alpha=0.8)

                ax.legend(fontsize=11, framealpha=0.9, loc='best')
                ax.set_title('All Radial Distribution Functions',
                           fontsize=14, fontweight='bold')
            else:
                # Plot specific pair
                g_r = df[pair_selection].values

                # Extract species indices
                idx = g_columns.index(pair_selection)
                i, j = species_pairs[idx]

                ax.plot(r, g_r, linewidth=lw, color='C0')
                ax.set_ylabel(f'$g_{{{i}{j}}}(r)$', fontsize=13)
                ax.set_title(f'Radial Distribution Function: Species {i}-{j}',
                           fontsize=14, fontweight='bold')

            # Reference line
            if show_ref:
                ax.axhline(y=1.0, color='gray', linestyle='--',
                          linewidth=1.5, alpha=0.6, label='g(r) = 1')

            ax.set_xlabel('r (Å)', fontsize=13)
            if pair_selection == 'all':
                ax.set_ylabel('g(r)', fontsize=13)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(left=0)

            # Set y-axis
            if pair_selection == 'all':
                y_min = df[g_columns].values.min()
            else:
                y_min = g_r.min()

            if y_min >= 0:
                ax.set_ylim(bottom=0)

            plt.tight_layout()
            plt.show()

    # Create interactive output
    interactive_plot = widgets.interactive(
        update_plot,
        pair_selection=pair_dropdown,
        show_ref=show_reference,
        lw=linewidth_slider
    )

    # Display widgets
    print(f"Loaded RDF data from {csv_file}")
    print(f"Found {len(species_pairs)} species pairs")
    print("\nUse the controls below to explore the data:\n")

    display(widgets.VBox([
        widgets.HBox([pair_dropdown]),
        widgets.HBox([show_reference, linewidth_slider]),
        output
    ]))

    # Trigger initial plot
    update_plot(pair_dropdown.value, show_reference.value, linewidth_slider.value)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Plot radial distribution functions from RDF CSV file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate static plots (default)
  ./plot_rdf.py rdf.csv

  # Use in Jupyter notebook for interactive plotting
  from plot_rdf import plot_interactive
  plot_interactive('rdf.csv')
        """)
    parser.add_argument('csv_file', nargs='?', default='rdf.csv',
                       help='Path to RDF CSV file (default: rdf.csv)')
    parser.add_argument('-o', '--output', default='rdf_plot.png',
                       help='Output plot filename (default: rdf_plot.png)')
    parser.add_argument('-i', '--interactive', action='store_true',
                       help='Launch interactive plot (requires ipywidgets)')

    args = parser.parse_args()

    if args.interactive:
        if not IPYWIDGETS_AVAILABLE:
            print("Error: ipywidgets is not available.")
            print("Install with: pip install ipywidgets")
            print("\nNote: Interactive mode is designed for Jupyter notebooks.")
            print("For command-line use, omit the -i flag for static plots.")
            sys.exit(1)
        plot_interactive(args.csv_file)
    else:
        plot_rdf(args.csv_file, args.output)
