#!/usr/bin/env python3
"""
Plot radial distribution functions (RDF) from bulk Monte Carlo simulation.

Reads rdf.csv and creates plots of g(r) for all species pairs.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from pathlib import Path


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


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Plot radial distribution functions from RDF CSV file')
    parser.add_argument('csv_file', nargs='?', default='rdf.csv',
                       help='Path to RDF CSV file (default: rdf.csv)')
    parser.add_argument('-o', '--output', default='rdf_plot.png',
                       help='Output plot filename (default: rdf_plot.png)')

    args = parser.parse_args()

    plot_rdf(args.csv_file, args.output)
