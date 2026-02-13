"""
Visualization utilities for geological map analysis.

This module contains functions for creating specialized plots including
ternary diagrams, RGB composition maps, and polygon-based visualizations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
import geopandas as gpd

# Default figure sizes and coordinate limits for British National Grid
DEFAULT_FIGSIZE = (7, 10)
BNG_XLIM = (-200000, 664000)
BNG_YLIM = (0, 1225000)


def plot_bedrock_polygons(bedrock_gdf, output_path="paperplots/bedrock_polygons.png"):
    """
    Plot bedrock geology polygons.
    
    Parameters:
    -----------
    bedrock_gdf : gpd.GeoDataFrame
        Bedrock polygons
    output_path : str
        Output file path
    """
    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    
    bedrock_gdf.plot(ax=ax, facecolor='lightgray', edgecolor='black', 
                     linewidth=0.2)
    
    ax.set_xlabel('Easting (metres BNG)')
    ax.set_ylabel('Northing (metres BNG)')
    ax.set_xlim(BNG_XLIM)
    ax.set_ylim(BNG_YLIM)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved bedrock polygon plot to {output_path}")


def plot_samples_with_bedrock(mapdat, bedrock_gdf, 
                               output_path="paperplots/samples_with_bedrock.png"):
    """
    Plot sample points colored by composition with bedrock overlay.
    
    Parameters:
    -----------
    mapdat : pd.DataFrame
        Sample data with 'col' column for colors
    bedrock_gdf : gpd.GeoDataFrame
        Bedrock polygons
    output_path : str
        Output file path
    """
    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    
    # Plot samples
    ax.scatter(mapdat['Easting_BNG'], mapdat['Northing_BNG'], 
               c=mapdat['col'], s=1, alpha=0.8, rasterized=True)
    
    # Overlay bedrock boundaries
    bedrock_gdf.boundary.plot(ax=ax, color='black', linewidth=0.2, alpha=0.3)
    
    ax.set_xlabel('Easting (metres BNG)')
    ax.set_ylabel('Northing (metres BNG)')
    ax.set_xlim(BNG_XLIM)
    ax.set_ylim(BNG_YLIM)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved samples with bedrock plot to {output_path}")


def plot_terrain_examples(imgs_norm, n_examples=3, seed=678,
                         output_dir="paperplots"):
    """
    Plot example terrain patches.
    
    Parameters:
    -----------
    imgs_norm : np.ndarray
        Normalized terrain images
    n_examples : int
        Number of examples to plot
    seed : int
        Random seed
    output_dir : str
        Output directory
    """
    np.random.seed(seed)
    indices = np.random.choice(len(imgs_norm), n_examples, replace=False)
    
    for i, idx in enumerate(indices):
        fig, ax = plt.subplots(figsize=(3, 3))
        
        im = ax.imshow(imgs_norm[idx], cmap='viridis')
        ax.set_xticks([])
        ax.set_yticks([])
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/terrain_example_{i+1}.png", 
                   dpi=300, bbox_inches='tight')
        plt.close()
    
    print(f"Saved {n_examples} terrain example plots to {output_dir}")


def plot_ternary_legend(rgb_data, element_names, 
                        output_path="paperplots/ternary_legend.png"):
    """
    Create a ternary plot as a legend for composition colors.
    
    Note: This is a simplified version. For full ternary plots,
    consider using the 'python-ternary' package.
    
    Parameters:
    -----------
    rgb_data : np.ndarray
        RGB normalized composition data (n_samples x 3)
    element_names : list
        Names of the three elements
    output_path : str
        Output file path
    """
    # For a proper ternary plot, would use python-ternary package
    # Here we create a simplified triangular color reference
    
    fig, ax = plt.subplots(figsize=(4, 4))
    
    # Sample points to create gradient
    n_points = 1000
    indices = np.random.choice(len(rgb_data), min(n_points, len(rgb_data)), 
                               replace=False)
    
    # Convert RGB to colors
    colors_hex = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) 
                  for r, g, b in rgb_data[indices]]
    
    # Create a simple scatter showing color distribution
    # (This is a placeholder - proper ternary plot would be more complex)
    x = rgb_data[indices, 0]  # Element 1
    y = rgb_data[indices, 1]  # Element 2
    
    ax.scatter(x, y, c=colors_hex, s=2, alpha=0.5)
    ax.set_xlabel(element_names[0])
    ax.set_ylabel(element_names[1])
    ax.set_title(f'Composition Color Reference\n({element_names[2]} shown by hue)')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved ternary legend to {output_path}")


def plot_prediction_map(preds_df, color_column='color',
                        output_path="paperplots/prediction_map.png"):
    """
    Plot spatial prediction map with RGB colors.
    
    Parameters:
    -----------
    preds_df : pd.DataFrame
        Predictions with x, y coordinates and color column
    color_column : str
        Column name containing color values
    output_path : str
        Output file path
    """
    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    
    ax.scatter(preds_df['x'], preds_df['y'], 
               c=preds_df[color_column], s=0.5, 
               marker='s', rasterized=True)
    
    ax.set_xlabel('Easting (metres BNG)')
    ax.set_ylabel('Northing (metres BNG)')
    ax.set_xlim(BNG_XLIM)
    ax.set_ylim(BNG_YLIM)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved prediction map to {output_path}")


def plot_polygon_composition_map(bedrock_gdf, color_column,
                                 output_path="paperplots/polygon_map.png"):
    """
    Plot bedrock polygons colored by average composition.
    
    Parameters:
    -----------
    bedrock_gdf : gpd.GeoDataFrame
        Bedrock polygons with composition colors
    color_column : str
        Column containing color hex codes
    output_path : str
        Output file path
    """
    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    
    bedrock_gdf.plot(ax=ax, color=bedrock_gdf[color_column], 
                     edgecolor='none')
    
    ax.set_xlabel('Easting (metres BNG)')
    ax.set_ylabel('Northing (metres BNG)')
    ax.set_xlim(BNG_XLIM)
    ax.set_ylim(BNG_YLIM)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved polygon composition map to {output_path}")


def create_comparison_plot(obs, pred, element_name, metrics=None,
                          output_path="paperplots/comparison.png"):
    """
    Create observed vs predicted scatter plot.
    
    Parameters:
    -----------
    obs : np.ndarray
        Observed values
    pred : np.ndarray
        Predicted values
    element_name : str
        Name of the element
    metrics : dict, optional
        Dictionary with 'r2' and 'rmse' keys
    output_path : str
        Output file path
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Scatter plot
    ax.scatter(obs, pred, alpha=0.3, s=10, rasterized=True)
    
    # 1:1 line
    lims = [min(obs.min(), pred.min()), max(obs.max(), pred.max())]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
    
    # Labels and title
    ax.set_xlabel('Observed (centered log-ratio)')
    ax.set_ylabel('Predicted (centered log-ratio)')
    
    if metrics:
        title = f"{element_name}: R² = {metrics['r2']:.3f}, RMSE = {metrics['rmse']:.3f}"
    else:
        title = element_name
    ax.set_title(title)
    
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_multi_element_comparison(obs_all, pred_all, element_names,
                                   output_path="paperplots/all_elements_comparison.png"):
    """
    Create a multi-panel comparison plot for all elements.
    
    Parameters:
    -----------
    obs_all : np.ndarray
        Observed values (n_samples x n_elements)
    pred_all : np.ndarray
        Predicted values (n_samples x n_elements)
    element_names : list
        Names of elements
    output_path : str
        Output file path
    """
    n_elements = len(element_names)
    fig, axes = plt.subplots(1, n_elements, figsize=(5*n_elements, 5))
    
    if n_elements == 1:
        axes = [axes]
    
    for i, (ax, elem) in enumerate(zip(axes, element_names)):
        obs = obs_all[:, i]
        pred = pred_all[:, i]
        
        # Calculate metrics
        r2 = np.corrcoef(obs, pred)[0, 1]**2
        rmse = np.sqrt(np.mean((obs - pred)**2))
        
        # Plot
        ax.scatter(obs, pred, alpha=0.3, s=10, rasterized=True)
        
        # 1:1 line
        lims = [min(obs.min(), pred.min()), max(obs.max(), pred.max())]
        ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
        
        ax.set_xlabel('Observed')
        ax.set_ylabel('Predicted')
        ax.set_title(f'{elem}\nR² = {r2:.3f}, RMSE = {rmse:.3f}')
        ax.grid(True, alpha=0.3)
        ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved multi-element comparison to {output_path}")


def assign_samples_to_polygons(samples_gdf, polygons_gdf):
    """
    Assign each sample to a bedrock polygon.
    
    Parameters:
    -----------
    samples_gdf : gpd.GeoDataFrame
        Sample points
    polygons_gdf : gpd.GeoDataFrame
        Bedrock polygons
        
    Returns:
    --------
    np.ndarray
        Array of polygon IDs for each sample
    """
    print("Assigning samples to polygons...")
    
    # Spatial join
    joined = gpd.sjoin(samples_gdf, polygons_gdf, how='left', predicate='within')
    
    polygon_ids = joined['polygon_ID'].values
    
    print(f"Assigned {(~pd.isna(polygon_ids)).sum()} samples to polygons")
    
    return polygon_ids


def calculate_polygon_averages(samples_df, polygon_ids, value_columns):
    """
    Calculate average values within each polygon.
    
    Parameters:
    -----------
    samples_df : pd.DataFrame
        Sample data
    polygon_ids : np.ndarray
        Polygon ID for each sample
    value_columns : list
        Columns to average
        
    Returns:
    --------
    pd.DataFrame
        Averages by polygon
    """
    df = samples_df.copy()
    df['polygon_id'] = polygon_ids
    
    averages = df.groupby('polygon_id')[value_columns].mean()
    
    return averages
