# Geological Map Explanatory Power

Python implementation of the geological map analysis code from the paper "Quantifying the limited explanatory power of traditional geological maps".

## Original Code

This repository contains Python code converted from the original R implementation:
- Original R code: https://github.com/charliekirkwood/geological-map-explanatory-power
- Paper code file: `Paper code after 1st revisions - includes ternary legends.R`

## Overview

This code analyzes the relationship between terrain features and geochemistry using Bayesian Neural Networks (BNN). The analysis includes:

1. **Data Loading**: Loading elevation rasters and geochemistry data
2. **Data Preprocessing**: CLR transformation of compositional data
3. **Terrain Feature Extraction**: Extracting terrain patches around sample locations
4. **Model Training**: Training a Bayesian Neural Network to predict geochemistry from terrain
5. **Visualization**: Creating maps and plots showing predictions vs observations
6. **Evaluation**: Comparing BNN predictions with traditional geological map averages

## Requirements

### Python Dependencies

Install required packages:

```bash
pip install -r requirements.txt
```

Main dependencies:
- `numpy`, `pandas` - Data manipulation
- `rasterio`, `geopandas` - Geospatial data handling
- `scikit-bio` - Compositional data analysis (CLR transformation)
- `tensorflow`, `tensorflow-probability` - Bayesian neural network
- `matplotlib`, `seaborn`, `plotnine` - Visualization

### Data Files

The code expects the following data files in the `data/` directory:

1. **Elevation data**: `British Isles 250m Copernicus DEM.tif`
   - Digital elevation model for British Isles at 250m resolution

2. **Geochemistry data**: `streamsedimentgeochemistryUKandIE.csv`
   - Stream sediment geochemistry samples for UK and Ireland
   - Must include columns: `Easting_BNG`, `Northing_BNG`, `K2O`, `Fe2O3`, `CaO`

3. **Bedrock geology**: `digmap625_bedrock_arc/625k_V5_BEDROCK_Geology_Polygons.shp`
   - Bedrock geology polygons shapefile

## Usage

### Basic Usage

Run the main analysis script:

```bash
python geological_analysis.py
```

### Custom Analysis

The main script provides modular functions that can be imported and used separately:

```python
from geological_analysis import (
    load_elevation_data,
    load_geochemistry_data,
    prepare_composition_data,
    build_bnn_model
)

# Load data
elevdat, transform, crs = load_elevation_data('data/British Isles 250m Copernicus DEM.tif')
data = load_geochemistry_data('data/streamsedimentgeochemistryUKandIE.csv')

# Prepare compositional data
mapdat = prepare_composition_data(data, ['K2O', 'Fe2O3', 'CaO'])

# Build and train model
model = build_bnn_model(image_dim=27)
# ... train model with your data
```

## Project Structure

```
.
├── README.md                                           # This file
├── requirements.txt                                    # Python dependencies
├── geological_analysis.py                              # Main analysis script
├── Paper code after 1st revisions - includes ternary legends.R  # Original R code
├── data/                                               # Data directory (not included)
│   ├── British Isles 250m Copernicus DEM.tif
│   ├── streamsedimentgeochemistryUKandIE.csv
│   └── digmap625_bedrock_arc/
│       └── 625k_V5_BEDROCK_Geology_Polygons.shp
├── paperplots/                                         # Output plots (created on run)
└── models/                                             # Saved models (created on run)
```

## Key Functions

### Data Loading and Preprocessing

- `load_elevation_data(filepath)`: Load elevation raster
- `load_geochemistry_data(filepath)`: Load geochemistry CSV
- `extract_elevation_at_points(data, elevation_path)`: Sample elevation at points
- `prepare_composition_data(data, comp_elements)`: Apply CLR transformation
- `rgb_to_color(mapdat, comp_elements)`: Convert compositions to RGB colors

### Model Building and Training

- `build_bnn_model(image_dim, dropout_spatial, dropout_dense)`: Build BNN architecture
- `prepare_train_test_data(mapdat, imgs_norm, loc_norm)`: Create train/val/test splits
- `negative_log_likelihood(y_true, y_pred)`: Loss function for probabilistic model

### Visualization and Evaluation

- `plot_uk_elevation(elevdat, output_path)`: Plot elevation map
- `plot_training_history(history, output_path)`: Plot training curves
- `evaluate_model(model, x_test, y_test, comp_elements)`: Evaluate and plot results

## Model Architecture

The Bayesian Neural Network consists of:

1. **Convolutional Branch**: Processes terrain images
   - 3 Conv2D layers (256 filters each)
   - Batch normalization and dropout for regularization
   - Extracts spatial terrain features

2. **Auxiliary Branch**: Processes location coordinates
   - Dense layer for coordinate encoding
   - Captures spatial trends

3. **Combined Processing**: 
   - Concatenates terrain and location features
   - 2 dense layers (2048, 1024 units)
   - Outputs distribution parameters

4. **Probabilistic Output**:
   - Multivariate Normal distribution
   - Predicts mean and uncertainty for K2O, Fe2O3, CaO

## Conversion Notes

### Key Differences from R Code

1. **Compositional Data Analysis**: Uses `scikit-bio` for CLR transformation instead of R's `compositions` package

2. **Geospatial Processing**: Uses `rasterio` and `geopandas` instead of R's `terra` and `sf`

3. **Neural Network**: TensorFlow/Keras with TensorFlow Probability instead of R's Keras wrapper

4. **Visualization**: Uses matplotlib/seaborn instead of ggplot2

5. **Parallel Processing**: Uses Python's `multiprocessing` instead of R's `parallel::mclapply`

### Maintained Functionality

- Same model architecture (convolutional + auxiliary branches)
- Same dropout rates and normalization approach
- Same loss function (negative log-likelihood)
- Same evaluation metrics (R², RMSE)
- Same visualization approach (RGB composition maps)

## Performance Considerations

- **Terrain Extraction**: Most computationally intensive step (~3-5 minutes for full dataset)
- **Model Training**: Depends on GPU availability (hours on CPU, minutes on GPU)
- **Memory**: Large raster files may require significant RAM
- **GPU Acceleration**: Strongly recommended for model training

## Output

The script creates the following outputs:

### Plots Directory (`paperplots/`)

- `UKelevation.png`: UK elevation map
- `BNN_training_progress.png`: Training/validation loss curves
- `K_BNN_mean_holdout.png`, `Fe_BNN_mean_holdout.png`, `Ca_BNN_mean_holdout.png`: Element-specific predictions
- RGB composition maps (if full analysis is run)

### Models Directory (`models/`)

- Saved model weights (if training is completed)

## Citation

If you use this code, please cite the original paper:

```
[Paper citation to be added]
```

## License

[License to be determined based on original repository]

## Authors

- Original R code: Charlie Kirkwood and collaborators
- Python conversion: Converted from R implementation

## Contributing

This is a conversion of published research code. For the original implementation, see:
https://github.com/charliekirkwood/geological-map-explanatory-power
