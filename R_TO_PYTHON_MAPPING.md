# R to Python Conversion Guide

This document provides a mapping between the original R code and the Python conversion.

## Library/Package Mappings

| R Package | Python Package | Purpose |
|-----------|---------------|---------|
| `data.table` | `pandas` | Data manipulation |
| `ggplot2` | `matplotlib`, `seaborn`, `plotnine` | Visualization |
| `terra` | `rasterio` | Raster data handling |
| `sf` | `geopandas` | Vector spatial data |
| `keras` | `tensorflow.keras` | Neural network framework |
| `tensorflow` | `tensorflow` | Deep learning |
| `tfprobability` | `tensorflow_probability` | Probabilistic models |
| `compositions` | `scikit-bio` | Compositional data (CLR) |
| `ranger` | `scikit-learn` | Random forest (not used in final code) |
| `parallel` | `multiprocessing` | Parallel processing |

## Function Mappings

### Data Loading

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `rast("file.tif")` | `rasterio.open("file.tif")` | Opens raster file |
| `fread("file.csv")` | `pd.read_csv("file.csv")` | Reads CSV |
| `vect("file.shp")` | `gpd.read_file("file.shp")` | Reads shapefile |

### Data Manipulation

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `data[, col := value]` | `data['col'] = value` | Assignment |
| `data[condition]` | `data[condition]` | Filtering (similar syntax) |
| `data[, ..cols]` | `data[cols]` | Column selection |
| `na.omit(data)` | `data.dropna()` | Remove NA values |
| `rowSums(data)` | `data.sum(axis=1)` | Row-wise sum |
| `apply(data, 2, fun)` | `data.apply(fun, axis=0)` | Column-wise apply |

### Compositional Data

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `acomp(data)` | Handled in `clr()` | Composition object |
| `clr(acomp(data))` | `clr(data.values)` from `skbio` | Centered log-ratio |

### Spatial Operations

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `extract(raster, points)` | `raster.sample(coords)` | Sample raster at points |
| `sf::st_as_sf(data)` | `gpd.GeoDataFrame(data, geometry=...)` | Create spatial dataframe |
| `sf::st_intersects()` | `gpd.sjoin(..., predicate='within')` | Spatial join |
| `sf::st_crs()` | `gdf.crs` | Get CRS |

### Neural Network

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `layer_input(shape = ...)` | `layers.Input(shape=...)` | Input layer |
| `layer_conv_2d(...)` | `layers.Conv2D(...)` | 2D convolution |
| `activation_relu()` | `layers.Activation('relu')` or `activation='relu'` | ReLU activation |
| `layer_batch_normalization()` | `layers.BatchNormalization()` | Batch norm |
| `layer_spatial_dropout_2d()` | `layers.SpatialDropout2D()` | Spatial dropout |
| `layer_dropout()` | `layers.Dropout()` | Dropout |
| `layer_flatten()` | `layers.Flatten()` | Flatten |
| `layer_dense()` | `layers.Dense()` | Dense layer |
| `layer_concatenate()` | `layers.Concatenate()` | Concatenate |
| `keras_model(inputs, outputs)` | `keras.Model(inputs, outputs)` | Create model |

### TensorFlow Probability

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `tfp$distributions$MultivariateNormalDiag(...)` | `tfp.distributions.MultivariateNormalDiag(...)` | Distribution |
| `tfd_log_prob(y)` | `distribution.log_prob(y)` | Log probability |
| `tf$math$softplus(x)` | `tf.nn.softplus(x)` | Softplus activation |
| `tfp$layers$DistributionLambda(fn)` | `tfp.layers.DistributionLambda(fn)` | Distribution layer |

### Training

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `optimizer_nadam(...)` | `keras.optimizers.Nadam(...)` | Nadam optimizer |
| `model %>% compile(...)` | `model.compile(...)` | Compile model |
| `model %>% fit(...)` | `model.fit(...)` | Train model |
| `predict(model, data)` | `model.predict(data)` or `model(data)` | Make predictions |

### Visualization

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `ggplot(...) + geom_point(...)` | `plt.scatter(...)` or `ax.scatter(...)` | Scatter plot |
| `ggplot(...) + geom_raster(...)` | `plt.imshow(...)` or `ax.imshow(...)` | Raster plot |
| `ggsave(...)` | `plt.savefig(...)` | Save figure |
| `scale_fill_viridis_c()` | `cmap='viridis'` | Viridis colormap |
| `theme_bw()` | `plt.style.use('seaborn-whitegrid')` | Plot theme |
| `coord_equal()` | `ax.set_aspect('equal')` | Equal aspect ratio |

### Utilities

| R Code | Python Code | Notes |
|--------|-------------|-------|
| `Sys.time()` | `time.time()` | Current time |
| `set.seed(n)` | `np.random.seed(n)`, `tf.random.set_seed(n)` | Set random seed |
| `paste0(...)` | `f"..."` or `"".join(...)` | String concatenation |
| `seq(from, to, by)` | `np.arange(from, to, by)` | Sequence |
| `which(condition)` | `np.where(condition)` | Find indices |
| `sample(x, n)` | `np.random.choice(x, n)` | Random sample |

## Key Structural Differences

### 1. Pipe Operator

**R:**
```r
data %>% 
  filter(x > 0) %>%
  select(a, b, c)
```

**Python:**
```python
data[data['x'] > 0][['a', 'b', 'c']]
# or using method chaining
(data
 .query('x > 0')
 .loc[:, ['a', 'b', 'c']])
```

### 2. Assignment

**R:**
```r
data[, new_col := value]  # data.table
data$new_col <- value     # base R
```

**Python:**
```python
data['new_col'] = value
```

### 3. Indexing

**R (1-indexed):**
```r
data[1]      # First element
data[1:5]    # Elements 1-5 (inclusive)
```

**Python (0-indexed):**
```python
data[0]      # First element
data[0:5]    # Elements 0-4 (5 is exclusive)
```

### 4. Functions

**R:**
```r
my_function <- function(x, y = 10) {
  result <- x + y
  return(result)
}
```

**Python:**
```python
def my_function(x, y=10):
    result = x + y
    return result
```

### 5. Parallel Processing

**R:**
```r
results <- mclapply(data, function(x) process(x), mc.cores = 8)
```

**Python:**
```python
from multiprocessing import Pool

with Pool(processes=8) as pool:
    results = pool.map(process, data)
```

## Code Structure Comparison

### Original R Script Structure
1. Library loading
2. Global options and functions
3. Data loading
4. Data preprocessing
5. Visualization setup
6. Model building
7. Model training
8. Prediction and evaluation
9. Additional visualizations

### Python Module Structure
**geological_analysis.py:**
- Imports and setup
- Individual functions for each task
- Main function orchestrating the workflow

**visualization_utils.py:**
- Specialized plotting functions
- Polygon analysis utilities

This modular structure makes the Python code more maintainable and testable.

## Running the Code

### R
```r
source("Paper code after 1st revisions - includes ternary legends.R")
```

### Python
```bash
python geological_analysis.py
```

Or import as a module:
```python
from geological_analysis import load_elevation_data, build_bnn_model
```

## Testing

### R
No automated tests in original code.

### Python
```bash
python test_conversion.py  # Validate syntax and structure
```

## Notable Implementation Differences

1. **Color Space Conversion**: The R code uses `convertColor()` for Lab color space. The Python version includes a simplified approach. For full Lab color space conversion, consider using the `colormath` library.

2. **Ternary Plots**: The R code uses `ggtern` for ternary diagrams. The Python version includes a placeholder. For proper ternary plots in Python, use the `python-ternary` library.

3. **Spatial Operations**: Some spatial operations may be slightly different due to different underlying algorithms in `sf` vs `geopandas`.

4. **Memory Management**: Python and R handle memory differently. Large raster operations may require different approaches.

5. **Mixed Precision**: Both versions support mixed precision training, but configuration differs slightly between R's Keras wrapper and native TensorFlow/Keras in Python.
