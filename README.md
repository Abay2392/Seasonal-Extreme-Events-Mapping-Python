# Seasonal-Extreme-Events-Mapping-Python
Python workflow for seasonal mapping of 95th and 99th percentile rainfall and temperature extremes using spatial interpolation, Gaussian smoothing, and geospatial visualization.
# Seasonal Mapping of Extreme Rainfall and Temperature Events

## Overview

This repository contains a Python workflow for producing seasonal spatial maps of extreme rainfall and temperature events based on 95th and 99th percentile exceedances.

The workflow processes event datasets, aggregates seasonal extreme statistics, interpolates station observations onto regular grids, and generates publication-quality maps.

The methodology was developed for Northeast England and can be adapted to any region with station-based extreme event datasets.

---

## Study Variables

The workflow maps:

### Rainfall Extremes

* 95th percentile rainfall exceedances
* 99th percentile rainfall exceedances

### Maximum Temperature Extremes

* 95th percentile Tmax exceedances
* 99th percentile Tmax exceedances

### Minimum Temperature Extremes

* 95th percentile Tmin cold extremes
* 99th percentile Tmin cold extremes

---

## Seasonal Classification

Extreme events are grouped into four climatological seasons:

| Season | Months                     |
| ------ | -------------------------- |
| DJF    | December–January–February  |
| MAM    | March–April–May            |
| JJA    | June–July–August           |
| SON    | September–October–November |

---

## Input Datasets

The workflow uses CSV files containing:

* Event date
* Latitude
* Longitude
* Observed value
* Threshold exceedance magnitude

Input files include:

```text
Rain_95pct_Events.csv
Exceedance_99pct_Events_corrected.csv
Tmax_95th_Exceedance_Events.csv
Tmax_99th_Exceedance_Events.csv
Tmin_ColdExtreme_95th.csv
Tmin_ColdExtreme_99th.csv
```

---

## Spatial Processing Workflow

### Step 1: Event Selection

Events are classified by season using event dates.

### Step 2: Seasonal Aggregation

For each station location:

#### Rainfall and Tmax

Maximum seasonal exceedance values are extracted.

#### Tmin

Minimum seasonal temperatures and maximum threshold departures are extracted.

---

### Step 3: Spatial Interpolation

Station observations are interpolated to a regular geographic grid using:

```python
scipy.interpolate.griddata()
```

Interpolation method:

```text
Linear interpolation
```

---

### Step 4: Spatial Smoothing

Interpolated surfaces are smoothed using:

```python
scipy.ndimage.gaussian_filter()
```

to improve map readability and reduce interpolation artefacts.

---

### Step 5: Geographic Masking

Interpolated surfaces are clipped to the study region using:

* GeoPandas
* Shapely geometry operations

Only values within the regional boundary are retained.

---

## Outputs

For each variable and season, two maps are generated:

### Main Extreme Variable

Examples:

* Rainfall amount
* Maximum temperature
* Minimum temperature

### Exceedance Magnitude

Examples:

* Rainfall above threshold
* Temperature above threshold
* Temperature below threshold

---

## Generated Maps

The workflow automatically creates:

### Rainfall

* DJF Rainfall 95th percentile

* MAM Rainfall 95th percentile

* JJA Rainfall 95th percentile

* SON Rainfall 95th percentile

* DJF Rainfall 99th percentile

* MAM Rainfall 99th percentile

* JJA Rainfall 99th percentile

* SON Rainfall 99th percentile

### Tmax

* Seasonal 95th percentile exceedance maps
* Seasonal 99th percentile exceedance maps

### Tmin

* Seasonal 95th percentile cold extreme maps
* Seasonal 99th percentile cold extreme maps

---

## Mapping Features

The maps include:

* Administrative boundaries
* High-resolution interpolation
* Publication-quality colour scales
* Geographic coordinate labels
* Automatic colour bar scaling
* Seasonal titles and metadata

---

## Python Libraries

Required packages:

```text
numpy
pandas
geopandas
matplotlib
scipy
shapely
```

---

## Output Resolution

Regular grid spacing:

```text
0.009°
```

Spatial smoothing:

```text
Gaussian filter (σ = 1.0)
```

---

## Applications

This workflow can be used for:

* Compound climate extremes
* Extreme rainfall analysis
* Heatwave assessment
* Cold spell assessment
* Climate hazard mapping
* Infrastructure risk studies
* Power outage impact assessments
* Climate change adaptation planning

---

## Author

Dr. Tagele Mossie Aschale

---

## Citation

If this code contributes to your research, please cite the associated publication and acknowledge this repository.

