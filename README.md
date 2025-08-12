# PyTransHelio
A Python-based Tool for Calculating and Visualizing Weighted Centroids of Full-Disk Solar Active Regions  

# Environment Configuration Guide

This system requires Python 3.11 or higher. The detailed configuration steps are as follows:

## Python Environment Preparation

### Version Check:`python --version`

If the version is lower than 3.11 or not installed, an upgrade/installation is required.

### Installation Methods:
- **macOS**: Use `brew install python` or download from the official website  
- **Windows**: Download the installer from the Python official website and check "Add Python to PATH" during installation  

## Dependency Installation

It is recommended to create a virtual environment to isolate dependencies:
`python -m venv venv`

## Activate the environment

```
source venv/bin/activate # macOS/Linux

venv\Scripts\activate # Windows
```
Install dependencies using pip:
`pip install -r requirements.txt`

## Verification

After installation, verify with the following command: 
`pip list`
After that,it should display output similar to the following:
```
Package           Version
----------------- ---------
asdf              4.1.0
asdf-astropy      0.7.1
astropy           7.0.1
astropy-healpix   1.1.2
dask              2025.2.0
matplotlib        3.10.0
numpy             2.2.3
opencv-python     4.11.0
pandas            2.2.3
reproject         0.14.1
scikit-image      0.25.2
sunpy             6.1.1
```

