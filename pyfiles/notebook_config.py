# ===================================================
# Notebook configuration settings
# By Shahab Fatemi
# ===================================================
#
import numpy as np
import matplotlib
import matplotlib.cm as cm
from cycler import cycler
import matplotlib.pyplot as plt
from IPython.display import clear_output

# High-resolution display for figures
matplotlib.rcParams['figure.dpi']       = 200
matplotlib.rcParams['figure.figsize']   = (6, 4)

# Automatically use tight layout for all figures
matplotlib.rcParams['figure.constrained_layout.use'] = True

# General font sizes for plots
matplotlib.rcParams['font.family']      = 'serif'
matplotlib.rcParams['axes.titlesize']   = 16  # Title font size
matplotlib.rcParams['axes.labelsize']   = 14  # Axis label font size
matplotlib.rcParams['axes.linewidth']   = 1.0
matplotlib.rcParams['xtick.labelsize']  = 12  # X-axis tick font size
matplotlib.rcParams['ytick.labelsize']  = 12  # Y-axis tick font size

# Grid settings
matplotlib.rcParams['axes.grid']      = True   # Always show grid
matplotlib.rcParams['grid.color']     = 'grey' # Color
matplotlib.rcParams['grid.linestyle'] = '--'   # Dashed lines
matplotlib.rcParams['grid.alpha']     = 0.5    # Transparency
matplotlib.rcParams['grid.linewidth'] = 0.5    # Optional for finer grids

# Line markers
matplotlib.rcParams['lines.linewidth']  = 2.0
matplotlib.rcParams['lines.markersize'] = 8

# Legend settings
matplotlib.rcParams['legend.fontsize']  = 12      # Legend font size
matplotlib.rcParams['legend.frameon']   = True    # Show legend frame
matplotlib.rcParams['legend.loc']       = 'best'  # Best location for legend

# Pick a categorical colormap
colors = matplotlib.colormaps['tab10'].colors   # gives you a list of 10 RGB tuples

# Set the color cycle in rcParams
matplotlib.rcParams['axes.prop_cycle'] = cycler(color=colors)

# Save figure settings
matplotlib.rcParams['savefig.dpi']      = 300     # DPI for saved figures
matplotlib.rcParams['savefig.format']   = 'png'   # Default format to save figures

# ===================================================
# END
# ===================================================