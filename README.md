# Repository for Nature Review Submission

This repository contains all the data and analysis scripts associated with our submission of **Mechanisms of transcription attenuation and condensation of RNA polymerase II by RECQ5** to *bioRxiv*. It is structured to facilitate easy access and review of the materials supporting our findings. Below is a guide to the repository's organization and contents.


### Jupyter Notebook with Analysis Methods and Plotting
- **Figures_and_Analysis.ipynb**: Analysis for Figure 4, and Extended Data Figures 7/8.

### Raw Data and Figures
- **raw/**: Contains all raw data needed to recreate figures in the manuscript
- **figures/**: Contains copies of the original figures generated from modelling data

### Run Files for Simulations

- **Protein Files/**: Text files listing protein sequences used in simulations.
- **To_Run/**: Contains Python scripts for running simulations (`run_example.py`). Different simulation styles can be run or edited in sim_urry.py, allowing the user to mix and match for their use case. Turning on the 'special_rule' can be used to bind RECQ5 to RNAPII, and add the extra potential used to math the known binding affinity between SRI and CTD.

## Important Notes

The following is available upon request, but not included due to the restrictions on filesize.

- **Condensates/**
  - **Final_Positions/**: Contains XML files detailing the final positions of particles in various conditions.
  - **Post_Compression_Positions/**: XML files of particle positions post-compression under different conditions.
  - **Restarts/**: Checkpoint files enabling simulation restarts, ensuring reproducibility.
  - **Simulation_Information/**: CSV files with metadata and parameters for each simulation run, and a `readme.txt` providing additional context.
  - **Starting_Configurations/**: Initial configuration files for simulations, including PDB and XML formats, with a `readme.txt` for descriptions.

## Packages Used

  - Python 3.11.6
  - OpenMM 8.2.0
  - MDAnalysis 2.5.0
  - MDTraj 1.9.9.
  - Analysis of Umbrella sampling data is performed using “WHAM: an implementation of the weighted histogram analysis method”, http://membrane.urmc.rochester.edu/content/wham/  - version 2.0.11
