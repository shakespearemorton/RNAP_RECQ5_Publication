# Repository for Nature Review Submission

This repository contains all the data and analysis scripts associated with our submission of **Mechanisms of transcription attenuation and condensation of RNA polymerase II by RECQ5** to *Nature*. It is structured to facilitate easy access and review of the materials supporting our findings. Below is a guide to the repository's organization and contents.


### Figure 4 - Detailed Analysis and Supporting Data

- **SubPlot_A**
  - **Analysis_Code/**: Contains Python scripts for analyzing radial distribution functions (`analyticalrdf.py`, `rad_experimental.py`, `rad_simulated.py`).
  - **cryo-ET_data/**: Cryo-electron tomography data files in `.cmm` format, used in the analysis.

### Supplementary Figures

- **FigureS-17** and **FigureS-18**
  - Includes PDF and PPTX files for the figures.
  - **SubPlot_A_B/**: Analysis scripts (`plot_fe.py`, `plot_umbrellas.py`) and raw data for free energy calculations.
  - **SubPlot_C/**: Simulation files and density calculations for slabs.
  - **SubPlot_D/**: Additional analysis results in PDF format.

### Analysis Scripts and Raw Data

- **Radius_of_Gyration/**
  - **Analysis_Codes/**: Contains a Python script (`radius.py`) for calculating the radius of gyration from simulation data.
  - **raw_data/**: CSV files with analysis results for different protein configurations.

### Run Files for Simulations

- **Generate_Initial_Config/**: Scripts for generating initial configurations for molecular dynamics simulations.
- **Protein Files/**: Text files listing protein sequences used in simulations.
- **To_Run/**: Contains Python scripts for running simulations (`run_large_condensate.py`), compressing condensate data, calculating the radius of gyration, and performing umbrella sampling.


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
  - OpenMM 8.0.0
  - MDAnalysis 2.5.0
  - MDTraj 1.9.9.
  - Analysis of Umbrella sampling data is performed using “WHAM: an implementation of the weighted histogram analysis method”, http://membrane.urmc.rochester.edu/content/wham/  - version 2.0.11
