# guenon_simulation_stats

MATLAB and R code to run statistical analyses of results from evolutionary simulations of face pattern diversification in guenons

Written by Sandra Winters, sandra.winters@nyu.edu

Please cite:
Winters S & Higham JP. 2022. Simulated evolution of mating signal diversification in a primate radiation. Proceedings of the Royal Society B: Biological Sciences. https://doi.org/10.1098/rspb.2022.0734

Usage:
* Run simulations using code available here: https://github.com/sandrawinters/guenon_evolutionary_simulations

* Move all auto-generated folders of results to a folder called 'simulation_results'. The following scripts will work if this folder is in the current directory for MATLAB & R; alternatively, the location of this folder can be set at the beginning of each script.

* Run the following scripts to compile and format simulation data:  
extract_results.m  
compile_sim_data.R  
cluster_populations.m 

* To analyze face pattern diversification (Results section 1), run:  
analysis_clusters.R 

* To analyze face pattern distinctiveness between populations (Results section 2) and face pattern variability within populations (Results section 3), run:  
analysis_distances.R

* To analyze female mating biases (Results section 4), run:  
analysis_preferences.R

Data files from simulations presented in Winters & Higham 2022 Proc B (used in the analysis_*.R files above) can be downloaded here: https://doi.org/10.5061/dryad.gf1vhhmsf
