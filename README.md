# Script for patch prediction for "Patchy Nanoparticles by Atomic Stencilling (2025)"

1. All scripts can run on MATLAB v2022a or higher
2. All relavant convenience scripts are found in the "polyhedron_function" folder. Add this folder to the MATLAB path directory
3. Installation Requirements
4. A standard MATLAB installation is needed, please follow the MATLAB installation instructions

## Scripts

Name: corona_predictor.m  
Description: this script generates corona prediction. Parameters to be changed are indicated in the comments under the label "PARAMETERS TO CHANGE" seciton

Name: corona_plotter.m  
Description: plots the predicted corona from corona_prediction.m script 

## Sample Dataset

Name: cuboctahedron_mesh.dat  
Description: mesh binary file to define surface of cuboctahedron for use in corona_predictor.m

Name: dipyramid5_mesh.dat  
Description: mesh binary file to define surface of dipyramid for use in corona_predictor.m 

Name: rhombicdodecahedron_mesh.dat  
Description: mesh binary file to define surface of rhombicdodecahedron for use in corona_predictor.m

Name: octahedron_mesh.dat  
Description: mesh binary file to define surface of octahedron for use in corona_predictor.m

Name: octahedron_mesh.dat  
Description: mesh binary file to define surface of octahedron for use in corona_predictor.m

Name: rhombic_sample_corona.dat  
Description: example output of corona data (can be open in MATLAB via corona_plottter.m script)
