# Supplementary Data Codes 

**Direct observation of 3D atomic packing in monatomic amorphous materials**

Coherent Imaging Group, UCLA

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

Liquids and solids are two fundamental states of matter. Although the structure of crystalline solids has long been solved by crystallography, our understanding of the 3D atomic structure of liquids and amorphous materials remained speculative due to the lack of direct experimental determination. Here we advance atomic electron tomography to resolve for the first time the 3D atomic positions in monatomic amorphous materials, including a Ta thin film and two Pd nanoparticles. The experimental data and source codes for the 3D image reconstruction and post analysis are provided here.

# System Requirements

## Hardware Requirements

We recommend a computer with 16G DRAM, standard i7 4-core CPU, and a GPU to run most data analysis source codes. But for the 3D reconstruction of the experimental data with RESIRE, atomic tracing and refinement, we recommend a computer with large memory (256G DRAM, 16-core CPU and 1 GPU).

## Software Requirements

### OS Requirements

This package has been tested on the following Operating System:

Linux: CentOS 6 2.6.32    
Windows: Windows 10 18368.778    
Mac OSX: We have not tested it on a Mac yet, but it should in principle work.     

### Matlab Version Requirements

This package has been tested with `Matlab` R2018b. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2018a or higher to test the data and source codes.

# Repositary Contents

### 1. Experiment Data

Folder: [1_Measured_data](./1_Measured_data)

This folder contains experimental projections after denoising and alignment as well as their corresponding angles for Ta film and two Pd nanoparticles.

### 2. The REal Space Iterative REconstruction (RESIRE) Package

Folder: [2_RESIRE_package](./2_RESIRE_package)

Run the sample codes `Main_RESIRE_film_example.m` to get the 3D reconstruction of a smaller test film object. Run the main codes `Main_RESIRE_Pd1.m`, `Main_RESIRE_Pd2.m` and `Main_RESIRE_Ta_film_sAET.m` to obtain the 3D reconstruction of the two Pd nanoparticles and the sAET reconstruction of the Ta film.

### 3. Reconstructed 3D Volume

Folder: [3_Final_reconstruction_volume](./3_Final_reconstruction_volume)

This folder contains the 3D volume of the two Pd nanoparticles and Ta film reconstructed from above step.

### 4. Atom Tracing

Folder: [4_Atom_tracing](./4_Atom_tracing)

Run the codes `Main_atom_tracing_Pd1_nanoparticle.m`, `Main_atom_tracing_Pd2_nanoparticle.m` and `Main_atom_tracing_Ta_film.m` to trace the initial atomic positions from the reconstructed 3D volume.

Run the codes `Main_remove_non_atom_peak_Pd1_nanoparticle.m`, `Main_remove_non_atom_peak_Pd2_nanoparticle.m` and `Main_remove_non_atom_peak_Ta_film.m` to remove non-atom peaks using K-mean method. The obtained atomic coordinates were manually inspected before sending to position refinement.

### 5. Atomic Position Refinement

Folder: [5_Position_refinement](./5_Position_refinement)

Run the codes `Main_position_refinement_Pd1_nanoparticle.m`, `Main_position_refinement_Pd2_nanoparticle.m` and `Main_position_refinement_Ta_film.m` to refine the 3D atomic coordinates of each sample.

### 6. Experimental Atomic Model

Folder: [6_Final_coordinates](./6_Final_coordinates)

The final 3D atomic models of the two Pd nanoparticles and the Ta film.

### 7. Post Data Analysis

Folder: [7_Data_analysis](./7_Data_analysis)

Run the codes `Main_1_rdf_and_boo_entire_sample_Pd1_nanoparticle.m`, `Main_1_rdf_and_boo_entire_sample_Pd2_nanoparticle.m` and `Main_1_rdf_and_boo_entire_sample_Ta_film.m` to calculate the radial distribution function and the bond orientation order parameter for all the atoms in the three samples.

Run the codes `Main_2_rdf_and_voronoi_amorphous_region_Pd1_nanoparticle.m`, `Main_2_rdf_and_voronoi_amorphous_region_Pd2_nanoparticle.m` and `Main_2_rdf_and_voronoi_amorphous_region_Ta_film.m` to compute the radial distribution function and Voronoi indices in the amorphous region for all three samples.

Run the codes `Main_3_polytetrahedral_analysis_Pd1_nanoparticle.m`, `Main_3_polytetrahedral_analysis_Pd2_nanoparticle.m` and `Main_3_polytetrahedral_analysis_Ta_film.m` to search for polytetrahedral motifs, i.e. triplets, quadrilateral, pentagonal and hexagonal bipyramids, in the three samples.

Run the codes `Main_4_motifs_connection_and_network_Pd1_nanoparticle.m`, `Main_4_motifs_connection_and_network_Pd2_nanoparticle.m` and `Main_4_motifs_connection_and_network_Ta_film.m` to analyze the edge-/vertex-sharing and icosahedral filling by pentagonal bipyramids, as well as the networks formed by quadrilateral/pentagonal/hexagonal bipyramids in each sample.
