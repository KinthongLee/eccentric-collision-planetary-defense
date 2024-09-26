# eccentric-collision-planetary-defense
This code primarily generates the original code used in the Icarus paper titled "Investigation of the Incremental Benefits of Eccentric Collisions in Kinetic Deflection of Potentially Hazardous Asteroids." It includes all computational models and result processing. Everyone is welcome to use and modify this code! Just remember to reference the paper! If you publish your work, please cite my article! I would greatly appreciate it!

Additionally, I am currently conducting research on asteroid defense at Tsinghua University in China and welcome any potential collaborations! You can reach me at: ktlee3819@gmail.com.

% ---------------------------------------------------------

The entire code is primarily divided into two parts: 1. Calculation 2. Plotting.

%---------------- Calculation------------------------------

All code will use the PHA_table.xlsx file as the central data storage. The code will read data from this file, and once the program completes, it will write the results back into the table and output a temporary_result.xlsx file. After confirming the results for the processed target are correct, you may delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.


This section will takes PHA row.34 (2023 DZ2) as an example to demonstrate the entire code execution process.

The data in PHA_table.xlsx comes from: https://cneos.jpl.nasa.gov/ca/ and is sorted in descending order based on Rarity. The original data is from columns A to I in the table (from 'Object' to 'Rarity'). The remaining columns are initially blank and will be filled in with results from subsequent code calculations. From the table, it can be seen that I have completed the PHA calculations up to Row 33, which are also the objects discussed in my paper.

1. Before begin simulation for any PHA, go to the NASA JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/app.html#/ . First, in the Ephemeris Type, select "Small-Body SPK File." In the Target Body, enter the name of the target PHA; for this example, use 2022 QX4. For the Time Specification, make sure to set the range to plus and minus 20 years from the PHA's close-approach year. There are no strict limits; you can choose a larger range as long as your computer has enough storage space. If the range is too narrow and exceeds the limit during calculations, an error will occur. For 2022 QX4, the close-approach year is 1977, so set Start to 1950-01-01 and Stop to 2000-01-01. Then, click "Generate Ephemeris" to download the corresponding .bsp file. Manually enter the corresponding BSP number for this PHA into Column J (BSP_file_name) of PHA_table.xlsx; in this case, it is 54297628. Next, place the 54297628.bsp file into the code folder: eccentric_collision_planetary_defense/mice/kernel.
   
2. Open the matlab code calculate_closest_approach_distance_and_relative_error.m and set line 87 to "p = 33:33". Here, p represents the PHA in row p+1 of PHA_table.xlsx. This example only calculates PHA 2022 QX4 (row 34), so set p to p = 33:33 and run the program. This program will use the High Precision Orbit Propagator that I have developed to perform recursion starting 10 years before the close approach of the PHA, calculating the closest approach (CA) distance to Earth and comparing it with the results from the SPICE model. If the relative error is acceptable, the process can continue, and the futher deflection distance will be calculated based on the CA distance obtained from this MATLAB model. After completion, a temporary_result.xlsx will be generated, and the previously blank Columns K and L will now contain results. Once you confirm that the results are correct, you can delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.xlsx.

3. Next, open the MATLAB file calculate_launch_window_and_best_transfer.m and change the code on line 96 to "for i = 33:33". Then run the program. This program will calculate all possible bi-impulsive transfer orbits, known as Lambert transfers, under the condition of a 10-year warning period before the close approach of 2022 QX4, as well as the deflection distance resulting from a collision with COG. The program will save all results in .mat file format in the directory eccentric_collision_planetary_defense\output_result\different_PHA\launch_window\matfile. These results will later be called upon in the plot section to create a launch window diagram. Additionally, the program will save the best Lambert transfer parameters to Columns M (Best_Launch Year) to AL (delta_t_by_T) in the table and output them to temporary_result.xlsx. As before, once you confirm that the results are correct, you can delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.xlsx.

4. Next, open the MATLAB file calculate_different_PHA_deflection_distance.m and change the code on line 117 to "for p = 33:33". Then run the program. This program will perform a Monte Carlo simulation with 100,000 samples (you can change the sample size on line 103 " sample_size = 100000;")to address uncertainties in the PHA's attitude and beta coefficient, the parameters will based on the optimal Lambert transfer parameters calculated earlier. The default 3D model used for the PHA is the light-curve 3D model of Apophis, but you can modify line 107 "OBJ = read_wobj('Apophis_Model.obj')" to replace it with any 3D model you need, make sure you store the .obj file in eccentric_collision_planetary_defense\code\3D_Model. The program will save all results in .mat file format in the directory eccentric_collision_planetary_defense\output_result\different_PHA\distribution\100k_result\matfile. These results will later be referenced in the plot section to create a distribution of deflection distance diagram. Additionally, the result parameters will be written to Columns AM (x_ast_impact) to BH (different_delta_v_h_relative) in the table and output to temporary_result.xlsx. As before, once you confirm that the results are correct, you can delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.xlsx.

5. Next, open the MATLAB file calculate_analytical_method and change line 10 to for ppp = 33:33. Then run the program. This program will use the analytical method to calculate the deflection distance, writing the results into Columns BI (delta_r_COG_theory) to BL (Gain_theory_percent) in the table, and output them to temporary_result.xlsx. As before, once you confirm that the results are correct, you can delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.xlsx.

At this point, all calculation processes are complete, and all important parameters are recorded in PHA_table.xlsx.

% ------------------------------------------------------------------------




%- --------------------- Plotting ---------------------------

Regarding the plotting aspect, the main figures produced by this model include:

1. Launch Window: Figures 7, 12, and 13 in the paper. Open the MATLAB file plot_launch_window.m. The resulting plots will be saved in eccentric_collision_planetary_defense\output_result\different_PHA\launch_window\pngfile.

2. Best Transfer Trajectory: Figures 8, 17, and 18 in the paper. Open the MATLAB file plot_transfer_trajectory.m. The resulting plots will be saved in eccentric_collision_planetary_defense\output_result\different_PHA\launch_window\transfer_orbit_png.

3. Impact Model Illustrator: Figure 6 in the paper. Open the MATLAB file plot_impact_model_illustrator.m.

4. Deflection Distance Distribution: Figures 10, 14, and 15 in the paper. Open the MATLAB file plot_deflection_distance_distribution.m. The resulting plots will be saved in eccentric_collision_planetary_defense\output_result\different_PHA\distribution\100k_result\pngfile.

The files I provided already contain the result plots for the first 32 PHAs in the corresponding locations, which you can check.

These plotting programs do not require significant modifications; just ensure that your target PHA is included in the for loop. For example, in plot_launch_window.m, line 21 reads for i = 1:32, which will generate and export the launch windows for PHAs in rows 2 to 33 of PHA_table to eccentric_collision_planetary_defense\output_result\different_PHA\launch_window\pngfile.

% -------------------------------------------------------------------
