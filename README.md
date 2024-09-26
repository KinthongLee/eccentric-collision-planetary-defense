# eccentric-collision-planetary-defense
This code primarily generates the original code used in the Icarus paper titled "Investigation of the Incremental Benefits of Eccentric Collisions in Kinetic Deflection of Potentially Hazardous Asteroids." It includes all computational models and result processing. Everyone is welcome to use and modify this code! Just remember to reference the paper! If you publish your work, please cite my article! I would greatly appreciate it!

Additionally, I am currently conducting research on asteroid defense at Tsinghua University in China and welcome any potential collaborations! You can reach me at: ktlee3819@gmail.com.

% --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



All code will use the PHA_table.xlsx file as the central data storage. The code will read data from this file, and once the program completes, it will write the results back into the table and output a temporary_result.xlsx file. After confirming the results for the processed target are correct, you may delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.


This section will takes PHA row.34 (2023 DZ2) as an example to demonstrate the entire code execution process.

The data in PHA_table.xlsx comes from: https://cneos.jpl.nasa.gov/ca/ and is sorted in descending order based on Rarity. The original data is from columns A to I in the table (from 'Object' to 'Rarity'). The remaining columns are initially blank and will be filled in with results from subsequent code calculations. From the table, it can be seen that I have completed the PHA calculations up to Row 33, which are also the objects discussed in my paper.

1. Before begin simulation for any PHA, go to the NASA JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/app.html#/ . First, in the Ephemeris Type, select "Small-Body SPK File." In the Target Body, enter the name of the target PHA; for this example, use 2022 QX4. For the Time Specification, make sure to set the range to plus and minus 20 years from the PHA's close-approach year. There are no strict limits; you can choose a larger range as long as your computer has enough storage space. If the range is too narrow and exceeds the limit during calculations, an error will occur. For 2022 QX4, the close-approach year is 1977, so set Start to 1950-01-01 and Stop to 2000-01-01. Then, click "Generate Ephemeris" to download the corresponding .bsp file. Manually enter the corresponding BSP number for this PHA into Column J (BSP_file_name) of PHA_table.xlsx; in this case, it is 54297628. Next, place the 54297628.bsp file into the code folder: High_Precision_Orbit_Propagator/mice/kernel.
   
2. Open the matlab code calculate_closest_approach_distance_and_relative_error.m and set line 87 to p = 33:33. Here, p represents the PHA in row p+1 of PHA_table.xlsx. This example only calculates PHA 2022 QX4 (row 34), so set p to p = 33:33 and run the program. After completion, a temporary_result.xlsx will be generated, and the previously blank Columns K and L will now contain results. Once you confirm that the results are correct, you can delete the original PHA_table.xlsx and rename temporary_result.xlsx to PHA_table.xlsx.

