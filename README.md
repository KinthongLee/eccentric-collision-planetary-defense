# eccentric-collision-planetary-defense
## Introduction
This code primarily generates the data and figures used in the Icarus paper titled "Investigation of the Incremental Benefits of Eccentric Collisions in Kinetic Deflection of Potentially Hazardous Asteroids." It includes all computational models and result processing. Everyone is welcome to use and modify this code! Just remember to reference the paper! If you publish your work, please cite my article! I would greatly appreciate it!

- Lee, Kinthong, Zhengqing Fang, and Zhaokui Wang. "Investigation of the incremental benefits of eccentric collisions in kinetic deflection of potentially hazardous asteroids." Icarus 425 (2025): 116312.

Additionally, I am currently conducting research on asteroid defense at Tsinghua University in China and welcome any potential collaborations! You can reach me at: ktlee3819@gmail.com.

---

## BEFORE YOU BEGIN

Make sure you download JPL planetary ephemerides ".bsp" files, specifically de441_part-1.bsp and de441_part-2.bsp (The lastest JPL planetary ephemerides on 26th Sept 2024) from https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/ or from BaiduNetdisk: : https://pan.baidu.com/s/1lEcd3QQUUZNuWDM-AmBgjg?pwd=6t5h passcode: 6t5h  (for chinese mainland users) and add it to location "eccentric_collision_planetary_defense/code/kernel". Otherwise the high precision orbit propagator will fails. Because both planetary ephemerides are way larger than github's limitation on files upload (100MB), both planetary ephemerides are leave for you to download separately.

---

## The entire code is primarily divided into two parts:  
1. Calculation  
2. Plotting  

---

## Calculation

All code will use the **`PHA_table.xlsx`** file as the central data storage. The code will read data from this file, and once the program completes, it will write the results back into the table and output a **`temporary_result.xlsx`** file. After confirming the results for the processed target are correct, you may delete the original **`PHA_table.xlsx`** and rename **`temporary_result.xlsx`** to **`PHA_table.xlsx`**.

---

### Example: Running the Code for PHA Row 34 (2023 DZ2)

This section will takes PHA row.34 (2023 DZ2) as an example to demonstrate the entire code execution process.

The data in **`PHA_table.xlsx`** comes from: https://cneos.jpl.nasa.gov/ca/ and is sorted in descending order based on Rarity. The original data is from columns A to I in the table (from 'Object' to 'Rarity'). The remaining columns are initially blank and will be filled in with results from subsequent code calculations. From the table, it can be seen that I have completed the PHA calculations up to Row 33, which are also the objects discussed in my paper.

---

#### Step 1. Generate SPK File from JPL Horizons

Before begin simulation for any PHA, go to the NASA JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/app.html#/ .  
- In the **Ephemeris Type**, select *Small-Body SPK File*.  
- In the **Target Body**, enter the name of the target PHA; for this example, use **2022 QX4**.  
- For the **Time Specification**, make sure to set the range to plus and minus 20 years from the PHA's close-approach year.  

For **2022 QX4**, the close-approach year is 1977, so set:  
- Start = **1950-01-01**  
- Stop = **2000-01-01**  

Click **Generate Ephemeris** to download the corresponding **.bsp** file.  
- Manually enter the corresponding BSP number for this PHA into Column J (**BSP_file_name**) of **PHA_table.xlsx**; in this case, it is **54297628**.  
- Next, place the **54297628.bsp** file into the folder:  
  **eccentric_collision_planetary_defense/code/kernel**

---

#### Step 2. Calculate Closest Approach Distance

Open the MATLAB code **`calculate_closest_approach_distance_and_relative_error.m`** and set line 87:

    p = 33:33

Here, `p` represents the PHA in row `p+1` of **PHA_table.xlsx**. This example only calculates PHA **2022 QX4** (row 34), so set `p = 33:33` and run the program.  

This program will:  
- Use the High Precision Orbit Propagator that I have developed to perform recursion starting 10 years before the close approach of the PHA.  
- Calculate the closest approach (CA) distance to Earth and compare it with the results from the SPICE model.  
- If the relative error is acceptable, the process can continue, and the further deflection distance will be calculated based on the CA distance obtained from this MATLAB model.  
- After completion, a **temporary_result.xlsx** will be generated, and the previously blank Columns K and L will now contain results.  

Once you confirm that the results are correct, you can delete the original **PHA_table.xlsx** and rename **temporary_result.xlsx** to **PHA_table.xlsx**.

---

#### Step 3. Calculate Launch Window and Best Transfer

Open the MATLAB file **`calculate_launch_window_and_best_transfer.m`** and change the code on line 96 to:

    for i = 33:33

Then run the program.  

This program will:  
- Calculate all possible two impulse transfer orbits, known as Lambert transfers, under the condition of a 10-year warning period before the close approach of **2022 QX4**.  
- Calculate the deflection distance resulting from a collision with COG.  
- Save all results in **.mat** file format in the directory:  
  **eccentric_collision_planetary_defense/output_result/different_PHA/launch_window/matfile**  
- Save the best Lambert transfer parameters to Columns M (**Best_Launch Year**) to AL (**delta_t_by_T**) in the table and output them to **temporary_result.xlsx**.  

As before, once you confirm that the results are correct, you can delete the original **PHA_table.xlsx** and rename **temporary_result.xlsx** to **PHA_table.xlsx**.

---

#### Step 4. Monte Carlo Deflection Simulation

Open the MATLAB file **`calculate_different_PHA_deflection_distance.m`** and change the code on line 117 to:

    for p = 33:33

You can also change the sample size on line 103:

    sample_size = 100000;

This program will:  
- Perform a Monte Carlo simulation with **100,000 samples** to address uncertainties in the PHA's attitude and beta coefficient, based on the optimal Lambert transfer parameters calculated earlier.  
- Use the default 3D model of Apophis (line 107):  

    OBJ = read_wobj('Apophis_Model.obj')

- You may replace it with any 3D model you need, just ensure the **.obj** file is stored in:  
  **eccentric_collision_planetary_defense/code/3D_Model**  
- Save all results in **.mat** file format in:  
  **eccentric_collision_planetary_defense/output_result/different_PHA/distribution/100k_result/matfile**  
- Write the result parameters to Columns AM (**x_ast_impact**) to BH (**different_delta_v_h_relative**) in the table and output them to **temporary_result.xlsx**.  

As before, once you confirm that the results are correct, you can delete the original **PHA_table.xlsx** and rename **temporary_result.xlsx** to **PHA_table.xlsx**.

---

#### Step 5. Analytical Method Calculation

Open the MATLAB file **`calculate_analytical_method.m`** and change line 10 to:

    for ppp = 33:33

Then run the program.  

This program will:  
- Use the analytical method to calculate the deflection distance.  
- Write the results into Columns BI (**delta_r_COG_theory**) to BL (**Gain_theory_percent**) in the table.  
- Output them to **temporary_result.xlsx**.  

As before, once you confirm that the results are correct, you can delete the original **PHA_table.xlsx** and rename **temporary_result.xlsx** to **PHA_table.xlsx**.

---

#### Final Step

At this point, all calculation processes are complete, and all important parameters are recorded in **PHA_table.xlsx**.

---

## Plotting

Regarding the plotting aspect, the main figures produced by this model include:

1. **Launch Window**  
   - Figures 7, 12, and 13 in the paper.  
   - Open the MATLAB file **`plot_launch_window.m`**.  
   - The resulting plots will be saved in:  
     **eccentric_collision_planetary_defense/output_result/different_PHA/launch_window/pngfile**

2. **Best Transfer Trajectory**  
   - Figures 8, 17, and 18 in the paper.  
   - Open the MATLAB file **`plot_transfer_trajectory.m`**.  
   - The resulting plots will be saved in:  
     **eccentric_collision_planetary_defense/output_result/different_PHA/launch_window/transfer_orbit_png**

3. **Impact Model Illustrator**  
   - Figure 6 in the paper.  
   - Open the MATLAB file **`plot_impact_model_illustrator.m`**.

4. **Deflection Distance Distribution**  
   - Figures 10, 14, and 15 in the paper.  
   - Open the MATLAB file **`plot_deflection_distance_distribution.m`**.  
   - The resulting plots will be saved in:  
     **eccentric_collision_planetary_defense/output_result/different_PHA/distribution/100k_result/pngfile**

---

The files I provided already contain the result plots for the first 32 PHAs in the corresponding locations, which you can check.

These plotting programs do not require significant modifications; just ensure that your target PHA is included in the `for` loop.  
For example, in **`plot_launch_window.m`**, line 21 reads:

    for i = 1:32

This will generate and export the launch windows for PHAs in rows 2 to 33 of **PHA_table.xlsx** to:  
**eccentric_collision_planetary_defense/output_result/different_PHA/launch_window/pngfile**
