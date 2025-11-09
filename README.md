# Normal Standardized LV Longitudinal Strain Curves
This is the official repository related to the paper *Normal Standardized Strain Curves Stratified by Age and Sex in Healthy Individuals: The Copenhagen City Heart Study* by Christensen and Simonsen et al. 

Here, we provide the age- and sex-appropriate reference ranges of the LV standardized longitudinal strain curves. In addition, we provide the official code implementations to standardize custom longitudinal strain curves and calculate EDS, LDS along with age- and sex-appropriate mean strain deviation and diastolic strain deviation as detailed in the paper. 
<img width="2952" height="3544" alt="Figure 4" src="https://github.com/user-attachments/assets/a99e84d6-a1fd-4580-a95f-ef7e4a8e404d" />

# Script to derive mean strain curves and calculate novel strain deviation measures
The file "Mean_curves_and_strain_deviation.R" provides a script to derive a mean curve for each apical chamber view and subsequently calculate mean strain deviation and diastolic strain deviation for a custom dataset.
The dataset is expected to follow the structure of a resulting dataset from the script "Standard_Strain_master_30_09_2025.m". 
