library(dplyr)
library(stringr)
library(tidyr)

# Read custom dataset from .csv
# Is expected to follow a wide format with the subject's age and sex, along with 60 measurement points per segment (ie. 60 point x 18 segments with 6 segments per view)
path <- "path/to/your/dataset"
df <- read.csv(path)

path_to_ref_ranges <- "path/to/directory/of/Reference_ranges.csv"

ref_ranges <- read.csv(path_to_ref_ranges)

#df <- df[,c(1:1080, which(str_detect(names(df), pattern = "K0ALDER5|KON")))]
names(df) <- str_remove_all(names(df), pattern = "pol_")

# Define the variable names for sex and age as str to be used to call for these later in the script
# Expects sex_var to be coded as binary with 0 = female and 1 = male
age_var <- "NAME_OF_AGE_VARIABLE"
sex_var <- "NAME_OF_SEX_VARIABLE"

# Remove subjects with missing data for sex and age
df <- df[-which(is.na(df[,age_var]) | is.na(df[,sex_var])),]

# Create age_group variable to sort subjects into age categories with same intervals as SCORE2
df$age_group <- NA
df[which(df[,age_var]<40),"age_group"] <- "<40"
df[which(between(df[,age_var],40,44.99)),"age_group"] <- "40-44"
df[which(between(df[,age_var],45,49.99)),"age_group"] <- "45-49"
df[which(between(df[,age_var],50,54.99)),"age_group"] <- "50-54"
df[which(between(df[,age_var],55,59.99)),"age_group"] <- "55-59"
df[which(between(df[,age_var],60,64.99)),"age_group"] <- "60-64"
df[which(df[,age_var]>=65),"age_group"] <- "\u226565"
age_groups <- c("<40", "40-44", "45-49", "50-54",
                "55-59", "60-64","\u226565")

# Define name of segments
# The script expects the name of the variables to follow this naming 
# Each of the 60 columns per segment should be named after the segment followed by an index, for instance apLat1, apLat2 etc
four_ch_names <- c('apLat', 'midLat', 'basLat',  # 4ch lateral wall
                    'apSept', 'midSept', 'basSept' # 4ch septal wall
              )

two_ch_names <- c('apAnt', 'midAnt', 'basAnt', # 2ch anterior wall
                  'apInf', 'midInf', 'basInf' # 2ch inferior wall
                  )

three_ch_names <- c('apAntSept', 'midAntSept', 'basAntSept', # 3ch/APLAX anterior septal wall
                    'apPost', 'midPost', 'basPost' # 3ch/APLAX posterior wall
                    )

# Create per-chamber mean curves
paste_collapse <- function(x){
  paste0(x, collapse = "")
}

df <- df %>% mutate(row_id = row_number())

cols <- c("row_id",
          paste0(four_ch_names, rep(1:60,each=6)),
          paste0(two_ch_names, rep(1:60,each=6)),
          paste0(three_ch_names, rep(1:60,each=6)))


# --- Step 1: Reshape the data from wide to long format ---
# This gathers all the measurement columns into a few key-value columns
long_data <- df[,cols] %>%
  pivot_longer(
    cols = -row_id, # Keep row_id as an identifier
    names_to = "column_name",
    values_to = "value"
  )

# --- Step 2: Extract the view (4ch, 2ch, 3ch) and group (j) from the column names ---
# Create a mapping data frame to join against
col_map <- data.frame(column_name = cols[-1]) %>%
  mutate(
    view = rep(c("4ch", "2ch", "3ch"), each = 60 * 6),
    group_j = rep(rep(1:60, each = 6), 3)
  )

# Join the mapping info to our long data
long_data_mapped <- left_join(long_data, col_map, by = "column_name")

# --- Step 3: Group by row, view, and j, then calculate the means ---
main_strain_long <- long_data_mapped %>%
  group_by(row_id, view, group_j) %>%
  summarise(
    mean_val = if_else(sum(is.na(value)) > 3, 
                       NA_real_, 
                       mean(value, na.rm = TRUE)),
    .groups = 'drop' # Drop grouping after summarising
  )

# --- Step 4: Pivot the data back to the desired wide format ---
main_strain <- main_strain_long %>%
  pivot_wider(
    id_cols = row_id,
    names_from = c(view, group_j),
    names_sep = "_",
    values_from = mean_val
  ) %>%
  select(-row_id) # Remove the helper row_id column

plot(1:60,main_strain[1,1:60])

# Columns are sorted by name with 2ch appearing first followed by 3ch and 4ch
# Change order of columns to 4ch, 2ch and 3ch
main_strain <- as.data.frame(main_strain[,c(which(str_detect(names(main_strain), pattern = "4ch")),
                                            which(str_detect(names(main_strain), pattern = "2ch")),
                                            which(str_detect(names(main_strain), pattern = "3ch")))])

# Ensure that the names of main_strain follow 4ch_1 - 4ch_60, then 2ch_1 - 2ch_60 etc.
names(main_strain)

# Calculate strain deviation measures
percent_within_normal_4ch <- c()
diastolic_percent_within_normal_4ch <- c()
diastolic_outfall_4ch <- c()
percent_within_normal_2ch <- c()
diastolic_percent_within_normal_2ch <- c()
diastolic_outfall_2ch <- c()
percent_within_normal_3ch <- c()
diastolic_percent_within_normal_3ch <- c()
diastolic_outfall_3ch <- c()
complete_strain_deviation_4ch <- c()
complete_strain_deviation_2ch <- c()
complete_strain_deviation_3ch <- c()
deviation_from_mean_4ch <- c()
deviation_from_mean_2ch <- c()
deviation_from_mean_3ch <- c()
i <- 1
for(i in 1:nrow(main_strain)){
  print(paste("Person",i))
  temp_4ch <- c()
  temp_diastolic_4ch <- c()
  temp_diastolic_size_4ch <- c()
  temp_complete_strain_deviation_4ch <- c()
  
  temp_2ch <- c()
  temp_diastolic_2ch <- c()
  temp_diastolic_size_2ch <- c()
  temp_complete_strain_deviation_2ch <- c()
  
  temp_3ch <- c()
  temp_diastolic_3ch <- c()
  temp_diastolic_size_3ch <- c()
  temp_complete_strain_deviation_3ch <- c()
  
  sex <- ifelse(df[i,sex_var]==0,"Female","Male")
  age_group <- df[i,"age_group"]
  rows <- which(str_detect(ref_ranges$age_group, pattern = age_group) & ref_ranges$sex==sex)
  gls_index_4ch <- which.min(main_strain[i,paste0("4ch_",1:60)])
  gls_index_2ch <- which.min(main_strain[i,paste0("2ch_",1:60)])
  gls_index_3ch <- which.min(main_strain[i,paste0("3ch_",1:60)])

  gls_4ch <- min(main_strain[i,paste0("4ch_",1:60)])
  gls_2ch <- min(main_strain[i,paste0("2ch_",1:60)])
  gls_3ch <- min(main_strain[i,paste0("3ch_",1:60)])
  
  gls_4ch_normal <- min(ref_ranges[rows,"mean_4ch"])
  gls_2ch_normal <- min(ref_ranges[rows,"mean_2ch"])
  gls_3ch_normal <- min(ref_ranges[rows,"mean_3ch"])
  
  factor_4ch <- gls_4ch/gls_4ch_normal
  factor_2ch <- gls_2ch/gls_2ch_normal
  factor_3ch <- gls_3ch/gls_3ch_normal
  
  for(j in 1:60){
    row <- rows[j]
    
    # 4ch
    lower_bound_4ch <- ref_ranges[row,"conf.low_4ch"]
    upper_bound_4ch <- ref_ranges[row,"conf.high_4ch"]
    mean_point_4ch <- ref_ranges[row,"mean_4ch"]
    temp_strain_point_4ch <- main_strain[i,paste0("4ch_",j)]
    
    within_normal_4ch <- ifelse(dplyr::between(temp_strain_point_4ch,
                                               left = lower_bound_4ch,
                                               right = upper_bound_4ch),1,0)
    
    temp_4ch <- c(temp_4ch,within_normal_4ch)
    
    # 2ch
    lower_bound_2ch <- ref_ranges[row,"conf.low_2ch"]
    upper_bound_2ch <- ref_ranges[row,"conf.high_2ch"]
    mean_point_2ch <- ref_ranges[row,"mean_2ch"]
    temp_strain_point_2ch <- main_strain[i,paste0("2ch_",j)]
    
    within_normal_2ch <- ifelse(dplyr::between(temp_strain_point_2ch,
                                               left = lower_bound_2ch,
                                               right = upper_bound_2ch),1,0)
    
    temp_2ch <- c(temp_2ch,within_normal_2ch)
    
    # 3ch
    lower_bound_3ch <- ref_ranges[row,"conf.low_3ch"]
    upper_bound_3ch <- ref_ranges[row,"conf.high_3ch"]
    mean_point_3ch <- ref_ranges[row,"mean_3ch"]
    temp_strain_point_3ch <- main_strain[i,paste0("3ch_",j)]
    
    within_normal_3ch <- ifelse(dplyr::between(temp_strain_point_3ch,
                                               left = lower_bound_3ch,
                                               right = upper_bound_3ch),1,0)
    
    temp_3ch <- c(temp_3ch,within_normal_3ch)
    
    if(is.na(temp_strain_point_4ch)){
      temp_4ch <- NA
      temp_diastolic_4ch <- NA
      temp_diastolic_size_4ch <- NA
    }else{
      # complete strain deviation
      if(dplyr::between(temp_strain_point_4ch/factor_4ch,
                        left = lower_bound_4ch,
                        right = upper_bound_4ch)){
        temp_complete_strain_deviation_4ch <- c(temp_complete_strain_deviation_4ch,0)
      }else{
        #how_far_off <- abs(temp_strain_point-mean_point)
        how_far_off_4ch <- ifelse((temp_strain_point_4ch/factor_4ch)<lower_bound_4ch, 
                                  abs((temp_strain_point_4ch/factor_4ch)-lower_bound_4ch),
                                  abs((temp_strain_point_4ch/factor_4ch)-upper_bound_4ch))
        temp_complete_strain_deviation_4ch <- c(temp_complete_strain_deviation_4ch,how_far_off_4ch)
      }
      
      if(j>=gls_index_4ch){
        within_normal_4ch <- ifelse(dplyr::between(temp_strain_point_4ch,
                                                   left = lower_bound_4ch,
                                                   right = upper_bound_4ch),1,0)
        temp_diastolic_4ch <- c(temp_diastolic_4ch, within_normal_4ch)
        if(dplyr::between(temp_strain_point_4ch,
                          left = lower_bound_4ch,
                          right = upper_bound_4ch)){
          temp_diastolic_size_4ch <- c(temp_diastolic_size_4ch,0)
        }else{
          #how_far_off <- abs(temp_strain_point-mean_point)
          how_far_off_4ch <- ifelse(temp_strain_point_4ch<lower_bound_4ch, 
                                    abs(temp_strain_point_4ch-lower_bound_4ch),
                                    abs(temp_strain_point_4ch-upper_bound_4ch))
          temp_diastolic_size_4ch <- c(temp_diastolic_size_4ch,how_far_off_4ch)
        }
      }
    }
    
    if(is.na(temp_strain_point_2ch)){
      temp_2ch <- NA
      temp_diastolic_2ch <- NA
      temp_diastolic_size_2ch <- NA
    }else{
      # complete strain deviation
      if(dplyr::between(temp_strain_point_2ch/factor_2ch,
                        left = lower_bound_2ch,
                        right = upper_bound_2ch)){
        temp_complete_strain_deviation_2ch <- c(temp_complete_strain_deviation_2ch,0)
      }else{
        #how_far_off <- abs(temp_strain_point-mean_point)
        how_far_off_2ch <- ifelse((temp_strain_point_2ch/factor_2ch)<lower_bound_2ch, 
                                  abs((temp_strain_point_2ch/factor_2ch)-lower_bound_2ch),
                                  abs((temp_strain_point_2ch/factor_2ch)-upper_bound_2ch))
        temp_complete_strain_deviation_2ch <- c(temp_complete_strain_deviation_2ch,how_far_off_2ch)
      }
      
      if(j>=gls_index_2ch){
        within_normal_2ch <- ifelse(dplyr::between(temp_strain_point_2ch,
                                                   left = lower_bound_2ch,
                                                   right = upper_bound_2ch),1,0)
        temp_diastolic_2ch <- c(temp_diastolic_2ch, within_normal_2ch)
        if(dplyr::between(temp_strain_point_2ch,
                          left = lower_bound_2ch,
                          right = upper_bound_2ch)){
          temp_diastolic_size_2ch <- c(temp_diastolic_size_2ch,0)
        }else{
          #how_far_off <- abs(temp_strain_point-mean_point)
          how_far_off_2ch <- ifelse(temp_strain_point_2ch<lower_bound_2ch, 
                                    abs(temp_strain_point_2ch-lower_bound_2ch),
                                    abs(temp_strain_point_2ch-upper_bound_2ch))
          temp_diastolic_size_2ch <- c(temp_diastolic_size_2ch,how_far_off_2ch)
        }
      }
    }
    
    if(is.na(temp_strain_point_3ch)){
      temp_3ch <- NA
      temp_diastolic_3ch <- NA
      temp_diastolic_size_3ch <- NA
    }else{
      # complete strain deviation
      if(dplyr::between(temp_strain_point_3ch/factor_3ch,
                        left = lower_bound_3ch,
                        right = upper_bound_3ch)){
        temp_complete_strain_deviation_3ch <- c(temp_complete_strain_deviation_3ch,0)
      }else{
        #how_far_off <- abs(temp_strain_point-mean_point)
        how_far_off_3ch <- ifelse((temp_strain_point_3ch/factor_3ch)<lower_bound_3ch, 
                                  abs((temp_strain_point_3ch/factor_3ch)-lower_bound_3ch),
                                  abs((temp_strain_point_3ch/factor_3ch)-upper_bound_3ch))
        temp_complete_strain_deviation_3ch <- c(temp_complete_strain_deviation_3ch,how_far_off_3ch)
      }
      
      if(j>=gls_index_3ch){
        within_normal_3ch <- ifelse(dplyr::between(temp_strain_point_3ch,
                                                   left = lower_bound_3ch,
                                                   right = upper_bound_3ch),1,0)
        temp_diastolic_3ch <- c(temp_diastolic_3ch, within_normal_3ch)
        if(dplyr::between(temp_strain_point_3ch,
                          left = lower_bound_3ch,
                          right = upper_bound_3ch)){
          temp_diastolic_size_3ch <- c(temp_diastolic_size_3ch,0)
        }else{
          #how_far_off <- abs(temp_strain_point-mean_point)
          how_far_off_3ch <- ifelse(temp_strain_point_3ch<lower_bound_3ch, 
                                    abs(temp_strain_point_3ch-lower_bound_3ch),
                                    abs(temp_strain_point_3ch-upper_bound_3ch))
          temp_diastolic_size_3ch <- c(temp_diastolic_size_3ch,how_far_off_3ch)
        }
      }
    }
  }
  
  deviation_from_mean_4ch <- c(deviation_from_mean_4ch,sum(abs(main_strain[i,paste0("4ch_",1:60)] - ref_ranges[rows,"mean_4ch"])))
  deviation_from_mean_2ch <- c(deviation_from_mean_2ch,sum(abs(main_strain[i,paste0("2ch_",1:60)] - ref_ranges[rows,"mean_2ch"])))
  deviation_from_mean_3ch <- c(deviation_from_mean_3ch,sum(abs(main_strain[i,paste0("3ch_",1:60)] - ref_ranges[rows,"mean_3ch"])))
  
  percent_4ch <- sum(temp_4ch)/60*100
  percent_within_normal_4ch <- c(percent_within_normal_4ch,percent_4ch)
  percent_diastolic_4ch <- sum(temp_diastolic_4ch)/length(temp_diastolic_4ch)*100
  diastolic_percent_within_normal_4ch <- c(diastolic_percent_within_normal_4ch,percent_diastolic_4ch)
  diastolic_outfall_4ch <- c(diastolic_outfall_4ch,sum(temp_diastolic_size_4ch))
  complete_strain_deviation_4ch <- c(complete_strain_deviation_4ch, sum(temp_complete_strain_deviation_4ch))
  
  percent_2ch <- sum(temp_2ch)/60*100
  percent_within_normal_2ch <- c(percent_within_normal_2ch,percent_2ch)
  percent_diastolic_2ch <- sum(temp_diastolic_2ch)/length(temp_diastolic_2ch)*100
  diastolic_percent_within_normal_2ch <- c(diastolic_percent_within_normal_2ch,percent_diastolic_2ch)
  diastolic_outfall_2ch <- c(diastolic_outfall_2ch,sum(temp_diastolic_size_2ch))
  complete_strain_deviation_2ch <- c(complete_strain_deviation_2ch, sum(temp_complete_strain_deviation_2ch))
  
  percent_3ch <- sum(temp_3ch)/60*100
  percent_within_normal_3ch <- c(percent_within_normal_3ch,percent_3ch)
  percent_diastolic_3ch <- sum(temp_diastolic_3ch)/length(temp_diastolic_3ch)*100
  diastolic_percent_within_normal_3ch <- c(diastolic_percent_within_normal_3ch,percent_diastolic_3ch)
  diastolic_outfall_3ch <- c(diastolic_outfall_3ch,sum(temp_diastolic_size_3ch))
  complete_strain_deviation_3ch <- c(complete_strain_deviation_3ch, sum(temp_complete_strain_deviation_3ch))
}

deviation_from_mean <- rowMeans(data.frame(deviation_from_mean_4ch,
                                           deviation_from_mean_2ch,
                                           deviation_from_mean_3ch),na.rm = T)

percent_within_normal <- rowMeans(data.frame(percent_within_normal_4ch,
                                             percent_within_normal_2ch,
                                             percent_within_normal_3ch),na.rm = T)

diastolic_percent_within_normal <- rowMeans(data.frame(diastolic_percent_within_normal_4ch,
                                                       diastolic_percent_within_normal_2ch,
                                                       diastolic_percent_within_normal_3ch),na.rm = T)

diastolic_outfall <- rowMeans(data.frame(diastolic_outfall_4ch,
                                         diastolic_outfall_2ch,
                                         diastolic_outfall_3ch),na.rm = T)

complete_strain_deviation <- rowMeans(data.frame(complete_strain_deviation_4ch,
                                                 complete_strain_deviation_2ch,
                                                 complete_strain_deviation_3ch),na.rm = T)

df$deviation_from_mean <- deviation_from_mean
df$strain_deviance <- percent_within_normal
df$diastolic_strain_deviation <- diastolic_outfall
df$indexed_strain_deviation <- complete_strain_deviation
