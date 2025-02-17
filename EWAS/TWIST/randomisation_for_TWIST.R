twist = read.csv("<datastore_home_dir>/Lab/Smoking/TWIST/2022-07-28-Twist_Sample_Info.csv")
data = twist[c("Sample_Name", "age", "pair", "Case", "sex")]

# Convert Group column to factor for easier filtering
data$Group <- as.factor(data$Case)

# Separate controls and cases based on Group column
controls <- data$age[data$Group == 0]
cases <- data$age[data$Group == 1]

# Find the minimum and maximum ages in each group
min_control_age <- min(controls)
max_control_age <- max(controls)
min_case_age <- min(cases)
max_case_age <- max(cases)

# Calculate the maximum difference
max_age_difference <- max(abs(max_control_age - min_case_age), abs(max_case_age - min_control_age))

# Print the result
cat("The maximum age difference between controls and cases is:", max_age_difference, "\n")

### within pair

# Filter pairs: select only those rows with labels that contain both "cont" and "case" in matched pairs
pairs <- data[grep("cont|case", data$pair), ]

# Sort data for pairing
pairs <- pairs[order(pairs$pair), ]

# Ensure the controls and cases alternate in the sorted data frame
controls <- pairs[pairs$Group == 0, ]
cases <- pairs[pairs$Group == 1, ]

# Check if we have an equal number of controls and cases to pair
if (nrow(controls) == nrow(cases)) {
  
  # Calculate the absolute age difference within each control-case pair
  age_diffs <- abs(controls$age - cases$age)
  
  # Combine into a data frame to display the pair differences
  pair_diffs <- data.frame(
    Control_Id = controls$Sample_Name, 
    Control = controls$pair,
    Control_Sex = controls$sex,
    Case_Id = cases$Sample_Name,
    Case = cases$pair,
    Case_Sex = cases$sex,
    Control_Age = controls$age,
    Case_Age = cases$age,
    Age_Difference = age_diffs
  )
  
  # Calculate the maximum difference within pairs
  max_pair_diff <- max(age_diffs)
  identical(pair_diffs$Control_Sex, pair_diffs$Case_Sex)
  # Output the pair differences and the maximum difference
  print(pair_diffs) # 6.67
  write.csv(pair_diffs, "<datastore_home_dir>/Lab/Smoking/TWIST/Pair_age_differences_TWIST.csv", row.names = F)
  
  
  wrong_ids = c(38010, 165495)
  cleaned_twist = subset(twist, !Sample_Name %in% wrong_ids)
  pair_diffs = subset(pair_diffs, !Case_Id %in% wrong_ids)
  
  
  identical(pair_diffs$Control_Sex, pair_diffs$Case_Sex) # TRUE
  table(pair_diffs$Age_Difference == 0) # TRUE
  
  write.csv(cleaned_twist, "<datastore_home_dir>/Lab/Smoking/TWIST/2024-11-14-Twist_Sample_Info_Cleaned.csv")
  
  cat("The maximum age difference within pairs is:", max_pair_diff, "\n")
  
} else {
  cat("The number of controls and cases does not match for pairing.\n")
}
