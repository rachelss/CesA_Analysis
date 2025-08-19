library(writexl)
library(dplyr)

# Add significance columns
branch_model_lrt_summary_with_omega$significant <- branch_model_lrt_summary_with_omega$p_value < 0.05
branch_model_lrt_summary_with_omega$FDR <- p.adjust(branch_model_lrt_summary_with_omega$p_value, method = "BH")
branch_model_lrt_summary_with_omega$FDR_significant <- branch_model_lrt_summary_with_omega$FDR < 0.05

# Count significant hits
sum(branch_model_lrt_summary_with_omega$FDR_significant)

# Subset only significant rows
significant_df <- branch_model_lrt_summary_with_omega %>%
  filter(FDR_significant == TRUE)

# Write full table
write_xlsx(branch_model_lrt_summary_with_omega, "Supplementary_branch_model_summary.xlsx")

# Save to full table to csv as well
write.csv(branch_model_lrt_summary_with_omega, "Supplementary_branch_model_summary.csv", row.names = FALSE)



