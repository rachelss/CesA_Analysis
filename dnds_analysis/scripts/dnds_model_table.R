library(tidyverse)   
library(flextable) # Required for table creation
library(officer) # Required for border formatting

# Data frame
dnds_table <- tibble::tibble(
  Model = c("M0 (one-ratio)", "M1a (nearly neutral)", "M1a (nearly neutral)", 
            "M2a (positive selection)", "M7 (beta)", "M8 (beta & ω)"),
  LogL = c("-59147.35", "-58502.65", "-58502.65", "-58502.65", "-57030.17", "-57030.18"),
  Free_parameters = c(98, 99, 99, 101, 99, 101),
  Comparison = c("M0 vs. M1a", "", "M1a vs. M2a", "", "M7 vs. M8", ""),
  df = c(1, NA, 2, NA, 2, NA),
  Delta_LogL = c(1289.4, NA, 0, NA, -0.02, NA),
  p_value = c("2.27 × 10⁻²⁸²", "", "1", "", "1", "")
)

# Create flextable
ft <- flextable(dnds_table) %>%
  set_header_labels(
    Model = "Model",
    LogL = "ℓ",
    Free_parameters = "Free parameters",
    Comparison = "Model comparison",
    df = "df",
    Delta_LogL = "LR (2Δℓ)",
    p_value = "p-value"
  ) %>%
  theme_vanilla() %>%
  bg(bg = "white", part = "all") %>%
  color(color = "black", part = "all") %>%
  align(align = "center", part = "all") %>%
  
  # Merge rows
  merge_at(j = "Comparison", i = 1:2) %>%
  merge_at(j = "Comparison", i = 3:4) %>%
  merge_at(j = "Comparison", i = 5:6) %>%
  merge_at(j = "df", i = 1:2) %>%
  merge_at(j = "df", i = 3:4) %>%
  merge_at(j = "df", i = 5:6) %>%
  merge_at(j = "Delta_LogL", i = 1:2) %>%
  merge_at(j = "Delta_LogL", i = 3:4) %>%
  merge_at(j = "Delta_LogL", i = 5:6) %>%
  merge_at(j = "p_value", i = 1:2) %>%
  merge_at(j = "p_value", i = 3:4) %>%
  merge_at(j = "p_value", i = 5:6) %>%
  
  # Font sizes & spacing
  fontsize(size = 17, part = "header") %>%
  fontsize(size = 16, part = "body") %>%
  line_spacing(space = 1.5, part = "all") %>%
  
  # Add column separators (fixed)
  border(j = seq_along(dnds_table), border = officer::fp_border(color = "grey10", width = 1)) %>%
  
  # Add dividing lines between header columns
  border(part = "header", 
         border.left = officer::fp_border(color = "grey10", width = 1), 
         border.right = officer::fp_border(color = "grey10", width = 1)) %>%
  
  # Set column widths
  width(j = 1, width = 2.5) %>%
  width(j = 2:7, width = 1.5)

# Save to PNG
save_as_image(ft, path = "dnds_site_model_table.png", webshot = "webshot2", zoom = 3, expand = 10)


