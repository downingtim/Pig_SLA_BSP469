library(tidyverse)
library(viridis)
setwd("/mnt/lustre/RDS-ephemeral/bioinformatics/proj/bsp/BSP469/blast_outputs/")
# 1. Load data
matrix_file <- "Filtered_SLA_Allele_Matrix.185.csv"
df_raw <- read_csv(matrix_file, show_col_types = FALSE)

colnames(df_raw)
# get row sums
df_raw[, colnames(df_raw)] <- lapply(df_raw[, colnames(df_raw)], as.numeric)
a2 <- sort ( rowSums(df_raw, na.rm=T)  , decreasing = TRUE)
a2 <- a2[a2 > 199]
summary(a2)
IQR(a2)
# 2. Reshape to long format
df_long <- df_raw %>%  pivot_longer(
    cols = -sample,    names_to = "Allele",
    values_to = "ReadCount"  ) %>%
  mutate(Allele_Clean = gsub("IPD-MHC:", "", Allele))

str(df_long)

# 3. Apply Minimum Depth Filter (> 200 reads per sample)
valid_samples <- df_long %>%  group_by(sample) %>%
  summarize(TotalSampleReads = sum(ReadCount), .groups = "drop") %>%
  filter(TotalSampleReads >= 200) %>%  pull(sample)

df_filtered <- df_long %>%  filter(sample %in% valid_samples) %>%
  group_by(sample) %>%
  mutate(    TotalSampleReads = sum(ReadCount),
    RelativeAbundance = (ReadCount / TotalSampleReads) * 100
  ) %>%  ungroup()

# 4. Apply Relative Abundance Threshold (< 1% masked to 0)
THRESHOLD_PCT <- 1

df_normalized <- df_filtered %>%
  mutate( FilteredAbundance = ifelse(RelativeAbundance >= THRESHOLD_PCT, RelativeAbundance, 0) )

# Drop alleles that have 2 across all remaining valid samples
active_alleles <- df_normalized %>%  group_by(Allele_Clean) %>%
  summarize(MaxAbundance = max(FilteredAbundance)) %>%
  filter(MaxAbundance > 2) %>%  pull(Allele_Clean)

df_plot <- df_normalized %>%  filter(Allele_Clean %in% active_alleles)
df_plot$sample <- factor(df_plot$sample)
df_plot$sample <- sub("_S\\d+_L\\d+", "", df_plot$sample)
# ==============================================================================
# Plot 1: Normalized Relative Abundance Stacked Bar Plot
# ==============================================================================
p_filtered_bar <- ggplot(df_plot, aes(x = sample, y = FilteredAbundance, fill = Allele_Clean)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_viridis_d(option = "plasma", name = "SLA allele") +
  theme_classic(base_size = 24) +
  theme(    axis.text.x = element_text(angle =90, hjust = 1, vjust = 1, 
    face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
legend.text = element_text(size = 20),
legend.key.size = unit(1.2, "cm") ) +
  labs(    x = "",    y = "Relative abundance (%)"  )
ggsave("../Normalised_SLA_Composition_Barplot.png", plot = p_filtered_bar,
   width =16, height =6.5, dpi = 300)

# ==============================================================================
# Plot 2: Normalized & Filtered Heatmap
# ==============================================================================
p_filtered_heatmap <- ggplot(df_plot, aes(x = Allele_Clean, y = sample, fill = FilteredAbundance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(    option = "magma",
    name = "Abundance (%)",    na.value = "grey95"  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle =90, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 19),
legend.key.size = unit(1.4, "cm")   ) +
  labs(    x = "",    y = ""  )
ggsave("../Normalised_SLA_Heatmap.png", plot = p_filtered_heatmap, width =6, height =11, dpi = 300)



df_pig_plot <- df_plot %>%   filter(grepl("Pig", sample)) %>% droplevels()
df_pig_plot$Allele
df_pig_plot <- subset(df_pig_plot,  Allele!="IPD-MHC:SLA09764_SLA-9*03:01_1086_bp" )
p_pig_filtered_bar <- ggplot(df_pig_plot, aes(x = sample, y = FilteredAbundance, fill = Allele_Clean)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_viridis_d(option = "plasma", name = "SLA allele") +
  theme_classic(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 19),
    legend.key.size = unit(1.5, "cm") ) +
  labs(x = "", y = "Relative abundance (%)")
ggsave("../Normalised_SLA_Composition_Barplot_Pigs.png", 
       plot = p_pig_filtered_bar, width = 10, height =7.5, dpi = 300)

# ==============================================================================
# Plot 2: Normalized & Filtered Heatmap (Pigs Only)
p_pig_filtered_heatmap <- ggplot(df_pig_plot, aes(x = Allele_Clean, y = sample, fill = FilteredAbundance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", name = "Abundance (%)", na.value = "grey95") +
  theme_minimal(base_size = 19) +
  theme(   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 18),
legend.key.size = unit(.72, "cm")  ) +
  labs(x = "", y = "")
ggsave("../Normalised_SLA_Heatmap_Pigs.png", 
       plot = p_pig_filtered_heatmap, width =6, height = 5.5, dpi = 300)

# ==============================================================================
# ==============================================================================

df_st_plot <- df_plot %>%  filter(grepl("ST", sample)) %>% droplevels()

# Step 1: Filter ST samples and strip the _SXX_L001 tag
df_st <- df_plot %>%  filter(grepl("ST", sample)) %>%
  mutate(sample = sub("_S\\d+_L\\d+", "", sample))

# Step 2: Order SAMPLES (x-axis) by abundance of the dominant allele (e.g., 14:02)
# Using grepl to prevent exact string matching failures
sample_ranks <- df_st %>% filter(grepl("14:02", Allele_Clean)) %>%
  arrange(desc(FilteredAbundance)) %>% pull(sample)

# Guarantee no samples are accidentally dropped if they lack 14:02
all_st_samples <- unique(df_st$sample)
final_sample_order <- c(sample_ranks, setdiff(all_st_samples, sample_ranks))

# Step 3: Order ALLELES (legend & stack colors) by total prevalence across all ST samples
allele_ranks <- df_st %>% group_by(Allele_Clean) %>%
  summarise(total_abundance = sum(FilteredAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>% pull(Allele_Clean)

# Step 4: Apply explicit factor levels
df_st_plot <- df_st %>%
  mutate( sample = factor(sample, levels = final_sample_order),
    Allele_Clean = factor(Allele_Clean, levels = allele_ranks)
  ) %>% droplevels()


# Summary table containing both detection prevalence and average abundance
a1 <- (df_st_plot %>%
  group_by(Allele_Clean) %>%
  summarise(total_abundance = sum(FilteredAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)))
  print(a1)
  a1$Allele_Clean[7:9]
df_st_plot2 <- subset(df_st_plot, Allele!="IPD-MHC:SLA06202_SLA-3*05:02_1086_bp" 
& Allele!="IPD-MHC:SLA09709_SLA-1*20:02_546_bp"
& Allele!="IPD-MHC:SLA08447_SLA-1*08:11_1086_bp" ) 

unique(df_st_plot2$Allele)

# Step 5: Plot
p_st_filtered_bar <- ggplot(df_st_plot2, aes(x = sample, y = FilteredAbundance, fill = Allele_Clean)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_viridis_d(option = "plasma", name = "SLA allele") +
  theme_classic(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 19),
    legend.key.size = unit(1.6, "cm") ) +
  labs(x = "", y = "Relative abundance (%)")
ggsave("../Normalised_SLA_Composition_Barplot_ST.png", 
       plot = p_st_filtered_bar, width = 15, height = 7, dpi = 300)

# Plot 2: Normalized & Filtered Heatmap (Pigs Only)
# ==============================================================================
p_st_filtered_heatmap <- ggplot(df_st_plot2, aes(x = Allele_Clean, y = sample, fill = FilteredAbundance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", name = "Abundance (%)", na.value = "grey95") +
  theme_minimal(base_size = 19) +
  theme(   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank(),
    legend.text = element_text(size = 18),
legend.key.size = unit(1.8, "cm")  ) +
  labs(x = "", y = "")
ggsave("../Normalised_SLA_Heatmap_ST.png", 
       plot = p_st_filtered_heatmap, width =7, height =11, dpi = 300)
