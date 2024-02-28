
# Importing data

library(dplyr)
data_associations = read.delim("gwas-association-downloaded_2024-02-27-EFO_1001870-withChildTraits.tsv")


############# Significant SNPs (p_value < 5E-8

# Deleting first row 
columns_deleted = c(1,2,3,4,5,6,7,8,9,10,32,33,34,35,36,37,38)
data_associations = data_associations[, -columns_deleted]


# Selecting p values that are significant
significant = data_associations[data_associations$P.VALUE < 5e-8, ]
significant

columns_deleted_2 = c(1,2,3,4,8,11,13,14,15,16,17,19,20,21,22)
significant_snps = significant[, -columns_deleted_2]


# Counting number of snps that are upstream, downstream and within genes 
filtered_data = significant_snps %>%
  mutate(
    DISTANCE = pmin(ifelse(is.na(UPSTREAM_GENE_DISTANCE), 0, abs(UPSTREAM_GENE_DISTANCE)),
                    ifelse(is.na(DOWNSTREAM_GENE_DISTANCE), 0, DOWNSTREAM_GENE_DISTANCE)),
    LOCATION = ifelse(DISTANCE == 0, 0,
                      ifelse(DISTANCE == DOWNSTREAM_GENE_DISTANCE, 2, 
                             ifelse(DISTANCE == UPSTREAM_GENE_DISTANCE, 1, NA)))
  )%>%
  select(MAPPED_GENE, DISTANCE, SNPS, P.VALUE, UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE, LOCATION)

filtered_data

sig_filtered_data = filtered_data[, -c(5,6,7)]
write.csv(sig_filtered_data, "C:/Users/Admin/Desktop/significant_genes.csv")


################ PLOTTING FOR UPSTREAM, DOWNSTREAM, WITHIN GENES

counts = filtered_data %>%
  group_by(LOCATION) %>%
  summarise(Count = n())

print(counts)

counts = counts %>%
  mutate(LOCATION = case_when(
    LOCATION == 0 ~ "Within gene",
    LOCATION == 1 ~ "Upstream",
    LOCATION == 2 ~ "Downstream",
  ))

# bar plot 
library(ggplot2)
library(data.table)

ggplot(counts, aes(x = LOCATION, y = Count, fill = LOCATION)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Counts of SNPs Distance",
       x = "Location",
       y = "Counts of SNP") + 
  theme_minimal() 
  

##################### CHROMOSOME COUNT WITH SNPS

chromosome = significant[, c(2, 12)]

# Subsetting rows based on chr values 
subsetted_rows = chromosome[apply(chromosome, 1, function(row) any(grepl("chr", row))), ]

# Separating the column values by setting : delimiter 
subsetted_rows = separate(subsetted_rows, SNPS, into = c("sep_chr", "sep_snp"), sep = ":", remove = FALSE)

# Removing chr string 
subsetted_rows$sep_chr = sub("chr", "", subsetted_rows$sep_chr)

# Removing first two columns 
subsetted_rows = subsetted_rows[, -c(1, 2)]

# Changing column names 
new_column_names = c("CHR_ID", "SNPS")
colnames(subsetted_rows) = new_column_names

# Removing empty rows 
chromosome = chromosome %>%
  mutate(CHR_ID = na_if(CHR_ID, ""))
chromosome

chromosome = na.omit(chromosome)

# Merging two dataframes to create new dataframe 
chromosome = rbind(chromosome, subsetted_rows[, 1:2])

# Creating empty list 
snps_chromosome = list()

# Loop all over 
for (i in 1:22){
  snps_chromosome[[paste0("Chromosome_", i)]] = chromosome %>%
    filter(`CHR_ID` == i) %>%
    select(`SNPS`)
  
}
print(snps_chromosome)

# Initialize an empty data frame
df <- data.frame(CHR_ID = integer(), Count = integer())

# Loop over the list
for (i in 1:length(snps_chromosome)) {
  # Get the chromosome ID
  chr_id <- paste0("Chromosome_", i)
  
  # Get the count of SNPs for this chromosome
  count <- nrow(snps_chromosome[[chr_id]])
  
  
  # Add this to the data frame
  df <- rbind(df, data.frame(CHR_ID = chr_id, Count = count))
}

print(df)

# Create a vector of colors
colors <- rainbow(length(unique(df$CHR_ID)))

#Plot the bar chart
question_3_plot = ggplot(df, aes(x = CHR_ID, y = Count, fill = CHR_ID)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "CHR_ID", y = "Count", title = "Chromosome and SNP")

question_3_plot


############### HISTOGRAM FOR SNPS AND DISTANCES

distance_filtered_data = subset(filtered_data, DISTANCE != 0)

# Plot the histogram
question_4_plot = ggplot(distance_filtered_data, aes(x = DISTANCE)) +
  geom_histogram(binwidth = 500, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(x = "Distance to Closest Gene", y = "Count", title = "Histogram of Distances from Significant SNPs to Closest Genes")

question_4_plot


############ BRAINEAC UPSTREAM 

filtered_data_upstream = significant

filtered_data_upstream <- significant %>%
  mutate(Selected = if_else(UPSTREAM_GENE_DISTANCE < DOWNSTREAM_GENE_DISTANCE, UPSTREAM_GENE_DISTANCE, DOWNSTREAM_GENE_DISTANCE),
         Label = if_else(UPSTREAM_GENE_DISTANCE < DOWNSTREAM_GENE_DISTANCE, "UPSTREAM", "DOWNSTREAM"))

filtered_data_upstream <- filtered_data_upstream %>%
  filter(Label == "UPSTREAM")

upstream_rsids = filtered_data_upstream$SNPS
write.csv(upstream_rsids, "C:/Users/Admin/Desktop/rsids.csv")
