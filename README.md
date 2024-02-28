Description:
This R script conducts an extensive eQTL analysis for Alzheimer's disease, utilizing data from the GWASCatalog and BRAINEAC database. 
The analysis includes the identification of significant SNPs associated with Alzheimer's disease, their mapping to genes, distance to closest genes, and downstream/upstream categorization. 
Subsequently, visualizations are generated, including bar plots for SNP distribution relative to genes, chromosome distribution, and a histogram illustrating the distance from significant SNPs to the closest genes.
Additionally, the script explores potential eQTL relationships by checking whether SNPs located upstream of genes are eQTLs in the BRAINEAC database.

The script generates the following output files:
Alzheimers_GWAS_Significant_SNPs.xlsx: Excel file containing significant SNPs, mapped genes, distances, and p-values.
SNP_Location_BarPlot.png: Bar plot illustrating the distribution of SNPs relative to genes.
Chromosome_Distribution_BarPlot.png: Bar plot displaying the chromosome distribution of significant SNPs.
Distance_to_Genes_Histogram.png: Histogram showing the distance from significant SNPs to the closest genes.
eQTL_Analysis_Results.xlsx: Table with eQTL information for SNPs located upstream of genes.
