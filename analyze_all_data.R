#R script to perform misc. analyses on Mash and pham data for Mavrich & Hatfull, Nature Micro, 2017
#Travis Mavrich
#Note: this code merges and analyzes various datasets and exports datasets for downstream analyses in Python and Excel
#Distinct analyses and code blocks separated by "###"





###Define functions

#The standard genomic similarity plot parameters
plot_genomic_similarity_standard <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$modified_mash_distance,table$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
  abline(0,2,lty=2,lwd=3,col="grey")
  
}
  
  
  


#Genomic similarities color-coded by type of intra- and inter-cluster comparisons
plot_cluster_specific_profiles <- function(table,cluster){
  
  table$cluster_specific_one_or_two <- ifelse(table$ref_phage_cluster == cluster | table$query_phage_cluster == cluster,TRUE,FALSE)
  table$cluster_specific_both <- ifelse(table$ref_phage_cluster == cluster & table$query_phage_cluster == cluster,TRUE,FALSE)
  table$cluster_specific_one <- ifelse(table$cluster_specific_one_or_two == TRUE & table$cluster_specific_both == FALSE,TRUE,FALSE)
  cluster_specific_one <- subset(table,table$cluster_specific_one == TRUE)
  cluster_specific_both <- subset(table,table$cluster_specific_both == TRUE)
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific_both$modified_mash_distance,cluster_specific_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
  par(new=TRUE)
  plot(cluster_specific_one$modified_mash_distance,cluster_specific_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
  abline(0,2,lty=2,lwd=3,col="grey")
  return(nrow(cluster_specific_one) + nrow(cluster_specific_both))
  
}


#Genomic similarities colored by evolutionary mode
plot_genomic_similarity_tricolored <- function(table1,table2,table3){
  
  par(mar=c(4,8,4,4))
  plot(table1$modified_mash_distance,table1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
  par(new=TRUE)
  plot(table2$modified_mash_distance,table2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
  par(new=TRUE)
  plot(table3$modified_mash_distance,table3$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
  
}
  


#Genomic similarities colored orange
plot_genomic_similarity_orange <- function(table){
  
  par(mar=c(4,8,4,4))
  plot(table$modified_mash_distance,table$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
  abline(0,2,lty=2,lwd=3,col="grey")
  
}
  

#Function to compute how many comparisons in each data subset are assigned to each sector
compute_sector_distribution <- function(table){
  
  sector_intra_subcluster_count <- sum(table$sector_intra_subcluster == TRUE)
  sector_intra_cluster_count <- sum(table$sector_intra_cluster == TRUE)
  sector_inter_cluster_distant_homology_count <- sum(table$sector_inter_cluster_distant_homology == TRUE)
  sector_inter_cluster_hgt_count <- sum(table$sector_inter_cluster_hgt == TRUE)
  sector_no_similarity_count <- sum(table$sector_no_similarity == TRUE)
  sector_non_cds_similarity_count <- sum(table$sector_non_cds_similarity == TRUE)
  sector_filtered_out_count <- sum(table$sector_filtered_out == TRUE)
  
  
  print(c("sector_intra_subcluster_count",
          "sector_intra_cluster_count",
          "sector_inter_cluster_distant_homology_count",
          "sector_inter_cluster_hgt_count",
          "sector_no_similarity_count",
          "sector_non_cds_similarity_count",
          "sector_filtered_out_count"))
  
  print(c(sector_intra_subcluster_count,
          sector_intra_cluster_count,
          sector_inter_cluster_distant_homology_count,
          sector_inter_cluster_hgt_count,
          sector_no_similarity_count,
          sector_non_cds_similarity_count,
          sector_filtered_out_count))
  
  print(c("sector_intra_subcluster_percent",
          "sector_intra_cluster_percent",
          "sector_inter_cluster_distant_homology_percent",
          "sector_inter_cluster_hgt_percent",
          "sector_no_similarity_percent",
          "sector_non_cds_similarity_percent",
          "sector_filtered_out_percent"))
  
  print(c(sector_intra_subcluster_count,
          sector_intra_cluster_count,
          sector_inter_cluster_distant_homology_count,
          sector_inter_cluster_hgt_count,
          sector_no_similarity_count,
          sector_non_cds_similarity_count,
          sector_filtered_out_count)/nrow(table))
  
}









###Import primary mash and gene content dissimilarity data

#Import mash dataset
#Format:
#0 = reference genome
#1 = query genome
#2 = mash distance 
#3 = mash p-value
#4 = mash kmer count
#5 = reference_query
mash_table <- read.csv("processed_mash_output.csv",sep=",",header=TRUE)
names(mash_table) <- c("mash_reference",
                       "mash_query",
                       "mash_distance",
                       "mash_pvalue",
                       "mash_count",
                       "mash_ref_query")



#Import host taxonomy data and phage metadata
#Format:
#0 = phage_identifier
#1 = header_source_info
#2 = database
#3 = host_superkingdom
#4 = host_phylum
#5 = host_class
#6 = host_order
#7 = host_family
#8 = host_genus
#9 = phage_superkingdom
#10 = phage_viral_type
#11 = phage_order
#12 = phage_family
#13 = phage_genus
#14 = phage_cluster
#15 = phage_subcluster
#16 = phage_cluster_source
#17 = phage_temperate
#18 = phage_size
#19 = phage_gene_count
#20 = phage_toxic
#21 = phage_predicted_temperate
phage_metadata_table <- read.csv("phage_host_data.csv",sep=",",header=TRUE)


#Convert all Unspecified fields to NA missing value
phage_metadata_table[phage_metadata_table == "Unspecified"] <- NA



#Modify column names of host data and merge with mash table
#Data for all phages in processed mash table should be present in phage metadata table.
#As a result, do not select all.x=TRUE option. This way, any missing rows indicates an error in the table.


#Match metadata for both the query and reference phages in each pairwise comparison
names(phage_metadata_table) <- c("phage_identifier","query_header_source_info","query_database",
                             "query_host_superkingdom","query_host_phylum","query_host_class","query_host_order","query_host_family","query_host_genus",
                             "query_phage_superkingdom","query_phage_viral_type","query_phage_order","query_phage_family","query_phage_genus",
                             "query_phage_cluster","query_phage_subcluster","query_phage_cluster_source",
                             "query_phage_temperate","query_size","query_gene_count",
                             "query_toxic","query_predicted_temperate")


main_data_table <- merge(mash_table,phage_metadata_table,by.x="mash_query",by.y="phage_identifier")

names(phage_metadata_table) <- c("phage_identifier","ref_header_source_info","ref_database",
                       "ref_host_superkingdom","ref_host_phylum","ref_host_class","ref_host_order","ref_host_family","ref_host_genus",
                       "ref_phage_superkingdom","ref_phage_viral_type","ref_phage_order","ref_phage_family","ref_phage_genus",
                       "ref_phage_cluster","ref_phage_subcluster","ref_phage_cluster_source",
                       "ref_phage_temperate","ref_size","ref_gene_count",
                       "ref_toxic","ref_predicted_temperate")



main_data_table <- merge(main_data_table,phage_metadata_table,by.x="mash_reference",by.y="phage_identifier")

names(phage_metadata_table) <- c("phage_identifier","header_source_info","database",
                       "host_superkingdom","host_phylum","host_class","host_order","host_family","host_genus",
                       "phage_superkingdom","phage_viral_type","phage_order","phage_family","phage_genus",
                       "phage_cluster","phage_subcluster","phage_cluster_source",
                       "phage_temperate","size","gene_count",
                       "toxic","predicted_temperate")



#Compare genome sizes
main_data_table$size_diff <- abs(main_data_table$ref_size - main_data_table$query_size)
main_data_table$size_diff_ref_percent <- main_data_table$size_diff / main_data_table$ref_size
main_data_table$size_diff_query_percent <- main_data_table$size_diff / main_data_table$query_size
main_data_table$size_diff_min_percent <- ifelse(main_data_table$size_diff_ref_percent < main_data_table$size_diff_query_percent,main_data_table$size_diff_ref_percent,main_data_table$size_diff_query_percent)
main_data_table$size_diff_max_percent <- ifelse(main_data_table$size_diff_ref_percent > main_data_table$size_diff_query_percent,main_data_table$size_diff_ref_percent,main_data_table$size_diff_query_percent)
main_data_table$size_diff_ave_percent <- (main_data_table$size_diff_ref_percent + main_data_table$size_diff_query_percent)/2






#Import ANI data to check Mash vs Pham of specific clusters
#This ANI data is different than the ANI data from 79 genomes used to optimize the Mash parameters
#Format:
#0 = ani_ref_query
#1 = ani_reference_genome
#2 = ani_query_genome
#3 = ani_ani
#4 = ani_distance
ani_data <- read.csv("ani_cluster_data.csv",sep=",",header=TRUE)
ani_data$ani_ref_query <- as.character(ani_data$ani_ref_query)
ani_data$ani_distance <- 1 - ani_data$ani_ani

#Merge with the main table
main_data_table <- merge(main_data_table,ani_data,by.x="mash_ref_query",by.y="ani_ref_query",all.x=TRUE)



#Import pham data and merge with mash table
#Format
#0 = phage1
#1 = phage1_number_unshared_phams
#2 = phage1_shared_proportion
#3 = phage2
#4 = phage2_number_unshared_phams
#5 = phage2_shared_proportion
#6 = number_shared_phams
#7 = average_shared_proportion
#8 = jaccard_similarity
#9 = shared_pham_distribution_mean
#10 = shared_pham_distribution_median
#11 = shared_pham_distribution_max
#12 = unshared_pham_distribution_mean
#13 = unshared_pham_distribution_median
#14 = unshared_pham_distribution_max
#15 = unshared_orpham_count
pham_table <- read.csv("pairwise_pham_proportions.csv",sep=",",header=TRUE)

#Compute gene content dissimilarity
pham_table$pham_dissimilarity <- 1 - pham_table$average_shared_proportion
pham_table$jaccard_dissimilarity <- 1 - pham_table$jaccard_similarity

names(pham_table) <- c("pham_phage1","pham_phage1_number_unshared_phams","pham_phage1_shared_proportion",
                                "pham_phage2","pham_phage2_number_unshared_phams","pham_phage2_shared_proportion",
                                "pham_number_shared_phams","pham_average_shared_proportion",
                                "pham_jaccard_similarity","pham_shared_pham_distribution_mean",
                                "pham_shared_pham_distribution_median","pham_shared_pham_distribution_max",
                                "pham_unshared_pham_distribution_mean","pham_unshared_pham_distribution_median",
                                "pham_unshared_pham_distribution_max","pham_unshared_orpham_count",
                                "pham_pham_dissimilarity","pham_jaccard_dissimilarity")


#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
pham_table$pham_phage1_phage2 <- paste(pham_table$pham_phage1,"_",pham_table$pham_phage2,sep="")
pham_table$pham_phage1_phage2 <- as.factor(pham_table$pham_phage1_phage2)

#To retain all rows, be sure to keep all.x=TRUE.
#But all rows don't need to be retained = making scatter plots or histograms can cause errors if some rows are missing data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
main_data_table <- merge(main_data_table,pham_table,by.x="mash_ref_query",by.y="pham_phage1_phage2")



#Assign filter status and change mash distance if data is not significant
main_data_table$filter <- ifelse(main_data_table$mash_pvalue < 1e-10 & main_data_table$size_diff_max_percent < 1,TRUE,FALSE)

#Alternatively, the max percent parameter can be omitted with minimal effects on final analysis.
#main_data_table$filter <- ifelse(main_data_table$mash_pvalue < 1e-10,TRUE,FALSE)

#At this point, the max mash distance of all filtered comparisons < 0.5. So set the distance of all comparisons that did not pass the filter = 0.5
main_data_table$modified_mash_distance <- ifelse(main_data_table$filter == TRUE,main_data_table$mash_distance,0.5)





#Assign evolutionary mode to each pairwise comparison.
#This assigns all data in the dataset, but evolutionary mode in this analysis is only relevant to dsDNA
main_data_table$gene_flux_part1 <- ifelse(main_data_table$modified_mash_distance < 0.1666667 & main_data_table$pham_pham_dissimilarity > (main_data_table$modified_mash_distance * 3.5),TRUE,FALSE)
main_data_table$gene_flux_part2 <- ifelse(main_data_table$modified_mash_distance > 0.1666667 & main_data_table$pham_pham_dissimilarity > (main_data_table$modified_mash_distance * 2 + 0.25),TRUE,FALSE)
main_data_table$gene_flux_category <- ifelse(main_data_table$gene_flux_part1 == TRUE | main_data_table$gene_flux_part2 ==TRUE,"high","low")
main_data_table$gene_flux_category <- as.factor(main_data_table$gene_flux_category)




#Compute % of comparisons that are positioned in each sector
#Assign each comparison to a sector:
#intra-subcluster
#intra-cluster
#inter-cluster with distant homology
#inter-cluster with horizontal gene transfer
#no similarity
#inter-cluster with non-CDS similarity
#data filtered out based on filtering parameters (p-value)

#Set up boundaries to define each sector
#The boundaries below represent subcluster or cluster inclusion boundaries
intra_subcluster_nucleotide_threshold_upper <- 0.20
intra_subcluster_gene_content_threshold_upper <- 0.62
intra_cluster_nucleotide_threshold_upper <- 0.42
intra_cluster_gene_content_threshold_upper <- 0.89

#Assign each comparison to each sector as TRUE or FALSE
#The last sector, "unsectored" is used to catch any errors
main_data_table$sector_intra_subcluster <- ifelse(main_data_table$modified_mash_distance < intra_subcluster_nucleotide_threshold_upper &
                                                main_data_table$modified_mash_distance >= 0 &
                                                main_data_table$pham_pham_dissimilarity < intra_subcluster_gene_content_threshold_upper &
                                                main_data_table$pham_pham_dissimilarity >= 0,TRUE,FALSE)

main_data_table$sector_intra_cluster <- ifelse(main_data_table$modified_mash_distance < intra_cluster_nucleotide_threshold_upper &
                                             main_data_table$pham_pham_dissimilarity < intra_cluster_gene_content_threshold_upper &
                                             (main_data_table$modified_mash_distance >= intra_subcluster_nucleotide_threshold_upper |
                                                (main_data_table$modified_mash_distance < intra_subcluster_nucleotide_threshold_upper &
                                                   main_data_table$pham_pham_dissimilarity >= intra_subcluster_gene_content_threshold_upper)),TRUE,FALSE)

main_data_table$sector_inter_cluster_distant_homology <- ifelse(main_data_table$modified_mash_distance < 0.5 &
                                                              main_data_table$modified_mash_distance >= intra_cluster_nucleotide_threshold_upper &
                                                              main_data_table$pham_pham_dissimilarity < intra_cluster_gene_content_threshold_upper &
                                                              main_data_table$pham_pham_dissimilarity >= 0,TRUE,FALSE)

main_data_table$sector_inter_cluster_hgt <- ifelse(main_data_table$modified_mash_distance < intra_cluster_nucleotide_threshold_upper &
                                                 main_data_table$modified_mash_distance >= 0 &
                                                 main_data_table$pham_pham_dissimilarity < 1 &
                                                 main_data_table$pham_pham_dissimilarity >= intra_cluster_gene_content_threshold_upper,TRUE,FALSE)

main_data_table$sector_no_similarity <- ifelse(main_data_table$modified_mash_distance < 0.5 &
                                             main_data_table$modified_mash_distance >= intra_cluster_nucleotide_threshold_upper &
                                             main_data_table$pham_pham_dissimilarity >= intra_cluster_gene_content_threshold_upper,TRUE,FALSE)

main_data_table$sector_non_cds_similarity <- ifelse(main_data_table$modified_mash_distance < intra_cluster_nucleotide_threshold_upper &
                                                  main_data_table$modified_mash_distance >= 0 &
                                                  main_data_table$pham_pham_dissimilarity == 1,TRUE,FALSE)

main_data_table$sector_filtered_out <- ifelse(main_data_table$modified_mash_distance == 0.5,TRUE,FALSE)

main_data_table$sector_unsectored <- ifelse(main_data_table$sector_intra_subcluster == FALSE &
                                          main_data_table$sector_intra_cluster == FALSE &
                                          main_data_table$sector_inter_cluster_distant_homology == FALSE &
                                          main_data_table$sector_inter_cluster_hgt == FALSE &
                                          main_data_table$sector_no_similarity == FALSE &
                                          main_data_table$sector_non_cds_similarity == FALSE &
                                          main_data_table$sector_filtered_out == FALSE,TRUE,FALSE)








#Compare ref and query host and phage metadata columns
main_data_table$host_superkingdom_compare <- ifelse(main_data_table$ref_host_superkingdom==main_data_table$query_host_superkingdom,as.character(main_data_table$ref_host_superkingdom),"different")
main_data_table$host_phylum_compare <- ifelse(main_data_table$ref_host_phylum==main_data_table$query_host_phylum,as.character(main_data_table$ref_host_phylum),"different")
main_data_table$host_class_compare <- ifelse(main_data_table$ref_host_class==main_data_table$query_host_class,as.character(main_data_table$ref_host_class),"different")
main_data_table$host_order_compare <- ifelse(main_data_table$ref_host_order==main_data_table$query_host_order,as.character(main_data_table$ref_host_order),"different")
main_data_table$host_family_compare <- ifelse(main_data_table$ref_host_family==main_data_table$query_host_family,as.character(main_data_table$ref_host_family),"different")
main_data_table$host_genus_compare <- ifelse(main_data_table$ref_host_genus==main_data_table$query_host_genus,as.character(main_data_table$ref_host_genus),"different")
main_data_table$phage_superkingdom_compare <- ifelse(main_data_table$ref_phage_superkingdom==main_data_table$query_phage_superkingdom,as.character(main_data_table$ref_phage_superkingdom),"different")
main_data_table$phage_viral_type_compare <- ifelse(main_data_table$ref_phage_viral_type==main_data_table$query_phage_viral_type,as.character(main_data_table$ref_phage_viral_type),"different")
main_data_table$phage_order_compare <- ifelse(main_data_table$ref_phage_order==main_data_table$query_phage_order,as.character(main_data_table$ref_phage_order),"different")
main_data_table$phage_family_compare <- ifelse(main_data_table$ref_phage_family==main_data_table$query_phage_family,as.character(main_data_table$ref_phage_family),"different")
main_data_table$phage_genus_compare <- ifelse(main_data_table$ref_phage_genus==main_data_table$query_phage_genus,as.character(main_data_table$ref_phage_genus),"different")
main_data_table$phage_cluster_compare <- ifelse(main_data_table$ref_phage_cluster==main_data_table$query_phage_cluster,as.character(main_data_table$ref_phage_cluster),"different")
main_data_table$phage_subcluster_compare <- ifelse(main_data_table$ref_phage_subcluster==main_data_table$query_phage_subcluster,as.character(main_data_table$ref_phage_subcluster),"different")
main_data_table$phage_cluster_source_compare <- ifelse(main_data_table$ref_phage_cluster_source==main_data_table$query_phage_cluster_source,as.character(main_data_table$ref_phage_cluster_source),"different")
main_data_table$phage_temperate_compare <- ifelse(main_data_table$ref_phage_temperate==main_data_table$query_phage_temperate,as.character(main_data_table$ref_phage_temperate),"different")
main_data_table$phage_toxic_compare <- ifelse(main_data_table$ref_toxic==main_data_table$query_toxic,as.character(main_data_table$ref_toxic),"different")
main_data_table$phage_predicted_temperate_compare <- ifelse(main_data_table$ref_predicted_temperate==main_data_table$query_predicted_temperate,as.character(main_data_table$ref_predicted_temperate),"different")





#Now set all new columns to factor class
main_data_table$host_superkingdom_compare <- as.factor(main_data_table$host_superkingdom_compare)
main_data_table$host_phylum_compare <- as.factor(main_data_table$host_phylum_compare)
main_data_table$host_class_compare <- as.factor(main_data_table$host_class_compare)
main_data_table$host_order_compare <- as.factor(main_data_table$host_order_compare)
main_data_table$host_family_compare <- as.factor(main_data_table$host_family_compare)
main_data_table$host_genus_compare <- as.factor(main_data_table$host_genus_compare)
main_data_table$phage_superkingdom_compare <- as.factor(main_data_table$phage_superkingdom_compare)
main_data_table$phage_viral_type_compare <- as.factor(main_data_table$phage_viral_type_compare)
main_data_table$phage_order_compare <- as.factor(main_data_table$phage_order_compare)
main_data_table$phage_family_compare <- as.factor(main_data_table$phage_family_compare)
main_data_table$phage_genus_compare <- as.factor(main_data_table$phage_genus_compare)
main_data_table$phage_cluster_compare <- as.factor(main_data_table$phage_cluster_compare)
main_data_table$phage_subcluster_compare <- as.factor(main_data_table$phage_subcluster_compare)
main_data_table$phage_cluster_source_compare <- as.factor(main_data_table$phage_cluster_source_compare)
main_data_table$phage_temperate_compare <- as.factor(main_data_table$phage_temperate_compare)
main_data_table$phage_toxic_compare <- as.factor(main_data_table$phage_toxic_compare)
main_data_table$phage_predicted_temperate_compare <- as.factor(main_data_table$phage_predicted_temperate_compare)


#At this point, the primary data for the analysis has been imported and processed.
#From here on out, commands are mostly analysis and data visualization.
#Although, for some specific analyses, additional datasets are either exported for downstream analysis with Python or Excel or imported
#to be paired up with the primary data.











###All dsDNA phages
bacteria_dsDNA <- subset(main_data_table,main_data_table$host_superkingdom_compare == 'Bacteria' & main_data_table$phage_viral_type_compare == 'dsDNA')



#Fig. 1a
plot_genomic_similarity_standard(bacteria_dsDNA)

#Fig. 1a
par(mar=c(4,8,15,4))
hist(bacteria_dsDNA$modified_mash_distance,breaks=((range(bacteria_dsDNA$modified_mash_distance)[2]-range(bacteria_dsDNA$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,3e4),col="black",cex.axis=2)

#Fig. 1a
par(mar=c(4,4,15,4))
hist(bacteria_dsDNA$pham_pham_dissimilarity,breaks=((range(bacteria_dsDNA$pham_pham_dissimilarity)[2]-range(bacteria_dsDNA$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)











###Empirical lifestyle analysis
bacteria_dsDNA <- subset(main_data_table,main_data_table$host_superkingdom_compare == 'Bacteria' & main_data_table$phage_viral_type_compare == 'dsDNA')
all_temperate <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == "yes")
all_lytic <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == "no")
all_different <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == "different")


#Fig. 1c
plot_genomic_similarity_standard(all_temperate)

#Fig. 1c
plot_genomic_similarity_standard(all_lytic)

#Fig. 1c
plot_genomic_similarity_standard(all_different)










###Cluster-specific analysis
cluster_actino <- subset(main_data_table,main_data_table$phage_cluster_source_compare == "actino")

#Fig. 2a
dev.off()
plot_cluster_specific_profiles(cluster_actino,"F")

#Fig. 2a
dev.off()
plot_cluster_specific_profiles(cluster_actino,"K")

#Fig. 2a
dev.off()
plot_cluster_specific_profiles(cluster_actino,"AO")

#Fig. 2a
dev.off()
plot_cluster_specific_profiles(cluster_actino,"B")

#Fig. 2a
dev.off()
plot_cluster_specific_profiles(cluster_actino,"BU")

#Fig. 2a
dev.off()
plot_cluster_specific_profiles(cluster_actino,"BD")

#Supp. Fig. 11c
dev.off()
plot_cluster_specific_profiles(cluster_actino,"N")







###Cluster A analysis
#Check A1 and non-A1 distances to other non-A phages.
#Since only actino phages that have been clustered are used, all rows should have cluster data.
#However, not all rows necessarily have subcluster data, so this must be taken into account.
cluster_actino <- subset(main_data_table,main_data_table$phage_cluster_source_compare == "actino")



cluster_actino$subcluster_A1_one_or_two <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_neither <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,TRUE,FALSE)
cluster_actino$subcluster_A1_both <- ifelse((cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_one <- ifelse(cluster_actino$subcluster_A1_one_or_two == TRUE & cluster_actino$subcluster_A1_both == FALSE,TRUE,FALSE)

cluster_actino$cluster_A_one_or_two <- ifelse(cluster_actino$ref_phage_cluster == "A" | cluster_actino$query_phage_cluster == "A",TRUE,FALSE)
cluster_actino$cluster_A_both <- ifelse(cluster_actino$ref_phage_cluster == "A" & cluster_actino$query_phage_cluster == "A",TRUE,FALSE)
cluster_actino$cluster_A_both_but_notA1 <- ifelse(cluster_actino$cluster_A_both == TRUE & cluster_actino$subcluster_A1_neither == TRUE,TRUE,FALSE)
cluster_actino$cluster_A_one <- ifelse(cluster_actino$cluster_A_one_or_two == TRUE & cluster_actino$cluster_A_both == FALSE,TRUE,FALSE)
cluster_actino$cluster_A_one_but_notA1 <- ifelse(cluster_actino$cluster_A_one == TRUE & cluster_actino$subcluster_A1_neither == TRUE,TRUE,FALSE)


#Comparisons with at least one A1 phage
subcluster_A1_one_or_two_all <- subset(cluster_actino,cluster_actino$subcluster_A1_one_or_two == TRUE)

#Comparisons between all A1 phages
subcluster_A1_both <- subset(subcluster_A1_one_or_two_all,subcluster_A1_one_or_two_all$phage_subcluster_compare == "A1")

#Comparisons of one and only one A1 phage
subcluster_A1_one_all <- subset(subcluster_A1_one_or_two_all,subcluster_A1_one_or_two_all$subcluster_A1_one == TRUE)

#Comparisons of one and only one A1 phage, and the other phage is a Cluster A phage
subcluster_A1_one_clusterA <- subset(subcluster_A1_one_all,subcluster_A1_one_all$phage_cluster_compare == "A")

#Comparisons of one and only one A1 phage, and the other phage is NOT a Cluster A phage
#The phage_cluster_compare should NOT be "A", but it can be "NA" (since some phages are not clustered), or it can be "different"
subcluster_A1_one_other_not_clusterA <- subset(subcluster_A1_one_all,
                                               is.na(subcluster_A1_one_all$phage_cluster_compare) | 
                                               subcluster_A1_one_all$phage_cluster_compare == "different")

#Comparisons of one and only non-A1, and no other Cluster A phages
cluster_A_one_but_notA1 <- subset(cluster_actino,cluster_actino$cluster_A_one_but_notA1 == TRUE)

#Comparisons of two Cluster A phages, but neither are A1 phages
cluster_A_both_but_notA1 <- subset(cluster_actino,cluster_actino$cluster_A_both_but_notA1 == TRUE)


#Fig. 3a (left)
par(mar=c(4,8,4,4))
plot(cluster_A_both_but_notA1$modified_mash_distance,cluster_A_both_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(subcluster_A1_both$modified_mash_distance,subcluster_A1_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(subcluster_A1_one_clusterA$modified_mash_distance,subcluster_A1_one_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")#422
par(new=TRUE)
plot(subcluster_A1_one_other_not_clusterA$modified_mash_distance,subcluster_A1_one_other_not_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(cluster_A_one_but_notA1$modified_mash_distance,cluster_A_one_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
abline(0,2,lty=2,lwd=3,col="grey")











###Misc clusters genometrics

#Format
#0 = phage
#1 = cluster
#2 = size
#3 = total_gene_count
#4 = Unspecified_gene_count
#5 = lysis_gene_count
#6 = lysogeny_gene_count 
#7 = recombination_replication_gene_count
#8 = structure_assembly_gene_count
misc_clusters_genometrics <- read.csv("misc_clusters_genometrics.csv",sep=",",header=TRUE)
misc_clusters_genometrics$phage_cluster <- factor(misc_clusters_genometrics$phage_cluster,c("A1","F","BD","K","non-A1","B"))


#Format
#0 = phage
#1 = cluster
#2 = GC
misc_clusters_gc <- read.csv("misc_clusters_gc.csv",sep=",",header=TRUE)
misc_clusters_gc$Cluster <- factor(misc_clusters_gc$Cluster,c("A1","F","BD","K","non-A1","B"))


#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$size ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$size ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$total_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$total_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$Unspecified_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$Unspecified_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$lysis_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$lysis_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$lysogeny_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$lysogeny_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$recombination_replication_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$recombination_replication_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$structure_assembly_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$structure_assembly_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

#Fig. 3f
par(mar=c(8,24,4,4))
boxplot(misc_clusters_gc$GC ~ misc_clusters_gc$Cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_gc$GC ~ misc_clusters_gc$Cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))











###Host phylum and empirical lifestyle analysis
#Only focus on dsDNA phages, so remove other phage types.
bacteria_dsDNA <- subset(main_data_table,main_data_table$host_superkingdom_compare == 'Bacteria' & main_data_table$phage_viral_type_compare == 'dsDNA')
host_phylum_diff <- subset(bacteria_dsDNA,bacteria_dsDNA$host_phylum_compare == "different")

actino <- subset(bacteria_dsDNA,bacteria_dsDNA$host_phylum_compare == "Actinobacteria")
bacter <- subset(bacteria_dsDNA,bacteria_dsDNA$host_phylum_compare == "Bacteroidetes")
cyano <- subset(bacteria_dsDNA,bacteria_dsDNA$host_phylum_compare == "Cyanobacteria")
firm <- subset(bacteria_dsDNA,bacteria_dsDNA$host_phylum_compare == "Firmicutes")
proteo <- subset(bacteria_dsDNA,bacteria_dsDNA$host_phylum_compare == "Proteobacteria")


#Same phylum

#Fig 4a
plot_genomic_similarity_standard(actino)

#Fig 4b
plot_genomic_similarity_standard(bacter)

#Fig 4c
plot_genomic_similarity_standard(cyano)

#Fig 4d
plot_genomic_similarity_standard(firm)

#Fig 4e
plot_genomic_similarity_standard(proteo)




#Truncated histogram distributions of mash distances

#Fig 4a
par(mar=c(4,8,15,4))
hist(actino$modified_mash_distance,breaks=((range(actino$modified_mash_distance)[2]-range(actino$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,3e4),col="black",cex.axis=2)

#Fig 4b
par(mar=c(4,8,15,4))
hist(bacter$modified_mash_distance,breaks=((range(bacter$modified_mash_distance)[2]-range(bacter$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,50),col="black",cex.axis=2)

#Fig 4c
par(mar=c(4,8,15,4))
hist(cyano$modified_mash_distance,breaks=((range(cyano$modified_mash_distance)[2]-range(cyano$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,100),col="black",cex.axis=2)

#Fig 4d
par(mar=c(4,8,15,4))
hist(firm$modified_mash_distance,breaks=((range(firm$modified_mash_distance)[2]-range(firm$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)

#Fig 4e
par(mar=c(4,8,15,4))
hist(proteo$modified_mash_distance,breaks=((range(proteo$modified_mash_distance)[2]-range(proteo$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)





#Truncated histogram distribution of gene content dissimilarities

#Fig 4a
par(mar=c(4,4,15,4))
hist(actino$pham_pham_dissimilarity,breaks=((range(actino$pham_pham_dissimilarity)[2]-range(actino$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

#Fig 4b
par(mar=c(4,4,15,4))
hist(bacter$pham_pham_dissimilarity,breaks=((range(bacter$pham_pham_dissimilarity)[2]-range(bacter$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,500),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

#Fig 4c
par(mar=c(4,4,15,4))
hist(cyano$pham_pham_dissimilarity,breaks=((range(cyano$pham_pham_dissimilarity)[2]-range(cyano$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,100),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

#Fig 4d
par(mar=c(4,4,15,4))
hist(firm$pham_pham_dissimilarity,breaks=((range(firm$pham_pham_dissimilarity)[2]-range(firm$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

#Fig 4e
par(mar=c(4,4,15,4))
hist(proteo$pham_pham_dissimilarity,breaks=((range(proteo$pham_pham_dissimilarity)[2]-range(proteo$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)




#Different phyla
#Supp. Fig. 5c
plot_genomic_similarity_standard(host_phylum_diff)










###Compare general vs jaccard gene content dissimilarity
bacteria_dsDNA <- subset(main_data_table,main_data_table$host_superkingdom_compare == 'Bacteria' & main_data_table$phage_viral_type_compare == 'dsDNA')



#Assess correlation of general dissimilarity and jaccard dissimilarity
#Supp. Fig. 2a
par(mar=c(4,8,4,4))
plot(main_data_table$pham_jaccard_dissimilarity,main_data_table$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")


#Genomic similarity (Mash vs Pham) plot using jaccard dissimilarity
#Supp. Fig. 2b
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###ANI vs Pham Dissimilarity data for the 79 genomes that were used to optimize Mash parameters
#Data contains complete matrix of 79 x 79 comparisons, including self comparisons and duplicate (reciprocal) comparisons
#Format
#0 = ref_query
#1 = ref_phage_identifier
#2 = query_phage_identifier
#3 = ani_distance
ani79_data <- read.csv("ani_79_data.csv",sep=",",header=TRUE)

names(ani79_data) <- c("ani79_ref_query","ani79_ref_phage_identifier","ani79_query_phage_identifier","ani79_ani_distance")

#Merge tables
#The main mash table contains non-redundant comparisona and no self comparisons
ani79_analysis <- merge(main_data_table,ani79_data,by.x="mash_ref_query",by.y="ani79_ref_query")

#Supp. Fig. 2c
plot_genomic_similarity_standard(ani79_analysis)

#Supp. Fig. 2c
par(mar=c(4,8,4,4))
plot(ani79_analysis$ani79_ani_distance,ani79_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###Compare Mash to ANI for specific clusters
#Compare patterns from pham vs mash and pham vs ani comparisons
ani_mash <- subset(main_data_table,complete.cases(main_data_table$ani_distance))

#Colored plot by assigned evolutionary mode 
ani_mash_hgcf <- subset(ani_mash,ani_mash$gene_flux_category == "high" & ani_mash$phage_predicted_temperate_compare == "yes")
ani_mash_lgcf <- subset(ani_mash,ani_mash$gene_flux_category == "low" & ani_mash$phage_predicted_temperate_compare == "yes")
ani_mash_lytic <- subset(ani_mash,ani_mash$phage_predicted_temperate_compare == "no")


#Supp. Fig. 2d
plot_genomic_similarity_tricolored(ani_mash_hgcf,ani_mash_lgcf,ani_mash_lytic)

#Supp. Fig. 2d
par(mar=c(4,8,4,4))
plot(ani_mash_hgcf$ani_distance,ani_mash_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(ani_mash_lgcf$ani_distance,ani_mash_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(ani_mash_lytic$ani_distance,ani_mash_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")









###VOG analysis
#Format
#0 = phage1_name
#1 = phage1_accession
#2 = phage1_number_of_genes
#3 = phage1_number_of_vogs
#4 = phage1_number_of_unshared_vogs
#5 = phage1_number_of_shared_vog_genes
#6 = phage1_number_of_unshared_vog_genes
#7 = phage1_number_of_unshared_other_genes
#8 = phage1_shared_vog_gene_proportion
#9 = phage2_name
#10 = phage2_accession
#11 = phage2_number_of_genes
#12 = phage2_number_of_vogs
#13 = phage2_number_of_unshared_vogs
#14 = phage2_number_of_shared_vog_genes
#15 = phage2_number_of_unshared_vog_genes
#16 = phage2_number_of_unshared_other_genes
#17 = phage2_shared_vog_gene_proportion
#18 = number_of_shared_vogs
#19 = average_shared_vog_gene_proportion
#20 = gene_content_dissimilarity
vog_table <- read.csv("pairwise_vog_proportions.csv",sep=",",header=TRUE)

names(vog_table) <- c("vog_reference","vog_ref_accession","vog_ref_number_of_genes","vog_ref_number_of_vogs",
                      "vog_ref_number_of_unshared_vogs","vog_ref_number_of_shared_vog_genes","vog_ref_number_of_unshared_vog_genes","vog_ref_number_of_unshared_other_genes",
                      "vog_ref_shared_vog_gene_proportion","vog_query","vog_query_accession","vog_query_number_of_genes",
                      "vog_query_number_of_vogs","vog_query_number_of_unshared_vogs","vog_query_number_of_shared_vog_genes","vog_query_number_of_unshared_vog_genes",
                      "vog_query_number_of_unshared_other_genes","vog_query_shared_vog_gene_proportion",
                      "vog_number_of_shared_vogs","vog_average_shared_vog_gene_proportion","vog_gene_content_dissimilarity")

vog_table$vog_ref_query <- paste(vog_table$vog_reference,vog_table$vog_query,sep="_")
vog_table$vog_ref_query <- as.factor(vog_table$vog_ref_query)




#VOG data is based on 1877 genomes that are present in the merged2333 dataset.
#But the VOG data contains redundant data rows, where each comparison is represented twice, with the ref and query reversed
#So when merged to main_data_table, no need to keep all rows in either table - it is expected there will be fewer rows than in both tables
#Also, there are 2 genomes (vb_paem_c1-14-ab28__NC_026600 and pv94__NC_027368) that contain no annotated genes. These are not in the pham data,
#so even though they are present in the VOG data, I am unable to compare these two genomes. This results in 1875 genomes.
main_data_table_vog <- merge(main_data_table,vog_table,by.x="mash_ref_query",by.y="vog_ref_query")
main_data_table_vog$mash_reference <- factor(main_data_table_vog$mash_reference)
main_data_table_vog$mash_query <- factor(main_data_table_vog$mash_query)

vog_bacteria_dsDNA <- subset(main_data_table_vog,main_data_table_vog$host_superkingdom_compare == 'Bacteria' & main_data_table_vog$phage_viral_type_compare == 'dsDNA')
vog_temperate_both <- subset(vog_bacteria_dsDNA,vog_bacteria_dsDNA$phage_temperate_compare == 'yes')
vog_temperate_neither <- subset(vog_bacteria_dsDNA,vog_bacteria_dsDNA$phage_temperate_compare == 'no')





#Compare pham-based gene content dissimilarity and vog-based gene content dissimilarity
#Supp. Fig. 2e
par(mar=c(4,8,4,4))
plot(vog_bacteria_dsDNA$pham_pham_dissimilarity,vog_bacteria_dsDNA$vog_gene_content_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")



#Compare pham-based and vog-based bacteria dsDNA phage empirical lifestyle plots 
#Temperate
plot_genomic_similarity_standard(vog_temperate_both)


#Supp Fig. 2f
par(mar=c(4,8,4,4))
plot(vog_temperate_both$modified_mash_distance,vog_temperate_both$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Lytic
plot_genomic_similarity_standard(vog_temperate_neither)

#Supp Fig. 2f
par(mar=c(4,8,4,4))
plot(vog_temperate_neither$modified_mash_distance,vog_temperate_neither$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###Phage nucleic acid type analysis

type_dsDNA <- subset(main_data_table,main_data_table$phage_viral_type_compare == "dsDNA")
type_dsRNA <- subset(main_data_table,main_data_table$phage_viral_type_compare == "dsRNA")
type_ssDNA <- subset(main_data_table,main_data_table$phage_viral_type_compare == "ssDNA")
type_ssRNA <- subset(main_data_table,main_data_table$phage_viral_type_compare == "ssRNA")
type_diff <- subset(main_data_table,main_data_table$phage_viral_type_compare == "different")


#Compare different nucleic acid types

#Supp Fig. 5b
plot_genomic_similarity_standard(type_diff)



#Asses each nucleic acid type

#Supp Fig. 4
plot_genomic_similarity_standard(type_dsDNA)

#Supp Fig. 4
plot_genomic_similarity_standard(type_dsRNA)

#Supp Fig. 4
plot_genomic_similarity_standard(type_ssDNA)

#Supp Fig. 4
plot_genomic_similarity_standard(type_ssRNA)









###Eukaryotic controls
euk_check <- subset(main_data_table,main_data_table$ref_host_superkingdom == "Eukaryota" | main_data_table$query_host_superkingdom == "Eukaryota")
euk_check <- subset(euk_check,euk_check$ref_host_superkingdom == "Bacteria" | euk_check$query_host_superkingdom == "Bacteria")

compute_sector_distribution(euk_check)

plot_genomic_similarity_standard(euk_check)


###Archaea controls
controls <- main_data_table
controls$archaea_one_or_two <- ifelse(controls$ref_host_superkingdom == "Archaea" | controls$query_host_superkingdom == "Archaea",TRUE,FALSE)
controls$archaea_both <- ifelse(controls$ref_host_superkingdom == "Archaea" & controls$query_host_superkingdom == "Archaea",TRUE,FALSE)
controls$archaea_one <- ifelse(controls$archaea_one_or_two == TRUE & controls$archaea_both == FALSE,TRUE,FALSE)
controls$bacteria_one_or_two <- ifelse(controls$ref_host_superkingdom == "Bacteria" | controls$query_host_superkingdom == "Bacteria",TRUE,FALSE)
controls$archaea_one_bacteria_one <- ifelse(controls$archaea_one == TRUE & controls$bacteria_one_or_two == TRUE,TRUE,FALSE)


archaea_check <- subset(controls,controls$archaea_one_bacteria_one == TRUE)
compute_sector_distribution(archaea_check)

#Supp Fig. 5a
plot_genomic_similarity_standard(archaea_check)








###Cluster and subcluster analyses
cluster_actino <- subset(main_data_table,main_data_table$phage_cluster_source_compare == "actino")
cluster_actino_same <- subset(cluster_actino,cluster_actino$phage_cluster_compare != "different")
cluster_actino_diff <- subset(cluster_actino,cluster_actino$phage_cluster_compare == "different")
subcluster_actino_same <- subset(cluster_actino_same,cluster_actino_same$phage_subcluster_compare != "different")
subcluster_actino_diff <- subset(cluster_actino_same,cluster_actino_same$phage_subcluster_compare == "different")

compute_sector_distribution(cluster_actino_same)
compute_sector_distribution(cluster_actino_diff)
compute_sector_distribution(subcluster_actino_same)
compute_sector_distribution(subcluster_actino_diff)

#Supp. Fig. 5d
plot_genomic_similarity_orange(cluster_actino_same)

#Supp. Fig. 5d
plot_genomic_similarity_standard(cluster_actino_diff)

#Supp. Fig. 5d
plot_genomic_similarity_orange(subcluster_actino_same)

#Supp. Fig. 5d
plot_genomic_similarity_orange(subcluster_actino_diff)







###Compare distances for phages that are not subclustered.
cluster_actino <- subset(main_data_table,main_data_table$phage_cluster_source_compare == "actino")

cluster_actino$query_subclustered <- ifelse(is.na(cluster_actino$query_phage_subcluster) == TRUE,FALSE,TRUE) 
cluster_actino$ref_subclustered <- ifelse(is.na(cluster_actino$ref_phage_subcluster) == TRUE,FALSE,TRUE) 
cluster_actino$subclustered_one_or_two <- ifelse(cluster_actino$query_subclustered == TRUE | cluster_actino$ref_subclustered == TRUE,TRUE,FALSE)
cluster_actino$subclustered_neither <- ifelse(cluster_actino$query_subclustered == FALSE & cluster_actino$ref_subclustered == FALSE,TRUE,FALSE)


cluster_actino_same_cluster_neither_subclustered <- subset(cluster_actino,cluster_actino$phage_cluster_compare != "different" &
                                                             cluster_actino$subclustered_neither == TRUE)


compute_sector_distribution(cluster_actino_same_cluster_neither_subclustered)

#Supp Fig. 5e
plot_genomic_similarity_orange(cluster_actino_same_cluster_neither_subclustered)






###Compute distance of all Actinobacteriophage Singletons from other Actinobacteriophages that have been clustered
cluster_actino <- subset(main_data_table,main_data_table$phage_cluster_source_compare == "actino")

cluster_actino_same <- subset(cluster_actino,cluster_actino$phage_cluster_compare != "different")
cluster_actino_diff <- subset(cluster_actino,cluster_actino$phage_cluster_compare == "different")
cluster_actino_diff$query_singleton <- grepl("^Singleton",cluster_actino_diff$query_phage_cluster)
cluster_actino_diff$ref_singleton <- grepl("^Singleton",cluster_actino_diff$ref_phage_cluster)
cluster_actino_diff$singleton_one_or_two <- ifelse(cluster_actino_diff$ref_singleton == TRUE | cluster_actino_diff$query_singleton == TRUE,TRUE,FALSE)

cluster_actino_singleton_one_or_two <- subset(cluster_actino_diff,cluster_actino_diff$singleton_one_or_two == TRUE)

compute_sector_distribution(cluster_actino_singleton_one_or_two)

#Supp Fig. 5f
plot_genomic_similarity_standard(cluster_actino_singleton_one_or_two)










###Predicted lifestyle analysis
bacteria_dsDNA <- subset(main_data_table,main_data_table$host_superkingdom_compare == 'Bacteria' & main_data_table$phage_viral_type_compare == 'dsDNA')
lifestyle_predicted_temperate <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'yes')
lifestyle_predicted_lytic <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'no')

#Supp. Fig. 6b
plot_genomic_similarity_standard(lifestyle_predicted_temperate)

#Supp. Fig. 6b
plot_genomic_similarity_standard(lifestyle_predicted_lytic)





###Predict evolutionary mode
data_for_mode_prediction <- subset(main_data_table,main_data_table$host_superkingdom_compare == "Bacteria" &
                                     main_data_table$phage_viral_type_compare == "dsDNA" &
                                     main_data_table$modified_mash_distance < 0.42 &
                                     main_data_table$pham_pham_dissimilarity < 0.89)
data_for_mode_prediction <- subset(data_for_mode_prediction,select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))
write.table(data_for_mode_prediction,"data_for_mode_prediction.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)



#Run the data through the analyze_mash_network.py script to predict the evolutionary mode.
#Then import the mode prediction back into R.
#Format
#0 = phage
#1 = hgcf_tally
#2 = lgcf_tally
#3 = out_of_range_tally
#4 = hgcf_percent
#5 = lgcf_percent"
#6 = mode_exact	
#7 = mode_approx_80_percent
mode_prediction_table <- read.csv("mode_prediction.csv",sep=",",header=TRUE)

names(mode_prediction_table) <- c("phage_identifier",
                                  "mode_prediction_hgcf_tally","mode_prediction_lgcf_tally","mode_prediction_out_of_range_tally",
                                  "mode_prediction_hgcf_percent","mode_prediction_lgcf_percent",
                                  "mode_prediction_exact","mode_prediction_approx_80")

#Supp. Fig. 6d
par(mar=c(4,8,4,4))
hist(mode_prediction_table$mode_prediction_hgcf_percent,col="black",breaks=25,cex.axis=2,ann=FALSE,las=1,ylim=c(0,1400))



















###Phage tail type analysis
bacteriophages <- subset(main_data_table,main_data_table$host_superkingdom_compare == "Bacteria")
caudovirales <- subset(bacteriophages,bacteriophages$phage_order_compare == "Caudovirales")
caudovirales_family <- subset(caudovirales,caudovirales$phage_family_compare != "different")
myo <- subset(caudovirales_family,caudovirales_family$phage_family_compare == "Myoviridae")
sipho <- subset(caudovirales_family,caudovirales_family$phage_family_compare == "Siphoviridae")
podo <- subset(caudovirales_family,caudovirales_family$phage_family_compare == "Podoviridae")

myo_temperate <- subset(myo,myo$phage_temperate_compare == "yes")
myo_lytic <- subset(myo,myo$phage_temperate_compare == "no")
myo_lifestyle_diff <- subset(myo,myo$phage_temperate_compare == "different")

sipho_temperate <- subset(sipho,sipho$phage_temperate_compare == "yes")
sipho_lytic <- subset(sipho,sipho$phage_temperate_compare == "no")
sipho_lifestyle_diff <- subset(sipho,sipho$phage_temperate_compare == "different")

podo_temperate <- subset(podo,podo$phage_temperate_compare == "yes")
podo_lytic <- subset(podo,podo$phage_temperate_compare == "no")
podo_lifestyle_diff <- subset(podo,podo$phage_temperate_compare == "different")


caudovirales_family_diff <- subset(caudovirales,caudovirales$phage_family_compare == "different")

caudovirales_family_diff$myo_one <- ifelse(caudovirales_family_diff$ref_phage_family == "Myoviridae" | caudovirales_family_diff$query_phage_family == "Myoviridae",TRUE,FALSE)
caudovirales_family_diff$podo_one <- ifelse(caudovirales_family_diff$ref_phage_family == "Podoviridae" | caudovirales_family_diff$query_phage_family == "Podoviridae",TRUE,FALSE)
caudovirales_family_diff$sipho_one <- ifelse(caudovirales_family_diff$ref_phage_family == "Siphoviridae" | caudovirales_family_diff$query_phage_family == "Siphoviridae",TRUE,FALSE)

caudovirales_family_diff$sipho_myo <- ifelse(caudovirales_family_diff$sipho_one == TRUE & caudovirales_family_diff$myo_one == TRUE,TRUE,FALSE)
caudovirales_family_diff$sipho_podo <- ifelse(caudovirales_family_diff$sipho_one == TRUE & caudovirales_family_diff$podo_one == TRUE,TRUE,FALSE)
caudovirales_family_diff$myo_podo <- ifelse(caudovirales_family_diff$myo_one == TRUE & caudovirales_family_diff$podo_one == TRUE,TRUE,FALSE)

caudovirales_family_diff_sipho_myo <- subset(caudovirales_family_diff,caudovirales_family_diff$sipho_myo == TRUE)
caudovirales_family_diff_sipho_podo <- subset(caudovirales_family_diff,caudovirales_family_diff$sipho_podo == TRUE)
caudovirales_family_diff_myo_podo <- subset(caudovirales_family_diff,caudovirales_family_diff$myo_podo == TRUE)


#Supp Fig. 7a
plot_genomic_similarity_standard(myo)

#Supp Fig. 7a
plot_genomic_similarity_standard(sipho)

#Supp Fig. 7a
plot_genomic_similarity_standard(podo)


#Supp Fig. 7a
plot_genomic_similarity_standard(podo_temperate)

#Supp Fig. 7a
plot_genomic_similarity_standard(podo_lytic)

#Supp Fig. 7a
plot_genomic_similarity_standard(podo_lifestyle_diff)


#Supp Fig. 7a
plot_genomic_similarity_standard(sipho_temperate)

#Supp Fig. 7a
plot_genomic_similarity_standard(sipho_lytic)

#Supp Fig. 7a
plot_genomic_similarity_standard(sipho_lifestyle_diff)



#Supp Fig. 7a
plot_genomic_similarity_standard(myo_temperate)

#Supp Fig. 7a
plot_genomic_similarity_standard(myo_lytic)

#Supp Fig. 7a
plot_genomic_similarity_standard(myo_lifestyle_diff)



#Supp Fig. 7b
plot_genomic_similarity_standard(caudovirales_family_diff_sipho_myo)

#Supp Fig. 7b
plot_genomic_similarity_standard(caudovirales_family_diff_sipho_podo)

#Supp Fig. 7b
plot_genomic_similarity_standard(caudovirales_family_diff_myo_podo)








###Export data for gene-specific mash analysis
#Only investigate comparisons that are within cluster boundaries
mash_filtered <- subset(main_data_table,main_data_table$filter == TRUE)
bacteria_dsDNA_filtered <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')

#Cluster-boundary dataset
bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA_filtered,bacteria_dsDNA_filtered$modified_mash_distance < 0.42 & bacteria_dsDNA_filtered$pham_pham_dissimilarity < 0.89)

#Now reduce the data table to only the columns needed for export, and export the data
bacteria_dsDNA_nuc042_gene089_reduced <- subset(bacteria_dsDNA_nuc042_gene089,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
write.table(bacteria_dsDNA_nuc042_gene089_reduced,"bacteria_dsDNA_nuc042_gene089_data.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)






###Gene-specific mash analysis
#After running gene-specific mash analysis, subset out only the comparisons used in the gsm analysis

mash_filtered <- subset(main_data_table,main_data_table$filter == TRUE)
bacteria_dsDNA_filtered <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')


#Cluster-boundary dataset
gene_specific_mash_table <- subset(bacteria_dsDNA_filtered,bacteria_dsDNA_filtered$modified_mash_distance < 0.42 & bacteria_dsDNA_filtered$pham_pham_dissimilarity < 0.89)


#Now import the gsm data
#Format:
#0 = phage1_phage2
#1 = phage1
#2 = phage2
#3 = phage1_phage2 # shared phams
#4 = phage1_phage2 gene content dissimilarity (general index)
#5 = phage1_phage2 gene content dissimilarity (jaccard index)
#6 = phage1_phage2 all genes mash distance
#7 = phage1_phage2 all genes mash p-value
#8 = phage1_phage2 all genes mash kmer count
#9 = phage1_phage2 shared genes mash distance
#10 = phage1_phage2 shared genes mash p-value
#11 = phage1_phage2 shared genes mash kmer count
#12 = phage1_phage2 unshared genes mash distance
#13 = phage1_phage2 unshared genes mash p-value
#14 = phage1_phage2 unshared genes mash kmer count
#15 = phage1_phage2 shared-unshared mash distance
#16 = phage1_phage2 shared-unshared mash p-value
#17 = phage1_phage2 shared-unshared mash kmer count
#18 = phage1_phage2 total length of combined shared sequence
#19 = phage1_phage2 total length of combined unshared sequence
#20 = phage1_phage2 combined shared gene GC content
#21 = phage1_phage2 combined unshared gene GC content
#22 = phage1 number unshared phams
#23 = phage1 number all phams
#24 = phage1 number all genes
#25 = phage1 number shared genes
#26 = phage1 number unshared genes
#27 = phage1 average length of all genes
#28 = phage1 average length of shared genes
#29 = phage1 average length of unshared genes
#30 = phage1 total length of all genes
#31 = phage1 total length of shared genes
#32 = phage1 total length of unshared genes
#33 = phage1 all genes GC content
#34 = phage1 shared genes GC content
#35 = phage1 unshared genes GC content
#36 = phage2 number unshared phams
#37 = phage2 number all phams
#38 = phage2 number all genes
#39 = phage2 number shared genes
#40 = phage2 number unshared genes
#41 = phage2 average length of all genes
#42 = phage2 average length of shared genes
#43 = phage2 average length of unshared genes
#44 = phage2 total length of all genes
#45 = phage2 total length of shared genes
#46 = phage2 total length of unshared genes
#47 = phage2 all genes GC content
#48 = phage2 shared genes GC content
#49 = phage2 unshared genes GC content
gene_specific_mash_data <- read.csv("gene_specific_mash_data.csv",sep=",",header=TRUE)

names(gene_specific_mash_data) = c("gsm_mash_ref_query","gsm_ref","gsm_query","gsm_num_shared_phams","gsm_pham_dissimilarity","gsm_pham_jaccard_dissimilarity",
                                    "gsm_all_mash_distance","gsm_all_mash_pvalue","gsm_all_mash_count",
                                    "gsm_shared_mash_distance","gsm_shared_mash_pvalue","gsm_shared_mash_count",
                                    "gsm_unshared_mash_distance","gsm_unshared_mash_pvalue","gsm_unshared_mash_count",
                                    "gsm_shared_unshared_mash_distance","gsm_shared_unshared_mash_pvalue","gsm_shared_unshared_mash_count",
                                    "gsm_shared_total_size","gsm_unshared_total_size","gsm_shared_GC","gsm_unshared_GC",
                                    "gsm_ref_num_unshared_phams","gsm_ref_num_all_phams",
                                    "gsm_ref_num_all_genes","gsm_ref_num_shared_genes","gsm_ref_num_unshared_genes",
                                    "gsm_ref_all_ave_size","gsm_ref_shared_ave_size","gsm_ref_unshared_ave_size",
                                    "gsm_ref_all_total_size","gsm_ref_shared_total_size","gsm_ref_unshared_total_size",
                                    "gsm_ref_all_GC","gsm_ref_shared_GC","gsm_ref_unshared_GC",
                                    "gsm_query_num_unshared_phams","gsm_query_num_all_phams",
                                    "gsm_query_num_all_genes","gsm_query_num_shared_genes","gsm_query_num_unshared_genes",
                                    "gsm_query_all_ave_size","gsm_query_shared_ave_size","gsm_query_unshared_ave_size",
                                    "gsm_query_all_total_size","gsm_query_shared_total_size","gsm_query_unshared_total_size",
                                    "gsm_query_all_GC","gsm_query_shared_GC","gsm_query_unshared_GC")


#Merge the gene-specific mash data to the main data table that has been reduced to contain
#only those comparisons analyzed by the gene-specific mash script.
#Since the ref_query comparison identifiers should be identical, there should be equal rows in both tables
gene_specific_mash_table <- merge(gene_specific_mash_table,gene_specific_mash_data,by.x="mash_ref_query",by.y="gsm_mash_ref_query")



#Choose which data should be retained. 
#Data not passing the filter get Mash distances re-assigned to 0.6.
gene_specific_mash_table$gsm_all_mash_filter <- ifelse(gene_specific_mash_table$gsm_all_mash_pvalue < 1e-10,TRUE,FALSE)
gene_specific_mash_table$gsm_shared_mash_filter <- ifelse(gene_specific_mash_table$gsm_shared_mash_pvalue < 1e-10,TRUE,FALSE)
gene_specific_mash_table$gsm_unshared_mash_filter <- ifelse(gene_specific_mash_table$gsm_unshared_mash_pvalue < 1e-10,TRUE,FALSE)
gene_specific_mash_table$gsm_shared_unshared_mash_filter <- ifelse(gene_specific_mash_table$gsm_shared_unshared_mash_pvalue < 1e-10,TRUE,FALSE)

gene_specific_mash_table$gsm_all_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_all_mash_filter == TRUE,gene_specific_mash_table$gsm_all_mash_distance,0.6)
gene_specific_mash_table$gsm_shared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_shared_mash_filter == TRUE,gene_specific_mash_table$gsm_shared_mash_distance,0.6)
gene_specific_mash_table$gsm_unshared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_unshared_mash_filter == TRUE,gene_specific_mash_table$gsm_unshared_mash_distance,0.6)
gene_specific_mash_table$gsm_shared_unshared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_shared_unshared_mash_filter == TRUE,gene_specific_mash_table$gsm_shared_unshared_mash_distance,0.6)


#Compute of the proportion of coding sequence per genome
#Note: the gene-specific sequence length used does not take into account overlapping CDS features, 
#so the sequences could have duplicate regions in the genome.
#Note: the 'all coding proportion' data is based on the real genome size, 
#but the 'shared/unshared coding proportion' data is based only on the gene-specific mash sizes
gene_specific_mash_table$gsm_all_total_size <- gene_specific_mash_table$gsm_ref_all_total_size + 
                                                gene_specific_mash_table$gsm_query_all_total_size

gene_specific_mash_table$gsm_ref_unshared_coding_proportion <- gene_specific_mash_table$gsm_ref_unshared_total_size/
                                                                gene_specific_mash_table$gsm_ref_all_total_size

gene_specific_mash_table$gsm_query_unshared_coding_proportion <- gene_specific_mash_table$gsm_query_unshared_total_size/
                                                                  gene_specific_mash_table$gsm_query_all_total_size

gene_specific_mash_table$gsm_unshared_coding_proportion <- gene_specific_mash_table$gsm_unshared_total_size/
                                                            gene_specific_mash_table$gsm_all_total_size


#Average gene sizes
#To compute average gene size, sum up all nucleotide sizes then divide by total number of genes
gene_specific_mash_table$gsm_ave_size_shared_shared <- ((gene_specific_mash_table$gsm_ref_shared_ave_size * gene_specific_mash_table$gsm_ref_num_shared_genes) +
                                                          (gene_specific_mash_table$gsm_query_shared_ave_size * gene_specific_mash_table$gsm_query_num_shared_genes))/
                                                          (gene_specific_mash_table$gsm_ref_num_shared_genes + gene_specific_mash_table$gsm_query_num_shared_genes)
  
gene_specific_mash_table$gsm_ave_size_unshared_unshared <- ((gene_specific_mash_table$gsm_ref_unshared_ave_size * gene_specific_mash_table$gsm_ref_num_unshared_genes) +
                                                          (gene_specific_mash_table$gsm_query_unshared_ave_size * gene_specific_mash_table$gsm_query_num_unshared_genes))/
                                                          (gene_specific_mash_table$gsm_ref_num_unshared_genes + gene_specific_mash_table$gsm_query_num_unshared_genes)

                                                          


#Convert all Unspecified fields to NA missing value
gene_specific_mash_table[gene_specific_mash_table == "Unspecified"] <- NA


#Split into HGCF and LGCF datasets based on position of data points on the plot
gene_specific_mash_table_temperate_hgcf <- subset(gene_specific_mash_table,gene_specific_mash_table$gene_flux_category == 'high' & gene_specific_mash_table$phage_temperate_compare == 'yes')
gene_specific_mash_table_temperate_lgcf <- subset(gene_specific_mash_table,gene_specific_mash_table$gene_flux_category == 'low' & gene_specific_mash_table$phage_temperate_compare == 'yes')
gene_specific_mash_table_lytic <- subset(gene_specific_mash_table,gene_specific_mash_table$phage_temperate_compare == 'no')






#Analysis after splitting into gene flux modes

#All data points used in gene-specific mash plots
#Supp. Fig. 8a
plot_genomic_similarity_tricolored(gene_specific_mash_table_temperate_hgcf,gene_specific_mash_table_temperate_lgcf,gene_specific_mash_table_lytic)





#Compare whole genome distance to shared and unshared distances
#Supp. Fig. 8b
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_shared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
abline(0,1,lty=2,lwd=3,col="grey")

#Supp. Fig. 8b
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_shared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")

#Supp. Fig. 8b
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_shared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")






#Coding proportion
#Supp. Fig. 8d
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_unshared_coding_proportion,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
abline(0,1,lty=2,lwd=3,col="grey")

#Supp. Fig. 8d
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$gsm_unshared_coding_proportion,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")

#Supp. Fig. 8d
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$gsm_unshared_coding_proportion,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")




#Compare average gene sizes
gsm_temperate_hgcf_shared <- subset(gene_specific_mash_table_temperate_hgcf,select=c("gsm_ave_size_shared_shared"))
gsm_temperate_lgcf_shared <- subset(gene_specific_mash_table_temperate_lgcf,select=c("gsm_ave_size_shared_shared"))
gsm_lytic_shared <- subset(gene_specific_mash_table_lytic,select=c("gsm_ave_size_shared_shared"))

gsm_temperate_hgcf_unshared <- subset(gene_specific_mash_table_temperate_hgcf,select=c("gsm_ave_size_unshared_unshared"))
gsm_temperate_lgcf_unshared <- subset(gene_specific_mash_table_temperate_lgcf,select=c("gsm_ave_size_unshared_unshared"))
gsm_lytic_unshared <- subset(gene_specific_mash_table_lytic,select=c("gsm_ave_size_unshared_unshared"))

names(gsm_temperate_hgcf_shared) <- c("ave_gene_size")
names(gsm_temperate_lgcf_shared) <- c("ave_gene_size")
names(gsm_lytic_shared) <- c("ave_gene_size")
names(gsm_temperate_hgcf_unshared) <- c("ave_gene_size")
names(gsm_temperate_lgcf_unshared) <- c("ave_gene_size")
names(gsm_lytic_unshared) <- c("ave_gene_size")

gsm_temperate_hgcf_shared$category <- "temperate_HGCF"
gsm_temperate_lgcf_shared$category <- "temperate_LGCF"
gsm_lytic_shared$category <- "lytic"
gsm_temperate_hgcf_unshared$category <- "temperate_HGCF"
gsm_temperate_lgcf_unshared$category <- "temperate_LGCF"
gsm_lytic_unshared$category <- "lytic"

gsm_temperate_hgcf_shared$gene_type <- "shared"
gsm_temperate_lgcf_shared$gene_type <- "shared"
gsm_lytic_shared$gene_type <- "shared"
gsm_temperate_hgcf_unshared$gene_type <- "unshared"
gsm_temperate_lgcf_unshared$gene_type <- "unshared"
gsm_lytic_unshared$gene_type <- "unshared"

gsm_plotting_table <- rbind(gsm_temperate_hgcf_shared,gsm_temperate_lgcf_shared,gsm_lytic_shared,gsm_temperate_hgcf_unshared,gsm_temperate_lgcf_unshared,gsm_lytic_unshared)

#Supp. Fig. 8e
par(mar=c(10,4,4,4))
boxplot(gsm_plotting_table$ave_gene_size ~ gsm_plotting_table$category*gsm_plotting_table$gene_type,las=2,col=c("dark red","dark blue","dark green","pink","light blue","light green"))











#Analysis after splitting into gene flux modes
#Sliding window average across ordered mash distances
library(caTools)

gene_specific_mash_table_temperate_hgcf_mmdsort <- gene_specific_mash_table_temperate_hgcf[order(gene_specific_mash_table_temperate_hgcf$modified_mash_distance),]
gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent,101)

gene_specific_mash_table_temperate_lgcf_mmdsort <- gene_specific_mash_table_temperate_lgcf[order(gene_specific_mash_table_temperate_lgcf$modified_mash_distance),]
gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent,101)

gene_specific_mash_table_lytic_mmdsort <- gene_specific_mash_table_lytic[order(gene_specific_mash_table_lytic$modified_mash_distance),]
gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent,101)


#Compare whole genome size disparity
#Supp. Fig. 8a
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#Sliding window analyis across ordered gene content dissimilarity instead of nucleotide distance
gene_specific_mash_table_temperate_hgcf_gcdsort <- gene_specific_mash_table_temperate_hgcf[order(gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity),]
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance,101)


gene_specific_mash_table_temperate_lgcf_gcdsort <- gene_specific_mash_table_temperate_lgcf[order(gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity),]
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance,101)


gene_specific_mash_table_lytic_gcdsort <- gene_specific_mash_table_lytic[order(gene_specific_mash_table_lytic$pham_pham_dissimilarity),]
gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance,101)


#Compare shared and unshared distances
#Supp. Fig. 8c
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance_runmean,gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance_runmean,gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance_runmean,gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance_runmean,gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance_runmean,gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance_runmean,gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")





###Phylogeny comparison 
#Note: Requires that the gene-specific mash data has been loaded and processed
#Format:
#0 = ref_query
#1 = phylogeny distance
phylogeny_data <- read.csv("phylogeny_data.csv",sep=",",header=TRUE)
phylogeny_analysis <- merge(gene_specific_mash_table,phylogeny_data,by.x="mash_ref_query",by.y="ref_query")





#Compare mash distance to phylogeny distance

#Split by Clusters to visualize by color
phylogeny_hgcf <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'F' |
                           (phylogeny_analysis$phage_cluster_compare == 'A' & phylogeny_analysis$phage_subcluster_compare == 'A1'))


phylogeny_lgcf <- subset(phylogeny_analysis,
                         phylogeny_analysis$phage_cluster_compare == 'K' |
                           phylogeny_analysis$phage_cluster_compare == 'BD' |
                           (phylogeny_analysis$phage_cluster_compare == 'A' & phylogeny_analysis$phage_subcluster_compare != 'A1'))


phylogeny_lytic <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'B')



#Supp. Fig. 9a
plot_genomic_similarity_tricolored(phylogeny_hgcf,phylogeny_lgcf,phylogeny_lytic)


#Supp. Fig 9a
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")




#Cluster A phylogenetic distance analysis
phylogeny_analysis_a <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'A')
phylogeny_analysis_a_both_a1 <- subset(phylogeny_analysis_a,phylogeny_analysis_a$phage_subcluster_compare == 'A1')
phylogeny_analysis_a_both_nonA1 <- subset(phylogeny_analysis_a,phylogeny_analysis_a$phage_subcluster_compare != 'A1')
phylogeny_analysis_a_one_a1 <- subset(phylogeny_analysis_a,(phylogeny_analysis_a$ref_phage_subcluster == 'A1' | phylogeny_analysis_a$query_phage_subcluster == 'A1') & phylogeny_analysis_a$phage_subcluster_compare == "different")


#Fig. 3a
par(mar=c(4,8,4,4))
plot(phylogeny_analysis_a_both_nonA1$phylogeny_distance,phylogeny_analysis_a_both_nonA1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_analysis_a_both_a1$phylogeny_distance,phylogeny_analysis_a_both_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(phylogeny_analysis_a_one_a1$phylogeny_distance,phylogeny_analysis_a_one_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")


#Fig. 3c
par(mar=c(4,8,4,4))
plot(phylogeny_analysis_a_both_a1$phylogeny_distance,phylogeny_analysis_a_both_a1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(phylogeny_analysis_a_one_a1$phylogeny_distance,phylogeny_analysis_a_one_a1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")
par(new=TRUE)
plot(phylogeny_analysis_a_both_nonA1$phylogeny_distance,phylogeny_analysis_a_both_nonA1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")

  






  
#Match up pham proportion data based only on clustered actinobacteriophages (785 genomes).
#Columns are designated as "pham2" since it is the second set of pham data that have been loaded so far. 
#This second pham data is similar, but not identical, to the first.
#Some metrics, such as pham distribution and orpham count, are impacted by the diversity in the dataset, so these values may
#not be the same as in the first pham dataset since it only contains comparisons from the actino785 set.

#Format
#0 = phage1_name
#1 = phage1_number_of_unshared_phams
#2 = phage1_shared_proportion
#3 = phage2_name
#4 = phage2_number_of_unshared_phams
#5 = phage2_shared_proportion
#6 = number_of_shared_phams
#7 = average_shared_proportion
#8 = jaccard_similarity
#9 = shared_pham_distribution_mean
#10 = shared_pham_distribution_median
#11 = shared_pham_distribution_max
#12 = unshared_pham_distribution_mean
#13 = unshared_pham_distribution_median
#14 = unshared_pham_distribution_max
#15 = unshared_orpham_count
actino_pham_data <- read.csv("actino_only_pairwise_pham_proportions.csv",sep=",",header=TRUE)

names(actino_pham_data) <- c("pham2_phage1","pham2_phage1_number_unshared_phams","pham2_phage1_shared_proportion","pham2_phage2",
                             "pham2_phage2_number_unshared_phams","pham2_phage2_shared_proportion","pham2_number_shared_phams","pham2_average_shared_proportion",
                             "pham2_jaccard_similarity","pham2_shared_pham_distribution_mean","pham2_shared_pham_distribution_median","pham2_shared_pham_distribution_max",
                             "pham2_unshared_pham_distribution_mean","pham2_unshared_pham_distribution_median",
                             "pham2_unshared_pham_distribution_max","pham2_unshared_orpham_count")

#Compute gene content dissimilarity
actino_pham_data$pham2_pham_dissimilarity <- 1 - actino_pham_data$pham2_average_shared_proportion
actino_pham_data$pham2_jaccard_dissimilarity <- 1 - actino_pham_data$pham2_jaccard_similarity


#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
actino_pham_data$pham2_phage1_phage2 <- paste(actino_pham_data$pham2_phage1,"_",actino_pham_data$pham2_phage2,sep="")
actino_pham_data$pham2_phage1_phage2 <- as.factor(actino_pham_data$pham2_phage1_phage2)






#All rows don't need to be retained = making scatter plots or histograms can cause errors if some rows are missing data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
phylogeny_analysis <- merge(phylogeny_analysis,actino_pham_data,by.x="mash_ref_query",by.y="pham2_phage1_phage2")


#Split by Clusters to visualize by color
phylogeny_hgcf2 <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'F' |
                           (phylogeny_analysis$phage_cluster_compare == 'A' & phylogeny_analysis$phage_subcluster_compare == 'A1'))

phylogeny_lgcf2 <- subset(phylogeny_analysis,
                         phylogeny_analysis$phage_cluster_compare == 'K' |
                           phylogeny_analysis$phage_cluster_compare == 'BD' |
                           (phylogeny_analysis$phage_cluster_compare == 'A' & phylogeny_analysis$phage_subcluster_compare != 'A1'))

phylogeny_lytic2 <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'B')




#Sliding window average across mash distance
#Only use data for sliding windows that are positioned within the scatter plot boundaries = phylogeny distance < 0.3.
library(caTools)

phylogeny_hgcf_phylosort <- phylogeny_hgcf2[order(phylogeny_hgcf2$phylogeny_distance),]
phylogeny_hgcf_phylosort <- phylogeny_hgcf_phylosort[phylogeny_hgcf_phylosort$phylogeny_distance < 0.3,]
phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_hgcf_phylosort$size_diff_ave_percent,101)
phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_hgcf_phylosort$gsm_unshared_coding_proportion_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_unshared_coding_proportion,101)
phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared,101)



phylogeny_lgcf_phylosort <- phylogeny_lgcf2[order(phylogeny_lgcf2$phylogeny_distance),]
phylogeny_lgcf_phylosort <- phylogeny_lgcf_phylosort[phylogeny_lgcf_phylosort$phylogeny_distance < 0.3,]
phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lgcf_phylosort$size_diff_ave_percent,101)
phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lgcf_phylosort$gsm_unshared_coding_proportion_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_unshared_coding_proportion,101)
phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared,101)


phylogeny_lytic_phylosort <- phylogeny_lytic2[order(phylogeny_lytic2$phylogeny_distance),]
phylogeny_lytic_phylosort <- phylogeny_lytic_phylosort[phylogeny_lytic_phylosort$phylogeny_distance < 0.3,]
phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lytic_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lytic_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lytic_phylosort$size_diff_ave_percent,101)
phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lytic_phylosort$gsm_unshared_coding_proportion_runmean <- runmean(phylogeny_lytic_phylosort$gsm_unshared_coding_proportion,101)
phylogeny_lytic_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_shared_shared,101)






#Gene content dissimilarity
#Supp. Fig. 9b
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Ave genome size difference
#Supp. Fig. 9b
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Unshared coding sequence proportion
#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_unshared_coding_proportion_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_unshared_coding_proportion_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_unshared_coding_proportion_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


#Shared/unshared ave gene size
#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")


#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



#Shared/unshared pham distribution mean
#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")

#Unshared orpham count
#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")








###Lambda analysis
#Phage identifier for lambda: lambda__nc_001416
#Since no self comparisons are present in the dataset, only need to compute if there's 'one' lambda, instead of 'one_or_two' or 'both'
lambda_figure <- main_data_table
lambda_figure$lambda_one <- ifelse(lambda_figure$mash_reference == 'lambda__nc_001416' | lambda_figure$mash_query == 'lambda__nc_001416',TRUE,FALSE)
lambda_comparisons <- subset(lambda_figure,lambda_figure$lambda_one == TRUE)

#Supp. Fig. 11a
plot_genomic_similarity_standard(lambda_comparisons)



###Toxic phage analysis
main_data_table_toxic <- main_data_table
main_data_table_toxic$toxic_one_or_two <- ifelse(main_data_table_toxic$ref_toxic == "yes" | main_data_table_toxic$query_toxic == "yes",TRUE,FALSE)
main_data_table_toxic$toxic_both <- ifelse(main_data_table_toxic$ref_toxic == "yes" & main_data_table_toxic$query_toxic == "yes",TRUE,FALSE)
toxic_one_or_two <- subset(main_data_table_toxic,main_data_table_toxic$toxic_one_or_two == TRUE)

#Supp. Fig. 11b
plot_genomic_similarity_standard(toxic_one_or_two)






