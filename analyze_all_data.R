#Script to analyze merged2333 phage classification mash data
#20161006
#Updated script to account for new phage and host data fields
#20161108
#Revised again
#20161121
#Revise again to include all data points, even non-significant ones
#20170221
#Revise again to include VOG data

#Function to process data
import_function <- function(kmer_file,kmer){
  
  kmer_table <- read.csv(kmer_file,sep=",",header=TRUE)
  
  #Convert fields from factor class (default) to character class
  kmer_table$reference <- as.character(kmer_table$reference)
  kmer_table$query <- as.character(kmer_table$query)
  kmer_table$ref_query <- as.character(kmer_table$ref_query)
  
  #Create column headers and assign them
  reference <- paste(kmer,"reference",sep="_")
  query <- paste(kmer,"query",sep="_")
  dist <- paste(kmer,"distance",sep="_")
  pvalue <- paste(kmer,"pvalue",sep="_")
  count <- paste(kmer,"count",sep="_")
  ref_query <- paste(kmer,"ref_query",sep="_")
  column_headers <- c(reference,query,dist,pvalue,count,ref_query)
  names(kmer_table) <- column_headers
  return(kmer_table)
  
  
  
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



#setwd("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/")


#Decide on which dataset to import
#Main dataset
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20161006_25000sketch_15kmer_mash_processed.csv"



##OR##
#Import alternative mash dataset
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170407_5000sketch_13kmer_mash_processed.csv"
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170407_5000sketch_17kmer_mash_processed.csv"
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170407_25000sketch_17kmer_mash_processed.csv"
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170410_25000sketch_13kmer_mash_processed.csv"
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170410_25000sketch_15kmer_mash_processed.csv"
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170411_50000sketch_15kmer_mash_processed.csv"
filename <- "/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/0_processed_mash_datasets/20170413_50000sketch_17kmer_mash_processed.csv"
##



#Import datasets
mash_table <- import_function(filename,"mash")






#Convert fields from character class (default) to factor class
mash_table$mash_reference <- as.factor(mash_table$mash_reference)
mash_table$mash_query <- as.factor(mash_table$mash_query)
mash_table$mash_ref_query <- as.factor(mash_table$mash_ref_query)




#Import host taxonomy data
host_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/merged2333_phage_host_data11.csv",sep=",",header=TRUE)


#Convert all Unspecified fields to NA missing value
host_table[host_table == "Unspecified"] <- NA



#Modify column names of host data and merge with kmer table
#Data for all phages in kmer15 table should be present in host table and size table.
#As a result, do not select all.x=TRUE option. This way, any missing rows indicates an error in the table.


#Original host table column header names:
#"phage_identifier","header_source_info","database","host_superkingdom","host_phylum","host_class","host_order","host_family","host_genus","phage_superkingdom","phage_viral_type","phage_order","phage_family","phage_genus","phage_cluster","phage_subcluster","phage_cluster_source","phage_temperate","size" 



names(host_table) <- c("phage_identifier","query_header_source_info","query_database",
                             "query_host_superkingdom","query_host_phylum","query_host_class","query_host_order","query_host_family","query_host_genus",
                             "query_phage_superkingdom","query_phage_viral_type","query_phage_order","query_phage_family","query_phage_genus",
                             "query_phage_cluster","query_phage_subcluster","query_phage_cluster_source",
                             "query_phage_temperate","query_size","query_rnrdb","query_phage_rnr","query_phamerator_status","query_gene_count",
                             "query_toxic","query_extra_chrome","query_predicted_temperate")




mash_table2 <- merge(mash_table,host_table,by.x="mash_query",by.y="phage_identifier")

names(host_table) <- c("phage_identifier","ref_header_source_info","ref_database",
                       "ref_host_superkingdom","ref_host_phylum","ref_host_class","ref_host_order","ref_host_family","ref_host_genus",
                       "ref_phage_superkingdom","ref_phage_viral_type","ref_phage_order","ref_phage_family","ref_phage_genus",
                       "ref_phage_cluster","ref_phage_subcluster","ref_phage_cluster_source",
                       "ref_phage_temperate","ref_size","ref_rnrdb","ref_phage_rnr","ref_phamerator_status","ref_gene_count",
                       "ref_toxic","ref_extra_chrome","ref_predicted_temperate")



mash_table2 <- merge(mash_table2,host_table,by.x="mash_reference",by.y="phage_identifier")

names(host_table) <- c("phage_identifier","header_source_info","database",
                       "host_superkingdom","host_phylum","host_class","host_order","host_family","host_genus",
                       "phage_superkingdom","phage_viral_type","phage_order","phage_family","phage_genus",
                       "phage_cluster","phage_subcluster","phage_cluster_source",
                       "phage_temperate","size","rnrdb","phage_rnr","phamerator_status","gene_count",
                       "toxic","extra_chrome","predicted_temperate")


#At this point, the ref and query cluster and taxonomic ranking columns are factor class






#Compare genome sizes
mash_table2$size_diff <- abs(mash_table2$ref_size - mash_table2$query_size)
mash_table2$size_diff_ref_percent <- mash_table2$size_diff / mash_table2$ref_size
mash_table2$size_diff_query_percent <- mash_table2$size_diff / mash_table2$query_size
mash_table2$size_diff_min_percent <- ifelse(mash_table2$size_diff_ref_percent < mash_table2$size_diff_query_percent,mash_table2$size_diff_ref_percent,mash_table2$size_diff_query_percent)
mash_table2$size_diff_max_percent <- ifelse(mash_table2$size_diff_ref_percent > mash_table2$size_diff_query_percent,mash_table2$size_diff_ref_percent,mash_table2$size_diff_query_percent)
mash_table2$size_diff_ave_percent <- (mash_table2$size_diff_ref_percent + mash_table2$size_diff_query_percent)/2



#Compare difference in gene abundance
mash_table2$gene_count_disparity <- abs(mash_table2$ref_gene_count - mash_table2$query_gene_count)
mash_table2$gene_count_disparity_ave_percent <- ((mash_table2$gene_count_disparity/mash_table2$ref_gene_count) + (mash_table2$gene_count_disparity/mash_table2$query_gene_count))/2



#Import ANI data to check Mash vs Pham of specific clusters
#Original file
#ani_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161107_ani_data.csv",sep=",",header=TRUE)
#New file containing cluster f data
ani_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170330_ani_data.csv",sep=",",header=TRUE)


ani_data$ani_ref_query <- as.character(ani_data$ani_ref_query)
ani_data$ani_distance <- 1 - ani_data$ani_ani

#Merge tables
mash_table2 <- merge(mash_table2,ani_data,by.x="mash_ref_query",by.y="ani_ref_query",all.x=TRUE)



#Import pham data and merge with kmer table
#Bacteriophages database pham table data:
#pham_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161121_bacteriophages_v5_pairwise_pham_proportions.csv",sep=",",header=TRUE)
##OR##
pham_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170505_bacteriophages_2333_v5_pairwise_pham_proportions.csv",sep=",",header=TRUE)


pham_table$pham_dissimilarity <- 1 - pham_table$average_shared_proportion
pham_table$jaccard_dissimilarity <- 1 - pham_table$jaccard_similarity

#20161121 pham table contains fewer output columns, since it was generated using an older version of the code
# names(pham_table) <- c("pham_phage1","pham_phage1_number_unshared_phams","pham_phage1_shared_proportion","pham_phage2",
#                        "pham_phage2_number_unshared_phams","pham_phage2_shared_proportion","pham_number_shared_phams","pham_averaged_shared_proportion",
#                        "pham_jaccard_similarity","pham_pham_dissimilarity","pham_jaccard_dissimilarity")



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

#Some genomes have no annotated genes, so no phams (e.g. vb_paem_c1-14_ab28__nc_026600).
#This amounts to 4663 rows.
#To retain all rows, be sure to keep all.x=TRUE, but you don't want to retain all rows = making scatter plots or histograms can cause errors if not all rows have data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
mash_table2 <- merge(mash_table2,pham_table,by.x="mash_ref_query",by.y="pham_phage1_phage2")








#Assign filter status and change mash distance if data is not significant
#Choose which data should be retained. Those not passing the filter getting Mash distances re-assigned to 0.5.


#Assess the distribution of mash values before and after p-value filter
mash_table2_pvalue <- subset(mash_table2,mash_table2$mash_pvalue < 1e-10)
mash_table2_pvalue_mash_above05 <- subset(mash_table2_pvalue,mash_table2_pvalue$mash_distance >= 0.5)
nrow(mash_table2_pvalue_mash_above05)
summary(mash_table2_pvalue$mash_distance)
filename

#Old filtering criteria
#mash_table2$filter <- ifelse(mash_table2$mash_pvalue < 1e-10 & mash_table2$size_diff_max_percent < 1,TRUE,FALSE)

#New filtering criteria
#If pvalue is the only filtering criteria, there will be a single comparison involving a dsDNA and dsRNA phage with distance of ~0.509. So set any value > 0.5 to 0.5.
mash_table2$filter <- ifelse(mash_table2$mash_pvalue < 1e-10 & mash_table2$mash_distance <= 0.5,TRUE,FALSE)



#At this point, the max mash distance of all filtered comparisons = 0.4903. So set the distance of all comparisons that did not pass the filter = 0.5
mash_table2$modified_mash_distance <- ifelse(mash_table2$filter == TRUE,mash_table2$mash_distance,0.5)

summary(mash_table2$mash_distance)
summary(mash_table2$modified_mash_distance)


#
#
#
##TEMP assess whether genome size disparity should really be accounted for
par(mar=c(4,8,4,4))
plot(mash_table2$modified_mash_distance,mash_table2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

temp <- subset(mash_table2,mash_table2$ref_host_superkingdom == 'Bacteria' &
                 mash_table2$query_host_superkingdom == 'Bacteria' &
                 mash_table2$ref_phage_viral_type == 'dsDNA' &
                 mash_table2$query_phage_viral_type == 'dsDNA')

par(mar=c(4,8,4,4))
plot(temp$modified_mash_distance,temp$pham_pham_dissimilarity,xlim=c(0,0.51),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

temp2 <- mash_table2
temp2$filter2 <- ifelse(temp2$mash_pvalue < 1e-10,TRUE,FALSE)
retained_by_pvalue <- subset(mash_table2,mash_table2$mash_pvalue < 1e-10)
retained_by_size <- subset(retained_by_pvalue,retained_by_pvalue$size_diff_max_percent < 1)
removed_by_size <- subset(retained_by_pvalue,retained_by_pvalue$size_diff_max_percent >= 1)

original_retained <- subset(mash_table2,mash_table2$mash_pvalue < 1e-10 & mash_table2$size_diff_max_percent < 1)

par(mar=c(4,8,4,4))
plot(original_retained$mash_distance,original_retained$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(retained_by_size$mash_distance,retained_by_size$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(removed_by_size$mash_distance,removed_by_size$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

compute_sector_distribution(removed_by_size)

retained_by_pvalue2 <- subset(retained_by_pvalue,retained_by_pvalue$mash_pvalue < 1e-11)
removed_by_pvalue2 <- subset(retained_by_pvalue,retained_by_pvalue$mash_pvalue >= 1e-11)

par(mar=c(4,8,4,4))
plot(retained_by_pvalue2$mash_distance,retained_by_pvalue2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(removed_by_pvalue2$mash_distance,removed_by_pvalue2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


retained_by_pvalue3 <- subset(retained_by_pvalue,retained_by_pvalue$mash_pvalue < 1e-15)


retained_by_value_viral_type <- subset(retained_by_pvalue,retained_by_pvalue$ref_phage_viral_type == retained_by_pvalue$query_phage_viral_type)
#
#
#




#Assign evolutionary mode
#This assigns all data in the dataset. However, evolutionary mode is only relevant to dsDNA
#hgf = high gene flux
#lgf = low gene flux
#Lines to demarcate groups
#abline(1.35,-2,lty=2,lwd=3,col="grey")
#abline(0.4,-2,lty=2,lwd=3,col="grey")
#abline(0,3.5,lty=2,lwd=3,col="grey")
#abline(0.25,2,lty=2,lwd=3,col="grey")
mash_table2$gene_flux_part1 <- ifelse(mash_table2$modified_mash_distance < 0.1666667 & mash_table2$pham_pham_dissimilarity > (mash_table2$modified_mash_distance * 3.5),TRUE,FALSE)
mash_table2$gene_flux_part2 <- ifelse(mash_table2$modified_mash_distance > 0.1666667 & mash_table2$pham_pham_dissimilarity > (mash_table2$modified_mash_distance * 2 + 0.25),TRUE,FALSE)
mash_table2$gene_flux_category <- ifelse(mash_table2$gene_flux_part1 == TRUE | mash_table2$gene_flux_part2 ==TRUE,"high","low")
mash_table2$gene_flux_category <- as.factor(mash_table2$gene_flux_category)




#Compute % of comparisons that are positioned in each similarity sector
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
#The last sector, "unsectored" is used to catch any that fall through the cracks
mash_table2$sector_intra_subcluster <- ifelse(mash_table2$modified_mash_distance < intra_subcluster_nucleotide_threshold_upper &
                                                mash_table2$modified_mash_distance >= 0 &
                                                mash_table2$pham_pham_dissimilarity < intra_subcluster_gene_content_threshold_upper &
                                                mash_table2$pham_pham_dissimilarity >= 0,TRUE,FALSE)

mash_table2$sector_intra_cluster <- ifelse(mash_table2$modified_mash_distance < intra_cluster_nucleotide_threshold_upper &
                                             mash_table2$pham_pham_dissimilarity < intra_cluster_gene_content_threshold_upper &
                                             (mash_table2$modified_mash_distance >= intra_subcluster_nucleotide_threshold_upper |
                                                (mash_table2$modified_mash_distance < intra_subcluster_nucleotide_threshold_upper &
                                                   mash_table2$pham_pham_dissimilarity >= intra_subcluster_gene_content_threshold_upper)),TRUE,FALSE)

mash_table2$sector_inter_cluster_distant_homology <- ifelse(mash_table2$modified_mash_distance < 0.5 &
                                                              mash_table2$modified_mash_distance >= intra_cluster_nucleotide_threshold_upper &
                                                              mash_table2$pham_pham_dissimilarity < intra_cluster_gene_content_threshold_upper &
                                                              mash_table2$pham_pham_dissimilarity >= 0,TRUE,FALSE)

mash_table2$sector_inter_cluster_hgt <- ifelse(mash_table2$modified_mash_distance < intra_cluster_nucleotide_threshold_upper &
                                                 mash_table2$modified_mash_distance >= 0 &
                                                 mash_table2$pham_pham_dissimilarity < 1 &
                                                 mash_table2$pham_pham_dissimilarity >= intra_cluster_gene_content_threshold_upper,TRUE,FALSE)

mash_table2$sector_no_similarity <- ifelse(mash_table2$modified_mash_distance < 0.5 &
                                             mash_table2$modified_mash_distance >= intra_cluster_nucleotide_threshold_upper &
                                             mash_table2$pham_pham_dissimilarity >= intra_cluster_gene_content_threshold_upper,TRUE,FALSE)

mash_table2$sector_non_cds_similarity <- ifelse(mash_table2$modified_mash_distance < intra_cluster_nucleotide_threshold_upper &
                                                  mash_table2$modified_mash_distance >= 0 &
                                                  mash_table2$pham_pham_dissimilarity == 1,TRUE,FALSE)

mash_table2$sector_filtered_out <- ifelse(mash_table2$modified_mash_distance == 0.5,TRUE,FALSE)

mash_table2$sector_unsectored <- ifelse(mash_table2$sector_intra_subcluster == FALSE &
                                          mash_table2$sector_intra_cluster == FALSE &
                                          mash_table2$sector_inter_cluster_distant_homology == FALSE &
                                          mash_table2$sector_inter_cluster_hgt == FALSE &
                                          mash_table2$sector_no_similarity == FALSE &
                                          mash_table2$sector_non_cds_similarity == FALSE &
                                          mash_table2$sector_filtered_out == FALSE,TRUE,FALSE)

#Subset based on assigned sector
sector_intra_subcluster <- subset(mash_table2,mash_table2$sector_intra_subcluster == TRUE)
sector_intra_cluster <- subset(mash_table2,mash_table2$sector_intra_cluster == TRUE)
sector_inter_cluster_distant_homology <- subset(mash_table2,mash_table2$sector_inter_cluster_distant_homology == TRUE)
sector_inter_cluster_hgt <- subset(mash_table2,mash_table2$sector_inter_cluster_hgt == TRUE)
sector_no_similarity <- subset(mash_table2,mash_table2$sector_no_similarity == TRUE)
sector_non_cds_similarity <- subset(mash_table2,mash_table2$sector_non_cds_similarity == TRUE)
sector_filtered_out <- subset(mash_table2,mash_table2$sector_filtered_out == TRUE)
sector_unsectored <- subset(mash_table2,mash_table2$sector_unsectored == TRUE)

#Check to verify that all comparisons have been sectored
nrow(mash_table2) -
  (nrow(sector_intra_subcluster) +
     nrow(sector_intra_cluster) +
     nrow(sector_inter_cluster_distant_homology) +
     nrow(sector_inter_cluster_hgt) +
     nrow(sector_no_similarity) +
     nrow(sector_non_cds_similarity) +
     nrow(sector_filtered_out))

nrow(sector_unsectored)

#Plot comparisons by sector  
par(mar=c(4,8,4,4))
plot(sector_intra_subcluster$modified_mash_distance,sector_intra_subcluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_intra_cluster$modified_mash_distance,sector_intra_cluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_inter_cluster_distant_homology$modified_mash_distance,sector_inter_cluster_distant_homology$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_inter_cluster_hgt$modified_mash_distance,sector_inter_cluster_hgt$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_no_similarity$modified_mash_distance,sector_no_similarity$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_non_cds_similarity$modified_mash_distance,sector_non_cds_similarity$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_filtered_out$modified_mash_distance,sector_filtered_out$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sector_unsectored$modified_mash_distance,sector_unsectored$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")











#Compare ref and query host data columns and convert columns to factor class
mash_table2$host_superkingdom_compare <- ifelse(mash_table2$ref_host_superkingdom==mash_table2$query_host_superkingdom,as.character(mash_table2$ref_host_superkingdom),"different")
mash_table2$host_phylum_compare <- ifelse(mash_table2$ref_host_phylum==mash_table2$query_host_phylum,as.character(mash_table2$ref_host_phylum),"different")
mash_table2$host_class_compare <- ifelse(mash_table2$ref_host_class==mash_table2$query_host_class,as.character(mash_table2$ref_host_class),"different")
mash_table2$host_order_compare <- ifelse(mash_table2$ref_host_order==mash_table2$query_host_order,as.character(mash_table2$ref_host_order),"different")
mash_table2$host_family_compare <- ifelse(mash_table2$ref_host_family==mash_table2$query_host_family,as.character(mash_table2$ref_host_family),"different")
mash_table2$host_genus_compare <- ifelse(mash_table2$ref_host_genus==mash_table2$query_host_genus,as.character(mash_table2$ref_host_genus),"different")
mash_table2$phage_superkingdom_compare <- ifelse(mash_table2$ref_phage_superkingdom==mash_table2$query_phage_superkingdom,as.character(mash_table2$ref_phage_superkingdom),"different")
mash_table2$phage_viral_type_compare <- ifelse(mash_table2$ref_phage_viral_type==mash_table2$query_phage_viral_type,as.character(mash_table2$ref_phage_viral_type),"different")
mash_table2$phage_order_compare <- ifelse(mash_table2$ref_phage_order==mash_table2$query_phage_order,as.character(mash_table2$ref_phage_order),"different")
mash_table2$phage_family_compare <- ifelse(mash_table2$ref_phage_family==mash_table2$query_phage_family,as.character(mash_table2$ref_phage_family),"different")
mash_table2$phage_genus_compare <- ifelse(mash_table2$ref_phage_genus==mash_table2$query_phage_genus,as.character(mash_table2$ref_phage_genus),"different")
mash_table2$phage_cluster_compare <- ifelse(mash_table2$ref_phage_cluster==mash_table2$query_phage_cluster,as.character(mash_table2$ref_phage_cluster),"different")
mash_table2$phage_subcluster_compare <- ifelse(mash_table2$ref_phage_subcluster==mash_table2$query_phage_subcluster,as.character(mash_table2$ref_phage_subcluster),"different")
mash_table2$phage_cluster_source_compare <- ifelse(mash_table2$ref_phage_cluster_source==mash_table2$query_phage_cluster_source,as.character(mash_table2$ref_phage_cluster_source),"different")
mash_table2$phage_temperate_compare <- ifelse(mash_table2$ref_phage_temperate==mash_table2$query_phage_temperate,as.character(mash_table2$ref_phage_temperate),"different")
mash_table2$phage_rnrdb_compare <- ifelse(mash_table2$ref_rnrdb==mash_table2$query_rnrdb,as.character(mash_table2$ref_rnrdb),"different")
mash_table2$phage_phage_rnr_compare <- ifelse(mash_table2$ref_phage_rnr==mash_table2$query_phage_rnr,as.character(mash_table2$ref_phage_rnr),"different")
mash_table2$phamerator_status_compare <- ifelse(mash_table2$ref_phamerator_status==mash_table2$query_phamerator_status,as.character(mash_table2$ref_phamerator_status),"different")
mash_table2$phage_toxic_compare <- ifelse(mash_table2$ref_toxic==mash_table2$query_toxic,as.character(mash_table2$ref_toxic),"different")
mash_table2$phage_extra_chrome_compare <- ifelse(mash_table2$ref_extra_chrome==mash_table2$query_extra_chrome,as.character(mash_table2$ref_extra_chrome),"different")
mash_table2$phage_predicted_temperate_compare <- ifelse(mash_table2$ref_predicted_temperate==mash_table2$query_predicted_temperate,as.character(mash_table2$ref_predicted_temperate),"different")





#Now set all new columns to factor class
mash_table2$host_superkingdom_compare <- as.factor(mash_table2$host_superkingdom_compare)
mash_table2$host_phylum_compare <- as.factor(mash_table2$host_phylum_compare)
mash_table2$host_class_compare <- as.factor(mash_table2$host_class_compare)
mash_table2$host_order_compare <- as.factor(mash_table2$host_order_compare)
mash_table2$host_family_compare <- as.factor(mash_table2$host_family_compare)
mash_table2$host_genus_compare <- as.factor(mash_table2$host_genus_compare)
mash_table2$phage_superkingdom_compare <- as.factor(mash_table2$phage_superkingdom_compare)
mash_table2$phage_viral_type_compare <- as.factor(mash_table2$phage_viral_type_compare)
mash_table2$phage_order_compare <- as.factor(mash_table2$phage_order_compare)
mash_table2$phage_family_compare <- as.factor(mash_table2$phage_family_compare)
mash_table2$phage_genus_compare <- as.factor(mash_table2$phage_genus_compare)
mash_table2$phage_cluster_compare <- as.factor(mash_table2$phage_cluster_compare)
mash_table2$phage_subcluster_compare <- as.factor(mash_table2$phage_subcluster_compare)
mash_table2$phage_cluster_source_compare <- as.factor(mash_table2$phage_cluster_source_compare)
mash_table2$phage_temperate_compare <- as.factor(mash_table2$phage_temperate_compare)
mash_table2$phage_rnrdb_compare <- as.factor(mash_table2$phage_rnrdb_compare)
mash_table2$phage_phage_rnr_compare <- as.factor(mash_table2$phage_phage_rnr_compare)
mash_table2$phamerator_status_compare <- as.factor(mash_table2$phamerator_status_compare)
mash_table2$phage_toxic_compare <- as.factor(mash_table2$phage_toxic_compare)
mash_table2$phage_extra_chrome_compare <- as.factor(mash_table2$phage_extra_chrome_compare)
mash_table2$phage_predicted_temperate_compare <- as.factor(mash_table2$phage_predicted_temperate_compare)









###At this point, all imported data has been processed. From here on out, commands are mostly analysis


bacteria_dsDNA <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria" & mash_table2$phage_viral_type_compare == "dsDNA")
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


















#If I want to delete certain columns to reduce size
# mash_table2 <- subset(mash_table2,select = -c(
#   pham_phage1,
#   pham_phage2,
#   pham_averaged_shared_proportion,
#   pham_jaccard_similarity,
#   ani_reference,
#   ani_query,
#   ani_ani,
#   query_header_source_info,
#   ref_header_source_info
# )
# )
#
#
#








#Plot histograms to see the range of genome size disparities
par(mar=c(4,8,4,4))
hist(mash_table2$size_diff_max_percent,breaks=1000,col="black",cex.axis=2,xlab=NULL,ylab=NULL,main=NULL,las=1)

par(mar=c(4,8,4,4))
hist(mash_table2$size_diff_max_percent,breaks=1000,col="black",cex.axis=2,xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,8))
abline(v=1,lty=2,lwd=3,col="grey")









#
#
#
###To subset based on the filtered results:
#mash_host_pham_pvalue_size <- subset(mash_host_pham,mash_host_pham$filter == TRUE)


#Subset based on pvalue
#mash_host_pham_pvalue <- subset(mash_host_pham,mash_host_pham$mash_pvalue < 1e-10)

#Subset based on genome size
#Size difference of 1 is equivalent to 100% genome size disparity. One of the two genomes is twice the size of the other genome.
#mash_host_pham_pvalue_size <- subset(mash_host_pham_pvalue,mash_host_pham_pvalue$size_diff_max_percent < 1)



#Compute percentage of comparisons that passed the filtering
passed_filter <- nrow(mash_host_pham_pvalue_size) / nrow(mash_table)
#
#
#













###Scatterplot
#xlab="Mash distance"
#ylab="Pham distance"

par(mar=c(4,8,4,4))
hist(mash_table2$modified_mash_distance,breaks=((range(mash_table2$modified_mash_distance)[2]-range(mash_table2$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,4e4),col="black",cex.axis=2)


par(mar=c(4,8,20,4))
hist(mash_table2$modified_mash_distance,breaks=((range(mash_table2$modified_mash_distance)[2]-range(mash_table2$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,4e4),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(mash_table2$modified_mash_distance,breaks=((range(mash_table2$modified_mash_distance)[2]-range(mash_table2$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),col="black",cex.axis=2)



par(mar=c(4,8,4,4))
hist(mash_table2$pham_pham_dissimilarity,breaks=((range(mash_table2$pham_pham_dissimilarity)[2]-range(mash_table2$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)



par(mar=c(4,8,4,4))
plot(mash_table2$modified_mash_distance,mash_table2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")














###Check Eukaryotic controls
euk_check <- subset(mash_table2,mash_table2$ref_host_superkingdom == "Eukaryota" | mash_table2$query_host_superkingdom == "Eukaryota")
euk_check <- subset(euk_check,euk_check$ref_host_superkingdom == "Bacteria" | euk_check$query_host_superkingdom == "Bacteria")

compute_sector_distribution(euk_check)

par(mar=c(4,8,4,4))
plot(euk_check$modified_mash_distance,euk_check$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Since there is such a narrow range of values, histograms don't work








###Check archaea controls
controls <- mash_table2
controls$archaea_one_or_two <- ifelse(controls$ref_host_superkingdom == "Archaea" | controls$query_host_superkingdom == "Archaea",TRUE,FALSE)
controls$archaea_both <- ifelse(controls$ref_host_superkingdom == "Archaea" & controls$query_host_superkingdom == "Archaea",TRUE,FALSE)
controls$archaea_one <- ifelse(controls$archaea_one_or_two == TRUE & controls$archaea_both == FALSE,TRUE,FALSE)
controls$bacteria_one_or_two <- ifelse(controls$ref_host_superkingdom == "Bacteria" | controls$query_host_superkingdom == "Bacteria",TRUE,FALSE)
controls$archaea_one_bacteria_one <- ifelse(controls$archaea_one == TRUE & controls$bacteria_one_or_two == TRUE,TRUE,FALSE)


archaea_check <- subset(controls,controls$archaea_one_bacteria_one == TRUE)
compute_sector_distribution(archaea_check)

par(mar=c(4,8,4,4))
plot(archaea_check$modified_mash_distance,archaea_check$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,20,4))
hist(archaea_check$modified_mash_distance,breaks=((range(archaea_check$modified_mash_distance)[2]-range(archaea_check$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,200),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(archaea_check$pham_pham_dissimilarity,breaks=((range(archaea_check$pham_pham_dissimilarity)[2]-range(archaea_check$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)















###By phage type
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_diff <- subset(mash_table2,mash_table2$phage_viral_type_compare == "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
type_dsRNA <- subset(type,type$phage_viral_type_compare == "dsRNA")
type_ssDNA <- subset(type,type$phage_viral_type_compare == "ssDNA")
type_ssRNA <- subset(type,type$phage_viral_type_compare == "ssRNA")


compute_sector_distribution(type_dsDNA)


#Compare similarity of different types
par(mar=c(4,8,4,4))
plot(type_diff$modified_mash_distance,type_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,20,4))
hist(type_diff$modified_mash_distance,breaks=((range(type_diff$modified_mash_distance)[2]-range(type_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,200),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(type_diff$pham_pham_dissimilarity,breaks=((range(type_diff$pham_pham_dissimilarity)[2]-range(type_diff$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)



#Now plot for each nucleic acid type
par(mar=c(4,8,4,4))
plot(type_dsDNA$modified_mash_distance,type_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(type_dsRNA$modified_mash_distance,type_dsRNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(type_ssDNA$modified_mash_distance,type_ssDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(type_ssRNA$modified_mash_distance,type_ssRNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")














###Compute distance of all Actinobacteriophage Singletons from other Actinobacteriophages
#Histogram = Mash distances for all Actino785 comparisons containing at least one Singleton
actino785_data <- subset(mash_table2,mash_table2$phage_cluster_source_compare == "actino")
actino785_data_same <- subset(actino785_data,actino785_data$phage_cluster_compare != "different")
actino785_data_diff <- subset(actino785_data,actino785_data$phage_cluster_compare == "different")
actino785_data_diff$query_singleton <- grepl("^Singleton",actino785_data_diff$query_phage_cluster)
actino785_data_diff$ref_singleton <- grepl("^Singleton",actino785_data_diff$ref_phage_cluster)
actino785_data_diff$singleton_one_or_two <- ifelse(actino785_data_diff$ref_singleton == TRUE | actino785_data_diff$query_singleton == TRUE,TRUE,FALSE)
#actino785_data_diff$singleton_neither <- ifelse(actino785_data_diff$ref_singleton == TRUE | actino785_data_diff$query_singleton == TRUE,FALSE,TRUE)
#actino785_data_diff$singleton_both <- ifelse(actino785_data_diff$ref_singleton == TRUE & actino785_data_diff$query_singleton == TRUE,TRUE,FALSE)
#actino785_data_diff$singleton_one <- ifelse(actino785_data_diff$singleton_both == TRUE | actino785_data_diff$singleton_neither == TRUE,FALSE,TRUE)

actino785_singleton_one_or_two <- subset(actino785_data_diff,actino785_data_diff$singleton_one_or_two == TRUE)
#actino785_singleton_one <- subset(actino785_data_diff,actino785_data_diff$singleton_one == TRUE)
#actino785_singleton_both <- subset(actino785_data_diff,actino785_data_diff$singleton_both == TRUE)
#actino785_singleton_neither <- subset(actino785_data_diff,actino785_data_diff$singleton_neither == TRUE)

compute_sector_distribution(actino785_singleton_one_or_two)

par(mar=c(4,8,4,4))
plot(actino785_singleton_one_or_two$modified_mash_distance,actino785_singleton_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,20,4))
hist(actino785_singleton_one_or_two$modified_mash_distance,breaks=((range(actino785_singleton_one_or_two$modified_mash_distance)[2]-range(actino785_singleton_one_or_two$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(actino785_singleton_one_or_two$pham_pham_dissimilarity,breaks=((range(actino785_singleton_one_or_two$pham_pham_dissimilarity)[2]-range(actino785_singleton_one_or_two$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)












###Compare distances for phages that are not subclustered.
#How similar are they to all other phages?
actino785_data <- subset(mash_table2,mash_table2$phage_cluster_source_compare == "actino")
actino785_data$query_subclustered <- ifelse(is.na(actino785_data$query_phage_subcluster) == TRUE,FALSE,TRUE) 
actino785_data$ref_subclustered <- ifelse(is.na(actino785_data$ref_phage_subcluster) == TRUE,FALSE,TRUE) 
actino785_data$subclustered_one_or_two <- ifelse(actino785_data$query_subclustered == TRUE | actino785_data$ref_subclustered == TRUE,TRUE,FALSE)
actino785_data$subclustered_neither <- ifelse(actino785_data$query_subclustered == FALSE & actino785_data$ref_subclustered == FALSE,TRUE,FALSE)
actino785_data$subclustered_both <- ifelse(actino785_data$query_subclustered == TRUE & actino785_data$ref_subclustered == TRUE,TRUE,FALSE)
actino785_data$subclustered_one <- ifelse(actino785_data$subclustered_one_or_two == TRUE & actino785_data$subclustered_both == FALSE,TRUE,FALSE)

actino785_data_same_cluster <- subset(actino785_data,actino785_data$phage_cluster_compare != "different")
actino785_data_diff_cluster <- subset(actino785_data,actino785_data$phage_cluster_compare == "different")
actino785_data_same_cluster_neither_subclustered <- subset(actino785_data_same_cluster,actino785_data_same_cluster$subclustered_neither == TRUE)
actino785_data_diff_cluster_one_subclustered <- subset(actino785_data_diff_cluster,actino785_data_diff_cluster$subclustered_one == TRUE)

compute_sector_distribution(actino785_data_same_cluster_neither_subclustered)

par(mar=c(4,8,4,4))
plot(actino785_data_same_cluster_neither_subclustered$modified_mash_distance,actino785_data_same_cluster_neither_subclustered$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino785_data_diff_cluster_one_subclustered$modified_mash_distance,actino785_data_diff_cluster_one_subclustered$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")
















###By phage family and tail type

bacteriophages <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria")
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







par(mar=c(4,8,4,4))
plot(caudovirales_family_diff$modified_mash_distance,caudovirales_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,15,4))
hist(caudovirales_family_diff$modified_mash_distance,breaks=((range(caudovirales_family_diff$modified_mash_distance)[2]-range(caudovirales_family_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)


par(mar=c(4,4,15,4))
hist(caudovirales_family_diff$pham_pham_dissimilarity,breaks=((range(caudovirales_family_diff$pham_pham_dissimilarity)[2]-range(caudovirales_family_diff$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)




par(mar=c(4,8,4,4))
plot(myo$modified_mash_distance,myo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sipho$modified_mash_distance,sipho$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(podo$modified_mash_distance,podo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





par(mar=c(4,8,4,4))
plot(caudovirales_family_diff_sipho_myo$modified_mash_distance,caudovirales_family_diff_sipho_myo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(caudovirales_family_diff_sipho_podo$modified_mash_distance,caudovirales_family_diff_sipho_podo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(caudovirales_family_diff_myo_podo$modified_mash_distance,caudovirales_family_diff_myo_podo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")








par(mar=c(4,8,4,4))
plot(podo_temperate$modified_mash_distance,podo_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(podo_lytic$modified_mash_distance,podo_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(podo_lifestyle_diff$modified_mash_distance,podo_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




par(mar=c(4,8,4,4))
plot(sipho_temperate$modified_mash_distance,sipho_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sipho_lytic$modified_mash_distance,sipho_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(sipho_lifestyle_diff$modified_mash_distance,sipho_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





par(mar=c(4,8,4,4))
plot(myo_temperate$modified_mash_distance,myo_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(myo_lytic$modified_mash_distance,myo_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(myo_lifestyle_diff$modified_mash_distance,myo_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






#
#
#
###By host superkingdom
superkingdom_bacteria <- subset(mash_host_pvalue,mash_host_pvalue$host_superkingdom_compare == "Bacteria")
hist(superkingdom_bacteria$mash_distance,breaks=((range(superkingdom_bacteria$mash_distance)[2]-range(superkingdom_bacteria$mash_distance)[1]) * 100),xlab=NULL,ylab="Number of comparisons",main="Merged2333 s25000k15size superkingdom bacteria",xlim=c(0,0.6),col="black")

###By phage order
phage_order <- subset(kmer14_host_pham_pvalue,kmer14_host_pham_pvalue$phage_order_compare != "different")
phage_order_caudo <- subset(phage_order,phage_order$phage_order_compare == "Caudovirales")
phage_order_diff <- subset(kmer14_host_pham_pvalue,kmer14_host_pham_pvalue$phage_order_compare == "different")

plot(phage_order_caudo$kmer14_distance,phage_order_caudo$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k10p1e10 Caudovirales",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.00001)
abline(0,2.3)
#
#
#








###By host phylum, then by lifestyle
#Since I'm mostly interested in dsDNA phage evolution, remove other phage types.
#This removes some background in the plot, particularly with Proteobacteria phages.

bacteria <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria")
type_dsDNA <- subset(bacteria,bacteria$phage_viral_type_compare == "dsDNA")

host_phylum <- subset(type_dsDNA,type_dsDNA$host_phylum_compare != "different")
host_phylum_diff <- subset(type_dsDNA,type_dsDNA$host_phylum_compare == "different")

actino <- subset(host_phylum,host_phylum$host_phylum_compare == "Actinobacteria")
bacter <- subset(host_phylum,host_phylum$host_phylum_compare == "Bacteroidetes")
cyano <- subset(host_phylum,host_phylum$host_phylum_compare == "Cyanobacteria")
firm <- subset(host_phylum,host_phylum$host_phylum_compare == "Firmicutes")
proteo <- subset(host_phylum,host_phylum$host_phylum_compare == "Proteobacteria")



compute_sector_distribution(type_dsDNA)
compute_sector_distribution(host_phylum_diff)






#Scatter plots of Mash vs Pham distances by host phyla


#all dsDNA phages
par(mar=c(4,8,4,4))
plot(type_dsDNA$modified_mash_distance,type_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(actino$modified_mash_distance,actino$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(bacter$modified_mash_distance,bacter$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano$modified_mash_distance,cyano$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm$modified_mash_distance,firm$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo$modified_mash_distance,proteo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(host_phylum_diff$modified_mash_distance,host_phylum_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


#Histogram distributions of mash distances


#par(mar=c(4,12,15,0))
#hist(type_dsDNA$modified_mash_distance,breaks=((range(type_dsDNA$modified_mash_distance)[2]-range(type_dsDNA$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,xlim=c(0,0.5),ylim=c(0,3e4),col="black",yaxt="n")
#axis(side=2,at=c(0,15000,30000),las=1,cex.axis=4)

#truncated
par(mar=c(4,8,15,4))
hist(type_dsDNA$modified_mash_distance,breaks=((range(type_dsDNA$modified_mash_distance)[2]-range(type_dsDNA$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,3e4),col="black",cex.axis=2)

par(mar=c(4,8,15,4))
hist(actino$modified_mash_distance,breaks=((range(actino$modified_mash_distance)[2]-range(actino$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,3e4),col="black",cex.axis=2)

par(mar=c(4,8,15,4))
hist(bacter$modified_mash_distance,breaks=((range(bacter$modified_mash_distance)[2]-range(bacter$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,50),col="black",cex.axis=2)

par(mar=c(4,8,15,4))
hist(cyano$modified_mash_distance,breaks=((range(cyano$modified_mash_distance)[2]-range(cyano$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,100),col="black",cex.axis=2)

par(mar=c(4,8,15,4))
hist(firm$modified_mash_distance,breaks=((range(firm$modified_mash_distance)[2]-range(firm$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)

par(mar=c(4,8,15,4))
hist(proteo$modified_mash_distance,breaks=((range(proteo$modified_mash_distance)[2]-range(proteo$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)





#Histogram distribution of pham distances


par(mar=c(4,4,15,4))
hist(type_dsDNA$pham_pham_dissimilarity,breaks=((range(type_dsDNA$pham_pham_dissimilarity)[2]-range(type_dsDNA$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)



par(mar=c(4,4,15,4))
hist(actino$pham_pham_dissimilarity,breaks=((range(actino$pham_pham_dissimilarity)[2]-range(actino$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,4,15,4))
hist(bacter$pham_pham_dissimilarity,breaks=((range(bacter$pham_pham_dissimilarity)[2]-range(bacter$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,500),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,4,15,4))
hist(cyano$pham_pham_dissimilarity,breaks=((range(cyano$pham_pham_dissimilarity)[2]-range(cyano$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,500),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,4,15,4))
hist(firm$pham_pham_dissimilarity,breaks=((range(firm$pham_pham_dissimilarity)[2]-range(firm$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,4,15,4))
hist(proteo$pham_pham_dissimilarity,breaks=((range(proteo$pham_pham_dissimilarity)[2]-range(proteo$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)










###Phylum then lifestyle
actino_temperate <- subset(actino,actino$phage_temperate_compare == "yes")
actino_lytic <- subset(actino,actino$phage_temperate_compare == "no")
actino_lytic_temperate <- subset(actino,actino$phage_temperate_compare == "different")
actino_lifestyle_unspecified <- subset(actino,is.na(actino$phage_temperate_compare))

bacter_temperate <- subset(bacter,bacter$phage_temperate_compare == "yes")
bacter_lytic <- subset(bacter,bacter$phage_temperate_compare == "no")
bacter_lytic_temperate <- subset(bacter,bacter$phage_temperate_compare == "different")
bacter_lifestyle_unspecified <- subset(bacter,is.na(bacter$phage_temperate_compare))



cyano_temperate <- subset(cyano,cyano$phage_temperate_compare == "yes")
cyano_lytic <- subset(cyano,cyano$phage_temperate_compare == "no")
cyano_lytic_temperate <- subset(cyano,cyano$phage_temperate_compare == "different")
cyano_lifestyle_unspecified <- subset(cyano,is.na(cyano$phage_temperate_compare))



firm_temperate <- subset(firm,firm$phage_temperate_compare == "yes")
firm_lytic <- subset(firm,firm$phage_temperate_compare == "no")
firm_lytic_temperate <- subset(firm,firm$phage_temperate_compare == "different")
firm_lifestyle_unspecified <- subset(firm,is.na(firm$phage_temperate_compare))


proteo_temperate <- subset(proteo,proteo$phage_temperate_compare == "yes")
proteo_lytic <- subset(proteo,proteo$phage_temperate_compare == "no")
proteo_lytic_temperate <- subset(proteo,proteo$phage_temperate_compare == "different")
proteo_lifestyle_unspecified <- subset(proteo,is.na(proteo$phage_temperate_compare))

#Check to verify all data was subsetted correctly...
nrow(actino) - (nrow(actino_temperate) +
                  nrow(actino_lytic) +
                  nrow(actino_lytic_temperate) +
                  nrow(actino_lifestyle_unspecified))

nrow(bacter) - (nrow(bacter_temperate) +
                  nrow(bacter_lytic) +
                  nrow(bacter_lytic_temperate) +
                  nrow(bacter_lifestyle_unspecified))

nrow(cyano) - (nrow(cyano_temperate) +
                 nrow(cyano_lytic) +
                 nrow(cyano_lytic_temperate) +
                 nrow(cyano_lifestyle_unspecified))

nrow(firm) - (nrow(firm_temperate) +
                nrow(firm_lytic) +
                nrow(firm_lytic_temperate) +
                nrow(firm_lifestyle_unspecified))

nrow(proteo) - (nrow(proteo_temperate) +
                  nrow(proteo_lytic) +
                  nrow(proteo_lytic_temperate) +
                  nrow(proteo_lifestyle_unspecified))






#Separate, uncolored scatters
par(mar=c(4,8,4,4))
plot(actino_temperate$modified_mash_distance,actino_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_lytic$modified_mash_distance,actino_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(bacter_temperate$modified_mash_distance,bacter_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(bacter_lytic$modified_mash_distance,bacter_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_temperate$modified_mash_distance,cyano_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_lytic$modified_mash_distance,cyano_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_temperate$modified_mash_distance,firm_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_lytic$modified_mash_distance,firm_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_temperate$modified_mash_distance,proteo_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_lytic$modified_mash_distance,proteo_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




#Combined, colored scatters
par(mar=c(4,8,4,4))
plot(actino_lifestyle_unspecified$modified_mash_distance,actino_lifestyle_unspecified$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(actino_temperate$modified_mash_distance,actino_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_lytic$modified_mash_distance,actino_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(actino_lytic_temperate$modified_mash_distance,actino_lytic_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(bacter_lifestyle_unspecified$modified_mash_distance,bacter_lifestyle_unspecified$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(bacter_temperate$modified_mash_distance,bacter_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(bacter_lytic$modified_mash_distance,bacter_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(bacter_lytic_temperate$modified_mash_distance,bacter_lytic_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_lifestyle_unspecified$modified_mash_distance,cyano_lifestyle_unspecified$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(cyano_temperate$modified_mash_distance,cyano_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(cyano_lytic$modified_mash_distance,cyano_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(cyano_lytic_temperate$modified_mash_distance,cyano_lytic_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(firm_lifestyle_unspecified$modified_mash_distance,firm_lifestyle_unspecified$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(firm_temperate$modified_mash_distance,firm_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(firm_lytic$modified_mash_distance,firm_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(firm_lytic_temperate$modified_mash_distance,firm_lytic_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_lifestyle_unspecified$modified_mash_distance,proteo_lifestyle_unspecified$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(proteo_temperate$modified_mash_distance,proteo_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(proteo_lytic$modified_mash_distance,proteo_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(proteo_lytic_temperate$modified_mash_distance,proteo_lytic_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,2,lty=2,lwd=3,col="grey")



















###By host phylum, then by taxonomy ranks
#Should be dsDNA phages, then same phylum, then subset from there.
#For comparing between different rankings, do a two-step subsetting.
#First, check if the first level is different, then check if the next level up is the same


actino_class_same <- subset(actino,actino$host_class_compare != "different")
actino_class_same <- subset(actino_class_same,actino_class_same$host_phylum_compare != "different")

actino_order_same <- subset(actino,actino$host_order_compare != "different")
actino_order_same <- subset(actino_order_same,actino_order_same$host_class_compare != "different")

actino_family_same <- subset(actino,actino$host_family_compare != "different")
actino_family_same <- subset(actino_family_same,actino_family_same$host_order_compare != "different")

actino_genus_same <- subset(actino,actino$host_genus_compare != "different")
actino_genus_same <- subset(actino_genus_same,actino_genus_same$host_family_compare != "different")


# bacter_class_same <- subset(bacter,bacter$host_class_compare != "different")
# bacter_order_same <- subset(bacter,bacter$host_order_compare != "different")
# bacter_family_same <- subset(bacter,bacter$host_family_compare != "different")
# bacter_genus_same <- subset(bacter,bacter$host_genus_compare != "different")



cyano_class_same <- subset(cyano,cyano$host_class_compare != "different")
cyano_class_same <- subset(cyano_class_same,cyano_class_same$host_phylum_compare != "different")

cyano_order_same <- subset(cyano,cyano$host_order_compare != "different")
cyano_order_same <- subset(cyano_order_same,cyano_order_same$host_class_compare != "different")

cyano_family_same <- subset(cyano,cyano$host_family_compare != "different")
cyano_family_same <- subset(cyano_family_same,cyano_family_same$host_order_compare != "different")

cyano_genus_same <- subset(cyano,cyano$host_genus_compare != "different")
cyano_genus_same <- subset(cyano_genus_same,cyano_genus_same$host_family_compare != "different")



firm_class_same <- subset(firm,firm$host_class_compare != "different")
firm_class_same <- subset(firm_class_same,firm_class_same$host_phylum_compare != "different")

firm_order_same <- subset(firm,firm$host_order_compare != "different")
firm_order_same <- subset(firm_order_same,firm_order_same$host_class_compare != "different")

firm_family_same <- subset(firm,firm$host_family_compare != "different")
firm_family_same <- subset(firm_family_same,firm_family_same$host_order_compare != "different")

firm_genus_same <- subset(firm,firm$host_genus_compare != "different")
firm_genus_same <- subset(firm_genus_same,firm_genus_same$host_family_compare != "different")


proteo_class_same <- subset(proteo,proteo$host_class_compare != "different")
proteo_class_same <- subset(proteo_class_same,proteo_class_same$host_phylum_compare != "different")

proteo_order_same <- subset(proteo,proteo$host_order_compare != "different")
proteo_order_same <- subset(proteo_order_same,proteo_order_same$host_class_compare != "different")

proteo_family_same <- subset(proteo,proteo$host_family_compare != "different")
proteo_family_same <- subset(proteo_family_same,proteo_family_same$host_order_compare != "different")

proteo_genus_same <- subset(proteo,proteo$host_genus_compare != "different")
proteo_genus_same <- subset(proteo_genus_same,proteo_genus_same$host_family_compare != "different")







actino_class_diff <- subset(actino,actino$host_class_compare == "different")
actino_class_diff <- subset(actino_class_diff,actino_class_diff$host_phylum_compare != "different")

actino_order_diff <- subset(actino,actino$host_order_compare == "different")
actino_order_diff <- subset(actino_order_diff,actino_order_diff$host_class_compare != "different")

actino_family_diff <- subset(actino,actino$host_family_compare == "different")
actino_family_diff <- subset(actino_family_diff,actino_family_diff$host_order_compare != "different")

actino_genus_diff <- subset(actino,actino$host_genus_compare == "different")
actino_genus_diff <- subset(actino_genus_diff,actino_genus_diff$host_family_compare != "different")



# bacter_class_diff <- subset(bacter,bacter$host_class_compare == "different")
# bacter_class_diff <- subset(bacter_class_diff,bacter_class_diff$host_phylum_compare != "different")
# 
# bacter_order_diff <- subset(bacter,bacter$host_order_compare == "different")
# bacter_order_diff <- subset(bacter_order_diff,bacter_order_diff$host_class_compare != "different")
# 
# bacter_family_diff <- subset(bacter,bacter$host_family_compare == "different")
# bacter_family_diff <- subset(bacter_family_diff,bacter_family_diff$host_order_compare != "different")
# 
# 
# bacter_genus_diff <- subset(bacter,bacter$host_genus_compare == "different")
# bacter_genus_diff <- subset(bacter_genus_diff,bacter_genus_diff$host_family_compare != "different")




cyano_class_diff <- subset(cyano,cyano$host_class_compare == "different")
cyano_class_diff <- subset(cyano_class_diff,cyano_class_diff$host_phylum_compare != "different")

cyano_order_diff <- subset(cyano,cyano$host_order_compare == "different")
cyano_order_diff <- subset(cyano_order_diff,cyano_order_diff$host_class_compare != "different")

cyano_family_diff <- subset(cyano,cyano$host_family_compare == "different")
cyano_family_diff <- subset(cyano_family_diff,cyano_family_diff$host_order_compare != "different")

cyano_genus_diff <- subset(cyano,cyano$host_genus_compare == "different")
cyano_genus_diff <- subset(cyano_genus_diff,cyano_genus_diff$host_family_compare != "different")


firm_class_diff <- subset(firm,firm$host_class_compare == "different")
firm_class_diff <- subset(firm_class_diff,firm_class_diff$host_phylum_compare != "different")

firm_order_diff <- subset(firm,firm$host_order_compare == "different")
firm_order_diff <- subset(firm_order_diff,firm_order_diff$host_class_compare != "different")

firm_family_diff <- subset(firm,firm$host_family_compare == "different")
firm_family_diff <- subset(firm_family_diff,firm_family_diff$host_order_compare != "different")

firm_genus_diff <- subset(firm,firm$host_genus_compare == "different")
firm_genus_diff <- subset(firm_genus_diff,firm_genus_diff$host_family_compare != "different")






proteo_class_diff <- subset(proteo,proteo$host_class_compare == "different")
proteo_class_diff <- subset(proteo_class_diff,proteo_class_diff$host_phylum_compare != "different")


proteo_order_diff <- subset(proteo,proteo$host_order_compare == "different")
proteo_order_diff <- subset(proteo_order_diff,proteo_order_diff$host_class_compare != "different")



proteo_family_diff <- subset(proteo,proteo$host_family_compare == "different")
proteo_family_diff <- subset(proteo_family_diff,proteo_family_diff$host_order_compare != "different")


proteo_genus_diff <- subset(proteo,proteo$host_genus_compare == "different")
proteo_genus_diff <- subset(proteo_genus_diff,proteo_genus_diff$host_family_compare != "different")




actino_class_same
actino_order_same
actino_family_same
actino_genus_same



par(mar=c(4,8,4,4))
plot(actino_class_same$modified_mash_distance,actino_class_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_order_same$modified_mash_distance,actino_order_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_family_same$modified_mash_distance,actino_family_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_genus_same$modified_mash_distance,actino_genus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




# bacter_class_same
# bacter_order_same
# bacter_family_same
# bacter_genus_same
# 
# par(mar=c(4,8,4,4))
# plot(bacter_class_same$modified_mash_distance,bacter_class_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacter_order_same$modified_mash_distance,bacter_order_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacter_family_same$modified_mash_distance,bacter_family_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacter_genus_same$modified_mash_distance,bacter_genus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")





cyano_class_same
cyano_order_same
cyano_family_same
cyano_genus_same

par(mar=c(4,8,4,4))
plot(cyano_class_same$modified_mash_distance,cyano_class_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_order_same$modified_mash_distance,cyano_order_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_family_same$modified_mash_distance,cyano_family_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_genus_same$modified_mash_distance,cyano_genus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






firm_class_same
firm_order_same
firm_family_same
firm_genus_same

par(mar=c(4,8,4,4))
plot(firm_class_same$modified_mash_distance,firm_class_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_order_same$modified_mash_distance,firm_order_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_family_same$modified_mash_distance,firm_family_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_genus_same$modified_mash_distance,firm_genus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






proteo_class_same
proteo_order_same
proteo_family_same
proteo_genus_same

par(mar=c(4,8,4,4))
plot(proteo_class_same$modified_mash_distance,proteo_class_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_order_same$modified_mash_distance,proteo_order_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_family_same$modified_mash_distance,proteo_family_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_genus_same$modified_mash_distance,proteo_genus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






actino_class_diff
actino_order_diff
actino_family_diff
actino_genus_diff

par(mar=c(4,8,4,4))
plot(actino_class_diff$modified_mash_distance,actino_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_order_diff$modified_mash_distance,actino_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_family_diff$modified_mash_distance,actino_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_genus_diff$modified_mash_distance,actino_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





# bacter_class_diff
# bacter_order_diff
# bacter_family_diff
# bacter_genus_diff
# 
# par(mar=c(4,8,4,4))
# plot(bacter_class_diff$modified_mash_distance,bacter_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacter_order_diff$modified_mash_distance,bacter_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacter_family_diff$modified_mash_distance,bacter_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacter_genus_diff$modified_mash_distance,bacter_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")


cyano_class_diff
cyano_order_diff
cyano_family_diff
cyano_genus_diff

par(mar=c(4,8,4,4))
plot(cyano_class_diff$modified_mash_distance,cyano_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_order_diff$modified_mash_distance,cyano_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_family_diff$modified_mash_distance,cyano_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cyano_genus_diff$modified_mash_distance,cyano_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


firm_class_diff
firm_order_diff
firm_family_diff
firm_genus_diff

par(mar=c(4,8,4,4))
plot(firm_class_diff$modified_mash_distance,firm_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_order_diff$modified_mash_distance,firm_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_family_diff$modified_mash_distance,firm_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(firm_genus_diff$modified_mash_distance,firm_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

proteo_class_diff
proteo_order_diff
proteo_family_diff
proteo_genus_diff

par(mar=c(4,8,4,4))
plot(proteo_class_diff$modified_mash_distance,proteo_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_order_diff$modified_mash_distance,proteo_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_family_diff$modified_mash_distance,proteo_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(proteo_genus_diff$modified_mash_distance,proteo_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")














###Check by lifestyle.
#First select only dsDNA phages
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")

all_temperate <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "yes")
all_lytic <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "no")
all_different <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "different")

compute_sector_distribution(all_temperate)
compute_sector_distribution(all_lytic)
compute_sector_distribution(all_different)


par(mar=c(4,8,4,4))
plot(all_temperate$modified_mash_distance,all_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(all_lytic$modified_mash_distance,all_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(all_different$modified_mash_distance,all_different$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





par(mar=c(4,8,4,4))
hist(all_temperate$modified_mash_distance,breaks=((range(all_temperate$modified_mash_distance)[2]-range(all_temperate$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(all_lytic$modified_mash_distance,breaks=((range(all_lytic$modified_mash_distance)[2]-range(all_lytic$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(all_different$modified_mash_distance,breaks=((range(all_different$modified_mash_distance)[2]-range(all_different$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,1e4),col="black",cex.axis=2)







par(mar=c(4,8,4,4))
hist(all_temperate$pham_pham_dissimilarity,breaks=((range(all_temperate$pham_pham_dissimilarity)[2]-range(all_temperate$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,4,4))
hist(all_lytic$pham_pham_dissimilarity,breaks=((range(all_lytic$pham_pham_dissimilarity)[2]-range(all_lytic$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,4,4))
hist(all_different$pham_pham_dissimilarity,breaks=((range(all_different$pham_pham_dissimilarity)[2]-range(all_different$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)







#
#
#Temperate Actino not A
all_temperate_actino <- subset(all_temperate,all_temperate$host_phylum_compare == "Actinobacteria")
all_temperate_actino_a <- subset(all_temperate_actino,all_temperate_actino$phage_cluster_compare == "A")
all_temperate_actino_not_a <- subset(all_temperate_actino,all_temperate_actino$phage_cluster_compare != "A")

par(mar=c(4,8,4,4))
plot(all_temperate_actino_a$modified_mash_distance,all_temperate_actino_a$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(all_temperate_actino_not_a$modified_mash_distance,all_temperate_actino_not_a$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")
#
#
#




#
#
#Method 2
#Alternative to original host tax analysis
#Account for NA rankings
#Only data with no missing values across all ranks is used. So each ranking is nested within the above rank.
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
mash_table2_hosttax_filtered <- subset(type_dsDNA,is.na(type_dsDNA$host_superkingdom_compare) == FALSE)
mash_table2_hosttax_filtered <- subset(mash_table2_hosttax_filtered,is.na(mash_table2_hosttax_filtered$host_phylum_compare) == FALSE)
mash_table2_hosttax_filtered <- subset(mash_table2_hosttax_filtered,is.na(mash_table2_hosttax_filtered$host_class_compare) == FALSE)
mash_table2_hosttax_filtered <- subset(mash_table2_hosttax_filtered,is.na(mash_table2_hosttax_filtered$host_order_compare) == FALSE)
mash_table2_hosttax_filtered <- subset(mash_table2_hosttax_filtered,is.na(mash_table2_hosttax_filtered$host_family_compare) == FALSE)
mash_table2_hosttax_filtered <- subset(mash_table2_hosttax_filtered,is.na(mash_table2_hosttax_filtered$host_genus_compare) == FALSE)



host_superkingdom_same <- subset(mash_table2_hosttax_filtered,mash_table2_hosttax_filtered$host_superkingdom_compare != "different")

#option 1 = to compare using all host phyla
host_superkingdom_diff <- subset(mash_table2_hosttax_filtered,mash_table2_hosttax_filtered$host_superkingdom_compare == "different")
host_superkingdom_diff_temperate <- subset(host_superkingdom_diff,host_superkingdom_diff$phage_temperate_compare == "yes")
host_superkingdom_diff_lytic <- subset(host_superkingdom_diff,host_superkingdom_diff$phage_temperate_compare == "no")


host_phylum_same <- subset(host_superkingdom_same,host_superkingdom_same$host_phylum_compare != "different")
host_phylum_diff <- subset(host_superkingdom_same,host_superkingdom_same$host_phylum_compare == "different")
host_phylum_diff_temperate <- subset(host_phylum_diff,host_phylum_diff$phage_temperate_compare == "yes")
host_phylum_diff_lytic <- subset(host_phylum_diff,host_phylum_diff$phage_temperate_compare == "no")
#



#option 2 = to compare using one particular phylum
host_phylum_same <- subset(host_superkingdom_same,host_superkingdom_same$host_phylum_compare == "Proteobacteria")
#




host_class_same <- subset(host_phylum_same,host_phylum_same$host_class_compare != "different")
host_class_diff <- subset(host_phylum_same,host_phylum_same$host_class_compare == "different")
host_class_diff_temperate <- subset(host_class_diff,host_class_diff$phage_temperate_compare == "yes")
host_class_diff_lytic <- subset(host_class_diff,host_class_diff$phage_temperate_compare == "no")
host_class_diff_lifestyle_diff <- subset(host_class_diff,host_class_diff$phage_temperate_compare == "different")

host_order_same <- subset(host_class_same,host_class_same$host_order_compare != "different")
host_order_diff <- subset(host_class_same,host_class_same$host_order_compare == "different")
host_order_diff_temperate <- subset(host_order_diff,host_order_diff$phage_temperate_compare == "yes")
host_order_diff_lytic <- subset(host_order_diff,host_order_diff$phage_temperate_compare == "no")
host_order_diff_lifestyle_diff <- subset(host_order_diff,host_order_diff$phage_temperate_compare == "different")

host_family_same <- subset(host_order_same,host_order_same$host_family_compare != "different")
host_family_diff <- subset(host_order_same,host_order_same$host_family_compare == "different")
host_family_diff_temperate <- subset(host_family_diff,host_family_diff$phage_temperate_compare == "yes")
host_family_diff_lytic <- subset(host_family_diff,host_family_diff$phage_temperate_compare == "no")
host_family_diff_lifestyle_diff <- subset(host_family_diff,host_family_diff$phage_temperate_compare == "different")



host_genus_same <- subset(host_family_same,host_family_same$host_genus_compare != "different")
host_genus_diff <- subset(host_family_same,host_family_same$host_genus_compare == "different")
host_genus_diff_temperate <- subset(host_genus_diff,host_genus_diff$phage_temperate_compare == "yes")
host_genus_diff_lytic <- subset(host_genus_diff,host_genus_diff$phage_temperate_compare == "no")
host_genus_diff_lifestyle_diff <- subset(host_genus_diff,host_genus_diff$phage_temperate_compare == "different")




#Scatter plot of Mash vs Pham distributions
par(mar=c(4,8,4,4))
plot(host_superkingdom_same$modified_mash_distance,host_superkingdom_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_superkingdom_diff$modified_mash_distance,host_superkingdom_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_phylum_same$modified_mash_distance,host_phylum_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_phylum_diff$modified_mash_distance,host_phylum_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_class_same$modified_mash_distance,host_class_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_class_diff$modified_mash_distance,host_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_order_same$modified_mash_distance,host_order_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_order_diff$modified_mash_distance,host_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_family_same$modified_mash_distance,host_family_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_family_diff$modified_mash_distance,host_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_genus_same$modified_mash_distance,host_genus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_genus_diff$modified_mash_distance,host_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")








#Scatter plots of lifestyle by different host rankings

par(mar=c(4,8,4,4))
plot(host_class_diff_temperate$modified_mash_distance,host_class_diff_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_class_diff_lytic$modified_mash_distance,host_class_diff_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_class_diff_lifestyle_diff$modified_mash_distance,host_class_diff_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




par(mar=c(4,8,4,4))
plot(host_order_diff_temperate$modified_mash_distance,host_order_diff_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_order_diff_lytic$modified_mash_distance,host_order_diff_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_order_diff_lifestyle_diff$modified_mash_distance,host_order_diff_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



par(mar=c(4,8,4,4))
plot(host_family_diff_temperate$modified_mash_distance,host_family_diff_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_family_diff_lytic$modified_mash_distance,host_family_diff_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_family_diff_lifestyle_diff$modified_mash_distance,host_family_diff_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



par(mar=c(4,8,4,4))
plot(host_genus_diff_temperate$modified_mash_distance,host_genus_diff_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_genus_diff_lytic$modified_mash_distance,host_genus_diff_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_genus_diff_lifestyle_diff$modified_mash_distance,host_genus_diff_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")























#Method 1
#Subset by taxonomy. For each diff taxonomy ranking, it's a two-step subsetting: first by the rank you're interested in.
#First subset by dsDNA phage type, then by each host taxonomic ranking
#Then, once you have a list of comparisons from different rank, subset by the next higher rank, and make sure they're within the same higher rank
#Not all phages are ranked at every taxonomic level, so forcing each comparison to have the same ranking at the level higher than the level of comparison will lead to data loss,
#but for all diff comparisons, you want to make sure you're comparing phages of different ranking, but nested within the same next hierarchical ranking.

#By host superkingdom
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")


host_superkingdom <- subset(type_dsDNA,type_dsDNA$host_superkingdom_compare != "different")
host_superkingdom_diff <- subset(type_dsDNA,type_dsDNA$host_superkingdom_compare == "different")


#By host phylum
host_phylum <- subset(type_dsDNA,type_dsDNA$host_phylum_compare != "different")
host_phylum <- subset(host_phylum,host_phylum$host_superkingdom_compare != "different")

host_phylum_diff <- subset(type_dsDNA,type_dsDNA$host_phylum_compare == "different")
host_phylum_diff <- subset(host_phylum_diff,host_phylum_diff$host_superkingdom_compare != "different")


#By host class
host_class <- subset(type_dsDNA,type_dsDNA$host_class_compare != "different")
host_class <- subset(host_class,host_class$host_phylum_compare != "different")

host_class_diff <- subset(type_dsDNA,type_dsDNA$host_class_compare == "different")
host_class_diff <- subset(host_class_diff,host_class_diff$host_phylum_compare != "different")

#By host order
host_order <- subset(type_dsDNA,type_dsDNA$host_order_compare != "different")
host_order <- subset(host_order,host_order$host_class_compare != "different")

host_order_diff <- subset(type_dsDNA,type_dsDNA$host_order_compare == "different")
host_order_diff <- subset(host_order_diff,host_order_diff$host_class_compare != "different")

#By host family
host_family <- subset(type_dsDNA,type_dsDNA$host_family_compare != "different")
host_family <- subset(host_family,host_family$host_order_compare != "different")

host_family_diff <- subset(type_dsDNA,type_dsDNA$host_family_compare == "different")
host_family_diff <- subset(host_family_diff,host_family_diff$host_order_compare != "different")

#By host genus
host_genus <- subset(type_dsDNA,type_dsDNA$host_genus_compare != "different")
host_genus <- subset(host_genus,host_genus$host_family_compare != "different")

host_genus_diff <- subset(type_dsDNA,type_dsDNA$host_genus_compare == "different")
host_genus_diff <- subset(host_genus_diff,host_genus_diff$host_family_compare != "different")






#Scatter plot of Mash vs Pham distributions
par(mar=c(4,8,4,4))
plot(host_superkingdom$modified_mash_distance,host_superkingdom$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_superkingdom_diff$modified_mash_distance,host_superkingdom_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_phylum$modified_mash_distance,host_phylum$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_phylum_diff$modified_mash_distance,host_phylum_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_class$modified_mash_distance,host_class$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_class_diff$modified_mash_distance,host_class_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_order$modified_mash_distance,host_order$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_order_diff$modified_mash_distance,host_order_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_family$modified_mash_distance,host_family$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_family_diff$modified_mash_distance,host_family_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_genus$modified_mash_distance,host_genus$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(host_genus_diff$modified_mash_distance,host_genus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





# par(mar=c(4,8,4,4))
# plot(host_genus_diff_actino$mash_distance,host_genus_diff_actino$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2)
# 
# par(mar=c(4,8,4,4))
# plot(host_genus_diff_bacter$mash_distance,host_genus_diff_bacter$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2)
# 
# par(mar=c(4,8,4,4))
# plot(host_genus_diff_cyano$mash_distance,host_genus_diff_cyano$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2)
# 
# par(mar=c(4,8,4,4))
# plot(host_genus_diff_firm$mash_distance,host_genus_diff_firm$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2)
# 
# par(mar=c(4,8,4,4))
# plot(host_genus_diff_proteo$mash_distance,host_genus_diff_proteo$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2)







#Histogram of Mash distance by host taxonomy
par(mar=c(4,8,4,4))
hist(host_superkingdom$modified_mash_distance,breaks=((range(host_superkingdom$modified_mash_distance)[2]-range(host_superkingdom$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_superkingdom_diff$modified_mash_distance,breaks=((range(host_superkingdom_diff$modified_mash_distance)[2]-range(host_superkingdom_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_phylum$modified_mash_distance,breaks=((range(host_phylum$modified_mash_distance)[2]-range(host_phylum$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_phylum_diff$modified_mash_distance,breaks=((range(host_phylum_diff$modified_mash_distance)[2]-range(host_phylum_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_class$modified_mash_distance,breaks=((range(host_class$modified_mash_distance)[2]-range(host_class$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_class_diff$modified_mash_distance,breaks=((range(host_class_diff$modified_mash_distance)[2]-range(host_class_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)







par(mar=c(4,8,4,4))
hist(host_order$modified_mash_distance,breaks=((range(host_order$modified_mash_distance)[2]-range(host_order$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_order_diff$modified_mash_distance,breaks=((range(host_order_diff$modified_mash_distance)[2]-range(host_order_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_family$modified_mash_distance,breaks=((range(host_family$modified_mash_distance)[2]-range(host_family$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_family_diff$modified_mash_distance,breaks=((range(host_family_diff$modified_mash_distance)[2]-range(host_family_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_genus$modified_mash_distance,breaks=((range(host_genus$modified_mash_distance)[2]-range(host_genus$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)

par(mar=c(4,8,4,4))
hist(host_genus_diff$modified_mash_distance,breaks=((range(host_genus_diff$modified_mash_distance)[2]-range(host_genus_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,6e3),col="black",cex.axis=2)















#Histogram distribution of pham distances by host phyla

# par(mar=c(4,8,4,4))
# hist(host_superkingdom$pham_pham_distance,breaks=((range(host_superkingdom$pham_pham_distance)[2]-range(host_superkingdom$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_superkingdom_diff$pham_pham_distance,xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_phylum$pham_pham_distance,breaks=((range(host_phylum$pham_pham_distance)[2]-range(host_phylum$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_phylum_diff$pham_pham_distance,breaks=((range(host_phylum_diff$pham_pham_distance)[2]-range(host_phylum_diff$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_class$pham_pham_distance,breaks=((range(host_class$pham_pham_distance)[2]-range(host_class$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_class_diff$pham_pham_distance,breaks=((range(host_class_diff$pham_pham_distance)[2]-range(host_class_diff$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_order$pham_pham_distance,breaks=((range(host_order$pham_pham_distance)[2]-range(host_order$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_order_diff$pham_pham_distance,breaks=((range(host_order_diff$pham_pham_distance)[2]-range(host_order_diff$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_family$pham_pham_distance,breaks=((range(host_family$pham_pham_distance)[2]-range(host_family$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_family_diff$pham_pham_distance,breaks=((range(host_family_diff$pham_pham_distance)[2]-range(host_family_diff$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_genus$pham_pham_distance,breaks=((range(host_genus$pham_pham_distance)[2]-range(host_genus$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)
# 
# par(mar=c(4,8,4,4))
# hist(host_genus_diff$pham_pham_distance,breaks=((range(host_genus_diff$pham_pham_distance)[2]-range(host_genus_diff$pham_pham_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,1),col="black",cex.axis=2)









###Check RNR data
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")

rnrdb_yes <- subset(type_dsDNA,type_dsDNA$phage_rnrdb_compare == "yes")
rnrdb_no <- subset(type_dsDNA,type_dsDNA$phage_rnrdb_compare == "no")
rnrdb_different <- subset(type_dsDNA,type_dsDNA$phage_rnrdb_compare == "different")


temperate <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "yes")
temperate_rnrdb_yes <- subset(temperate,temperate$phage_rnrdb_compare == "yes")
temperate_rnrdb_no <- subset(temperate,temperate$phage_rnrdb_compare == "no")

lytic <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "no")
lytic_rnrdb_yes <- subset(lytic,lytic$phage_rnrdb_compare == "yes")
lytic_rnrdb_no <- subset(lytic,lytic$phage_rnrdb_compare == "no")


par(mar=c(4,8,4,4))
plot(rnrdb_yes$modified_mash_distance,rnrdb_yes$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(rnrdb_no$modified_mash_distance,rnrdb_no$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(rnrdb_different$modified_mash_distance,rnrdb_different$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(temperate_rnrdb_yes$modified_mash_distance,temperate_rnrdb_yes$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_rnrdb_no$modified_mash_distance,temperate_rnrdb_no$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lytic_rnrdb_yes$modified_mash_distance,lytic_rnrdb_yes$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lytic_rnrdb_no$modified_mash_distance,lytic_rnrdb_no$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###By cluster and subcluster
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")
cluster_actino_same <- subset(cluster_actino,cluster_actino$phage_cluster_compare != "different")
cluster_actino_diff <- subset(cluster_actino,cluster_actino$phage_cluster_compare == "different")
subcluster_actino_same <- subset(cluster_actino_same,cluster_actino_same$phage_subcluster_compare != "different")
subcluster_actino_diff <- subset(cluster_actino_same,cluster_actino_same$phage_subcluster_compare == "different")

compute_sector_distribution(cluster_actino_same)
compute_sector_distribution(cluster_actino_diff)
compute_sector_distribution(subcluster_actino_same)
compute_sector_distribution(subcluster_actino_diff)


par(mar=c(4,8,4,4))
plot(cluster_actino_same$modified_mash_distance,cluster_actino_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_actino_diff$modified_mash_distance,cluster_actino_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(subcluster_actino_same$modified_mash_distance,subcluster_actino_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_actino_diff$modified_mash_distance,subcluster_actino_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")







par(mar=c(4,8,20,4))
hist(cluster_actino_same$modified_mash_distance,breaks=((range(cluster_actino_same$modified_mash_distance)[2]-range(cluster_actino_same$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_actino_diff$modified_mash_distance,breaks=((range(cluster_actino_diff$modified_mash_distance)[2]-range(cluster_actino_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e4),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(subcluster_actino_same$modified_mash_distance,breaks=((range(subcluster_actino_same$modified_mash_distance)[2]-range(subcluster_actino_same$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(subcluster_actino_diff$modified_mash_distance,breaks=((range(subcluster_actino_diff$modified_mash_distance)[2]-range(subcluster_actino_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,2e3),col="black",cex.axis=2)



par(mar=c(4,8,20,4))
hist(cluster_actino_same$pham_pham_dissimilarity,breaks=((range(cluster_actino_same$pham_pham_dissimilarity)[2]-range(cluster_actino_same$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)


par(mar=c(4,8,20,4))
hist(cluster_actino_diff$pham_pham_dissimilarity,breaks=((range(cluster_actino_diff$pham_pham_dissimilarity)[2]-range(cluster_actino_diff$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e4),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,20,4))
hist(subcluster_actino_same$pham_pham_dissimilarity,breaks=((range(subcluster_actino_same$pham_pham_dissimilarity)[2]-range(subcluster_actino_same$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,600),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,20,4))
hist(subcluster_actino_diff$pham_pham_dissimilarity,breaks=((range(subcluster_actino_diff$pham_pham_dissimilarity)[2]-range(subcluster_actino_diff$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)















###Cluster A Subcluster plots
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")
cluster_A <- subset(cluster_actino,cluster_actino$phage_cluster_compare == "A")
cluster_A$subcluster_A1_one_or_two <- ifelse(cluster_A$ref_phage_subcluster == "A1" | cluster_A$query_phage_subcluster == "A1",TRUE,FALSE)
cluster_A$subcluster_A1_neither <- ifelse(cluster_A$ref_phage_subcluster == "A1" | cluster_A$query_phage_subcluster == "A1",FALSE,TRUE)
cluster_A$subcluster_A1_both <- ifelse(cluster_A$ref_phage_subcluster == "A1" & cluster_A$query_phage_subcluster == "A1",TRUE,FALSE)
cluster_A$subcluster_A1_one <- ifelse(cluster_A$subcluster_A1_both == TRUE | cluster_A$subcluster_A1_neither == TRUE,FALSE,TRUE)

subcluster_A1 <- subset(cluster_A,cluster_A$phage_subcluster_compare == "A1")
subcluster_A1_one <- subset(cluster_A,cluster_A$subcluster_A1_one == TRUE)
subcluster_A1_both <- subset(cluster_A,cluster_A$subcluster_A1_both == TRUE)
subcluster_A1_neither <- subset(cluster_A,cluster_A$subcluster_A1_neither == TRUE)





par(mar=c(4,8,4,4))
plot(cluster_A$modified_mash_distance,cluster_A$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_A1_neither$modified_mash_distance,subcluster_A1_neither$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_A1_both$modified_mash_distance,subcluster_A1_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_A1_one$modified_mash_distance,subcluster_A1_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






#
#
#
###Check A1 and non-A1 distances to other non-A phages
#Since only actino phages that have been clustered are used, all rows should have cluster data.
#However, not all rows necessarily have subcluster data, so you have to take this into account.
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")



cluster_actino$subcluster_A1_one_or_two <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_neither <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,TRUE,FALSE)
cluster_actino$subcluster_A1_both <- ifelse((cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_one <- ifelse(cluster_actino$subcluster_A1_one_or_two == TRUE & cluster_actino$subcluster_A1_both == FALSE,TRUE,FALSE)

cluster_actino$cluster_A_one_or_two <- ifelse(cluster_actino$ref_phage_cluster == "A" | cluster_actino$query_phage_cluster == "A",TRUE,FALSE)
cluster_actino$cluster_A_both <- ifelse(cluster_actino$ref_phage_cluster == "A" & cluster_actino$query_phage_cluster == "A",TRUE,FALSE)
cluster_actino$cluster_A_both_but_notA1 <- ifelse(cluster_actino$cluster_A_both == TRUE & cluster_actino$subcluster_A1_neither == TRUE,TRUE,FALSE)
cluster_actino$cluster_A_one <- ifelse(cluster_actino$cluster_A_one_or_two == TRUE & cluster_actino$cluster_A_both == FALSE,TRUE,FALSE)
cluster_actino$cluster_A_one_but_notA1 <- ifelse(cluster_actino$cluster_A_one == TRUE & cluster_actino$subcluster_A1_neither == TRUE,TRUE,FALSE)










#Subset of all comparisons with at least one A1 phage
subcluster_A1_one_or_two_all <- subset(cluster_actino,cluster_actino$subcluster_A1_one_or_two == TRUE)

#Subset of all A1 comparisons
subcluster_A1_both <- subset(subcluster_A1_one_or_two_all,subcluster_A1_one_or_two_all$phage_subcluster_compare == "A1")

#Subset of all comparisons with one and only one A1 phage
subcluster_A1_one_all <- subset(subcluster_A1_one_or_two_all,subcluster_A1_one_or_two_all$subcluster_A1_one == TRUE)

#Subset of all comparisons with one and only one A1 phage, and the other phage is a Cluster A phage
subcluster_A1_one_clusterA <- subset(subcluster_A1_one_all,subcluster_A1_one_all$phage_cluster_compare == "A")

#Subset of all comparisons with one and only one A1 phage, and the other phage is NOT a Cluster A phage
#The phage_cluster_compare should NOT be "A", but it can be "NA" (since some phages are not clustered), or it can be "different"
subcluster_A1_one_other_not_clusterA <- subset(subcluster_A1_one_all,is.na(subcluster_A1_one_all$phage_cluster_compare) | subcluster_A1_one_all$phage_cluster_compare == "different")



#Subset of all comparisons with one and only non-A1, and no other cluster A
cluster_A_one_but_notA1 <- subset(cluster_actino,cluster_actino$cluster_A_one_but_notA1 == TRUE)

#Subset of all comparisons with both A, but no A1
cluster_A_both_but_notA1 <- subset(cluster_actino,cluster_actino$cluster_A_both_but_notA1 == TRUE)

nrow(subcluster_A1_both) +
  nrow(subcluster_A1_one_clusterA) +
  nrow(subcluster_A1_one_other_not_clusterA) +
  nrow(cluster_A_one_but_notA1) +
  nrow(cluster_A_both_but_notA1)


compute_sector_distribution(subcluster_A1_both)
compute_sector_distribution(subcluster_A1_one_clusterA)
compute_sector_distribution(subcluster_A1_one_other_not_clusterA)
compute_sector_distribution(cluster_A_both_but_notA1)
compute_sector_distribution(cluster_A_one_but_notA1)
compute_sector_distribution(type_dsDNA)

#Individual plots
par(mar=c(4,8,4,4))
plot(subcluster_A1_both$modified_mash_distance,subcluster_A1_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_A1_one_clusterA$modified_mash_distance,subcluster_A1_one_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_A1_one_other_not_clusterA$modified_mash_distance,subcluster_A1_one_other_not_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(cluster_A_both_but_notA1$modified_mash_distance,cluster_A_both_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(cluster_A_one_but_notA1$modified_mash_distance,cluster_A_one_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




#Colored on same plot
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








#Figure split up

par(mar=c(4,8,4,4))
plot(cluster_A_both_but_notA1$modified_mash_distance,cluster_A_both_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
par(new=TRUE)
plot(subcluster_A1_both$modified_mash_distance,subcluster_A1_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
par(new=TRUE)
plot(subcluster_A1_one_clusterA$modified_mash_distance,subcluster_A1_one_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")
par(new=TRUE)
plot(subcluster_A1_one_other_not_clusterA$modified_mash_distance,subcluster_A1_one_other_not_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(cluster_A_one_but_notA1$modified_mash_distance,cluster_A_one_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")






par(mar=c(4,8,4,4))
plot(cluster_A_both_but_notA1$modified_mash_distance,cluster_A_both_but_notA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(subcluster_A1_both$modified_mash_distance,subcluster_A1_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(subcluster_A1_one_clusterA$modified_mash_distance,subcluster_A1_one_clusterA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")
abline(0,2,lty=2,lwd=3,col="grey")













###Check clusters with Cluster A like architecture
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")

cluster_actino$cluster_BD_one_or_two <- ifelse(cluster_actino$ref_phage_cluster == "BD" | cluster_actino$query_phage_cluster == "BD",TRUE,FALSE)
cluster_actino$cluster_BD_both <- ifelse(cluster_actino$ref_phage_cluster == "BD" & cluster_actino$query_phage_cluster == "BD",TRUE,FALSE)
cluster_actino$cluster_BD_one <- ifelse(cluster_actino$cluster_BD_one_or_two == TRUE & cluster_actino$cluster_BD_both == FALSE,TRUE,FALSE)


cluster_actino$cluster_CA_one_or_two <- ifelse(cluster_actino$ref_phage_cluster == "CA" | cluster_actino$query_phage_cluster == "CA",TRUE,FALSE)
cluster_actino$cluster_CA_both <- ifelse(cluster_actino$ref_phage_cluster == "CA" & cluster_actino$query_phage_cluster == "CA",TRUE,FALSE)
cluster_actino$cluster_CA_one <- ifelse(cluster_actino$cluster_CA_one_or_two == TRUE & cluster_actino$cluster_CA_both == FALSE,TRUE,FALSE)

sum(cluster_actino$cluster_CA_one_or_two == TRUE) - (sum(cluster_actino$cluster_CA_both == TRUE)+ sum(cluster_actino$cluster_CA_one == TRUE))


cluster_actino$subcluster_A1_one_or_two <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_neither <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,TRUE,FALSE)
cluster_actino$subcluster_A1_both <- ifelse((cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_one <- ifelse(cluster_actino$subcluster_A1_one_or_two == TRUE & cluster_actino$subcluster_A1_both == FALSE,TRUE,FALSE)

cluster_actino$cluster_A_one_or_two <- ifelse(cluster_actino$ref_phage_cluster == "A" | cluster_actino$query_phage_cluster == "A",TRUE,FALSE)
cluster_actino$cluster_A_both <- ifelse(cluster_actino$ref_phage_cluster == "A" & cluster_actino$query_phage_cluster == "A",TRUE,FALSE)
cluster_actino$cluster_A_both_but_notA1 <- ifelse(cluster_actino$cluster_A_both == TRUE & cluster_actino$subcluster_A1_neither == TRUE,TRUE,FALSE)
cluster_actino$cluster_A_one <- ifelse(cluster_actino$cluster_A_one_or_two == TRUE & cluster_actino$cluster_A_both == FALSE,TRUE,FALSE)
cluster_actino$cluster_A_one_but_notA1 <- ifelse(cluster_actino$cluster_A_one == TRUE & cluster_actino$subcluster_A1_neither == TRUE,TRUE,FALSE)

cluster_actino$cluster_BD_one_A1_one <- ifelse(cluster_actino$cluster_BD_one == TRUE & cluster_actino$subcluster_A1_one == TRUE,TRUE,FALSE)
cluster_actino$cluster_BD_one_AnonA1_one <- ifelse(cluster_actino$cluster_BD_one == TRUE & cluster_actino$cluster_A_one_but_notA1 == TRUE,TRUE,FALSE)
cluster_actino$cluster_BD_one_nonA_one <- ifelse(cluster_actino$cluster_BD_one == TRUE & cluster_actino$cluster_A_one_or_two == FALSE,TRUE,FALSE)

sum(cluster_actino$cluster_BD_one_or_two == TRUE) - (sum(cluster_actino$cluster_BD_one_A1_one == TRUE) + sum(cluster_actino$cluster_BD_one_AnonA1_one == TRUE) + sum(cluster_actino$cluster_BD_one_nonA_one == TRUE) + sum(cluster_actino$cluster_BD_both == TRUE))


cluster_actino$cluster_CA_one_A1_one <- ifelse(cluster_actino$cluster_CA_one == TRUE & cluster_actino$subcluster_A1_one == TRUE,TRUE,FALSE)
cluster_actino$cluster_CA_one_AnonA1_one <- ifelse(cluster_actino$cluster_CA_one == TRUE & cluster_actino$cluster_A_one_but_notA1 == TRUE,TRUE,FALSE)
cluster_actino$cluster_CA_one_nonA_one <- ifelse(cluster_actino$cluster_CA_one == TRUE & cluster_actino$cluster_A_one_or_two == FALSE,TRUE,FALSE)

sum(cluster_actino$cluster_CA_one_or_two == TRUE) - (sum(cluster_actino$cluster_CA_one_A1_one == TRUE) + sum(cluster_actino$cluster_CA_one_AnonA1_one == TRUE) + sum(cluster_actino$cluster_CA_one_nonA_one == TRUE) + sum(cluster_actino$cluster_CA_both == TRUE))






cluster_actino$cluster_a_like_one_or_two <- ifelse(cluster_actino$ref_phage_cluster == "A" | cluster_actino$ref_phage_cluster == "BD" | cluster_actino$ref_phage_cluster == "CA" | cluster_actino$query_phage_cluster == "A" | cluster_actino$query_phage_cluster == "BD" | cluster_actino$query_phage_cluster == "CA",TRUE,FALSE)
cluster_actino$cluster_a_like_both <- ifelse((cluster_actino$ref_phage_cluster == "A" | cluster_actino$ref_phage_cluster == "BD" | cluster_actino$ref_phage_cluster == "CA") & (cluster_actino$query_phage_cluster == "A" | cluster_actino$query_phage_cluster == "BD" | cluster_actino$query_phage_cluster == "CA"),TRUE,FALSE)
cluster_actino$cluster_a_like_one <- ifelse(cluster_actino$cluster_a_like_one_or_two == TRUE & cluster_actino$cluster_a_like_both == FALSE,TRUE,FALSE)
cluster_actino$cluster_a_like_neither <- ifelse(cluster_actino$cluster_a_like_one_or_two == TRUE,FALSE,TRUE)


cluster_actino_diff_cluster_a_like_both <- subset(cluster_actino,cluster_actino$phage_cluster_compare == "different" & cluster_actino$cluster_a_like_both == TRUE)
cluster_actino_diff_cluster_a_like_neither <- subset(cluster_actino,cluster_actino$phage_cluster_compare == "different" & cluster_actino$cluster_a_like_neither == TRUE)

cluster_BD_both <- subset(cluster_actino,cluster_actino$cluster_BD_both == TRUE)
cluster_BD_one_A1_one <- subset(cluster_actino,cluster_actino$cluster_BD_one_A1_one == TRUE)
cluster_BD_one_AnonA1_one <- subset(cluster_actino,cluster_actino$cluster_BD_one_AnonA1_one == TRUE)
cluster_BD_one_nonA_one <- subset(cluster_actino,cluster_actino$cluster_BD_one_nonA_one == TRUE)



cluster_CA_both <- subset(cluster_actino,cluster_actino$cluster_CA_both == TRUE)
cluster_CA_one_A1_one <- subset(cluster_actino,cluster_actino$cluster_CA_one_A1_one == TRUE)
cluster_CA_one_AnonA1_one <- subset(cluster_actino,cluster_actino$cluster_CA_one_AnonA1_one == TRUE)
cluster_CA_one_nonA_one <- subset(cluster_actino,cluster_actino$cluster_CA_one_nonA_one == TRUE)



par(mar=c(4,8,4,4))
plot(cluster_CA_one_nonA_one$modified_mash_distance,cluster_CA_one_nonA_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(cluster_CA_one_A1_one$modified_mash_distance,cluster_CA_one_A1_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(cluster_CA_one_AnonA1_one$modified_mash_distance,cluster_CA_one_AnonA1_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(cluster_CA_both$modified_mash_distance,cluster_CA_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_BD_one_AnonA1_one$modified_mash_distance,cluster_BD_one_AnonA1_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(cluster_BD_one_nonA_one$modified_mash_distance,cluster_BD_one_nonA_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(cluster_BD_one_A1_one$modified_mash_distance,cluster_BD_one_A1_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(cluster_BD_both$modified_mash_distance,cluster_BD_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,2,lty=2,lwd=3,col="grey")











par(mar=c(4,8,4,4))
plot(cluster_actino_diff_cluster_a_like_both$modified_mash_distance,cluster_actino_diff_cluster_a_like_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_actino_diff_cluster_a_like_neither$modified_mash_distance,cluster_actino_diff_cluster_a_like_neither$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")










#
#
#
###Make intra-Cluster and intra-Subcluster boxplot ranked by median
#I WILL NEED TO USE MODIFIED MASH DATA

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")
cluster_actino_same <- subset(cluster_actino,cluster_actino$phage_cluster_compare != "different")
subcluster_actino_same <- subset(cluster_actino_same,cluster_actino_same$phage_subcluster_compare != "different")

#First re-factor the pertinent columns to remove unused levels that are present after subsetting data
cluster_actino_same$phage_cluster_compare <- factor(cluster_actino_same$phage_cluster_compare)
subcluster_actino_same$phage_subcluster_compare <- factor(subcluster_actino_same$phage_subcluster_compare)



#Now plot data

par(mar=c(4,2,4,15))
stripchart(cluster_actino_same$mash_distance ~ cluster_actino_same$phage_cluster_compare,at=rank(tapply(cluster_actino_same$mash_distance,cluster_actino_same$phage_cluster_compare,median)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)

par(mar=c(4,2,4,15))
stripchart(cluster_actino_same$pham_pham_dissimilarity ~ cluster_actino_same$phage_cluster_compare,at=rank(tapply(cluster_actino_same$mash_distance,cluster_actino_same$phage_cluster_compare,median)),las=1,col="black",xlim=c(0,1),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)


par(mar=c(4,2,4,15))
stripchart(subcluster_actino_same$mash_distance ~ subcluster_actino_same$phage_subcluster_compare,at=rank(tapply(subcluster_actino_same$mash_distance,subcluster_actino_same$phage_subcluster_compare,median)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)

par(mar=c(4,2,4,15))
stripchart(subcluster_actino_same$pham_pham_dissimilarity ~ subcluster_actino_same$phage_subcluster_compare,at=rank(tapply(subcluster_actino_same$mash_distance,subcluster_actino_same$phage_subcluster_compare,median)),las=1,col="black",xlim=c(0,1),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)


#
#
###Histogram of number of phages per cluster or subcluster
#STILL IN DEVELOPMENT


cluster_compare_frequency_table <- as.data.frame(summary(cluster_actino_same$phage_cluster_compare))
setDT(cluster_compare_frequency_table,keep.rownames = TRUE)
names(cluster_compare_frequency_table) = c("cluster","frequency")

cluster_compare_rank_by_median <- rank(tapply(cluster_actino_same$mash_distance,cluster_actino_same$phage_cluster_compare,median))
cluster_compare_rank_by_median_table <- as.data.frame(cluster_compare_rank_by_median)
setDT(cluster_compare_rank_by_median_table,keep.rownames = TRUE)
names(cluster_compare_rank_by_median_table) <- c("cluster","rank_by_median")

cluster_compare_frequency_table <- merge(cluster_compare_frequency_table,cluster_compare_rank_by_median_table,by.x="cluster",by.y="cluster")

par(mar=c(4,8,4,4))
barplot(cluster_compare_frequency_table[order(cluster_compare_frequency_table2$rank_by_median)]$frequency,horiz=TRUE,las=1,cex.axis=2)




subcluster_compare_frequency_table <- as.data.frame(summary(subcluster_actino_same$phage_subcluster_compare))
setDT(subcluster_compare_frequency_table,keep.rownames = TRUE)
names(subcluster_compare_frequency_table) = c("subcluster","frequency")

subcluster_compare_rank_by_median <- rank(tapply(subcluster_actino_same$mash_distance,subcluster_actino_same$phage_subcluster_compare,median))
subcluster_compare_rank_by_median_table <- as.data.frame(subcluster_compare_rank_by_median)
setDT(subcluster_compare_rank_by_median_table,keep.rownames = TRUE)
names(subcluster_compare_rank_by_median_table) <- c("subcluster","rank_by_median")

subcluster_compare_frequency_table <- merge(subcluster_compare_frequency_table,subcluster_compare_rank_by_median_table,by.x="subcluster",by.y="subcluster")

par(mar=c(4,8,4,4))
barplot(subcluster_compare_frequency_table[order(subcluster_compare_frequency_table$rank_by_median)]$frequency,horiz=TRUE,las=1,cex.axis=2)
#
#
#
#




















###Plots of Entero and Bacillus Clusters
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")

cluster_entero <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "entero")
cluster_entero_same <- subset(cluster_entero,cluster_entero$phage_cluster_compare != "different")
cluster_entero_diff <- subset(cluster_entero,cluster_entero$phage_cluster_compare == "different")

cluster_bacillus <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "bacillus")
cluster_bacillus_same <- subset(cluster_bacillus,cluster_bacillus$phage_cluster_compare != "different")
cluster_bacillus_diff <- subset(cluster_bacillus,cluster_bacillus$phage_cluster_compare == "different")



par(mar=c(4,8,4,4))
plot(cluster_entero_same$modified_mash_distance,cluster_entero_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_entero_diff$modified_mash_distance,cluster_entero_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(cluster_bacillus_same$modified_mash_distance,cluster_bacillus_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_bacillus_diff$modified_mash_distance,cluster_bacillus_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





par(mar=c(4,8,20,4))
hist(cluster_entero_same$modified_mash_distance,breaks=((range(cluster_entero_same$modified_mash_distance)[2]-range(cluster_entero_same$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,200),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_entero_diff$modified_mash_distance,breaks=((range(cluster_entero_diff$modified_mash_distance)[2]-range(cluster_entero_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,200),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_bacillus_same$modified_mash_distance,breaks=((range(cluster_bacillus_same$modified_mash_distance)[2]-range(cluster_bacillus_same$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,40),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_bacillus_diff$modified_mash_distance,breaks=((range(cluster_bacillus_diff$modified_mash_distance)[2]-range(cluster_bacillus_diff$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,40),col="black",cex.axis=2)







par(mar=c(4,8,20,4))
hist(cluster_entero_same$pham_pham_dissimilarity,breaks=((range(cluster_entero_same$pham_pham_dissimilarity)[2]-range(cluster_entero_same$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,200),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_entero_diff$pham_pham_dissimilarity,breaks=((range(cluster_entero_diff$pham_pham_dissimilarity)[2]-range(cluster_entero_diff$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,1e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_bacillus_same$pham_pham_dissimilarity,breaks=((range(cluster_bacillus_same$pham_pham_dissimilarity)[2]-range(cluster_bacillus_same$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,40),yaxt="n")
axis(side=4,pos=0,cex.axis=2)

par(mar=c(4,8,20,4))
hist(cluster_bacillus_diff$pham_pham_dissimilarity,breaks=((range(cluster_bacillus_diff$pham_pham_dissimilarity)[2]-range(cluster_bacillus_diff$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,40),yaxt="n")
axis(side=4,pos=0,cex.axis=2)


















###Specific clusters
#List of Clusters to investigate
# A
# B
# F
# ZB = T4
# BU
# C
# E
# ZE = T7
# K
# BD
# YV
# YX
# G
# L
# ZC
# XG
# ZZ
# 
# AK
# YE
# ZA
# J
# N
# ZN
# AN
# YG
# YJ
# ZI
# ZR
# D
# CA
# P
# ZF
# 
# CV
# YS
# ZD
# ZM
# YB
# AO
# CZ
# O
# YD
# YQ
# ZK
# ZO
# ZL
# AQ
# CQ
# CS
# H
# I
# Q
# XE
# 
# BG
# BW
# DA
# M
# XA
# XC
# XD
# YH
# YM
# YU
# Singleton60
# ZG
# BB
# BV
# CR
# CT
# CX
# DE
# DF
# R
# XN
# YA
# YC
# 
# YN
# YW
# YZ
# ZH
# ZP
# ZS
# ZT
# ZV
# ZW = Lambda
# ZX
# YK
# ZY

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster <- subset(type_dsDNA,type_dsDNA$phage_cluster_compare != "different")
#cluster_diff <- subset(mash_host_pham_pvalue_size,mash_host_pham_pvalue_size$phage_cluster_compare == "different")

cluster_A <- subset(cluster,cluster$phage_cluster_compare == "A")
cluster_F <- subset(cluster,cluster$phage_cluster_compare == "F")
cluster_BD <- subset(cluster,cluster$phage_cluster_compare == "BD")
cluster_AO <- subset(cluster,cluster$phage_cluster_compare == "AO")
cluster_N <- subset(cluster,cluster$phage_cluster_compare == "N")



cluster_ZW <- subset(cluster,cluster$phage_cluster_compare == "ZW")
cluster_ZE <- subset(cluster,cluster$phage_cluster_compare == "ZE")
cluster_ZB <- subset(cluster,cluster$phage_cluster_compare == "ZB")


cluster_B <- subset(cluster,cluster$phage_cluster_compare == "B")
cluster_BU <- subset(cluster,cluster$phage_cluster_compare == "BU")
cluster_K <- subset(cluster,cluster$phage_cluster_compare == "K")
cluster_CZ <- subset(cluster,cluster$phage_cluster_compare == "CZ")
cluster_YG <- subset(cluster,cluster$phage_cluster_compare == "YG")
cluster_ZZ <- subset(cluster,cluster$phage_cluster_compare == "ZZ")

cluster_J <- subset(cluster,cluster$phage_cluster_compare == "J")
cluster_P <- subset(cluster,cluster$phage_cluster_compare == "P")
cluster_I <- subset(cluster,cluster$phage_cluster_compare == "I")
cluster_E <- subset(cluster,cluster$phage_cluster_compare == "E")
cluster_O <- subset(cluster,cluster$phage_cluster_compare == "O")
cluster_G <- subset(cluster,cluster$phage_cluster_compare == "G")



# cluster_A <- subset(cluster,cluster$phage_cluster_compare == "A")
# cluster_B <- subset(cluster,cluster$phage_cluster_compare == "B")
# cluster_F <- subset(cluster,cluster$phage_cluster_compare == "F")
# cluster_ZB <- subset(cluster,cluster$phage_cluster_compare == "ZB")
# cluster_BU <- subset(cluster,cluster$phage_cluster_compare == "BU")
# cluster_C <- subset(cluster,cluster$phage_cluster_compare == "C")
# cluster_E <- subset(cluster,cluster$phage_cluster_compare == "E")
# cluster_ZE <- subset(cluster,cluster$phage_cluster_compare == "ZE")
# cluster_K <- subset(cluster,cluster$phage_cluster_compare == "K")
# cluster_BD <- subset(cluster,cluster$phage_cluster_compare == "BD")
# cluster_YV <- subset(cluster,cluster$phage_cluster_compare == "YV")
# cluster_YX <- subset(cluster,cluster$phage_cluster_compare == "YX")
# cluster_G <- subset(cluster,cluster$phage_cluster_compare == "G")
# cluster_L <- subset(cluster,cluster$phage_cluster_compare == "L")
# cluster_ZC <- subset(cluster,cluster$phage_cluster_compare == "ZC")
# cluster_XG <- subset(cluster,cluster$phage_cluster_compare == "XG")
# cluster_ZZ <- subset(cluster,cluster$phage_cluster_compare == "ZZ")
# cluster_AK <- subset(cluster,cluster$phage_cluster_compare == "AK")
# cluster_YE <- subset(cluster,cluster$phage_cluster_compare == "YE")
# cluster_ZA <- subset(cluster,cluster$phage_cluster_compare == "ZA")
# cluster_J <- subset(cluster,cluster$phage_cluster_compare == "J")
# cluster_N <- subset(cluster,cluster$phage_cluster_compare == "N")
# cluster_ZN <- subset(cluster,cluster$phage_cluster_compare == "ZN")
# cluster_AN <- subset(cluster,cluster$phage_cluster_compare == "AN")
# cluster_YG <- subset(cluster,cluster$phage_cluster_compare == "YG")
# cluster_YJ <- subset(cluster,cluster$phage_cluster_compare == "YJ")
# cluster_ZI <- subset(cluster,cluster$phage_cluster_compare == "ZI")
# cluster_ZR <- subset(cluster,cluster$phage_cluster_compare == "ZR")
# cluster_D <- subset(cluster,cluster$phage_cluster_compare == "D")
# cluster_CA <- subset(cluster,cluster$phage_cluster_compare == "CA")
# cluster_P <- subset(cluster,cluster$phage_cluster_compare == "P")
# cluster_ZF <- subset(cluster,cluster$phage_cluster_compare == "ZF")
# cluster_CV <- subset(cluster,cluster$phage_cluster_compare == "CV")
# cluster_YS <- subset(cluster,cluster$phage_cluster_compare == "YS")
# cluster_ZD <- subset(cluster,cluster$phage_cluster_compare == "ZD")
# cluster_ZM <- subset(cluster,cluster$phage_cluster_compare == "ZM")
# cluster_YB <- subset(cluster,cluster$phage_cluster_compare == "YB")
# cluster_AO <- subset(cluster,cluster$phage_cluster_compare == "AO")
# cluster_CZ <- subset(cluster,cluster$phage_cluster_compare == "CZ")
# cluster_O <- subset(cluster,cluster$phage_cluster_compare == "O")
# cluster_YD <- subset(cluster,cluster$phage_cluster_compare == "YD")
# cluster_YQ <- subset(cluster,cluster$phage_cluster_compare == "YQ")
# cluster_ZK <- subset(cluster,cluster$phage_cluster_compare == "ZK")
# cluster_ZO <- subset(cluster,cluster$phage_cluster_compare == "ZO")
# cluster_ZL <- subset(cluster,cluster$phage_cluster_compare == "ZL")
# cluster_AQ <- subset(cluster,cluster$phage_cluster_compare == "AQ")
# cluster_CQ <- subset(cluster,cluster$phage_cluster_compare == "CQ")
# cluster_CS <- subset(cluster,cluster$phage_cluster_compare == "CS")
# cluster_H <- subset(cluster,cluster$phage_cluster_compare == "H")
# cluster_I <- subset(cluster,cluster$phage_cluster_compare == "I")
# cluster_Q <- subset(cluster,cluster$phage_cluster_compare == "Q")
# cluster_XE <- subset(cluster,cluster$phage_cluster_compare == "XE")
# cluster_BG <- subset(cluster,cluster$phage_cluster_compare == "BG")
# cluster_BW <- subset(cluster,cluster$phage_cluster_compare == "BW")
# cluster_DA <- subset(cluster,cluster$phage_cluster_compare == "DA")
# cluster_M <- subset(cluster,cluster$phage_cluster_compare == "M")
# cluster_XA <- subset(cluster,cluster$phage_cluster_compare == "XA")
# cluster_XC <- subset(cluster,cluster$phage_cluster_compare == "XC")
# cluster_XD <- subset(cluster,cluster$phage_cluster_compare == "XD")
# cluster_YH <- subset(cluster,cluster$phage_cluster_compare == "YH")
# cluster_YM <- subset(cluster,cluster$phage_cluster_compare == "YM")
# cluster_YU <- subset(cluster,cluster$phage_cluster_compare == "YU")
# cluster_Singleton60 <- subset(cluster,cluster$phage_cluster_compare == "Singleton60")
# cluster_ZG <- subset(cluster,cluster$phage_cluster_compare == "ZG")
# cluster_BB <- subset(cluster,cluster$phage_cluster_compare == "BB")
# cluster_BV <- subset(cluster,cluster$phage_cluster_compare == "BV")
# cluster_CR <- subset(cluster,cluster$phage_cluster_compare == "CR")
# cluster_CT <- subset(cluster,cluster$phage_cluster_compare == "CT")
# cluster_CX <- subset(cluster,cluster$phage_cluster_compare == "CX")
# cluster_DE <- subset(cluster,cluster$phage_cluster_compare == "DE")
# cluster_DF <- subset(cluster,cluster$phage_cluster_compare == "DF")
# cluster_R <- subset(cluster,cluster$phage_cluster_compare == "R")
# cluster_XN <- subset(cluster,cluster$phage_cluster_compare == "XN")
# cluster_YA <- subset(cluster,cluster$phage_cluster_compare == "YA")
# cluster_YC <- subset(cluster,cluster$phage_cluster_compare == "YC")
# cluster_YN <- subset(cluster,cluster$phage_cluster_compare == "YN")
# cluster_YW <- subset(cluster,cluster$phage_cluster_compare == "YW")
# cluster_YZ <- subset(cluster,cluster$phage_cluster_compare == "YZ")
# cluster_ZH <- subset(cluster,cluster$phage_cluster_compare == "ZH")
# cluster_ZP <- subset(cluster,cluster$phage_cluster_compare == "ZP")
# cluster_ZS <- subset(cluster,cluster$phage_cluster_compare == "ZS")
# cluster_ZT <- subset(cluster,cluster$phage_cluster_compare == "ZT")
# cluster_ZV <- subset(cluster,cluster$phage_cluster_compare == "ZV")
# cluster_ZW <- subset(cluster,cluster$phage_cluster_compare == "ZW")
# cluster_ZX <- subset(cluster,cluster$phage_cluster_compare == "ZX")
# cluster_YK <- subset(cluster,cluster$phage_cluster_compare == "YK")
# cluster_ZY <- subset(cluster,cluster$phage_cluster_compare == "ZY")




#Lines to demarcate groups
abline(1.35,-2,lty=2,lwd=3,col="grey")
abline(0.4,-2,lty=2,lwd=3,col="grey")
abline(0,3.5,lty=2,lwd=3,col="grey")
abline(0.25,2,lty=2,lwd=3,col="grey")



par(mar=c(4,8,4,4))
plot(cluster_F$modified_mash_distance,cluster_F$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(cluster_BD$modified_mash_distance,cluster_BD$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_AO$modified_mash_distance,cluster_AO$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_N$modified_mash_distance,cluster_N$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





par(mar=c(4,8,4,4))
plot(cluster_ZW$modified_mash_distance,cluster_ZW$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_ZE$modified_mash_distance,cluster_ZE$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_ZB$modified_mash_distance,cluster_ZB$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






par(mar=c(4,8,4,4))
plot(cluster_B$modified_mash_distance,cluster_B$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_BU$modified_mash_distance,cluster_BU$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_K$modified_mash_distance,cluster_K$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_CZ$modified_mash_distance,cluster_CZ$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_YG$modified_mash_distance,cluster_YG$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_ZZ$modified_mash_distance,cluster_ZZ$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







par(mar=c(4,8,4,4))
plot(cluster_J$modified_mash_distance,cluster_J$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_P$modified_mash_distance,cluster_P$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_I$modified_mash_distance,cluster_I$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_E$modified_mash_distance,cluster_E$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_O$modified_mash_distance,cluster_O$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_G$modified_mash_distance,cluster_G$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_A$modified_mash_distance,cluster_A$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")












###Compare Mash to ANI
#These commands are to verify that the evolutionary patterns we see with pham-mash comparisons are also present in pham-ani comparisons


ani_mash <- subset(mash_table2,complete.cases(mash_table2$ani_distance))

ani_mash_clusterg <- subset(ani_mash,ani_mash$phage_cluster_compare=="G")
ani_mash_clusterj <- subset(ani_mash,ani_mash$phage_cluster_compare=="J")
ani_mash_clusterl <- subset(ani_mash,ani_mash$phage_cluster_compare=="L")
ani_mash_clustern <- subset(ani_mash,ani_mash$phage_cluster_compare=="N")
ani_mash_clusteryv <- subset(ani_mash,ani_mash$phage_cluster_compare=="YV")
ani_mash_clusterf <- subset(ani_mash,ani_mash$phage_cluster_compare=="F")

#Compare ani to mash
par(mar=c(4,8,4,4))
plot(ani_mash$ani_distance,ani_mash$modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")


#Now check specific clusters using mash
par(mar=c(4,8,4,4))
plot(ani_mash_clusterg$modified_mash_distance,ani_mash_clusterg$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusterj$modified_mash_distance,ani_mash_clusterj$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusterl$modified_mash_distance,ani_mash_clusterl$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clustern$modified_mash_distance,ani_mash_clustern$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusteryv$modified_mash_distance,ani_mash_clusteryv$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusterf$modified_mash_distance,ani_mash_clusterf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


#Now check specific clusters using ani
par(mar=c(4,8,4,4))
plot(ani_mash_clusterg$ani_distance,ani_mash_clusterg$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusterj$ani_distance,ani_mash_clusterj$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusterl$ani_distance,ani_mash_clusterl$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clustern$ani_distance,ani_mash_clustern$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusteryv$ani_distance,ani_mash_clusteryv$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani_mash_clusterf$ani_distance,ani_mash_clusterf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




#Colored plot by assigne evolutionary mode 
ani_mash_hgcf <- subset(ani_mash,ani_mash$gene_flux_category == "high" & ani_mash$phage_predicted_temperate_compare == "yes")
ani_mash_lgcf <- subset(ani_mash,ani_mash$gene_flux_category == "low" & ani_mash$phage_predicted_temperate_compare == "yes")
ani_mash_lytic <- subset(ani_mash,ani_mash$phage_predicted_temperate_compare == "no")



par(mar=c(4,8,4,4))
plot(ani_mash_hgcf$modified_mash_distance,ani_mash_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(ani_mash_lgcf$modified_mash_distance,ani_mash_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(ani_mash_lytic$modified_mash_distance,ani_mash_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


par(mar=c(4,8,4,4))
plot(ani_mash_hgcf$ani_distance,ani_mash_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(ani_mash_lgcf$ani_distance,ani_mash_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(ani_mash_lytic$ani_distance,ani_mash_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


###CURRENT
#temp
ani_mash_hgcf_subset <- subset(ani_mash_hgcf,
                               select=c('mash_reference','mash_query','mash_count','modified_mash_distance','pham_pham_dissimilarity','ref_size','query_size','ani_distance'))

ani_mash_lgcf_subset <- subset(ani_mash_lgcf,
                     select=c('mash_reference','mash_query','mash_count','modified_mash_distance','pham_pham_dissimilarity','ref_size','query_size','ani_distance'))
#





















###pham dissimilarity vs jaccard dissimilarity


mash_table2_filtered <- subset(mash_table2,mash_table2$filter == TRUE)

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
bacteria_dsDNA <- subset(type_dsDNA,type_dsDNA$host_superkingdom_compare == 'Bacteria')
host_phylum <- subset(type_dsDNA,type_dsDNA$host_phylum_compare != "different")
actino <- subset(host_phylum,host_phylum$host_phylum_compare == "Actinobacteria")



#Check how correlated pham dissimilarity and jaccard dissimilarity are
par(mar=c(4,8,4,4))
plot(mash_table2$pham_jaccard_dissimilarity,mash_table2$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_table2_filtered$pham_jaccard_dissimilarity,mash_table2_filtered$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")



#Mash vs Pham plot using jaccard
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


cluster_g <- subset(actino,actino$phage_cluster_compare == "G")
cluster_j <- subset(actino,actino$phage_cluster_compare == "J")
cluster_l <- subset(actino,actino$phage_cluster_compare == "L")
cluster_n <- subset(actino,actino$phage_cluster_compare == "N")

par(mar=c(4,8,4,4))
plot(actino$modified_mash_distance,actino$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




par(mar=c(4,8,4,4))
plot(cluster_g$ani_distance,cluster_g$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_j$ani_distance,cluster_j$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_l$ani_distance,cluster_l$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_n$ani_distance,cluster_n$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



















#
#
#
###Evolutionary analysis of Mash Groups
#Subset by mash groups


mash_groups_d030 <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d030_compare!="different")
mash_groups_d020 <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d020_compare!="different")
mash_groups_d010 <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d010_compare!="different")


mash_groups_d030_diff <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d030_compare=="different")
mash_groups_d020_diff <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d020_compare=="different")
mash_groups_d010_diff <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d010_compare=="different")




plot(mash_groups_d010$mash_distance,mash_groups_d010$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 same",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.000001)
abline(0,2.3)
plot(mash_groups_d010_diff$mash_distance,mash_groups_d010_diff$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 diff",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.000001)
abline(0,2.3)





plot(mash_groups_d020$mash_distance,mash_groups_d020$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 same",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.000001)
abline(0,2.3)
plot(mash_groups_d020_diff$mash_distance,mash_groups_d020_diff$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 diff",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.000001)
abline(0,2.3)




plot(mash_groups_d030$mash_distance,mash_groups_d030$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d030 same",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.000001)
abline(0,2.3)
plot(mash_groups_d030_diff$mash_distance,mash_groups_d030_diff$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d030 diff",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.000001)
abline(0,2.3)




mash_groups_d010_actino <- subset(mash_groups_d010,mash_groups_d010$host_phylum_compare=="Actinobacteria")
mash_groups_d010_bacter <- subset(mash_groups_d010,mash_groups_d010$host_phylum_compare=="Bacteroidetes")
mash_groups_d010_cyano <- subset(mash_groups_d010,mash_groups_d010$host_phylum_compare=="Cyanobacteria")
mash_groups_d010_firm <- subset(mash_groups_d010,mash_groups_d010$host_phylum_compare=="Firmicutes")
mash_groups_d010_proteo <- subset(mash_groups_d010,mash_groups_d010$host_phylum_compare=="Proteobacteria")


mash_groups_d020_actino <- subset(mash_groups_d020,mash_groups_d020$host_phylum_compare=="Actinobacteria")
mash_groups_d020_bacter <- subset(mash_groups_d020,mash_groups_d020$host_phylum_compare=="Bacteroidetes")
mash_groups_d020_cyano <- subset(mash_groups_d020,mash_groups_d020$host_phylum_compare=="Cyanobacteria")
mash_groups_d020_firm <- subset(mash_groups_d020,mash_groups_d020$host_phylum_compare=="Firmicutes")
mash_groups_d020_proteo <- subset(mash_groups_d020,mash_groups_d020$host_phylum_compare=="Proteobacteria")


plot(actino$mash_distance,actino$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size actino",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d020_actino$mash_distance,mash_groups_d020_actino$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 actino",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_actino$mash_distance,mash_groups_d010_actino$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 actino",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)


plot(bacter$mash_distance,bacter$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size bacter",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d020_bacter$mash_distance,mash_groups_d020_bacter$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 bacter",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_bacter$mash_distance,mash_groups_d010_bacter$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 bacter",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)


plot(cyano$mash_distance,cyano$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size cyano",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d020_cyano$mash_distance,mash_groups_d020_cyano$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 cyano",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_cyano$mash_distance,mash_groups_d010_cyano$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 cyano",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)


plot(firm$mash_distance,firm$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size firm",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d020_firm$mash_distance,mash_groups_d020_firm$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 firm",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_firm$mash_distance,mash_groups_d010_firm$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 firm",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)


plot(proteo$mash_distance,proteo$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size proteo",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d020_proteo$mash_distance,mash_groups_d020_proteo$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d020 proteo",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_proteo$mash_distance,mash_groups_d010_proteo$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
#
#
#



#Specific Mash Groups

#
#
#
mash_groups_d010_firm_mg14 <- subset(mash_groups_d010_firm,mash_groups_d010_firm$mash_group_d010_compare=="Mash_Group_14")
mash_groups_d010_firm_mg129 <- subset(mash_groups_d010_firm,mash_groups_d010_firm$mash_group_d010_compare=="Mash_Group_129")
mash_groups_d010_firm_mg224 <- subset(mash_groups_d010_firm,mash_groups_d010_firm$mash_group_d010_compare=="Mash_Group_224")
mash_groups_d010_firm_mg367 <- subset(mash_groups_d010_firm,mash_groups_d010_firm$mash_group_d010_compare=="Mash_Group_367")

plot(mash_groups_d010_firm_mg14$mash_distance,mash_groups_d010_firm_mg14$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 firm MG14",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)

plot(mash_groups_d010_firm_mg129$mash_distance,mash_groups_d010_firm_mg129$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 firm MG129",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)

plot(mash_groups_d010_firm_mg224$mash_distance,mash_groups_d010_firm_mg224$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 firm MG224",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)

plot(mash_groups_d010_firm_mg367$mash_distance,mash_groups_d010_firm_mg367$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 firm MG367",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
#
#
#




#Lambda = MG_d010_210
#Only 3 comparisons
#mash_groups_d010_210 <- subset(mash_table2,mash_table2$mash_group_d010_compare=="Mash_Group_210")


#Lambda = MG_d030_391
mash_groups_d030_391 <- subset(mash_table2,mash_table2$mash_group_d030_compare=="Mash_Group_391")

par(mar=c(4,8,4,4))
plot(mash_groups_d030_391$modified_mash_distance,mash_groups_d030_391$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#T7 = MG_d030_392
mash_groups_d030_392 <- subset(mash_table2,mash_table2$mash_group_d030_compare=="Mash_Group_392")

par(mar=c(4,8,4,4))
plot(mash_groups_d030_392$modified_mash_distance,mash_groups_d030_392$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#T4 = MG_d030_191
mash_groups_d030_191 <- subset(mash_table2,mash_table2$mash_group_d030_compare=="Mash_Group_191")

par(mar=c(4,8,4,4))
plot(mash_groups_d030_191$modified_mash_distance,mash_groups_d030_191$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Mycobacterium phage large mash group
mash_groups_d030_196 <- subset(mash_table2,mash_table2$mash_group_d030_compare=="Mash_Group_196")

par(mar=c(4,8,4,4))
plot(mash_groups_d030_196$modified_mash_distance,mash_groups_d030_196$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




#
#
#
mash_groups_d010_proteo_mg19 <- subset(mash_groups_d010_proteo,mash_groups_d010_proteo$mash_group_d010_compare=="Mash_Group_19")
mash_groups_d010_proteo_mg227 <- subset(mash_groups_d010_proteo,mash_groups_d010_proteo$mash_group_d010_compare=="Mash_Group_227")
mash_groups_d010_proteo_mg155 <- subset(mash_groups_d010_proteo,mash_groups_d010_proteo$mash_group_d010_compare=="Mash_Group_155")
mash_groups_d010_proteo_mg252 <- subset(mash_groups_d010_proteo,mash_groups_d010_proteo$mash_group_d010_compare=="Mash_Group_252")
mash_groups_d010_proteo_mg93 <- subset(mash_groups_d010_proteo,mash_groups_d010_proteo$mash_group_d010_compare=="Mash_Group_93")
mash_groups_d010_proteo_mg3 <- subset(mash_groups_d010_proteo,mash_groups_d010_proteo$mash_group_d010_compare=="Mash_Group_3")



plot(mash_groups_d010_proteo_mg19$mash_distance,mash_groups_d010_proteo_mg19$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo MG19",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_proteo_mg227$mash_distance,mash_groups_d010_proteo_mg227$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo MG227",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_proteo_mg155$mash_distance,mash_groups_d010_proteo_mg155$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo MG155",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_proteo_mg252$mash_distance,mash_groups_d010_proteo_mg252$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo MG252",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_proteo_mg93$mash_distance,mash_groups_d010_proteo_mg93$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo MG93",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
plot(mash_groups_d010_proteo_mg3$mash_distance,mash_groups_d010_proteo_mg3$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size mash_groups_d010 proteo MG3",xlim=c(0,0.6),ylim=c(0,1),pch=20,cex=0.1)
abline(0,2.3)
#
#
#









###Histograms of Mash Group Membership
library(plyr)



mash_group_d005_count <- count(mash_group_table,'mash_group_d005')
mash_group_d010_count <- count(mash_group_table,'mash_group_d010')
mash_group_d020_count <- count(mash_group_table,'mash_group_d020')
mash_group_d030_count <- count(mash_group_table,'mash_group_d030')
mash_group_d040_count <- count(mash_group_table,'mash_group_d040')

mash_group_d005_count$freq_levels <- as.factor(mash_group_d005_count$freq)
mash_group_d010_count$freq_levels <- as.factor(mash_group_d010_count$freq)
mash_group_d020_count$freq_levels <- as.factor(mash_group_d020_count$freq)
mash_group_d030_count$freq_levels <- as.factor(mash_group_d030_count$freq)
mash_group_d040_count$freq_levels <- as.factor(mash_group_d040_count$freq)

par(mar=c(4,8,4,4))
barplot(summary(mash_group_d005_count$freq_levels),las=2,xlab="",ylab="",cex.axis=2)

par(mar=c(4,8,4,4))
barplot(summary(mash_group_d010_count$freq_levels),las=2,xlab="",ylab="",cex.axis=2)

par(mar=c(4,8,4,4))
barplot(summary(mash_group_d020_count$freq_levels),las=2,xlab="",ylab="",cex.axis=2)

par(mar=c(4,8,4,4))
barplot(summary(mash_group_d030_count$freq_levels),las=2,xlab="",ylab="",cex.axis=2)

par(mar=c(4,8,4,4))
barplot(summary(mash_group_d040_count$freq_levels),las=2,xlab="",ylab="",cex.axis=2)
axis(1,las=2,cex.axis=0.5)






#Max Mash Group Distance
d005_distances <- subset(mash_table2,mash_table2$mash_group_d005_compare != "different")
d010_distances <- subset(mash_table2,mash_table2$mash_group_d010_compare != "different")
d015_distances <- subset(mash_table2,mash_table2$mash_group_d015_compare != "different")
d020_distances <- subset(mash_table2,mash_table2$mash_group_d020_compare != "different")
d030_distances <- subset(mash_table2,mash_table2$mash_group_d030_compare != "different")


par(mar=c(4,2,4,15))
stripchart(d005_distances$modified_mash_distance ~ d005_distances$mash_group_d005_compare,at=rank(tapply(d005_distances$modified_mash_distance,d005_distances$mash_group_d005_compare,max)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)
abline(v=0.05,lty=2,lwd=3,col="grey")

par(mar=c(4,2,4,15))
stripchart(d010_distances$modified_mash_distance ~ d010_distances$mash_group_d010_compare,at=rank(tapply(d010_distances$modified_mash_distance,d010_distances$mash_group_d010_compare,max)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)
abline(v=0.1,lty=2,lwd=3,col="grey")

par(mar=c(4,2,4,15))
stripchart(d015_distances$modified_mash_distance ~ d015_distances$mash_group_d015_compare,at=rank(tapply(d015_distances$modified_mash_distance,d015_distances$mash_group_d015_compare,max)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)
abline(v=0.15,lty=2,lwd=3,col="grey")



par(mar=c(4,2,4,15))
stripchart(d020_distances$modified_mash_distance ~ d020_distances$mash_group_d020_compare,at=rank(tapply(d020_distances$modified_mash_distance,d020_distances$mash_group_d020_compare,max)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)
abline(v=0.2,lty=2,lwd=3,col="grey")

par(mar=c(4,2,4,15))
stripchart(d030_distances$modified_mash_distance ~ d030_distances$mash_group_d030_compare,at=rank(tapply(d030_distances$modified_mash_distance,d030_distances$mash_group_d030_compare,max)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)
abline(v=0.3,lty=2,lwd=3,col="grey")






d010_mash_group_max_distances <- tapply(d010_distances$modified_mash_distance,d010_distances$mash_group_d010_compare,max)
d015_mash_group_max_distances <- tapply(d015_distances$modified_mash_distance,d015_distances$mash_group_d015_compare,max)
sort(d015_mash_group_max_distances)

mash_groups_d010_323 <- subset(mash_table2,mash_table2$mash_group_d010_compare=="Mash_Group_323")
mash_groups_d010_311 <- subset(mash_table2,mash_table2$mash_group_d010_compare=="Mash_Group_311")
mash_groups_d010_3 <- subset(mash_table2,mash_table2$mash_group_d010_compare=="Mash_Group_3")



mash_groups_d015_443 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_443")
mash_groups_d015_373 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_373")
mash_groups_d015_331 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_331")
mash_groups_d015_30 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_30")
mash_groups_d015_3 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_3")





par(mar=c(4,8,4,4))
plot(mash_groups_d010_323$modified_mash_distance,mash_groups_d010_323$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_groups_d010_311$modified_mash_distance,mash_groups_d010_311$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_groups_d010_3$modified_mash_distance,mash_groups_d010_3$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






mash_groups_d015_443 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_443")
mash_groups_d015_373 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_373")
mash_groups_d015_331 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_331")
mash_groups_d015_30 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_30")
mash_groups_d015_3 <- subset(mash_table2,mash_table2$mash_group_d015_compare=="Mash_Group_3")

par(mar=c(4,8,4,4))
plot(mash_groups_d015_443$modified_mash_distance,mash_groups_d015_443$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_groups_d015_373$modified_mash_distance,mash_groups_d015_373$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_groups_d015_331$modified_mash_distance,mash_groups_d015_331$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_groups_d015_30$modified_mash_distance,mash_groups_d015_30$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mash_groups_d015_3$modified_mash_distance,mash_groups_d015_3$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




###Output updated mash data for cytoscape, with no significant data retained. 
#This is an alternative to the filtering script, because here it takes into account other factors such as genome size disparity.
#Also, this data file will contain gene content dissimilarity data as well.
#All modified_mash_distance data takes into account the filtering parameters and pvalue

mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
bacteria_dsDNA <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')

bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.42 & bacteria_dsDNA$pham_pham_dissimilarity < 0.89)


#Plot the data to make sure it contains only what I want
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA_nuc042_gene089$modified_mash_distance,bacteria_dsDNA_nuc042_gene089$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Now reduce the data table to only the columns needed for export, and export the data
bacteria_dsDNA_nuc042_gene089_reduced <- subset(bacteria_dsDNA_nuc042_gene089,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
write.table(bacteria_dsDNA_nuc042_gene089_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_bacteria_dsDNA_nuc042_gene089_data.csv",sep=",",row.names = FALSE,quote=FALSE)




#Below: old script used for cytoscape maps, before I took into account filtering out comparisons using gene content dissimilarity as well...
mash_filtered_d030 <- subset(mash_filtered,mash_filtered$mash_distance < 0.3)
mash_filtered_d050_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.5,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d045_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.45,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d040_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.4,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d035_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.35,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d030_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.3,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d025_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.25,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d020_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.2,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d015_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.15,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d010_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.1,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d005_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.05,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
mash_filtered_d001_reduced <- subset(mash_filtered,mash_filtered$mash_distance < 0.001,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))


write.table(mash_filtered_d050_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d050_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d045_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d045_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d040_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d040_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d035_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d035_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d030_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d030_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d025_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d025_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d020_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d020_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d015_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d015_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d010_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d010_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d005_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d005_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
write.table(mash_filtered_d001_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_d001_pham_data.csv",sep=",",row.names = FALSE,quote=FALSE)
#
#
#
mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
mash_filtered_d010 <- subset(mash_filtered,mash_filtered$mash_group_d010_compare != "different")

par(mar=c(4,2,4,15))
stripchart(mash_filtered$modified_mash_distance ~ mash_filtered$mash_group_d010_compare,at=rank(tapply(mash_filtered$modified_mash_distance,mash_filtered$mash_group_d010_compare,median)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)

par(mar=c(4,2,4,15))
stripchart(mash_filtered$modified_mash_distance ~ mash_filtered$mash_group_d005_compare,at=rank(tapply(mash_filtered$modified_mash_distance,mash_filtered$mash_group_d005_compare,median)),las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)


mash_filtered$mash_group_d010_compare <- as.factor(mash_filtered$mash_group_d010_compare)


mash_filtered_d010_max_distances <- tapply(mash_filtered$mash_distance,mash_filtered$mash_group_d010_compare,max)


mash_filtered_d010 <- subset(mash_filtered,mash_filtered$mash_distance < 0.1)
mash_filtered_d010 <- subset(mash_filtered_d010,mash_filtered_d010$mash_group_d010_compare != "different")

mash_filtered_d010$mash_group_d010_compare <- as.factor(mash_filtered_d010$mash_group_d010_compare)
mash_filtered_d010_max_distances <- tapply(mash_filtered_d010$mash_distance,mash_filtered_d010$mash_group_d010_compare,max)
#
#
#












###Gene count disparity analysis

#Match gene counts per genome to pairwise comparison data
gene_count_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161213_phage_gene_count.csv",sep=",",header=TRUE)
names(gene_count_table) = c('phage_identifier','query_gene_count')
mash_table2_gene_count <- merge(mash_table2,gene_count_table,by.x="mash_query",by.y="phage_identifier")
names(gene_count_table) = c('phage_identifier','ref_gene_count')
mash_table2_gene_count <- merge(mash_table2_gene_count,gene_count_table,by.x="mash_reference",by.y="phage_identifier")
mash_table2_gene_count[mash_table2_gene_count == "Unspecified"] <- NA









###Compute gene flux differences
#Before subsetting data, be sure to cut all mash distance data into bins, so that bins are not affected by the size or distribution of each subset
#For the analysis, use only:
#1. dsDNA phages
#2. intra-cluster comparisons, so that all gene flux rates are among the same type of phages
#3. phages that have been annotated by SEA-PHAGES (phamerator_status = 'final')


type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type$mash_distance_bin <- cut(type$modified_mash_distance,breaks=10)

dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
dsDNA_cluster <- subset(dsDNA,dsDNA$phage_cluster_compare != "different")
dsDNA_finals <- subset(dsDNA_cluster,dsDNA_cluster$phamerator_status_compare == "final")
dsDNA_finals_reduced <- subset(dsDNA_finals,dsDNA_finals$phage_temperate_compare != "different")
dsDNA_finals_reduced$phage_temperate_compare <- factor(dsDNA_finals_reduced$phage_temperate_compare)

#BELOW: THESE COMMANDS NEED UPDATED. MASH_TABLE2 NOW CONTAINS THESE FIELDS
dsDNA_finals_reduced$size_diff_ave_percent <- (dsDNA_finals_reduced$size_diff_ref_percent + dsDNA_finals_reduced$size_diff_query_percent)/2
dsDNA_finals_reduced$gene_count_disparity_ave_percent <- ((dsDNA_finals_reduced$gene_count_disparity/dsDNA_finals_reduced$ref_gene_count) + (dsDNA_finals_reduced$gene_count_disparity/dsDNA_finals_reduced$query_gene_count))/2
#



#Subset data to plot
dsDNA_finals_temperate <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$phage_temperate_compare == "yes")
dsDNA_finals_lytic <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$phage_temperate_compare == "no")

dsDNA_finals_temperate_hgf <- subset(dsDNA_finals_temperate,dsDNA_finals_temperate$gene_flux_category == "high")
dsDNA_finals_temperate_lgf <- subset(dsDNA_finals_temperate,dsDNA_finals_temperate$gene_flux_category == "low")

dsDNA_finals_lytic_hgf <- subset(dsDNA_finals_lytic,dsDNA_finals_lytic$gene_flux_category == "high")
dsDNA_finals_lytic_lgf <- subset(dsDNA_finals_lytic,dsDNA_finals_lytic$gene_flux_category == "low")


#Scatter plot of each subset of data
par(mar=c(4,8,4,4))
plot(dsDNA_finals_temperate_hgf$modified_mash_distance,dsDNA_finals_temperate_hgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(dsDNA_finals_temperate_lgf$modified_mash_distance,dsDNA_finals_temperate_lgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(dsDNA_finals_lytic_hgf$modified_mash_distance,dsDNA_finals_lytic_hgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(dsDNA_finals_lytic_lgf$modified_mash_distance,dsDNA_finals_lytic_lgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






par(mar=c(4,8,4,4))
plot(dsDNA_finals_temperate_hgf$modified_mash_distance,dsDNA_finals_temperate_hgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,2,lty=2,lwd=3,col="grey")
par(new=TRUE)
plot(dsDNA_finals_temperate_lgf$modified_mash_distance,dsDNA_finals_temperate_lgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_finals_lytic$modified_mash_distance,dsDNA_finals_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")





#sliding window average
library(caTools)
dsDNA_finals_temp_hgf_mmdsort <- dsDNA_finals_temperate_hgf[order(dsDNA_finals_temperate_hgf$modified_mash_distance),]
dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity_runmean <- runmean(dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity,101)
dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity_runmed <- runmed(dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity,501,algorithm = "Stuetzle",endrule = "keep")
dsDNA_finals_temp_hgf_mmdsort$size_diff_runmean <- runmean(dsDNA_finals_temp_hgf_mmdsort$size_diff,101)
dsDNA_finals_temp_hgf_mmdsort$size_diff_ave_percent_runmean <- runmean(dsDNA_finals_temp_hgf_mmdsort$size_diff_ave_percent,101)
dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity_ave_percent_runmean <- runmean(dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity_ave_percent,101)


dsDNA_finals_temp_lgf_mmdsort <- dsDNA_finals_temperate_lgf[order(dsDNA_finals_temperate_lgf$modified_mash_distance),]
dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity_runmean <- runmean(dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity,101)
dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity_runmed <- runmed(dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity,501,algorithm = "Stuetzle",endrule = "keep")
dsDNA_finals_temp_lgf_mmdsort$size_diff_runmean <- runmean(dsDNA_finals_temp_lgf_mmdsort$size_diff,101)
dsDNA_finals_temp_lgf_mmdsort$size_diff_ave_percent_runmean <- runmean(dsDNA_finals_temp_lgf_mmdsort$size_diff_ave_percent,101)
dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity_ave_percent_runmean <- runmean(dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity_ave_percent,101)


dsDNA_finals_lytic_mmdsort <- dsDNA_finals_lytic[order(dsDNA_finals_lytic$modified_mash_distance),]
dsDNA_finals_lytic_mmdsort$gene_count_disparity_runmean <- runmean(dsDNA_finals_lytic_mmdsort$gene_count_disparity,101)
dsDNA_finals_lytic_mmdsort$gene_count_disparity_runmed <- runmed(dsDNA_finals_lytic_mmdsort$gene_count_disparity,501,algorithm = "Stuetzle",endrule = "keep")
dsDNA_finals_lytic_mmdsort$size_diff_runmean <- runmean(dsDNA_finals_lytic_mmdsort$size_diff,101)
dsDNA_finals_lytic_mmdsort$size_diff_ave_percent_runmean <- runmean(dsDNA_finals_lytic_mmdsort$size_diff_ave_percent,101)
dsDNA_finals_lytic_mmdsort$gene_count_disparity_ave_percent_runmean <- runmean(dsDNA_finals_lytic_mmdsort$gene_count_disparity_ave_percent,101)




par(mar=c(4,8,4,4))
plot(dsDNA_finals_lytic_mmdsort$modified_mash_distance,dsDNA_finals_lytic_mmdsort$gene_count_disparity_runmean,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_finals_temp_lgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity_runmean,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_finals_temp_hgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity_runmean,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(dsDNA_finals_lytic_mmdsort$modified_mash_distance,dsDNA_finals_lytic_mmdsort$size_diff_runmean,xlim=c(0,0.5),ylim=c(0,3500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_finals_temp_lgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_lgf_mmdsort$size_diff_runmean,xlim=c(0,0.5),ylim=c(0,3500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_finals_temp_hgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_hgf_mmdsort$size_diff_runmean,xlim=c(0,0.5),ylim=c(0,3500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(dsDNA_finals_lytic_mmdsort$modified_mash_distance,dsDNA_finals_lytic_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.08),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_finals_temp_lgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_lgf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.08),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_finals_temp_hgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_hgf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.08),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(dsDNA_finals_lytic_mmdsort$modified_mash_distance,dsDNA_finals_lytic_mmdsort$gene_count_disparity_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.17),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_finals_temp_lgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_lgf_mmdsort$gene_count_disparity_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.17),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_finals_temp_hgf_mmdsort$modified_mash_distance,dsDNA_finals_temp_hgf_mmdsort$gene_count_disparity_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.17),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")






#Now subset data by bins
dsDNA_finals_bin1 <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$modified_mash_distance > 0 & (dsDNA_finals_reduced$modified_mash_distance < 0.05 | dsDNA_finals_reduced$modified_mash_distance == 0.05))
dsDNA_finals_bin2 <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$modified_mash_distance > 0.05 & (dsDNA_finals_reduced$modified_mash_distance < 0.10 | dsDNA_finals_reduced$modified_mash_distance == 0.10))
dsDNA_finals_bin3 <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$modified_mash_distance > 0.1 & (dsDNA_finals_reduced$modified_mash_distance < 0.15 | dsDNA_finals_reduced$modified_mash_distance == 0.15))
dsDNA_finals_bin4 <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$modified_mash_distance > 0.15 & (dsDNA_finals_reduced$modified_mash_distance < 0.20 | dsDNA_finals_reduced$modified_mash_distance == 0.20))
dsDNA_finals_bin5 <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$modified_mash_distance > 0.2 & (dsDNA_finals_reduced$modified_mash_distance < 0.25 | dsDNA_finals_reduced$modified_mash_distance == 0.25))
dsDNA_finals_bin6 <- subset(dsDNA_finals_reduced,dsDNA_finals_reduced$modified_mash_distance > 0.25 & (dsDNA_finals_reduced$modified_mash_distance < 0.30 | dsDNA_finals_reduced$modified_mash_distance == 0.30))


par(mar=c(4,20,4,4))
boxplot(dsDNA_finals_bin1$gene_count_disparity ~ dsDNA_finals_bin1$gene_flux_category*dsDNA_finals_bin1$phage_temperate_compare,las=2)
par(mar=c(4,20,4,4))
boxplot(dsDNA_finals_bin2$gene_count_disparity ~ dsDNA_finals_bin2$gene_flux_category*dsDNA_finals_bin2$phage_temperate_compare,las=2)
par(mar=c(4,20,4,4))
boxplot(dsDNA_finals_bin3$gene_count_disparity ~ dsDNA_finals_bin3$gene_flux_category*dsDNA_finals_bin3$phage_temperate_compare,las=2)
par(mar=c(4,20,4,4))
boxplot(dsDNA_finals_bin4$gene_count_disparity ~ dsDNA_finals_bin4$gene_flux_category*dsDNA_finals_bin4$phage_temperate_compare,las=2)
par(mar=c(4,20,4,4))
boxplot(dsDNA_finals_bin5$gene_count_disparity ~ dsDNA_finals_bin5$gene_flux_category*dsDNA_finals_bin5$phage_temperate_compare,las=2)
par(mar=c(4,20,4,4))
boxplot(dsDNA_finals_bin6$gene_count_disparity ~ dsDNA_finals_bin6$gene_flux_category*dsDNA_finals_bin6$phage_temperate_compare,las=2)













###Gene flux analysis on all data, not just SEA-PHAGES final status phages
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type$mash_distance_bin <- cut(type$modified_mash_distance,breaks=10)
bacteria <- subset(type,type$host_superkingdom_compare == "Bacteria")
dsDNA <- subset(bacteria,bacteria$phage_viral_type_compare == "dsDNA")
dsDNA$phage_temperate_compare <- factor(dsDNA$phage_temperate_compare)
dsDNA$size_diff_ave_percent <- (dsDNA$size_diff_ref_percent + dsDNA$size_diff_query_percent)/2
dsDNA$gene_count_disparity_ave_percent <- ((dsDNA$gene_count_disparity/dsDNA$ref_gene_count) + (dsDNA$gene_count_disparity/dsDNA$query_gene_count))/2


#Subset data to plot
dsDNA_temperate <- subset(dsDNA,dsDNA$phage_temperate_compare == "yes")
dsDNA_lytic <- subset(dsDNA,dsDNA$phage_temperate_compare == "no")

dsDNA_temperate_hgf <- subset(dsDNA_temperate,dsDNA_temperate$gene_flux_category == "high")
dsDNA_temperate_lgf <- subset(dsDNA_temperate,dsDNA_temperate$gene_flux_category == "low")

dsDNA_lytic_hgf <- subset(dsDNA_lytic,dsDNA_lytic$gene_flux_category == "high")
dsDNA_lytic_lgf <- subset(dsDNA_lytic,dsDNA_lytic$gene_flux_category == "low")


#Scatter plot of each subset of data
par(mar=c(4,8,4,4))
plot(dsDNA_temperate_hgf$modified_mash_distance,dsDNA_temperate_hgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(dsDNA_temperate_lgf$modified_mash_distance,dsDNA_temperate_lgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(dsDNA_lytic_hgf$modified_mash_distance,dsDNA_lytic_hgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(dsDNA_lytic_lgf$modified_mash_distance,dsDNA_lytic_lgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






par(mar=c(4,8,4,4))
plot(dsDNA_temperate_hgf$modified_mash_distance,dsDNA_temperate_hgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,2,lty=2,lwd=3,col="grey")
par(new=TRUE)
plot(dsDNA_temperate_lgf$modified_mash_distance,dsDNA_temperate_lgf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_lytic$modified_mash_distance,dsDNA_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")





#sliding window average
library(caTools)
dsDNA_temp_hgf_mmdsort <- dsDNA_temperate_hgf[order(dsDNA_temperate_hgf$modified_mash_distance),]
dsDNA_temp_hgf_mmdsort$gene_count_disparity_runmean <- runmean(dsDNA_temp_hgf_mmdsort$gene_count_disparity,101)
dsDNA_temp_hgf_mmdsort$size_diff_runmean <- runmean(dsDNA_temp_hgf_mmdsort$size_diff,101)
dsDNA_temp_hgf_mmdsort$size_diff_ave_percent_runmean <- runmean(dsDNA_temp_hgf_mmdsort$size_diff_ave_percent,101)
dsDNA_temp_hgf_mmdsort$gene_count_disparity_ave_percent_runmean <- runmean(dsDNA_temp_hgf_mmdsort$gene_count_disparity_ave_percent,101)
dsDNA_temp_hgf_mmdsort$size_diff_max_percent_runmean <- runmean(dsDNA_temp_hgf_mmdsort$size_diff_max_percent,101)

dsDNA_temp_lgf_mmdsort <- dsDNA_temperate_lgf[order(dsDNA_temperate_lgf$modified_mash_distance),]
dsDNA_temp_lgf_mmdsort$gene_count_disparity_runmean <- runmean(dsDNA_temp_lgf_mmdsort$gene_count_disparity,101)
dsDNA_temp_lgf_mmdsort$size_diff_runmean <- runmean(dsDNA_temp_lgf_mmdsort$size_diff,101)
dsDNA_temp_lgf_mmdsort$size_diff_ave_percent_runmean <- runmean(dsDNA_temp_lgf_mmdsort$size_diff_ave_percent,101)
dsDNA_temp_lgf_mmdsort$gene_count_disparity_ave_percent_runmean <- runmean(dsDNA_temp_lgf_mmdsort$gene_count_disparity_ave_percent,101)
dsDNA_temp_lgf_mmdsort$size_diff_max_percent_runmean <- runmean(dsDNA_temp_lgf_mmdsort$size_diff_max_percent,101)


dsDNA_lytic_mmdsort <- dsDNA_lytic[order(dsDNA_lytic$modified_mash_distance),]
dsDNA_lytic_mmdsort$gene_count_disparity_runmean <- runmean(dsDNA_lytic_mmdsort$gene_count_disparity,101)
dsDNA_lytic_mmdsort$size_diff_runmean <- runmean(dsDNA_lytic_mmdsort$size_diff,101)
dsDNA_lytic_mmdsort$size_diff_ave_percent_runmean <- runmean(dsDNA_lytic_mmdsort$size_diff_ave_percent,101)
dsDNA_lytic_mmdsort$gene_count_disparity_ave_percent_runmean <- runmean(dsDNA_lytic_mmdsort$gene_count_disparity_ave_percent,101)
dsDNA_lytic_mmdsort$size_diff_max_percent_runmean <- runmean(dsDNA_lytic_mmdsort$size_diff_max_percent,101)



par(mar=c(4,8,4,4))
plot(dsDNA_lytic_mmdsort$modified_mash_distance,dsDNA_lytic_mmdsort$gene_count_disparity_runmean,xlim=c(0,0.5),ylim=c(0,60),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_temp_lgf_mmdsort$modified_mash_distance,dsDNA_temp_lgf_mmdsort$gene_count_disparity_runmean,xlim=c(0,0.5),ylim=c(0,60),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_temp_hgf_mmdsort$modified_mash_distance,dsDNA_temp_hgf_mmdsort$gene_count_disparity_runmean,xlim=c(0,0.5),ylim=c(0,60),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(dsDNA_lytic_mmdsort$modified_mash_distance,dsDNA_lytic_mmdsort$size_diff_runmean,xlim=c(0,0.5),ylim=c(0,25000),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_temp_lgf_mmdsort$modified_mash_distance,dsDNA_temp_lgf_mmdsort$size_diff_runmean,xlim=c(0,0.5),ylim=c(0,25000),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_temp_hgf_mmdsort$modified_mash_distance,dsDNA_temp_hgf_mmdsort$size_diff_runmean,xlim=c(0,0.5),ylim=c(0,25000),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(dsDNA_lytic_mmdsort$modified_mash_distance,dsDNA_lytic_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_temp_lgf_mmdsort$modified_mash_distance,dsDNA_temp_lgf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_temp_hgf_mmdsort$modified_mash_distance,dsDNA_temp_hgf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(dsDNA_lytic_mmdsort$modified_mash_distance,dsDNA_lytic_mmdsort$size_diff_max_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_temp_lgf_mmdsort$modified_mash_distance,dsDNA_temp_lgf_mmdsort$size_diff_max_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_temp_hgf_mmdsort$modified_mash_distance,dsDNA_temp_hgf_mmdsort$size_diff_max_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")


par(mar=c(4,8,4,4))
plot(dsDNA_lytic_mmdsort$modified_mash_distance,dsDNA_lytic_mmdsort$gene_count_disparity_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(dsDNA_temp_lgf_mmdsort$modified_mash_distance,dsDNA_temp_lgf_mmdsort$gene_count_disparity_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(dsDNA_temp_hgf_mmdsort$modified_mash_distance,dsDNA_temp_hgf_mmdsort$gene_count_disparity_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")






















###statistics

par(mar=c(4,8,4,4))
plot(mash_table2$modified_mash_distance,mash_table2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Lines to demarcate groups
abline(1.35,-2,lty=2,lwd=3,col="grey")
abline(0.4,-2,lty=2,lwd=3,col="grey")
abline(0,3.5,lty=2,lwd=3,col="grey")
abline(0.25,2,lty=2,lwd=3,col="grey")
#abline(v=0.166667,lty=2,lwd=3,col="grey")



lifestyle_data <- subset(mash_table2,is.na(mash_table2$phage_temperate_compare) == FALSE)
lifestyle_data$retain_data <- ifelse(lifestyle_data$pham_pham_dissimilarity > (lifestyle_data$modified_mash_distance * (-2) + 0.4) & lifestyle_data$pham_pham_dissimilarity < (lifestyle_data$modified_mash_distance * (-2) + 1.35),TRUE,FALSE) 
lifestyle_data$mode_part1 <- ifelse(lifestyle_data$modified_mash_distance < 0.1666667 & lifestyle_data$pham_pham_dissimilarity > (lifestyle_data$modified_mash_distance * 3.5),TRUE,FALSE)
lifestyle_data$mode_part2 <- ifelse(lifestyle_data$modified_mash_distance > 0.1666667 & lifestyle_data$pham_pham_dissimilarity > (lifestyle_data$modified_mash_distance * 2 + 0.25),TRUE,FALSE)
lifestyle_data$mode_part3 <- ifelse(lifestyle_data$mode_part1 == TRUE | lifestyle_data$mode_part2 ==TRUE,TRUE,FALSE)

lifestyle_data_retained <- subset(lifestyle_data,lifestyle_data$retain_data == TRUE)
lifestyle_data_retained$mash_reference <- factor(lifestyle_data_retained$mash_reference)
lifestyle_data_retained$mash_query <- factor(lifestyle_data_retained$mash_query)
lifestyle_data_retained_phages <- unique(c(levels(lifestyle_data_retained$mash_reference),levels(lifestyle_data_retained$mash_query)))
lifestyle_data_retained_phages <- factor(lifestyle_data_retained_phages)


lifestyle_data_retained_phages <- as.data.frame(unique(c(levels(lifestyle_data_retained$mash_reference),levels(lifestyle_data_retained$mash_query))))
names(lifestyle_data_retained_phages) = c("phage_identifier")
lifestyle_data_retained_phages <- merge(lifestyle_data_retained_phages,host_table,by.x="phage_identifier",by.y ="phage_identifier")









mode_log <- subset(lifestyle_data_retained,lifestyle_data_retained$mode_part3 == TRUE)
mode_linear <- subset(lifestyle_data_retained,lifestyle_data_retained$mode_part3 == FALSE)

nrow(mode_log)
nrow(mode_linear)

mode_log$mash_reference <- factor(mode_log$mash_reference)
mode_log$mash_query <- factor(mode_log$mash_query)
mode_log_phages <- unique(c(levels(mode_log$mash_reference),levels(mode_log$mash_query)))
mode_log_phages <- factor(mode_log_phages)

mode_log_phages_df <- as.data.frame(mode_log_phages)
names(mode_log_phages_df) = c("phage_identifier")
mode_log_phages_df <- merge(mode_log_phages_df,host_table,by.x="phage_identifier",by.y ="phage_identifier")








mode_linear$mash_reference <- factor(mode_linear$mash_reference)
mode_linear$mash_query <- factor(mode_linear$mash_query)
mode_linear_phages <- unique(c(levels(mode_linear$mash_reference),levels(mode_linear$mash_query)))
mode_linear_phages <- factor(mode_linear_phages)


mode_linear_phages_df <- as.data.frame(mode_linear_phages)
names(mode_linear_phages_df) = c("phage_identifier")
mode_linear_phages_df <- merge(mode_linear_phages_df,host_table,by.x="phage_identifier",by.y ="phage_identifier")



mode_log_temperate <- subset(mode_log,mode_log$phage_temperate_compare == "yes")
mode_log_temperate$mash_reference <- factor(mode_log_temperate$mash_reference)
mode_log_temperate$mash_query <- factor(mode_log_temperate$mash_query)
mode_log_temperate_phages <- unique(c(levels(mode_log_temperate$mash_reference),levels(mode_log_temperate$mash_query)))
mode_log_temperate_phages <- factor(mode_log_temperate_phages)

mode_log_lytic <- subset(mode_log,mode_log$phage_temperate_compare == "no")
mode_log_lytic$mash_reference <- factor(mode_log_lytic$mash_reference)
mode_log_lytic$mash_query <- factor(mode_log_lytic$mash_query)
mode_log_lytic_phages <- unique(c(levels(mode_log_lytic$mash_reference),levels(mode_log_lytic$mash_query)))
mode_log_lytic_phages <- factor(mode_log_lytic_phages)

mode_log_mixed <- subset(mode_log,mode_log$phage_temperate_compare == "different")
mode_log_mixed$mash_reference <- factor(mode_log_mixed$mash_reference)
mode_log_mixed$mash_query <- factor(mode_log_mixed$mash_query)
mode_log_mixed_phages <- unique(c(levels(mode_log_mixed$mash_reference),levels(mode_log_mixed$mash_query)))
mode_log_mixed_phages <- factor(mode_log_mixed_phages)




mode_linear_temperate <- subset(mode_linear,mode_linear$phage_temperate_compare == "yes")
mode_linear_temperate$mash_reference <- factor(mode_linear_temperate$mash_reference)
mode_linear_temperate$mash_query <- factor(mode_linear_temperate$mash_query)
mode_linear_temperate_phages <- unique(c(levels(mode_linear_temperate$mash_reference),levels(mode_linear_temperate$mash_query)))
mode_linear_temperate_phages <- factor(mode_linear_temperate_phages)

mode_linear_lytic <- subset(mode_linear,mode_linear$phage_temperate_compare == "no")
mode_linear_lytic$mash_reference <- factor(mode_linear_lytic$mash_reference)
mode_linear_lytic$mash_query <- factor(mode_linear_lytic$mash_query)
mode_linear_lytic_phages <- unique(c(levels(mode_linear_lytic$mash_reference),levels(mode_linear_lytic$mash_query)))
mode_linear_lytic_phages <- factor(mode_linear_lytic_phages)

mode_linear_mixed <- subset(mode_linear,mode_linear$phage_temperate_compare == "different")
mode_linear_mixed$mash_reference <- factor(mode_linear_mixed$mash_reference)
mode_linear_mixed$mash_query <- factor(mode_linear_mixed$mash_query)
mode_linear_mixed_phages <- unique(c(levels(mode_linear_mixed$mash_reference),levels(mode_linear_mixed$mash_query)))
mode_linear_mixed_phages <- factor(mode_linear_mixed_phages)




all_temperate <- subset(lifestyle_data_retained,lifestyle_data_retained$phage_temperate_compare == "yes")
all_temperate$mash_reference <- factor(all_temperate$mash_reference)
all_temperate$mash_query <- factor(all_temperate$mash_query)
all_temperate_phages <- unique(c(levels(all_temperate$mash_reference),levels(all_temperate$mash_query)))
all_temperate_phages <- factor(all_temperate_phages)


all_lytic <- subset(lifestyle_data_retained,lifestyle_data_retained$phage_temperate_compare == "no")
all_lytic$mash_reference <- factor(all_lytic$mash_reference)
all_lytic$mash_query <- factor(all_lytic$mash_query)
all_lytic_phages <- unique(c(levels(all_lytic$mash_reference),levels(all_lytic$mash_query)))
all_lytic_phages <- factor(all_lytic_phages)

all_mixed <- subset(lifestyle_data_retained,lifestyle_data_retained$phage_temperate_compare == "different")
all_mixed$mash_reference <- factor(all_mixed$mash_reference)
all_mixed$mash_query <- factor(all_mixed$mash_query)
all_mixed_phages <- unique(c(levels(all_mixed$mash_reference),levels(all_mixed$mash_query)))
all_mixed_phages <- factor(all_mixed_phages)

host_table$phage_temperate
all_mixed_phages_table <- subset(host_table,select = c("phage_identifier","phage_temperate"))






length(union(mode_log_temperate_phages,mode_linear_temperate_phages))
length(intersect(mode_log_temperate_phages,mode_linear_temperate_phages))
length(setdiff(mode_log_temperate_phages,mode_linear_temperate_phages))
length(setdiff(mode_linear_temperate_phages,mode_log_temperate_phages))


length(union(mode_log_lytic_phages,mode_linear_lytic_phages))
length(intersect(mode_log_lytic_phages,mode_linear_lytic_phages))
length(setdiff(mode_log_lytic_phages,mode_linear_lytic_phages))
length(setdiff(mode_linear_lytic_phages,mode_log_lytic_phages))


#based on # of comparisons
binom.test(3478,11961,p = 0.243,alternative = "two.sided")
binom.test(69,2639,p = 0.243,alternative = "two.sided")
binom.test(18,22,p = 0.243,alternative = "two.sided")


#based on # of unique phages involved
binom.test(318,398,p = 0.61,alternative = "two.sided")
binom.test(301,631,p = 0.61,alternative = "two.sided")

#scatter plots
par(mar=c(4,8,4,4))
plot(temperate$modified_mash_distance,temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lytic$modified_mash_distance,lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(lifestyle_data_retained$modified_mash_distance,lifestyle_data_retained$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mode_log$modified_mash_distance,mode_log$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mode_linear$modified_mash_distance,mode_linear$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




par(mar=c(4,8,4,4))
plot(mode_log_temperate$modified_mash_distance,mode_log_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(mode_log_lytic$modified_mash_distance,mode_log_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Fit best trend to high and low gene content flux patterns
bacteria <- subset(mash_table2,mash_table2$host_superkingdom_compare == 'Bacteria')
dsDNA <- subset(bacteria,bacteria$phage_viral_type_compare == 'dsDNA')
filtered <- subset(dsDNA,dsDNA$filter == TRUE)
lifestyle_data <- subset(filtered,is.na(filtered$phage_temperate_compare) == FALSE)
#lifestyle_data$retain_data <- ifelse(lifestyle_data$pham_pham_dissimilarity > (lifestyle_data$modified_mash_distance * (-2) + 0.4) & lifestyle_data$pham_pham_dissimilarity < (lifestyle_data$modified_mash_distance * (-2) + 1.35),TRUE,FALSE) 
lifestyle_data$mode_part1 <- ifelse(lifestyle_data$modified_mash_distance < 0.1666667 & lifestyle_data$pham_pham_dissimilarity > (lifestyle_data$modified_mash_distance * 3.5),TRUE,FALSE)
lifestyle_data$mode_part2 <- ifelse(lifestyle_data$modified_mash_distance > 0.1666667 & lifestyle_data$pham_pham_dissimilarity > (lifestyle_data$modified_mash_distance * 2 + 0.25),TRUE,FALSE)
lifestyle_data$mode_part3 <- ifelse(lifestyle_data$mode_part1 == TRUE | lifestyle_data$mode_part2 ==TRUE,TRUE,FALSE)

hgcf <- subset(lifestyle_data,lifestyle_data$mode_part3 == TRUE)
lgcf <- subset(lifestyle_data,lifestyle_data$mode_part3 == FALSE)

temperate_hgcf <- subset(hgcf,hgcf$phage_temperate_compare == "yes")
temperate_lgcf <- subset(lgcf,lgcf$phage_temperate_compare == "yes")
lytic <- subset(lifestyle_data,lifestyle_data$phage_temperate_compare == "no")

par(mar=c(4,8,4,4))
plot(hgcf$modified_mash_distance,hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(lgcf$modified_mash_distance,lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_hgcf$modified_mash_distance,temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(temperate_lgcf$modified_mash_distance,temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(lytic$modified_mash_distance,lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")
abline(lm(lytic$pham_pham_dissimilarity ~ lytic$modified_mash_distance))
abline(lm(temperate_lgcf$pham_pham_dissimilarity ~ temperate_lgcf$modified_mash_distance))
abline(lm(temperate_hgcf$pham_pham_dissimilarity ~ temperate_hgcf$modified_mash_distance))
lines(predict(temperate_hgcf_pm))

lytic_lm <- lm(lytic$pham_pham_dissimilarity ~ lytic$modified_mash_distance,data = lytic)
temperate_lgcf_lm <- lm(temperate_lgcf$pham_pham_dissimilarity ~ temperate_lgcf$modified_mash_distance,data = lytic)
temperate_hgcf_lm <- lm(temperate_hgcf$pham_pham_dissimilarity ~ temperate_hgcf$modified_mash_distance,data = lytic)

summary(lytic_lm)
summary(temperate_lgcf_lm)
summary(temperate_hgcf_lm)



temperate_hgcf_pm <- lm(temperate_hgcf$pham_pham_dissimilarity ~ poly(temperate_hgcf$modified_mash_distance,3))
temperate_hgcf_pm_predict <- predict(temperate_hgcf_pm,temperate_hgcf$modified_mash_distance)



temperate_hgcf_sorted <- temperate_hgcf[order(temperate_hgcf$modified_mash_distance),]

lines(temperate_hgcf[order(temperate_hgcf),],fitted(temperate_hgcf[order(temperate_hgcf),]))

lines(temperate_hgcf_sorted$modified_mash_distance,fitted(temperate_hgcf_pm))












###Investigate enterobacteria myoviridae

proteo <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$host_phylum_compare == "Proteobacteria")
proteo_sipho <- subset(proteo,proteo$phage_family_compare == "Siphoviridae")
proteo_myo <- subset(proteo,proteo$phage_family_compare == "Myoviridae")
proteo_d030_MG391 <- subset(proteo,proteo$mash_group_d030_compare == "Mash_Group_391")

proteo_d020_MG348 <- subset(proteo,proteo$mash_group_d020_compare == "Mash_Group_348")

Mash_Group_348

mash_groups_d030 <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$mash_group_d030_compare!="different")


par(mar=c(4,8,4,4))
plot(proteo_sipho$mash_distance,proteo_sipho$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=0.1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2)


par(mar=c(4,8,4,4))
plot(proteo_myo$mash_distance,proteo_myo$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=0.1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2)

par(mar=c(4,8,4,4))
plot(proteo_d020_MG348$mash_distance,proteo_d020_MG348$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=0.1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2)


par(mar=c(4,8,4,4))
hist(proteo_myo$mash_distance,breaks=((range(proteo_myo$mash_distance)[2]-range(proteo_myo$mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),col="black",cex.axis=2)








###Check Cluster M

cluster_m <- mash_table2
cluster_m$cluster_m_one_or_two <- ifelse(cluster_m$ref_phage_cluster == "M" | cluster_m$query_phage_cluster == "M",TRUE,FALSE)
#cluster_m$cluster_m_two <- ifelse(cluster_m$ref_phage_cluster == "M" & cluster_m$query_phage_cluster == "M",TRUE,FALSE)
cluster_m$cluster_m_two <- ifelse(cluster_m$phage_cluster_compare == "M",TRUE,FALSE)
cluster_m$cluster_m_one <- ifelse(cluster_m$cluster_m_one_or_two == TRUE & cluster_m$cluster_m_two == FALSE,TRUE,FALSE)

cluster_m_one_or_two <- subset(cluster_m,cluster_m$cluster_m_one_or_two == TRUE)
cluster_m_one <- subset(cluster_m,cluster_m$cluster_m_one == TRUE)
cluster_m_two <- subset(cluster_m,cluster_m$cluster_m_two == TRUE)


par(mar=c(4,8,4,4))
plot(cluster_m_one_or_two$modified_mash_distance,cluster_m_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_m_one$modified_mash_distance,cluster_m_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_m_two$modified_mash_distance,cluster_m_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")











###Check Patience data


mash_host_pham_pvalue_size_ani_groups$patience_one <- ifelse(mash_host_pham_pvalue_size_ani_groups$mash_reference == "patience__actino785" | mash_host_pham_pvalue_size_ani_groups$mash_query == "patience__actino785",TRUE,FALSE)



patience_one <- subset(mash_host_pham_pvalue_size_ani_groups,mash_host_pham_pvalue_size_ani_groups$patience_one == TRUE)

par(mar=c(4,8,4,4))
plot(patience_one$mash_distance,patience_one$pham_pham_distance,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2)





###Pie chart by host taxonomic rankings

pie(summary(mash_host_pham_pvalue_size_ani$host_phylum_compare),main="Phage diversity by host phylum",clockwise=TRUE)

pie(summary(host_table$host_phylum),main="Phage diversity by host phylum",clockwise=TRUE)
pie(summary(host_table$host_class),main="Phage diversity by host class",clockwise=TRUE)
pie(summary(host_table$host_order),main="Phage diversity by host order",clockwise=TRUE)
pie(summary(host_table$host_family),main="Phage diversity by host family",clockwise=TRUE)
pie(summary(host_table$host_genus),main="Phage diversity by host genus",clockwise=TRUE)

pie(c(retained_all,1-retained_all),labels=NA,col=c("black","white"),main="All comparisons retained",clockwise=TRUE)












#
#
#
###Plot number of interactions as p-value data is progressively filtered
filter_data_by_pvalue <- function(x){
  sum(kmer14_host_pham$kmer14_pvalue < x)
}

pvalue_vector <- seq(0,1,0.0001)
filtered_interactions_by_pvalue <- sapply(pvalue_vector,filter_data_by_pvalue)
plot(pvalue_vector,filtered_interactions_by_pvalue,pch=20,cex=0.1,xlim=c(0,1),xlab="pvalue distance",ylab="# of interactions",main="Total # of interactions vs pvalue",type="l")
plot(pvalue_vector,filtered_interactions_by_pvalue,pch=20,cex=0.1,xlim=c(0,1e-2),xlab="pvalue distance",ylab="# of interactions",main="Total # of interactions vs pvalue",type="l")



#Plot number of interactions as data is progressively filtered
filter_data_by_distance <- function(x){
  sum(kmer14_host_pham_pvalue$kmer14_distance < x)
}

distance_vector <- seq(0,1,0.0001)
filtered_interactions <- sapply(distance_vector,filter_data_by_distance)
plot(distance_vector,filtered_interactions,pch=20,cex=0.1,xlim=c(0,1),xlab="Mash distance",ylab="# of interactions",main="Total # of interactions vs Mash distance",type="l")
abline(v=0.2,lty=2)
plot(distance_vector,filtered_interactions,pch=20,cex=0.1,ylim=c(0,8e4),xlim=c(0,0.35),xlab="Mash distance",ylab="# of interactions",main="Total # of interactions vs Mash distance",type="l")
abline(v=0.2,lty=2)
#
#
#






###Host taxa diversity
library(vegan)
library(data.table)

#Rows = margin 1 = host phyla
#Columns = margin 2 = class, order, family, genus
host_class_diversity_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_class_frequency.csv",sep=",",row.names=1,header=TRUE)
host_order_diversity_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_order_frequency.csv",sep=",",row.names=1,header=TRUE)
host_family_diversity_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_family_frequency.csv",sep=",",row.names=1,header=TRUE)
host_genus_diversity_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_genus_frequency.csv",sep=",",row.names=1,header=TRUE)


host_class_richness <- specnumber(host_class_diversity_table)
host_order_richness <- specnumber(host_order_diversity_table)
host_family_richness <- specnumber(host_family_diversity_table)
host_genus_richness <- specnumber(host_genus_diversity_table)


host_class_diversity <- diversity(host_class_diversity_table)
host_order_diversity <- diversity(host_order_diversity_table)
host_family_diversity <- diversity(host_family_diversity_table)
host_genus_diversity <- diversity(host_genus_diversity_table)


host_class_evenness <- host_class_diversity/log(host_class_richness)
host_order_evenness <- host_order_diversity/log(host_order_richness)
host_family_evenness <- host_family_diversity/log(host_family_richness)
host_genus_evenness <- host_genus_diversity/log(host_genus_richness)




host_class_richness_df <- as.data.frame(host_class_richness)
host_order_richness_df <- as.data.frame(host_order_richness)
host_family_richness_df <- as.data.frame(host_family_richness)
host_genus_richness_df <- as.data.frame(host_genus_richness)

host_class_diversity_df <- as.data.frame(host_class_diversity)
host_order_diversity_df <- as.data.frame(host_order_diversity)
host_family_diversity_df <- as.data.frame(host_family_diversity)
host_genus_diversity_df <- as.data.frame(host_genus_diversity)


host_class_evenness_df <- as.data.frame(host_class_evenness)
host_order_evenness_df <- as.data.frame(host_order_evenness)
host_family_evenness_df <- as.data.frame(host_family_evenness)
host_genus_evenness_df <- as.data.frame(host_genus_evenness)








setDT(host_class_richness_df,keep.rownames = TRUE)
setDT(host_order_richness_df,keep.rownames = TRUE)
setDT(host_family_richness_df,keep.rownames = TRUE)
setDT(host_genus_richness_df,keep.rownames = TRUE)

setDT(host_class_diversity_df,keep.rownames = TRUE)
setDT(host_order_diversity_df,keep.rownames = TRUE)
setDT(host_family_diversity_df,keep.rownames = TRUE)
setDT(host_genus_diversity_df,keep.rownames = TRUE)

setDT(host_class_evenness_df,keep.rownames = TRUE)
setDT(host_order_evenness_df,keep.rownames = TRUE)
setDT(host_family_evenness_df,keep.rownames = TRUE)
setDT(host_genus_evenness_df,keep.rownames = TRUE)


host_richness_summary <- merge(host_class_richness_df,host_order_richness_df,by.x="rn",by.y="rn")
host_richness_summary <- merge(host_richness_summary,host_family_richness_df,by.x="rn",by.y="rn")
host_richness_summary <- merge(host_richness_summary,host_genus_richness_df,by.x="rn",by.y="rn")

host_diversity_summary <- merge(host_class_diversity_df,host_order_diversity_df,by.x="rn",by.y="rn")
host_diversity_summary <- merge(host_diversity_summary,host_family_diversity_df,by.x="rn",by.y="rn")
host_diversity_summary <- merge(host_diversity_summary,host_genus_diversity_df,by.x="rn",by.y="rn")

host_evenness_summary <- merge(host_class_evenness_df,host_order_evenness_df,by.x="rn",by.y="rn")
host_evenness_summary <- merge(host_evenness_summary,host_family_evenness_df,by.x="rn",by.y="rn")
host_evenness_summary <- merge(host_evenness_summary,host_genus_evenness_df,by.x="rn",by.y="rn")


write.table(host_richness_summary,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_rank_richness.csv",sep=",")
write.table(host_diversity_summary,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_rank_diversity.csv",sep=",")
write.table(host_evenness_summary,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_rank_evenness.csv",sep=",")



#write.table(host_class_richness,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_class_richness.csv",sep=",")
#write.table(host_order_richness,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_order_richness.csv",sep=",")
#write.table(host_family_richness,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_family_richness.csv",sep=",")
#write.table(host_genus_richness,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/host_diversity/host_genus_richness.csv",sep=",")






host_table_actino <- subset(host_table,host_table$host_phylum == "Actinobacteria")
host_table_bacter <- subset(host_table,host_table$host_phylum == "Bacteroidetes")
host_table_cyano <- subset(host_table,host_table$host_phylum == "Cyanobacteria")
host_table_firm <- subset(host_table,host_table$host_phylum == "Firmicutes")
host_table_proteo <- subset(host_table,host_table$host_phylum == "Proteobacteria")



barplot(host_table$host_genus ~ host_table$host_phylum)












#
#
#
###Code to analyze biodiversity

library(vegan)
library(data.table)





#Rows = margin 1 = host phyla
#Columns = margin 2 = Mash Groups (species)

community_100 <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_matrix_d100.csv",sep=",",row.names=1,header=TRUE)
community_040 <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_matrix_d040.csv",sep=",",row.names=1,header=TRUE)
community_030 <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_matrix_d030.csv",sep=",",row.names=1,header=TRUE)
community_020 <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_matrix_d020.csv",sep=",",row.names=1,header=TRUE)
community_010 <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_matrix_d010.csv",sep=",",row.names=1,header=TRUE)
community_005 <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_matrix_d005.csv",sep=",",row.names=1,header=TRUE)

community_richness_100 <- specnumber(community_100)
community_richness_040 <- specnumber(community_040)
community_richness_030 <- specnumber(community_030)
community_richness_020 <- specnumber(community_020)
community_richness_010 <- specnumber(community_010)
community_richness_005 <- specnumber(community_005)




community_diversity_100 <- diversity(community_100)
community_diversity_040 <- diversity(community_040)
community_diversity_030 <- diversity(community_030)
community_diversity_020 <- diversity(community_020)
community_diversity_010 <- diversity(community_010)
community_diversity_005 <- diversity(community_005)


community_evenness_100 <- community_diversity_100/log(community_richness_100)
community_evenness_040 <- community_diversity_040/log(community_richness_040)
community_evenness_030 <- community_diversity_030/log(community_richness_030)
community_evenness_020 <- community_diversity_020/log(community_richness_020)
community_evenness_010 <- community_diversity_010/log(community_richness_010)
community_evenness_005 <- community_diversity_005/log(community_richness_005)





community_richness_df_100 <- as.data.frame(community_richness_100)
community_richness_df_040 <- as.data.frame(community_richness_040)
community_richness_df_030 <- as.data.frame(community_richness_030)
community_richness_df_020 <- as.data.frame(community_richness_020)
community_richness_df_010 <- as.data.frame(community_richness_010)
community_richness_df_005 <- as.data.frame(community_richness_005)


community_diversity_df_100 <- as.data.frame(community_diversity_100)
community_diversity_df_040 <- as.data.frame(community_diversity_040)
community_diversity_df_030 <- as.data.frame(community_diversity_030)
community_diversity_df_020 <- as.data.frame(community_diversity_020)
community_diversity_df_010 <- as.data.frame(community_diversity_010)
community_diversity_df_005 <- as.data.frame(community_diversity_005)


community_evenness_df_100 <- as.data.frame(community_evenness_100)
community_evenness_df_040 <- as.data.frame(community_evenness_040)
community_evenness_df_030 <- as.data.frame(community_evenness_030)
community_evenness_df_020 <- as.data.frame(community_evenness_020)
community_evenness_df_010 <- as.data.frame(community_evenness_010)
community_evenness_df_005 <- as.data.frame(community_evenness_005)




setDT(community_richness_df_100,keep.rownames = TRUE)
setDT(community_richness_df_040,keep.rownames = TRUE)
setDT(community_richness_df_030,keep.rownames = TRUE)
setDT(community_richness_df_020,keep.rownames = TRUE)
setDT(community_richness_df_010,keep.rownames = TRUE)
setDT(community_richness_df_005,keep.rownames = TRUE)


setDT(community_diversity_df_100,keep.rownames = TRUE)
setDT(community_diversity_df_040,keep.rownames = TRUE)
setDT(community_diversity_df_030,keep.rownames = TRUE)
setDT(community_diversity_df_020,keep.rownames = TRUE)
setDT(community_diversity_df_010,keep.rownames = TRUE)
setDT(community_diversity_df_005,keep.rownames = TRUE)


setDT(community_evenness_df_100,keep.rownames = TRUE)
setDT(community_evenness_df_040,keep.rownames = TRUE)
setDT(community_evenness_df_030,keep.rownames = TRUE)
setDT(community_evenness_df_020,keep.rownames = TRUE)
setDT(community_evenness_df_010,keep.rownames = TRUE)
setDT(community_evenness_df_005,keep.rownames = TRUE)



community_metrics_df_ALL <- merge(community_diversity_df_005,community_diversity_df_010,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_diversity_df_020,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_diversity_df_030,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_diversity_df_040,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_diversity_df_100,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_richness_df_005,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_richness_df_010,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_richness_df_020,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_richness_df_030,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_richness_df_040,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_richness_df_100,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_evenness_df_005,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_evenness_df_010,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_evenness_df_020,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_evenness_df_030,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_evenness_df_040,by.x="rn",by.y="rn")
community_metrics_df_ALL <- merge(community_metrics_df_ALL,community_evenness_df_100,by.x="rn",by.y="rn")



write.table(community_metrics_df_ALL,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/community_diversity/community_metrics.csv",sep=",")













###Cluster-specific colored analysis
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
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")

dev.off()
plot_cluster_specific_profiles(cluster_actino,"F")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"J")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"N")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"I")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"P")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"K")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"AO")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"B")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"A")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"BU")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"CZ")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"E")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"G")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"L")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"AO")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"H")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"BD")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"AN")


cluster_g_ani <- subset(cluster_actino,cluster_actino$phage_cluster_compare == 'G')
cluster_j_ani <- subset(cluster_actino,cluster_actino$phage_cluster_compare == 'J')
cluster_l_ani <- subset(cluster_actino,cluster_actino$phage_cluster_compare == 'L')
cluster_n_ani <- subset(cluster_actino,cluster_actino$phage_cluster_compare == 'N')

par(mar=c(4,8,4,4))
plot(cluster_g_ani$ani_distance,cluster_g_ani$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_j_ani$ani_distance,cluster_j_ani$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_l_ani$ani_distance,cluster_l_ani$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_n_ani$ani_distance,cluster_n_ani$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")






###Check to see how related high gcf temperate, low gcf temperate, and lytic phages are between actino hosts

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster_actino <- subset(type_dsDNA,type_dsDNA$phage_cluster_source_compare == "actino")
cluster_actino_diff_host <- subset(cluster_actino,cluster_actino$host_genus_compare == "different")
cluster_actino_diff_host_temperate_high_gcf <- subset(cluster_actino_diff_host,cluster_actino_diff_host$gene_flux_category == "high" & cluster_actino_diff_host$phage_temperate_compare == "yes")
cluster_actino_diff_host_temperate_low_gcf <- subset(cluster_actino_diff_host,cluster_actino_diff_host$gene_flux_category == "low" & cluster_actino_diff_host$phage_temperate_compare == "yes")
cluster_actino_diff_host_lytic <- subset(cluster_actino_diff_host,cluster_actino_diff_host$phage_temperate_compare == "no")
cluster_actino_diff_host_temperate_lytic <- subset(cluster_actino_diff_host,cluster_actino_diff_host$phage_temperate_compare == "different")
cluster_actino_diff_host_lifestyle_unspecified <- subset(cluster_actino_diff_host,is.na(cluster_actino_diff_host$phage_temperate_compare))

nrow(cluster_actino_diff_host) -
  (nrow(cluster_actino_diff_host_temperate_high_gcf) +
     nrow(cluster_actino_diff_host_temperate_low_gcf) +
     nrow(cluster_actino_diff_host_lytic) +
     nrow(cluster_actino_diff_host_temperate_lytic) +
     nrow(cluster_actino_diff_host_lifestyle_unspecified))




par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host_temperate_low_gcf$modified_mash_distance,cluster_actino_diff_host_temperate_low_gcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(cluster_actino_diff_host_temperate_lytic$modified_mash_distance,cluster_actino_diff_host_temperate_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(cluster_actino_diff_host_lytic$modified_mash_distance,cluster_actino_diff_host_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host_lifestyle_unspecified$modified_mash_distance,cluster_actino_diff_host_lifestyle_unspecified$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(cluster_actino_diff_host_temperate_low_gcf$modified_mash_distance,cluster_actino_diff_host_temperate_low_gcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(cluster_actino_diff_host_temperate_lytic$modified_mash_distance,cluster_actino_diff_host_temperate_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(cluster_actino_diff_host_lytic$modified_mash_distance,cluster_actino_diff_host_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host_temperate_high_gcf$modified_mash_distance,cluster_actino_diff_host_temperate_high_gcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host_temperate_lytic$modified_mash_distance,cluster_actino_diff_host_temperate_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host_lytic$modified_mash_distance,cluster_actino_diff_host_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host_temperate_low_gcf$modified_mash_distance,cluster_actino_diff_host_temperate_low_gcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")




par(mar=c(4,8,4,4))
plot(cluster_actino_diff_host$modified_mash_distance,cluster_actino_diff_host$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")





















#
#
#
#
#
#
###Code to analyze Arthrobacter phages
#20161108
#Revised 20170120



#If needed, reduce the metadata table to only those phages that have been clustered by us 
actino_phages <- subset(host_group_table,host_group_table$phage_cluster_source == 'actino')
actino_phages$host_genus <- factor(actino_phages$host_genus)


arthro_analysis <- subset(mash_table2,mash_table2$phage_cluster_source_compare=="actino")

arthro_analysis$arthro_one_or_two <- ifelse(arthro_analysis$ref_host_genus == "Arthrobacter" | arthro_analysis$query_host_genus == "Arthrobacter",TRUE,FALSE)
arthro_analysis$arthro_neither <- ifelse(arthro_analysis$ref_host_genus == "Arthrobacter" | arthro_analysis$query_host_genus == "Arthrobacter",FALSE,TRUE)
arthro_analysis$arthro_both <- ifelse(arthro_analysis$ref_host_genus == "Arthrobacter" & arthro_analysis$query_host_genus == "Arthrobacter",TRUE,FALSE)
arthro_analysis$arthro_one <- ifelse(arthro_analysis$arthro_both == TRUE | arthro_analysis$arthro_neither == TRUE,FALSE,TRUE)

arthro_analysis_one <- subset(arthro_analysis,arthro_analysis$arthro_one == TRUE)
arthro_analysis_both <- subset(arthro_analysis,arthro_analysis$arthro_both == TRUE)
arthro_analysis_neither <- subset(arthro_analysis,arthro_analysis$arthro_neither == TRUE)

arthro_analysis_intercluster <- subset(arthro_analysis_both,arthro_analysis_both$phage_cluster_compare == 'different')
arthro_analysis_intracluster <- subset(arthro_analysis_both,arthro_analysis_both$phage_cluster_compare != 'different')


par(mar=c(4,8,4,4))
plot(arthro_analysis_one$modified_mash_distance,arthro_analysis_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="grey")
par(new=TRUE)
plot(arthro_analysis_intercluster$modified_mash_distance,arthro_analysis_intercluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(arthro_analysis_intracluster$modified_mash_distance,arthro_analysis_intracluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=2,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")



#hist(arthro_analysis_one$mash_distance,breaks=((range(arthro_analysis_one$mash_distance)[2]-range(arthro_analysis_one$mash_distance)[1]) * 100),xlab=NULL,ylab="Number of same comparisons",main="Merged2333 s25000k15p1e10size arthro one",xlim=c(0,0.5),col="black")
#hist(arthro_analysis_both$mash_distance,breaks=((range(arthro_analysis_both$mash_distance)[2]-range(arthro_analysis_both$mash_distance)[1]) * 100),xlab=NULL,ylab="Number of same comparisons",main="Merged2333 s25000k15p1e10size arthro both",xlim=c(0,0.5),col="black")
#plot(arthro_analysis_one$mash_distance,arthro_analysis_one$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size arthro one",xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1)
#plot(arthro_analysis_both$mash_distance,arthro_analysis_both$pham_pham_distance,xlab="Mash distance",ylab="Pham distance",main="Phams Shared vs Mash s25000k15p1e10size arthro both",xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1)
#
#
#













#
#
###Cluster-level pham conservation analysis

myco_cluster_pham_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161209_myco_cluster_pham_conservation.csv",sep=",",header=TRUE)
myco_cluster_pham_matrix <- as.matrix(myco_cluster_pham_table)

myco_cluster_pham_matrix_small <- subset(myco_cluster_pham_matrix,select = c('A','A1','F'))

library(gplots)
library(Heatplus)
library(ggplot2)

cscale <- colorRampPalette(c("blue","green"),space="rgb")(100)

heatmap.2(myco_cluster_pham_matrix,col=cscale,trace=c("none"))
#
#
#


#
#
#
###Analysis of shared phams between clusters
#Of the phams from a cluster that are found to be shared between other clusters, what percentage of those shared phams are found in each specific other cluster?

myco_cluster_pham_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161213_myco_cluster_shared_phams2.csv",sep=",",row.names=1,header=TRUE)
myco_cluster_pham_matrix <- as.matrix(myco_cluster_pham_table)

myco_cluster_pham_matrix_small <- subset(myco_cluster_pham_matrix,select = c('A','A1','F'))

library(gplots)
library(Heatplus)
library(ggplot2)

cscale <- colorRampPalette(c("blue","green"),space="rgb")(100)

cscale <- colorRampPalette(c("white","blue"),space="rgb")(100)

heatmap.2(myco_cluster_pham_matrix,col=cscale,trace=c("none"))
#
#
#








###Toxic phage analysis

mash_table2_toxic <- mash_table2
mash_table2_toxic$toxic_one_or_two <- ifelse(mash_table2_toxic$ref_toxic == "yes" | mash_table2_toxic$query_toxic == "yes",TRUE,FALSE)
mash_table2_toxic$toxic_both <- ifelse(mash_table2_toxic$ref_toxic == "yes" & mash_table2_toxic$query_toxic == "yes",TRUE,FALSE)


toxic_one_or_two <- subset(mash_table2_toxic,mash_table2_toxic$toxic_one_or_two == TRUE)
toxic_both <- subset(mash_table2_toxic,mash_table2_toxic$toxic_both == TRUE)
#toxic_neither <- subset(mash_table2_toxic,mash_table2_toxic$toxic_neither == TRUE)
#toxic_one <- subset(mash_table2_toxic,mash_table2_toxic$toxic_one == TRUE)


par(mar=c(4,8,4,4))
plot(toxic_both$modified_mash_distance,toxic_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(toxic_one_or_two$modified_mash_distance,toxic_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(mash_table2_tener$modified_mash_distance,mash_table2_tener$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


















#
#
#
###Code to analyze shared pham diversity using Shannon Index

library(vegan)
library(data.table)





#Rows = margin 1 = host phyla
#Columns = margin 2 = Mash Groups (species)

shared_pham_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161213_myco_cluster_shared_pham_diversity_input.csv",sep=",",row.names=1,header=TRUE)

shared_pham_table_diversity <- diversity(shared_pham_table)
shared_pham_table_diversity_df <- as.data.frame(shared_pham_table_diversity)
setDT(shared_pham_table_diversity_df,keep.rownames = TRUE)

write.table(shared_pham_table_diversity_df,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161213_myco_cluster_shared_pham_diversity_output.csv",sep=",")









###Empty plot
empty_data <- subset(mash_table2,mash_table2$phage_cluster_compare == "empty")

par(mar=c(4,8,4,4))
plot(empty_data$modified_mash_distance,empty_data$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)




###Gordonia analysis


actino <- subset(mash_table2,mash_table2$phage_cluster_source_compare == 'actino')
actino$gordonia_one_or_two <- ifelse(actino$ref_host_genus == 'Gordonia' | actino$query_host_genus == 'Gordonia',TRUE,FALSE)
actino$gordonia_both <- ifelse(actino$host_genus_compare == 'Gordonia',TRUE,FALSE)
actino$gordonia_one <- ifelse(actino$gordonia_one_or_two == TRUE & actino$gordonia_both == FALSE,TRUE,FALSE)


gordonia_both_df <- subset(actino,actino$gordonia_both == TRUE)
gordonia_one_df <- subset(actino,actino$gordonia_one == TRUE)
gordonia_both_intra_cluster <- subset(gordonia_both_df,gordonia_both_df$phage_cluster_compare != 'different')
gordonia_both_inter_cluster <- subset(gordonia_both_df,gordonia_both_df$phage_cluster_compare == 'different')
gordonia_both_intra_subcluster <- subset(gordonia_both_df,gordonia_both_df$phage_subcluster_compare != 'different')
gordonia_both_inter_subcluster <- subset(gordonia_both_df,gordonia_both_df$phage_subcluster_compare == 'different')


write.table(gordonia_both_df,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/gordonia_comparisons.csv",sep=",",row.names = FALSE,quote=FALSE)





par(mar=c(4,8,4,4))
plot(gordonia_both_df$modified_mash_distance,gordonia_both_df$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gordonia_one_df$modified_mash_distance,gordonia_one_df$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gordonia_both_intra_cluster$modified_mash_distance,gordonia_both_intra_cluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gordonia_both_inter_cluster$modified_mash_distance,gordonia_both_inter_cluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gordonia_both_intra_subcluster$modified_mash_distance,gordonia_both_intra_subcluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gordonia_both_inter_subcluster$modified_mash_distance,gordonia_both_inter_subcluster$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






actino785_data <- subset(mash_table2,mash_table2$phage_cluster_source_compare == "actino")
actino785_data_same <- subset(actino785_data,actino785_data$phage_cluster_compare != "different")
actino785_data_diff <- subset(actino785_data,actino785_data$phage_cluster_compare == "different")
actino785_data_diff$query_singleton <- grepl("^Singleton",actino785_data_diff$query_phage_cluster)
actino785_data_diff$ref_singleton <- grepl("^Singleton",actino785_data_diff$ref_phage_cluster)
actino785_data_diff$singleton_one_or_two <- ifelse(actino785_data_diff$ref_singleton == TRUE | actino785_data_diff$query_singleton == TRUE,TRUE,FALSE)
actino785_data_diff$singleton_one <- ifelse(actino785_data_diff$singleton_both == TRUE | actino785_data_diff$singleton_neither == TRUE,FALSE,TRUE)
actino785_singleton_one_or_two <- subset(actino785_data_diff,actino785_data_diff$singleton_one_or_two == TRUE)
gordonia_singleton <- subset(actino785_singleton_one_or_two,actino785_singleton_one_or_two$host_genus_compare == 'Gordonia')

par(mar=c(4,8,4,4))
plot(gordonia_singleton$modified_mash_distance,gordonia_singleton$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






#reduced gene-specific data for welkin
#from table generated in the gene-specific analysis section
welkin <- gene_specific_mash_table
welkin <- subset(welkin,welkin$phage_cluster_source_compare == 'actino')
welkin <- subset(welkin,select=c("mash_reference","mash_query","modified_mash_distance",
                                 "query_host_genus","query_phage_order","query_phage_family","query_phage_cluster","query_phage_subcluster","query_size","query_gene_count",
                                 "ref_host_genus","ref_phage_order","ref_phage_family","ref_phage_cluster","ref_phage_subcluster","ref_size","ref_gene_count",
                                 "pham_phage1_number_unshared_phams","pham_phage2_number_unshared_phams","pham_number_shared_phams","pham_pham_dissimilarity",
                                 "gsm_ref_shared_total_size","gsm_ref_unshared_total_size","gsm_ref_shared_ave_size","gsm_ref_unshared_ave_size",
                                 "gsm_query_shared_total_size","gsm_query_unshared_total_size","gsm_query_shared_ave_size","gsm_query_unshared_ave_size"
))


write.table(welkin,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/welkin_actino_comparisons.csv",sep=",",row.names = FALSE,quote=FALSE)




































###ANI vs Pham Dissimilarity data
#Import ANI data from the 79 genome optimization test to show ANI vs Pham Distance
#Data contains complete matrix of 79 x 79 comparisons, including self comparisons and duplicate (reciprocal) comparisons
ani79_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161006_set2_79_ani_data.csv",sep=",",header=TRUE)

names(ani79_data) <- c("ani79_ref_query","ani79_ref_phage_identifier","ani79_query_phage_identifier","ani79_ani_distance")

#Merge tables
#The main mash table contains non-redundant comparisona and no self comparisons
ani79_analysis <- merge(mash_table2,ani79_data,by.x="mash_ref_query",by.y="ani79_ref_query")

par(mar=c(4,8,4,4))
plot(ani79_analysis$modified_mash_distance,ani79_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(ani79_analysis$ani79_ani_distance,ani79_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###Scatter plot of all comparisons involving lambda

#Phage identifier for lambda:lambda__nc_001416
#Copy main data table then determine which comparisons involve lambda. Since no self comparisons are present in the dataset, only need to compute if there's 'one' lambda, instead of 'one_or_two' or 'both'
lambda_figure <- mash_table2
lambda_figure$lambda_one <- ifelse(lambda_figure$mash_reference == 'lambda__nc_001416' | lambda_figure$mash_query == 'lambda__nc_001416',TRUE,FALSE)
lambda_comparisons <- subset(lambda_figure,lambda_figure$lambda_one == TRUE)


par(mar=c(4,8,4,4))
plot(lambda_comparisons$modified_mash_distance,lambda_comparisons$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")














###Assess genome size disparity filter

genome_size <- mash_table2
genome_size_exceeds_pvalue <- subset(genome_size,genome_size$mash_pvalue >= 1e-10)
genome_size_within_pvalue <- subset(genome_size,genome_size$mash_pvalue < 1e-10)
genome_size_within_pvalue_exceeds_size_diff <- subset(genome_size_within_pvalue,genome_size_within_pvalue$size_diff_max_percent >= 1)
genome_size_within_pvalue_exceeds_size_diff_bacteria_dsDNA <- subset(genome_size_within_pvalue_exceeds_size_diff,genome_size_within_pvalue_exceeds_size_diff$host_superkingdom_compare == 'Bacteria')
genome_size_within_pvalue_exceeds_size_diff_bacteria_dsDNA <- subset(genome_size_within_pvalue_exceeds_size_diff_bacteria_dsDNA,genome_size_within_pvalue_exceeds_size_diff_bacteria_dsDNA$phage_viral_type_compare == 'dsDNA')

#mash_table2$filter <- ifelse(mash_table2$mash_pvalue < 1e-10 & mash_table2$size_diff_max_percent < 1,TRUE,FALSE)

par(mar=c(4,8,4,4))
plot(genome_size_within_pvalue_exceeds_size_diff_bacteria_dsDNA$mash_distance,genome_size_within_pvalue_exceeds_size_diff_bacteria_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temp$mash_distance,temp$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

















###Export data for gene-specific mash analysis and/or cytoscape networks
#Only investigate comparisons that are within cluster boundaries
mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
bacteria_dsDNA <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')

#Cluster-boundary dataset
bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.42 & bacteria_dsDNA$pham_pham_dissimilarity < 0.89)

#Subcluster-boundary dataset
bacteria_dsDNA_nuc020_gene062 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.20 & bacteria_dsDNA$pham_pham_dissimilarity < 0.62)


#Plot the data to make sure it contains only what I want
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA_nuc042_gene089$modified_mash_distance,bacteria_dsDNA_nuc042_gene089$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Now reduce the data table to only the columns needed for export, and export the data
bacteria_dsDNA_nuc042_gene089_reduced <- subset(bacteria_dsDNA_nuc042_gene089,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
write.table(bacteria_dsDNA_nuc042_gene089_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_bacteria_dsDNA_nuc042_gene089_data.csv",sep=",",row.names = FALSE,quote=FALSE)

bacteria_dsDNA_nuc020_gene062_reduced <- subset(bacteria_dsDNA_nuc020_gene062,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
write.table(bacteria_dsDNA_nuc020_gene062_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/s25k_k15_bacteria_dsDNA_nuc020_gene062_data.csv",sep=",",row.names = FALSE,quote=FALSE)





###Gene-specific mash analysis
#After running gene-specific mash analysis, subset out only the comparisons used in the gsm analysis

mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
bacteria_dsDNA <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')
bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.42 & bacteria_dsDNA$pham_pham_dissimilarity < 0.89)

#Rename the reduced data to the generic gsm analysis table name
gene_specific_mash_table <- bacteria_dsDNA_nuc042_gene089


#Now import the gsm data
gene_specific_mash_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170109_gene_specific_mash_analysis.csv",sep=",",header=TRUE)

#Original column headers after import, but the names got truncated:
#names(gene_specific_mash_data) = c("phage1_phage2","phage1","phage2","phage1_phage2...shared.phams","phage1_phage2.gene.content.dissimilarity..general.index.","phage1_phage2.gene.content.dissimilarity..jaccard.index.","phage1_phage2.all.genes.mash.distance","phage1_phage2.all.genes.mash.p.value",
#                                    "phage1_phage2.all.genes.kmer.count","phage1_phage2.shared.genes.mash.distance","phage1_phage2.shared.genes.mash.p.value","phage1_phage2.shared.genes.kmer.count","phage1_phage2.unshared.genes.mash.distance","phage1_phage2.unshared.genes.mash.p.value","phage1_phage2.unshared.genes.kmer.count","phage1...unshared.phams",
#                                    "phage1...all.phams","phage1...shared.genes","phage1...unshared.genes","phage1...all.genes","phage1...average.length.of.all.genes","phage1.average.length.of.shared.genes","phage1.average.length.of.unshared.genes","phage1.total.length.of.all.genes",
#                                    "phage1.total.length.of.shared.genes","phage1.total.length.of.unshared.genes","phage1...all.genes.GC.content","phage1...shared.genes.GC.content","phage1...unshared.genes.GC.content","phage2...unshared.phams","phage2...all.phams","phage2...shared.genes",
#                                    "phage2...unshared.genes","phage2...all.genes","phage2...average.length.of.all.genes","phage2.average.length.of.shared.genes","phage2.average.length.of.unshared.genes","phage2.total.length.of.all.genes","phage2.total.length.of.shared.genes","phage2.total.length.of.unshared.genes",
#                                    "phage2...all.genes.GC.content","phage2...shared.genes.GC.content","phage2...unshared.genes.GC.content")

#Original column headers copied from Excel, with no truncation:
#"phage1_phage2","phage1","phage2","phage1_phage2","# shared phams","phage1_phage2 gene content dissimilarity (general index)","phage1_phage2 gene content dissimilarity (jaccard index)","phage1_phage2 all genes mash distance","phage1_phage2 all genes mash p-value","phage1_phage2 all genes kmer count","phage1_phage2 shared genes mash distance","phage1_phage2 shared genes mash p-value","phage1_phage2 shared genes kmer count","phage1_phage2 unshared genes mash distance","phage1_phage2 unshared genes mash p-value","phage1_phage2 unshared genes kmer count","phage1 # unshared phams","phage1 # all phams","phage1 # shared genes","phage1 # unshared genes","phage1 # all genes","phage1 - average length of all genes","phage1 average length of shared genes","phage1 average length of unshared genes","phage1 total length of all genes","phage1 total length of shared genes","phage1 total length of unshared genes","phage1 - all genes GC content","phage1 - shared genes GC content","ref - unshared genes GC content","phage2 # unshared phams","phage2 # all phams","phage2 # shared genes","phage2 # unshared genes","phage2 # all genes","phage2 - average length of all genes","phage2 average length of shared genes","phage2 average length of unshared genes","phage2 total length of all genes","phage2 total length of shared genes","phage2 total length of unshared genes","phage2 - all genes GC content","phage2 - shared genes GC content","phage2 - unshared genes GC content"

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







#Now merge the gene-specific mash data to the original data used.
#Since the ref_query comparison identifiers should be identical, there should be equal rows in both tables
gene_specific_mash_table <- merge(gene_specific_mash_table,gene_specific_mash_data,by.x="mash_ref_query",by.y="gsm_mash_ref_query")






#There several types of misleading data in the gsm dataset, since some comparisons have no shared or unshared genes. The program can't compute some values, and so must output default data. This data needs filtered out.
#1. For failed mash computations, the script outputs mash distance = 0, pvalue = 1, and kmer count = 0/0. If I filter by pvalue, these are successfully removed.
#2. For failed gene ave size computations, the script outputs ave gene size = 0. These should be converted to "NA".
#3. For failed GC content computations, the script outputs GC content = 0. These should be converted to "NA".

  
#Choose which data should be retained. Those not passing the filter getting Mash distances re-assigned to 0.5.
gene_specific_mash_table$gsm_all_mash_filter <- ifelse(gene_specific_mash_table$gsm_all_mash_pvalue < 1e-10,TRUE,FALSE)
gene_specific_mash_table$gsm_shared_mash_filter <- ifelse(gene_specific_mash_table$gsm_shared_mash_pvalue < 1e-10,TRUE,FALSE)
gene_specific_mash_table$gsm_unshared_mash_filter <- ifelse(gene_specific_mash_table$gsm_unshared_mash_pvalue < 1e-10,TRUE,FALSE)
gene_specific_mash_table$gsm_shared_unshared_mash_filter <- ifelse(gene_specific_mash_table$gsm_shared_unshared_mash_pvalue < 1e-10,TRUE,FALSE)




#Compute modified values
#Remove gsm mash distances with high pvalues, to determine the range of distance values, then created modified distance data
#gene_specific_mash_table_temp1 <- subset(gene_specific_mash_table,gene_specific_mash_table$gsm_all_mash_filter == TRUE)
#gene_specific_mash_table_temp2 <- subset(gene_specific_mash_table,gene_specific_mash_table$gsm_shared_mash_filter == TRUE)
#gene_specific_mash_table_temp3 <- subset(gene_specific_mash_table,gene_specific_mash_table$gsm_unshared_mash_filter == TRUE)

#At this point, the max mash distance of all filtered gsm_all_mash_distance comparisons = 0.5095.
#At this point, the max mash distance of all filtered gsm_shared_mash_distance comparisons = 0.52160
#At this point, the max mash distance of all filtered gsm_unshared_mash_distance comparisons = 0.5365

#Create modified data for mash
gene_specific_mash_table$gsm_all_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_all_mash_filter == TRUE,gene_specific_mash_table$gsm_all_mash_distance,0.6)
gene_specific_mash_table$gsm_shared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_shared_mash_filter == TRUE,gene_specific_mash_table$gsm_shared_mash_distance,0.6)
gene_specific_mash_table$gsm_unshared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_unshared_mash_filter == TRUE,gene_specific_mash_table$gsm_unshared_mash_distance,0.6)
gene_specific_mash_table$gsm_shared_unshared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_shared_unshared_mash_filter == TRUE,gene_specific_mash_table$gsm_shared_unshared_mash_distance,0.6)

#Create modified data for average gene size
gene_specific_mash_table$gsm_ref_all_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_ref_num_all_genes > 0,gene_specific_mash_table$gsm_ref_all_ave_size,NA)
gene_specific_mash_table$gsm_ref_shared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_ref_num_shared_genes > 0,gene_specific_mash_table$gsm_ref_shared_ave_size,NA)
gene_specific_mash_table$gsm_ref_unshared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_ref_num_unshared_genes > 0,gene_specific_mash_table$gsm_ref_unshared_ave_size,NA)
gene_specific_mash_table$gsm_query_all_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_query_num_all_genes > 0,gene_specific_mash_table$gsm_query_all_ave_size,NA)
gene_specific_mash_table$gsm_query_shared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_query_num_shared_genes > 0,gene_specific_mash_table$gsm_query_shared_ave_size,NA)
gene_specific_mash_table$gsm_query_unshared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_query_num_unshared_genes > 0,gene_specific_mash_table$gsm_query_unshared_ave_size,NA)




#Create modified data for GC content
gene_specific_mash_table$gsm_ref_all_modified_GC <- ifelse(gene_specific_mash_table$gsm_ref_all_total_size > 0,gene_specific_mash_table$gsm_ref_all_GC,NA)
gene_specific_mash_table$gsm_ref_shared_modified_GC <- ifelse(gene_specific_mash_table$gsm_ref_shared_total_size > 0,gene_specific_mash_table$gsm_ref_shared_GC,NA)
gene_specific_mash_table$gsm_ref_unshared_modified_GC <- ifelse(gene_specific_mash_table$gsm_ref_unshared_total_size > 0,gene_specific_mash_table$gsm_ref_unshared_GC,NA)
gene_specific_mash_table$gsm_query_all_modified_GC <- ifelse(gene_specific_mash_table$gsm_query_all_total_size > 0,gene_specific_mash_table$gsm_query_all_GC,NA)
gene_specific_mash_table$gsm_query_shared_modified_GC <- ifelse(gene_specific_mash_table$gsm_query_shared_total_size > 0,gene_specific_mash_table$gsm_query_shared_GC,NA)
gene_specific_mash_table$gsm_query_unshared_modified_GC <- ifelse(gene_specific_mash_table$gsm_query_unshared_total_size > 0,gene_specific_mash_table$gsm_query_unshared_GC,NA)
gene_specific_mash_table$gsm_shared_modified_GC <- ifelse(gene_specific_mash_table$gsm_shared_total_size > 0,gene_specific_mash_table$gsm_shared_GC,NA)
gene_specific_mash_table$gsm_unshared_modified_GC <- ifelse(gene_specific_mash_table$gsm_unshared_total_size > 0,gene_specific_mash_table$gsm_unshared_GC,NA)




#Plots to check the modifications
# plot(gene_specific_mash_table$gsm_all_mash_distance,gene_specific_mash_table$gsm_all_modified_mash_distance,xlim = c(0,0.6),ylim = c(0,0.6))
# plot(gene_specific_mash_table$gsm_shared_mash_distance,gene_specific_mash_table$gsm_shared_modified_mash_distance,xlim = c(0,0.6),ylim = c(0,0.6))
# plot(gene_specific_mash_table$gsm_unshared_mash_distance,gene_specific_mash_table$gsm_unshared_modified_mash_distance,xlim = c(0,0.6),ylim = c(0,0.6))
# 
# plot(gene_specific_mash_table$gsm_ref_all_ave_size,gene_specific_mash_table$gsm_ref_all_modified_ave_size)
# plot(gene_specific_mash_table$gsm_ref_shared_ave_size,gene_specific_mash_table$gsm_ref_shared_modified_ave_size)
# plot(gene_specific_mash_table$gsm_ref_unshared_ave_size,gene_specific_mash_table$gsm_ref_unshared_modified_ave_size)
# 
# plot(gene_specific_mash_table$gsm_query_all_ave_size,gene_specific_mash_table$gsm_query_all_modified_ave_size)
# plot(gene_specific_mash_table$gsm_query_shared_ave_size,gene_specific_mash_table$gsm_query_shared_modified_ave_size)
# plot(gene_specific_mash_table$gsm_query_unshared_ave_size,gene_specific_mash_table$gsm_query_unshared_modified_ave_size)
# 
# plot(gene_specific_mash_table$gsm_ref_all_GC,gene_specific_mash_table$gsm_ref_all_modified_GC)
# plot(gene_specific_mash_table$gsm_ref_shared_GC,gene_specific_mash_table$gsm_ref_shared_modified_GC)
# plot(gene_specific_mash_table$gsm_ref_unshared_GC,gene_specific_mash_table$gsm_ref_unshared_modified_GC)
# 
# plot(gene_specific_mash_table$gsm_query_all_GC,gene_specific_mash_table$gsm_query_all_modified_GC)
# plot(gene_specific_mash_table$gsm_query_shared_GC,gene_specific_mash_table$gsm_query_shared_modified_GC)
# plot(gene_specific_mash_table$gsm_query_unshared_GC,gene_specific_mash_table$gsm_query_unshared_modified_GC)

  


#QC and Analysis values to compute, now that the gene-specific mash data has been added
#Misc
gene_specific_mash_table$gsm_num_unshared_phams <- gene_specific_mash_table$gsm_ref_num_unshared_phams + gene_specific_mash_table$gsm_query_num_unshared_phams

gene_specific_mash_table$gsm_num_all_genes <- gene_specific_mash_table$gsm_ref_num_all_genes + gene_specific_mash_table$gsm_query_num_all_genes
gene_specific_mash_table$gsm_num_shared_genes <- gene_specific_mash_table$gsm_ref_num_shared_genes + gene_specific_mash_table$gsm_query_num_shared_genes
gene_specific_mash_table$gsm_num_unshared_genes <- gene_specific_mash_table$gsm_ref_num_unshared_genes + gene_specific_mash_table$gsm_query_num_unshared_genes

#Difference between number of genes and number of phams
gene_specific_mash_table$gsm_ref_all_pham_gene_disparity <- gene_specific_mash_table$gsm_ref_num_all_genes - gene_specific_mash_table$gsm_ref_num_all_phams
gene_specific_mash_table$gsm_ref_shared_pham_gene_disparity <- gene_specific_mash_table$gsm_ref_num_shared_genes - gene_specific_mash_table$gsm_num_shared_phams
gene_specific_mash_table$gsm_ref_unshared_pham_gene_disparity <- gene_specific_mash_table$gsm_ref_num_unshared_genes - gene_specific_mash_table$gsm_ref_num_unshared_phams

gene_specific_mash_table$gsm_query_all_pham_gene_disparity <- gene_specific_mash_table$gsm_query_num_all_genes - gene_specific_mash_table$gsm_query_num_all_phams
gene_specific_mash_table$gsm_query_shared_pham_gene_disparity <- gene_specific_mash_table$gsm_query_num_shared_genes - gene_specific_mash_table$gsm_num_shared_phams
gene_specific_mash_table$gsm_query_unshared_pham_gene_disparity <- gene_specific_mash_table$gsm_query_num_unshared_genes - gene_specific_mash_table$gsm_query_num_unshared_phams


#Estimate of the proportion of coding sequence per genome
#Note: the gene-specific sequence length used does not take into account overlapping CDS features, so the sequences could have duplicate regions in the genome
#Note: the 'all coding potential' data is based on the real genome size, but the 'shared/unshared coding potential' data is based only on the gene-specific mash sizes

gene_specific_mash_table$gsm_all_total_size <- gene_specific_mash_table$gsm_ref_all_total_size + gene_specific_mash_table$gsm_query_all_total_size


gene_specific_mash_table$gsm_ref_all_coding_potential <- gene_specific_mash_table$gsm_ref_all_total_size/gene_specific_mash_table$ref_size
gene_specific_mash_table$gsm_query_all_coding_potential <- gene_specific_mash_table$gsm_query_all_total_size/gene_specific_mash_table$query_size

gene_specific_mash_table$gsm_all_coding_potential <- (gene_specific_mash_table$gsm_ref_all_coding_potential + gene_specific_mash_table$gsm_query_all_coding_potential)/2


gene_specific_mash_table$gsm_ref_shared_coding_potential <- gene_specific_mash_table$gsm_ref_shared_total_size/gene_specific_mash_table$gsm_ref_all_total_size
gene_specific_mash_table$gsm_query_shared_coding_potential <- gene_specific_mash_table$gsm_query_shared_total_size/gene_specific_mash_table$gsm_query_all_total_size

gene_specific_mash_table$gsm_ref_unshared_coding_potential <- gene_specific_mash_table$gsm_ref_unshared_total_size/gene_specific_mash_table$gsm_ref_all_total_size
gene_specific_mash_table$gsm_query_unshared_coding_potential <- gene_specific_mash_table$gsm_query_unshared_total_size/gene_specific_mash_table$gsm_query_all_total_size

gene_specific_mash_table$gsm_shared_coding_potential <- gene_specific_mash_table$gsm_shared_total_size / gene_specific_mash_table$gsm_all_total_size
gene_specific_mash_table$gsm_unshared_coding_potential <- gene_specific_mash_table$gsm_unshared_total_size / gene_specific_mash_table$gsm_all_total_size





#Differences in GC content
gene_specific_mash_table$gsm_modified_gc_diff_all_all <- abs(gene_specific_mash_table$gsm_ref_all_modified_GC - gene_specific_mash_table$gsm_query_all_modified_GC)
gene_specific_mash_table$gsm_modified_gc_diff_shared_shared <- abs(gene_specific_mash_table$gsm_ref_shared_modified_GC - gene_specific_mash_table$gsm_query_shared_modified_GC)
gene_specific_mash_table$gsm_modified_gc_diff_unshared_unshared <- abs(gene_specific_mash_table$gsm_ref_unshared_modified_GC - gene_specific_mash_table$gsm_query_unshared_modified_GC)
gene_specific_mash_table$gsm_modified_gc_diff_shared_unshared <- abs(gene_specific_mash_table$gsm_shared_modified_GC - gene_specific_mash_table$gsm_unshared_modified_GC)


#Old code that is probably unnecessary
# gene_specific_mash_table$gc_content_diff_all_all <- abs(gene_specific_mash_table$gsm_ref_all_GC - gene_specific_mash_table$gsm_query_all_GC)
# gene_specific_mash_table$ref_gc_content_diff_shared_all <- abs(gene_specific_mash_table$gsm_ref_all_GC - gene_specific_mash_table$gsm_ref_shared_GC)
# gene_specific_mash_table$ref_gc_content_diff_unshared_all <- abs(gene_specific_mash_table$gsm_ref_all_GC - gene_specific_mash_table$gsm_ref_unshared_GC)
# gene_specific_mash_table$query_gc_content_diff_shared_all <- abs(gene_specific_mash_table$gsm_query_all_GC - gene_specific_mash_table$gsm_query_shared_GC)
# gene_specific_mash_table$query_gc_content_diff_unshared_all <- abs(gene_specific_mash_table$gsm_query_all_GC - gene_specific_mash_table$gsm_query_unshared_GC)
# 
# gene_specific_mash_table$ref_gc_content_diff_shared_unshared <- abs(gene_specific_mash_table$gsm_ref_unshared_GC - gene_specific_mash_table$gsm_ref_shared_GC)
# gene_specific_mash_table$query_gc_content_diff_shared_unshared <- abs(gene_specific_mash_table$gsm_query_unshared_GC - gene_specific_mash_table$gsm_query_shared_GC)




#Difference in shared and unshared mash distances
gene_specific_mash_table$gsm_mash_distance_diff_unshared_shared <- gene_specific_mash_table$gsm_unshared_modified_mash_distance - gene_specific_mash_table$gsm_shared_modified_mash_distance


#Compute gene content dissimilarity by shared genes instead of shared phams
gene_specific_mash_table$gsm_gene_dissimilarity <- 1 - ((gene_specific_mash_table$gsm_ref_num_shared_genes/gene_specific_mash_table$gsm_ref_num_all_genes) + (gene_specific_mash_table$gsm_query_num_shared_genes/gene_specific_mash_table$gsm_query_num_all_genes))/2




#Average gene sizes
#This data does not use the modified gene size data, because if there are no shared/unshared genes, the default ave size is 0, and this is fine to use to compute averages
gene_specific_mash_table$gsm_ave_size_all_all <- (gene_specific_mash_table$gsm_ref_all_ave_size + gene_specific_mash_table$gsm_query_all_ave_size)/2
gene_specific_mash_table$gsm_ave_size_shared_shared <- (gene_specific_mash_table$gsm_ref_shared_ave_size + gene_specific_mash_table$gsm_query_shared_ave_size)/2
gene_specific_mash_table$gsm_ave_size_unshared_unshared <- (gene_specific_mash_table$gsm_ref_unshared_ave_size + gene_specific_mash_table$gsm_query_unshared_ave_size)/2


#Average GC
#This data uses the raw and modified GC content. The script default if there are no shared/unshared genes, is to set GC content to 0, but this could be confused with a biologically relevant result.
#Get around this with a two-step process:
#First, check if EITHER modified GC content = NA. If it is, simply use only the other genome's data. If neither are NA, then average them.
#Second, check if BOTH modified GC contents = NA. If they do, set to NA.
gene_specific_mash_table$gsm_ave_gc_all_all <- ifelse(is.na(gene_specific_mash_table$gsm_ref_all_modified_GC) == TRUE | is.na(gene_specific_mash_table$gsm_query_all_modified_GC) == TRUE,
                                                                 (gene_specific_mash_table$gsm_ref_all_GC + gene_specific_mash_table$gsm_query_all_GC),
                                                                 (gene_specific_mash_table$gsm_ref_all_GC + gene_specific_mash_table$gsm_query_all_GC)/2)
gene_specific_mash_table$gsm_ave_gc_all_all <- ifelse(is.na(gene_specific_mash_table$gsm_ref_all_modified_GC) == TRUE & is.na(gene_specific_mash_table$gsm_query_all_modified_GC) == TRUE,
                                                           NA,gene_specific_mash_table$gsm_ave_gc_all_all)

gene_specific_mash_table$gsm_ave_gc_shared_shared <- ifelse(is.na(gene_specific_mash_table$gsm_ref_shared_modified_GC) == TRUE | is.na(gene_specific_mash_table$gsm_query_shared_modified_GC) == TRUE,
                                                                 (gene_specific_mash_table$gsm_ref_shared_GC + gene_specific_mash_table$gsm_query_shared_GC),
                                                                 (gene_specific_mash_table$gsm_ref_shared_GC + gene_specific_mash_table$gsm_query_shared_GC)/2)
gene_specific_mash_table$gsm_ave_gc_shared_shared <- ifelse(is.na(gene_specific_mash_table$gsm_ref_shared_modified_GC) == TRUE & is.na(gene_specific_mash_table$gsm_query_shared_modified_GC) == TRUE,
                                                           NA,gene_specific_mash_table$gsm_ave_gc_shared_shared)

gene_specific_mash_table$gsm_ave_gc_unshared_unshared <- ifelse(is.na(gene_specific_mash_table$gsm_ref_unshared_modified_GC) == TRUE | is.na(gene_specific_mash_table$gsm_query_unshared_modified_GC) == TRUE,
                                                                 (gene_specific_mash_table$gsm_ref_unshared_GC + gene_specific_mash_table$gsm_query_unshared_GC),
                                                                 (gene_specific_mash_table$gsm_ref_unshared_GC + gene_specific_mash_table$gsm_query_unshared_GC)/2)
gene_specific_mash_table$gsm_ave_gc_unshared_unshared <- ifelse(is.na(gene_specific_mash_table$gsm_ref_unshared_modified_GC) == TRUE & is.na(gene_specific_mash_table$gsm_query_unshared_modified_GC) == TRUE,
                                                           NA,gene_specific_mash_table$gsm_ave_gc_unshared_unshared)

#Subset to verify GC averaging worked okay
#temp <- subset(gene_specific_mash_table,is.na(gene_specific_mash_table$gsm_ref_unshared_modified_GC) == TRUE,select=c("gsm_ref_unshared_modified_GC","gsm_query_unshared_modified_GC","gsm_ref_unshared_GC","gsm_query_unshared_GC","gsm_ave_gc_unshared_unshared"))
#temp <- subset(gene_specific_mash_table,is.na(gene_specific_mash_table$gsm_query_unshared_modified_GC) == TRUE,select=c("gsm_ref_unshared_modified_GC","gsm_query_unshared_modified_GC","gsm_ref_unshared_GC","gsm_query_unshared_GC","gsm_ave_gc_unshared_unshared"))
#temp <- subset(gene_specific_mash_table,is.na(gene_specific_mash_table$gsm_ref_unshared_modified_GC) == FALSE & is.na(gene_specific_mash_table$gsm_query_unshared_modified_GC) == FALSE,select=c("gsm_ref_unshared_modified_GC","gsm_query_unshared_modified_GC","gsm_ref_unshared_GC","gsm_query_unshared_GC","gsm_ave_gc_unshared_unshared"))









#Convert all Unspecified fields to NA missing value
gene_specific_mash_table[gene_specific_mash_table == "Unspecified"] <- NA

















#Analysis before splitting into gene flux mode

#Compare whole genome mash distance to gene-specific distance
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table$modified_mash_distance,gene_specific_mash_table$gsm_all_modified_mash_distance,xlim=c(0,0.6),ylim=c(0,0.6))
abline(0,1,lty=2,lwd=3,col="grey")

#Compare gene content dissimilarity based on pham counts or gene counts
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table$pham_pham_dissimilarity,gene_specific_mash_table$gsm_gene_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")


#Compare shared and unshared distances
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table$gsm_shared_modified_mash_distance,gene_specific_mash_table$gsm_unshared_modified_mash_distance,xlim=c(0,0.6),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")

















#Make backup table before subsetting
gene_specific_mash_table_complete <- gene_specific_mash_table

gene_specific_mash_table <- gene_specific_mash_table_complete

#Now decide to retain only SEA-PHAGES annotated genomes or not
gene_specific_mash_table <- subset(gene_specific_mash_table,gene_specific_mash_table$phage_cluster_source_compare == 'actino' & gene_specific_mash_table$phage_cluster_compare != 'different')





#Analysis after splitting into gene flux modes

#Split into HGCF and LGCF datasets

gene_specific_mash_table_temperate_hgcf <- subset(gene_specific_mash_table,gene_specific_mash_table$gene_flux_category == 'high' & gene_specific_mash_table$phage_temperate_compare == 'yes')
gene_specific_mash_table_temperate_lgcf <- subset(gene_specific_mash_table,gene_specific_mash_table$gene_flux_category == 'low' & gene_specific_mash_table$phage_temperate_compare == 'yes')
gene_specific_mash_table_lytic <- subset(gene_specific_mash_table,gene_specific_mash_table$phage_temperate_compare == 'no')







#Compare average gene sizes
temp_temperate_hgcf_shared <- subset(gene_specific_mash_table_temperate_hgcf,select=c("gsm_ave_size_shared_shared"))
temp_temperate_lgcf_shared <- subset(gene_specific_mash_table_temperate_lgcf,select=c("gsm_ave_size_shared_shared"))
temp_lytic_shared <- subset(gene_specific_mash_table_lytic,select=c("gsm_ave_size_shared_shared"))

temp_temperate_hgcf_unshared <- subset(gene_specific_mash_table_temperate_hgcf,select=c("gsm_ave_size_unshared_unshared"))
temp_temperate_lgcf_unshared <- subset(gene_specific_mash_table_temperate_lgcf,select=c("gsm_ave_size_unshared_unshared"))
temp_lytic_unshared <- subset(gene_specific_mash_table_lytic,select=c("gsm_ave_size_unshared_unshared"))

names(temp_temperate_hgcf_shared) <- c("ave_gene_size")
names(temp_temperate_lgcf_shared) <- c("ave_gene_size")
names(temp_lytic_shared) <- c("ave_gene_size")
names(temp_temperate_hgcf_unshared) <- c("ave_gene_size")
names(temp_temperate_lgcf_unshared) <- c("ave_gene_size")
names(temp_lytic_unshared) <- c("ave_gene_size")

temp_temperate_hgcf_shared$category <- "temperate_HGCF"
temp_temperate_lgcf_shared$category <- "temperate_LGCF"
temp_lytic_shared$category <- "lytic"
temp_temperate_hgcf_unshared$category <- "temperate_HGCF"
temp_temperate_lgcf_unshared$category <- "temperate_LGCF"
temp_lytic_unshared$category <- "lytic"

temp_temperate_hgcf_shared$gene_type <- "shared"
temp_temperate_lgcf_shared$gene_type <- "shared"
temp_lytic_shared$gene_type <- "shared"
temp_temperate_hgcf_unshared$gene_type <- "unshared"
temp_temperate_lgcf_unshared$gene_type <- "unshared"
temp_lytic_unshared$gene_type <- "unshared"

temp_plotting_table <- rbind(temp_temperate_hgcf_shared,temp_temperate_lgcf_shared,temp_lytic_shared,temp_temperate_hgcf_unshared,temp_temperate_lgcf_unshared,temp_lytic_unshared)

par(mar=c(10,4,4,4))
boxplot(temp_plotting_table$ave_gene_size ~ temp_plotting_table$category*temp_plotting_table$gene_type,las=2,col=c("dark red","dark blue","dark green","pink","light blue","light green"))













#How distant are shared and unshared sequences?



#Compare shared to unshared
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_shared_modified_mash_distance,gene_specific_mash_table_lytic$gsm_unshared_modified_mash_distance,xlim=c(0,0.6),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$gsm_shared_modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.6),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$gsm_shared_modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.6),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,1,lty=2,lwd=3,col="grey")
# 
# 
#Compare all to shared and unshared by mash distance
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_shared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
abline(0,1,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_shared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_shared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")










par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$size_diff,xlim=c(0,0.5),ylim=c(0,1e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$size_diff,xlim=c(0,0.5),ylim=c(0,1e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")








#Compare un-averaged genome size disparity
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$size_diff,xlim=c(0,0.5),ylim=c(0,10000),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$size_diff,xlim=c(0,0.5),ylim=c(0,10000),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$size_diff,xlim=c(0,0.5),ylim=c(0,10000),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")









#Compare all to shared and unshared by gene content dissimilarity
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_shared_modified_mash_distance,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$gsm_unshared_modified_mash_distance,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
abline(0,1,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$gsm_shared_modified_mash_distance,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$gsm_unshared_modified_mash_distance,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$gsm_shared_modified_mash_distance,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$gsm_unshared_modified_mash_distance,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,0.6),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")




# 
# 
# 
# 
# 
# par(new=FALSE)
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_mash_distance_diff_unshared_shared,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_mash_distance_diff_unshared_shared,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# 
# 
# 
# 
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_shared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_shared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark grey")
# abline(0,1,lty=2,lwd=3,col="grey")
# 
# 
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0,2.5e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark grey")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0,2.5e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
# 
# 
# 
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_shared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_shared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_total_size,xlim=c(0,0.5),ylim=c(0,2e5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# abline(0,1,lty=2,lwd=3,col="grey")




#Coding potential

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_shared_coding_potential,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$gsm_shared_coding_potential,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$gsm_shared_coding_potential,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
abline(0,1,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_unshared_coding_potential,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$gsm_unshared_coding_potential,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$gsm_unshared_coding_potential,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")



par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_unshared_coding_potential,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$gsm_unshared_coding_potential,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$gsm_unshared_coding_potential,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")







#
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")







par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")




#Compare unshared and shared ave gene sizes
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_ave_size_shared_shared,xlim=c(0,0.5),ylim=c(0,1500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")


par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_ave_size_shared_shared,xlim=c(0,0.5),ylim=c(0,1500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")



par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_ave_size_shared_shared,xlim=c(0,0.5),ylim=c(0,1500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1500),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")





#Compare proportion of shared and unshared coding proportion
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")















#Analysis after splitting into gene flux modes and after computing sliding windows
#sliding window average
library(caTools)

gene_specific_mash_table_temperate_hgcf_mmdsort <- gene_specific_mash_table_temperate_hgcf[order(gene_specific_mash_table_temperate_hgcf$modified_mash_distance),]
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_modified_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_ave_size_shared_shared_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_ave_size_shared_shared,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_ave_size_unshared_unshared_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_ave_size_unshared_unshared,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_coding_potential,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_coding_potential,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_shared_genes_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_shared_genes,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_unshared_genes_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_unshared_genes,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_shared_phams_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_shared_phams,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_unshared_phams_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_unshared_phams,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_all_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_all_coding_potential,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_mash_distance_diff_unshared_shared_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_mash_distance_diff_unshared_shared,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_unshared_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_unshared_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_all_total_size_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_all_total_size,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_total_size,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_total_size,101)
gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff,101)

#gene_specific_mash_table_temperate_hgcf_mmdsort$ <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$,101)













gene_specific_mash_table_temperate_lgcf_mmdsort <- gene_specific_mash_table_temperate_lgcf[order(gene_specific_mash_table_temperate_lgcf$modified_mash_distance),]
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_modified_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_ave_size_shared_shared_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_ave_size_shared_shared,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_ave_size_unshared_unshared_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_ave_size_unshared_unshared,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_coding_potential,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_coding_potential,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_shared_genes_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_shared_genes,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_unshared_genes_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_unshared_genes,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_shared_phams_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_shared_phams,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_unshared_phams_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_unshared_phams,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_all_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_all_coding_potential,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_mash_distance_diff_unshared_shared_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_mash_distance_diff_unshared_shared,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_unshared_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_unshared_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_all_total_size_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_all_total_size,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_total_size,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_total_size,101)
gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff,101)



#gene_specific_mash_table_temperate_lgcf_mmdsort$ <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$,101)










gene_specific_mash_table_lytic_mmdsort <- gene_specific_mash_table_lytic[order(gene_specific_mash_table_lytic$modified_mash_distance),]
gene_specific_mash_table_lytic_mmdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_lytic_mmdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_unshared_modified_mash_distance,101)
gene_specific_mash_table_lytic_mmdsort$gsm_ave_size_shared_shared_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_ave_size_shared_shared,101)
gene_specific_mash_table_lytic_mmdsort$gsm_ave_size_unshared_unshared_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_ave_size_unshared_unshared,101)
gene_specific_mash_table_lytic_mmdsort$gsm_shared_coding_potential_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_shared_coding_potential,101)
gene_specific_mash_table_lytic_mmdsort$gsm_unshared_coding_potential_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_unshared_coding_potential,101)
gene_specific_mash_table_lytic_mmdsort$gsm_num_shared_genes_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_num_shared_genes,101)
gene_specific_mash_table_lytic_mmdsort$gsm_num_unshared_genes_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_num_unshared_genes,101)
gene_specific_mash_table_lytic_mmdsort$gsm_num_shared_phams_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_num_shared_phams,101)
gene_specific_mash_table_lytic_mmdsort$gsm_num_unshared_phams_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_num_unshared_phams,101)
gene_specific_mash_table_lytic_mmdsort$gsm_all_coding_potential_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_all_coding_potential,101)
gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent,101)
gene_specific_mash_table_lytic_mmdsort$gsm_mash_distance_diff_unshared_shared_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_mash_distance_diff_unshared_shared,101)
gene_specific_mash_table_lytic_mmdsort$gsm_shared_unshared_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_shared_unshared_mash_distance,101)
gene_specific_mash_table_lytic_mmdsort$gsm_all_total_size_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_all_total_size,101)
gene_specific_mash_table_lytic_mmdsort$gsm_shared_total_size_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_shared_total_size,101)
gene_specific_mash_table_lytic_mmdsort$gsm_unshared_total_size_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$gsm_unshared_total_size,101)
gene_specific_mash_table_lytic_mmdsort$size_diff_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$size_diff,101)


#gene_specific_mash_table_lytic_mmdsort$ <- runmean(gene_specific_mash_table_lytic_mmdsort$,101)








#Standard plot as a reference
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")





# #Compare all coding potential
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_all_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0.9,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_all_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0.9,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_all_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0.9,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

#Compare all genome size disparity
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

# #Compare difference between shared and unshared mash distances
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_mash_distance_diff_unshared_shared_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_mash_distance_diff_unshared_shared_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_mash_distance_diff_unshared_shared_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# 
# #Compare mash distance between shared and unshared genes
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_unshared_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_unshared_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_shared_unshared_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# 
# 
# #Compare all total size
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_all_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1.5e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_all_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1.5e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_all_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1.5e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# 
# 
# 
# 
# 
# 
# 
# 
#Compare shared and unshared distances
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_shared_modified_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_modified_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_modified_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_unshared_modified_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_modified_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_modified_mash_distance_runmean,xlim=c(0,0.5),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")

#Compare shared and unshared coding potential
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_shared_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")





# #Compare shared and unshared average gene sizes
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.5),ylim=c(0,1100),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.5),ylim=c(0,1100),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.5),ylim=c(0,1100),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.5),ylim=c(0,1100),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.5),ylim=c(0,1100),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.5),ylim=c(0,1100),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
# 
# 
# #Compare shared and unshared number of genes
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_num_shared_genes_runmean,xlim=c(0,0.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_shared_genes_runmean,xlim=c(0,0.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_shared_genes_runmean,xlim=c(0,0.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_num_unshared_genes_runmean,xlim=c(0,0.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_unshared_genes_runmean,xlim=c(0,0.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_unshared_genes_runmean,xlim=c(0,0.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
# 
# 
# #Compare shared and unshared number of phams
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_num_shared_phams_runmean,xlim=c(0,0.5),ylim=c(0,200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_shared_phams_runmean,xlim=c(0,0.5),ylim=c(0,200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_shared_phams_runmean,xlim=c(0,0.5),ylim=c(0,200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$gsm_num_unshared_phams_runmean,xlim=c(0,0.5),ylim=c(0,200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_num_unshared_phams_runmean,xlim=c(0,0.5),ylim=c(0,200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_num_unshared_phams_runmean,xlim=c(0,0.5),ylim=c(0,200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
# 
# 
# 
# 
# 
# #Compare shared and unshared total size
# par(mar=c(4,8,4,4))
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_shared_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_shared_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$gsm_unshared_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$gsm_unshared_total_size_runmean,xlim=c(0,0.5),ylim=c(0,1e5),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")

























#sliding window analyis based on gene content and NOT nucleotide distance
gene_specific_mash_table_temperate_hgcf_gcdsort <- gene_specific_mash_table_temperate_hgcf[order(gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity),]
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_ave_size_shared_shared_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_ave_size_shared_shared,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_ave_size_unshared_unshared_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_ave_size_unshared_unshared,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_coding_potential,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_coding_potential,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_shared_genes_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_shared_genes,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_unshared_genes_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_unshared_genes,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_shared_phams_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_shared_phams,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_unshared_phams_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_num_unshared_phams,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_all_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_all_coding_potential,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$size_diff_ave_percent,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_mash_distance_diff_unshared_shared_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_mash_distance_diff_unshared_shared,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_unshared_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_unshared_mash_distance,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_all_total_size_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_all_total_size,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_total_size,101)
gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_total_size,101)














gene_specific_mash_table_temperate_lgcf_gcdsort <- gene_specific_mash_table_temperate_lgcf[order(gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity),]
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_ave_size_shared_shared_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_ave_size_shared_shared,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_ave_size_unshared_unshared_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_ave_size_unshared_unshared,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_coding_potential,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_coding_potential,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_shared_genes_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_shared_genes,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_unshared_genes_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_unshared_genes,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_shared_phams_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_shared_phams,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_unshared_phams_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_num_unshared_phams,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_all_coding_potential_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_all_coding_potential,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$size_diff_ave_percent,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_mash_distance_diff_unshared_shared_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_mash_distance_diff_unshared_shared,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_unshared_mash_distance_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_unshared_mash_distance,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_all_total_size_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_all_total_size,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_total_size,101)
gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_total_size_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_total_size,101)












gene_specific_mash_table_lytic_gcdsort <- gene_specific_mash_table_lytic[order(gene_specific_mash_table_lytic$pham_pham_dissimilarity),]
gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance,101)
gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance,101)
gene_specific_mash_table_lytic_gcdsort$gsm_ave_size_shared_shared_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_ave_size_shared_shared,101)
gene_specific_mash_table_lytic_gcdsort$gsm_ave_size_unshared_unshared_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_ave_size_unshared_unshared,101)
gene_specific_mash_table_lytic_gcdsort$gsm_shared_coding_potential_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_shared_coding_potential,101)
gene_specific_mash_table_lytic_gcdsort$gsm_unshared_coding_potential_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_unshared_coding_potential,101)
gene_specific_mash_table_lytic_gcdsort$gsm_num_shared_genes_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_num_shared_genes,101)
gene_specific_mash_table_lytic_gcdsort$gsm_num_unshared_genes_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_num_unshared_genes,101)
gene_specific_mash_table_lytic_gcdsort$gsm_num_shared_phams_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_num_shared_phams,101)
gene_specific_mash_table_lytic_gcdsort$gsm_num_unshared_phams_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_num_unshared_phams,101)
gene_specific_mash_table_lytic_gcdsort$gsm_all_coding_potential_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_all_coding_potential,101)
gene_specific_mash_table_lytic_gcdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$size_diff_ave_percent,101)
gene_specific_mash_table_lytic_gcdsort$gsm_mash_distance_diff_unshared_shared_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_mash_distance_diff_unshared_shared,101)
gene_specific_mash_table_lytic_gcdsort$gsm_shared_unshared_mash_distance_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_shared_unshared_mash_distance,101)
gene_specific_mash_table_lytic_gcdsort$gsm_all_total_size_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_all_total_size,101)
gene_specific_mash_table_lytic_gcdsort$gsm_shared_total_size_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_shared_total_size,101)
gene_specific_mash_table_lytic_gcdsort$gsm_unshared_total_size_runmean <- runmean(gene_specific_mash_table_lytic_gcdsort$gsm_unshared_total_size,101)







#Compare shared and unshared distances
#par(mar=c(4,8,4,4))
#plot(gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
#par(new=TRUE)
#plot(gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
#par(new=TRUE)
#plot(gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
#par(new=TRUE)
#plot(gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
#par(new=TRUE)
# plot(gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
# par(new=TRUE)
# plot(gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")



par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_lytic_gcdsort$gsm_shared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_lytic_gcdsort$gsm_unshared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")

par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_modified_mash_distance,xlim=c(0,1),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")









#Standard plot for reference
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_lgcf_gcdsort$modified_mash_distance,xlim=c(0,1),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_temperate_hgcf_gcdsort$modified_mash_distance,xlim=c(0,1),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,gene_specific_mash_table_lytic_gcdsort$modified_mash_distance,xlim=c(0,1),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Compare shared and unshared distances
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



#Compare shared and unshared coding potential
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic_gcdsort$gsm_unshared_coding_potential_runmean,gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_unshared_coding_potential_runmean,gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_unshared_coding_potential_runmean,gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_gcdsort$gsm_shared_coding_potential_runmean,gene_specific_mash_table_lytic_gcdsort$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf_gcdsort$gsm_shared_coding_potential_runmean,gene_specific_mash_table_temperate_lgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_gcdsort$gsm_shared_coding_potential_runmean,gene_specific_mash_table_temperate_hgcf_gcdsort$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")























###Export HGCF and LGCF data for analysis in Excel and Phamerator

hgcf_lgcf_data <- mash_table2
hgcf_lgcf_data <- subset(hgcf_lgcf_data,hgcf_lgcf_data$phage_cluster_source_compare == 'actino' & hgcf_lgcf_data$modified_mash_distance < 0.42 & hgcf_lgcf_data$pham_pham_dissimilarity < 0.89)

temperate_hgcf <- subset(hgcf_lgcf_data,hgcf_lgcf_data$phage_temperate_compare == 'yes' & hgcf_lgcf_data$gene_flux_category == 'high',select=c('mash_reference','mash_query','modified_mash_distance','pham_pham_dissimilarity','ref_phage_cluster','query_phage_cluster'))
temperate_lgcf <- subset(hgcf_lgcf_data,hgcf_lgcf_data$phage_temperate_compare == 'yes' & hgcf_lgcf_data$gene_flux_category == 'low',select=c('mash_reference','mash_query','modified_mash_distance','pham_pham_dissimilarity','ref_phage_cluster','query_phage_cluster'))
lytic_hgcf_lgcf <- subset(hgcf_lgcf_data,hgcf_lgcf_data$phage_temperate_compare == 'no',select=c('mash_reference','mash_query','modified_mash_distance','pham_pham_dissimilarity','ref_phage_cluster','query_phage_cluster'))


write.table(temperate_hgcf,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170106_actino_temperate_hgcf_output.csv",sep=",")
write.table(temperate_lgcf,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170106_actino_temperate_lgcf_output.csv",sep=",")
write.table(lytic_hgcf_lgcf,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170106_actino_lytic_hgcf_lgcf_output.csv",sep=",")


par(mar=c(4,8,4,4))
plot(hgcf_lgcf_data$modified_mash_distance,hgcf_lgcf_data$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_hgcf$modified_mash_distance,temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_lgcf$modified_mash_distance,temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lytic_hgcf_lgcf$modified_mash_distance,lytic_hgcf_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")











###Misc old code that may be useful to have quick access to

#par(mar=c(4,8,4,4))
#plot(density(mash_table2$modified_mash_distance),xlim=c(0,0.5),ylim=c(0,2))
#polygon(density(mash_table2$modified_mash_distance),col="black",border="black")




###Pseudomonas and streptococcus phages

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
bacteria_dsDNA <- subset(type_dsDNA,type_dsDNA$host_superkingdom_compare == 'Bacteria')

pseudomonas <- subset(bacteria_dsDNA,bacteria_dsDNA$host_genus_compare == "Pseudomonas")
streptococcus <- subset(bacteria_dsDNA,bacteria_dsDNA$host_genus_compare == "Streptococcus")

par(mar=c(4,8,4,4))
plot(pseudomonas$modified_mash_distance,pseudomonas$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(streptococcus$modified_mash_distance,streptococcus$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



###Output mash and gene content data for publication
#Merge the ani79 data as before, but retain all rows
#Merge the gene-specific mash data as before, but retain all rows
#Change the column names to be more readable
output_table <- mash_table2

ani79_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20161006_set2_79_ani_data.csv",sep=",",header=TRUE)
names(ani79_data) <- c("ani79_ref_query","ani79_ref_phage_identifier","ani79_query_phage_identifier","ani79_ani_distance")
output_table <- merge(output_table,ani79_data,by.x="mash_ref_query",by.y="ani79_ref_query",all.x=TRUE)

gene_specific_mash_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170109_gene_specific_mash_analysis.csv",sep=",",header=TRUE)
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

output_table <- merge(output_table,gene_specific_mash_data,by.x="mash_ref_query",by.y="gsm_mash_ref_query",all.x=TRUE)
output_table_reduced <- subset(output_table,select=c("mash_reference","mash_query","mash_distance","mash_pvalue","pham_pham_dissimilarity",
                                                     "gsm_shared_mash_distance","gsm_shared_mash_pvalue",
                                                     "gsm_unshared_mash_distance","gsm_unshared_mash_pvalue"))

#If data truncation is needed to reduce filesize:
# output_table_reduced$mash_distance <- round(output_table_reduced$mash_distance, digits = 2)
# output_table_reduced$pham_pham_dissimilarity <- round(output_table_reduced$pham_pham_dissimilarity, digits = 2)
# output_table_reduced$ani79_ani_distance <- round(output_table_reduced$ani79_ani_distance, digits = 2)
# output_table_reduced$gsm_shared_mash_distance <- round(output_table_reduced$gsm_shared_mash_distance, digits = 2)
# output_table_reduced$gsm_unshared_mash_distance <- round(output_table_reduced$gsm_unshared_mash_distance, digits = 2)

names(output_table_reduced) = c("Virus 1","Virus 2",
                                "Whole genome mash raw distance","Whole genome mash pvalue",
                                "Gene content dissimilarity",
                                "Shared gene mash raw distance","Shared gene mash pvalue",
                                "Unshared gene mash raw distance","Unshared gene mash pvalue")
write.table(output_table_reduced,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/phage_evolution_data.csv",sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)








###Analysis for committee meeting

mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
bacteria_dsDNA <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')

#Cluster-boundary dataset
bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.42 & bacteria_dsDNA$pham_pham_dissimilarity < 0.89)

#Subcluster-boundary dataset
bacteria_dsDNA_nuc020_gene062 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.20 & bacteria_dsDNA$pham_pham_dissimilarity < 0.62)


cluster_archernm <- subset(bacteria_dsDNA_nuc042_gene089,bacteria_dsDNA_nuc042_gene089$mash_reference == 'archernm__actino785' | bacteria_dsDNA_nuc042_gene089$mash_query == 'archernm__actino785')
subcluster_archernm <- subset(bacteria_dsDNA_nuc020_gene062,bacteria_dsDNA_nuc020_gene062$mash_reference == 'archernm__actino785' | bacteria_dsDNA_nuc020_gene062$mash_query == 'archernm__actino785')

cluster_eagleeye <- subset(bacteria_dsDNA_nuc042_gene089,bacteria_dsDNA_nuc042_gene089$mash_reference == 'eagleeye__actino785' | bacteria_dsDNA_nuc042_gene089$mash_query == 'eagleeye__actino785')
subcluster_eagleeye <- subset(bacteria_dsDNA_nuc020_gene062,bacteria_dsDNA_nuc020_gene062$mash_reference == 'eagleeye__actino785' | bacteria_dsDNA_nuc020_gene062$mash_query == 'eagleeye__actino785')



#Immunity analysis
#Import immunity data that has already been filtered for immunity tests involving phages not present in the merged2333 mash data
#Both the mash table and immunity data table contain non-redundant rows, and the order of the reference and query phage names may not match the order used for the immunity data. 
#To solve this, create a new mash data table that contains duplicated data. Each duplicate will contain the reference and query phage names in reversed order.
#Then match this new table to the immunity data, and only keep matched rows. If successfull, there should be no loss of immunity data rows

immunity_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170117_immunity_data.csv",sep=",",header=TRUE)
names(immunity_data) <- c("immunity_lysogen_lysate","immunity_lysogen","immunity_lysate","immunity_data")
immunity_data$immunity_data <- as.factor(immunity_data$immunity_data)

temp <- mash_table2
temp$mash_query_ref <- paste(temp$mash_query,temp$mash_reference,sep="_")
temp$mash_query_ref <- as.factor(temp$mash_query_ref)
mash_match_a <- subset(temp,select = c('mash_ref_query','modified_mash_distance','pham_pham_dissimilarity'))
mash_match_b <- subset(temp,select = c('mash_query_ref','modified_mash_distance','pham_pham_dissimilarity'))
names(mash_match_a) <- c('phage_comparison','modified_mash_distance','pham_pham_dissimilarity')
names(mash_match_b) <- c('phage_comparison','modified_mash_distance','pham_pham_dissimilarity')
mash_match <- rbind(mash_match_a,mash_match_b)

immunity_analysis <- merge(immunity_data,mash_match,by.x="immunity_lysogen_lysate",by.y="phage_comparison")

immunity_analysis_1 <- subset(immunity_analysis,immunity_analysis$immunity_data == '1')
immunity_analysis_2 <- subset(immunity_analysis,immunity_analysis$immunity_data == '2')
immunity_analysis_3 <- subset(immunity_analysis,immunity_analysis$immunity_data == '3')
immunity_analysis_4 <- subset(immunity_analysis,immunity_analysis$immunity_data == '4')
immunity_analysis_5 <- subset(immunity_analysis,immunity_analysis$immunity_data == '5')


par(mar=c(4,8,4,4))
plot(immunity_analysis_1$modified_mash_distance,immunity_analysis_1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
par(new=TRUE)
plot(immunity_analysis_2$modified_mash_distance,immunity_analysis_2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light grey")
par(new=TRUE)
plot(immunity_analysis_4$modified_mash_distance,immunity_analysis_4$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(immunity_analysis_5$modified_mash_distance,immunity_analysis_5$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(immunity_analysis_3$modified_mash_distance,immunity_analysis_3$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,2,4,15))
stripchart(immunity_analysis$modified_mash_distance ~ immunity_analysis$immunity_data,las=1,col="black",xlim=c(0,0.5),pch=20,cex=0.5,xlab="",ylab="",xaxt="n",cex.axis=0.5)
axis(1,cex.axis=2)

#output matched data

write.table(immunity_analysis,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/r_filtered_networks_for_cytoscape/20170117_matched_immunity_data.csv",sep=",",row.names = FALSE,quote=FALSE)






###Investigate Cluster A subcluster diversities

subcluster_a1 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A1')
subcluster_a2 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A2')
subcluster_a3 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A3')
subcluster_a4 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A4')
subcluster_a5 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A5')
subcluster_a6 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A6')
subcluster_a7 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A7')
subcluster_a8 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A8')
subcluster_a9 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A9')
subcluster_a10 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A10')
subcluster_a11 <- subset(mash_table2,mash_table2$phage_subcluster_compare == 'A11')

par(mar=c(4,8,4,4))
plot(subcluster_a1$modified_mash_distance,subcluster_a1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a2$modified_mash_distance,subcluster_a2$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a3$modified_mash_distance,subcluster_a3$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a4$modified_mash_distance,subcluster_a4$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a5$modified_mash_distance,subcluster_a5$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a6$modified_mash_distance,subcluster_a6$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a7$modified_mash_distance,subcluster_a7$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a8$modified_mash_distance,subcluster_a8$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a9$modified_mash_distance,subcluster_a9$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a10$modified_mash_distance,subcluster_a10$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(subcluster_a11$modified_mash_distance,subcluster_a11$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





###Extrachromosomal phages
extrachromosomal_data <- mash_table2


cluster_actino$subcluster_A1_one_or_two <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_neither <- ifelse((cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" | cluster_actino$query_phage_subcluster == "A1") == TRUE,TRUE,FALSE)
cluster_actino$subcluster_A1_both <- ifelse((cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == FALSE | is.na(cluster_actino$ref_phage_subcluster == "A1" & cluster_actino$query_phage_subcluster == "A1") == TRUE,FALSE,TRUE)
cluster_actino$subcluster_A1_one <- ifelse(cluster_actino$subcluster_A1_one_or_two == TRUE & cluster_actino$subcluster_A1_both == FALSE,TRUE,FALSE)





extrachromosomal_data$extra_chrome_one_or_two <- ifelse((extrachromosomal_data$ref_extra_chrome == "yes" | extrachromosomal_data$query_extra_chrome == "yes") == FALSE | is.na(extrachromosomal_data$ref_extra_chrome == "yes" | extrachromosomal_data$query_extra_chrome == "yes") == TRUE,FALSE,TRUE)
extrachromosomal_data$extra_chrome_both <- ifelse((extrachromosomal_data$ref_extra_chrome == "yes" & extrachromosomal_data$query_extra_chrome == "yes") == FALSE | is.na(extrachromosomal_data$ref_extra_chrome == "yes" & extrachromosomal_data$query_extra_chrome == "yes") == TRUE,FALSE,TRUE)
extrachromosomal_data$extra_chrome_one <- ifelse(extrachromosomal_data$extra_chrome_one_or_two == TRUE & extrachromosomal_data$extra_chrome_both == FALSE,TRUE,FALSE)







#list of all non-Cluster A extrachromosomal phages in this dataset
#n15__nc_001901
#p1__nc_005856
#phihap-1__nc_010342
#py54__nc_005069
#vp58-5__nc_027981
#pzl12__actino785

extrachromosomal_data$nonClusterA_extra_chrome_one_or_two <- ifelse(extrachromosomal_data$mash_reference == 'n15__nc_001901' | 
                                                                    extrachromosomal_data$mash_reference == 'p1__nc_005856' | 
                                                                    extrachromosomal_data$mash_reference == 'phihap-1__nc_010342' | 
                                                                    extrachromosomal_data$mash_reference == 'py54__nc_005069' | 
                                                                    extrachromosomal_data$mash_reference == 'vp58-5__nc_027981' | 
                                                                    extrachromosomal_data$mash_reference == 'pzl12__actino785' | 
                                                                    extrachromosomal_data$mash_query == 'n15__nc_001901' | 
                                                                    extrachromosomal_data$mash_query == 'p1__nc_005856' | 
                                                                    extrachromosomal_data$mash_query == 'phihap-1__nc_010342' | 
                                                                    extrachromosomal_data$mash_query == 'py54__nc_005069' | 
                                                                    extrachromosomal_data$mash_query == 'vp58-5__nc_027981' |
                                                                    extrachromosomal_data$mash_query == 'pzl12__actino785',TRUE,FALSE)



extra_chrome_one_or_two_table <- subset(extrachromosomal_data,extrachromosomal_data$extra_chrome_one_or_two == TRUE)
extra_chrome_both_table <- subset(extrachromosomal_data,extrachromosomal_data$extra_chrome_both == TRUE)
extra_chrome_one_table <- subset(extrachromosomal_data,extrachromosomal_data$extra_chrome_one == TRUE)



non_clusterA_extra_chrome_one_or_two_table <- subset(extra_chrome_one_or_two_table,extra_chrome_one_or_two_table$phage_cluster_compare != "A" | is.na(extra_chrome_one_or_two_table$phage_cluster_compare))
#OR
non_clusterA_extra_chrome_one_or_two_table <- subset(extrachromosomal_data,extrachromosomal_data$nonClusterA_extra_chrome_one_or_two == TRUE)



par(mar=c(4,8,4,4))
plot(extra_chrome_one_or_two$modified_mash_distance,extra_chrome_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(extra_chrome_both$modified_mash_distance,extra_chrome_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(extra_chrome_one$modified_mash_distance,extra_chrome_one$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(non_clusterA_extra_chrome_one_or_two_table$modified_mash_distance,non_clusterA_extra_chrome_one_or_two_table$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")










###Noon seminar

#phage abundance by host
pie(summary(host_group_table$host_phylum),main="Phage diversity by host phylum",clockwise=TRUE)
pie(summary(host_group_table$host_class),main="Phage diversity by host class",clockwise=TRUE)
pie(summary(host_group_table$host_order),main="Phage diversity by host order",clockwise=TRUE)
pie(summary(host_group_table$host_family),main="Phage diversity by host family",clockwise=TRUE)
pie(summary(host_group_table$host_genus),main="Phage diversity by host genus",clockwise=TRUE)






#Phage identifier for lambda:lambda__nc_001416
#Copy main data table then determine which comparisons involve lambda. Since no self comparisons are present in the dataset, only need to compute if there's 'one' lambda, instead of 'one_or_two' or 'both'

#bobi__actino785
#cerasum__actino785
#rhyno__actino785
#bactobuster__actino785
#lambda__nc_001416
#kubed__actino785, cluster = BU
#shedlockholmes__actino785, cluster = K


plot_specific_phage_comparisons <- function(table,phage,cluster){
  
  
  phage_figure <- table
  phage_figure$phage_one <- ifelse(phage_figure$mash_reference == phage | phage_figure$mash_query == phage,TRUE,FALSE)
  phage_figure$cluster_same <- ifelse(phage_figure$phage_cluster_compare == cluster,TRUE,FALSE)
  phage_figure$cluster_diff <- ifelse(phage_figure$phage_cluster_compare != cluster,TRUE,FALSE)
  
  phage_comparisons_all <- subset(phage_figure,phage_figure$phage_one == TRUE)
  
  #Now use only those comparisons that have been clustered by SEA-PHAGES
  phage_comparisons_clustered <- subset(phage_comparisons_all,phage_comparisons_all$phage_cluster_source_compare == 'actino')
  phage_comparisons_cluster_same <- subset(phage_comparisons_clustered,phage_comparisons_clustered$cluster_same == TRUE)
  phage_comparisons_cluster_diff <- subset(phage_comparisons_clustered,phage_comparisons_clustered$cluster_diff == TRUE)  
  
  print('All comparisons')
  print(nrow(phage_comparisons_all))
  
  print('All clustered comparisons')
  print(nrow(phage_comparisons_clustered))

  print('Cluster same comparisons')
  print(nrow(phage_comparisons_cluster_same))
  
  print('Cluster diff comparisons')
  print(nrow(phage_comparisons_cluster_diff))
  
  par(mar=c(4,8,4,4))
  plot(phage_comparisons_all$modified_mash_distance,phage_comparisons_all$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(phage_comparisons_cluster_same$modified_mash_distance,phage_comparisons_cluster_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
  par(new=TRUE)
  plot(phage_comparisons_cluster_diff$modified_mash_distance,phage_comparisons_cluster_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
  abline(0,2,lty=2,lwd=3,col="grey")
  
}

dev.off()
plot_specific_phage_comparisons(mash_table2,'bobi__actino785','F')

dev.off()
plot_specific_phage_comparisons(mash_table2,'bactobuster__actino785','A')

dev.off()
plot_specific_phage_comparisons(mash_table2,'rhyno__actino785','A')

dev.off()
plot_specific_phage_comparisons(mash_table2,'lambda__nc_001416','')

dev.off()
plot_specific_phage_comparisons(mash_table2,'kubed__actino785','BU')

dev.off()
plot_specific_phage_comparisons(mash_table2,'shedlockholmes__actino785','K')





#split up SEA-PHAGES and other data

bac_dsDNA <- subset(mash_table2,mash_table2$host_superkingdom_compare == 'Bacteria' & mash_table2$phage_viral_type_compare == 'dsDNA')
sea_phages <- subset(bac_dsDNA,bac_dsDNA$phage_cluster_source_compare == 'actino')
non_sea_phages <- subset(bac_dsDNA,bac_dsDNA$phage_cluster_source_compare != 'actino' | is.na(bac_dsDNA$phage_cluster_source_compare) == TRUE)


par(mar=c(4,8,4,4))
plot(sea_phages$modified_mash_distance,sea_phages$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(non_sea_phages$modified_mash_distance,non_sea_phages$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





#Tail type diff
#Bacterial dsDNA phages are myo, podo, sipho, or tecti
tail <- subset(mash_table2,mash_table2$host_superkingdom_compare == 'Bacteria' & mash_table2$phage_viral_type_compare == 'dsDNA')
tail <- subset(tail,tail$phage_order_compare == 'Caudovirales')
tail <- subset(tail,tail$phage_family_compare == 'different')


par(mar=c(4,8,4,4))
plot(tail$modified_mash_distance,tail$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




















###VOG analysis





vog_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170221_shared_vog_proportion_data.csv",sep=",",header=TRUE)

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
#So when merged to mash_table2, no need to keep all rows in either table - it is expected there will be fewer rows than in both tables
#Also, there are 2 genomes (vb_paem_c1-14-ab28__NC_026600 and pv94__NC_027368) that contain no annotated genes. These are not in the pham data,
#so even though they are present in the VOG data, I am unable to compare these two genomes. This results in 1875 genomes.
mash_table2_vog <- merge(mash_table2,vog_table,by.x="mash_ref_query",by.y="vog_ref_query")
mash_table2_vog$mash_reference <- factor(mash_table2_vog$mash_reference)
mash_table2_vog$mash_query <- factor(mash_table2_vog$mash_query)



bacteria_dsDNA <- subset(mash_table2_vog,mash_table2_vog$host_superkingdom_compare == 'Bacteria' & mash_table2_vog$phage_viral_type_compare == 'dsDNA')
temperate_both <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == 'yes')
temperate_neither <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == 'no')





#How well do pham-based gcd and vog-based gcd correlate?
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$pham_pham_dissimilarity,bacteria_dsDNA$vog_gene_content_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")


#Compare pham-based and vog-based bacteria dsDNA phage plots 
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Compare pham-based and vog-based bacteria dsDNA phage lifestyle plots 
par(mar=c(4,8,4,4))
plot(temperate_both$modified_mash_distance,temperate_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_both$modified_mash_distance,temperate_both$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(temperate_neither$modified_mash_distance,temperate_neither$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_neither$modified_mash_distance,temperate_neither$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




# #Try to track down what causes the discrepancies between VOG and pham gene content dissimilarity
# #Filter out comparisons involving phages with fewer than X comparisons below modified mash distance < Y
# library(plyr)
# 
# #Obtain list of all unique phages in the combined VOG-Pham dataset
# vog_phages <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170221_vog_phage_identifiers.csv",sep=",",header=TRUE)
# 
# 
# 
# 
# #Create unique lists of phages at different stages of the analysis
# mash_table2_vog_phages <- rbind(as.character(mash_table2_vog$mash_reference),as.character(mash_table2_vog$mash_query))
# mash_table2_vog_phages <- factor(mash_table2_vog_phages)
# 
# bacteria_dsDNA_phages <- rbind(as.character(bacteria_dsDNA$mash_reference),as.character(bacteria_dsDNA$mash_query))
# bacteria_dsDNA_phages <- factor(bacteria_dsDNA_phages)
# 
# 
# vog_analysis_nuc010 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.1)
# 
# vog_analysis_nuc010_phages <- rbind(as.character(vog_analysis_nuc010$mash_reference),as.character(vog_analysis_nuc010$mash_query))
# vog_analysis_nuc010_phages <- factor(vog_analysis_nuc010_phages)
# 
# vog_analysis_nuc010_phages_count <- count(vog_analysis_nuc010_phages)
# names(vog_analysis_nuc010_phages_count) <- c('phage_identifier','frequency')
# vog_analysis_nuc010_phages_3relatives <- subset(vog_analysis_nuc010_phages_count,vog_analysis_nuc010_phages_count$freq > 2)
# 
# bacteria_dsDNA_nuc010relatives <- subset(bacteria_dsDNA,bacteria_dsDNA$mash_reference %in% vog_analysis_nuc010_phages_3relatives$phage_identifier & bacteria_dsDNA$mash_query %in% vog_analysis_nuc010_phages_3relatives$phage_identifier)
# 
# par(mar=c(4,8,4,4))
# plot(bacteria_dsDNA_nuc010relatives$modified_mash_distance,bacteria_dsDNA_nuc010relatives$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacteria_dsDNA_nuc010relatives$modified_mash_distance,bacteria_dsDNA_nuc010relatives$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# vog_analysis_nuc010_phages_10relatives <- subset(vog_analysis_nuc010_phages_count,vog_analysis_nuc010_phages_count$freq > 9)
# 
# bacteria_dsDNA_nuc010_10relatives <- subset(bacteria_dsDNA,bacteria_dsDNA$mash_reference %in% vog_analysis_nuc010_phages_10relatives$phage_identifier &
#                                                            bacteria_dsDNA$mash_query %in% vog_analysis_nuc010_phages_10relatives$phage_identifier)
# 
# par(mar=c(4,8,4,4))
# plot(bacteria_dsDNA_nuc010_10relatives$modified_mash_distance,bacteria_dsDNA_nuc010_10relatives$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacteria_dsDNA_nuc010_10relatives$modified_mash_distance,bacteria_dsDNA_nuc010_10relatives$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# vog_analysis_nuc005 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.05)
# vog_analysis_nuc005_phages <- rbind(as.character(vog_analysis_nuc005$mash_reference),as.character(vog_analysis_nuc005$mash_query))
# vog_analysis_nuc005_phages <- factor(vog_analysis_nuc005_phages)
# 
# vog_analysis_nuc005_phages_count <- count(vog_analysis_nuc005_phages)
# names(vog_analysis_nuc005_phages_count) <- c('phage_identifier','frequency')
# 
# 
# 
# vog_analysis_nuc005_phages_10relatives <- subset(vog_analysis_nuc005_phages_count,vog_analysis_nuc005_phages_count$freq > 9)
# 
# bacteria_dsDNA_nuc005_10relatives <- subset(bacteria_dsDNA,bacteria_dsDNA$mash_reference %in% vog_analysis_nuc005_phages_10relatives$phage_identifier &
#                                               bacteria_dsDNA$mash_query %in% vog_analysis_nuc005_phages_10relatives$phage_identifier)
# 
# par(mar=c(4,8,4,4))
# plot(bacteria_dsDNA_nuc005_10relatives$modified_mash_distance,bacteria_dsDNA_nuc005_10relatives$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(bacteria_dsDNA_nuc005_10relatives$modified_mash_distance,bacteria_dsDNA_nuc005_10relatives$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")



#Add phylogenetic data if needed
phylogeny_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170303_phylogeny_data.csv",sep=",",header=TRUE)


phylogeny_vog_analysis <- merge(bacteria_dsDNA,phylogeny_data,by.x="mash_ref_query",by.y="ref_query")

phylogeny_vog_hgcf <- subset(phylogeny_vog_analysis,phylogeny_vog_analysis$phage_cluster_compare == 'F' |
                           (phylogeny_vog_analysis$phage_cluster_compare == 'A' & phylogeny_vog_analysis$phage_subcluster_compare == 'A1'))


phylogeny_vog_lgcf <- subset(phylogeny_vog_analysis,
                             phylogeny_vog_analysis$phage_cluster_compare == 'K' |
                               phylogeny_vog_analysis$phage_cluster_compare == 'BD' |
                           (phylogeny_vog_analysis$phage_cluster_compare == 'A' & phylogeny_vog_analysis$phage_subcluster_compare != 'A1'))


phylogeny_vog_lytic <- subset(phylogeny_vog_analysis,phylogeny_vog_analysis$phage_cluster_compare == 'B')


par(mar=c(4,8,4,4))
plot(phylogeny_vog_hgcf$modified_mash_distance,phylogeny_vog_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_vog_lgcf$modified_mash_distance,phylogeny_vog_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_vog_lytic$modified_mash_distance,phylogeny_vog_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(phylogeny_vog_hgcf$phylogeny_distance,phylogeny_vog_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_vog_lgcf$phylogeny_distance,phylogeny_vog_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_vog_lytic$phylogeny_distance,phylogeny_vog_lytic$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


par(mar=c(4,8,4,4))
plot(phylogeny_vog_hgcf$phylogeny_distance,phylogeny_vog_hgcf$vog_gene_content_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_vog_lgcf$phylogeny_distance,phylogeny_vog_lgcf$vog_gene_content_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_vog_lytic$phylogeny_distance,phylogeny_vog_lytic$vog_gene_content_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")











###Analyze predicted lifestyle data


lifestyle_analysis <- mash_table2


bacteria_dsDNA <- subset(lifestyle_analysis,lifestyle_analysis$host_superkingdom_compare == 'Bacteria' & lifestyle_analysis$phage_viral_type_compare == 'dsDNA')


lifestyle_predicted_temperate <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'yes')
lifestyle_predicted_lytic <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'no')
lifestyle_predicted_both <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'different')


par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate$modified_mash_distance,lifestyle_predicted_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic$modified_mash_distance,lifestyle_predicted_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(lifestyle_predicted_both$modified_mash_distance,lifestyle_predicted_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



# #Directly compare to actual lifestyle data
# lifestyle_actual_temperate <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == 'yes')
# lifestyle_actual_lytic <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == 'no')
# lifestyle_actual_both <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_temperate_compare == 'different')
# 
# par(mar=c(4,8,4,4))
# plot(lifestyle_actual_temperate$modified_mash_distance,lifestyle_actual_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(lifestyle_actual_lytic$modified_mash_distance,lifestyle_actual_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")
# 
# par(mar=c(4,8,4,4))
# plot(lifestyle_actual_both$modified_mash_distance,lifestyle_actual_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
# abline(0,2,lty=2,lwd=3,col="grey")











lifestyle_predicted_temperate_actino <- subset(lifestyle_predicted_temperate,lifestyle_predicted_temperate$host_phylum_compare == 'Actinobacteria')
lifestyle_predicted_temperate_cyano <- subset(lifestyle_predicted_temperate,lifestyle_predicted_temperate$host_phylum_compare == 'Cyanobacteria')
lifestyle_predicted_temperate_firm <- subset(lifestyle_predicted_temperate,lifestyle_predicted_temperate$host_phylum_compare == 'Firmicutes')
lifestyle_predicted_temperate_proteo <- subset(lifestyle_predicted_temperate,lifestyle_predicted_temperate$host_phylum_compare == 'Proteobacteria')

lifestyle_predicted_lytic_actino <- subset(lifestyle_predicted_lytic,lifestyle_predicted_lytic$host_phylum_compare == 'Actinobacteria')
lifestyle_predicted_lytic_cyano <- subset(lifestyle_predicted_lytic,lifestyle_predicted_lytic$host_phylum_compare == 'Cyanobacteria')
lifestyle_predicted_lytic_firm <- subset(lifestyle_predicted_lytic,lifestyle_predicted_lytic$host_phylum_compare == 'Firmicutes')
lifestyle_predicted_lytic_proteo <- subset(lifestyle_predicted_lytic,lifestyle_predicted_lytic$host_phylum_compare == 'Proteobacteria')


par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate_actino$modified_mash_distance,lifestyle_predicted_temperate_actino$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate_cyano$modified_mash_distance,lifestyle_predicted_temperate_cyano$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate_firm$modified_mash_distance,lifestyle_predicted_temperate_firm$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate_proteo$modified_mash_distance,lifestyle_predicted_temperate_proteo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic_actino$modified_mash_distance,lifestyle_predicted_lytic_actino$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic_cyano$modified_mash_distance,lifestyle_predicted_lytic_cyano$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic_firm$modified_mash_distance,lifestyle_predicted_lytic_firm$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic_proteo$modified_mash_distance,lifestyle_predicted_lytic_proteo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")









#To investigate conflicting lifestyle data
lifestyle_analysis$predicted_lifestyle_conflict <- ifelse(as.character(lifestyle_analysis$ref_predicted_temperate) == as.character(lifestyle_analysis$ref_phage_temperate) & as.character(lifestyle_analysis$query_predicted_temperate) == as.character(lifestyle_analysis$query_phage_temperate) | (is.na(lifestyle_analysis$ref_phage_temperate) == TRUE | is.na(lifestyle_analysis$query_phage_temperate) == TRUE),"same_lifestyle","conflicting_lifestyle")
lifestyle_analysis$predicted_lifestyle_conflict <- as.factor(lifestyle_analysis$predicted_lifestyle_conflict)

bacteria_dsDNA <- subset(lifestyle_analysis,lifestyle_analysis$host_superkingdom_compare == 'Bacteria' & lifestyle_analysis$phage_viral_type_compare == 'dsDNA')

lifestyle_no_conflict <- subset(bacteria_dsDNA,bacteria_dsDNA$predicted_lifestyle_conflict == 'same_lifestyle')

lifestyle_no_conflict_predicted_temperate <- subset(lifestyle_no_conflict,lifestyle_no_conflict$predicted_lifestyle_compare == 'yes')
lifestyle_no_conflict_predicted_lytic <- subset(lifestyle_no_conflict,lifestyle_no_conflict$predicted_lifestyle_compare == 'no')
lifestyle_no_conflict_predicted_both <- subset(lifestyle_no_conflict,lifestyle_no_conflict$predicted_lifestyle_compare == 'different')



par(mar=c(4,8,4,4))
plot(lifestyle_no_conflict_predicted_temperate$modified_mash_distance,lifestyle_no_conflict_predicted_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_no_conflict_predicted_lytic$modified_mash_distance,lifestyle_no_conflict_predicted_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(lifestyle_no_conflict_predicted_both$modified_mash_distance,lifestyle_no_conflict_predicted_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


















###Phylogeny comparison
#phylogeny_data_round1 <- phylogeny_data
phylogeny_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170317_phylogeny_data.csv",sep=",",header=TRUE)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Decide which table to merge with...
phylogeny_analysis <- merge(mash_table2,phylogeny_data,by.x="mash_ref_query",by.y="ref_query")
#OR...
phylogeny_analysis <- merge(gene_specific_mash_table,phylogeny_data,by.x="mash_ref_query",by.y="ref_query")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#Compare mash to phylogeny
par(mar=c(4,8,4,4))
plot(phylogeny_analysis$modified_mash_distance,phylogeny_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(phylogeny_analysis$phylogeny_distance,phylogeny_analysis$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


#Split into modes to visualize by color
phylogeny_analysis_hgcf <- subset(phylogeny_analysis,phylogeny_analysis$gene_flux_category == "high" & phylogeny_analysis$phage_predicted_temperate_compare == "yes")
phylogeny_analysis_lgcf <- subset(phylogeny_analysis,phylogeny_analysis$gene_flux_category == "low" & phylogeny_analysis$phage_predicted_temperate_compare == "yes")
phylogeny_analysis_lytic <- subset(phylogeny_analysis,phylogeny_analysis$phage_predicted_temperate_compare == "no")

par(mar=c(4,8,4,4))
plot(phylogeny_analysis_hgcf$modified_mash_distance,phylogeny_analysis_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_analysis_lgcf$modified_mash_distance,phylogeny_analysis_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_analysis_lytic$modified_mash_distance,phylogeny_analysis_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

par(mar=c(4,8,4,4))
plot(phylogeny_analysis_hgcf$phylogeny_distance,phylogeny_analysis_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_analysis_lgcf$phylogeny_distance,phylogeny_analysis_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_analysis_lytic$phylogeny_distance,phylogeny_analysis_lytic$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")








#Assess Cluster A data
phylogeny_analysis_a <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'A')
phylogeny_analysis_a_both_a1 <- subset(phylogeny_analysis_a,phylogeny_analysis_a$phage_subcluster_compare == 'A1')
phylogeny_analysis_a_both_nonA1 <- subset(phylogeny_analysis_a,phylogeny_analysis_a$phage_subcluster_compare != 'A1')
phylogeny_analysis_a_one_a1 <- subset(phylogeny_analysis_a,(phylogeny_analysis_a$ref_phage_subcluster == 'A1' | phylogeny_analysis_a$query_phage_subcluster == 'A1') & phylogeny_analysis_a$phage_subcluster_compare == "different")


par(mar=c(4,8,4,4))
plot(phylogeny_analysis_a_both_nonA1$modified_mash_distance,phylogeny_analysis_a_both_nonA1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_analysis_a_both_a1$modified_mash_distance,phylogeny_analysis_a_both_a1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(phylogeny_analysis_a_one_a1$modified_mash_distance,phylogeny_analysis_a_one_a1$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")

par(mar=c(4,8,4,4))
plot(phylogeny_analysis_a_both_nonA1$phylogeny_distance,phylogeny_analysis_a_both_nonA1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_analysis_a_both_a1$phylogeny_distance,phylogeny_analysis_a_both_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(phylogeny_analysis_a_one_a1$phylogeny_distance,phylogeny_analysis_a_one_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")


par(mar=c(4,8,4,4))
plot(phylogeny_analysis_a_both_nonA1$phylogeny_distance,phylogeny_analysis_a_both_nonA1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_analysis_a_both_a1$phylogeny_distance,phylogeny_analysis_a_both_a1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(phylogeny_analysis_a_one_a1$phylogeny_distance,phylogeny_analysis_a_one_a1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")

par(mar=c(4,8,4,4))
plot(phylogeny_analysis_a_both_a1$phylogeny_distance,phylogeny_analysis_a_both_a1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="cyan")
par(new=TRUE)
plot(phylogeny_analysis_a_one_a1$phylogeny_distance,phylogeny_analysis_a_one_a1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="purple")
par(new=TRUE)
plot(phylogeny_analysis_a_both_nonA1$phylogeny_distance,phylogeny_analysis_a_both_nonA1$pham_pham_dissimilarity,xlim=c(0,0.2),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")

  






  
#match up pham proportion data that contains orpham count, pham distribution, etc.
#I designate these columns as "pham2" since it is the second set of pham data that I have loaded so far. 
#This second pham data overlaps the first set, but contains additional columns. Also, the extra columns are impacted by the data subset - it only contains comparisons from the actino782 set. So it is not a replacement for the first pham data file.
#If the analysis works well, the updated pham data can be loaded at the beginning, and this section removed.


#Actino pham data import code IF the data does not have pham-function data
# actino_pham_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170315_bac2333_actino_only_pairwise_pham_proportions.csv",sep=",",header=TRUE)
# actino_pham_data$pham_dissimilarity <- 1 - actino_pham_data$average_shared_proportion
# actino_pham_data$jaccard_dissimilarity <- 1 - actino_pham_data$jaccard_similarity
# names(actino_pham_data) <- c("pham2_phage1","pham2_phage1_number_unshared_phams","pham2_phage1_shared_proportion","pham2_phage2",
#                              "pham2_phage2_number_unshared_phams","pham2_phage2_shared_proportion","pham2_number_shared_phams","pham2_averaged_shared_proportion",
#                              "pham2_jaccard_similarity","pham2_shared_pham_distribution_mean","pham2_shared_pham_distribution_median","pham2_shared_pham_distribution_max",
#                              "pham2_unshared_pham_distribution_mean","pham2_unshared_pham_distribution_median",
#                              "pham2_unshared_pham_distribution_max","pham2_unshared_orpham_count","pham2_pham_dissimilarity","pham2_jaccard_dissimilarity")  
# #Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
# actino_pham_data$pham2_phage1_phage2 <- paste(actino_pham_data$pham2_phage1,"_",actino_pham_data$pham2_phage2,sep="")
# actino_pham_data$pham2_phage1_phage2 <- as.factor(actino_pham_data$pham2_phage1_phage2)



actino_pham_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170315_bac2333_actino_only_pairwise_pham_proportions.csv",sep=",",header=TRUE)

names(actino_pham_data) <- c("pham2_phage1","pham2_phage1_number_unshared_phams","pham2_phage1_shared_proportion","pham2_phage2",
                             "pham2_phage2_number_unshared_phams","pham2_phage2_shared_proportion","pham2_number_shared_phams","pham2_average_shared_proportion",
                             "pham2_jaccard_similarity","pham2_shared_pham_distribution_mean","pham2_shared_pham_distribution_median","pham2_shared_pham_distribution_max",
                             "pham2_unshared_pham_distribution_mean","pham2_unshared_pham_distribution_median",
                             "pham2_unshared_pham_distribution_max","pham2_unshared_orpham_count",
                             "pham2_Unspecified_phage1_number_of_unshared_phams","pham2_Unspecified_phage2_number_of_unshared_phams","pham2_Unspecified_number_of_shared_phams","pham2_Unspecified_average_shared_proportion",
                             "pham2_defense_phage1_number_of_unshared_phams","pham2_defense_phage2_number_of_unshared_phams","pham2_defense_number_of_shared_phams","pham2_defense_average_shared_proportion",
                             "pham2_dna_metabolism_phage1_number_of_unshared_phams","pham2_dna_metabolism_phage2_number_of_unshared_phams","pham2_dna_metabolism_number_of_shared_phams","pham2_dna_metabolism_average_shared_proportion",
                             "pham2_lysis_phage1_number_of_unshared_phams","pham2_lysis_phage2_number_of_unshared_phams","pham2_lysis_number_of_shared_phams","pham2_lysis_average_shared_proportion",
                             "pham2_lysogeny_phage1_number_of_unshared_phams","pham2_lysogeny_phage2_number_of_unshared_phams","pham2_lysogeny_number_of_shared_phams","pham2_lysogeny_average_shared_proportion",
                             "pham2_mobile_phage1_number_of_unshared_phams","pham2_mobile_phage2_number_of_unshared_phams","pham2_mobile_number_of_shared_phams","pham2_mobile_average_shared_proportion",
                             "pham2_other_phage1_number_of_unshared_phams","pham2_other_phage2_number_of_unshared_phams","pham2_other_number_of_shared_phams","pham2_other_average_shared_proportion",
                             "pham2_recombination_replication_phage1_number_of_unshared_phams","pham2_recombination_replication_phage2_number_of_unshared_phams","pham2_recombination_replication_number_of_shared_phams","pham2_recombination_replication_average_shared_proportion",
                             "pham2_structure_assembly_phage1_number_of_unshared_phams","pham2_structure_assembly_phage2_number_of_unshared_phams","pham2_structure_assembly_number_of_shared_phams","pham2_structure_assembly_average_shared_proportion")



#The function-specific average shared proportion data is not exactly accurate. If there were no shared or unshared phams in the function category, technically you can't compute a shared proportion
#(since dividing by zero), so a "-1" was inserted in these cases. But really, they should be NA.
actino_pham_data$pham2_Unspecified_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_Unspecified_average_shared_proportion != -1,actino_pham_data$pham2_Unspecified_average_shared_proportion,NA)
actino_pham_data$pham2_defense_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_defense_average_shared_proportion != -1,actino_pham_data$pham2_defense_average_shared_proportion,NA)
actino_pham_data$pham2_dna_metabolism_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_dna_metabolism_average_shared_proportion != -1,actino_pham_data$pham2_dna_metabolism_average_shared_proportion,NA)
actino_pham_data$pham2_lysis_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_lysis_average_shared_proportion != -1,actino_pham_data$pham2_lysis_average_shared_proportion,NA)
actino_pham_data$pham2_lysogeny_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_lysogeny_average_shared_proportion != -1,actino_pham_data$pham2_lysogeny_average_shared_proportion,NA)
actino_pham_data$pham2_mobile_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_mobile_average_shared_proportion != -1,actino_pham_data$pham2_mobile_average_shared_proportion,NA)
actino_pham_data$pham2_other_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_other_average_shared_proportion != -1,actino_pham_data$pham2_other_average_shared_proportion,NA)
actino_pham_data$pham2_recombination_replication_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_recombination_replication_average_shared_proportion != -1,actino_pham_data$pham2_recombination_replication_average_shared_proportion,NA)
actino_pham_data$pham2_structure_assembly_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_structure_assembly_average_shared_proportion != -1,actino_pham_data$pham2_structure_assembly_average_shared_proportion,NA)




#Compute gene content dissimilarity
actino_pham_data$pham2_pham_dissimilarity <- 1 - actino_pham_data$pham2_average_shared_proportion
actino_pham_data$pham2_jaccard_dissimilarity <- 1 - actino_pham_data$pham2_jaccard_similarity
actino_pham_data$pham2_pham_dissimilarity_Unspecified <- 1 - actino_pham_data$pham2_Unspecified_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_defense <- 1 - actino_pham_data$pham2_defense_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_dna_metabolism <- 1 - actino_pham_data$pham2_dna_metabolism_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_lysis <- 1 - actino_pham_data$pham2_lysis_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_lysogeny <- 1 - actino_pham_data$pham2_lysogeny_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_mobile <- 1 - actino_pham_data$pham2_mobile_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_other <- 1 - actino_pham_data$pham2_other_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_recomb_rep <- 1 - actino_pham_data$pham2_recombination_replication_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_structural <- 1 - actino_pham_data$pham2_structure_assembly_average_shared_proportion_mod


#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
actino_pham_data$pham2_phage1_phage2 <- paste(actino_pham_data$pham2_phage1,"_",actino_pham_data$pham2_phage2,sep="")
actino_pham_data$pham2_phage1_phage2 <- as.factor(actino_pham_data$pham2_phage1_phage2)





#To retain all rows, be sure to keep all.x=TRUE, but you don't want to retain all rows = making scatter plots or histograms can cause errors if not all rows have data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
phylogeny_analysis <- merge(phylogeny_analysis,actino_pham_data,by.x="mash_ref_query",by.y="pham2_phage1_phage2")





#Verify both sets of uploaded pham data match
# plot(phylogeny_analysis$pham_phage1_number_unshared_phams,phylogeny_analysis$pham2_phage1_number_unshared_phams)
# plot(phylogeny_analysis$pham_phage1_shared_proportion,phylogeny_analysis$pham2_phage1_shared_proportion)
# plot(phylogeny_analysis$pham_phage2_number_unshared_phams,phylogeny_analysis$pham2_phage2_number_unshared_phams)
# plot(phylogeny_analysis$pham_phage2_shared_proportion,phylogeny_analysis$pham2_phage2_shared_proportion)
# plot(phylogeny_analysis$pham_number_shared_phams,phylogeny_analysis$pham2_number_shared_phams)
# plot(phylogeny_analysis$pham_averaged_shared_proportion,phylogeny_analysis$pham2_average_shared_proportion)
# plot(phylogeny_analysis$pham_jaccard_similarity,phylogeny_analysis$pham2_jaccard_similarity)
# plot(phylogeny_analysis$pham_pham_dissimilarity,phylogeny_analysis$pham2_pham_dissimilarity)
# plot(phylogeny_analysis$pham_jaccard_dissimilarity,phylogeny_analysis$pham2_jaccard_dissimilarity)







#Import Count cumulative pham gain/loss data
count_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170321_count_wagner_analysis.csv",sep=",",header=TRUE)

names(count_data) <- c("ref_query","wagner_distance","wagner_pham_gain_sum","wagner_pham_loss_sum")


#Count data contains redundant data, so only half the rows should be retained. No need to specificy all.x or all.y.
phylogeny_analysis <- merge(phylogeny_analysis,count_data,by.x="mash_ref_query",by.y="ref_query")

#Verify branch distances from count data match those from original phylogeny data
plot(phylogeny_analysis$phylogeny_distance,phylogeny_analysis$wagner_distance)





#Now split by evolutionary mode
phylogeny_hgcf <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'F' |
                           (phylogeny_analysis$phage_cluster_compare == 'A' & phylogeny_analysis$phage_subcluster_compare == 'A1'))


phylogeny_lgcf <- subset(phylogeny_analysis,
                           phylogeny_analysis$phage_cluster_compare == 'K' |
                           phylogeny_analysis$phage_cluster_compare == 'BD' |
                           (phylogeny_analysis$phage_cluster_compare == 'A' & phylogeny_analysis$phage_subcluster_compare != 'A1'))


phylogeny_lytic <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'B')


par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$modified_mash_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$modified_mash_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$modified_mash_distance,phylogeny_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

# par(mar=c(4,8,4,4))
# plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")


par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$wagner_pham_gain_sum,xlim=c(0,2.5),ylim=c(0,250),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$wagner_pham_gain_sum,xlim=c(0,2.5),ylim=c(0,250),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$wagner_pham_gain_sum,xlim=c(0,2.5),ylim=c(0,250),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,100),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,100),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,100),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$wagner_pham_gain_sum + phylogeny_hgcf$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$wagner_pham_gain_sum + phylogeny_lgcf$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$wagner_pham_gain_sum + phylogeny_lytic$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,300),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$wagner_pham_gain_sum / phylogeny_hgcf$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,150),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$wagner_pham_gain_sum / phylogeny_lgcf$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,150),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$wagner_pham_gain_sum / phylogeny_lytic$wagner_pham_loss_sum,xlim=c(0,2.5),ylim=c(0,150),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,log(phylogeny_lytic$wagner_pham_gain_sum / phylogeny_lytic$wagner_pham_loss_sum,2),xlim=c(0,2.5),ylim=c(-4,7),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(phylogeny_hgcf$phylogeny_distance,log(phylogeny_hgcf$wagner_pham_gain_sum / phylogeny_hgcf$wagner_pham_loss_sum,2),xlim=c(0,2.5),ylim=c(-4,7),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,log(phylogeny_lgcf$wagner_pham_gain_sum / phylogeny_lgcf$wagner_pham_loss_sum,2),xlim=c(0,2.5),ylim=c(-4,7),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$wagner_pham_gain_sum,phylogeny_hgcf$wagner_pham_loss_sum,xlim=c(0,250),ylim=c(0,100),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$wagner_pham_gain_sum,phylogeny_lgcf$wagner_pham_loss_sum,xlim=c(0,250),ylim=c(0,100),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$wagner_pham_gain_sum,phylogeny_lytic$wagner_pham_loss_sum,xlim=c(0,250),ylim=c(0,100),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,0.4,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$wagner_pham_gain_sum,phylogeny_lgcf$wagner_pham_loss_sum,xlim=c(0,250),ylim=c(0,250),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$wagner_pham_gain_sum,phylogeny_lytic$wagner_pham_loss_sum,xlim=c(0,250),ylim=c(0,250),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
par(new=TRUE)
plot(phylogeny_hgcf$wagner_pham_gain_sum,phylogeny_hgcf$wagner_pham_loss_sum,xlim=c(0,250),ylim=c(0,250),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,1,lty=2,lwd=3,col="grey")




#Density distributions

hgcf_ratio_data <- log(phylogeny_hgcf$wagner_pham_gain_sum / phylogeny_hgcf$wagner_pham_loss_sum,2)
lgcf_ratio_data <- log(phylogeny_lgcf$wagner_pham_gain_sum / phylogeny_lgcf$wagner_pham_loss_sum,2)
lytic_ratio_data <- log(phylogeny_lytic$wagner_pham_gain_sum / phylogeny_lytic$wagner_pham_loss_sum,2)

hgcf_ratio_data <- hgcf_ratio_data[!is.na(hgcf_ratio_data)]
lgcf_ratio_data <- lgcf_ratio_data[!is.na(lgcf_ratio_data)]
lytic_ratio_data <- lytic_ratio_data[!is.na(lytic_ratio_data)]

par(mar=c(4,8,4,4))
plot(density(hgcf_ratio_data),xlim=c(-3,6),ylim=c(0,1),ann=FALSE,col="blue")
par(new=TRUE)
plot(density(lgcf_ratio_data),xlim=c(-3,6),ylim=c(0,1),ann=FALSE,,col="green")
par(new=TRUE)
plot(density(lytic_ratio_data),xlim=c(-3,6),ylim=c(0,1),ann=FALSE,,col="red")


hgcf_gain_sum_dist_ratio_data <- log(phylogeny_hgcf$wagner_pham_gain_sum / phylogeny_hgcf$wagner_pham_loss_sum,2) / phylogeny_hgcf$phylogeny_distance
lgcf_gain_sum_dist_ratio_data <- log(phylogeny_lgcf$wagner_pham_gain_sum / phylogeny_lgcf$wagner_pham_loss_sum,2) / phylogeny_lgcf$phylogeny_distance
lytic_gain_sum_dist_ratio_data <- log(phylogeny_lytic$wagner_pham_gain_sum / phylogeny_lytic$wagner_pham_loss_sum,2) / phylogeny_lytic$phylogeny_distance

hgcf_gain_sum_dist_ratio_data <- hgcf_gain_sum_dist_ratio_data[!is.na(hgcf_gain_sum_dist_ratio_data)]
lgcf_gain_sum_dist_ratio_data <- lgcf_gain_sum_dist_ratio_data[!is.na(lgcf_gain_sum_dist_ratio_data)]
lytic_gain_sum_dist_ratio_data <- lytic_gain_sum_dist_ratio_data[!is.na(lytic_gain_sum_dist_ratio_data)]


par(mar=c(4,8,4,4))
plot(density(hgcf_gain_sum_dist_ratio_data),ylim=c(0,1),ann=FALSE,col="blue")
par(new=TRUE)
plot(density(lgcf_gain_sum_dist_ratio_data),ylim=c(0,1),ann=FALSE,,col="green")
par(new=TRUE)
plot(density(lytic_gain_sum_dist_ratio_data),ylim=c(0,1),ann=FALSE,,col="red")







#Ratios

phylogeny_hgcf_ratios <- phylogeny_hgcf
phylogeny_lgcf_ratios <- phylogeny_lgcf
phylogeny_lytic_ratios <- phylogeny_lytic

phylogeny_hgcf_ratios$wagner_gain_rate <- phylogeny_hgcf_ratios$wagner_pham_gain_sum / phylogeny_hgcf_ratios$phylogeny_distance
phylogeny_lgcf_ratios$wagner_gain_rate <- phylogeny_lgcf_ratios$wagner_pham_gain_sum / phylogeny_lgcf_ratios$phylogeny_distance
phylogeny_lytic_ratios$wagner_gain_rate <- phylogeny_lytic_ratios$wagner_pham_gain_sum / phylogeny_lytic_ratios$phylogeny_distance

phylogeny_hgcf_ratios$wagner_loss_rate <- phylogeny_hgcf_ratios$wagner_pham_loss_sum / phylogeny_hgcf_ratios$phylogeny_distance
phylogeny_lgcf_ratios$wagner_loss_rate <- phylogeny_lgcf_ratios$wagner_pham_loss_sum / phylogeny_lgcf_ratios$phylogeny_distance
phylogeny_lytic_ratios$wagner_loss_rate <- phylogeny_lytic_ratios$wagner_pham_loss_sum / phylogeny_lytic_ratios$phylogeny_distance


phylogeny_hgcf_ratios <- phylogeny_hgcf_ratios[!is.na(phylogeny_hgcf_ratios$wagner_gain_rate) & !is.na(phylogeny_hgcf_ratios$wagner_loss_rate),]
phylogeny_lgcf_ratios <- phylogeny_lgcf_ratios[!is.na(phylogeny_lgcf_ratios$wagner_gain_rate) & !is.na(phylogeny_lgcf_ratios$wagner_loss_rate),]
phylogeny_lytic_ratios <- phylogeny_lytic_ratios[!is.na(phylogeny_lytic_ratios$wagner_gain_rate) & !is.na(phylogeny_lytic_ratios$wagner_loss_rate),]

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_ratios$wagner_gain_rate,phylogeny_hgcf_ratios$wagner_loss_rate,xlim=c(0,2e4),ylim=c(0,2e4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_ratios$wagner_gain_rate,phylogeny_lgcf_ratios$wagner_loss_rate,xlim=c(0,2e4),ylim=c(0,2e4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_ratios$wagner_gain_rate,phylogeny_lytic_ratios$wagner_loss_rate,xlim=c(0,2e4),ylim=c(0,2e4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_ratios$wagner_gain_rate,phylogeny_hgcf_ratios$wagner_loss_rate,xlim=c(0,5e4),ylim=c(0,5e4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_ratios$wagner_gain_rate,phylogeny_lgcf_ratios$wagner_loss_rate,xlim=c(0,5e4),ylim=c(0,5e4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_ratios$wagner_gain_rate,phylogeny_lytic_ratios$wagner_loss_rate,xlim=c(0,5e4),ylim=c(0,5e4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
abline(0,1,lty=2,lwd=3,col="grey")

hist(phylogeny_hgcf_ratios$wagner_gain_rate,xlim=c(0,1e4),breaks=100000,col="blue")
hist(phylogeny_lgcf_ratios$wagner_gain_rate,xlim=c(0,1e4),breaks=100000,col="green")
hist(phylogeny_lytic_ratios$wagner_gain_rate,xlim=c(0,1e4),breaks=100000,col="red")

plot(density(phylogeny_hgcf_ratios$wagner_gain_rate),xlim=c(0,2e5))
plot(density(phylogeny_lgcf_ratios$wagner_gain_rate),xlim=c(0,2e5))
plot(density(phylogeny_lytic_ratios$wagner_gain_rate),xlim=c(0,2e5))


plot(density(phylogeny_hgcf_ratios$wagner_gain_rate))


phylogeny_hgcf_ratios_gains <- subset(phylogeny_hgcf_ratios,select=c("wagner_gain_rate"))
phylogeny_hgcf_ratios_gains$mode <- "hgcf"

phylogeny_lgcf_ratios_gains <- subset(phylogeny_lgcf_ratios,select=c("wagner_gain_rate"))
phylogeny_lgcf_ratios_gains$mode <- "lgcf"

phylogeny_lytic_ratios_gains <- subset(phylogeny_lytic_ratios,select=c("wagner_gain_rate"))
phylogeny_lytic_ratios_gains$mode <- "lytic"

mode_gains <- rbind(phylogeny_hgcf_ratios_gains,phylogeny_lgcf_ratios_gains,phylogeny_lytic_ratios_gains)
boxplot(mode_gains$wagner_gain_rate ~ mode_gains$mode,ylim=c(0,1e4))



phylogeny_hgcf_ratios_losses <- subset(phylogeny_hgcf_ratios,select=c("wagner_loss_rate"))
phylogeny_hgcf_ratios_losses$mode <- "hgcf"

phylogeny_lgcf_ratios_losses <- subset(phylogeny_lgcf_ratios,select=c("wagner_loss_rate"))
phylogeny_lgcf_ratios_losses$mode <- "lgcf"

phylogeny_lytic_ratios_losses <- subset(phylogeny_lytic_ratios,select=c("wagner_loss_rate"))
phylogeny_lytic_ratios_losses$mode <- "lytic"

mode_losses <- rbind(phylogeny_hgcf_ratios_losses,phylogeny_lgcf_ratios_losses,phylogeny_lytic_ratios_losses)
boxplot(mode_losses$wagner_loss_rate ~ mode_losses$mode,ylim=c(0,6e3))


phylogeny_hgcf_ratios$wag
phylogeny_hgcf_ratios_gains_losses <- subset(phylogeny_hgcf_ratios,select=c("wagner_loss_rate"))
phylogeny_hgcf_ratios_gains_losses$mode <- "hgcf"

phylogeny_lgcf_ratios_gains_losses <- subset(phylogeny_lgcf_ratios,select=c("wagner_loss_rate"))
phylogeny_lgcf_ratios_gains_losses$mode <- "lgcf"

phylogeny_lytic_ratios_gains_losses <- subset(phylogeny_lytic_ratios,select=c("wagner_loss_rate"))
phylogeny_lytic_ratios_gains_losses$mode <- "lytic"

mode_losses <- rbind(phylogeny_hgcf_ratios_losses,phylogeny_lgcf_ratios_losses,phylogeny_lytic_ratios_losses)
boxplot(mode_losses$wagner_loss_rate ~ mode_losses$mode,ylim=c(0,6e3))








# hgcf_gain_rate <- phylogeny_hgcf$wagner_pham_gain_sum / phylogeny_hgcf$phylogeny_distance
# lgcf_gain_rate <- phylogeny_lgcf$wagner_pham_gain_sum / phylogeny_lgcf$phylogeny_distance
# lytic_gain_rate <- phylogeny_lytic$wagner_pham_gain_sum / phylogeny_lytic$phylogeny_distance
# 
# hgcf_gain_rate <- hgcf_gain_rate[!is.na(hgcf_gain_rate)]
# lgcf_gain_rate <- lgcf_gain_rate[!is.na(lgcf_gain_rate)]
# lytic_gain_rate <- lytic_gain_rate[!is.na(lytic_gain_rate)]
# 
# hgcf_loss_rate <- phylogeny_hgcf$wagner_pham_loss_sum / phylogeny_hgcf$phylogeny_distance
# lgcf_loss_rate <- phylogeny_lgcf$wagner_pham_loss_sum / phylogeny_lgcf$phylogeny_distance
# lytic_loss_rate <- phylogeny_lytic$wagner_pham_loss_sum / phylogeny_lytic$phylogeny_distance
# 
# hgcf_loss_rate <- hgcf_loss_rate[!is.na(hgcf_loss_rate)]
# lgcf_loss_rate <- lgcf_loss_rate[!is.na(lgcf_loss_rate)]
# lytic_loss_rate <- lytic_loss_rate[!is.na(lytic_loss_rate)]
# 
# par(mar=c(4,8,4,4))
# plot(hgcf_gain_rate,hgcf_loss_rate,pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# par(new=TRUE)
# plot(lgcf_gain_rate,lgcf_loss_rate,pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(lytic_gain_rate,lytic_loss_rate,pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# abline(0,1,lty=2,lwd=3,col="grey")
# 
# 
# 
# 
# par(mar=c(4,8,4,4))
# plot(phylogeny_lgcf$wagner_pham_gain_sum / phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$wagner_pham_loss_sum / phylogeny_lgcf$phylogeny_distance,xlim=c(0,),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
# par(new=TRUE)
# plot(phylogeny_lytic$wagner_pham_gain_sum / phylogeny_lytic$phylogeny_distance,phylogeny_lytic$wagner_pham_loss_sum / phylogeny_lytic$phylogeny_distance,pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")
# par(new=TRUE)
# plot(phylogeny_hgcf$wagner_pham_gain_sum / phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$wagner_pham_loss_sum / phylogeny_hgcf$phylogeny_distance,pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
# abline(0,1,lty=2,lwd=3,col="grey")














par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")



plot(phylogeny_hgcf$modified_mash_distance,phylogeny_hgcf$phylogeny_distance)
plot(phylogeny_lgcf$modified_mash_distance,phylogeny_lgcf$phylogeny_distance)

temp <- subset(phylogeny_hgcf,phylogeny_hgcf$phylogeny_distance > 0.5,select=c('mash_ref_query','mash_reference','mash_query','modified_mash_distance','pham_pham_dissimilarity','phylogeny_distance','phage_cluster_compare','phage_subcluster_compare'))








phylogeny_a <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'A')
phylogeny_a1 <- subset(phylogeny_a,phylogeny_a$phage_subcluster_compare == 'A1')
phylogeny_a_not_a1 <- subset(phylogeny_a,phylogeny_a$phage_subcluster_compare != 'A1')

phylogeny_b <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'B')
phylogeny_f <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'F')
phylogeny_k <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'K')
phylogeny_bd <- subset(phylogeny_analysis,phylogeny_analysis$phage_cluster_compare == 'BD')




par(mar=c(4,8,4,4))
plot(phylogeny_f$phylogeny_distance,phylogeny_f$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_b$phylogeny_distance,phylogeny_b$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

par(mar=c(4,8,4,4))
plot(phylogeny_k$phylogeny_distance,phylogeny_k$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_bd$phylogeny_distance,phylogeny_bd$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")


par(mar=c(4,8,4,4))
plot(phylogeny_a_not_a1$phylogeny_distance,phylogeny_a_not_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_a1$phylogeny_distance,phylogeny_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_a_not_a1$phylogeny_distance,phylogeny_a_not_a1$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")






#Compare types, sizes, GC content of shared and unshared genes
par(mar=c(4,8,4,4))
plot(phylogeny_a1$phylogeny_distance,phylogeny_a1$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")

par(mar=c(4,8,4,4))
plot(phylogeny_a_not_a1$phylogeny_distance,phylogeny_a_not_a1$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_a1$phylogeny_distance,phylogeny_a1$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")


par(mar=c(4,8,4,4))
plot(phylogeny_f$phylogeny_distance,phylogeny_f$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")



par(mar=c(4,8,4,4))
plot(phylogeny_a_not_a1$phylogeny_distance,phylogeny_a_not_a1$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_a1$phylogeny_distance,phylogeny_a1$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")



par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_ave_size_unshared_unshared,xlim=c(0,2.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_ave_size_unshared_unshared,xlim=c(0,0.5),ylim=c(0,1200),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")




#Standard plot
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")





#Coding proportion
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_unshared_coding_potential,xlim=c(0,0.5),ylim=c(0,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_unshared_modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.6),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



#total genome sequence size
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0.8E5,1.4E5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0.8E5,1.4E5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_all_total_size,xlim=c(0,0.5),ylim=c(0.8E5,1.4E5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#change in genome size
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,0.2),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,0.2),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$size_diff_ave_percent,xlim=c(0,0.5),ylim=c(0,0.2),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")




#GC content
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_ave_gc_all_all,xlim=c(0,0.5),ylim=c(0.5,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_ave_gc_all_all,xlim=c(0,0.5),ylim=c(0.5,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_ave_gc_all_all,xlim=c(0,0.5),ylim=c(0.5,0.8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Change in GC content
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_modified_gc_diff_all_all,xlim=c(0,0.5),ylim=c(0,0.05),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_modified_gc_diff_all_all,xlim=c(0,0.5),ylim=c(0,0.05),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_modified_gc_diff_all_all,xlim=c(0,0.5),ylim=c(0,0.05),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Change in unshared gene GC content
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_modified_gc_diff_unshared_unshared,xlim=c(0,0.5),ylim=c(0,0.2),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_modified_gc_diff_unshared_unshared,xlim=c(0,0.5),ylim=c(0,0.2),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_modified_gc_diff_unshared_unshared,xlim=c(0,0.5),ylim=c(0,0.2),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#Change in shared gene GC content

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_modified_gc_diff_unshared_unshared,xlim=c(0,0.5),ylim=c(0,0.15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$gsm_modified_gc_diff_shared_shared,xlim=c(0,0.5),ylim=c(0,0.15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")



par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_modified_gc_diff_unshared_unshared,xlim=c(0,0.5),ylim=c(0,0.15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$gsm_modified_gc_diff_shared_shared,xlim=c(0,0.5),ylim=c(0,0.15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")


par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_modified_gc_diff_unshared_unshared,xlim=c(0,0.5),ylim=c(0,0.15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$gsm_modified_gc_diff_shared_shared,xlim=c(0,0.5),ylim=c(0,0.15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")








#Shared pham distribution mean
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_shared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_shared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_shared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_shared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_shared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_shared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")






#Shared pham distribution median
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_shared_pham_distribution_median,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_shared_pham_distribution_median,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_shared_pham_distribution_median,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#Shared pham distribution max
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_shared_pham_distribution_max,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_shared_pham_distribution_max,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_shared_pham_distribution_max,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")





#Unshared pham distribution mean
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_unshared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,10),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_unshared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,10),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_unshared_pham_distribution_mean,xlim=c(0,0.5),ylim=c(0,10),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#Unshared pham distribution median
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_unshared_pham_distribution_median,xlim=c(0,0.5),ylim=c(0,10),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_unshared_pham_distribution_median,xlim=c(0,0.5),ylim=c(0,10),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_unshared_pham_distribution_median,xlim=c(0,0.5),ylim=c(0,10),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Unshared pham distribution max
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_unshared_pham_distribution_max,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_unshared_pham_distribution_max,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_unshared_pham_distribution_max,xlim=c(0,0.5),ylim=c(0,15),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Unshared orpham count
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham2_unshared_orpham_count,xlim=c(0,0.5),ylim=c(0,30),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham2_unshared_orpham_count,xlim=c(0,0.5),ylim=c(0,30),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")

par(mar=c(4,8,4,4))
plot(phylogeny_lytic$phylogeny_distance,phylogeny_lytic$pham2_unshared_orpham_count,xlim=c(0,0.5),ylim=c(0,30),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")




###REVIEW PROGRESS



###CURRENT
#sliding window average
library(caTools)

phylogeny_hgcf_phylosort <- phylogeny_hgcf[order(phylogeny_hgcf$phylogeny_distance),]
phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_median_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_median,101)
phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_max_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_max,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_median_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_median,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_max_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_max,101)
phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_hgcf_phylosort$size_diff_ave_percent,101)
phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_hgcf_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_unshared_coding_potential,101)




phylogeny_lgcf_phylosort <- phylogeny_lgcf[order(phylogeny_lgcf$phylogeny_distance),]
phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_median_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_median,101)
phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_max_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_max,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_median_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_median,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_max_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_max,101)
phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lgcf_phylosort$size_diff_ave_percent,101)
phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lgcf_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_unshared_coding_potential,101)



phylogeny_lytic_phylosort <- phylogeny_lytic[order(phylogeny_lytic$phylogeny_distance),]
phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_shared_pham_distribution_median_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_median,101)
phylogeny_lytic_phylosort$pham2_shared_pham_distribution_max_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_max,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_median_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_median,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_max_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_max,101)
phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lytic_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lytic_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lytic_phylosort$size_diff_ave_percent,101)
phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lytic_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_lytic_phylosort$gsm_unshared_coding_potential,101)


#Shared pham distribution mean
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(0.5,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(0.5,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(0.5,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

#Unshared pham distribution mean
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(0.5,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(0.5,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(0.5,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


#Shared pham distribution median
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

#Unshared pham distribution median
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



#Shared pham distribution max
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

#Unshared pham distribution max
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")





#Unshared orpham count
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),ylim=c(0,8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),ylim=c(0,8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),ylim=c(0,8),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(mar=c(4,8,4,4))
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(mar=c(4,8,4,4))
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")








#Shared/unshared pham distribution mean
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(1,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(1,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(1,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(1,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(1,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.5),ylim=c(1,3.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")





#Shared/unshared pham distribution median
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_median_runmean,xlim=c(0,0.5),ylim=c(0.5,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")




#Shared/unshared pham distribution max
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")
par(new=TRUE)
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_max_runmean,xlim=c(0,0.5),ylim=c(2,11),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


#Unshared orpham count
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.5),ylim=c(0,5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")




#sliding window average
#X axis = c(0,0.3)
#only use data for sliding windows that fit within the scatter plot boundaries
library(caTools)

phylogeny_hgcf_phylosort <- phylogeny_hgcf[order(phylogeny_hgcf$phylogeny_distance),]
phylogeny_hgcf_phylosort <- phylogeny_hgcf_phylosort[phylogeny_hgcf_phylosort$phylogeny_distance < 0.3,]

phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_median_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_median,101)
phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_max_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_max,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_median_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_median,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_max_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_max,101)
phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_hgcf_phylosort$size_diff_ave_percent,101)
phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_hgcf_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_unshared_coding_potential,101)
phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared,101)



phylogeny_lgcf_phylosort <- phylogeny_lgcf[order(phylogeny_lgcf$phylogeny_distance),]
phylogeny_lgcf_phylosort <- phylogeny_lgcf_phylosort[phylogeny_lgcf_phylosort$phylogeny_distance < 0.3,]

phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_median_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_median,101)
phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_max_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_max,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_median_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_median,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_max_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_max,101)
phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lgcf_phylosort$size_diff_ave_percent,101)
phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lgcf_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_unshared_coding_potential,101)
phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared,101)


phylogeny_lytic_phylosort <- phylogeny_lytic[order(phylogeny_lytic$phylogeny_distance),]
phylogeny_lytic_phylosort <- phylogeny_lytic_phylosort[phylogeny_lytic_phylosort$phylogeny_distance < 0.3,]

phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_shared_pham_distribution_median_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_median,101)
phylogeny_lytic_phylosort$pham2_shared_pham_distribution_max_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_max,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_median_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_median,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_max_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_max,101)
phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lytic_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lytic_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lytic_phylosort$size_diff_ave_percent,101)
phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lytic_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_lytic_phylosort$gsm_unshared_coding_potential,101)
phylogeny_lytic_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_shared_shared,101)






#Gene content dissimilarity
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Ave genome size difference
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")






###NEED TO CHANGE COLORS
#Unshared coding sequence proportion
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


#Shared/unshared ave gene size
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")




#Shared/unshared pham distribution mean
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")




#Unshared orpham count
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
















###Function-specific pham analysis
#Actino782 analysis
#similar code to the phylogeny section
#Designate pham data as set #2, since it is for a subset of all comparisons (actino782), and contains additional columns 
#One version of the pham analysis output file now contains two rows at the top of the file containing data on the when the analysis was run. These need to first be removed in bash ($ tail -n +3 'file.name' > 'new_file.name'), so the file is now labeled as 'mod'


actino_pham_data <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170315_bac2333_actino_only_pairwise_pham_proportions.csv",sep=",",header=TRUE)


#Old renaming. Invalid now that the script output has changed.
# names(actino_pham_data) <- c("pham2_phage1","pham2_phage1_number_unshared_phams","pham2_phage1_shared_proportion","pham2_phage2",
#                              "pham2_phage2_number_unshared_phams","pham2_phage2_shared_proportion","pham2_number_shared_phams","pham2_average_shared_proportion",
#                              "pham2_jaccard_similarity","pham2_shared_pham_distribution_mean","pham2_shared_pham_distribution_median","pham2_shared_pham_distribution_max",
#                              "pham2_unshared_pham_distribution_mean","pham2_unshared_pham_distribution_median",
#                              "pham2_unshared_pham_distribution_max","pham2_unshared_orpham_count",
#                              "pham2_shared_tally_Unspecified","pham2_unshared_tally_Unspecified","pham2_unshared_percentage_Unspecified",
#                              "pham2_shared_tally_defense","pham2_unshared_tally_defense","pham2_unshared_percentage_defense",
#                              "pham2_shared_tally_dna_metabolism","pham2_unshared_tally_dna_metabolism","pham2_unshared_percentage_dna_metabolism",
#                              "pham2_shared_tally_lysis","pham2_unshared_tally_lysis","pham2_unshared_percentage_lysis",
#                              "pham2_shared_tally_lysogeny","pham2_unshared_tally_lysogeny","pham2_unshared_percentage_lysogeny",
#                              "pham2_shared_tally_mobile","pham2_unshared_tally_mobile","pham2_unshared_percentage_mobile",
#                              "pham2_shared_tally_other","pham2_unshared_tally_other","pham2_unshared_percentage_other",
#                              "pham2_shared_tally_recombination_replication","pham2_unshared_tally_recombination_replication","pham2_unshared_percentage_recombination_replication",
#                              "pham2_shared_tally_structure_assembly","pham2_unshared_tally_structure_assembly","pham2_unshared_percentage_structure_assembly")
# #The percentage data is not exactly accurate. If there were no shared or unshared phams in the function category, technically you can't compute a percentage
# #(since dividing by zero), so a Zero was inserted in these cases. But really, they should be NA.
# actino_pham_data$pham2_unshared_percentage_Unspecified_mod <- ifelse(actino_pham_data$pham2_shared_tally_Unspecified + actino_pham_data$pham2_unshared_tally_Unspecified > 0,actino_pham_data$pham2_unshared_percentage_Unspecified,NA)
# actino_pham_data$pham2_unshared_percentage_defense_mod <- ifelse(actino_pham_data$pham2_shared_tally_defense + actino_pham_data$pham2_unshared_tally_defense > 0,actino_pham_data$pham2_unshared_percentage_defense,NA)
# actino_pham_data$pham2_unshared_percentage_dna_metabolism_mod <- ifelse(actino_pham_data$pham2_shared_tally_dna_metabolism + actino_pham_data$pham2_unshared_tally_dna_metabolism > 0,actino_pham_data$pham2_unshared_percentage_dna_metabolism,NA)
# actino_pham_data$pham2_unshared_percentage_lysis_mod <- ifelse(actino_pham_data$pham2_shared_tally_lysis + actino_pham_data$pham2_unshared_tally_lysis > 0,actino_pham_data$pham2_unshared_percentage_lysis,NA)
# actino_pham_data$pham2_unshared_percentage_lysogeny_mod <- ifelse(actino_pham_data$pham2_shared_tally_lysogeny + actino_pham_data$pham2_unshared_tally_lysogeny > 0,actino_pham_data$pham2_unshared_percentage_lysogeny,NA)
# actino_pham_data$pham2_unshared_percentage_mobile_mod <- ifelse(actino_pham_data$pham2_shared_tally_mobile + actino_pham_data$pham2_unshared_tally_mobile > 0,actino_pham_data$pham2_unshared_percentage_mobile,NA)
# actino_pham_data$pham2_unshared_percentage_other_mod <- ifelse(actino_pham_data$pham2_shared_tally_other + actino_pham_data$pham2_unshared_tally_other > 0,actino_pham_data$pham2_unshared_percentage_other,NA)
# actino_pham_data$pham2_unshared_percentage_recombination_replication_mod <- ifelse(actino_pham_data$pham2_shared_tally_recombination_replication + actino_pham_data$pham2_unshared_tally_recombination_replication > 0,actino_pham_data$pham2_unshared_percentage_recombination_replication,NA)
# actino_pham_data$pham2_unshared_percentage_structure_assembly_mod <- ifelse(actino_pham_data$pham2_shared_tally_structure_assembly + actino_pham_data$pham2_unshared_tally_structure_assembly > 0,actino_pham_data$pham2_unshared_percentage_structure_assembly,NA)




names(actino_pham_data) <- c("pham2_phage1","pham2_phage1_number_unshared_phams","pham2_phage1_shared_proportion","pham2_phage2",
                             "pham2_phage2_number_unshared_phams","pham2_phage2_shared_proportion","pham2_number_shared_phams","pham2_average_shared_proportion",
                             "pham2_jaccard_similarity","pham2_shared_pham_distribution_mean","pham2_shared_pham_distribution_median","pham2_shared_pham_distribution_max",
                             "pham2_unshared_pham_distribution_mean","pham2_unshared_pham_distribution_median",
                             "pham2_unshared_pham_distribution_max","pham2_unshared_orpham_count",
                             "pham2_Unspecified_phage1_number_of_unshared_phams","pham2_Unspecified_phage2_number_of_unshared_phams","pham2_Unspecified_number_of_shared_phams","pham2_Unspecified_average_shared_proportion",
                             "pham2_defense_phage1_number_of_unshared_phams","pham2_defense_phage2_number_of_unshared_phams","pham2_defense_number_of_shared_phams","pham2_defense_average_shared_proportion",
                             "pham2_dna_metabolism_phage1_number_of_unshared_phams","pham2_dna_metabolism_phage2_number_of_unshared_phams","pham2_dna_metabolism_number_of_shared_phams","pham2_dna_metabolism_average_shared_proportion",
                             "pham2_lysis_phage1_number_of_unshared_phams","pham2_lysis_phage2_number_of_unshared_phams","pham2_lysis_number_of_shared_phams","pham2_lysis_average_shared_proportion",
                             "pham2_lysogeny_phage1_number_of_unshared_phams","pham2_lysogeny_phage2_number_of_unshared_phams","pham2_lysogeny_number_of_shared_phams","pham2_lysogeny_average_shared_proportion",
                             "pham2_mobile_phage1_number_of_unshared_phams","pham2_mobile_phage2_number_of_unshared_phams","pham2_mobile_number_of_shared_phams","pham2_mobile_average_shared_proportion",
                             "pham2_other_phage1_number_of_unshared_phams","pham2_other_phage2_number_of_unshared_phams","pham2_other_number_of_shared_phams","pham2_other_average_shared_proportion",
                             "pham2_recombination_replication_phage1_number_of_unshared_phams","pham2_recombination_replication_phage2_number_of_unshared_phams","pham2_recombination_replication_number_of_shared_phams","pham2_recombination_replication_average_shared_proportion",
                             "pham2_structure_assembly_phage1_number_of_unshared_phams","pham2_structure_assembly_phage2_number_of_unshared_phams","pham2_structure_assembly_number_of_shared_phams","pham2_structure_assembly_average_shared_proportion")



#The function-specific average shared proportion data is not exactly accurate. If there were no shared or unshared phams in the function category, technically you can't compute a shared proportion
#(since dividing by zero), so a "-1" was inserted in these cases. But really, they should be NA.
actino_pham_data$pham2_Unspecified_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_Unspecified_average_shared_proportion != -1,actino_pham_data$pham2_Unspecified_average_shared_proportion,NA)
actino_pham_data$pham2_defense_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_defense_average_shared_proportion != -1,actino_pham_data$pham2_defense_average_shared_proportion,NA)
actino_pham_data$pham2_dna_metabolism_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_dna_metabolism_average_shared_proportion != -1,actino_pham_data$pham2_dna_metabolism_average_shared_proportion,NA)
actino_pham_data$pham2_lysis_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_lysis_average_shared_proportion != -1,actino_pham_data$pham2_lysis_average_shared_proportion,NA)
actino_pham_data$pham2_lysogeny_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_lysogeny_average_shared_proportion != -1,actino_pham_data$pham2_lysogeny_average_shared_proportion,NA)
actino_pham_data$pham2_mobile_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_mobile_average_shared_proportion != -1,actino_pham_data$pham2_mobile_average_shared_proportion,NA)
actino_pham_data$pham2_other_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_other_average_shared_proportion != -1,actino_pham_data$pham2_other_average_shared_proportion,NA)
actino_pham_data$pham2_recombination_replication_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_recombination_replication_average_shared_proportion != -1,actino_pham_data$pham2_recombination_replication_average_shared_proportion,NA)
actino_pham_data$pham2_structure_assembly_average_shared_proportion_mod <- ifelse(actino_pham_data$pham2_structure_assembly_average_shared_proportion != -1,actino_pham_data$pham2_structure_assembly_average_shared_proportion,NA)




#Compute gene content dissimilarity
actino_pham_data$pham2_pham_dissimilarity <- 1 - actino_pham_data$pham2_average_shared_proportion
actino_pham_data$pham2_jaccard_dissimilarity <- 1 - actino_pham_data$pham2_jaccard_similarity
actino_pham_data$pham2_pham_dissimilarity_Unspecified <- 1 - actino_pham_data$pham2_Unspecified_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_defense <- 1 - actino_pham_data$pham2_defense_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_dna_metabolism <- 1 - actino_pham_data$pham2_dna_metabolism_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_lysis <- 1 - actino_pham_data$pham2_lysis_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_lysogeny <- 1 - actino_pham_data$pham2_lysogeny_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_mobile <- 1 - actino_pham_data$pham2_mobile_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_other <- 1 - actino_pham_data$pham2_other_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_recomb_rep <- 1 - actino_pham_data$pham2_recombination_replication_average_shared_proportion_mod
actino_pham_data$pham2_pham_dissimilarity_structural <- 1 - actino_pham_data$pham2_structure_assembly_average_shared_proportion_mod








#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
actino_pham_data$pham2_phage1_phage2 <- paste(actino_pham_data$pham2_phage1,"_",actino_pham_data$pham2_phage2,sep="")
actino_pham_data$pham2_phage1_phage2 <- as.factor(actino_pham_data$pham2_phage1_phage2)






#To retain all rows, be sure to keep all.x=TRUE, but you don't want to retain all rows = making scatter plots or histograms can cause errors if not all rows have data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
actino_pham_analysis <- merge(mash_table2,actino_pham_data,by.x="mash_ref_query",by.y="pham2_phage1_phage2")
#OR merge with comparisons from gene-specific mash analysis
actino_pham_analysis <- merge(gene_specific_mash_table,actino_pham_data,by.x="mash_ref_query",by.y="pham2_phage1_phage2")
#OR intra-cluster comparisons
actino_intra_cluster_comparisons <- subset(mash_table2,mash_table2$phage_cluster_source_compare == "actino" &
                                      mash_table2$phage_cluster_compare != "different")
actino_pham_analysis <- merge(actino_intra_cluster_comparisons,actino_pham_data,by.x="mash_ref_query",by.y="pham2_phage1_phage2")


#Check for correct data before proceeding
par(mar=c(4,8,4,4))
plot(actino_pham_analysis$modified_mash_distance,actino_pham_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)




#Define modes
actino_pham_analysis_hgcf <- subset(actino_pham_analysis,actino_pham_analysis$gene_flux_category == 'high' & actino_pham_analysis$phage_temperate_compare == 'yes')
actino_pham_analysis_lgcf <- subset(actino_pham_analysis,actino_pham_analysis$gene_flux_category == 'low' & actino_pham_analysis$phage_temperate_compare == 'yes')
actino_pham_analysis_lytic <- subset(actino_pham_analysis,actino_pham_analysis$phage_temperate_compare == 'no')


par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_recomb_rep,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_recomb_rep,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_recomb_rep,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_structural,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_structural,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_structural,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_Unspecified,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_Unspecified,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_Unspecified,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_defense,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_defense,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_defense,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_dna_metabolism,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_dna_metabolism,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_dna_metabolism,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_lysis,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_lysis,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_lysis,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_lysogeny,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_lysogeny,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_lysogeny,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_mobile,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_mobile,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_mobile,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf$modified_mash_distance,actino_pham_analysis_hgcf$pham2_pham_dissimilarity_other,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf$modified_mash_distance,actino_pham_analysis_lgcf$pham2_pham_dissimilarity_other,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic$modified_mash_distance,actino_pham_analysis_lytic$pham2_pham_dissimilarity_other,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col="red")
abline(0,2,lty=2,lwd=3,col="grey")




#Cluster-specific, function-specific plots
plot_function_specific_data <- function(cluster_specific,color){
  

  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_recomb_rep,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_structural,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_Unspecified,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_defense,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_dna_metabolism,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_lysis,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_lysogeny,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")

  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_mobile,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")

  par(mar=c(4,8,4,4))
  plot(cluster_specific$modified_mash_distance,cluster_specific$pham2_pham_dissimilarity_other,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,las=1,main=NULL,col=color)
  abline(0,2,lty=2,lwd=3,col="grey")
  
}


#hgcf clusters
dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'F')
plot_function_specific_data(cluster_specific,'blue')

dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_subcluster_compare == 'A1')
plot_function_specific_data(cluster_specific,'blue')


dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'I')
plot_function_specific_data(cluster_specific,'blue')


dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'CZ')
plot_function_specific_data(cluster_specific,'blue')


dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'CV')
plot_function_specific_data(cluster_specific,'blue')

dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'BW')
plot_function_specific_data(cluster_specific,'blue')



#lgcf
dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'K')
plot_function_specific_data(cluster_specific,'green')


#lytic
dev.off()
cluster_specific <- subset(actino_pham_analysis,actino_pham_analysis$phage_cluster_compare == 'B')
plot_function_specific_data(cluster_specific,'red')






#sliding window average
#Note, there are NA values in many of these vectors. The default algorithm for runmean is "C", which is able to handle NAs.
#However, if a different algorithm needs to be used, it may/may not be able to handle the NA values.
library(caTools)

actino_pham_analysis_hgcf_mmdsort <- actino_pham_analysis_hgcf[order(actino_pham_analysis_hgcf$modified_mash_distance),]
actino_pham_analysis_hgcf_mmdsort$pham2_shared_pham_distribution_mean_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_shared_pham_distribution_mean,101)
actino_pham_analysis_hgcf_mmdsort$pham2_shared_pham_distribution_median_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_shared_pham_distribution_median,101)
actino_pham_analysis_hgcf_mmdsort$pham2_shared_pham_distribution_max_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_shared_pham_distribution_max,101)
actino_pham_analysis_hgcf_mmdsort$pham2_unshared_pham_distribution_mean_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_unshared_pham_distribution_mean,101)
actino_pham_analysis_hgcf_mmdsort$pham2_unshared_pham_distribution_median_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_unshared_pham_distribution_median,101)
actino_pham_analysis_hgcf_mmdsort$pham2_unshared_pham_distribution_max_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_unshared_pham_distribution_max,101)
actino_pham_analysis_hgcf_mmdsort$pham2_unshared_orpham_count_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_unshared_orpham_count,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_Unspecified,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_defense_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_defense,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_dna_metabolism_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_dna_metabolism,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_lysis_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_lysis,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_lysogeny_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_lysogeny,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_mobile_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_mobile,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_other_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_other,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_recomb_rep_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_recomb_rep,101)
actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_structural_runmean <- runmean(actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_structural,101)




actino_pham_analysis_lgcf_mmdsort <- actino_pham_analysis_lgcf[order(actino_pham_analysis_lgcf$modified_mash_distance),]
actino_pham_analysis_lgcf_mmdsort$pham2_shared_pham_distribution_mean_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_shared_pham_distribution_mean,101)
actino_pham_analysis_lgcf_mmdsort$pham2_shared_pham_distribution_median_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_shared_pham_distribution_median,101)
actino_pham_analysis_lgcf_mmdsort$pham2_shared_pham_distribution_max_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_shared_pham_distribution_max,101)
actino_pham_analysis_lgcf_mmdsort$pham2_unshared_pham_distribution_mean_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_unshared_pham_distribution_mean,101)
actino_pham_analysis_lgcf_mmdsort$pham2_unshared_pham_distribution_median_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_unshared_pham_distribution_median,101)
actino_pham_analysis_lgcf_mmdsort$pham2_unshared_pham_distribution_max_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_unshared_pham_distribution_max,101)
actino_pham_analysis_lgcf_mmdsort$pham2_unshared_orpham_count_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_unshared_orpham_count,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_Unspecified,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_defense_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_defense,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_dna_metabolism_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_dna_metabolism,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_lysis_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_lysis,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_lysogeny_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_lysogeny,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_mobile_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_mobile,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_other_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_other,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_recomb_rep_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_recomb_rep,101)
actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_structural_runmean <- runmean(actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_structural,101)





actino_pham_analysis_lytic_mmdsort <- actino_pham_analysis_lytic[order(actino_pham_analysis_lytic$modified_mash_distance),]
actino_pham_analysis_lytic_mmdsort$pham2_shared_pham_distribution_mean_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_shared_pham_distribution_mean,101)
actino_pham_analysis_lytic_mmdsort$pham2_shared_pham_distribution_median_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_shared_pham_distribution_median,101)
actino_pham_analysis_lytic_mmdsort$pham2_shared_pham_distribution_max_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_shared_pham_distribution_max,101)
actino_pham_analysis_lytic_mmdsort$pham2_unshared_pham_distribution_mean_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_unshared_pham_distribution_mean,101)
actino_pham_analysis_lytic_mmdsort$pham2_unshared_pham_distribution_median_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_unshared_pham_distribution_median,101)
actino_pham_analysis_lytic_mmdsort$pham2_unshared_pham_distribution_max_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_unshared_pham_distribution_max,101)
actino_pham_analysis_lytic_mmdsort$pham2_unshared_orpham_count_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_unshared_orpham_count,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_Unspecified,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_defense_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_defense,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_dna_metabolism_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_dna_metabolism,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_lysis_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_lysis,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_lysogeny_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_lysogeny,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_mobile_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_mobile,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_other_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_other,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_recomb_rep_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_recomb_rep,101)
actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_structural_runmean <- runmean(actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_structural,101)














# #List of runmean variables to plot
# pham2_shared_pham_distribution_mean_runmean
# pham2_shared_pham_distribution_median_runmean
# pham2_shared_pham_distribution_max_runmean
# pham2_unshared_pham_distribution_mean_runmean
# pham2_unshared_pham_distribution_median_runmean
# pham2_unshared_pham_distribution_max_runmean
# pham2_unshared_orpham_count_runmean
# pham2_pham_dissimilarity_runmean
# pham2_pham_dissimilarity_Unspecified_runmean
# pham2_pham_dissimilarity_defense_runmean
# pham2_pham_dissimilarity_dna_metabolism_runmean
# pham2_pham_dissimilarity_lysis_runmean
# pham2_pham_dissimilarity_lysogeny_runmean
# pham2_pham_dissimilarity_mobile_runmean
# pham2_pham_dissimilarity_other_runmean
# pham2_pham_dissimilarity_recomb_rep_runmean
# pham2_pham_dissimilarity_structural_runmean














#Plot all runmean data


#Standard plot
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#function = ALL
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")




#function = structure/assembly
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_structural_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_structural_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_structural_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#function = defense
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_defense_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_defense_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_defense_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#function = dna metabolism
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_dna_metabolism_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_dna_metabolism_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_dna_metabolism_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#function = lysis
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_lysis_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_lysis_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_lysis_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#function = lysogeny
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_lysogeny_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_lysogeny_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_lysogeny_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#function = mobile
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_mobile_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_mobile_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_mobile_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")





#function = recombination/replication
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_recomb_rep_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_recomb_rep_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_recomb_rep_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")




#function = other
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_other_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_other_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_other_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#function = unspecified
par(mar=c(4,8,4,4))
plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(actino_pham_analysis_lgcf_mmdsort$modified_mash_distance,actino_pham_analysis_lgcf_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(actino_pham_analysis_lytic_mmdsort$modified_mash_distance,actino_pham_analysis_lytic_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")













###Count analysis - pham gains/losses
actino_pham_function_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170309_bac2333_actino_only_pham_function_mod.csv",sep=",",header=TRUE)
actino_pham_function_table$pham <- factor(actino_pham_function_table$pham)

names(actino_pham_function_table) <- c("pham","count_pham_function")


#count_cluster_f_pham_gains_losses <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170321_count_cluster_f_pham_gains_losses.csv",sep=",",header=TRUE)
# names(count_cluster_f_pham_gains_losses) <- c("count_cluster","count_pham","count_gains","count_losses")
# count_cluster_f_pham_gains_losses$count_pham <- factor(count_cluster_f_pham_gains_losses$count_pham)
# count_cluster_f_pham_gains_losses <- merge(count_cluster_f_pham_gains_losses,actino_pham_function_table,by.x="count_pham",by.y="pham",all.x=TRUE)


count_pham_gains_losses <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170321_count_pham_gains_losses.csv",sep=",",header=TRUE)

names(count_pham_gains_losses) <- c("count_cluster","count_pham","count_gains","count_losses")
count_pham_gains_losses$count_pham <- factor(count_pham_gains_losses$count_pham)

count_pham_gains_losses <- merge(count_pham_gains_losses,actino_pham_function_table,by.x="count_pham",by.y="pham",all.x=TRUE)

# par(mar=c(4,8,4,4))
# plot(actino_pham_analysis_hgcf_mmdsort$modified_mash_distance,actino_pham_analysis_hgcf_mmdsort$pham2_pham_dissimilarity_Unspecified_runmean,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")


par(mar=c(4,12,4,4))
boxplot(count_pham_gains_losses$count_gains ~ count_pham_gains_losses$count_cluster * count_pham_gains_losses$count_pham_function,las=1,horizontal=TRUE)

par(mar=c(4,12,4,4))
boxplot(count_pham_gains_losses$count_gains ~ count_pham_gains_losses$count_cluster,las=1,horizontal=TRUE)

par(mar=c(4,12,4,4))
boxplot(count_pham_gains_losses$count_losses ~ count_pham_gains_losses$count_cluster,ylim=c(0,5),las=1,horizontal=TRUE)























###Predict evolutionary mode

data_for_mode_prediction <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria" &
                                     mash_table2$phage_viral_type_compare == "dsDNA" &
                                     mash_table2$modified_mash_distance < 0.42 &
                                     mash_table2$pham_pham_dissimilarity < 0.89)
data_for_mode_prediction <- subset(data_for_mode_prediction,select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))
write.table(data_for_mode_prediction,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/20170315_data_for_mode_prediction.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)



#Run the data through the analyze_mash_network_script to predict the evolutionary mode

#Then import the mode prediction back into R
mode_prediction_table <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170315_mode_prediction.csv",sep=",",header=TRUE)

names(mode_prediction_table) <- c("phage_identifier",
                                  "mode_prediction_hgcf_tally","mode_prediction_lgcf_tally","mode_prediction_out_of_range_tally",
                                  "mode_prediction_hgcf_percent","mode_prediction_lgcf_percent","mode_prediction_mode")



host_mode_table <- merge(host_table,mode_prediction_table,by.x="phage_identifier",by.y="phage_identifier",all.x=TRUE)


par(mar=c(4,8,4,4))
hist(host_mode_table$mode_prediction_hgcf_percent,col="black",breaks=25,cex.axis=2,ann=FALSE,las=1,ylim=c(0,1400))



#Mode prediction script assigns evolutionary mode to each phage if > 50% comparisons are in the HGCF or LGCF mode
#Now, assign mode based on "exact" (100% of comparisons are in either mode) or "approximate" (> 80% of comparisons are in either mode)
host_mode_table$mode_prediction_mode_exact <- ifelse(is.na(host_mode_table$mode_prediction_mode) == FALSE & (host_mode_table$mode_prediction_hgcf_percent == 1 | host_mode_table$mode_prediction_hgcf_percent == 0),as.character(host_mode_table$mode_prediction_mode),'unknown')
host_mode_table$mode_prediction_mode_exact <- factor(host_mode_table$mode_prediction_mode_exact)

host_mode_table$mode_prediction_mode_approx <- ifelse(is.na(host_mode_table$mode_prediction_mode) == FALSE & (host_mode_table$mode_prediction_hgcf_percent > 0.8 | host_mode_table$mode_prediction_hgcf_percent < 0.2),as.character(host_mode_table$mode_prediction_mode),'unknown')
host_mode_table$mode_prediction_mode_approx <- factor(host_mode_table$mode_prediction_mode_approx)


#Export table for analysis in Excel
write.table(host_mode_table,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/20170315_host_mode_table.csv",sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)

#Convert all 'NA' values to 'unknown'




#Analysis based on known lifestyle
actual_temperate_mode <- subset(host_mode_table,host_mode_table$phage_temperate == "yes")
actual_lytic_mode <- subset(host_mode_table,host_mode_table$phage_temperate == "no")

actual_temperate_mode_exact_actino <- subset(actual_temperate_mode,actual_temperate_mode$host_phylum == 'Actinobacteria')
actual_temperate_mode_exact_cyano <- subset(actual_temperate_mode,actual_temperate_mode$host_phylum == 'Cyanobacteria')
actual_temperate_mode_exact_firm <- subset(actual_temperate_mode,actual_temperate_mode$host_phylum == 'Firmicutes')
actual_temperate_mode_exact_proteo <- subset(actual_temperate_mode,actual_temperate_mode$host_phylum == 'Proteobacteria')

actual_lytic_mode_exact_actino <- subset(actual_lytic_mode,actual_lytic_mode$host_phylum == 'Actinobacteria')
actual_lytic_mode_exact_cyano <- subset(actual_lytic_mode,actual_lytic_mode$host_phylum == 'Cyanobacteria')
actual_lytic_mode_exact_firm <- subset(actual_lytic_mode,actual_lytic_mode$host_phylum == 'Firmicutes')
actual_lytic_mode_exact_proteo <- subset(actual_lytic_mode,actual_lytic_mode$host_phylum == 'Proteobacteria')



pie(summary(actual_temperate_mode_exact_actino$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(actual_temperate_mode_exact_cyano$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(actual_temperate_mode_exact_firm$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(actual_temperate_mode_exact_proteo$mode_prediction_mode_exact),col=c("blue","green","white"))


pie(summary(actual_lytic_mode_exact_actino$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(actual_lytic_mode_exact_cyano$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(actual_lytic_mode_exact_firm$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(actual_lytic_mode_exact_proteo$mode_prediction_mode_exact),col=c("blue","green","white"))






#Analysis based on predicted lifestyle
predicted_temperate_mode <- subset(host_mode_table,host_mode_table$predicted_temperate == "yes")
predicted_lytic_mode <- subset(host_mode_table,host_mode_table$predicted_temperate == "no")




hist(predicted_temperate_mode$mode_prediction_hgcf_percent,breaks=100,ylim=c(0,20))
hist(predicted_lytic_mode$mode_prediction_hgcf_percent,breaks=100,ylim=c(0,20))


# predicted_temperate_mode_exact <- subset(predicted_temperate_mode,(predicted_temperate_mode$mode_prediction_hgcf_percent == 1 | predicted_temperate_mode$mode_prediction_hgcf_percent == 0) & predicted_temperate_mode$mode_prediction_mode != "unknown")
# predicted_lytic_mode_exact <- subset(predicted_lytic_mode,(predicted_lytic_mode$mode_prediction_hgcf_percent == 1 | predicted_lytic_mode$mode_prediction_hgcf_percent == 0) & predicted_lytic_mode$mode_prediction_mode != "unknown")
# predicted_temperate_mode_exact$mode_prediction_mode <- factor(predicted_temperate_mode_exact$mode_prediction_mode)
# predicted_lytic_mode_exact$mode_prediction_mode <- factor(predicted_lytic_mode_exact$mode_prediction_mode)
# predicted_temperate_mode_exact_actino <- subset(predicted_temperate_mode_exact,predicted_temperate_mode_exact$host_phylum == 'Actinobacteria')
# predicted_temperate_mode_exact_cyano <- subset(predicted_temperate_mode_exact,predicted_temperate_mode_exact$host_phylum == 'Cyanobacteria')
# predicted_temperate_mode_exact_firm <- subset(predicted_temperate_mode_exact,predicted_temperate_mode_exact$host_phylum == 'Firmicutes')
# predicted_temperate_mode_exact_proteo <- subset(predicted_temperate_mode_exact,predicted_temperate_mode_exact$host_phylum == 'Proteobacteria')
# 
# predicted_lytic_mode_exact_actino <- subset(predicted_lytic_mode_exact,predicted_lytic_mode_exact$host_phylum == 'Actinobacteria')
# predicted_lytic_mode_exact_cyano <- subset(predicted_lytic_mode_exact,predicted_lytic_mode_exact$host_phylum == 'Cyanobacteria')
# predicted_lytic_mode_exact_firm <- subset(predicted_lytic_mode_exact,predicted_lytic_mode_exact$host_phylum == 'Firmicutes')
# predicted_lytic_mode_exact_proteo <- subset(predicted_lytic_mode_exact,predicted_lytic_mode_exact$host_phylum == 'Proteobacteria')

predicted_temperate_mode_exact_actino <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Actinobacteria')
predicted_temperate_mode_exact_cyano <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Cyanobacteria')
predicted_temperate_mode_exact_firm <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Firmicutes')
predicted_temperate_mode_exact_proteo <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Proteobacteria')

predicted_lytic_mode_exact_actino <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Actinobacteria')
predicted_lytic_mode_exact_cyano <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Cyanobacteria')
predicted_lytic_mode_exact_firm <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Firmicutes')
predicted_lytic_mode_exact_proteo <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Proteobacteria')


pie(summary(predicted_temperate_mode_exact_actino$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_temperate_mode_exact_cyano$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_temperate_mode_exact_firm$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_temperate_mode_exact_proteo$mode_prediction_mode_exact),col=c("blue","green","white"))


pie(summary(predicted_lytic_mode_exact_actino$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_lytic_mode_exact_cyano$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_lytic_mode_exact_firm$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_lytic_mode_exact_proteo$mode_prediction_mode_exact),col=c("blue","green","white"))





predicted_temperate_actino <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Actinobacteria')
predicted_temperate_bacter <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Bacteroidetes')
predicted_temperate_cyano <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Cyanobacteria')
predicted_temperate_firm <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Firmicutes')
predicted_temperate_proteo <- subset(predicted_temperate_mode,predicted_temperate_mode$host_phylum == 'Proteobacteria')



predicted_lytic_actino <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Actinobacteria')
predicted_lytic_bacter <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Bacteroidetes')
predicted_lytic_cyano <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Cyanobacteria')
predicted_lytic_firm <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Firmicutes')
predicted_lytic_proteo <- subset(predicted_lytic_mode,predicted_lytic_mode$host_phylum == 'Proteobacteria')


pie(summary(predicted_temperate_actino$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_temperate_cyano$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_temperate_firm$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_temperate_proteo$mode_prediction_mode_exact),col=c("blue","green","white"))


pie(summary(predicted_lytic_mode_exact_actino$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_lytic_mode_exact_cyano$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_lytic_mode_exact_firm$mode_prediction_mode_exact),col=c("blue","green","white"))
pie(summary(predicted_lytic_mode_exact_proteo$mode_prediction_mode_exact),col=c("blue","green","white"))







#Match the predicted mode to the main comparison data, to be able to subset certain types of phages


host_mode_table_reduced <- subset(host_mode_table,select= c("phage_identifier","mode_prediction_mode","mode_prediction_mode_exact","mode_prediction_mode_approx"))



names(host_mode_table_reduced) <- c("phage_identifier","query_mode_prediction_mode","query_mode_prediction_mode_exact","query_mode_prediction_mode_approx")


predicted_mode_analysis <- merge(mash_table2,host_mode_table_reduced,by.x="mash_query",by.y="phage_identifier")
names(host_mode_table_reduced) <- c("phage_identifier","ref_mode_prediction_mode","ref_mode_prediction_mode_exact","ref_mode_prediction_mode_approx")
predicted_mode_analysis <- merge(predicted_mode_analysis,host_mode_table_reduced,by.x="mash_reference",by.y="phage_identifier")



predicted_mode_analysis$mode_approx_compare <- ifelse(predicted_mode_analysis$ref_mode_prediction_mode_approx == predicted_mode_analysis$query_mode_prediction_mode_approx,as.character(predicted_mode_analysis$ref_mode_prediction_mode_approx),"different")
predicted_mode_analysis$mode_approx_compare <- factor(predicted_mode_analysis$mode_approx_compare)


predicted_mode_analysis_lytic_hgcf <- subset(predicted_mode_analysis,(predicted_mode_analysis$ref_predicted_temperate == "no" &
                                               predicted_mode_analysis$ref_mode_prediction_mode_approx == "hgcf") |
                                               (predicted_mode_analysis$query_predicted_temperate == "no" &
                                               predicted_mode_analysis$query_mode_prediction_mode_approx == "hgcf")
                                               )



predicted_mode_analysis_actino_lytic_hgcf <- subset(predicted_mode_analysis,predicted_mode_analysis$ref_host_phylum == "Actinobacteria" &
                                                      predicted_mode_analysis$query_host_phylum == "Actinobacteria" &
                                                      ((predicted_mode_analysis$ref_predicted_temperate == "no" &
                                                      predicted_mode_analysis$ref_mode_prediction_mode_approx == "hgcf") |
                                                      (predicted_mode_analysis$query_predicted_temperate == "no" &
                                                      predicted_mode_analysis$query_mode_prediction_mode_approx == "hgcf")))


predicted_mode_analysis_bacter_lytic_hgcf <- subset(predicted_mode_analysis,predicted_mode_analysis$ref_host_phylum == "Bacteroidetes" &
                                                      predicted_mode_analysis$query_host_phylum == "Bacteroidetes" &
                                                      ((predicted_mode_analysis$ref_predicted_temperate == "no" &
                                                          predicted_mode_analysis$ref_mode_prediction_mode_approx == "hgcf") |
                                                         (predicted_mode_analysis$query_predicted_temperate == "no" &
                                                            predicted_mode_analysis$query_mode_prediction_mode_approx == "hgcf")))

predicted_mode_analysis_cyano_lytic_hgcf <- subset(predicted_mode_analysis,predicted_mode_analysis$ref_host_phylum == "Cyanobacteria" &
                                                      predicted_mode_analysis$query_host_phylum == "Cyanobacteria" &
                                                      ((predicted_mode_analysis$ref_predicted_temperate == "no" &
                                                          predicted_mode_analysis$ref_mode_prediction_mode_approx == "hgcf") |
                                                         (predicted_mode_analysis$query_predicted_temperate == "no" &
                                                            predicted_mode_analysis$query_mode_prediction_mode_approx == "hgcf")))


predicted_mode_analysis_firm_lytic_hgcf <- subset(predicted_mode_analysis,predicted_mode_analysis$ref_host_phylum == "Firmicutes" &
                                                     predicted_mode_analysis$query_host_phylum == "Firmicutes" &
                                                     ((predicted_mode_analysis$ref_predicted_temperate == "no" &
                                                         predicted_mode_analysis$ref_mode_prediction_mode_approx == "hgcf") |
                                                        (predicted_mode_analysis$query_predicted_temperate == "no" &
                                                           predicted_mode_analysis$query_mode_prediction_mode_approx == "hgcf")))

predicted_mode_analysis_proteo_lytic_hgcf <- subset(predicted_mode_analysis,predicted_mode_analysis$ref_host_phylum == "Proteobacteria" &
                                                     predicted_mode_analysis$query_host_phylum == "Proteobacteria" &
                                                     ((predicted_mode_analysis$ref_predicted_temperate == "no" &
                                                         predicted_mode_analysis$ref_mode_prediction_mode_approx == "hgcf") |
                                                        (predicted_mode_analysis$query_predicted_temperate == "no" &
                                                           predicted_mode_analysis$query_mode_prediction_mode_approx == "hgcf")))

  
par(mar=c(4,8,4,4))
plot(predicted_mode_analysis_lytic_hgcf$modified_mash_distance,predicted_mode_analysis_lytic_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


par(mar=c(4,8,4,4))
plot(predicted_mode_analysis_bacter_lytic_hgcf$modified_mash_distance,predicted_mode_analysis_bacter_lytic_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(predicted_mode_analysis_firm_lytic_hgcf$modified_mash_distance,predicted_mode_analysis_firm_lytic_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(predicted_mode_analysis_proteo_lytic_hgcf$modified_mash_distance,predicted_mode_analysis_proteo_lytic_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Check distribution by host phylum and predicted lifestyle
predicted_mode_analysis_actino_predicted_lytic <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Actinobacteria" &
                                                      predicted_mode_analysis$phage_predicted_temperate_compare == "no")
predicted_mode_analysis_actino_predicted_temperate <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Actinobacteria" &
                                                           predicted_mode_analysis$phage_predicted_temperate_compare == "yes")

predicted_mode_analysis_bacter_predicted_lytic <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Bacteroidetes" &
                                                           predicted_mode_analysis$phage_predicted_temperate_compare == "no")
predicted_mode_analysis_bacter_predicted_temperate <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Bacteroidetes" &
                                                               predicted_mode_analysis$phage_predicted_temperate_compare == "yes")

predicted_mode_analysis_cyano_predicted_lytic <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Cyanobacteria" &
                                                           predicted_mode_analysis$phage_predicted_temperate_compare == "no")
predicted_mode_analysis_cyano_predicted_temperate <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Cyanobacteria" &
                                                               predicted_mode_analysis$phage_predicted_temperate_compare == "yes")

predicted_mode_analysis_firm_predicted_lytic <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Firmicutes" &
                                                           predicted_mode_analysis$phage_predicted_temperate_compare == "no")
predicted_mode_analysis_firm_predicted_temperate <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Firmicutes" &
                                                               predicted_mode_analysis$phage_predicted_temperate_compare == "yes")

predicted_mode_analysis_proteo_predicted_lytic <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Proteobacteria" &
                                                           predicted_mode_analysis$phage_predicted_temperate_compare == "no")
predicted_mode_analysis_proteo_predicted_temperate <- subset(predicted_mode_analysis,predicted_mode_analysis$host_phylum_compare == "Proteobacteria" &
                                                               predicted_mode_analysis$phage_predicted_temperate_compare == "yes")

par(mar=c(4,8,4,4))
plot(predicted_mode_analysis_cyano_predicted_lytic$modified_mash_distance,predicted_mode_analysis_cyano_predicted_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(predicted_mode_analysis_cyano_predicted_temperate$modified_mash_distance,predicted_mode_analysis_cyano_predicted_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
abline(0,2,lty=2,lwd=3,col="grey")










#Estimate how many different clusters are hgcf and lgcf
library(plyr)
library(data.table)

#Create list of unique clusters
host_mode_table_actino <- subset(host_mode_table,host_mode_table$phage_cluster_source == 'actino')

host_mode_table_actino$phage_cluster <- factor(host_mode_table_actino$phage_cluster)
actino782_cluster_frequency <- as.data.frame(levels(host_mode_table_actino$phage_cluster))
names(actino782_cluster_frequency) = c('cluster')




#Original version to compute, but only quantifies for temperate lifestyle
# predicted_temp_hgcf_approx <- subset(predicted_temperate_mode,predicted_temperate_mode$mode_prediction_mode_approx == 'hgcf' & is.na(predicted_temperate_mode$phage_cluster) == FALSE)
# predicted_temp_lgcf_approx <- subset(predicted_temperate_mode,predicted_temperate_mode$mode_prediction_mode_approx == 'lgcf'& is.na(predicted_temperate_mode$phage_cluster) == FALSE)
# predicted_temp_hgcf_approx$phage_cluster <- factor(predicted_temp_hgcf_approx$phage_cluster)
# predicted_temp_lgcf_approx$phage_cluster <- factor(predicted_temp_lgcf_approx$phage_cluster)
# 
# #Count frequency of each cluster for each mode
# total_freq <- count(host_mode_table_actino,c('phage_cluster'))
# names(total_freq) <- c('cluster','total_frequency')
# predicted_temp_hgcf_approx_cluster_freq <- count(predicted_temp_hgcf_approx,c('phage_cluster'))
# names(predicted_temp_hgcf_approx_cluster_freq) <- c('cluster','predicted_temp_hgcf_approx_frequency')
# predicted_temp_lgcf_approx_cluster_freq <- count(predicted_temp_lgcf_approx,c('phage_cluster'))
# names(predicted_temp_lgcf_approx_cluster_freq) <- c('cluster','predicted_temp_lgcf_approx_frequency')
# 
# #Now combine data
# actino782_cluster_frequency <- merge(actino782_cluster_frequency,total_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
# actino782_cluster_frequency <- merge(actino782_cluster_frequency,predicted_temp_hgcf_approx_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
# actino782_cluster_frequency <- merge(actino782_cluster_frequency,predicted_temp_lgcf_approx_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
# 
# 
# #Convert NA values to 0
# actino782_cluster_frequency[is.na(actino782_cluster_frequency)] <- 0



predicted_hgcf_approx <- subset(host_mode_table_actino,host_mode_table_actino$mode_prediction_mode_approx == 'hgcf')
predicted_lgcf_approx <- subset(host_mode_table_actino,host_mode_table_actino$mode_prediction_mode_approx == 'lgcf')
predicted_unknown_approx <- subset(host_mode_table_actino,host_mode_table_actino$mode_prediction_mode_approx == 'unknown')
predicted_hgcf_approx$phage_cluster <- factor(predicted_hgcf_approx$phage_cluster)
predicted_lgcf_approx$phage_cluster <- factor(predicted_lgcf_approx$phage_cluster)
predicted_unknown_approx$phage_cluster <- factor(predicted_unknown_approx$phage_cluster)



#Count frequency of each cluster for each mode
total_freq <- count(host_mode_table_actino,c('phage_cluster'))
names(total_freq) <- c('cluster','total_frequency')
predicted_hgcf_approx_cluster_freq <- count(predicted_hgcf_approx,c('phage_cluster'))
names(predicted_hgcf_approx_cluster_freq) <- c('cluster','predicted_hgcf_approx_frequency')
predicted_lgcf_approx_cluster_freq <- count(predicted_lgcf_approx,c('phage_cluster'))
names(predicted_lgcf_approx_cluster_freq) <- c('cluster','predicted_lgcf_approx_frequency')
predicted_unknown_approx_cluster_freq <- count(predicted_unknown_approx,c('phage_cluster'))
names(predicted_unknown_approx_cluster_freq) <- c('cluster','predicted_unknown_approx_frequency')




#Now combine data
actino782_cluster_frequency <- merge(actino782_cluster_frequency,total_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
actino782_cluster_frequency <- merge(actino782_cluster_frequency,predicted_hgcf_approx_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
actino782_cluster_frequency <- merge(actino782_cluster_frequency,predicted_lgcf_approx_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
actino782_cluster_frequency <- merge(actino782_cluster_frequency,predicted_unknown_approx_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)

#Convert NA values to 0
actino782_cluster_frequency[is.na(actino782_cluster_frequency)] <- 0


#Compute the % of HGCF phages per cluster
actino782_cluster_frequency$predicted_hgcf_approx_cluster_percent <- actino782_cluster_frequency$predicted_hgcf_approx_frequency / actino782_cluster_frequency$total_frequency
actino782_cluster_frequency$predicted_lgcf_approx_cluster_percent <- actino782_cluster_frequency$predicted_lgcf_approx_frequency / actino782_cluster_frequency$total_frequency
actino782_cluster_frequency$predicted_unknown_approx_cluster_percent <- actino782_cluster_frequency$predicted_unknown_approx_frequency / actino782_cluster_frequency$total_frequency



#Now compute the % phages/cluster that are predicted to be temperate or lytic


host_mode_table_actino_predicted_temp <- subset(host_mode_table_actino,host_mode_table_actino$predicted_temperate == "yes")
host_mode_table_actino_predicted_lytic <- subset(host_mode_table_actino,host_mode_table_actino$predicted_temperate == "no")
host_mode_table_actino_predicted_temp$phage_cluster <- factor(host_mode_table_actino_predicted_temp$phage_cluster)
host_mode_table_actino_predicted_lytic$phage_cluster <- factor(host_mode_table_actino_predicted_lytic$phage_cluster)

host_mode_table_actino_predicted_temp_cluster_freq <- count(host_mode_table_actino_predicted_temp,c('phage_cluster'))
names(host_mode_table_actino_predicted_temp_cluster_freq) <- c('cluster','predicted_temperate_frequency')
host_mode_table_actino_predicted_lytic_cluster_freq <- count(host_mode_table_actino_predicted_lytic,c('phage_cluster'))
names(host_mode_table_actino_predicted_lytic_cluster_freq) <- c('cluster','predicted_lytic_frequency')

actino782_cluster_frequency <- merge(actino782_cluster_frequency,host_mode_table_actino_predicted_temp_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)
actino782_cluster_frequency <- merge(actino782_cluster_frequency,host_mode_table_actino_predicted_lytic_cluster_freq,by.x='cluster',by.y='cluster',all.x=TRUE)

actino782_cluster_frequency[is.na(actino782_cluster_frequency)] <- 0
actino782_cluster_frequency$predicted_temperate_cluster_percent <- actino782_cluster_frequency$predicted_temperate_frequency / actino782_cluster_frequency$total_frequency
actino782_cluster_frequency$predicted_lytic_cluster_percent <- actino782_cluster_frequency$predicted_lytic_frequency / actino782_cluster_frequency$total_frequency







#Export to file for analysis in Excel...
write.table(actino782_cluster_frequency,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170316_cluster_mode_data.csv",sep=",",row.names = FALSE,quote=FALSE)


#Create bar plot
actino782_cluster_frequency_sorted <- actino782_cluster_frequency[order(actino782_cluster_frequency$predicted_temp_hgcf_approx_cluster_percent,actino782_cluster_frequency$predicted_temp_lgcf_approx_cluster_percent),]
actino782_cluster_frequency_matrix <- as.matrix(actino782_cluster_frequency_sorted[,6:8])
rownames(actino782_cluster_frequency_matrix) <- actino782_cluster_frequency_sorted$cluster
actino782_cluster_frequency_matrix_trans <- t(actino782_cluster_frequency_matrix)
barplot(actino782_cluster_frequency_matrix_trans,col=c('blue','green','white'),horiz=TRUE,las=1,cex.names=0.5)


#Specifically check a few clusters
cluster_p <- subset(mash_table2,mash_table2$ref_phage_cluster == 'P' | mash_table2$query_phage_cluster == 'P')
cluster_p_intra <- subset(mash_table2,mash_table2$phage_cluster_compare == 'P')


par(mar=c(4,8,4,4))
plot(cluster_p$modified_mash_distance,cluster_p$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(cluster_p_intra$modified_mash_distance,cluster_p_intra$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")


cluster_j <- subset(mash_table2,mash_table2$ref_phage_cluster == 'J' | mash_table2$query_phage_cluster == 'J')
cluster_j_intra <- subset(mash_table2,mash_table2$phage_cluster_compare == 'J')


par(mar=c(4,8,4,4))
plot(cluster_j$modified_mash_distance,cluster_j$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(cluster_j_intra$modified_mash_distance,cluster_j_intra$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")












###Cluster A genometrics



cluster_a_genometrics <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170328_cluster_a_genometrics.csv",sep=",",header=TRUE)
















par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$size ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)
par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$total_gene_count ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)
par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$Unspecified_gene_count ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)
par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$lysis_gene_count ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)
par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$lysogeny_gene_count ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)
par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$recombination_replication_gene_count ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)
par(mar=c(8,20,4,4))
boxplot(cluster_a_genometrics$structure_assembly_gene_count ~ cluster_a_genometrics$subcluster,las=2,col=c("cyan","dark green"),cex.axis=2)












myco_genometrics <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170328_myco_genometrics.csv",sep=",",header=TRUE)
myco_genometrics$subcluster <- factor(myco_genometrics$subcluster,c("A1","non-A1","myco"))

myco_gc <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170328_myco_gc.csv",sep=",",header=TRUE)
myco_gc$group <- factor(myco_gc$group,c("A1","non-A1","myco"))





par(mar=c(8,24,4,4))
boxplot(myco_genometrics$size ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$size ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_genometrics$total_gene_count ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$total_gene_count ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_genometrics$Unspecified_gene_count ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$Unspecified_gene_count ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_genometrics$lysis_gene_count ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$lysis_gene_count ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_genometrics$lysogeny_gene_count ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$lysogeny_gene_count ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_genometrics$recombination_replication_gene_count ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$recombination_replication_gene_count ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_genometrics$structure_assembly_gene_count ~ myco_genometrics$subcluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_genometrics$structure_assembly_gene_count ~ myco_genometrics$subcluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))

par(mar=c(8,24,4,4))
boxplot(myco_gc$GC ~ myco_gc$group,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(myco_gc$GC ~ myco_gc$group,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("cyan","dark green","dark grey"))





par(mar=c(8,20,4,4))
boxplot(myco_gc$GC ~ myco_gc$group,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)




# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$size ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$total_gene_count ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$Unspecified_gene_count ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$lysis_gene_count ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$lysogeny_gene_count ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$recombination_replication_gene_count ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(myco_genometrics$structure_assembly_gene_count ~ myco_genometrics$subcluster,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)


# par(mar=c(8,20,4,4))
# boxplot(myco_gc$GC ~ myco_gc$group,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)




















###Actino782 genometrics



actino782_genometrics <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170329_actino782_genometrics.csv",sep=",",header=TRUE)
actino782_genometrics$category_for_graphing <- factor(actino782_genometrics$category_for_graphing,c("hgcf","lgcf","lytic"))


par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))

par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$total_gene_count ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$total_gene_count ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))

par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$Unspecified_gene_count ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$Unspecified_gene_count ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))

par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$lysis_gene_count ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$lysis_gene_count ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))

par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$lysogeny_gene_count ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$lysogeny_gene_count ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))

par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$recombination_replication_gene_count ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$recombination_replication_gene_count ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))

par(mar=c(8,24,4,4))
boxplot(actino782_genometrics$structure_assembly_gene_count ~ actino782_genometrics$category_for_graphing,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(actino782_genometrics$structure_assembly_gene_count ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))



















# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$total_gene_count ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$Unspecified_gene_count ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$lysis_gene_count ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$lysogeny_gene_count ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$recombination_replication_gene_count ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$structure_assembly_gene_count ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2)





# library(vioplot)
# 
# actino782_genometrics_hgcf <- subset(actino782_genometrics,actino782_genometrics$category_for_graphing == "hgcf")
# actino782_genometrics_lgcf <- subset(actino782_genometrics,actino782_genometrics$category_for_graphing == "lgcf")
# actino782_genometrics_lytic <- subset(actino782_genometrics,actino782_genometrics$category_for_graphing == "lytic")
# 
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$size,actino782_genometrics_lgcf$size,actino782_genometrics_lytic$size,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$total_gene_count,actino782_genometrics_lgcf$total_gene_count,actino782_genometrics_lytic$total_gene_count,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$Unspecified_gene_count,actino782_genometrics_lgcf$Unspecified_gene_count,actino782_genometrics_lytic$Unspecified_gene_count,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$lysis_gene_count,actino782_genometrics_lgcf$lysis_gene_count,actino782_genometrics_lytic$lysis_gene_count,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$lysogeny_gene_count,actino782_genometrics_lgcf$lysogeny_gene_count,actino782_genometrics_lytic$lysogeny_gene_count,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$recombination_replication_gene_count,actino782_genometrics_lgcf$recombination_replication_gene_count,actino782_genometrics_lytic$recombination_replication_gene_count,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))
# par(mar=c(8,20,4,4))
# vioplot(actino782_genometrics_hgcf$structure_assembly_gene_count,actino782_genometrics_lgcf$structure_assembly_gene_count,actino782_genometrics_lytic$structure_assembly_gene_count,names=c("hgcf","lgcf","lytic"),col=c("blue","green","red"))





# par(mar=c(8,20,4,4))
# stripchart(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","green","red"))
# boxplot(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,add=TRUE,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0)
# 
# par(mar=c(8,20,4,4))
# boxplot(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,las=2,col=c("blue","green","red"),cex.axis=2,whisklty=0,staplelty=0,range=0)
# stripchart(actino782_genometrics$size ~ actino782_genometrics$category_for_graphing,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5)




###Misc clusters genometrics


misc_clusters_genometrics <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170403_misc_clusters_genometrics.csv",sep=",",header=TRUE)
misc_clusters_genometrics$phage_cluster <- factor(misc_clusters_genometrics$phage_cluster,c("A1","F","BD","K","non-A1","B"))

misc_clusters_gc <- read.csv("/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/current/20170403_misc_clusters_gc.csv",sep=",",header=TRUE)
misc_clusters_gc$Cluster <- factor(misc_clusters_gc$Cluster,c("A1","F","BD","K","non-A1","B"))



par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$size ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$size ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$total_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$total_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$Unspecified_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$Unspecified_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$lysis_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$lysis_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$lysogeny_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$lysogeny_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$recombination_replication_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$recombination_replication_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_genometrics$structure_assembly_gene_count ~ misc_clusters_genometrics$phage_cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_genometrics$structure_assembly_gene_count ~ misc_clusters_genometrics$phage_cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))

par(mar=c(8,24,4,4))
boxplot(misc_clusters_gc$GC ~ misc_clusters_gc$Cluster,las=2,cex.axis=2,whisklty=0,staplelty=0,range=0,varwidth=FALSE)
stripchart(misc_clusters_gc$GC ~ misc_clusters_gc$Cluster,add=TRUE,vertical=TRUE,method="jitter",pch=19,cex=0.5,col=c("blue","blue","green","green","green","red"))



par(mar=c(8,20,4,4))
boxplot(myco_gc$GC ~ myco_gc$group,las=2,col=c("cyan","dark green","dark grey"),cex.axis=2)



















###scratch

actino_only <- subset(mash_table2,mash_table2$host_phylum_compare == 'Actinobacteria')

brusacoram <- subset(actino_only,actino_only$mash_reference == 'brusacoram__actino785' | actino_only$mash_query == 'brusacoram__actino785',
                     select=c('mash_reference','mash_query','mash_count','modified_mash_distance','pham_pham_dissimilarity','ref_size','query_size'))


iracema64 <- subset(actino_only,actino_only$mash_reference == 'iracema64__actino785' | actino_only$mash_query == 'iracema64__actino785',
                     select=c('mash_reference','mash_query','mash_count','modified_mash_distance','pham_pham_dissimilarity','ref_size','query_size'))


alvin <- subset(mash_table2,mash_table2$mash_reference == 'alvin__actino785' | mash_table2$mash_query == 'alvin__actino785',
                    select=c('mash_reference','mash_query','mash_count','modified_mash_distance','pham_pham_dissimilarity','ref_size','query_size'))

alvin_subset <- subset(alvin,alvin$modified_mash_distance < 0.42 &
                         alvin$pham_pham_dissimilarity < 0.89)

alvin_subset <- subset(alvin,alvin$modified_mash_distance < 0.42 &
                         alvin$modified_mash_distance >= 0.06 &
                         alvin$modified_mash_distance < 0.28 &
                         alvin$pham_pham_dissimilarity < 0.89 &
                         alvin$pham_pham_dissimilarity > 0.22 &
                         alvin$pham_pham_dissimilarity < 0.79)


par(mar=c(4,8,4,4))
plot(alvin$modified_mash_distance,alvin$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)

par(mar=c(4,8,4,4))
plot(alvin_subset$modified_mash_distance,alvin_subset$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)


abline(v=0.42,lty=2,lwd=3,col="grey")
abline(v=0.06,lty=2,lwd=3,col="grey")
abline(v=0.28,lty=2,lwd=3,col="grey")

abline(h=0.89,lty=2,lwd=3,col="grey")
abline(h=0.22,lty=2,lwd=3,col="grey")
abline(h=0.79,lty=2,lwd=3,col="grey")


abline(0,3.5,lty=2,lwd=3,col="grey")
abline(0.25,2,lty=2,lwd=3,col="grey")



abline(0,2,lty=2,lwd=3,col="grey")

abline(v=0.16,lty=2,lwd=3,col="grey")
abline(v=0.42,lty=2,lwd=3,col="grey")





















