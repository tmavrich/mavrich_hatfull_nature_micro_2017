#R script to perform misc. analyses on Mash and pham data
#Travis Mavrich
#Note: this code is not intended to be run from start to finish at once.
#Instead, several sections require exporting data for other code/software or
#importing additional data for analysis.


#Function to import Mash data
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






#Import mash dataset
#Format:
#0 = reference
#1 = query
#2 = mash distance 
#3 = mash p-value
#4 = mash kmer count
mash_table <- import_function("processed_mash_output.csv","mash")


#Convert fields from character class (default) to factor class
mash_table$mash_reference <- as.factor(mash_table$mash_reference)
mash_table$mash_query <- as.factor(mash_table$mash_query)
mash_table$mash_ref_query <- as.factor(mash_table$mash_ref_query)




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
#18 = size
#19 = gene_count
#20 = toxic
#21 = predicted_temperate
host_table <- read.csv("phage_host_data.csv",sep=",",header=TRUE)


#Convert all Unspecified fields to NA missing value
host_table[host_table == "Unspecified"] <- NA



#Modify column names of host data and merge with kmer table
#Data for all phages in processed mash table should be present in host and phage metadata table.
#As a result, do not select all.x=TRUE option. This way, any missing rows indicates an error in the table.


#Match metadata for both the query and reference phages in each pairwise comparison
names(host_table) <- c("phage_identifier","query_header_source_info","query_database",
                             "query_host_superkingdom","query_host_phylum","query_host_class","query_host_order","query_host_family","query_host_genus",
                             "query_phage_superkingdom","query_phage_viral_type","query_phage_order","query_phage_family","query_phage_genus",
                             "query_phage_cluster","query_phage_subcluster","query_phage_cluster_source",
                             "query_phage_temperate","query_size","query_gene_count",
                             "query_toxic","query_predicted_temperate")


mash_table2 <- merge(mash_table,host_table,by.x="mash_query",by.y="phage_identifier")

names(host_table) <- c("phage_identifier","ref_header_source_info","ref_database",
                       "ref_host_superkingdom","ref_host_phylum","ref_host_class","ref_host_order","ref_host_family","ref_host_genus",
                       "ref_phage_superkingdom","ref_phage_viral_type","ref_phage_order","ref_phage_family","ref_phage_genus",
                       "ref_phage_cluster","ref_phage_subcluster","ref_phage_cluster_source",
                       "ref_phage_temperate","ref_size","ref_gene_count",
                       "ref_toxic","ref_predicted_temperate")



mash_table2 <- merge(mash_table2,host_table,by.x="mash_reference",by.y="phage_identifier")

names(host_table) <- c("phage_identifier","header_source_info","database",
                       "host_superkingdom","host_phylum","host_class","host_order","host_family","host_genus",
                       "phage_superkingdom","phage_viral_type","phage_order","phage_family","phage_genus",
                       "phage_cluster","phage_subcluster","phage_cluster_source",
                       "phage_temperate","size","gene_count",
                       "toxic","predicted_temperate")



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
#Format:
#0 = ani_ref_query
#1 = ani_reference
#2 = ani_query
#3 = ani_ani
#4 = ani_distance
ani_data <- read.csv("ani_data.csv",sep=",",header=TRUE)
ani_data$ani_ref_query <- as.character(ani_data$ani_ref_query)
ani_data$ani_distance <- 1 - ani_data$ani_ani

#Merge with the main table
mash_table2 <- merge(mash_table2,ani_data,by.x="mash_ref_query",by.y="ani_ref_query",all.x=TRUE)



#Import pham data and merge with kmer table
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
#But you don't want to retain all rows = making scatter plots or histograms can cause errors if not all rows have data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
mash_table2 <- merge(mash_table2,pham_table,by.x="mash_ref_query",by.y="pham_phage1_phage2")



#Assign filter status and change mash distance if data is not significant
#Alternatively, the max percent parameter can be omitted with minimal change to final analysis.
mash_table2$filter <- ifelse(mash_table2$mash_pvalue < 1e-10 & mash_table2$size_diff_max_percent < 1,TRUE,FALSE)

#At this point, the max mash distance of all filtered comparisons = 0.4903. So set the distance of all comparisons that did not pass the filter = 0.5
mash_table2$modified_mash_distance <- ifelse(mash_table2$filter == TRUE,mash_table2$mash_distance,0.5)





#Assign evolutionary mode to each pairwise comparison.
#This assigns all data in the dataset. However, evolutionary mode is only relevant to dsDNA
mash_table2$gene_flux_part1 <- ifelse(mash_table2$modified_mash_distance < 0.1666667 & mash_table2$pham_pham_dissimilarity > (mash_table2$modified_mash_distance * 3.5),TRUE,FALSE)
mash_table2$gene_flux_part2 <- ifelse(mash_table2$modified_mash_distance > 0.1666667 & mash_table2$pham_pham_dissimilarity > (mash_table2$modified_mash_distance * 2 + 0.25),TRUE,FALSE)
mash_table2$gene_flux_category <- ifelse(mash_table2$gene_flux_part1 == TRUE | mash_table2$gene_flux_part2 ==TRUE,"high","low")
mash_table2$gene_flux_category <- as.factor(mash_table2$gene_flux_category)




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








#Compare ref and query host and phage metadata columns
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
mash_table2$phage_toxic_compare <- ifelse(mash_table2$ref_toxic==mash_table2$query_toxic,as.character(mash_table2$ref_toxic),"different")
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
mash_table2$phage_toxic_compare <- as.factor(mash_table2$phage_toxic_compare)
mash_table2$phage_predicted_temperate_compare <- as.factor(mash_table2$phage_predicted_temperate_compare)


#At this point, all imported data has been processed. From here on out, commands are mostly analysis
#based on subsets of the data, although for some specific analyses dataset are imported or exported. 











###All dsDNA phages
bacteria <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria")
type_dsDNA <- subset(bacteria,bacteria$phage_viral_type_compare == "dsDNA")

#Fig. 1a
par(mar=c(4,8,4,4))
plot(type_dsDNA$modified_mash_distance,type_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig. 1a
par(mar=c(4,8,15,4))
hist(type_dsDNA$modified_mash_distance,breaks=((range(type_dsDNA$modified_mash_distance)[2]-range(type_dsDNA$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,3e4),col="black",cex.axis=2)

#Fig. 1a
par(mar=c(4,4,15,4))
hist(type_dsDNA$pham_pham_dissimilarity,breaks=((range(type_dsDNA$pham_pham_dissimilarity)[2]-range(type_dsDNA$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)











###Check by empirical lifestyle.
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")

all_temperate <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "yes")
all_lytic <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "no")
all_different <- subset(type_dsDNA,type_dsDNA$phage_temperate_compare == "different")


#Fig. 1c
par(mar=c(4,8,4,4))
plot(all_temperate$modified_mash_distance,all_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig. 1c
par(mar=c(4,8,4,4))
plot(all_lytic$modified_mash_distance,all_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig. 1c
par(mar=c(4,8,4,4))
plot(all_different$modified_mash_distance,all_different$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")









###Cluster-specific analysis
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
#1 = phage_cluster
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
#0 = phageName
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


































###By host phylum, then by lifestyle
#Only focus on dsDNA phages, so remove other phage types.
bacteria <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria")
type_dsDNA <- subset(bacteria,bacteria$phage_viral_type_compare == "dsDNA")

host_phylum <- subset(type_dsDNA,type_dsDNA$host_phylum_compare != "different")
host_phylum_diff <- subset(type_dsDNA,type_dsDNA$host_phylum_compare == "different")

actino <- subset(host_phylum,host_phylum$host_phylum_compare == "Actinobacteria")
bacter <- subset(host_phylum,host_phylum$host_phylum_compare == "Bacteroidetes")
cyano <- subset(host_phylum,host_phylum$host_phylum_compare == "Cyanobacteria")
firm <- subset(host_phylum,host_phylum$host_phylum_compare == "Firmicutes")
proteo <- subset(host_phylum,host_phylum$host_phylum_compare == "Proteobacteria")








#Same phylum
#Scatter plots of Mash vs Pham distances by host phyla


#Fig 4a
par(mar=c(4,8,4,4))
plot(actino$modified_mash_distance,actino$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig 4b
par(mar=c(4,8,4,4))
plot(bacter$modified_mash_distance,bacter$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig 4c
par(mar=c(4,8,4,4))
plot(cyano$modified_mash_distance,cyano$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig 4d
par(mar=c(4,8,4,4))
plot(firm$modified_mash_distance,firm$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Fig 4e
par(mar=c(4,8,4,4))
plot(proteo$modified_mash_distance,proteo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




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
par(mar=c(4,8,4,4))
plot(host_phylum_diff$modified_mash_distance,host_phylum_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


















































###pham general vs jaccard dissimilarity

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
bacteria_dsDNA <- subset(type_dsDNA,type_dsDNA$host_superkingdom_compare == 'Bacteria')


#Check how correlated pham dissimilarity and jaccard dissimilarity are
#Supp. Fig. 2a
par(mar=c(4,8,4,4))
plot(mash_table2$pham_jaccard_dissimilarity,mash_table2$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")


#Mash vs Pham plot using jaccard
#Supp. Fig. 2b
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###ANI vs Pham Dissimilarity data
#Import ANI data from the 79 genome optimization test to show ANI vs Pham Distance
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
ani79_analysis <- merge(mash_table2,ani79_data,by.x="mash_ref_query",by.y="ani79_ref_query")

#Supp. Fig. 2c
par(mar=c(4,8,4,4))
plot(ani79_analysis$modified_mash_distance,ani79_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp. Fig. 2c
par(mar=c(4,8,4,4))
plot(ani79_analysis$ani79_ani_distance,ani79_analysis$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")







###Compare Mash to ANI
#Compare patterns from pham-mash and pham-ani comparisons
ani_mash <- subset(mash_table2,complete.cases(mash_table2$ani_distance))

#Colored plot by assigned evolutionary mode 
ani_mash_hgcf <- subset(ani_mash,ani_mash$gene_flux_category == "high" & ani_mash$phage_predicted_temperate_compare == "yes")
ani_mash_lgcf <- subset(ani_mash,ani_mash$gene_flux_category == "low" & ani_mash$phage_predicted_temperate_compare == "yes")
ani_mash_lytic <- subset(ani_mash,ani_mash$phage_predicted_temperate_compare == "no")


#Supp. Fig. 2d
par(mar=c(4,8,4,4))
plot(ani_mash_hgcf$modified_mash_distance,ani_mash_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(ani_mash_lgcf$modified_mash_distance,ani_mash_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(ani_mash_lytic$modified_mash_distance,ani_mash_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")

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
#Supp. Fig. 2e
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$pham_pham_dissimilarity,bacteria_dsDNA$vog_gene_content_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")



#Compare pham-based and vog-based bacteria dsDNA phage lifestyle plots 
par(mar=c(4,8,4,4))
plot(temperate_both$modified_mash_distance,temperate_both$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 2f
par(mar=c(4,8,4,4))
plot(temperate_both$modified_mash_distance,temperate_both$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(temperate_neither$modified_mash_distance,temperate_neither$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 2f
par(mar=c(4,8,4,4))
plot(temperate_neither$modified_mash_distance,temperate_neither$vog_gene_content_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")








###Plot by phage type
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_diff <- subset(mash_table2,mash_table2$phage_viral_type_compare == "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
type_dsRNA <- subset(type,type$phage_viral_type_compare == "dsRNA")
type_ssDNA <- subset(type,type$phage_viral_type_compare == "ssDNA")
type_ssRNA <- subset(type,type$phage_viral_type_compare == "ssRNA")


#Compare similarity of different types

#Supp Fig. 5b
par(mar=c(4,8,4,4))
plot(type_diff$modified_mash_distance,type_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Now plot for each nucleic acid type

#Supp Fig. 4
par(mar=c(4,8,4,4))
plot(type_dsDNA$modified_mash_distance,type_dsDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 4
par(mar=c(4,8,4,4))
plot(type_dsRNA$modified_mash_distance,type_dsRNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 4
par(mar=c(4,8,4,4))
plot(type_ssDNA$modified_mash_distance,type_ssDNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


#Supp Fig. 4
par(mar=c(4,8,4,4))
plot(type_ssRNA$modified_mash_distance,type_ssRNA$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")









###Check Eukaryotic controls
euk_check <- subset(mash_table2,mash_table2$ref_host_superkingdom == "Eukaryota" | mash_table2$query_host_superkingdom == "Eukaryota")
euk_check <- subset(euk_check,euk_check$ref_host_superkingdom == "Bacteria" | euk_check$query_host_superkingdom == "Bacteria")

compute_sector_distribution(euk_check)

#No figure
par(mar=c(4,8,4,4))
plot(euk_check$modified_mash_distance,euk_check$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



###Check archaea controls
controls <- mash_table2
controls$archaea_one_or_two <- ifelse(controls$ref_host_superkingdom == "Archaea" | controls$query_host_superkingdom == "Archaea",TRUE,FALSE)
controls$archaea_both <- ifelse(controls$ref_host_superkingdom == "Archaea" & controls$query_host_superkingdom == "Archaea",TRUE,FALSE)
controls$archaea_one <- ifelse(controls$archaea_one_or_two == TRUE & controls$archaea_both == FALSE,TRUE,FALSE)
controls$bacteria_one_or_two <- ifelse(controls$ref_host_superkingdom == "Bacteria" | controls$query_host_superkingdom == "Bacteria",TRUE,FALSE)
controls$archaea_one_bacteria_one <- ifelse(controls$archaea_one == TRUE & controls$bacteria_one_or_two == TRUE,TRUE,FALSE)


archaea_check <- subset(controls,controls$archaea_one_bacteria_one == TRUE)
compute_sector_distribution(archaea_check)

#Supp Fig. 5a
par(mar=c(4,8,4,4))
plot(archaea_check$modified_mash_distance,archaea_check$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
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

#Supp. Fig. 5d
par(mar=c(4,8,4,4))
plot(cluster_actino_same$modified_mash_distance,cluster_actino_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

#Supp. Fig. 5d
par(mar=c(4,8,4,4))
plot(cluster_actino_diff$modified_mash_distance,cluster_actino_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="black")
abline(0,2,lty=2,lwd=3,col="grey")

#Supp. Fig. 5d
par(mar=c(4,8,4,4))
plot(subcluster_actino_same$modified_mash_distance,subcluster_actino_same$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")

#Supp. Fig. 5d
par(mar=c(4,8,4,4))
plot(subcluster_actino_diff$modified_mash_distance,subcluster_actino_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")







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

#Supp Fig. 5e
par(mar=c(4,8,4,4))
plot(actino785_data_same_cluster_neither_subclustered$modified_mash_distance,actino785_data_same_cluster_neither_subclustered$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="orange")
abline(0,2,lty=2,lwd=3,col="grey")






###Compute distance of all Actinobacteriophage Singletons from other Actinobacteriophages
#Histogram = Mash distances for all Actino785 comparisons containing at least one Singleton
actino785_data <- subset(mash_table2,mash_table2$phage_cluster_source_compare == "actino")
actino785_data_same <- subset(actino785_data,actino785_data$phage_cluster_compare != "different")
actino785_data_diff <- subset(actino785_data,actino785_data$phage_cluster_compare == "different")
actino785_data_diff$query_singleton <- grepl("^Singleton",actino785_data_diff$query_phage_cluster)
actino785_data_diff$ref_singleton <- grepl("^Singleton",actino785_data_diff$ref_phage_cluster)
actino785_data_diff$singleton_one_or_two <- ifelse(actino785_data_diff$ref_singleton == TRUE | actino785_data_diff$query_singleton == TRUE,TRUE,FALSE)

actino785_singleton_one_or_two <- subset(actino785_data_diff,actino785_data_diff$singleton_one_or_two == TRUE)

compute_sector_distribution(actino785_singleton_one_or_two)

#Supp Fig. 5f
par(mar=c(4,8,4,4))
plot(actino785_singleton_one_or_two$modified_mash_distance,actino785_singleton_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")










###Analyze predicted lifestyle data
lifestyle_analysis <- mash_table2

bacteria_dsDNA <- subset(lifestyle_analysis,lifestyle_analysis$host_superkingdom_compare == 'Bacteria' & lifestyle_analysis$phage_viral_type_compare == 'dsDNA')

lifestyle_predicted_temperate <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'yes')
lifestyle_predicted_lytic <- subset(bacteria_dsDNA,bacteria_dsDNA$phage_predicted_temperate_compare == 'no')

#Supp. Fig. 6b
par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate$modified_mash_distance,lifestyle_predicted_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp. Fig. 6b
par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic$modified_mash_distance,lifestyle_predicted_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





###Predict evolutionary mode

data_for_mode_prediction <- subset(mash_table2,mash_table2$host_superkingdom_compare == "Bacteria" &
                                     mash_table2$phage_viral_type_compare == "dsDNA" &
                                     mash_table2$modified_mash_distance < 0.42 &
                                     mash_table2$pham_pham_dissimilarity < 0.89)
data_for_mode_prediction <- subset(data_for_mode_prediction,select=c("mash_reference","mash_query","modified_mash_distance","pham_pham_dissimilarity"))
write.table(data_for_mode_prediction,"/Users/Hatfull_Lab/Desktop/Project_phage_classification/7_merged2333_analysis/20170315_data_for_mode_prediction.csv",sep=",",row.names = FALSE,col.names = FALSE,quote=FALSE)

#Run the data through the analyze_mash_network_script to predict the evolutionary mode



#Then import the mode prediction back into R
#Format
#0 = phage
#1 = hgcf_tally
#2 = lgcf_tally
#3 = out_of_range_tally
#4 = hgcf_percent
#5 = lgcf_percent"
#6 = mode"
mode_prediction_table <- read.csv("mode_prediction.csv",sep=",",header=TRUE)

names(mode_prediction_table) <- c("phage_identifier",
                                  "mode_prediction_hgcf_tally","mode_prediction_lgcf_tally","mode_prediction_out_of_range_tally",
                                  "mode_prediction_hgcf_percent","mode_prediction_lgcf_percent","mode_prediction_mode")

#Supp. Fig. 6d
par(mar=c(4,8,4,4))
hist(mode_prediction_table$mode_prediction_hgcf_percent,col="black",breaks=25,cex.axis=2,ann=FALSE,las=1,ylim=c(0,1400))



















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


#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(myo$modified_mash_distance,myo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(sipho$modified_mash_distance,sipho$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(podo$modified_mash_distance,podo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(podo_temperate$modified_mash_distance,podo_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(podo_lytic$modified_mash_distance,podo_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(podo_lifestyle_diff$modified_mash_distance,podo_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(sipho_temperate$modified_mash_distance,sipho_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(sipho_lytic$modified_mash_distance,sipho_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(sipho_lifestyle_diff$modified_mash_distance,sipho_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(myo_temperate$modified_mash_distance,myo_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(myo_lytic$modified_mash_distance,myo_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7a
par(mar=c(4,8,4,4))
plot(myo_lifestyle_diff$modified_mash_distance,myo_lifestyle_diff$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





#Supp Fig. 7b
par(mar=c(4,8,4,4))
plot(caudovirales_family_diff_sipho_myo$modified_mash_distance,caudovirales_family_diff_sipho_myo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7b
par(mar=c(4,8,4,4))
plot(caudovirales_family_diff_sipho_podo$modified_mash_distance,caudovirales_family_diff_sipho_podo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

#Supp Fig. 7b
par(mar=c(4,8,4,4))
plot(caudovirales_family_diff_myo_podo$modified_mash_distance,caudovirales_family_diff_myo_podo$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")








































































###Export data for gene-specific mash analysis
#Only investigate comparisons that are within cluster boundaries
mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
bacteria_dsDNA <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')

#Cluster-boundary dataset
bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.42 & bacteria_dsDNA$pham_pham_dissimilarity < 0.89)

#Now reduce the data table to only the columns needed for export, and export the data
bacteria_dsDNA_nuc042_gene089_reduced <- subset(bacteria_dsDNA_nuc042_gene089,select = c('mash_reference','mash_query','mash_distance','pham_pham_dissimilarity'))
write.table(bacteria_dsDNA_nuc042_gene089_reduced,"bacteria_dsDNA_nuc042_gene089_data.csv",sep=",",row.names = FALSE,quote=FALSE)







###Gene-specific mash analysis
#After running gene-specific mash analysis, subset out only the comparisons used in the gsm analysis

mash_filtered <- subset(mash_table2,mash_table2$filter == TRUE)
bacteria_dsDNA <- subset(mash_filtered,mash_filtered$host_superkingdom_compare == 'Bacteria' & mash_filtered$phage_viral_type_compare == 'dsDNA')
bacteria_dsDNA_nuc042_gene089 <- subset(bacteria_dsDNA,bacteria_dsDNA$modified_mash_distance < 0.42 & bacteria_dsDNA$pham_pham_dissimilarity < 0.89)

#Rename the reduced data to the generic gsm analysis table name
gene_specific_mash_table <- bacteria_dsDNA_nuc042_gene089


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
#22 = phage1 # unshared phams
#23 = phage1 # all phams
#24 = phage1 # all genes
#25 = phage1 # shared genes
#26 = phage1 # unshared genes
#27 = phage1 average length of all genes
#28 = phage1 average length of shared genes
#29 = phage1 average length of unshared genes
#30 = phage1 total length of all genes
#31 = phage1 total length of shared genes
#32 = phage1 total length of unshared genes
#33 = phage1 all genes GC content
#34 = phage1 shared genes GC content
#35 = phage1 unshared genes GC content
#36 = phage2 # unshared phams
#37 = phage2 # all phams
#38 = phage2 # all genes
#39 = phage2 # shared genes
#40 = phage2 # unshared genes
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

#Create modified data for mash
gene_specific_mash_table$gsm_all_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_all_mash_filter == TRUE,gene_specific_mash_table$gsm_all_mash_distance,0.6)
gene_specific_mash_table$gsm_shared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_shared_mash_filter == TRUE,gene_specific_mash_table$gsm_shared_mash_distance,0.6)
gene_specific_mash_table$gsm_unshared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_unshared_mash_filter == TRUE,gene_specific_mash_table$gsm_unshared_mash_distance,0.6)
gene_specific_mash_table$gsm_shared_unshared_modified_mash_distance <- ifelse(gene_specific_mash_table$gsm_shared_unshared_mash_filter == TRUE,gene_specific_mash_table$gsm_shared_unshared_mash_distance,0.6)


####PROBABLY CAN DELETE THIS AFTER CHECKING AVE GENE SIZE MATH
#Create modified data for average gene size
gene_specific_mash_table$gsm_ref_all_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_ref_num_all_genes > 0,gene_specific_mash_table$gsm_ref_all_ave_size,NA)
gene_specific_mash_table$gsm_ref_shared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_ref_num_shared_genes > 0,gene_specific_mash_table$gsm_ref_shared_ave_size,NA)
gene_specific_mash_table$gsm_ref_unshared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_ref_num_unshared_genes > 0,gene_specific_mash_table$gsm_ref_unshared_ave_size,NA)
gene_specific_mash_table$gsm_query_all_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_query_num_all_genes > 0,gene_specific_mash_table$gsm_query_all_ave_size,NA)
gene_specific_mash_table$gsm_query_shared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_query_num_shared_genes > 0,gene_specific_mash_table$gsm_query_shared_ave_size,NA)
gene_specific_mash_table$gsm_query_unshared_modified_ave_size <- ifelse(gene_specific_mash_table$gsm_query_num_unshared_genes > 0,gene_specific_mash_table$gsm_query_unshared_ave_size,NA)
####


#Estimate of the proportion of coding sequence per genome
#Note: the gene-specific sequence length used does not take into account overlapping CDS features, so the sequences could have duplicate regions in the genome
#Note: the 'all coding potential' data is based on the real genome size, but the 'shared/unshared coding potential' data is based only on the gene-specific mash sizes

gene_specific_mash_table$gsm_all_total_size <- gene_specific_mash_table$gsm_ref_all_total_size + gene_specific_mash_table$gsm_query_all_total_size

gene_specific_mash_table$gsm_ref_unshared_coding_potential <- gene_specific_mash_table$gsm_ref_unshared_total_size/gene_specific_mash_table$gsm_ref_all_total_size
gene_specific_mash_table$gsm_query_unshared_coding_potential <- gene_specific_mash_table$gsm_query_unshared_total_size/gene_specific_mash_table$gsm_query_all_total_size

gene_specific_mash_table$gsm_unshared_coding_potential <- gene_specific_mash_table$gsm_unshared_total_size / gene_specific_mash_table$gsm_all_total_size



####VERIFY THIS MATH IS CORRECT!
#Average gene sizes
#This data does not use the modified gene size data, because if there are no shared/unshared genes, the default ave size is 0, and this is fine to use to compute averages
gene_specific_mash_table$gsm_ave_size_shared_shared <- (gene_specific_mash_table$gsm_ref_shared_ave_size + gene_specific_mash_table$gsm_query_shared_ave_size)/2
gene_specific_mash_table$gsm_ave_size_unshared_unshared <- (gene_specific_mash_table$gsm_ref_unshared_ave_size + gene_specific_mash_table$gsm_query_unshared_ave_size)/2






















#Convert all Unspecified fields to NA missing value
gene_specific_mash_table[gene_specific_mash_table == "Unspecified"] <- NA



#Make backup table before subsetting
gene_specific_mash_table_complete <- gene_specific_mash_table





#Split into HGCF and LGCF datasets

gene_specific_mash_table_temperate_hgcf <- subset(gene_specific_mash_table,gene_specific_mash_table$gene_flux_category == 'high' & gene_specific_mash_table$phage_temperate_compare == 'yes')
gene_specific_mash_table_temperate_lgcf <- subset(gene_specific_mash_table,gene_specific_mash_table$gene_flux_category == 'low' & gene_specific_mash_table$phage_temperate_compare == 'yes')
gene_specific_mash_table_lytic <- subset(gene_specific_mash_table,gene_specific_mash_table$phage_temperate_compare == 'no')






#Analysis after splitting into gene flux modes

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













#How distant are shared and unshared sequences?


#Compare all to shared and unshared by mash distance


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






#Coding potential

#Supp. Fig. 8d
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_lytic$gsm_unshared_coding_potential,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")
abline(0,1,lty=2,lwd=3,col="grey")

#Supp. Fig. 8d
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf$gsm_unshared_coding_potential,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
abline(0,1,lty=2,lwd=3,col="grey")

#Supp. Fig. 8d
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$gsm_unshared_coding_potential,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
abline(0,1,lty=2,lwd=3,col="grey")







#All data points used in gsm plots
#Supp. Fig. 8a
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_hgcf$modified_mash_distance,gene_specific_mash_table_temperate_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_lgcf$modified_mash_distance,gene_specific_mash_table_temperate_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_lytic$modified_mash_distance,gene_specific_mash_table_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")








#Analysis after splitting into gene flux modes and after computing sliding windows
#sliding window average
library(caTools)

gene_specific_mash_table_temperate_hgcf_mmdsort <- gene_specific_mash_table_temperate_hgcf[order(gene_specific_mash_table_temperate_hgcf$modified_mash_distance),]
gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent,101)

gene_specific_mash_table_temperate_lgcf_mmdsort <- gene_specific_mash_table_temperate_lgcf[order(gene_specific_mash_table_temperate_lgcf$modified_mash_distance),]
gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent,101)

gene_specific_mash_table_lytic_mmdsort <- gene_specific_mash_table_lytic[order(gene_specific_mash_table_lytic$modified_mash_distance),]
gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent_runmean <- runmean(gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent,101)


#Compare all genome size disparity
#Supp. Fig. 8a
par(mar=c(4,8,4,4))
plot(gene_specific_mash_table_temperate_lgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_lgcf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(gene_specific_mash_table_temperate_hgcf_mmdsort$modified_mash_distance,gene_specific_mash_table_temperate_hgcf_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(gene_specific_mash_table_lytic_mmdsort$modified_mash_distance,gene_specific_mash_table_lytic_mmdsort$size_diff_ave_percent_runmean,xlim=c(0,0.5),ylim=c(0,0.2),pch=20,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#sliding window analyis based on gene content instead of nucleotide distance
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





####STILL NEED TO FINISH REVIEWING ALL GENE-SPECIFIC-MASH DATA ABOVE





















###Phylogeny comparison 
#Format:
#0 = ref_query
#1 = phylogeny distance
phylogeny_data <- read.csv("phylogeny_data.csv",sep=",",header=TRUE)
phylogeny_analysis <- merge(gene_specific_mash_table,phylogeny_data,by.x="mash_ref_query",by.y="ref_query")



#Compare mash distance to phylogeny distance

#Split into modes to visualize by color
phylogeny_analysis_hgcf <- subset(phylogeny_analysis,phylogeny_analysis$gene_flux_category == "high" & phylogeny_analysis$phage_predicted_temperate_compare == "yes")
phylogeny_analysis_lgcf <- subset(phylogeny_analysis,phylogeny_analysis$gene_flux_category == "low" & phylogeny_analysis$phage_predicted_temperate_compare == "yes")
phylogeny_analysis_lytic <- subset(phylogeny_analysis,phylogeny_analysis$phage_predicted_temperate_compare == "no")



#Supp. Fig. 9a (left)? Not quite right, some datapoints wrong = intra-cluster a points incorrect
par(mar=c(4,8,4,4))
plot(phylogeny_analysis_hgcf$modified_mash_distance,phylogeny_analysis_hgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_analysis_lgcf$modified_mash_distance,phylogeny_analysis_lgcf$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_analysis_lytic$modified_mash_distance,phylogeny_analysis_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")



#Supp. Fig 9a (right) Not quite right, some datapoints wrong = intra-cluster a points incorrect
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


#Fig. 3a (right)
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

  






  
#match up pham proportion data that contains orpham count, pham distribution, etc.
#Columns designated as "pham2" since it is the second set of pham data that have been loaded so far. 
#This second pham data overlaps the first set, but contains additional columns. 
#Also, the extra columns are impacted by the data subset - it only contains comparisons from the actino785 set. 
#So it is not a replacement for the first pham data file.

#Format
#0 = phage1_name"                                              
#1 = phage1_number_of_unshared_phams"                          
#2 = phage1_shared_proportion"                                 
#3 = phage2_name"                                              
#4 = phage2_number_of_unshared_phams"                          
#5 = phage2_shared_proportion"                                 
#6 = number_of_shared_phams"                                   
#7 = average_shared_proportion"                                
#8 = jaccard_similarity"                                       
#9 = shared_pham_distribution_mean"                            
#10 = shared_pham_distribution_median"                          
#11 = shared_pham_distribution_max"                             
#12 = unshared_pham_distribution_mean"                          
#13 = unshared_pham_distribution_median"                        
#14 = unshared_pham_distribution_max"                           
#15 = unshared_orpham_count"                                    
#16 = Unspecified_phage1_number_of_unshared_phams"              
#17 = Unspecified_phage2_number_of_unshared_phams"              
#18 = Unspecified_number_of_shared_phams"                       
#19 = Unspecified_average_shared_proportion"                    
#20 = defense_phage1_number_of_unshared_phams"                  
#21 = defense_phage2_number_of_unshared_phams"                  
#22 = defense_number_of_shared_phams"                           
#23 = defense_average_shared_proportion"                        
#24 = dna_metabolism_phage1_number_of_unshared_phams"           
#25 = dna_metabolism_phage2_number_of_unshared_phams"           
#26 = dna_metabolism_number_of_shared_phams"                    
#27 = dna_metabolism_average_shared_proportion"                 
#28 = lysis_phage1_number_of_unshared_phams"                    
#29 = lysis_phage2_number_of_unshared_phams"                    
#30 = lysis_number_of_shared_phams"                             
#31 = lysis_average_shared_proportion"                          
#32 = lysogeny_phage1_number_of_unshared_phams"                 
#33 = lysogeny_phage2_number_of_unshared_phams"                 
#34 = lysogeny_number_of_shared_phams"                          
#35 = lysogeny_average_shared_proportion"                       
#36 = mobile_phage1_number_of_unshared_phams"                   
#37 = mobile_phage2_number_of_unshared_phams"                   
#38 = mobile_number_of_shared_phams"                            
#39 = mobile_average_shared_proportion"                         
#40 = other_phage1_number_of_unshared_phams"                    
#41 = other_phage2_number_of_unshared_phams"                    
#42 = other_number_of_shared_phams"                             
#43 = other_average_shared_proportion"                          
#44 = recombination_replication_phage1_number_of_unshared_phams"
#45 = recombination_replication_phage2_number_of_unshared_phams"
#46 = recombination_replication_number_of_shared_phams"         
#47 = recombination_replication_average_shared_proportion"      
#48 = structure_assembly_phage1_number_of_unshared_phams"       
#49 = structure_assembly_phage2_number_of_unshared_phams"       
#50 = structure_assembly_number_of_shared_phams"                
#51 = structure_assembly_average_shared_proportion
actino_pham_data <- read.csv("actino_only_pairwise_pham_proportions.csv",sep=",",header=TRUE)

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




#Since pham data contains pairwise duplicates, no need to worry about which phage is which when creating ref_query match column
actino_pham_data$pham2_phage1_phage2 <- paste(actino_pham_data$pham2_phage1,"_",actino_pham_data$pham2_phage2,sep="")
actino_pham_data$pham2_phage1_phage2 <- as.factor(actino_pham_data$pham2_phage1_phage2)





#To retain all rows, be sure to keep all.x=TRUE, but you don't want to retain all rows = making scatter plots or histograms can cause errors if not all rows have data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
phylogeny_analysis <- merge(phylogeny_analysis,actino_pham_data,by.x="mash_ref_query",by.y="pham2_phage1_phage2")


#Split into modes to visualize by color
phylogeny_analysis_hgcf2 <- subset(phylogeny_analysis,phylogeny_analysis$gene_flux_category == "high" & phylogeny_analysis$phage_predicted_temperate_compare == "yes")
phylogeny_analysis_lgcf2 <- subset(phylogeny_analysis,phylogeny_analysis$gene_flux_category == "low" & phylogeny_analysis$phage_predicted_temperate_compare == "yes")
phylogeny_analysis_lytic2 <- subset(phylogeny_analysis,phylogeny_analysis$phage_predicted_temperate_compare == "no")





#sliding window average
#only use data for sliding windows that fit within the scatter plot boundaries = phylogeny distance < 0.3
library(caTools)


phylogeny_hgcf_phylosort <- phylogeny_analysis_hgcf[order(phylogeny_analysis_hgcf2$phylogeny_distance),]
phylogeny_hgcf_phylosort <- phylogeny_hgcf_phylosort[phylogeny_hgcf_phylosort$phylogeny_distance < 0.3,]

phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_hgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_hgcf_phylosort$size_diff_ave_percent,101)
phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_hgcf_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_unshared_coding_potential,101)
phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared,101)



phylogeny_lgcf_phylosort <- phylogeny_analysis_lgcf[order(phylogeny_analysis_lgcf2$phylogeny_distance),]
phylogeny_lgcf_phylosort <- phylogeny_lgcf_phylosort[phylogeny_lgcf_phylosort$phylogeny_distance < 0.3,]

phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lgcf_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lgcf_phylosort$size_diff_ave_percent,101)
phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lgcf_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_unshared_coding_potential,101)
phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared,101)


phylogeny_lytic_phylosort <- phylogeny_analysis_lytic[order(phylogeny_analysis_lytic2$phylogeny_distance),]
phylogeny_lytic_phylosort <- phylogeny_lytic_phylosort[phylogeny_lytic_phylosort$phylogeny_distance < 0.3,]

phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean,101)
phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean <- runmean(phylogeny_lytic_phylosort$pham2_unshared_orpham_count,101)
phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean <- runmean(phylogeny_lytic_phylosort$pham2_pham_dissimilarity,101)
phylogeny_lytic_phylosort$size_diff_ave_percent_runmean <- runmean(phylogeny_lytic_phylosort$size_diff_ave_percent,101)
phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared,101)
phylogeny_lytic_phylosort$gsm_unshared_coding_potential_runmean <- runmean(phylogeny_lytic_phylosort$gsm_unshared_coding_potential,101)
phylogeny_lytic_phylosort$gsm_ave_size_shared_shared_runmean <- runmean(phylogeny_lytic_phylosort$gsm_ave_size_shared_shared,101)






#Gene content dissimilarity
#Supp. Fig. 9b (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_pham_dissimilarity_runmean,xlim=c(0,0.3),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Ave genome size difference
#Supp. Fig. 9b (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$size_diff_ave_percent_runmean,xlim=c(0,0.3),ylim=c(0,0.06),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="red")


#Unshared coding sequence proportion
#Supp. Fig. 9c (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_unshared_coding_potential_runmean,xlim=c(0,0.3),ylim=c(0,0.4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


#Shared/unshared ave gene size
#Supp. Fig. 9c
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_ave_size_shared_shared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")


#Supp. Fig. 9c (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$gsm_ave_size_unshared_unshared_runmean,xlim=c(0,0.3),ylim=c(0,800),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")



#Shared/unshared pham distribution mean
#Supp. Fig. 9c (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_shared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="dark red")

#Supp. Fig. 9c (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_pham_distribution_mean_runmean,xlim=c(0,0.3),ylim=c(0,3),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")

#Unshared orpham count
#Supp. Fig. 9c (but doesn't quite look right)
par(mar=c(4,8,4,4))
plot(phylogeny_hgcf_phylosort$phylogeny_distance,phylogeny_hgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light blue")
par(new=TRUE)
plot(phylogeny_lgcf_phylosort$phylogeny_distance,phylogeny_lgcf_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="light green")
par(new=TRUE)
plot(phylogeny_lytic_phylosort$phylogeny_distance,phylogeny_lytic_phylosort$pham2_unshared_orpham_count_runmean,xlim=c(0,0.3),ylim=c(0,4),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="pink")


















###Scatter plot of all comparisons involving lambda
#Phage identifier for lambda:lambda__nc_001416
#Copy main data table then determine which comparisons involve lambda. Since no self comparisons are present in the dataset, only need to compute if there's 'one' lambda, instead of 'one_or_two' or 'both'
lambda_figure <- mash_table2
lambda_figure$lambda_one <- ifelse(lambda_figure$mash_reference == 'lambda__nc_001416' | lambda_figure$mash_query == 'lambda__nc_001416',TRUE,FALSE)
lambda_comparisons <- subset(lambda_figure,lambda_figure$lambda_one == TRUE)

#Supp. Fig. 11a
par(mar=c(4,8,4,4))
plot(lambda_comparisons$modified_mash_distance,lambda_comparisons$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")





###Toxic phage analysis
mash_table2_toxic <- mash_table2
mash_table2_toxic$toxic_one_or_two <- ifelse(mash_table2_toxic$ref_toxic == "yes" | mash_table2_toxic$query_toxic == "yes",TRUE,FALSE)
mash_table2_toxic$toxic_both <- ifelse(mash_table2_toxic$ref_toxic == "yes" & mash_table2_toxic$query_toxic == "yes",TRUE,FALSE)
toxic_one_or_two <- subset(mash_table2_toxic,mash_table2_toxic$toxic_one_or_two == TRUE)

#Supp. Fig. 11b
par(mar=c(4,8,4,4))
plot(toxic_one_or_two$modified_mash_distance,toxic_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")






