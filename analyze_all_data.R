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





#Decide on which dataset to import
#Main dataset
filename <- "processed_mash_output.csv"


#Import datasets
mash_table <- import_function(filename,"mash")


#Convert fields from character class (default) to factor class
mash_table$mash_reference <- as.factor(mash_table$mash_reference)
mash_table$mash_query <- as.factor(mash_table$mash_query)
mash_table$mash_ref_query <- as.factor(mash_table$mash_ref_query)




#Import host taxonomy data and phage metadata
host_table <- read.csv("phage_host_data.csv",sep=",",header=TRUE)


#Convert all Unspecified fields to NA missing value
host_table[host_table == "Unspecified"] <- NA



#Modify column names of host data and merge with kmer table
#Data for all phages in processed mash table should be present in host and phage metadata table.
#As a result, do not select all.x=TRUE option. This way, any missing rows indicates an error in the table.


#Original host table column header names:
#"phage_identifier","header_source_info","database","host_superkingdom","host_phylum","host_class","host_order",
#"host_family","host_genus","phage_superkingdom","phage_viral_type","phage_order","phage_family","phage_genus",
#"phage_cluster","phage_subcluster","phage_cluster_source","phage_temperate","size" 


#Match metadata for both the query and reference phages in each pairwise comparison
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
ani_data <- read.csv("ani_data.csv",sep=",",header=TRUE)
ani_data$ani_ref_query <- as.character(ani_data$ani_ref_query)
ani_data$ani_distance <- 1 - ani_data$ani_ani

#Merge with the main table
mash_table2 <- merge(mash_table2,ani_data,by.x="mash_ref_query",by.y="ani_ref_query",all.x=TRUE)



#Import pham data and merge with kmer table
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

#To retain all rows, be sure to keep all.x=TRUE, but you don't want to retain all rows = making scatter plots or histograms can cause errors if not all rows have data.
#Omitting all.x, all rows with no matching pham data are removed, so no errors are encountered when making scatterplots
mash_table2 <- merge(mash_table2,pham_table,by.x="mash_ref_query",by.y="pham_phage1_phage2")



#Assign filter status and change mash distance if data is not significant
mash_table2$filter <- ifelse(mash_table2$mash_pvalue < 1e-10 & mash_table2$size_diff_max_percent < 1,TRUE,FALSE)

#At this point, the max mash distance of all filtered comparisons = 0.4903. So set the distance of all comparisons that did not pass the filter = 0.5
mash_table2$modified_mash_distance <- ifelse(mash_table2$filter == TRUE,mash_table2$mash_distance,0.5)





#Assign evolutionary mode to each pairwise comparison.
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










#At this point, all imported data has been processed. From here on out, commands are mostly analysis
#based on subsets of the data






####REVIEW: The section below might not be needed
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

####












###Check Eukaryotic controls
euk_check <- subset(mash_table2,mash_table2$ref_host_superkingdom == "Eukaryota" | mash_table2$query_host_superkingdom == "Eukaryota")
euk_check <- subset(euk_check,euk_check$ref_host_superkingdom == "Bacteria" | euk_check$query_host_superkingdom == "Bacteria")

compute_sector_distribution(euk_check)

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

par(mar=c(4,8,4,4))
plot(archaea_check$modified_mash_distance,archaea_check$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,20,4))
hist(archaea_check$modified_mash_distance,breaks=((range(archaea_check$modified_mash_distance)[2]-range(archaea_check$modified_mash_distance)[1]) * 100),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(0,0.5),ylim=c(0,200),col="black",cex.axis=2)

par(mar=c(4,8,20,4))
hist(archaea_check$pham_pham_dissimilarity,breaks=((range(archaea_check$pham_pham_dissimilarity)[2]-range(archaea_check$pham_pham_dissimilarity)[1]) * 100 +1),xlab=NULL,ylab=NULL,main=NULL,las=1,xlim=c(1,0),col="black",cex.axis=2,ylim=c(0,2e3),yaxt="n")
axis(side=4,pos=0,cex.axis=2)












###Plot by phage type
type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_diff <- subset(mash_table2,mash_table2$phage_viral_type_compare == "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
type_dsRNA <- subset(type,type$phage_viral_type_compare == "dsRNA")
type_ssDNA <- subset(type,type$phage_viral_type_compare == "ssDNA")
type_ssRNA <- subset(type,type$phage_viral_type_compare == "ssRNA")


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

actino785_singleton_one_or_two <- subset(actino785_data_diff,actino785_data_diff$singleton_one_or_two == TRUE)

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



#Phylum then lifestyle
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

















































###Specific clusters

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
cluster <- subset(type_dsDNA,type_dsDNA$phage_cluster_compare != "different")

cluster_F <- subset(cluster,cluster$phage_cluster_compare == "F")
cluster_BD <- subset(cluster,cluster$phage_cluster_compare == "BD")
cluster_AO <- subset(cluster,cluster$phage_cluster_compare == "AO")
cluster_B <- subset(cluster,cluster$phage_cluster_compare == "B")
cluster_BU <- subset(cluster,cluster$phage_cluster_compare == "BU")
cluster_K <- subset(cluster,cluster$phage_cluster_compare == "K")



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
plot(cluster_B$modified_mash_distance,cluster_B$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_BU$modified_mash_distance,cluster_BU$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(cluster_K$modified_mash_distance,cluster_K$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")















###Compare Mash to ANI
#Compare patterns from pham-mash and pham-ani comparisons


ani_mash <- subset(mash_table2,complete.cases(mash_table2$ani_distance))

ani_mash_clusterg <- subset(ani_mash,ani_mash$phage_cluster_compare=="G")
ani_mash_clusterj <- subset(ani_mash,ani_mash$phage_cluster_compare=="J")
ani_mash_clusterl <- subset(ani_mash,ani_mash$phage_cluster_compare=="L")
ani_mash_clustern <- subset(ani_mash,ani_mash$phage_cluster_compare=="N")
ani_mash_clusterf <- subset(ani_mash,ani_mash$phage_cluster_compare=="F")

#Compare ani to mash
par(mar=c(4,8,4,4))
plot(ani_mash$ani_distance,ani_mash$modified_mash_distance,xlim=c(0,0.5),ylim=c(0,0.5),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")



#Colored plot by assigned evolutionary mode 
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














###pham general vs jaccard dissimilarity

type <- subset(mash_table2,mash_table2$phage_viral_type_compare != "different")
type_dsDNA <- subset(type,type$phage_viral_type_compare == "dsDNA")
bacteria_dsDNA <- subset(type_dsDNA,type_dsDNA$host_superkingdom_compare == 'Bacteria')


#Check how correlated pham dissimilarity and jaccard dissimilarity are
par(mar=c(4,8,4,4))
plot(mash_table2$pham_jaccard_dissimilarity,mash_table2$pham_pham_dissimilarity,xlim=c(0,1),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,1,lty=2,lwd=3,col="grey")


#Mash vs Pham plot using jaccard
par(mar=c(4,8,4,4))
plot(bacteria_dsDNA$modified_mash_distance,bacteria_dsDNA$pham_jaccard_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")




































####NOT SURE IF I NEED THIS. IT INCLUDES CODE FOR GENOME SIZE DISPARITY RUNMEAN
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












####I MAY NOT NEED THE CLUSTER-SPECIFIC PROFILES ABOVE


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
plot_cluster_specific_profiles(cluster_actino,"K")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"AO")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"B")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"BU")
dev.off()
plot_cluster_specific_profiles(cluster_actino,"BD")






































###Toxic phage analysis

mash_table2_toxic <- mash_table2
mash_table2_toxic$toxic_one_or_two <- ifelse(mash_table2_toxic$ref_toxic == "yes" | mash_table2_toxic$query_toxic == "yes",TRUE,FALSE)
mash_table2_toxic$toxic_both <- ifelse(mash_table2_toxic$ref_toxic == "yes" & mash_table2_toxic$query_toxic == "yes",TRUE,FALSE)
toxic_one_or_two <- subset(mash_table2_toxic,mash_table2_toxic$toxic_one_or_two == TRUE)

par(mar=c(4,8,4,4))
plot(toxic_one_or_two$modified_mash_distance,toxic_one_or_two$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")


























###ANI vs Pham Dissimilarity data
#Import ANI data from the 79 genome optimization test to show ANI vs Pham Distance
#Data contains complete matrix of 79 x 79 comparisons, including self comparisons and duplicate (reciprocal) comparisons
ani79_data <- read.csv("ani_79_data.csv",sep=",",header=TRUE)

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
gene_specific_mash_data <- read.csv("gene_specific_mash_analysis.csv",sep=",",header=TRUE)

#Original column headers
#"phage1_phage2","phage1","phage2","phage1_phage2","# shared phams",
#"phage1_phage2 gene content dissimilarity (general index)","phage1_phage2 gene content dissimilarity (jaccard index)",
#"phage1_phage2 all genes mash distance","phage1_phage2 all genes mash p-value","phage1_phage2 all genes kmer count",
#"phage1_phage2 shared genes mash distance","phage1_phage2 shared genes mash p-value","phage1_phage2 shared genes kmer count",
#"phage1_phage2 unshared genes mash distance","phage1_phage2 unshared genes mash p-value","phage1_phage2 unshared genes kmer count",
#"phage1 # unshared phams","phage1 # all phams","phage1 # shared genes","phage1 # unshared genes","phage1 # all genes",
#"phage1 - average length of all genes","phage1 average length of shared genes","phage1 average length of unshared genes",
#"phage1 total length of all genes","phage1 total length of shared genes","phage1 total length of unshared genes",
#"phage1 - all genes GC content","phage1 - shared genes GC content","ref - unshared genes GC content","phage2 # unshared phams",
#"phage2 # all phams","phage2 # shared genes","phage2 # unshared genes","phage2 # all genes","phage2 - average length of all genes",
#"phage2 average length of shared genes","phage2 average length of unshared genes","phage2 total length of all genes",
#"phage2 total length of shared genes","phage2 total length of unshared genes","phage2 - all genes GC content",
#"phage2 - shared genes GC content","phage2 - unshared genes GC content"

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














###VOG analysis





vog_table <- read.csv("shared_vog_proportion_data.csv",sep=",",header=TRUE)

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








#Add phylogenetic data if needed
phylogeny_data <- read.csv("phylogeny_data.csv",sep=",",header=TRUE)


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


par(mar=c(4,8,4,4))
plot(lifestyle_predicted_temperate$modified_mash_distance,lifestyle_predicted_temperate$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")

par(mar=c(4,8,4,4))
plot(lifestyle_predicted_lytic$modified_mash_distance,lifestyle_predicted_lytic$pham_pham_dissimilarity,xlim=c(0,0.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1)
abline(0,2,lty=2,lwd=3,col="grey")



















###Phylogeny comparison
#phylogeny_data_round1 <- phylogeny_data
phylogeny_data <- read.csv("phylogeny_data.csv",sep=",",header=TRUE)


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












par(mar=c(4,8,4,4))
plot(phylogeny_hgcf$phylogeny_distance,phylogeny_hgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="blue")

par(mar=c(4,8,4,4))
plot(phylogeny_lgcf$phylogeny_distance,phylogeny_lgcf$pham_pham_dissimilarity,xlim=c(0,2.5),ylim=c(0,1),pch=1,cex=1,cex.axis=2,ann=FALSE,main=NULL,las=1,col="green")



plot(phylogeny_hgcf$modified_mash_distance,phylogeny_hgcf$phylogeny_distance)
plot(phylogeny_lgcf$modified_mash_distance,phylogeny_lgcf$phylogeny_distance)






####CURRENT REVIEW STOP



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

