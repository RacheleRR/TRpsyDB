
to do list

# things to add in version 2 
-TF
-enhancers from https://lcbb.swjtu.edu.cn/EnhancerDB/ maybe as well from fantome 







okey so i was thinking we can start 


1) for rexprt 
for the established variants we can simply use the table that is here https://pmc.ncbi.nlm.nih.gov/articles/PMC10832122/table/Tab1/ and for the infos about range etc of those variants we can go here https://stripy.org/database . 

for the polymorphic once we can get the catalog we can either go this route 
https://www.nature.com/articles/s41467-025-66153-5#data-availability

or https://github.com/broadinstitute/trexplorer-catalog (https://github.com/broadinstitute/trexplorer-catalog/releases/tag/v2.0) 

then i think we could add the data from the qtls like the tandem repeats for the qtls that are here https://wlcb.oit.uci.edu/TRxQTL/ they can be donwloaded quiete easly should do the rexptr on them ? 

2) and then for the encode we should add this data https://screen.wenglab.org/downloads


3) not sure how i should get the Brain chromatin state
4) and for the Pli etc should we get it from gnomad https://gnomad.broadinstitute.org/data#v4-constraint

https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=gnomadPLI


5) for the TSS distance and splicing junction distance in my code i used this ( refflat and approns and exons  )  should i calculate it for all the trs but i think that would be hell or can i just see where the tss is etc and then see where that tr is and put see it this way ? 
 see code below 

#===========================================
#! TSS proximity (1st type of analysis)
#===========================================
#load data
        refflat <- read.delim(file.path(resource_dir, "refFlat.txt"), stringsAsFactors = F, header = F)
        known.expansion <- read.delim(file.path(resource_dir, "UCSC_simple_repeats_hg38.period_lte20.txt"), stringsAsFactors = F, header = F)
        detected.expansion <- read.delim(file.path(combined_count_dir, "EHDn_DBSCAN.combinedCounts.bed"), header = F)

# Create a list for distance analysis
distance_list <- list()

# Add dynamic groups
for (g in group_names) {
  group_data <- filtered_data %>% filter(outlier_label == g)
  gR  <- prepare_granges(group_data, reduce = TRUE)$gr
  distance_list[[g]] <- calculate_tss_distance(gR, refflat, tss_window)
}

# Add the fixed special groups
known_strs <- prepare_granges(known.expansion,reduce =TRUE)$gr
all_expansions <- prepare_granges(detected.expansion,reduce=TRUE)$gr
distance_list$KnownSTRs <- calculate_tss_distance(known_strs, refflat, tss_window)
distance_list$all_expansions <- calculate_tss_distance(all_expansions, refflat, tss_window)

# Add privacy prefix if private
for (name in names(distance_list)) {
  if (privacy == "private") {
    distance_list[[name]]$rarity <- paste0("Private_", name)
  } else {
    distance_list[[name]]$rarity <- name
  }
  distance_list[[name]]$analysis <- "TSS"
}

# # Define comparisons dynamically (all pairs)
# comparison_groups <- c(group_names, "KnownSTRs", "all_expansions")
# comparisons <- combn(comparison_groups, 2, simplify = FALSE)
# #!should i just make it so that it needs to be specified in the configuration ?? 
# #!it either does the automatic comaprisons or it does the one we tell him in the configuration 

# Perform tests
tss_results <- perform_distance_tests(distance_list, comparisons)

# Combine for plotting
tss_combined <- do.call(rbind, distance_list)

# Prefix based on privacy
privacy_label <- ifelse(privacy == "private", "Private ", "")

# TSS title
tss_title <- sprintf("%sTSS Distance Analysis: %s", 
                     privacy_label, 
                     paste(comparison_groups, collapse = " vs "))

# Plot
tss_plot <- plot_distance_density(tss_combined, analysis_type ="TSS", analysis_name= tss_title)

fname <- title_to_filename(tss_title)

write.table(
  tss_results,
  file.path(output_dir, paste0(fname, "_results.tsv")),
  sep = "\t",
  row.names = FALSE
)

ggsave(
  file.path(output_dir, paste0(fname, ".png")),
  tss_plot,
  width = 10, height = 6, dpi = 300
)

cat("TSS analysis results and plot saved:", file.path(output_dir, fname), "\n")


#===========================================
#! splicing junctions proximity (2nd type of analysis)
#===========================================
#load data
  refflat <- read.delim(file.path(resource_dir, "refFlat.txt"), stringsAsFactors = F, header = F)
  known.expansion <- read.delim(file.path(resource_dir, "UCSC_simple_repeats_hg38.period_lte20.txt"), stringsAsFactors = F, header = F)
  introns <- read.delim(file.path(resource_dir, "hg38_intron_refFlat 1.txt"), stringsAsFactors = F)
  appris <- rbind(read.delim(file.path(resource_dir, "appris_data.principal.refseq108.hg38.txt"), stringsAsFactors = F, header = F),
                        read.delim(file.path(resource_dir, "appris_data.principal.refseq109.hg38.txt"), stringsAsFactors = F, header = F))
  exons <- read.delim(file.path(resource_dir, "hg38_exon_refFlat.txt"), stringsAsFactors = F)
  detected.expansion <- read.delim(file.path(combined_count_dir, "EHDn_DBSCAN.combinedCounts.bed"), header = F)


# Create a list for distance analysis
distance_list_SJ <- list()

# Add dynamic groups
for (g in group_names) {
  group_data <- filtered_data %>% filter(outlier_label == g)
  sj  <- prepare_for_splicing(group_data, max_size = sj_window)
  gR  <- prepare_granges(sj,reduce=FALSE)$gr
  distance_list_SJ[[g]] <- getsplicedistance(gR, introns, sj, exons)
}

appris$V3 <- sapply(sapply(appris$V3, strsplit, "\\."), "[", 1)
introns <- introns[introns$isoform %in% appris$V3, ]

all_expansions_sj <- prepare_for_splicing(detected.expansion, max_size = sj_window)
known_strs_sj <- prepare_for_splicing(known.expansion, max_size = sj_window)

known_strs_prep <- prepare_granges(known_strs_sj, reduce = TRUE)
known_strs_gr <- known_strs_prep$gr
known_strs_sj <- known_strs_prep$df  # Use reduced data frame
all_expansions_gr <- prepare_granges(all_expansions_sj, reduce = FALSE)$gr

distance_list_SJ$KnownSTRs <- getsplicedistance(known_strs_gr,  introns, known_strs_sj, exons) 
distance_list_SJ$all_expansions <- getsplicedistance(all_expansions_gr,  introns, all_expansions_sj, exons)

# Add privacy prefix if private
for (name in names(distance_list_SJ)) {
  if (privacy == "private") {
    distance_list_SJ[[name]]$rarity <- paste0("Private_", name)
  } else {
    distance_list_SJ[[name]]$rarity <- name
  }
  distance_list_SJ[[name]]$analysis <- "SJ"
}

6) for finding the TRs associated with psychiatry i will mine pubmed etc for them and then give the articles or infos to you and then create a json or tsv file with those infos ? 

7) probably should create a resource site in the website like they did in here 

https://strchive.org/resources/

8) add constraint from dazar after 
