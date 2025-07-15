library(tidyverse)
source("../../scripts/full_voclano_plots_unix.R")

load("../../data/ProcessedMetabolomeMicrobiome.RData")
load("../../data/tissue_processed.RData")

volcano_pipeline <- function(
    data, age_vec, group_vec, line_vec,
    age_, groups, line_, metabo_names,
    ontology_data, folder, stage="genotype",
    return_sample_numbers=FALSE,
    x_lim=NULL
) {
    # masks
    geno_mask <- group_vec[colnames(data)] == groups[1] | group_vec[colnames(data)] == groups[2]
    age_mask <- age_vec[colnames(data)] == age_
    line_mask <- line_vec[colnames(data)] == line_
    # subsetting
    # browser()
    sub_data <- data[, geno_mask & age_mask & line_mask]
    
    if (length(unique(group_vec[colnames(sub_data)])) < 2) {
        cat("Less than two sample groups left for",
            paste(groups, age_, line_, sep=", "), "! Returning...\n")
        return(NULL)
    }
    
    sample_numbers <- sapply(
        groups,
        function(g) sum(group_vec[colnames(sub_data)] == g)
    ) 
    names(sample_numbers) <- groups
    
    adj_p <- get_padj(sub_data, group_vec[colnames(sub_data)],
                      test.method=t.test, p.adjust.method="BH")
    fc <- fold.change(data=sub_data,
                      groups=group_vec[colnames(sub_data)],
                      group.order=groups)

    name_plot <- plot_volcano_names(adj_p, fc$fold.changes,
                                    metabo_names, min.p=.01)
    onto_plot <- plot_multi_volcano(adj_p, fc$fold.changes,
                                    ontology_data, min.p=.01,
                                    "Ontology")
    # browser()
    if (!is.null(x_lim)) {
        name_plot <- name_plot + xlim(-x_lim, x_lim)
        # onto_plot <- onto_plot + xlim(-x_lim, x_lim)
    }
    # saving plots
    comp <- str_remove_all(fc$comparison, "/")
    file_ <- paste(line_, age_, comp, sep="_")

    ggsave(filename=paste(folder, paste0(file_, ".pdf"), sep="/"),
           plot=name_plot, device="pdf", width=16, height=9)
    ggsave(filename=paste(folder, paste0(file_, "_by_ontology.pdf"), sep="/"),
           plot=onto_plot$plot, device="pdf", width=16, height=9)

    # export table of: metabolite & logFC & FDR
    name_plot$data$Ontology <- ontology_data[rownames(name_plot$data), "Ontology"]
    name_plot$data %>%
        as_tibble %>%
        dplyr::arrange(desc(log_pvalue), desc(abs(logFC))) %>%
        write_csv(file=paste(folder, paste0(file_, ".csv"), sep="/"))
    
    if (return_sample_numbers) return(sample_numbers)
}

# Preparing data
rownames(fc.scaled) <- rownames(cc.scaled)
rownames(fc.log.normed) <- rownames(cc.log.normed)

total_samples <- c(colnames(fc.scaled), colnames(cc.scaled))

line <- sample.info[total_samples, "Mouse_Line"]
names(line) <- total_samples

age <- sample.info[total_samples,"Age"]
names(age) <- total_samples
age <- factor(age, levels=c("5wk", "12wk", "20wk"))

genotype <- sample.info[total_samples,"Genotype"]
names(genotype) <- total_samples
genotype <- factor(genotype, levels=c("fl/fl", "tg/wt", "tg/tg",
                                      "fl/fl;-/-", "tg/wt;-/-"))

phenotype <- sample.info[total_samples,"Phenotype"]
names(phenotype) <- total_samples
phenotype <- factor(phenotype, levels=c("NT", "T"))

responder <- sample.info[total_samples,"Responder_Status"]
names(responder) <- total_samples

cc.av.mask <- line[colnames(cc.scaled)] == "AV"
cc.avi.mask <- line[colnames(cc.scaled)] == "AVI"
fc.av.mask <- line[colnames(fc.scaled)] == "AV"
fc.avi.mask <- line[colnames(fc.scaled)] == "AVI"


add_features <- metabolite.info[,c("Ontology", "Metabolite name")] %>%
                    add_column(Feature=rownames(metabolite.info))

fn <- metabolite.info$`Metabolite name`
names(fn) <- rownames(metabolite.info)


# tissue processing
tissue_line <- tissue_meta[, "Mouse_Line"]
names(tissue_line) <- rownames(tissue_meta)

tissue_age <- tissue_meta[,"Age"]
names(tissue_age) <- rownames(tissue_meta)
tissue_age <- factor(tissue_age, levels=c("5wk", "12wk", "20wk"))

tissue_genotype <- tissue_meta[,"Genotype"]
names(tissue_genotype) <- rownames(tissue_meta)
tissue_genotype <- factor(tissue_genotype, levels=c("fl/fl", "tg/wt", "tg/tg",
                                                    "fl/fl;-/-", "tg/wt;-/-"))

tissue_phenotype <- tissue_meta[,"Phenotype"]
names(tissue_phenotype) <- rownames(tissue_meta)
tissue_phenotype <- factor(tissue_phenotype, levels=c("NT", "T"))

tissue_status <- tissue_meta[,"Tissue"]
names(tissue_status) <- rownames(tissue_meta)
tissue_status <- factor(tissue_status, levels=c("Healthy", "Susceptible",
                                                "TA", "T"))


tissue_features <- tissue_data$info[,c("Ontology", "Metabolite name")] %>%
                    add_column(Feature=rownames(tissue_data$info))

tissue_fn <- tissue_data$info$`Metabolite name`
names(tissue_fn) <- rownames(tissue_data$info)


#### MAIN ####
if (sys.nframe() == 0) {
    lipids_only <- FALSE
    if (lipids_only) {
        lipid_ontologies <- c(
            "MG", "CAR", "EtherLPE", "EtherDG",
            "DHSph",   "Sph", "LPE", "DG", "SM", "EtherLPC",  
            "EtherPC", "LPC", "SL", "EtherPE", "PI",      
            "BASulfate", "CE", "NAGlySer", "TG", "PE", "STSE", "NAOrn", "FA",
            "LPG", "PA", "LPS", "FAHFA", "PG", "BileAcid", "LPI", "EtherLPG",
            "SSulfate", "PS", "EtherPG", "PC", "EtherSMGDG", "SHex", "PEtOH", "PMeOH"
        )
        lipid_features <- rownames(metabolite.info)[metabolite.info$Ontology %in% lipid_ontologies]
        cc.log.normed <- cc.log.normed[lipid_features,]
        fc.log.normed <- fc.log.normed[lipid_features,]
        tissue_log <- tissue_log[rownames(tissue_data$info)[tissue_data$info$Ontology %in% lipid_ontologies],]
    }
    # ======== #
    #### Q1 ####
    # ======== #
    # caecal
    # pipeline over all ages
    lapply(
        list(c("tg/tg", "fl/fl"), c("tg/wt", "fl/fl")),
        function(comp) {
            lapply(
                unique(age),
                function(x) {
                    volcano_pipeline(
                        cc.log.normed, age, genotype, line,
                        x, comp, "AV",
                        fn, add_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q1/caecal/VolcanoPlots"
                    )
                }
            )
        }
    )

    # faecal
    lapply(
        list(c("tg/tg", "fl/fl"), c("tg/wt", "fl/fl")),
        function(comp) {
            lapply(
                unique(age),
                function(x) {
                    volcano_pipeline(
                        fc.log.normed, age, genotype, line,
                        x, comp, "AV",
                        fn, add_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q1/faecal/VolcanoPlots",
                        x_lim=6.5
                    )
                }
            )
        }
    )

    # tissue
    lapply(
        list(c("tg/tg", "fl/fl"), c("tg/wt", "fl/fl")),
        function(comp) {
            lapply(
                unique(age),
                function(x) {
                    volcano_pipeline(
                        tissue_log, tissue_age, tissue_genotype, tissue_line,
                        x, comp, "AV",
                        tissue_fn, tissue_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q1/tissue/VolcanoPlots"
                    )
                }
            )
        }
    )
    
    # ======== #
    #### Q2 ####
    # ======== #
    phenotype_inferred <- ifelse(grepl("tg/tg", genotype), "T", "NT")
    names(phenotype_inferred) <- names(genotype)
    phenotype_inferred[age[names(phenotype_inferred)] == "5wk" & phenotype_inferred == "T"] <- "Susceptible"
    # caecal
    lapply(
        list(c("T", "NT")),
        function(comp) {
            lapply(
                unique(age),
                function(x) {
                    volcano_pipeline(
                        cc.log.normed, age, phenotype_inferred, line,
                        x, comp, "AV",
                        fn, add_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q2/caecal/VolcanoPlots"
                    )
                }
            )
        }
    )
    
    # faecal
    lapply(
        list(c("T", "NT")),
        function(comp) {
            lapply(
                unique(age),
                function(x) {
                    volcano_pipeline(
                        fc.log.normed, age, phenotype_inferred, line,
                        x, comp, "AV",
                        fn, add_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q2/faecal/VolcanoPlots",
                        x_lim=6.5
                    )
                }
            )
        }
    )
    
    # tissue 
    tissue_pheno_inferred <- ifelse(grepl("tg/tg", tissue_genotype), "T", "NT")
    names(tissue_pheno_inferred) <- names(tissue_genotype)
    tissue_pheno_inferred[age[names(tissue_pheno_inferred)] == "5wk" & tissue_pheno_inferred == "T"] <- "Susceptible"
    lapply(
        list(c("Susceptible", "NT")),
        function(comp) {
            lapply(
                unique(age),
                function(x) {
                    volcano_pipeline(
                        tissue_log, tissue_age, tissue_pheno_inferred, tissue_line,
                        x, comp, "AV",
                        tissue_fn, tissue_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q2/tissue/VolcanoPlots"
                    )
                }
            )
        }
    )
    
    
    # ======== #
    #### Q3 ####
    # ======== #
    time_points <- unique(age)
    comps <- list(c("T", "NT"), c("T", "TA"), c("Healthy", "T"), c("Healthy", "TA"))
    sample_numbers <- lapply(
        comps,
        function(comp) {
            comp_nrs <- sapply(
                time_points,
                function(x) {
                    volcano_pipeline(
                        tissue_log, tissue_age, tissue_status, tissue_line,
                        x, comp, "AV",
                        tissue_fn, tissue_features,
                        "/home/adamsorbie/user/Dropbox/Haller_Lab/Data/Metabolomics/analysis_Nikolai/atf6/analysis/FiguresSandra/Q3/tissue/VolcanoPlots",
                        return_sample_numbers=TRUE
                    )
                }
            )
            names(comp_nrs) <- time_points
            return(comp_nrs)
        }
    )
    names(sample_numbers) <- sapply(comps, paste, collapse="_")
}
