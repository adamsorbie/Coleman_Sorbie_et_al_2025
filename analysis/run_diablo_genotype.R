library(mixOmics)

# data processing 

load("../data/integration/integration_processed.RData")

source('../scripts/DIABLO.R')
source('../scripts/plsda_utils.R')

# # remove tg/wt
# tgwt_samples <- sample.info.av.luminal %>%
#   filter(Genotype == "tg/wt") %>%
#   rownames()
# 
# remove <- c(tgwt_samples)
# #
# luminal_clr <- luminal_clr[,!colnames(luminal_clr) %in% remove]
# luminal_metabolomics <- luminal_metabolomics[,!colnames(luminal_metabolomics) %in% remove]
# mucosal_clr <- mucosal_clr[,!colnames(mucosal_clr) %in% remove]
# sample.info.av.luminal <- sample.info.av.luminal[!sample.info.av.luminal$SampleID %in% remove,]

# ## remove RIKEN
ontologies <- metabolite.info[,"Ontology"] %>% unique()
rem <- c("Unknown", "null")
int_ontologies <- ontologies[!(ontologies %in% rem)]

luminal_metabolomics <- luminal_metabolomics[metabolite.info[rownames(luminal_metabolomics), "Ontology"] %in% int_ontologies,]

# order columns according to luminal metabolomics
luminal_clr <- luminal_clr[,colnames(luminal_metabolomics)]
mucosal_clr <- mucosal_clr[,colnames(luminal_metabolomics)]

# redundant since they are already matching but good to include just as a sanity check 
luminal_samples <- intersect(colnames(luminal_clr), colnames(luminal_metabolomics))
mucosal_samples <- intersect(colnames(mucosal_clr), colnames(luminal_metabolomics))

### replace Phenotype with your target vector ###
# here we change to Genotype to model 5+12 tgtg versus the rest
luminal_phenotype <- sample.info.av.luminal[luminal_samples, 'Genotype']
names(luminal_phenotype) <- luminal_samples

mucosal_phenotype <- sample.info.av.luminal[mucosal_samples, 'Genotype']
names(mucosal_phenotype) <- mucosal_samples

### Luminal - Luminal
combined_data_luminal <- list(Metabolome=t(luminal_metabolomics),
                      ASV=t(luminal_clr))



### All variables not generated in this file are coming from DIABLO.R ###
luminal_splsda_geno <- tune.run.block.splsda(keepXlist=X.list,
                                X=combined_data_luminal,
                                Y=luminal_phenotype,
                                design=design,
                                ncomp=args$ncomp,
                                folds=args$folds,
                                cpus=args$cpus,
                                nrepeat=5,
                                verbose=FALSE)

luminal_perf_geno <- perf(luminal_splsda_geno, folds=args$folds, nrepeat=10)

save(luminal_splsda_geno, luminal_perf_geno, file="../data/integration/luminal_integration_genotype.RData")

### Luminal - Mucosal
combined_data_mucosal <- list(Metabolome=t(luminal_metabolomics),
                         ASV=t(mucosal_clr))

### All variables not generated in this file are coming from DIABLO.R ###
mucosal_splsda_geno <- tune.run.block.splsda(keepXlist=X.list,
                                X=combined_data_mucosal,
                                Y=mucosal_phenotype,
                                design=design,
                                ncomp=args$ncomp,
                                folds=args$folds,
                                cpus=args$cpus,
                                nrepeat=5,
                                verbose=FALSE)

mucosal_perf_geno <- perf(mucosal_splsda_geno, folds=args$folds, nrepeat=10)

save(mucosal_splsda_geno, mucosal_perf_geno, file="../data/integration/mucosal_integration_genotype.RData")
