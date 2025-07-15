library(mixOmics)

# data processing 

load("../data/integration/integration_processed.RData")

source('../scripts/DIABLO.R')
source('../scripts/plsda_utils.R')


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
luminal_phenotype <- sample.info.av.luminal[luminal_samples, 'Phenotype']
names(luminal_phenotype) <- luminal_samples

mucosal_phenotype <- sample.info.av.luminal[mucosal_samples, 'Phenotype']
names(mucosal_phenotype) <- mucosal_samples


combined_data <- list(Metabolome=t(luminal_metabolomics),
                      ASV_luminal=t(luminal_clr),
                      ASV_mucosal=t(mucosal_clr))



### All variables not generated in this file are coming from DIABLO.R ###
phenotype_splsda <- tune.run.block.splsda(keepXlist=X.list,
                                X=combined_data,
                                Y=luminal_phenotype,
                                design=design,
                                ncomp=args$ncomp,
                                folds=args$folds,
                                cpus=args$cpus,
                                nrepeat=5,
                                verbose=FALSE)

phenotype_perf <- perf(phenotype_splsda, folds=args$folds, nrepeat=10)

save(phenotype_splsda, phenotype_perf, file="../data/integration/luminal_mucosal_integration_Phenotype.RData")

