#! /usr/bin/env R
library("CNOGpro")

hit_files <- Sys.glob("*.hits")
bootstrap_results <- data.frame()
for (hit_file in hit_files) {
    
    cnogpro_experiment <- CNOGpro(hitsfile = hit_file,
                                 gbkfile = "genomic.gbff",
                                 windowlength = 100,
                                 name = hit_file)
    cnogpro_normalized <- normalizeGC(cnogpro_experiment)
    bootstpared <- runBootstrap(cnogpro_experiment, replicates = 1000, quantiles = c(0.025, 0.975))
    df_x <- bootstpared$genes
    df_x$Sample <- hit_file
    bootstrap_results <- rbind(bootstrap_results, df_x)
}
write.csv(bootstrap_results, "gene_cp_nr.csv", row.names = FALSE)