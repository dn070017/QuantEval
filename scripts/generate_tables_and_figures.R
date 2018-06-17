library('plyr')
library('ggplot2')
library('reshape2')

source('/Users/Hsieh/Google 雲端硬碟/utilities.R')
source('utilities.R')

cat('- setup parameters\n')
{
    setwd('/Users/Hsieh/Downloads/Temporary Files/')
    dir.create('tables', showWarnings=F, recursive=T)
    dir.create('figures', showWarnings=F, recursive=T)
    dataset_list <- c('simulation', 'real')
    dataset_label_list <- c('Simulated', 'Experimental')
    species_list <- c('yeast', 'dog', 'mouse')
    species_label_list <- c('Yeast', 'Dog', 'Mouse')
    assembly_list <- c('rnaspades', 'transabyss', 'trinity')
    assembly_label_list <- c('rnaSPAdes', 'Trans-ABySS', 'Trinity')
    xprs_list <- c('kallisto', 'rsem', 'salmon')
    xprs_label_list <- c('Kallisto', 'RSEM', 'Salmon')
    num_mark <- ','
    max_range <- 15
    score_threshold <- 90
    stage <- 1
    no_table <- 1
    no_figure <- 1
    width <- 1650
    height <- 1500
    res <- 150
}
cat('- stage', stage, 'reading data\n')
{
    ref_dataset <- data.frame()
    contig_dataset <- data.frame()
    match_dataset <- data.frame()
    for(i in seq(1, 2)){
        dataset <- dataset_list[i]
        dataset_label <- dataset_label_list[i]
        for(j in seq(1, 3)){
            species <- species_list[j]
            species_label <- species_label_list[j]
            for(k in seq(1, 3)){
                assembly <- assembly_list[k]
                assembly_label <- assembly_label_list[k]
                
                cat('    - reading', dataset_label, species_label, assembly_label, '\n')
                folder <- paste(dataset, species, assembly, 'QuantEval', sep='/')
                contig <- read.table(paste0(folder, '/contig.tsv'), sep='\t', header=T, stringsAsFactors=F)
                rownames(contig) <- contig$contig_name
                
                contig$dataset <- factor(dataset_label)
                contig$species <- factor(species_label)
                contig$assembly <- factor(assembly_label)
               
                ref <- read.table(paste0(folder, '/ref.tsv'), sep='\t', header=T, stringsAsFactors=F)
                rownames(ref) <- ref$ref_name
                
                ref$dataset <- factor(dataset_label)
                ref$species <- factor(species_label)
                
                match <- read.table(paste0(folder, '/match.tsv'), sep='\t', header=T, stringsAsFactors=F)
                
                match <- subset(match, match$accuracy > score_threshold | match$recovery > score_threshold)
                #match <- calc_xprs_error(match)
                match <- categorize_length(match)
                
                match$alignment_score <- sqrt(match$accuracy * match$recovery)
                
                match$dataset <- factor(dataset_label)
                match$species <- factor(species_label)
                match$assembly <- factor(assembly_label)
                
                contig_dataset <- rbind(contig_dataset, contig)
                match_dataset <- rbind(match_dataset, match)
    
                if(k == 1){
                    ref_dataset <- rbind(ref_dataset, ref)
                }
            }
        }
    }
    write.table(match_dataset, paste0('./tables/match_dataset.tsv'), sep='\t', row.names=F, quote=F)
    write.table(contig_dataset, paste0('./tables/contig_dataset.tsv'), sep='\t', row.names=F, quote=F)
    write.table(ref_dataset, paste0('./tables/ref_dataset.tsv'), sep='\t', row.names=F, quote=F)
    rm(ref, contig, match)
    stage <- stage + 1
}
match_dataset <- read.table('./tables/match_dataset.tsv', sep='\t', header=T, stringsAsFactors=F)
contig_dataset <- read.table('./tables/contig_dataset.tsv', sep='\t', header=T, stringsAsFactors=F)
ref_dataset <- read.table('./tables/ref_dataset.tsv', sep='\t', header=T, stringsAsFactors=F)

cat('- stage', stage, 'generate transcript summary\n')
{
    transcript_summary = data.frame()
    transcript_complexity = data.frame()
    for(species_label in species_label_list){
        ref <- fetch_data(ref_dataset, 'Simulation', species_label)
        gene <- ref[, c('ref_gene', 'ref_gene_size')]
        gene <- gene[!duplicated(gene), ]
        component <- ref[, c('ref_component', 'ref_component_size')]
        component <- component[!duplicated(component), ]
        
        attach(ref)
        X <- data.frame(species=species_label, genes=prettyNum(length(unique(ref_gene)), big.mark=num_mark), 
                        transcripts=prettyNum(nrow(ref), big.mark=num_mark),
                        transcripts_per_gene=format(round(nrow(ref)/length(unique(ref_gene)), 3), nsmall=3),
                        transcript_1000=prettyNum(length(which(ref_length>=1000)), big.mark=num_mark), 
                        transcript_5000=prettyNum(length(which(ref_length>=5000)), big.mark=num_mark),
                        largest_length=prettyNum(max(ref_length), big.mark=num_mark),
                        average_length=prettyNum(format(round(mean(ref_length), 3), nsmall=3), big.mark=num_mark),
                        N50_length=prettyNum(calc_N50(ref_length), big.mark=num_mark), 
                        total_nucleotide=prettyNum(max(total_nucleotide=sum(ref_length), big.mark=num_mark)))
        transcript_summary <- rbind(transcript_summary, X)
        detach(ref)
        X <- data.frame(species=species_label, genes=prettyNum(nrow(gene), big.mark=num_mark),
                        transcripts=prettyNum(nrow(ref), big.mark=num_mark),
                        no_gene_with_single_transcript=prettyNum(length(which(gene$ref_gene_size==1)), big.mark=num_mark),
                        average_no_transcripts_in_gene=prettyNum(format(round(nrow(ref)/nrow(gene), 3), nsmall=3), big.mark=num_mark), 
                        maximum_no_transcripts_in_gene=prettyNum(max(gene$ref_gene_size), big.mark=num_mark),
                        unique_transcript=prettyNum(length(which(ref$ref_component_size==1)), big.mark=num_mark),
                        no_cluster=prettyNum(nrow(component), big.mark=num_mark), 
                        average_no_trasncripts_in_cluster=prettyNum(format(round(mean(component$ref_component_size), 3), nsmall=3), big.mark=num_mark),
                        maximum_no_trasncripts_in_cluster=prettyNum(max(component$ref_component_size), big.mark=num_mark))
        transcript_complexity <- rbind(transcript_complexity, X)
    }
    write.table(transcript_summary, paste0('./tables/table_', no_table, '_transcript_summary.tsv'), sep='\t', row.names=F, quote=F)
    no_table <- no_table + 1
    write.table(transcript_complexity, paste0('./tables/table_', no_table, '_transcript_complexity.tsv'), sep='\t', row.names=F, quote=F)
    no_table <- no_table + 1
    stage <- stage + 1
    rm(species_label, ref, gene, component, X)
}
cat('- stage', stage, 'generate contig summary\n')
{
    figure_list <- list()
    contig_summary <- data.frame()
    contig_ambiguity <- data.frame()
    contig_category <- data.frame()
    alignment_summary <- data.frame()
    transrate <- data.frame()
    for(dataset_label in dataset_label_list){
        for(species_label in species_label_list){
            cat('    - start analyzing', dataset_label, species_label, '\n')
            for(assembly_label in assembly_label_list){
                match <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                contig <- fetch_data(contig_dataset, dataset_label, species_label, assembly_label)
                ref <- fetch_data(ref_dataset, 'Simulation', species_label)
                
                label <- paste(dataset_label, species_label)
                
                component <- contig[, c('contig_component', 'contig_component_size')]
                component <- component[!duplicated(component), ]
                
                full_length <- categorize_match(match, 'full-length')
                over_extension <- categorize_match(match, 'over-extension')
                incompleteness <- categorize_match(match, 'incompleteness')
                family_collapse <- categorize_match(match, 'family-collapse')
                duplication <- categorize_match(match, 'duplication')
                
                contig_best_hit <- best_hit_filter(match, 'accuracy', on='contig_name', func=max)
                transcript_best_hit <- best_hit_filter(match, 'recovery', on='ref_name', func=max)
                correct <- subset(contig_best_hit, contig_best_hit$accuracy >= score_threshold)
                recover <- subset(transcript_best_hit, transcript_best_hit$recovery >= score_threshold)
                
                X <- data.frame(dataset=dataset_label, species=species_label, transcript=prettyNum(nrow(ref), big.mark=num_mark),
                                assembly=assembly_label, contig=prettyNum(nrow(contig), big.mark=num_mark),
                                unique_contig=prettyNum(length(which(contig$contig_component_size==1)), big.mark=num_mark), 
                                no_cluster=prettyNum(nrow(component), big.mark=num_mark),
                                average_no_contigs_in_cluster=prettyNum(format(round(mean(component$contig_component_size), 3), nsmall=3), big.mark=num_mark),
                                maximum_no_contigs_in_cluster=prettyNum(max(component$contig_component_size), big.mark=num_mark),
                                N50_contig_length=prettyNum(calc_N50(contig$contig_length), big.mark=num_mark),
                                aligned_contig=prettyNum(length(unique(match$contig_name)), big.mark=num_mark), 
                                correct_contig=prettyNum(length(unique(correct$contig_name)), big.mark=num_mark), 
                                aligned_transcript=prettyNum(length(unique(match$ref_name)), big.mark=num_mark), 
                                recover_transcript=prettyNum(length(unique(recover$ref_name)), big.mark=num_mark))
                contig_summary <- rbind(contig_summary, X)
                
                X <- data.frame(dataset=dataset_label, species=species_label, contig=prettyNum(nrow(contig), big.mark=num_mark),
                                full_length=prettyNum(length(unique(full_length$contig_name)), big.mark=num_mark), 
                                over_extension=prettyNum(length(unique(over_extension$contig_name)), big.mark=num_mark), 
                                incompleteness=prettyNum(length(unique(incompleteness$contig_name)), big.mark=num_mark), 
                                family_collaspe=prettyNum(length(unique(family_collapse$contig_name)), big.mark=num_mark), 
                                duplication=prettyNum(length(unique(duplication$contig_name)), big.mark=num_mark))
                contig_category <- rbind(contig_category, X)
                
                X <- data.frame(dataset=dataset_label, species=label, assembly=assembly_label, 
                                correct=round(nrow(correct)/nrow(contig), 3), recover=round(nrow(recover)/nrow(ref), 3))
                alignment_summary <- rbind(alignment_summary, X)
                
                X <- data.frame(dataset=dataset_label, species=label, assembly=assembly_label,
                                score=round(median(contig$contig_tr_score), 3),
                                bases_covered=round(median(contig$contig_tr_bases_covered), 3),
                                good_mapping=round(median(contig$contig_tr_good), 3),
                                not_segmented=round(median(contig$contig_tr_not_segmented), 3),
                                seq_true=round(median(contig$contig_tr_seq_true), 3))
                transrate <- rbind(transrate, X)
            }
        }
    }
    figure_list[[1]] <- return_barplot(alignment_summary, 'species', 'correct', 'assembly', y_label='Proportion of Accurate Contigs', legend_title='Assembly')
    figure_list[[2]] <- return_barplot(alignment_summary, 'species', 'recover', 'assembly', y_label='Proportion of Recovered Transcripts', legend_title='Assembly')
    figure_list[[3]] <- return_barplot(transrate, 'species', 'score', 'assembly', y_label='Overall Score', legend_title='Assembly')
    figure_list[[4]] <- return_barplot(transrate, 'species', 'bases_covered', 'assembly', range=c(0.75, 1), y_label='Bases Covered', legend_title='Assembly')
    figure_list[[5]] <- return_barplot(transrate, 'species', 'good_mapping', 'assembly', y_label='Good Mapping', legend_title='Assembly')
    figure_list[[6]] <- return_barplot(transrate, 'species', 'not_segmented', 'assembly', range=c(0.85, 1), y_label='Not Segmented', legend_title='Assembly')
    
    jpeg(paste0('figures/figure_', no_figure, '_proportion_recover_correct_barplot.jpeg'), width=width, height=height/2, res=res)
    grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], ncol=2, nrow=1)
    dev.off()
    no_figure <- no_figure + 1
    jpeg(paste0('figures/figure_', no_figure, '_transrate_score_barplot.jpeg'), width=width, height=height, res=res)
    grid_arrange_shared_legend(figure_list[[3]], figure_list[[4]], figure_list[[5]], figure_list[[6]], ncol=2, nrow=2)
    dev.off()
    no_figure <- no_figure + 1
    
    write.table(contig_summary, paste0('./tables/table_', no_table, '_contig_summary.tsv'), sep='\t', row.names=F, quote=F)
    no_table <- no_table + 1
    write.table(contig_category, paste0('./tables/table_', no_table, '_contig_category.tsv'), sep='\t', row.names=F, quote=F)
    no_table <- no_table + 1
    
    rm(ref, match, contig, component, X, contig_best_hit, transcript_best_hit, correct, recover)
    rm(full_length, over_extension, incompleteness, family_collapse, duplication)
    rm(dataset_label, species_label, assembly_label, label, figure_list)
    stage <- stage + 1
}
cat('- stage', stage, 'generate contig correlation matrix\n')
{
    for(metric in c('tpm', 'count')){
        for(method in c('pearson', 'spearman')){
            for(dataset_label in dataset_label_list){
                i <- 1
                figure_list = list()
                for(species_label in species_label_list){
                    for(assembly_label in assembly_label_list){
                        if(method == 'pearson'){
                            legend_title <- "Pearson's r"
                        } else {
                            legend_title <- "Spearman's r"
                        }
                        
                        if(assembly_label == 'rnaSPAdes'){
                            title = species_label
                        } else {
                            title = ''
                        }
                        
                        contig <- fetch_data(contig_dataset, dataset_label, species_label, assembly_label)
                        cor_matrix <- round(cor(contig[, c(paste0('contig_xprs_', metric, '_kallisto'), paste0('contig_xprs_', metric, '_rsem'), 
                                                           paste0('contig_xprs_', metric, '_salmon'))], method=method), 3)
                        
                        rownames(cor_matrix) <- c('Kallisto', 'RSEM', 'Salmon')
                        colnames(cor_matrix) <- c('Kallisto', 'RSEM', 'Salmon')
                        
                        cor_matrix[lower.tri(cor_matrix)] <- NA
                        
                        figure_list[[i]] <- return_cormatrix(cor_matrix, title=title, subtitle=assembly_label, legend_title=legend_title)
                        i <- i + 1
                    }
                }
                jpeg(paste0('figures/figure_', no_figure, '_contigs_cormatrix_', dataset_label, '_', metric, '_', method, '.jpeg'), width=width, height=height, res=res)
                grid_arrange_shared_legend(figure_list[[1]], figure_list[[4]], figure_list[[7]],
                                           figure_list[[2]], figure_list[[5]], figure_list[[8]],
                                           figure_list[[3]], figure_list[[6]], figure_list[[9]], ncol=3, nrow=3)
                dev.off()
            }
        }
    }
    rm(contig, cor_matrix, i)
    rm(method, dataset_label, species_label, assembly_label, figure_list)
    stage <- stage + 1
    no_figure <- no_figure + 1
}
cat('- stage', stage, 'generate scatter plots\n')
{
    for(metric in c('count', 'tpm')){
        for(category in c('full-length', 'incompleteness', 'over-extension', 'family-collapse-max-xprs', 'family-collapse-tot-xprs', 'family-collapse-max-score', 'duplication-tot-xprs', 'duplication-max-xprs', 'duplication-max-score')){
            for(dataset_label in dataset_label_list){
                cat('    - start analyzing', category, dataset_label, '\n')
                for(species_label in species_label_list){
                    i <- 1
                    figure_list <- list()
                    for(assembly_label in assembly_label_list){
                        for(j in seq(1:3)){
                            answer <- paste0('ref_xprs_', metric, '_answer')
                            xprs_label <- xprs_label_list[j]
                            xprs <- paste0('contig_xprs_', metric, '_', xprs_list[j])
                            xprs_error <- paste0('xprs_', metric, '_error_', xprs_list[j])
                            
                            X <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                            
                            if(category == 'family-collapse-max-xprs'){
                                X <- categorize_match(X, 'family-collapse')
                                X <- best_hit_filter(X, paste0('ref_component_contribute_xprs_', metric, '_answer'), on='contig_name', func=max)
                            } else if (category == 'family-collapse-tot-xprs'){
                                X <- categorize_match(X, 'family-collapse')
                                answer <- paste0('ref_component_tot_xprs_', metric, '_answer')
                                xprs_error <- paste0('tot_component_xprs_', metric, '_error_', xprs_list[j])
                                X[, xprs_error] <- (X[, xprs] - X[, answer]) / (X[, xprs] + X[, answer]) * 100 
                                X[is.na(X[, xprs_error]), xprs_error] <- 0 
                            } else if (category == 'family-collapse-max-score'){
                                X <- categorize_match(X, 'family-collapse')
                                X <- best_hit_filter(X, 'alignment_score', on='contig_name', func=max)
                            } else if (category == 'duplication-max-xprs'){
                                X <- categorize_match(X, 'duplication')
                                X <- best_hit_filter(X, paste0('contig_component_contribute_xprs_', metric, '_', xprs_list[j]), on='ref_name', func=max)
                            } else if (category == 'duplication-tot-xprs'){
                                X <- categorize_match(X, 'duplication')
                                xprs <- paste0('contig_component_tot_xprs_', metric, '_', xprs_list[j])
                                xprs_error <- paste0('tot_component_xprs_', metric, '_error_', xprs_list[j])
                                X[, xprs_error] <- (X[, xprs] - X[, answer]) / (X[, xprs] + X[, answer]) * 100 
                                X[is.na(X[, xprs_error]), xprs_error] <- 0 
                            } else if (category == 'duplication-max-score'){
                                X <- categorize_match(X, 'duplication')
                                X <- best_hit_filter(X, 'alignment_score', on='ref_name', func=max)
                            } else {
                                X <- categorize_match(X, category)
                            }
                            if(metric == 'tpm'){
                                max_range <- 15
                                x_label <- 'log2(Ground Truth TPM + 1)'
                                y_label <- 'log2(Estimated TPM + 1)'
                            } else {
                                max_range <- 22
                                x_label <- 'log2(Ground Truth Count + 1)'
                                y_label <- 'log2(Estimated Count + 1)'
                            }
                            figure_list[[i]] <- return_scatterplot(X, answer, xprs, xprs_error, x_label=x_label, y_label=y_label, max_range=max_range, title=xprs_label, subtitle=assembly_label)
                            i <- i + 1
                        }
                    }
                    jpeg(paste0('figures/figure_', no_figure, '_', category, '_scatterplot_', dataset_label, '_', species_label, '_', metric, '.jpeg'), width=width, height=height, res=res)
                    grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]],
                                               figure_list[[4]], figure_list[[5]], figure_list[[6]],
                                               figure_list[[7]], figure_list[[8]], figure_list[[9]], ncol=3, nrow=3)
                    dev.off()
                }
            }
        }
    }
    rm(answer, xprs_label, xprs, xprs_error, max_range, y_label, X, i, j)
    rm(category, dataset_label, species_label, assembly_label, figure_list)
    stage <- stage + 1
    no_figure <- no_figure + 1
}
cat('- stage', stage, 'generate figure error boxplot (full-length, incompleteness, over-extension)\n')
{
    for(metric in c('tpm', 'count')){
        for(dataset_label in dataset_label_list){
            for(species_label in species_label_list){
                cat('    - start analyzing', dataset_label, species_label, '\n')
                figure_list <- list()
                i <- 1
                for(assembly_label in assembly_label_list){
                    for(j in seq(1:3)){
                        answer <- paste0('ref_xprs_', metric, '_answer')
                        xprs <- paste0('contig_xprs_', metric, '_', xprs_list[j])
                        xprs_error <- paste0('xprs_', metric, '_error_', xprs_list[j])
                        xprs_label <- xprs_label_list[j]
                        
                        X <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                        Y <- rbind(categorize_match(X, 'full-length'), categorize_match(X, 'over-extension'), categorize_match(X, 'incompleteness'))
                        
                        figure_list[[i]] <- return_boxplot(Y, 'length_difference', xprs_error, 'length_label', cut=T, breaks=seq(-100, 40, by=10), x_label='Difference in Length (%)', y_label='Relative Error (%)', title=xprs_label, subtitle=assembly_label)
                        i <- i + 1
                    }
                }
                jpeg(paste0('figures/figure_', no_figure, '_usc_error_boxplot_', dataset_label, '_', species_label, '_', metric, '.jpeg'), width=width, height=height, res=res)
                grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]],
                                           figure_list[[4]], figure_list[[5]], figure_list[[6]],
                                           figure_list[[7]], figure_list[[8]], figure_list[[9]], ncol=3, nrow=3)
                dev.off()
            }
        }
    }
    
    rm(row_index, col_index, col_name, row_name)
    rm(dataset_label, species_label, assembly_label, figure_list)
    rm(answer, xprs_label, xprs, xprs_error, X, Y, i, j)
    no_figure <- no_figure + 2
    stage <- stage + 1
}
cat('- stage', stage, 'generate figure correlation matrix (full-length, incompleteness, over-extension)\n')
{
    for(metric in c('tpm', 'count')){
        i <- 1
        figure_list <- list()
        pdf_list <- list()
        for(dataset_label in dataset_label_list){
            pearson <- matrix(nrow=9, ncol=9)
            spearman <- matrix(nrow=9, ncol=9)
            col_name <- vector(length=9)
            row_name <- vector(length=9)
            row_index <- 1
            for(species_label in species_label_list){
                cat('    - start analyzing', dataset_label, species_label, '\n')
                for(assembly_label in assembly_label_list){
                    row_name[row_index] <- paste(dataset_label, species_label, assembly_label)
                    col_index <- 1
                    for(category in c('full-length', 'incompleteness', 'over-extension')){
                        for(j in seq(1:3)){
                            answer <- paste0('ref_xprs_', metric, '_answer')
                            xprs <- paste0('contig_xprs_', metric, '_', xprs_list[j])
                            xprs_error <- paste0('xprs_', metric, '_error_', xprs_list[j])
                            xprs_label <- xprs_label_list[j]
                            X <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                            X <- categorize_match(X, category)
                            if(category == 'full-length'){
                                category_label <- 'Full-length'
                            } else if (category == 'incompleteness'){
                                category_label <- 'Incompleteness'
                            } else {
                                category_label <- 'Over-extension'
                            }
                            pearson[row_index, col_index] <- round(cor(log2(X[, answer]+1), log2(X[, xprs]+1)), 3)
                            spearman[row_index, col_index] <- round(cor(log2(X[, answer]+1), log2(X[, xprs]+1), method='spearman'), 3)
                            col_name[col_index] <- paste(category_label, xprs_label)
                            col_index <- col_index + 1
                        }
                    }
                    row_index <- row_index + 1
                }
            }
            colnames(pearson) <- col_name
            rownames(pearson) <- row_name
            colnames(spearman) <- col_name
            rownames(spearman) <- row_name
            figure_list[[i]] <- return_cormatrix(pearson, title=dataset_label, subtitle="Pearson's R", legend_title="Correlation Coefficient", grid_step_h=3, grid_step_v=3)
            figure_list[[i + 1]] <- return_cormatrix(spearman, title=dataset_label, subtitle="Spearman's R", legend_title="Correlation Coefficient", grid_step_h=3, grid_step_v=3)
            pdf_list[[i]] <- return_cormatrix(pearson, pdf_font=T, title=dataset_label, subtitle="Pearson's R", legend_title="Correlation Coefficient", grid_step_h=3, grid_step_v=3)
            pdf_list[[i + 1]] <- return_cormatrix(spearman, pdf_font=T, title=dataset_label, subtitle="Spearman's R", legend_title="Correlation Coefficient", grid_step_h=3, grid_step_v=3)
            i <- i + 2
        }
        jpeg(paste0('figures/figure_', no_figure + 1, '_unique_alignment_cormatrix_', metric, '.jpeg'), width=width*1.5, height=height*1.5, res=res)
        grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]], figure_list[[4]], ncol=2, nrow=2)
        dev.off()
        if(metric == 'tpm'){
            ggsave(paste0('./figures/Figure 5.pdf'), grid_arrange_shared_legend(pdf_list[[1]], pdf_list[[2]], pdf_list[[3]], pdf_list[[4]], ncol=2, nrow=2), width=170, height=180, units='mm')
        }
    }
    
    rm(row_index, col_index, col_name, row_name, spearman, pearson)
    rm(dataset_label, species_label, assembly_label, figure_list, pdf_list)
    rm(answer, xprs_label, xprs, xprs_error, X, Y, i, j)
    no_figure <- no_figure + 2
    stage <- stage + 1
}
cat('- stage', stage, 'generate average alignment number for family-collapse and duplication\n')
{
    reflect_xprs = data.frame()
    possible_target = data.frame()
    for(dataset_label in dataset_label_list){
        for(species_label in species_label_list){
            cat('    - start analyzing', dataset_label, species_label, '\n')
            for(assembly_label in assembly_label_list){
                i <- 1
                mean_target <- vector(length = 2)
                for(category in c('family-collapse', 'duplication')){
                    match <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                    match <- categorize_match(match, category)
                    match$match_count <- 1
                    if(category == 'family-collapse'){
                        match <- aggregate_column(match, 'match_count', 'contig_name', 'possible_target', sum)
                        X <- match[, c('contig_name', 'possible_target')]
                        colnames(X) <- c('name', 'possible_target')
                    } else if (category == 'duplication') {
                        match <- aggregate_column(match, 'match_count', 'ref_name', 'possible_target', sum)
                        X <- match[, c('ref_name', 'possible_target')]
                        colnames(X) <- c('name', 'possible_target')
                    }
                    X <- X[!duplicated(X),]
                    mean_target[i] <- round(mean(X$possible_target), 3)
                    i <- i + 1
                }
                possible_target = rbind(possible_target, data.frame(dataset=dataset_label, species=species_label, assembly=assembly_label,
                                                                    family_collapse=mean_target[1], duplication=mean_target[2]))
            }
        }
    }
    write.table(possible_target, paste0('./tables/table_', no_table, '_possible_target.tsv'), sep='\t', row.names=F)
    no_table <- no_table + 1
    stage <- stage + 1
}
cat('- stage', stage, 'generate manuscript figure error boxplot (full-length, incompleteness, over-extension)\n')
{
    figure_list <- list()
    i <- 1
    for(species_label in species_label_list){
        answer <- paste0('ref_xprs_tpm_answer')
        xprs <- 'contig_xprs_tpm_kallisto'
        xprs_error <- 'xprs_tpm_error_kallisto'
        xprs_label <- 'Kallisto'
        
        X <- fetch_data(match_dataset, 'Simulated', species_label, 'Trinity')
        Y <- rbind(categorize_match(X, 'full-length'), categorize_match(X, 'over-extension'), categorize_match(X, 'incompleteness'))
        
        figure_list[[i]] <- return_boxplot(Y, 'length_difference', xprs_error, 'length_label', cut=T, breaks=seq(-100, 40, by=10), x_label='Difference in Length (%)', y_label='Relative Error (%)', pdf_font=T, title=species_label)
        i <- i + 1
    }
    #pdf(paste0('./figures/Figure 3.pdf'), width=6.69, height=2.5)
    ggsave(paste0('./figures/Figure 3.pdf'), grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]], ncol=3, nrow=1), width=170, height=60, units='mm')
    #dev.off()
    
    rm(row_index, col_index, col_name, row_name, spearman, pearson, pearson_figure, spearman_figure)
    rm(dataset_label, species_label, assembly_label, figure_list)
    rm(answer, xprs_label, xprs, xprs_error, X, Y, i, j)
    no_figure <- no_figure + 2
    stage <- stage + 1
}
cat('- stage', stage, 'generate manuscript figure scatterplot boxplot (full-length, incompleteness, over-extension)\n')
{
    i <- 1
    figure_list <- list()
    for(category in c('full-length', 'incompleteness', 'over-extension')){
        answer <- paste0('ref_xprs_tpm_answer')
        xprs_label <- 'Kallisto'
        xprs <- 'contig_xprs_tpm_kallisto'
        xprs_error <- 'xprs_tpm_error_kallisto'
        
        X <- fetch_data(match_dataset, 'Simulated', 'Dog', 'Trinity')
        X <- categorize_match(X, category)

        if(category == 'full-length'){
            category_label <- 'Full-length'
        } else if (category == 'incompleteness'){
            category_label <- 'Incompleteness'
        } else {
            category_label <- 'Over-extension'
        }
        
        max_range <- 15
        x_label <- 'log2(Ground Truth TPM + 1)'
        y_label <- 'log2(Estimated TPM + 1)'
        
        figure_list[[i]] <- return_scatterplot(X, answer, xprs, xprs_error, pdf_font=T, x_label=x_label, y_label=y_label, max_range=max_range, title=category_label)
        i <- i + 1            
    }
    #jpeg(paste0('figures/figure_', no_figure, '_', category, '_scatterplot_', dataset_label, '_', species_label, '_', metric, '.jpeg'), width=width, height=height, res=res)
    ggsave(paste0('./figures/Figure 4.pdf'), grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]], ncol=3, nrow=1), width=170, height=65, units='mm')
    #dev.off()
    rm(answer, xprs_label, xprs, xprs_error, max_range, y_label, X, i, j)
    rm(category, dataset_label, species_label, assembly_label, figure_list)
    stage <- stage + 1
    no_figure <- no_figure + 1
}
cat('- stage', stage, 'generate manuscript figure scatterplot boxplot (famiyl-collapse)\n')
{
    i <- 1
    figure_list <- list()
    for(category in c('Best Alignment Score', 'Highest Expression', 'Others')){
        answer <- paste0('ref_xprs_tpm_answer')
        xprs_label <- 'Kallisto'
        xprs <- 'contig_xprs_tpm_kallisto'
        xprs_error <- 'xprs_tpm_error_kallisto'
        
        X <- fetch_data(match_dataset, 'Simulated', 'Dog', 'Trinity')
        X <- categorize_match(X, 'family-collapse')
        
        max_range <- 15
        x_label <- 'log2(Ground Truth TPM + 1)'
        y_label <- 'log2(Estimated TPM + 1)'
        
        Y <- best_hit_filter(X, paste0('ref_component_contribute_xprs_', metric, '_answer'), on='contig_name', func=max)
        Z <- best_hit_filter(X, 'alignment_score', on='contig_name', func=max)
        
        if (category == 'Highest Expression'){
            X <- subset(X, X$ref_name %in% Y$ref_name)
        } else if (category == 'Best Alignment Score'){
            X <- subset(X, X$ref_name %in% Z$ref_name)
        } else {
            X <- subset(X, ! X$ref_name %in% Y$ref_name & ! X$ref_name %in% Z$ref_name)
        }
        
        figure_list[[i]] <- return_scatterplot(X, answer, xprs, xprs_error, pdf_font=T, x_label=x_label, y_label=y_label, max_range=max_range, title=category)
        i <- i + 1            
    }
    #jpeg(paste0('figures/figure_', no_figure, '_', category, '_scatterplot_', dataset_label, '_', species_label, '_', metric, '.jpeg'), width=width, height=height, res=res)
    ggsave(paste0('./figures/Figure 6.pdf'), grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]], ncol=3, nrow=1), width=170, height=65, units='mm')
    #dev.off()
    rm(answer, xprs_label, xprs, xprs_error, max_range, y_label, X, i, j)
    rm(category, dataset_label, species_label, assembly_label, figure_list)
    stage <- stage + 1
    no_figure <- no_figure + 1
}
cat('- stage', stage, 'generate figure error boxplot and correlation matrix (family-collapse, duplication)\n')
{
    for(metric in c('tpm', 'count')){
        for(category in c('family-collapse')){
            cormatrix_index <- 1
            pearson_figure_list <- list()
            spearman_figure_list <- list()
            for(dataset_label in dataset_label_list){
                error_index <- 1
                error_figure_list <- list()
                pearson <- matrix(nrow=9, ncol=9)
                spearman <- matrix(nrow=9, ncol=9)
                row_name <- vector(length=9)
                col_name <- vector(length=9)
                col_index <- 1
                for(species_label in species_label_list){
                    cat('    - start analyzing', category, dataset_label, species_label, '\n')
                    for(assembly_label in assembly_label_list){
                        contig <- fetch_data(contig_dataset, dataset_label, species_label, assembly_label)
                        ref <- fetch_data(ref_dataset, dataset_label, species_label)
                        error_df <- data.frame()
                        row_index <- 1
                        for(subcategory in c('Maximum Expression and Score', 'Maximum Expression Only', 'Maximum Alignment Score Only', 'Others')){
                             for(i in seq(1:3)){
                                
                                X <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                                X <- categorize_match(X, category)
                                
                                answer <- paste0('ref_xprs_', metric, '_answer')
                                xprs <- paste0('contig_xprs_', metric, '_', xprs_list[i])
                                xprs_error <- paste0('xprs_', metric, '_error_', xprs_list[i])
                                xprs_label <- xprs_label_list[i]
                           
                                if(category == 'family-collapse'){
                                    Y <- best_hit_filter(X, paste0('ref_component_contribute_xprs_', metric, '_answer'), on='contig_name', func=max)
                                    Z <- best_hit_filter(X, 'alignment_score', on='contig_name', func=max)
                                    if (subcategory == 'Maximum Expression and Score'){
                                        X <- subset(X, X$ref_name %in% Y$ref_name & X$ref_name %in% Z$ref_name)
                                    } else if (subcategory == 'Maximum Expression Only'){
                                        X <- subset(X, X$ref_name %in% Y$ref_name & ! X$ref_name %in% Z$ref_name)
                                    } else if (subcategory == 'Maximum Alignment Score Only'){
                                        X <- subset(X, ! X$ref_name %in% Y$ref_name & X$ref_name %in% Z$ref_name)
                                    } else {
                                        X <- subset(X, ! X$ref_name %in% Y$ref_name & ! X$ref_name %in% Z$ref_name)
                                    }
                                } else {
                                    Y <- best_hit_filter(X, paste0('ref_component_contribute_xprs_', metric, '_', xprs_list[i]), on='ref_name', func=max)
                                    Z <- best_hit_filter(X, 'alignment_score', on='ref_name', func=max)
                                    if (subcategory == 'Maximum Expression and Score'){
                                        X <- subset(X, X$contig_name %in% Y$contig_name & X$contig_name %in% Z$contig_name)
                                    } else if (subcategory == 'Maximum Expression Only'){
                                        X <- subset(X, X$contig_name %in% Y$contig_name & ! X$contig_name %in% Z$contig_name)
                                    } else if (subcategory == 'Maximum Alignment Score Only') {
                                        X <- subset(X, ! X$contig_name %in% Y$contig_name & X$contig_name %in% Z$contig_name)
                                    } else {
                                        X <- subset(X, ! X$contig_name %in% Y$contig_name & ! X$contig_name %in% Z$contig_name)
                                    }
                                }
                                if(nrow(X) > 0){
                                    X$error <- X[, xprs_error]
                                    X$xprs_label <- xprs_label
                                    X$subcategory_label <- subcategory
                                    error_df <- rbind(error_df, X[, c('xprs_label', 'error', 'subcategory_label')])
                                    #pearson[row_index, col_index] <- round(cor(log2(X[, answer]+1), log2(X[, xprs]+1)), 3)
                                    #spearman[row_index, col_index] <- round(cor(log2(X[, answer]+1), log2(X[, xprs]+1), method='spearman'), 3)
                                    #row_name[row_index] <- paste(subcategory, xprs_label)
                                    #row_index <- row_index + 1
                                }
                            }
                        }
                        error_df$xprs_label = factor(error_df$xprs_label)
                        error_figure_list[[error_index]] <- return_boxplot(error_df, 'xprs_label', 'error', 'subcategory_label', cut=F, guide_name='Selection Method', x_label='', y_label='Relative Error (%)', title=species_label, subtitle=assembly_label)
                        error_index <- error_index + 1
                        #col_name[col_index] <- paste(species_label, assembly_label)
                        #col_index <- col_index + 1
                    }
                }
                #rownames(pearson) <- row_name
                #colnames(pearson) <- col_name
                #rownames(spearman) <- row_name
                #colnames(spearman) <- col_name
                jpeg(paste0('figures/figure_', no_figure, '_', category, '_error_boxplot_', dataset_label, '_', metric, '.jpeg'), width=width, height=height, res=res)
                grid_arrange_shared_legend(error_figure_list[[1]], error_figure_list[[2]], error_figure_list[[3]],
                                           error_figure_list[[4]], error_figure_list[[5]], error_figure_list[[6]],
                                           error_figure_list[[7]], error_figure_list[[8]], error_figure_list[[9]], ncol=3, nrow=3)
                dev.off()
                #pearson_figure_list[[cormatrix_index]] <- return_cormatrix(pearson, title=dataset_label, subtitle="Pearson's R", legend_title="Correlation Coefficient", grid_step_v=3, grid_step_h=3)
                #spearman_figure_list[[cormatrix_index]] <- return_cormatrix(spearman, title=dataset_label, subtitle="Spearmans's R", legend_title="Correlation Coefficient", grid_step_v=3, grid_step_h=3)
                #cormatrix_index <- cormatrix_index + 1
            }
            #jpeg(paste0('figures/figure_', no_figure + 1, '_', category, '_cormatrix.jpeg'), width=width*1.5, height=height*1.5, res=res)
            #grid_arrange_shared_legend(pearson_figure_list[[1]], pearson_figure_list[[2]], spearman_figure_list[[1]], spearman_figure_list[[2]], ncol=2, nrow=2)
            #dev.off()
        }
    }
    rm(answer, xprs_label, xprs, xprs_error, X, i, contig, ref)
    rm(error_index, cormatrix_index, row_index, col_index, col_name, row_name)
    rm(spearman, pearson, pearson_figure_list, spearman_figure_list)
    rm(dataset_label, species_label, assembly_label)
    no_figure <- no_figure + 2
    stage <- stage + 1
}
cat('- stage', stage, 'generate figure duplication proportion of allocated reads line graph\n')
{
    for(category in c('duplication')){
        for(dataset_label in dataset_label_list){
            error_index <- 1
            error_figure_list <- list()
            for(species_label in species_label_list){
                cat('    - start analyzing', category, dataset_label, species_label, '\n')
                for(assembly_label in assembly_label_list){
                    match <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                    match <- categorize_match(match, category)
                    
                    cumulative <- data.frame()
                    for(i in seq(1:3)){
                        X <- match
                        answer <- 'ref_xprs_tpm_answer'
                        xprs <- paste0('contig_xprs_tpm_', xprs_list[i])
                        xprs_error <- paste0('xprs_tpm_error_', xprs_list[i])
                        xprs_label <- xprs_label_list[i]
                        xprs_contribute <- paste0('contig_component_contribute_xprs_count_', xprs_list[i])
                        
                        X <- aggregate_column(X, xprs_contribute, 'contig_component', 'max_contribute', max)
                        
                        X$xprs_label <- xprs_label
                        X <- X[, c('contig_component', 'max_contribute', 'xprs_label')]
                        
                        cumulative <- rbind(cumulative, X)
                    }
                    error_figure_list[[error_index]] <- return_distribution(cumulative, 'max_contribute', 'xprs_label', title=species_label, subtitle=assembly_label)
                    error_index <- error_index + 1
                }
            }
            jpeg(paste0('figures/figure_', no_figure, '_duplication_maximum_contribution_', dataset_label, '.jpeg'), width=width, height=height, res=res)
            grid_arrange_shared_legend(error_figure_list[[1]], error_figure_list[[2]], error_figure_list[[3]],
                                       error_figure_list[[4]], error_figure_list[[5]], error_figure_list[[6]],
                                       error_figure_list[[7]], error_figure_list[[8]], error_figure_list[[9]], ncol=3, nrow=3)
            dev.off()
        }
    }
    no_figure <- no_figure + 1
    stage <- stage + 1
}
cat('- stage', stage, 'generate manuscript figure duplication error plot\n')
{
    error_index <- 1
    figure_list <- list()
    for(category in c('Best Alignment Score', 'Highest Estimation', 'Total Expression')){
        match <- fetch_data(match_dataset, 'Experimental', 'Mouse', 'rnaSPAdes')
        match <- categorize_match(match, 'duplication')
        
        X <- match
        answer <- 'ref_xprs_tpm_answer'
        xprs <- 'contig_xprs_tpm_kallisto'
        xprs_error <- 'xprs_tpm_error_kallisto'
        xprs_label <- 'Kallisto'
        xprs_contribute <- paste0('contig_component_contribute_xprs_tpm_kallisto')
        X <- aggregate_column(X, xprs_contribute, 'contig_component', 'max_contribute', max)
        
        if(category == 'Best Alignment Score'){
            X <- best_hit_filter(X, 'alignment_score', on='ref_name', func=max)
        } else if (category == 'Highest Estimation'){
            X <- best_hit_filter(X, xprs_contribute, on='ref_name', func=max)
        } else {
            X <- subset(X, X[, xprs] == X[, 'contig_component_max_xprs_tpm_kallisto'])
            xprs <- 'contig_component_tot_xprs_tpm_kallisto'
            X[, xprs_error] <- (X[, xprs] - X[, answer]) / (X[, xprs] + X[, answer]) * 100 
            X[is.na(X[, xprs_error]), xprs_error] <- 0 
        }
        figure_list[[error_index]] <- return_boxplot(X, 'max_contribute', xprs_error, 'assembly', pdf_font=T, group_color=T, cut=T, breaks=seq(0, 1, by=0.05), x_label='Proportion of Maximum Estimated Abundance', y_label='Relative Error (%)', show_guide=F, title=category)
        error_index <- error_index + 1
    }
        
    #jpeg(paste0('figures/figure_', no_figure, '_duplication_error_boxplot_on_contribution_max_xprs_', dataset_label, '_', species_label, '.jpeg'), width=width, height=height, res=res)
    #grid.arrange(error_figure_list[[1]], error_figure_list[[2]], error_figure_list[[3]], ncol=3, nrow=3)
    #dev.off()
    ggsave(paste0('./figures/Figure 7.pdf'), grid.arrange(figure_list[[1]], figure_list[[2]], figure_list[[3]], ncol=3, nrow=1), width=170, height=65, units='mm')
    stage <- stage + 1
}
cat('- stage', stage, 'generate figure component scatter plots\n')
{
    for(category in c('Best Alignment Score', 'Highest Estimation', 'Total Expression')){
        for(dataset_label in dataset_label_list){
            for(species_label in species_label_list){
                i <- 1
                figure_list <- list()
                for(assembly_label in assembly_label_list){
                    for(j in seq(1, 3)){
                        xprs_label <- xprs_label_list[j]
                        answer <- 'ref_component_max_xprs_tpm_answer'
                        xprs_error <- paste0('xprs_tpm_error_', xprs_list[j])
                        match <- fetch_data(match_dataset, dataset_label, species_label, assembly_label)
                        
                        X <- match
                        if(category == 'Best Alignment Score'){
                            xprs <- paste0('contig_xprs_tpm_', xprs_list[j])
                            X <- best_hit_filter(X, 'alignment_score', on='ref_name', func=max)
                            X <- X[, c('ref_component', 'contig_component', answer, xprs, xprs_error)]
                            X <- X[!duplicated(X), ]
                        } else if (category == 'Highest Estimation'){
                            xprs <- paste0('contig_xprs_tpm_', xprs_list[j])
                            X <- best_hit_filter(X, xprs, on='contig_component', func=max)
                            X <- X[, c('ref_component', 'contig_component', answer, xprs, xprs_error)]
                            X <- X[!duplicated(X), ]
                        } else {
                            xprs <- paste0('contig_component_tot_xprs_tpm_', xprs_list[j])
                            X <- X[, c('ref_component', 'contig_component', answer, xprs)]
                            X <- X[!duplicated(X), ]
                            X[, xprs_error] <- (X[, xprs] - X[, answer]) / (X[, xprs] + X[, answer]) * 100
                            X[is.na(X[, xprs_error]), xprs_error] <- 0
                        }
                        figure_list[[i]] <- return_scatterplot(X, answer, xprs, xprs_error, title=assembly_label, subtitle=xprs_label)
                        i <- i + 1
                    }
                }
                jpeg(paste0('figures/figure_', no_figure, '_test2_component_scatterplot_', dataset_label, '_', species_label, '_', category, '.jpeg'), width=width, height=height, res=res)
                grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]],
                                           figure_list[[4]], figure_list[[5]], figure_list[[6]],
                                           figure_list[[7]], figure_list[[8]], figure_list[[9]], ncol=3, nrow=3)
                dev.off()
            }
        }
    }
}
cat('- stage', stage, 'generate manuscript figure component scatter plots\n')
{
    i <- 1
    figure_list <- list()
    for(category in c('ATAC', 'ATHC', 'HTAC', 'HTHC')){
        xprs_label <- 'Kallisto'
        xprs_error <- 'xprs_tpm_error_kallisto'
        match <- fetch_data(match_dataset, 'Simulated', 'Dog', 'Trinity')
        
        X <- match
        if(category == 'ATAC'){
            xprs <- 'contig_xprs_tpm_kallisto'
            answer <- 'ref_xprs_tpm_kallisto'
            X <- best_hit_filter(X, 'alignment_score', on='contig_name', fun=max)
            X <- best_hit_filter(X, 'alignment_score', on='ref_name', func=max)
            #X <- X[, c('ref_component', 'contig_component', answer, xprs, xprs_error)]
            #X <- X[!duplicated(X), ]
        } else if (category == 'ATHC'){
            xprs <- 'contig_xprs_tpm_kallisto'
            answer <- 'ref_xprs_tpm_kallisto'
            X <- best_hit_filter(X, 'alignment_score', on='contig_name', func=max)
            X <- best_hit_filter(X, xprs, on='contig_component', func=max)
            #X <- X[, c('ref_component', 'contig_component', answer, xprs, xprs_error)]
            #X <- X[!duplicated(X), ]
        } else if (category == 'HTAC'){
            xprs <- 'contig_xprs_tpm_kallisto'
            answer <- 'ref_xprs_tpm_kallisto'
            X <- best_hit_filter(X, answer, on='ref_component', func=max)
            X <- best_hit_filter(X, 'alignment_score', on='ref_name', fun=max)
            #X[, xprs_error] <- (X[, xprs] - X[, answer]) / (X[, xprs] + X[, answer]) * 100
            #X[is.na(X[, xprs_error]), xprs_error] <- 0
        } else {
            xprs <- 'contig_xprs_tpm_kallisto'
            answer <- 'ref_xprs_tpm_kallisto'
            X <- best_hit_filter(X, answer, on='ref_component', func=max)
            X <- best_hit_filter(X, xprs, on='contig_component', func=max)
        }
        figure_list[[i]] <- return_scatterplot(X, answer, xprs, xprs_error, pdf_font=T, title=category)
        i <- i + 1
    }
    #jpeg(paste0('figures/figure_', no_figure, '_test2_component_scatterplot_', dataset_label, '_', species_label, '_', category, '.jpeg'), width=width, height=height, res=res)
    #grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]],
    #                           figure_list[[4]], figure_list[[5]], figure_list[[6]],
    #                           figure_list[[7]], figure_list[[8]], figure_list[[9]], ncol=3, nrow=3)
    #dev.off()
    ggsave(paste0('./figures/Figure 8.pdf'), grid_arrange_shared_legend(figure_list[[1]], figure_list[[2]], figure_list[[3]], figure_list[[4]], ncol=2, nrow=2), width=120, height=130, units='mm')
}