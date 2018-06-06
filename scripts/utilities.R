library('gridExtra')
library('grid')

aggregate_column <- function(X, target, on_which_column, new_name, fun){
    X$target <- X[, target]
    X$on_which_column <- X[, on_which_column]
    tmp <- sapply(split(X$target, X$on_which_column), fun)
    X[, new_name] <- tmp[X$on_which_column]
    X$target <- NA
    X$on_which_column <- NA
    return(X)
}
best_hit_filter <- function(X, target_1='', target_2='', on='', func){
    
    X$on = X[, on]
    
    if(target_1 != ''){
        X$target_1 <- X[, target_1]
        tmp <- sapply(split(X$target_1, X$on), func)
        X <- subset(X, X$target_1 == tmp[X$on])
    }
    if(target_2 != ''){
        X$target_2 <- X[, target_2]
        tmp <- sapply(split(X$target_2, X$on), func)
        X <- subset(X, X$target_2 == tmp[X$on])
    }
    
    tmp <- sapply(split(X$alignment_score, X$on), max)
    return(subset(X, X$alignment_score == tmp[X$on]))
}
calc_component_match <- function(X){
    Y <- X[ , c('contig_component', 'ref_component', 'contig_component_size', 'ref_component_size')]
    Y <- Y[!duplicated(Y), ]
    Y$match_count = 1
    Y <- aggregate_column(Y, 'match_count', 'contig_component', 'contig_component_match_count', sum)
    Y$match_count = 1
    Y <- aggregate_column(Y, 'match_count', 'ref_component', 'ref_component_match_count', sum)
    return(Y)
}
calc_N50 <- function(X){
    X = sort(X, decreasing=T)
    half_total_length = sum(X)/2
    total = 0
    for(i in X){
        total = total + i
        if(total >= half_total_length){
            return(i)
        }
    }
}
calc_xprs_error <- function(X, answer_label='ref_xprs'){
    for(xprs_label in c('kallisto', 'rsem', 'salmon')){
        target_xprs <- X[, paste0('contig_xprs_', xprs_label)]
        answer_xprs <- X[, paste0(answer_label, '_answer')]
        error_label <- paste0('xprs_error_', xprs_label)
        X[, error_label] <- (target_xprs - answer_xprs) / (target_xprs + answer_xprs) * 100 
        X[is.na(X[, error_label]), error_label] <- 0 
        X[, paste0('abs_', error_label)] <- abs(X[, error_label])
    }
    return(X)
}
categorize_length <- function(X, len_threshold=10, score_threshold=90){
    X$length_difference <- (X$contig_length - X$ref_length) / (X$contig_length + X$ref_length) * 100
    condition <- X$length_difference >= -1 * len_threshold & 
                 X$length_difference <= len_threshold &
                 X$recovery >= score_threshold & 
                 X$accuracy >= score_threshold
    X[condition, 'length_label'] <- 'Full-length'
    
    condition <- X$length_difference > len_threshold &
                 X$recovery >= score_threshold
    X[condition, 'length_label'] <- 'Over-extension'
    
    condition <- X$length_difference < -1 * len_threshold &
                 X$accuracy >= score_threshold
    X[condition, 'length_label'] <- 'Incompleteness'
    
    X[is.na(X$length_label), 'length_label'] <- 'Others'
    X$length_label <- factor(X$length_label)
    
    return(X)
}
categorize_match <- function(X, category, len_threshold=10, score_threshold=90){
    if(category == 'full-length'){
        Y <- categorize_match(X, 'usc')
        return(subset(Y, Y$length_label == 'Full-length')) 
    } else if (category == 'over-extension'){
        Y <- categorize_match(X, 'usc')
        return(subset(Y, Y$length_label == 'Over-extension'))
    } else if (category == 'incompleteness'){
        Y <- categorize_match(X, 'usc')
        return(subset(Y, Y$length_label == 'Incompleteness'))
    } else {
        component_match <- calc_component_match(X)
        component_match <- subset(component_match, component_match$ref_component_match_count == 1 & component_match$contig_component_match_count == 1)
        if(category == 'family-collapse'){
            component_match <- subset(component_match, component_match$ref_component_size > 1 & component_match$contig_component_size == 1)
        } else if (category == 'duplication'){
            component_match <- subset(component_match, component_match$ref_component_size == 1 & component_match$contig_component_size > 1)
        } else if (category == 'multiple-alignment'){
            component_match <- subset(component_match, component_match$ref_component_size > 1 & component_match$contig_component_size > 1)
        } else if (category == 'usc'){
            component_match <- subset(component_match, component_match$ref_component_size == 1 & component_match$contig_component_size == 1)
            return(subset(X, X$contig_component %in% component_match$contig_component &
                             X$ref_component %in% component_match$ref_component))
        }
        X <- subset(X, X$contig_component %in% component_match$contig_component &
                       X$ref_component %in% component_match$ref_component)
        valid_match = subset(X, X$length_label == 'Full-length')
        component_match = calc_component_match(valid_match)
        return(subset(X, X$contig_component %in% component_match$contig_component &
                         X$ref_component %in% component_match$ref_component))
    } 
}
fetch_data <- function(X, dataset_name, species_name, assembly_name=''){
    if(assembly_name == ''){
        return(subset(X, X$dataset==dataset_name & X$species==species_name))   
    }
    return(subset(X, X$dataset==dataset_name & X$species==species_name & X$assembly == assembly_name))
}

grid_arrange_shared_legend <- function(..., ncol=length(list(...)), nrow = 1, position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    
    grid.newpage()
    grid.draw(combined)
    invisible(combined)
}
return_boxplot <- function(X, x, y, color, group_color=F, cut=T, breaks=seq(-100, 100, by=10), max_range=c(-100, 100), guide_name='Categories', show_guide=T, x_label='', y_label='', title='', subtitle=''){
    X$x <- X[, x]
    X$y <- X[, y]
    X$color <- factor(X[, color])
    
    if(cut==T){
        X$group <- cut(X$x, breaks=breaks, include.lowest=T, right=F)
        levels(X$group) <- append(levels(X$group), 'outlier')
        X[which(is.na(X$group)), 'group'] <- 'outlier'
        levels(X$group) <- breaks
    } else {
        X$group <- factor(X$x)
        X$color <- factor(X$color)
    }

    colors <- colorRampPalette(c("#0091ff", "#f0650e"))(length(levels(X$color)))
    names(colors) <- levels(X$color)
    
    if(group_color == T){
        colors <- colorRampPalette(c("#0091ff", "#f0650e"))(length(levels(X$group)))
        names(colors) <- levels(X$group)
        X$color <- X$group
    }
    
    figure <- ggplot(X, aes(x=group, y=y, fill=color, color=color)) + xlab(x_label) + ylab(y_label) +
              geom_hline(color='lightpink', yintercept=0, linetype=2) + geom_boxplot(notch=FALSE, outlier.shape=16, outlier.alpha=0.5, alpha=0.5) +
              scale_color_manual(values=colors, guide=F) +
              labs(title=title, subtitle=subtitle) + guides(fill=guides(title.position='top', title.hjust=0.5))+
              theme(axis.line=element_line(colour='black'), panel.grid.major=element_blank(),
                    panel.border=element_blank(), panel.background=element_blank(),
                    axis.text.x=element_text(angle = 45, hjust = 1)) +
              coord_cartesian(ylim=max_range)
    if(show_guide){
        figure <- figure + scale_fill_manual(values=colors, name=guide_name)
    } else {
        figure <- figure + scale_fill_manual(values=colors, name=guide_name, guide=F)
    }
    
    return(figure)
}
return_cormatrix <- function(X, legend_title="Pearson's r", limit=c(0, 1), title='', subtitle='', grid_step_v=1, grid_step_h=1){
    M <- X
    X <- melt(X, na.rm=T)
    figure <- ggplot(X, aes(Var2, Var1, fill=value)) + labs(title=title, subtitle=subtitle) +
              geom_tile(alpha=0.9) + geom_text(aes(Var2, Var1, label=value), color = "white", size = 4) +
              scale_fill_gradient(low="#0091ff", high = "#f0650e", limit=limit, space = "Lab", name=legend_title) +
              guides(fill=guide_colorbar(title.position='top', title.hjust=0.5, ticks=F, raster=F, barheight=0.5, barwidth=15, override.aes=list(alpha=0.5))) +
              theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank(),
                    panel.grid.major=element_line(colour='gainsboro'), panel.border=element_blank(), panel.background=element_blank()) +
              coord_fixed() + 
              geom_vline(xintercept=seq(1, ncol(M), grid_step_v)-0.5, colour="white", size=1.5) +
              geom_hline(yintercept=seq(1, ncol(M), grid_step_h)-0.5, colour="white", size=1.5) + 
              geom_vline(xintercept=seq(1, ncol(M))-0.5, colour="white", size=0.5) +
              geom_hline(yintercept=seq(1, ncol(M))-0.5, colour="white", size=0.5)  
    return(figure)
}
return_scatterplot <- function(X, x='x', y='y', z='', shape=1, correlation=T, legend_title='Relative Error (%)', axis_breaks=c(0, 5, 10, 15, 20), diagonal=T, x_label='log2(Ground Truth TPM + 1)', y_label='log2(Estimated TPM + 1)', max_range=15, title='', subtitle=''){
    
    X$x <- log2(X[, x] + 1)
    X$y <- log2(X[, y] + 1)
    X$z <- X[,z]
    
    if(correlation==T){
        pearsons = round(cor(X$x, X$y), 3)
        spearman = round(cor(X$x, X$y, method='spearman'), 3)
        subtitle = paste0(subtitle, "\nPearsons'r: ", pearsons, ", Spearman's r: ", spearman)
    }
    
    figure <- ggplot(X, aes(x=x, y=y, colour=z))
    figure <- figure + geom_point(alpha=0.5, shape=shape, size=3) + xlab(x_label) + ylab(y_label) +
              labs(title=title, subtitle=subtitle) + geom_rug(alpha=0.5, show.legend=F) +
              scale_colour_gradientn(colours=colorRampPalette(c("#0091ff", "#f0650e"))(100), limits=c(-100, 100), name=legend_title) +
              scale_x_continuous(breaks=axis_breaks) + scale_y_continuous(breaks=axis_breaks) +
              guides(colour=guide_colorbar(title.position='top', title.hjust=0.5, ticks=F, raster=F, barheight=0.5, barwidth=15, override.aes=list(alpha=0.5)),
              shape=guide_legend(title.position='top', title.hjust=0.5)) +
              theme(axis.line=element_line(colour='black'), panel.grid.major=element_line(colour='gainsboro'), 
                    panel.border=element_blank(), panel.background=element_blank())
    if(diagonal==T){
        figure <- figure + coord_cartesian(xlim=c(0, max_range), ylim=c(0, max_range)) + 
                  geom_abline(intercept=0, slope=1, col='lightpink', linetype=1)
    }  
    
    return(figure)
}
return_barplot <- function(X, x, y, color, range=c(0, 1), x_label='', y_label='', title='', subtitle='', legend_title){
    X$x <- factor(X[, x])
    X$y <- X[, y]
    X$color <- X[, color]
    
    colors <- colorRampPalette(c("#0091ff", "#f0650e"))(length(levels(X$color)))
    names(colors) <- levels(X$color)
    
    figure <- ggplot(X, aes(x=x, y=y, fill=color, color=color)) + xlab(x_label) + ylab(y_label) +
              geom_bar(alpha=0.5, stat="identity", position=position_dodge(), width=0.5) + labs(title=title, subtitle=subtitle) +
              scale_fill_manual(values=colors, name=legend_title) + scale_color_manual(values=colors, guide=F)  +
              guides(fill=guides(title.position='top', title.hjust=0.5)) +
              theme(axis.line=element_line(colour='black'), panel.grid.major=element_blank(),
                    panel.border=element_blank(), panel.background=element_blank(),
                    axis.text.x=element_text(angle=45, hjust=1)) +
              coord_cartesian(ylim=range)
    return(figure)
    
}
return_cumplot <- function(X, x, color, x_label='Proportion of Max. Estimated Abundance', y_label='Cumulative Frequency', range=c(0, 1), title='', subtitle=''){
    X$x <- X[, x]
    X$color <- factor(X[, color])
    figure <- ggplot(X, aes(x=x, colour=color)) + stat_ecdf(geom='step', pad=F)  +
              #geom_line(aes(y=1 - ..y..), stat='ecdf') +
              xlab(x_label) + ylab(y_label) + labs(title=title, subtitle=subtitle) +
              theme(axis.line=element_line(colour='black'), panel.grid.major=element_blank(),
                    panel.border=element_blank(), panel.background=element_blank(),
                    axis.text.x=element_text(angle=45, hjust=1)) + 
              coord_cartesian(xlim=range)
    return(figure)
}
