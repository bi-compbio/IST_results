## plotting params
## # add arial font, default font sizes
gg_def <- ggplot2::theme_bw() + 
    ggplot2::theme(text = ggplot2::element_text(size = 6, family = "Helvetica"))


# limma doc:
# 
# 9.3 Several Groups
# The above approaches for two groups extend easily to any number of groups. Suppose that three
# RNA targets to be compared. Suppose that the three targets are called “RNA1”, “RNA2” and
# “RNA3” and that the column targets$Target indicates which one was hybridized to each array. An
# appropriate design matrix can be created using
# > f <- factor(targets$Target, levels=c("RNA1","RNA2","RNA3"))
# > design <- model.matrix(~0+f)
# > colnames(design) <- c("RNA1","RNA2","RNA3")
# To make all pair-wise comparisons between the three groups one could proceed
# > fit <- lmFit(eset, design)
# > contrast.matrix <- makeContrasts(RNA2-RNA1, RNA3-RNA2, RNA3-RNA1,
#                                    + levels=design)
# > fit2 <- contrasts.fit(fit, contrast.matrix)
# > fit2 <- eBayes(fit2)
# A list of top genes for RNA2 versus RNA1 can be obtained from
# > topTable(fit2, coef=1, adjust="BH")
# The outcome of each hypothesis test can be assigned using
# > results <- decideTests(fit2)
# A Venn diagram showing numbers of genes significant in each comparison can be obtained from
# > vennDiagram(results)
# The statistic fit2$F and the corresponding fit2$F.p.value combine the three pair-wise comparisons into one F-test. This is equivalent to a one-way ANOVA for each gene except that the residual
# mean squares have been moderated between genes. To find genes which vary between the three RNA
# targets in any way, look for genes with small p-values. To find the top 30 genes:


# obtain KEGG annotations for an organism
get_kegg_pathway_table <- function(org = "hsa", entity = "hsa") {
    message("Retrieving files using KEGGREST...")
    
    # species-specific pathways
    path.name <- KEGGREST::keggList("pathway", org) 
    
    # pathways with names (to have the short pathway name)
    path.name.general <- KEGGREST::keggList("pathway") %>%
        stats::setNames(gsub("map", org, names(.)))
    
    stopifnot(all(names(path.name) %in% names(path.name.general)))
    v.name <- path.name.general[names(path.name)]
    df.name <- data.frame(
        path.id = names(v.name),
        path.name = v.name,
        stringsAsFactors = FALSE)
    
    # entities within pathways
    path.ids <- KEGGREST::keggLink(entity, "pathway")
    
    df.path <- data.frame(
        path.id = gsub("path:[[:alpha:]]+", org, names(path.ids)),
        entity = gsub(paste0(entity, ":"), "", path.ids),
        stringsAsFactors = FALSE
    )
    
    message("Joining tables...")
    
    df.full <- plyr::join(df.name, df.path, by = "path.id", type = "inner")
    
    # db metadata
    v.meta <- capture.output(
        cat(KEGGREST::keggInfo(org), sep = "\n"))
    
    list(
        data = df.full, 
        meta = v.meta
    )
}

# NOTE: designs using this function are limited to two groups
# NOTE: features are in rows
get_DE_table = function(X, y, ...){
    # design
    design = stats::model.matrix(~group, data = data.frame(group = y)) %>%
        magrittr::set_colnames(gsub('group', '', colnames(.)))
    # contrast
    cnt = limma::makeContrasts(
        contrasts=colnames(design)[-1], levels=design)
    # limma chain
    res.limma = X %>%
        limma::lmFit(design=design) %>%
        limma::contrasts.fit(contrasts=cnt) %>%
        limma::eBayes(robust=T)
    # toptable
    limma::topTreat(res.limma, n = Inf, adjust.method = 'fdr') %>%
        tibble::rownames_to_column('ensemblID') %>%
        data.table::data.table()
    # use dplyr to filter for p or fdr
}


# genes that change at all
# y is a factor
# logFC does not make sense here
get_1w_anova_genes <- function(X, y, p.value = 0.05, ...){
    
    design <- stats::model.matrix(~y, data = data.frame(y = y))
    fit <- limma::lmFit(t(X), design)
    cnt <- limma::makeContrasts(contrasts = colnames(design)[-1], levels = design)
    fit.cnt <- limma::contrasts.fit(fit, cnt)
    ebayes <- limma::eBayes(fit.cnt, robust = TRUE)
    
    # # this would be the 1w-ANOVA according to limma's doc
    # # coincides with the topTable without specifying contrasts
    # ans <- data.frame(
    #     gene.id = colnames(X), 
    #     F.stat = ebayes$F, 
    #     P.Value = ebayes$F.p.value
    # ) %>% dplyr::arrange(P.Value)
    
    tt <- limma::topTable(ebayes, number = Inf, adjust.method = "fdr")
    deg <- subset(tt, adj.P.Val < p.value)
    
    list(deg = rownames(deg), all = rownames(tt), topTable = tt)
}

get_up_and_down_genes <- function(X, y, lfc = 0.25, p.value = 0.05, ...){
    
    design <- stats::model.matrix(~y, data = data.frame(y = y))
    fit <- limma::lmFit(t(X), design)
    cnt <- limma::makeContrasts(contrasts = colnames(design)[2], levels = design)
    
    fit.cnt <- limma::contrasts.fit(fit, cnt)
    ebayes <- limma::eBayes(fit.cnt, robust = TRUE)
    
    tt <- limma::topTable(ebayes, number = Inf, adjust.method = "fdr")
    dt <- limma::decideTests(ebayes, lfc = lfc, p.value = p.value, adjust.method = "fdr")
    
    up <- rownames(dt)[dt[, 1] == 1L]
    down <- rownames(dt)[dt[, 1] == -1L]
    all <- rownames(dt)
    
    list(up = up, down = down, deg = c(up, down), all = all, topTable = tt)
}


# ORA from gene list
ora_from_genelist <- function(list.genes, list.bkgd, df.path, min.size = 50, max.size = 500) {
    enr <- clusterProfiler::enricher(
        gene = list.genes, pvalueCutoff = Inf, 
        pAdjustMethod = "fdr", qvalueCutoff = Inf, 
        universe = list.bkgd, TERM2GENE = df.path, 
        minGSSize = min.size, maxGSSize = max.size)
    
    head(enr, Inf)
}


# GSEA from toptable
gsea_from_genevalues <- function(v.genes, df.path, min.size = 50, max.size = 500) {
    v.genes.sort <- sort(v.genes, decreasing = TRUE)
    
    enr <- clusterProfiler::GSEA(
        geneList = v.genes.sort, pvalueCutoff = Inf, 
        pAdjustMethod = "fdr", TERM2GENE = df.path, 
        minGSSize = min.size, maxGSSize = max.size, seed = TRUE)
    
    head(enr, Inf)
}

## Adapted IST functions
save.pheatmap.pdf <- function (x, filename, dev = grDevices::pdf, width = 7, height = 7, 
                               dev.args = list()) 
{
    if (inherits(x, "list")) {
        message("List provided in x. Using x <- x$plot.obj")
        x <- x$plot.obj
    }
    checkmate::assertClass(x, "pheatmap")
    checkmate::qassert(filename, "S1")
    checkmate::qassert(width, "N1")
    checkmate::qassert(height, "N1")
    checkmate::assertClass(dev, "function")
    checkmate::assertClass(dev.args, "list")
    dev.args <- c(dev.args, list(file = filename, width = width, 
                                 height = height))
    do.call(dev, dev.args)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    grDevices::dev.off()
}


save.genemaps.pdf <- function (ist.results, dir.out, id.paths = getPathways(ist.results), 
                               cex.width = 1, cex.height = 1, ...) 
{
    checkmate::qassert(dir.out, "S1")
    checkmate::qassert(id.paths, "S+")
    checkmate::qassert(cex.width, "N1")
    checkmate::qassert(cex.height, "N1")
    if (!dir.exists(dir.out)) 
        dir.create(dir.out, recursive = TRUE)
    list.genemap <- lapply(stats::setNames(id.paths, id.paths), 
                           plot.ist.genemaps, x = ist.results, type = "pheatmap", 
                           ...)
    list.genemap <- list.genemap[!sapply(list.genemap, is.null)]
    if (length(list.genemap) == 0) 
        return(invisible())
    list.size <- lapply(list.genemap, function(x) dim(x$plot.data$data.wide))
    v.pathsize <- sapply(list.size, tail, 1)
    v.sigsize <- sapply(list.size, head, 1)
    v.files <- paste0(dir.out, "/", gsub(" ", "_", names(list.genemap)), 
                      ".pdf")
    width.min <- 7
    height.min <- 3
    size.max <- 50
    v.width <- pmax(pmin(20 * v.pathsize/100 * cex.width, size.max), 
                    width.min)
    v.height <- pmax(pmin(3 + v.sigsize/5 * cex.height, size.max), 
                     height.min)
    mapply(save.pheatmap.pdf, x = list.genemap, filename = v.files, 
           width = v.width, height = v.height)
    invisible(list.genemap)
}


# ## example code
# 
# library(SummarizedExperiment)
# library(magrittr)
# library(limma)
# 
# se.nash = readRDS("1_preprocessing/3_NASH_human_output/sumexp_NASH_hvd.rds")
# map.hsa = read.csv("1_preprocessing/0_genemapping_output/tab_genes_hsa_map.csv")
# df.path = read.csv("3_NASH/1_descriptive_output/df.export.pathways.csv")
# 
# 
# ids = se.nash$fibrosis_stage != -1
# X = assay(se.nash[, ids]) %>% t
# 
# y.fact.multilevel = se.nash[, ids]$fibrosis_stage %>% factor(ordered = FALSE)
# y.fact.twolevel = factor(ifelse(se.nash[, ids]$fibrosis_stage == "0", "Control", "NASH"))
# y.cont = se.nash[, ids]$fibrosis_stage %>% as.character %>% as.numeric
# 
# lfc = 0.25
# p.value = 0.05
# 
# df.symbol2entrez <- unique(na.omit(map.hsa[c("ensembl_gene_id", "hgnc_symbol")]))
# dt.1w <- get_1w_anova_genes(X, y.fact.multilevel)
# 
# length(dt.1w$deg)
# dim(dt.1w$topTable)
# subset(dt.1w$topTable, adj.P.Val < .05) %>% nrow
# 
# # y factor
# dt.updown <- get_up_and_down_genes(X, y.fact.twolevel)
# 
# lapply(dt.updown, length)
# intersect(dt.updown$up, dt.updown$down)
# dim(dt.updown$topTable)
# 
# # y continuous
# dt.cont <- get_up_and_down_genes(X, y.cont)
# 
# lapply(dt.cont, length)
# intersect(dt.cont$up, dt.cont$down)
# dim(dt.cont$topTable)
# 
# ora.up <- ora_from_genelist(dt.cont$up, list.bkgd = dt.cont$all, df.path = df.path)
# dim(ora.up)
# ora.dn <- ora_from_genelist(dt.cont$up, list.bkgd = dt.cont$all, df.path = df.path)
# dim(ora.dn)
# 
# gsea.cont <- gsea_from_genevalues(setNames(dt.cont$topTable$t, rownames(dt.cont$topTable)), df.path = df.path)
# dim(gsea.cont)


summarise_dataset <- function(
        mae, name.se, quote.aes = NULL,
        perc = 0.95, .transform = identity,
        comps = c(1, 2), return.interim = FALSE) {
    
    if (inherits(mae, "MultiAssayExperiment")) {
        message("Summarising a MultiAssayExperiment object")
        
        # transpose to get the usual rows = samples, cols = features
        mat <- do.call(.transform,
                       list(t(MultiAssayExperiment::assays(mae)[[name.se]])))
        
        # for PCA plot
        df.join <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
        join.by <- "primary"
        
        df.coldata <- mae %>%
            MultiAssayExperiment::colData() %>%
            as.data.frame %>%
            dplyr::mutate(primary = rownames(.))
        
    } else if (inherits(mae, "SummarizedExperiment")) {
        message("Summarising a SummarizedExperiment object")
        
        mat <- do.call(.transform, list(t(SummarizedExperiment::assay(mae))))
        
        df.join <- data.frame(colname = rownames(mat))
        join.by <- "colname"
        
        df.coldata <- mae %>%
            SummarizedExperiment::colData() %>%
            as.data.frame %>%
            dplyr::mutate(primary = rownames(.),
                          colname = primary)
    } else {
        stop("'mae' must be of type MultiAssayExperiment or SummarizedExperiment")
    }
    
    message("Dimension")
    show(dim(mat))
    
    message("General properties")
    v <- as.vector(mat)
    print(summary(v))
    
    message("Missings")
    message("- Total")
    print(sum(is.na(v)))
    mat.na <- is.na(mat)
    sample.na <- rowSums(mat.na)
    message("- Number of samples with missings")
    print(sum(sample.na > 0))
    message("- Number of missings per sample")
    print(summary(sample.na))
    message("- Number of missings per feature")
    print(summary(colSums(mat.na)))
    
    message("Summary of the features")
    col.means <- apply(mat, 2, mean)
    col.sds <- apply(mat, 2, sd)
    message("- Means")
    print(summary(col.means))
    message("- Standard deviations")
    print(summary(col.sds))
    
    message("Histograms")
    graphics::par(mfrow = c(1, 2))
    graphics::hist(v,
                   xlab = "Numeric value",
                   main = "Whole matrix")
    
    q3 <- stats::quantile(v, perc, na.rm = TRUE)
    graphics::hist(v[v < q3],
                   xlab = "Numeric value",
                   main = paste0("Values under perc. ", perc))
    
    message("PCA")
    nPcs <- max(comps, 2)
    pca <- pcaMethods::pca(mat, method = "nipals",
                           scale = "uv", center = TRUE,
                           nPcs = nPcs)
    show(pca)
    
    name.comp <- paste0("PC", comps)
    name.compvar <- paste0(name.comp, " (", round(pca@R2[comps]*100, 2), "%)")
    
    if (inherits(mae, "MultiAssayExperiment")) {
        df.join <- as.data.frame(MultiAssayExperiment::sampleMap(mae))
    } else {
        
    }
    
    # df.join is a dummy data.frame with one single column
    # if mae is a SE object, otherwise it contains the mapping to the
    # primary sample IDs in the MAE
    df.pca <- pcaMethods::scores(pca) %>% as.data.frame %>%
        dplyr::mutate(colname = rownames(.)) %>%
        plyr::join(df.join)
    
    # check if PCx columns already exist, overwrite if so
    col.existing <- intersect(colnames(df.coldata), name.comp)
    if (length(col.existing) > 0) {
        message("Overwriting already existing PCx columns: ", col.existing)
        df.coldata <- df.coldata[setdiff(colnames(df.coldata), col.existing)]
    }
    
    df.plot <- plyr::join(df.pca, df.coldata, by = join.by, type = "left")
    
    if (length(comps) == 2) {
        # two component: scatterplot
        gg.obj <- ggplot(df.plot, aes_string(x = name.comp[1], y = name.comp[2])) +
            geom_hline(yintercept = 0, lty = 2) +
            geom_vline(xintercept = 0, lty = 2) +
            geom_point(quote.aes) +
            coord_fixed() +
            xlab(name.compvar[1]) +
            ylab(name.compvar[2])
    } else if (length(comps) > 2) {
        # more than two: pairs plot
        gg.obj <- GGally::ggpairs(
            df.plot,
            mapping = quote.aes,
            columns = name.comp,
            columnLabels = name.compvar)
    } else {
        message("ncomps must have length 2 or more to plot the PCA components")
    }
    
    if (return.interim) {
        list(gg.obj = gg.obj, df.plot = df.plot, pca = pca)
    } else {
        gg.obj
    }
}


# upper triangle of matrix to dataframe
getUpper = function(m, value.name = 'value'){
    rowCol = expand.grid(rownames(m), colnames(m))
    df = rowCol[as.vector(upper.tri(m,diag=F)),]
    df[, value.name] = m[upper.tri(m,diag=F)]
    data.table(df)
}

#' KB - #' Compare List of Genesets
#'
#' This function uses simple hypergeomtric tests to determine the
#' overlap between sets in a list of sets.
#'  
#' @param geneSetsA  list of gene sets A
#' @param geneSetsB  list of gene sets B (replaced by A if NULL)
#' @param bkgSet set of background genes
#' 
#' @return object$setInfo, object$setTable, object$matrix$...
#' 
compareSets.new = function(geneSetsA, geneSetsB=NULL, bkgSet, setInfoIn = NULL, minSize=10){
    # merge genesets
    geneSets = c(geneSetsA, geneSetsB)
    # remove all genes not in background
    geneSets = lapply(geneSets, function(geneSet)
        intersect(geneSet, bkgSet))
    # remove sets with size < minSize
    geneSets = geneSets[!(lapply(geneSets, length) < minSize)]
    # statistics on individual gene sets
    setInfoOut = data.table(
        setName=names(geneSets),
        setSize=sapply(geneSets, function(x) length(x))
    )
    # merge setInfoOut with setInfoIn
    if (!is.null(setInfoIn)){
        setInfoOut = merge(setInfoIn, setInfoOut, by='setName')
    }
    # if no second list provided, test all in A against A
    geneSetsA = geneSets[intersect(names(geneSetsA), names(geneSets))]
    if (is.null(geneSetsB)){
        geneSetsB = geneSetsA
    } else {
        geneSetsB = geneSets[intersect(names(geneSetsB), names(geneSets))]
    }
    # prepare matrices
    mat = matrix(nrow=length(geneSetsA), ncol=length(geneSetsB))
    rownames(mat) = names(geneSetsA)
    colnames(mat) = names(geneSetsB)
    mat.p = mat.n = mat.exp = mat.fold = mat
    # comparison of genesets
    # TODO: parallelize
    # TODO: only calculate upper triangle in case A = B
    for (i in 1:length(geneSetsA)){
        for (j in 1:length(geneSetsB)){
            nameA = names(geneSetsA)[i]
            nameB = names(geneSetsB)[j]
            setA = geneSetsA[[i]]
            setB = geneSetsB[[j]]
            setAB = intersect(setA, setB)
            p = phyper(
                length(setAB), 
                length(setA), 
                length(bkgSet) - length(setA), 
                length(setB), lower.tail = F)
            # write matrices
            mat.p[i, j] = p
            mat.n[i, j] = length(setAB)
            mat.exp[i, j] = length(setA) / length(bkgSet) * length(setB)
            mat.fold = mat.n / mat.exp
        }
    }
    # combine matrices to data table
    df.p = getUpper(mat.p, value.name = 'pVal')
    df.n = getUpper(mat.n, value.name = 'observed')
    df.exp = getUpper(mat.exp, value.name = 'expected')
    df.fold = getUpper(mat.fold, value.name = 'fold')
    setTable = Reduce(function(x, y) merge(x, y, by=c('Var1', 'Var2')), list(df.exp, df.n, df.fold, df.p))
    setTable$setSize1 = setInfoOut[match(setTable$Var1, setInfoOut$setName), setSize]
    setTable$setSize2 = setInfoOut[match(setTable$Var2, setInfoOut$setName), setSize]
    # multiple testing
    setTable$fdr = p.adjust(setTable$pVal, method = 'BH')
    # return
    out = list()
    out$setInfo = setInfoOut
    out$setTable = setTable
    out$matrix$pVal = mat.p
    out$matrix$observed = mat.n
    out$matrix$expected = mat.exp
    out$matrix$fold = mat.fold
    out
}