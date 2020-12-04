## libraries
library(rgl)
library(TCGAbiolinks)
library(TCGA2STAT)
library(caret)
library(ggplot2)
library(MASS)
library(heatmap.plus)
library(reshape2)
library(RColorBrewer)
library(ConsensusClusterPlus)
library(sigclust)
library(pheatmap)
library(tsne)
library(gplots)
library(ggradar)
library(doMC)
registerDoMC(cores = 7)
source("/home/ye/R/publication-theme-ggplot/ggplot_pub_theme.R")

## color
col = c(
    "#386cb0",
    "#fdb462",
    "#7fc97f",
    "#ef3b2c",
    "#662506",
    "#a6cee3",
    "#fb9a99",
    "#984ea3",
    "#ffff33"
)

## working directory
CC_dir <- file.path('~/Documents/data/')
work_dir <- file.path(CC_dir, 'gbm/')
setwd(work_dir)

## GT genes
gt_file <- file.path(CC_dir, 'GT_Genes2.txt')
gt <- read.table(gt_file, header = F, stringsAsFactors = F)[,1]

## get expression data
#gbm <- getTCGA(disease = 'GBM', clinical = T)
#lgg <- getTCGA(disease = 'LGG', clinical = T)
#save(list = c('gbm', 'lgg'), file = 'glioma.rda')
load('glioma.rda')

## subset expression data on the GT genes
gbm_gt <- as.data.frame(t(gbm$dat[rownames(gbm$dat) %in% gt,]))
lgg_gt <- as.data.frame(t(lgg$dat[rownames(lgg$dat) %in% gt,]))
glioma_gt <- rbind(gbm_gt, lgg_gt)
glioma_gt$patient <- sapply(
    rownames(glioma_gt),
    function (x) paste0(strsplit(x, "-")[[1]][1:3], collapse = "-")
)

## get subtype data
gbm_subtype <- TCGAbiolinks::TCGAquery_subtype(tumor = 'GBM')
lgg_subtype <- TCGAbiolinks::TCGAquery_subtype(tumor = 'LGG')

## subset the subtype data on Original.Subtype, ATRX.status and IDH.status
grep_term <- "Original.Subtype|IDH.status|ATRX.status|patient"
glioma_subtype <- rbind(
    gbm_subtype[, grep(grep_term, colnames(gbm_subtype))],
    lgg_subtype[, grep(grep_term, colnames(lgg_subtype))]
)

## merge expression data with subtype data
glioma <- merge(glioma_subtype, glioma_gt)

## preporcess the data for exporation
set.seed(61186)
pp_explore <- preProcess(glioma, method = c('nzv', 'center', 'scale', 'YeoJohnson'))
ppglioma <- predict(pp_explore, glioma)

## explore: PCA
pca_glioma <- prcomp(ppglioma[,-c(1:4)])
idh.idx <- !is.na(glioma$IDH.status)
atrx.idx <- !is.na(glioma$ATRX.status)
original.idx <- !is.na(glioma$Original.Subtype)
IDH.status <- data.frame(
    IDH.status = ifelse(
        ppglioma$Original.Subtype == 'IDHmut-codel',
        'IDHmut-codel',
        ifelse(
            ppglioma$Original.Subtype == 'IDHmut-non-codel',
            'IDHmut-non-codel',
            'IDH-wt'
        )
    )
)
p1 <- ggplot(
    as.data.frame(pca_glioma$x[idh.idx,]),
    aes(
        PC1,
        PC2,
        color = glioma$IDH.status[idh.idx],
        fill = glioma$IDH.status[idh.idx]))+
    geom_point()+
    stat_ellipse(geom="polygon", alpha=0.2, type = "t", linetype = 2)
p1
p2 <- ggplot(
    as.data.frame(pca_glioma$x[atrx.idx,]),
    aes(
        PC1,
        PC2,
        color = glioma$ATRX.status[atrx.idx],
        fill = glioma$ATRX.status[atrx.idx]))+
    geom_point()+
    stat_ellipse(geom="polygon", alpha=0.2, type = "t", linetype = 2)
p2
p3 <- ggplot(
    as.data.frame(pca_glioma$x[original.idx,]),
    aes(
        PC1,
        PC2,
        color = glioma$Original.Subtype[original.idx],
        fill = glioma$Original.Subtype[original.idx]))+
    geom_point()+
    stat_ellipse(geom="polygon", alpha=0.2, type = "t", linetype = 2)
p3
pltdt <- pca_glioma$x[original.idx,]
plot3d(pltdt[,1], pltdt[,2], pltdt[,3], col = col[glioma$Original.Subtype[original.idx]], type = 's', size = 1)
p4 <- ggplot(
    as.data.frame(pca_glioma$x[original.idx,]),
    aes(
        PC1,
        PC2,
        color = IDH.status[original.idx,1],
        fill = IDH.status[original.idx,1]))+
    geom_point()+
    stat_ellipse(geom="polygon", alpha=0.2, type = "t", linetype = 2)+
    theme_Publication()+
    scale_colour_Publication(name="IDH status")+
    scale_fill_Publication(name="IDH status")
p4
ggsave("gbm-pc12.pdf", width = 5, height = 5, device = cairo_pdf)
p5 <- ggplot(
    as.data.frame(pca_glioma$x[original.idx,]),
    aes(
        PC1,
        PC3,
        color = IDH.status[original.idx,1],
        fill = IDH.status[original.idx,1]))+
    geom_point()+
    stat_ellipse(geom="polygon", alpha=0.2, type = "t", linetype = 2)+
    theme_Publication()+
    scale_colour_Publication(name="IDH status")+
    scale_fill_Publication(name="IDH status")
p5
ggsave("gbm-pc13.pdf", width = 5, height = 5, device = cairo_pdf)
p6 <- ggplot(
    as.data.frame(pca_glioma$x[original.idx,]),
    aes(
        PC2,
        PC3,
        color = IDH.status[original.idx,1],
        fill = IDH.status[original.idx,1]))+
    geom_point()+
    stat_ellipse(geom="polygon", alpha=0.2, type = "t", linetype = 2)+
    theme_Publication()+
    scale_colour_Publication(name="IDH status")+
    scale_fill_Publication(name="IDH status")
p6
ggsave("gbm-pc23.pdf", width = 5, height = 5, device = cairo_pdf)
plot3d(
    pltdt[,1], pltdt[,2], pltdt[,3], 
    col = col[IDH.status[original.idx,1]],
    xlab = "PC1 (23.6%)",
    ylab = "PC2 (8.6%)",
    zlab = "PC3 (7.5%)",
    type = 's', 
    size = 1
    #axes = F
)
snapshot3d("pca-3d.png")

## tsne
#tt3 <- tsne(ppglioma[,-c(1:4)], k=3, perplexity = 100)
#plot3d(tt3[,1], tt3[,2], tt3[,3], col = col[IDH.status$IDH.status], type='s', size=0.7, axes=F)

## explore: LD for the original.subtype variable
ld_original <- lda(Original.Subtype ~ ., data = glioma[,-c(1:3)])
ld_projecions <- as.data.frame(predict(ld_original, glioma)$x)
ld_subtype <- cbind(
    data.frame(subtype = glioma$Original.Subtype),
    ld_projecions
)
write.csv(ld_subtype, 'ld_subtype.csv', row.names = T, quote = F)
cc <- data.frame(group=rownames(confmat$byClass), confmat$byClass[,-c(6,8:10)])
ld_idh <- cbind(
    data.frame(idh = glioma$IDH.status),
    ld_projecions
)
write.csv(ld_idh, 'ld_idh.csv', row.names = T, quote = F)

plot3d(
    ld_subtype[]
)

## explore: differential expression: heatmap
idh.colors <- c("black", "gray")
subtype.colors <- palette()
color.mat <- as.matrix(
    data.frame(
        idh = idh.colors[factor(glioma$IDH.status)],
        subtype = subtype.colors[factor(glioma$Original.Subtype)],
        atrx = ifelse(
            ppglioma[,4] == 'IDHmut-codel',
            "red",
            ifelse(
                ppglioma[,4] == 'IDHmut-non-codel',
                "blue",
                "white"
            )
        )
    )
)
d <- dist(x = ppglioma[,-c(1:4)], method = 'max')
hfun <- function (x) {hclust(
    d = dist(x, method = 'euc'),
    method = 'ward.D2'
)}
numeric_mat <- as.matrix(ppglioma[,-c(1:4)])
numeric_mat <- numeric_mat[,order(pca_glioma$rotation[,1])]
numeric_mat[numeric_mat> 2] <- 2
numeric_mat[numeric_mat<  -2] <- -2
pdf("hc-gbm.pdf", width = 8, height = 16)
heatmap.2(
    as.matrix(numeric_mat), 
    Colv = T,
    RowSideColors = col[factor(color.mat[,3])],
    trace = 'none',
    col = colorRampPalette(c(rep("red", 1), "black", rep("green", 1)))(256),
    hclustfun = hfun, 
    labRow = F, key = F
)
dev.off()

## explore: differential expression: box plot
numeric_mat <- as.matrix(ppglioma[,-c(1:4)])
numeric_mat <- numeric_mat[,order(pca_glioma$rotation[,1])]
plotDF <- melt(cbind(IDH.status, numeric_mat))
plotDF <- na.omit(plotDF)
ggplot(plotDF, aes(variable, value, fill = IDH.status))+
    geom_boxplot()+
    theme_Publication()+
    scale_fill_Publication(name="Subtype")+
    theme(axis.text.x = element_text(angle = 90))
ggsave(width = 40, height = 20, file = 'test.pdf')

## expression profile
plotDF <- melt(
    cbind(
        data.frame(patient = ppglioma$patient),
        IDH.status,
        numeric_mat
    ),
    id.vars = c('patient', 'IDH.status')
)
pltDF_mean <- aggregate(value ~ IDH.status + variable, data = plotDF, FUN = mean)
pltDF_sd <- aggregate(value ~ IDH.status + variable, data = plotDF, FUN = sd)
ggplot()+
    #geom_line(aes(x = plotDF$variable, y = plotDF$value, color = plotDF$IDH.status, group = plotDF$patient), alpha = 0.01)+
    geom_errorbar(aes(x = pltDF_mean$variable, ymin = pltDF_mean$value-pltDF_sd$value, ymax = pltDF_mean$value+pltDF_sd$value, color = pltDF_sd$IDH.status), width = 0.3, size = 0.5)+
    geom_line(aes(x = pltDF_mean$variable, y = pltDF_mean$value, color = pltDF_mean$IDH.status, group = pltDF_mean$IDH.status), size = 0.5)+
    geom_point(aes(x = pltDF_mean$variable, y = pltDF_mean$value, color = pltDF_mean$IDH.status, group = pltDF_mean$IDH.status))+
    #facet_grid(IDH.status~.)+
    geom_hline(yintercept = 0, size = 0.2)+
    ylab("Normalized expression")+
    xlab("")+
    theme_Publication()+
    scale_colour_Publication(name="Subtype")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3))
ggsave("gbm-exp.pdf", width = 12, height = 4, device = cairo_pdf)

## expression profile: heatmap
meanProfile <- dcast(pltDF_mean, IDH.status ~ variable)
rownames(meanProfile) <- meanProfile[,1]
colors <- colorRampPalette(c("darkblue", "white", "darkred"))(7)

pdf("gbm-meanexp.pdf", width = 12, height = 2.3)
pheatmap::pheatmap(as.matrix(meanProfile[,-1]), color = colors)
dev.off()

## Survival
## subset expression data on the GT genes
subset_vars <- c('bcr', 'status', 'OS', gt)
gbm_gt <- as.data.frame(gbm$merged.dat[, colnames(gbm$merged.dat) %in% subset_vars])
lgg_gt <- as.data.frame(lgg$merged.dat[, colnames(lgg$merged.dat) %in% subset_vars])
glioma_gt <- rbind(gbm_gt, lgg_gt)
glioma_gt$patient <- sapply(
    glioma_gt$bcr,
    function (x) paste0(strsplit(x, "-")[[1]][1:3], collapse = "-")
)

## get subtype data
gbm_subtype <- TCGAquery_subtype(tumor = 'GBM')
lgg_subtype <- TCGAquery_subtype(tumor = 'LGG')

## subset the subtype data on Original.Subtype, ATRX.status and IDH.status
grep_term <- "Original.Subtype|IDH.status|patient"
glioma_subtype <- rbind(
    gbm_subtype[, grep(grep_term, colnames(gbm_subtype))],
    lgg_subtype[, grep(grep_term, colnames(lgg_subtype))]
)

## merge expression data with subtype data
glioma <- merge(glioma_subtype, glioma_gt)

## survival plotting function based on survival package
## adopted from http://www.liuzlab.org/TCGA2STAT/viz.misc.R
mypamr.plotsurvival <- function (group, survival.time, censoring.status, cols, lty, lwd)
{
    require(survival)
    n.class <- length(unique(group))
    junk <- survfit(Surv(survival.time, censoring.status) ~ as.factor(group))
    junk2 <- coxph(Surv(survival.time, censoring.status) ~ as.factor(group))
    pv <- 1 - pchisq(2 * (junk2$loglik[2] - junk2$loglik[1]),
                     df = n.class - 1)
    plot(junk, col = cols, xlab = "Years", ylab = "Probability of survival", lty=lty, lwd=lwd)
    legend(x=15, y=1, col = cols, lty = lty, lwd=lwd, legend = as.character(1:n.class), box.lwd=1)
    text(0.1 * max(survival.time), 0.25, paste("p =", as.character(round(pv, 4))), cex=0.8)
    print(pv)
    return()
}

## Plot survival data based on the IDH.status
plotdat <- glioma[,c(2,5,6)]
mypamr.plotsurvival(plotdat[,1], plotdat[,3]/365, plotdat[,2], cols = c('blue', 'red'), lty = c(1,2), lwd = c(2,2))

## Plot survival data based on IDH.statu and codeletion
idh.codel <- ifelse(glioma$Original.Subtype == 'IDHmut-codel', 'IDHmut-codel', ifelse(glioma$Original.Subtype == 'IDHmut-non-codel', 'IDHmut-non-codel', 'IDHwt'))
plotdat$idh.codel <- factor(idh.codel)
pdf("gbm-survival.pdf", width = 5, height = 5)
mypamr.plotsurvival(plotdat[,4], plotdat[,3]/365, plotdat[,2], cols = col, lty = c(1,1,1), lwd = c(2,2,2))
dev.off()

## plot survival data based on de novo clustering
d <- glioma[,-c(1:6)]
pp <- preProcess(d, method = c('nzv', 'center', 'scale', 'YeoJohnson'))
ppd <- predict(pp, d)
tppd <- t(ppd)
set.seed(61186)
cc <- ConsensusClusterPlus(d = tppd, maxK = 10, plot = 'pdf', seed = 1, reps = 1000)
denovo2 <- cc[[2]]$consensusClass
denovo3 <- cc[[5]]$consensusClass
pdf("gbm-survival-dn.pdf", width = 5, height = 5)
mypamr.plotsurvival(denovo3, plotdat[,3]/365, plotdat[,2], cols = brewer.pal(9, "Dark2"), lwd = c(2,2,2), lty = rep(1,5))
dev.off()

## heatmap for denovo clusters
# color.mat <- cbind(
#     denovo2 <- brewer.pal(9, 'Dark2')[factor(denovo2)],
#     idh <- c('black', 'white')[glioma$IDH.status],
#     denovo3 <- brewer.pal(9, 'Dark2')[factor(denovo3)],
#     idh.codel <- c('red', 'blue', 'gray')[factor(idh.codel)]
# )
# colnames(color.mat) <- c('denovo2', 'idh', 'denovo3', 'idh.codel')
# d <- dist(x = ppd, method = 'euc')
# hfun <- hclust(
#     d = d,
#     method = 'ward.D2'
# )
# ppd <- ppd[,order(pca_glioma$rotation[,1])]
# heatmap.plus(
#     as.matrix(ppd),
#     RowSideColors = color.mat,
#     hclustfun = function (x){hfun},
#     Colv = NA,
#     col = colorRampPalette(c('red', 'red', 'black', 'green', 'green'))(255)
# )

# Modelling (IDH) ---------------------------------------------------------
glioma$idh.codel <- idh.codel
attach(glioma)
## split the data into training and testing
set.seed(61186)
intrain <- createDataPartition(glioma$IDH.status[!is.na(glioma$IDH.status)], p = 0.7, list = F)
training <- glioma[intrain,]
testing <- glioma[-intrain,]

## preprocessing
pp <- preProcess(training[,-c(1:6,62)], method = c('nzv', 'center', 'scale', 'YeoJohnson'))
pptraining <- predict(pp, training)

## training control
hold_out <- lapply(unique(idh.codel), function(x) grep(x, pptraining))
steps <- seq(1, nrow(pptraining), ceiling(nrow(pptraining)/10))
kfolds <- lapply(
    steps,
    function(x) {
        seq(x, ifelse(x==steps[length(steps)], nrow(pptraining), x+steps[2]-steps[1]-1), 1)
    }
)
indecies <- c(hold_out, kfolds)
ctrl <- trainControl(
    method = 'cv',
    number = 10,
    #summaryFunction = twoClassSummary,
    savePredictions = T,
    classProbs = T,
    indexOut = indecies,
    verboseIter = T
)
rdaGrid <- expand.grid(gamma = 1:4/4, lambda = 0:4/4)

## training the RDA model
rdafit.idh <- train(
    IDH.status ~ .,
    data = pptraining[!is.na(pptraining$IDH.status),-c(1:4, 65)[-2]],
    method = 'rda',
    trControl = ctrl,
    tuneGrid = rdaGrid,
    metirc = 'ROC'
)
confmat <- confusionMatrix(rdafit.idh)
pheatmap(
    confmat$table,
    color = brewer.pal(10, 'Blues'),
    display_numbers = T,
    number_color = 'black',
    number_format = '%.2f',
    cluster_rows = F,
    cluster_cols = F
)

## blind internal testing
pptesting <- predict(pp, testing)
preds <- predict(rdafit.idh, pptesting[,-c(1:4, 65)])
confmat <- confusionMatrix(preds, pptesting$IDH.status)
pheatmap(
    confmat$table,
    color = brewer.pal(10, 'Blues'),
    display_numbers = T,
    number_color = 'black',
    number_format = '%d',
    cluster_rows = F,
    cluster_cols = F
)

# Modelling (IDH + codel) -------------------------------------------------

## split the data into training and testing
set.seed(61186)
glioma$idh.codel <- make.names(glioma$idh.codel)
glioma <- na.omit(glioma)
intrain <- createDataPartition(glioma$idh.codel[!is.na(glioma$idh.codel)], p = 0.7, list = F)
training <- glioma[intrain,]
testing <- glioma[-intrain,]

## preprocessing
pp <- preProcess(training[,-c(1:6,62)], method = c('nzv', 'center', 'scale', 'YeoJohnson'))
pptraining <- predict(pp, training)

## training control
hold_out <- lapply(unique(idh.codel), function(x) grep(x, pptraining))
steps <- seq(1, nrow(pptraining), ceiling(nrow(pptraining)/10))
kfolds <- lapply(
    steps,
    function(x) {
        seq(x, ifelse(x==steps[length(steps)], nrow(glioma), x+steps[2]-steps[1]-1), 1)
    }
)
indecies <- c(hold_out, kfolds)
ctrl <- trainControl(
    method = 'repeatedcv',
    number = 10, repeats = 10,
    #summaryFunction = twoClassSummary,
    savePredictions = T,
    classProbs = T,
    #indexOut = indecies,
    verboseIter = T
)
rdaGrid <- expand.grid(gamma = 1:4/4, lambda = 0:4/4)

## training the RDA model
rdafit.idh <- train(
    idh.codel ~ .,
    data = pptraining[!is.na(pptraining$idh.codel),-c(1:6)],
    method = 'rda',
    trControl = ctrl,
    tuneGrid = rdaGrid,
    metirc = 'ROC'
)
confmat <- confusionMatrix(rdafit.idh)
pdf("confmat-cv.pdf", width=3, height = 2.5)
pheatmap(
    confmat$table,
    color = brewer.pal(10, 'Blues'),
    display_numbers = T,
    number_color = 'black',
    number_format = '%.2f',
    cluster_rows = F,
    cluster_cols = F
)
dev.off()
tt <- rdafit.idh$pred[rdafit.idh$pred$gamma == 0.25 & rdafit.idh$pred$lambda == 0.75,]
confmat <- confusionMatrix(tt$pred, tt$obs)
confmat
pheatmap(
    confmat$table,
    color = brewer.pal(10, 'Blues'),
    display_numbers = T,
    number_color = 'black',
    number_format = '%.2f',
    cluster_rows = F,
    cluster_cols = F
)
cc <- data.frame(group=rownames(confmat$byClass), confmat$byClass[,-c(6,8:10)])
ggradar(cc)
ggsave("confmat-cv-radar.pdf", width = 12, height = 12, device = cairo_pdf)

write.confmat.tex <- function(df, outfile) {
    write.table(
        round(df, 3),
        file = outfile,
        quote = F,
        sep = "\t&\t",
        eol = "\t\t\\\\\n",
        row.names = T,
        col.names = T
    )
}
write.confmat.tex(confmat$table, "confmat-cv.tex")
write.confmat.tex(confmat$byClass[,c(1:4,7,11)], "confmat-byclass-cv.tex")

## blind internal testing
pptesting <- predict(pp, testing)
preds <- predict(rdafit.idh, pptesting[,-c(1:6, 62)])
confmat <- confusionMatrix(preds, pptesting$idh.codel)
pdf("confmat-testing.pdf", width = 3, height = 2.5)
pheatmap(
    confmat$table,
    color = brewer.pal(10, 'Blues'),
    display_numbers = T,
    number_color = 'black',
    number_format = '%d',
    cluster_rows = F,
    cluster_cols = F
)
dev.off()
cc <- data.frame(group=rownames(confmat$byClass), confmat$byClass[,-c(6,8:10)])
ggradar(cc)
ggsave("confmat-testing-radar.pdf", width = 12, height = 12, device = cairo_pdf)
write.confmat.tex(confmat$table, "confmat-testing.tex")
write.confmat.tex(confmat$byClass[,c(1:4,7,11)], "confmat-byclass-testing.tex")

## feature importance
vi <- filterVarImp(x = pptraining[,-c(1:6,ncol(pptraining))], y = factor(pptraining$idh.codel))
vi <- vi[order(apply(vi,1,sum)),]
vi$gene <- rownames(vi)
vi$gene <- factor(vi$gene, levels = vi$gene)
mvi <- melt(vi, id.vars = "gene") 
ggplot(mvi, aes(gene, value, color=variable))+
    #geom_hline(yintercept = seq(0.5,1.0,0.1), linetype=2, size=0.1)+
    geom_point(shape=15, size=2, alpha=0.8)+
    ylim(c(0.5,1))+
    #coord_flip()+
    theme_Publication()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.3, size=8))+
    xlab("Gene")+
    ylab("AUROC")+
    scale_colour_Publication(name="Subtype")
ggsave("featureImp.pdf", width = 8, height = 3, device = cairo_pdf)

## circular heatmap
## CV
preds.cv <- data.frame(
    subtype = tt$obs,
    tt[,3:5]
)
preds.cv <- preds.cv[order(preds.cv$subtype),]
preds.cv$n <- factor(rownames(preds.cv), levels = as.numeric(rownames(preds.cv)))
m.preds.cv <- melt(preds.cv, id.vars = c("subtype", "n"))
m.preds.cv$var2 <- as.numeric(m.preds.cv$variable) + 5
ggplot()+
    geom_tile(aes(m.preds.cv$n, m.preds.cv$var2, fill=m.preds.cv$value))+
    scale_fill_gradient(low = "grey95", high = "midnightblue", name="Probability")+
    coord_polar(theta = "x")+
    geom_tile(aes(m.preds.cv$n, 5), fill=col[as.numeric(as.factor(m.preds.cv$subtype))])+
    ylim(c(0,max(m.preds.cv$var2)+0.5))+
    theme_classic()+
    theme(
        panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(), 
        axis.line = element_blank()
    )
ggsave("circ-heatmap-cv.pdf", width = 4, height = 4, device = cairo_pdf)

## testing
preds.testing <- data.frame(
    subtype = pptesting$idh.codel,
    predict(rdafit.idh, pptesting, type = "prob")
)
preds.testing <- preds.testing[order(preds.testing$subtype),]
preds.testing$n <- factor(rownames(preds.testing), levels = as.numeric(rownames(preds.testing)))
m.preds.testing <- melt(preds.testing, id.vars = c("subtype", "n"))
m.preds.testing$var2 <- as.numeric(m.preds.testing$variable) + 5
ggplot()+
    geom_tile(aes(m.preds.testing$n, m.preds.testing$var2, fill=m.preds.testing$value))+
    scale_fill_gradient(low = "grey95", high = "midnightblue", name="Probability")+
    coord_polar(theta = "x")+
    geom_tile(aes(m.preds.testing$n, 5), fill=col[as.numeric(as.factor(m.preds.testing$subtype))])+
    ylim(c(0,max(m.preds.testing$var2)+0.5))+
    theme_classic()+
    theme(
        panel.background=element_blank(),
        axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(), 
        axis.line = element_blank()
    )
ggsave("circ-heatmap-testing.pdf", width = 4, height = 4, device = cairo_pdf)

