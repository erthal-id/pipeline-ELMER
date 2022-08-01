#Baixar os dados de metilação das duas amostras
query.met <- GDCquery(project = "TCGA-LIHC",
                  data.category = "DNA Methylation",
                  legacy = FALSE,
                  platform = c("Illumina Human Methylation 450"),
                  sample.type = c("Primary Tumor", "Solid Tissue Normal") 
)

#Baixar os dados de expressão gênica das duas amostras 
#Nesta função, legacy = FALSE também (é default)
query.exp <- GDCquery(project = "TCGA-GBM",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts",
)



#Tabela com os dados do query.exp
#N de amostras -> Primary Tumor: 371
#                 Solid Tissue Normal: 50
#                 Total: 421
datatable(getResults(query.exp, cols = c("data_type","cases","sample_type")), 
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE) 

GDCdownload(query.met, method = "api", files.per.chunk = 10)
GDCdownload(query.exp), method = "api", files.per.chunk = 10)
data.meth <- GDCprepare(query.met)
data.exp <- GDCprepare(query.exp)
assay(data.meth)[1:5,1:5]
assay(data.exp)[1:5,1:5]

#Agora que eu já tenho os dados de metilação e expressão gênica,
#Baixo os distal probes para fazer o MAE
distal.probes <- get.feature.probe(genome = "hg38",
                                   met.platform = "450K",
                                   rm.chr = paste0("chr", c("X", "Y")))
#Antes "rm.chr = paste0("chr", c(2:22, "X", "Y")))" -> 15142 distal probes
#Depois "rm.chr = paste0("chr", c("X", "Y")))" -> 158803 distal probes

#Criando o MultiAssayExperiment (MAE) 
LIHC.mae <- createMAE(exp = data.exp, 
                      met = data.meth,
                      filter.probes = distal.probes,
                      linearize.exp = TRUE,
                      met.platform = "450K",
                      genome = "hg38",
                      save = TRUE,
                      TCGA = TRUE) 

group1 <- "Primary solid Tumor"  #Grupo1: tumorais
group2 <- "Solid Tissue Normal"  #Grupo2: normais

############################ Hypo analysis ###################################


#Volcano plot 
TCGAVisualize_volcano(x = as.data.frame(sig.diff.hypo)[,grep("Minus",colnames(sig.diff.hypo),value = T)],
                      y = sig.diff.hypo$adjust.p,
                      title = paste0("Volcano plot - Probes ", "hypomethylated in ", group1, " vs ", group2, "\n"),
                      label =  c("Not Significant",
                                 paste0("Hypermethylated in ",group1),
                                 paste0("Hypomethylated in ",group1)),
                      ylab =  expression(paste(-Log[10],
                                               " (FDR corrected P-values) [one tailed test]")),
                      xlab =  expression(paste(
                        "DNA Methylation difference (",beta,"-values)")
                      ),
                      x.cut = 0.4, 
                      y.cut = 0.01
)

#Busca os 20 genes mais proximos à probe, 10 acima e 10 abaixo
neargenes.hypo <- GetNearGenes(data=LIHC.mae, probes = sig.diff.hypo$probe, numFlankingGenes = 20)
#total: 952138 genes

Hypo.pair <- get.pair(data = LIHC.mae, 
                      group.col = "definition", 
                      group1 = group1,
                      group2 = group2, 
                      nearGenes = neargenes, 
                      mode = "unsupervised",
                      permu.size = 35482, 
                      raw.pvalue = 0.001,
                      Pe = 0.001,
                      filter.probes = TRUE, 
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      label = "hypo", 
                      cores = 1)
#Para pvalue 0.01 e sig.dif 0.4 -> 35482: 765, 574 não se repetem
#length(unique(Hypo.pair$Probe))

for (i in 10:15) {
  scatter.plot(data = LIHC.mae,
               byPair = list(probe = Hypo.pair$Probe[i], gene = Hypo.pair$GeneID[i]), 
               category = "definition", save = TRUE, lm_line = TRUE)
}

enriched.motif <- get.enriched.motif(data = LIHC.mae,
                                     probes = Hypo.pair$Probe,
                                     label = "hypo",
                                     min.incidence = 1,
                                     lower.OR = 1.1)
#Para pvalue 0.01 e sig.dif 0.4 -> 16 enriched motifs

motif.enrichment.plot(motif.enrichment = "/media/hd/isabela/getMotif.hypo.motif.enrichment.hypo.csv", 
                      significant = list(OR = 1.5,lowerOR = 1.3), 
                      label = "hypo", 
                      summary = TRUE,
                      save = FALSE)

TF.hypo <- get.TFs(data = LIHC.mae, mode = "unsupervised", group.col = "definition", 
              group1 = group1, group2 = group2, enriched.motif = enriched.motif.hypo,
              cores = 1, label = "hypo")

load("/media/hd/isabela/getTF.hypo.TFs.with.motif.pvalue.rda")

motif1 <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif1,
             save = TRUE)

motif2 <- colnames(TF.meth.cor)[2]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif2,
             save = TRUE) 

motif16 <- colnames(TF.meth.cor)[16]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif16,
             save = TRUE) 

scatter.plot(data = LIHC.mae,
             byTF = list(TF = c("ZNF716", "PRDM9", "CENPX"),
                         probe = c(enriched.motif$HNF1B_HUMAN.H11MO.0.A)), 
             category = "definition",
             save = TRUE, 
             lm_line = TRUE)

scatter.plot(data = LIHC.mae,
             byTF = list(TF = c("CENPX", "ZNF716", "PRDM9"),
                         probe = c(enriched.motif$ZFHX3_HUMAN.H11MO.0.D)), 
             category = "definition",
             save = TRUE, 
             lm_line = TRUE)

scatter.plot(data = LIHC.mae,
             byTF = list(TF = c("MXD3", "ZNF716", "MYBL2" ),
                         probe = c(enriched.motif$PO5F1_HUMAN.H11MO.1.A)), 
             category = "definition",
             save = TRUE, 
             lm_line = TRUE) 

for (i in 1:2) {
schematic.plot(pair = Hypo.pair, 
               data = LIHC.mae,
               group.col = "definition",
               byProbe = Hypo.pair$Probe[2],
               save = FALSE)
}

heatmapPairs(data = LIHC.mae, 
             group.col = "definition",
             group1 = "Primary solid Tumor",
             group2 = "Solid Tissue Normal",
             pairs = Hypo.pair,
             filename =  NULL)

#####################################################################################

# Filtrando as probes em apenas as que apareceram na tabela 'enriched.motif'
sig.diff.em <- read.xlsx("file_show.xlsx") 
#le a tabela xlsx com as probes e a armazena em sig.diff.em

neargenes.em <- GetNearGenes(data=LIHC.mae, probes = sig.diff.em$probe,
                             numFlankingGenes = 20)
#busca os 20 genes próximos (10 US e 10 DS) das probes da nova tabela

Hypo.pair.em <- get.pair(data = LIHC.mae, 
                      group.col = "definition", 
                      group1 = group1,
                      group2 = group2, 
                      nearGenes = neargenes.em, 
                      mode = "unsupervised",
                      permu.size = 5995, 
                      raw.pvalue = 0.001,
                      Pe = 0.001,
                      filter.probes = TRUE, 
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      label = "hypo", 
                      cores = 1)
#procura pares entre as probes selecionadas e os genes com uma correlação negativa

for (i in 1:11) {
  schematic.plot(pair = Hypo.pair.em, 
                 data = LIHC.mae,
                 group.col = "definition",
                 byProbe = Hypo.pair.em$Probe[i],
                 save = FALSE)
}

heatmapPairs(data = LIHC.mae, 
             group.col = "definition",
             group1 = "Primary solid Tumor",
             group2 = "Solid Tissue Normal",
             annotation.col = c("")
             pairs = Hypo.pair.em,
             filename =  NULL)

scatter.plot(data = LIHC.mae,
             byPair = list(probe = Hypo.pair.em$Probe[5], gene = Hypo.pair.em$GeneID[5]), 
             category = "definition", save = FALSE, lm_line = TRUE)

#De 412 amostras no MAE, 371 são tumorais. No .csv tenho informações sobre 196 e 
#nem todas estão presentes do experimento.
#lihc.metadata
#209-258 - alcohol_history_documented, frequency_of_alcohol_consumption
#history_hepato_carcinoma_risk_factor 
#500 - alcoholic_exposure_category

################################## Hyper analysis ################################
sig.diff.hyper <- get.diff.meth(data = LIHC.mae, group.col = "definition",
                          group1 = group1, group2 = group2,
                          minSubgroupFrac = 0.2, sig.dif = 0.4,
                          diff.dir = "hyper", cores = 1,
                          pvalue = 0.01)
#sig.dif 0.4 e pvalue 0.01: 2409 probes

#Volcano plot
TCGAVisualize_volcano(x = as.data.frame(sig.diff.hyper)[,grep("Minus",colnames(sig.diff.hyper),value = T)],
                      y = sig.diff.hyper$adjust.p,
                      title = paste0("Volcano plot - Probes ", "hypermethylated in ", group1, " vs ", group2, "\n"),
                      label =  c("Not Significant",
                                 paste0("Hypermethylated in ",group1),
                                 paste0("Hypomethylated in ",group1)),
                      ylab =  expression(paste(-Log[10],
                                               " (FDR corrected P-values) [one tailed test]")),
                      xlab =  expression(paste(
                        "DNA Methylation difference (",beta,"-values)")
                      ),
                      x.cut = 0.4, 
                      y.cut = 0.01,
                      filename = "volcanoHyper.png"
)

neargenes.hyper <- GetNearGenes(data=LIHC.mae, probes = sig.diff.hyper$probe, numFlankingGenes = 20)
#total: 48175 genes

Hyper.pair <- get.pair(data = LIHC.mae, 
                      group.col = "definition", 
                      group1 = group1,
                      group2 = group2, 
                      nearGenes = neargenes.hyper, 
                      mode = "unsupervised",
                      permu.size = 61488, 
                      raw.pvalue = 0.001,
                      Pe = 0.001,
                      filter.probes = TRUE, 
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      label = "hyper", 
                      cores = 4)
#Para pvalue 0.01 e sig.dif 0.4 -> 61488: 121, 7 únicas
#length(unique(Hypo.pair$Probe))

for (i in 6:10) {
  scatter.plot(data = LIHC.mae,
               byPair = list(probe = Hyper.pair$Probe[i], gene = Hyper.pair$GeneID[i]), 
               category = "definition", save = TRUE, lm_line = TRUE)
}

heatmapPairs(data = LIHC.mae, 
             group.col = "definition",
             group1 = "Primary solid Tumor",
             group2 = "Solid Tissue Normal",
             pairs = Hyper.pair,
             filename =  NULL)


enriched.motif.hyper <- get.enriched.motif(data = LIHC.mae,
                                     probes = Hyper.pair$Probe,
                                     label = "hyper",
                                     min.incidence = 1,
                                     lower.OR = 1.1)
#Para pvalue 0.01 e sig.dif 0.4 -> 58 enriched motifs

TF.hyper <- get.TFs(data = LIHC.mae, mode = "unsupervised", group.col = "definition", 
                   group1 = group1, group2 = group2, enriched.motif = enriched.motif.hyper,
                   cores = 1, label = "hyper")


########################## Supervisionada #####################################
#save(data.exp, data.meth, iClusters12, file = "tathi.rda")

#Alterando barcodes
# data.meth <- assay(data.meth) #converte o RangedSummarizedExperiment para um dataframe (vc pode usar o proprio dataframe para fazer o MAE)
# colnames(data.meth) <- substr(colnames(data.meth), 1, 16)
# 
# data.exp <- assay(data.exp) #converte o RangedSummarizedExperiment para um dataframe (vc pode usar o proprio dataframe para fazer o MAE)
# colnames(data.exp) <- substr(colnames(data.exp), 1, 16)
# 
# coldata <- subset(iClusters12, iClusters12$primary %in% colnames(data.meth))
# 
# iClusters12$iClusters <- as.character(iClusters12$iClusters)

#pra larissa
#Alterando barcodes
# data.meth <- assay(data.meth) #converte o RangedSummarizedExperiment para um dataframe (vc pode usar o proprio dataframe para fazer o MAE)
# colnames(data.meth) <- substr(colnames(data.meth), 1, 16)
# 
# data.exp <- assay(data.exp) #converte o RangedSummarizedExperiment para um dataframe (vc pode usar o proprio dataframe para fazer o MAE)
# colnames(data.exp) <- substr(colnames(data.exp), 1, 16)
# 
# coldata <- subset(hpv-status, hpv-status$Barcode %in% colnames(data.meth))
# 
# hpv-status$HPV-Status <- as.character(hpv-status$HPV-Status)
#elmer.test.R
class(iClusters12$iClusters)
iClusters12 <- as.data.frame(iClusters12)
iClusters12$iClusters <- as.factor(iClusters12$iClusters) #convert to factor to change the names
levels(iClusters12$iClusters) <- c("Cluster.1", "Cluster.2") #rename the levels
iClusters12$iClusters <- as.character(iClusters12$iClusters) #convert back to character 
length(intersect(colnames(data.exp), iClusters12$primary)) #120
length(intersect(colnames(data.meth), iClusters12$primary)) #120
table(iClusters12$iClusters)

S.LIHC.mae <- createMAE(exp = data.exp, 
                            met = data.meth,
                            colData = iClusters12,
                            filter.probes = distal.probes,
                            linearize.exp = TRUE,
                            met.platform = "450K",
                            genome = "hg38",
                            save = TRUE,
                            TCGA = FALSE) # A logical. FALSE indicate data is not from TCGA (FALSE is default). 
                                         #TRUE indicates data is from TCGA and sample section will automatically filled in )

#table(S.LIHC.mae@colData@listData$iClusters)
# 1  2  3 
#65 55 63

# Análise entre iClusters 1 e 2
Sgroup1.12 <- "Cluster.2"  
Sgroup2.12 <- "Cluster.1"  
Sgroup.col.12 <- "iClusters" 

S.sig.diff12 <- get.diff.meth(data = S.LIHC.mae,
                               group.col = "iClusters",
                               group1 =  "Cluster.2",
                               group2 = "Cluster.1",
                               sig.dif = 0.25,
                               diff.dir = "hypo", 
                               cores = 1,
                               mode ="supervised",
                               minSubgroupFrac = 1,
                               pvalue = 0.05,
                               test = wilcox.test,
                               save = TRUE)

S.neargenes12 <- GetNearGenes(data=S.LIHC.mae,
                     probes= S.sig.diff12$probe,
                     numFlankingGenes = 20)

S.hypo.pair12 <- get.pair(data = S.LIHC.mae,
                    group.col = "iClusters",
                    group1 =  "Cluster.2",
                    group2 = "Cluster.1",
                    nearGenes = S.neargenes12,
                    mode = "supervised",
                    permu.size = 100000,
                    raw.pvalue = 0.05,
                    diff.dir = "hypo",
                    Pe = 0.001, 
                    filter.probes = FALSE, # See preAssociationProbeFiltering function
                    #filter.percentage = 0.05,
                    filter.portion = 0.3,
                    dir.out = "enhancer_result",
                    cores = 1,
                    label = "hypo")

S.Motif12 <- get.enriched.motif(data = S.LIHC.mae,
                                 probes = S.hypo.pair12$Probe,
                                 #dir.out = "result",
                                 label = "hypo",
                                 min.incidence = 10,
                                 lower.OR = 1.1)

S.TF12 <- get.TFs(data = S.LIHC.mae, mode = "supervised", group.col = Sgroup.col.12, 
                        group1 = Sgroup1.12, group2 = Sgroup2.12, enriched.motif = S.Motif12,
                        cores = 1, diff.dir = "hypo", label = "hypo")

# heatmapPairs(data = S.LIHC.mae, 
#              group.col = Sgroup.col.12,
#              group1 = Sgroup1.12,
#              group2 = Sgroup2.12,
#              pairs = S.HypoPair12,
#              filename =  NULL)
# 
# heatmapPairs(data = S.LIHC.mae, 
#              group.col = Sgroup.col.12,
#              group1 = Sgroup1.12,
#              group2 = Sgroup2.12,
#              pairs = S.HyperPair12,
#              filename =  NULL)
#Genes que mais se repetem nos pares hipermetilados: TPRG1L (92 vezes)

heatmapPairs(data = S.LIHC.mae,
             group.col = "iClusters",
             group1 =  "Cluster.2",
             group2 = "Cluster.1",
             #annotation.col = "iClusters",
             #filename =  "Supervised_heatmap_enhancer.png",
             filename =  NULL,
             pairs = S.hypo.pair12) 

for (i in 1:5) {
  scatter.plot(data = S.LIHC.mae,
               byPair = list(probe = S.pair.both$Probe[i], gene = S.pair.both$GeneID[i]), 
               category = "iClusters", save = TRUE, lm_line = TRUE)
}

load("/media/hd/isabela/getTF.hypo.TFs.with.motif.pvalue.rda")
motif12 <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif12,
             save = FALSE)

scatter.plot(data = S.LIHC.mae,
             byTF = list(TF = "NR1I2",
                         probe = S.Motif12$HXA2_HUMAN.H11MO.0.D), 
             category = "iClusters",
             save = FALSE, 
             lm_line = TRUE) 
############################iClusters13#########################

iClusters13 <- read_csv("iClusters13.csv")

class(iClusters13$iClusters)
iClusters13 <- as.data.frame(iClusters13)
iClusters13$iClusters <- as.factor(iClusters13$iClusters) #convert to factor to change the names
levels(iClusters13$iClusters) <- c("Cluster.1", "Cluster.3") #rename the levels
iClusters13$iClusters <- as.character(iClusters13$iClusters) #convert back to character 
length(intersect(colnames(data.exp), iClusters13$primary)) #128
length(intersect(colnames(data.meth), iClusters13$primary)) #128
table(iClusters13$iClusters)

S.LIHC.mae13 <- createMAE(exp = data.exp, 
                        met = data.meth,
                        colData = iClusters13,
                        filter.probes = distal.probes,
                        linearize.exp = TRUE,
                        met.platform = "450K",
                        genome = "hg38",
                        save = TRUE,
                        TCGA = FALSE)

Sgroup1.13 <- "Cluster.3"  
Sgroup2.13 <- "Cluster.1"  
Sgroup.col.13 <- "iClusters"

S.sig.diff13 <- get.diff.meth(data = S.LIHC.mae13,
                                 group.col = Sgroup.col.13,
                                 group1 =  Sgroup1.13,
                                 group2 = Sgroup2.13,
                                 sig.dif = 0.4,
                                 diff.dir = "hypo", 
                                 cores = 1,
                                 mode ="supervised",
                                 minSubgroupFrac = 1,
                                 pvalue = 0.05,
                                 test = wilcox.test,
                                 save = TRUE)

S.neargenes13 <- GetNearGenes(data=S.LIHC.mae13,
                                 probes= S.sig.diff13$probe,
                                 numFlankingGenes = 20) 

S.hypo.pair13 <- get.pair(data = S.LIHC.mae13,
                        group.col = "iClusters",
                        group1 =  "Cluster.3",
                        group2 = "Cluster.1",
                        nearGenes = S.neargenes13,
                        mode = "supervised",
                        permu.size = 100000,
                        raw.pvalue = 0.05,
                        diff.dir = "hypo",
                        Pe = 0.001, 
                        filter.probes = FALSE, # See preAssociationProbeFiltering function
                        #filter.percentage = 0.05,
                        filter.portion = 0.3,
                        dir.out = "enhancer_result",
                        cores = 1,
                        label = "hypo")

heatmapPairs(data = S.LIHC.mae13,
             group.col = "iClusters",
             group1 =  "Cluster.3",
             group2 = "Cluster.1",
             #annotation.col = "iClusters",
             #filename =  "Supervised_heatmap_enhancer.png",
             filename =  NULL,
             pairs = S.hypo.pair13)

S.Motif13 <- get.enriched.motif(data = S.LIHC.mae13,
                                     probes = S.hypo.pair13$Probe,
                                     #dir.out = "result",
                                     label = "hypo",
                                     min.incidence = 10,
                                     lower.OR = 1.1)

S.TF13 <- get.TFs(data = S.LIHC.mae13, mode = "supervised", group.col = Sgroup.col.13, 
                       group1 = Sgroup1.13, group2 = Sgroup2.13, enriched.motif = S.Motif13,
                       cores = 1, diff.dir = "hypo", label = "hypo")

for (i in 2000:2005) {
  scatter.plot(data = S.LIHC.mae13,
               byPair = list(probe = S.pair.both13$Probe[i], gene = S.pair.both13$GeneID[i]), 
               category = "iClusters", save = TRUE, lm_line = TRUE)
}

motif.enrichment.plot(motif.enrichment = "/media/hd/isabela/getMotif.both.motif.enrichment.csv", 
                      significant = list(OR = 1.5,lowerOR = 1.3), 
                      label = "hyper", 
                      summary = TRUE,
                      save = FALSE)

load("/media/hd/isabela/getTF.hypo.TFs.with.motif.pvalue.rda")
motif13 <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif13,
             save = FALSE)
############################iClusters23#########################

iClusters23 <- read_csv("iClusters23.csv")

class(iClusters23$iClusters) #numeric
iClusters23 <- as.data.frame(iClusters23)
iClusters23$iClusters <- as.factor(iClusters23$iClusters) #convert to factor to change the names
levels(iClusters23$iClusters) <- c("Cluster.2", "Cluster.3") #rename the levels
iClusters23$iClusters <- as.character(iClusters23$iClusters) #convert back to character 
length(intersect(colnames(data.exp), iClusters23$primary)) #128
length(intersect(colnames(data.meth), iClusters23$primary)) #128
table(iClusters23$iClusters)

S.LIHC.mae23 <- createMAE(exp = data.exp, 
                          met = data.meth,
                          colData = iClusters23,
                          filter.probes = distal.probes,
                          linearize.exp = TRUE,
                          met.platform = "450K",
                          genome = "hg38",
                          save = TRUE,
                          TCGA = FALSE)

Sgroup1.23 <- "Cluster.3"  
Sgroup2.23 <- "Cluster.2"  
Sgroup.col.23 <- "iClusters"

S.sig.diff23 <- get.diff.meth(data = S.LIHC.mae23,
                                   group.col = Sgroup.col.23,
                                   group1 =  Sgroup1.23,
                                   group2 = Sgroup2.23,
                                   sig.dif = 0.3,
                                   diff.dir = "hypo", 
                                   cores = 1,
                                   mode ="supervised",
                                   minSubgroupFrac = 1,
                                   pvalue = 0.05,
                                   test = wilcox.test,
                                   save = TRUE)

S.neargenes23 <- GetNearGenes(data=S.LIHC.mae23,
                              probes= S.sig.diff23$probe,
                              numFlankingGenes = 20) 

S.hypo.pair23 <- get.pair(data = S.LIHC.mae23,
                          group.col = "iClusters",
                          group1 =  "Cluster.3",
                          group2 = "Cluster.2",
                          nearGenes = S.neargenes23,
                          mode = "supervised",
                          permu.size = 100000,
                          raw.pvalue = 0.05,
                          diff.dir = "hypo",
                          Pe = 0.001, 
                          filter.probes = FALSE, # See preAssociationProbeFiltering function
                          #filter.percentage = 0.05,
                          filter.portion = 0.3,
                          dir.out = "enhancer_result",
                          cores = 1,
                          label = "hypo")

heatmapPairs(data = S.LIHC.mae23,
             group.col = "iClusters",
             group1 =  "Cluster.3",
             group2 = "Cluster.2",
             #annotation.col = "iClusters",
             #filename =  "Supervised_heatmap_enhancer.png",
             filename =  NULL,
             pairs = S.hypo.pair23)

S.Motif23 <- get.enriched.motif(data = S.LIHC.mae23,
                                     probes = S.hypo.pair23$Probe,
                                     #dir.out = "result",
                                     label = "hypo",
                                     min.incidence = 10,
                                     lower.OR = 1.1)

S.TF23 <- get.TFs(data = S.LIHC.mae23, mode = "supervised", group.col = Sgroup.col.23, 
                       group1 = Sgroup1.23, group2 = Sgroup2.23, enriched.motif = S.Motif23,
                       cores = 1, diff.dir = "hypo", label = "hypo")

for (i in 245:250) {
  scatter.plot(data = S.LIHC.mae23,
               byPair = list(probe = S.pair.both23$Probe[i], gene = S.pair.both23$GeneID[i]), 
               category = "iClusters", save = TRUE, lm_line = TRUE)
}

motif.enrichment.plot(motif.enrichment = "/media/hd/isabela/getMotif.both.motif.enrichment.csv", 
                      significant = list(OR = 1.5,lowerOR = 1.3), 
                      label = "both", 
                      summary = TRUE,
                      save = FALSE)

load("/media/hd/isabela/getTF.hypo.TFs.with.motif.pvalue.rda")
motif23 <- colnames(TF.meth.cor)[1]
TF.rank.plot(motif.pvalue = TF.meth.cor, 
             motif = motif23,
             save = FALSE)

filter.motif.probes <- function(input){ #input é o obj retornado da função get.enriched.motif
  result <- NULL
  for (i in 1:length(input)){
    motif.df <- data.frame(a = input[i]) #transforma a lista de probes do primeiro motif em um dataframe
    motif.df <- rename(motif.df, probes = colnames(motif.df[1])) #renomeia a coluna para "probes"
    result <- bind_rows(result, motif.df) #assim é possível juntar as colunas dos dataframes
  }
  result <- unique(result, by="probes") #remove as duplicatas
  return(result)
}
