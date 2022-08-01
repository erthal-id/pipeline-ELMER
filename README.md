# ELMER-pipeline
### Esse repositório tem como objetivo disponibilizar o script que foi utilizado para análise através do pacote [ELMER](https://bioconductor.org/packages/ELMER/). 
___
O ELMER (Enhancer Linking by Methylation-Expression Relationship) (Silva TC et al., 2018) é um pacote que correlaciona a metilação de regiões distais ou promotoras com a expressão de genes próximos, a fim de identificar sítios de ligação de fatores de transcrição e elucidar os mecanismos de regulação epigenética presentes no objeto de estudo. 
<br>

## ELMER Workflow
<img src="[https://github.com/erthal-id/ELMER-pipeline/blob/794896b76414bc72d45c9c2400a86d45d68dca61/images/ELMER-workflow.png](https://github.com/erthal-id/pipeline-ELMER/blob/main/images/ELMER-workflow.png)" width ='700px' height='450px'>
<br>

## Para instalar o pacote, basta rodar no R o seguinte código:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ELMER")
```

## Citação
Silva TC, Coetzee SG, Gull N, Yao L, Hazelett DJ, Noushmehr H, Lin D, Berman BP (2018). “ELMER v.2: An R/Bioconductor package to reconstruct gene regulatory networks from DNA methylation and transcriptome profiles.” Bioinformatics. doi: 10.1093/bioinformatics/bty902.
