# Introduction to Study

The following code was used for a study that analyzed the effects of vitamin D on preeclampsia for pregnant women (Mirzakhani et al. 2016). Study participants were divided into two groups: one group was given a daily dose 400 IU of vitamin D and the other group was given a daily does of 4,400 IU of vitamin D. Blood samples of patients were collected and expression of mRNA and miRNA genes were investigated. Whether the women were diagonosed with preeclampsia further divided the groups. 

# Outline of Code

Before continuing, it is important to note that patient data is not included because of patient privacy. 

The genes used in the analysis come from the following sources:

- LCC stands for Largest Connected Components. A regulatory network was inferred from expression data and these genes were in the largest connected component of the network graph.
- JCI stands for Journal of Clinical Investigation. The Mirzakhani et al. 2016 paper identified that these genes are differentially expressed.

The purpose of each R script are as follows:

- Figure2.R - This code predicted the targets of differentially expressed miRNA. The miRNA which targeted differentially expressed mRNA are stored in a list.
- PCA.R - a PCA plot is drawn which plots the genomic expression of different patients
- SignificantGeneCorrelations.R - heatmaps of the LCC and JCI genes were drawn which are in the heatmaps folder
- generateFigures.R - The predicted targets of differentially expressed miRNA genes. The miRNA which targeted differentially expressed mRNA genes were shown in the images/plot.png figure.
- venn.R - a Venn Diagram showing the number of genes found in each group

Citation

Mirzakhani, H., Litonjua, A. A., McElrath, T. F., O’Connor, G., Lee-Parritz, A., Iverson, R., Macones, G., Strunk, R. C., Bacharier, L. B., Zeiger, R.,       
  Hollis, B. W., Handy, D. E., Sharma, A., Laranjo, N., Carey, V., Qiu, W., Santolini, M., Liu, S., Chhabra, D., … Weiss, S. T. (2016). Early pregnancy 
  vitamin D status and risk of preeclampsia. Journal of Clinical Investigation, 126(12), 4702–4715. https://doi.org/10.1172/jci89031 
