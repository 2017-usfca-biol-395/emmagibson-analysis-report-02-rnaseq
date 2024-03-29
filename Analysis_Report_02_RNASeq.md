Analysis Report 2: Gender differences in lung cancer
================
Emma Gibson
November 20, 2017

Introduction
============

In today's society, cancer in general and lung cancer specifically is one of the most common causes of death. There have been numerous studies on the prevalance and risk factors of lung cancer, many of which agree that smoking presents a serious increase in risk for developing lung cancer. (Peto *et al.*, 2000) However, there have been conflicting results as to how biological sex affects risk of lung cancer. Some studies indicate that there there are no significant differences in the overall lung cancer rates between men and women, althought here might be differences in the rates of certain types of lung cancer. (Bain *et al.*, 2004) However, other studies indicate that women are at a higher risk of developing all types of lung cancer. (Zang and Wynder, 1996) However, both of the cited studies agree with each other and other studies that adenocarcinoma rates are higher in females. (Devesa *et al.*, 2005) This type of cancer is characterized by regulation of several distinct genes, including those known to be associated with lung sulfricant proteins. (Taguchi *et al.*, 2011)

Methods
=======

Sample origin and sequencing
----------------------------

The original data that was analyzed in this report was gathered by Seo et al. in their paper studying the transcriptional landscape of lung cancer patients. (Seo *et al.*, 2012) Original tissue samples were taken surgically from cancer patients, and were extracted from both healthy and cancerous tissue. cDNA libraries were sequenced using an Illumina HiSeq machine. The goal of this initial study was to identify potential point mutations and alternative splicing in the genes of cancer patients.

After this initial study, Li et al. re-analyzed the sequence data looking for differences in gene expression between smoking and non-smoking cancer patients. (Li *et al.*, 2015) Using various statistical analysis from R, such as Bionconductor and EdgeR, they compared overall gene transcription levels in smoking and non-smoking cancer patients. Their analysis excluded patients over age 75, and included additional data from 6 nonsmoker patients from another study.

Computational
-------------

The data from the initial study was downloaded off of NCBI, then re-proscessed and re-mapped using several programs. First, trimmomatic was used to properly clean and proscess the sequences from the initial RNAseq dataset. (Bolger *et al.*, 2014) BioMart was used to download the annotated reference genomes that the RNA sequences would be mapped to. (Durinck *et al.*, 2005) This RNAseq information mapped to these reference genomes usingsailfish. (Patro *et al.*, 2014) Finally, the cleaned RNAseq data and its associated metadata was analyzed in R using extensions such as ggplott

Results
=======

| gender | genename |  mean\_count|
|:-------|:---------|------------:|
| male   | EEF1A1   |    183911.87|
| female | EEF1A1   |    149021.94|
| female | SFTPB    |    135305.42|
| male   | SFTPB    |    105009.56|
| male   | MIR6723  |     91501.45|
| male   | FTL      |     77787.92|
| male   | RNA28S5  |     71265.16|
| female | SFTPA2   |     66520.77|
| male   | FN1      |     65454.33|
| male   | CD74     |     63047.61|
| male   | SFTPA2   |     60627.62|
| female | CD74     |     58994.05|
| female | MIR6723  |     58916.57|
| female | RNA28S5  |     58079.55|
| male   | FTH1     |     52390.78|

**Table 1.** The most highly expressed genes in both genders included *SFTPB* and *EEF1A1*.

![](Analysis_Report_02_RNASeq_files/figure-markdown_github-ascii_identifiers/make-barplot-of-highly-expressed-genes-1.png)

**Figure 1.** Fifteen most common genes in male and female patients.

![](Analysis_Report_02_RNASeq_files/figure-markdown_github-ascii_identifiers/make-boxplot-of-highly-expressed-genes-1.png)

**Figure 2.** Expression counts in top 15 genes, by gender and smoking status

``` r
# first, this code finds various factors for each patient
gender_and_smoking_status <- final_table %>%
  group_by(smoking_status, cancer_stage, gender, sample_name_s) %>%
  count() %>%
  arrange(desc(n))

#next, this code interprets it as a dataframe and reads it into ggplot
with(gender_and_smoking_status, table(gender, smoking_status, cancer_stage)) %>%
  as.data.frame() %>% #to make it a ggplot-friendly dataframe
  ggplot(aes(x = smoking_status,
             y = Freq,
             fill = gender)) +
  facet_wrap(~cancer_stage) +
  geom_col(position = "dodge") +
  ggtitle("Smoking status of patients with each cancer stage") +
  xlab("Smoking status") +
  ylab("Number of patients") +
  theme_bw() +
  theme(axis.text.x =
            element_text(angle = 90,
                         hjust = 1))
```

![](Analysis_Report_02_RNASeq_files/figure-markdown_github-ascii_identifiers/gender-and-smoking-1.png)

**Figure 3.** Smoking statues of patients with each cancer stage by gender

``` r
# this code finds the top 9 genes in female patients only
# I chose 9 instead of 15 because the original top 15 had
# 9 genes because some were common in both genders
top_9_f <- final_table %>%
  filter(gender == "female") %>%
  group_by(genename) %>%
  summarize(mean_count = mean(counts_lengthscaledtpm)) %>%
  arrange(desc(mean_count)) %>%
  head(n = 9)

#next, this line puts it into a nice markdown table
kable(top_9_f)
```

| genename |  mean\_count|
|:---------|------------:|
| EEF1A1   |    149021.94|
| SFTPB    |    135305.42|
| SFTPA2   |     66520.77|
| CD74     |     58994.05|
| MIR6723  |     58916.57|
| RNA28S5  |     58079.55|
| FTL      |     51200.27|
| SFTPA1   |     50997.97|
| TPT1     |     40372.54|

**Table 2.** Top 9 genes in female patients

``` r
# this puts the top 9 female genes into a facet-friendly format
top_genes_f <- top_9_f %>%
  ungroup() %>%
  select(genename) %>%
  unique() %>%
  pull()

# this code makes a figure similar to Figure 2, but with top female genes
# it doesn't plot top_genes_f because I wanted to include male patients
final_table %>%
  filter(genename %in% top_genes_f) %>%
  ggplot(aes(x = genename,
             y = counts_lengthscaledtpm,
             fill = gender)) +
    geom_boxplot() +
    facet_wrap(~smoking_status) +
    xlab("Gene name") +
    ylab("Scaled read counts per gene") +
    ggtitle("Read counts per gene of genes common to female smokers") +
  theme_bw() +
  theme(axis.text.x =
            element_text(angle = 90,
                         hjust = 1))
```

![](Analysis_Report_02_RNASeq_files/figure-markdown_github-ascii_identifiers/female-gene-graph-1.png)

**Figure 4.** Expression counts in top 9 most common female genes, by gender and smoking status

``` r
# this gets the appropriate metadata to make a figure faceted by smoking status
top_f_in_all_status <- final_table %>%
  filter(genename %in% top_genes_f) %>%
  group_by(gender, genename, smoking_status) %>%
  summarize(mean_count = mean(counts_lengthscaledtpm))

top_f_in_all_status %>%
  ggplot(aes(x = genename,
             y = mean_count,
             fill = gender)) +
  facet_wrap(~smoking_status) +
  geom_col(position = "dodge") +
  ggtitle("A. Expression of most common genes for female
          patients in all patients by smoking status") +
  xlab("Gene name") +
  ylab("Mean read counts per gene") +
  theme_bw() +
  theme(axis.text.x =
            element_text(angle = 90,
                         hjust = 1))
```

![](Analysis_Report_02_RNASeq_files/figure-markdown_github-ascii_identifiers/smoking-in-female-genes-1.png)

``` r
# this gets the appropriate metadata to make a figure faceted by cancer stage
top_f_in_all_stage <- final_table %>%
  filter(genename %in% top_genes_f) %>%
  group_by(gender, genename, cancer_stage) %>%
  summarize(mean_count = mean(counts_lengthscaledtpm))

top_f_in_all_stage %>%
  ggplot(aes(x = genename,
             y = mean_count,
             fill = gender)) +
  facet_wrap(~cancer_stage) +
  geom_col(position = "dodge") +
  ggtitle("B. Expression of most common genes for female
          patients in all patients by cancer stage") +
  xlab("Gene name") +
  ylab("Mean read counts per gene") +
  theme_bw() +
  theme(axis.text.x =
            element_text(angle = 90,
                         hjust = 1))
```

![](Analysis_Report_02_RNASeq_files/figure-markdown_github-ascii_identifiers/smoking-in-female-genes-2.png)

**Figure 5.** Abundance of top 9 female genes by gender, faceted by A. smoking status and B. cancer stage

Discussion
==========

Initially, I noticed that gene expression was relatively even in both genders for the most common genes, witht he exception of EEF1A1 and SFTB (Figure 1). After breaking further breaking down this data by smoking status, I found a similar correlation between gene expression between the genders, but I noticed that there seemed to be less data on female smokers (Figure 2). When I looked into this further, I found that most of the female patients had never smoked and were in early stages of cancer. (Figure 3) Given this information, I decided to investigate whether the most common genes in female patients might be different than the most common genes in all patients, and found that they were (Table 2). Most of the genes tht were most commom in female patients were the same as the most common in general patients, and they showed similar expression patterns as well (Figure 4) Lastly, I decided to examine the expression of these genes compared to smoking status and cancer stage, and found that the genes that were more common in females were also most common in nonsmokers with early cancer stages, which makes sense given the other analyses (Figure 5).

Of the few genes that are consistently more common in female patients than males are genes like SFTBP, which seems to be especially common in cancer stage 1B (Figure 5B). Although it is quite common in both genders at this stage, SFTBP is the most common gene in females at this stage, but not in males. Although this could be due to the gene simply being expressed in early stages of cancer, the difference could potentially hint that a different type of cancer common to females at this stage. Studies have shown that the SFTBP gene is associated with adenocarcinoma, the type of lung cancer that is especially common in females. (Taguchi *et al.*, 2011) Given that Li et al. showed that there are differences in gene expression among lung cancer patients who do and do not smoke, it is possible that adenocarcinomaadenocarcinoma might be more common in people who do not smoke, rather than females. (Li *et al.*, 2015) However, this gene is still more common in females who have never smoked than males who have never smoked, indicating that SFTBP's abundance in female patients is more likely due to gender than smoking status.

Sources Cited
=============

Bain,C. *et al.* (2004) Lung cancer rates in men and women with comparable histories of smoking. *Journal of the National Cancer Institute*, **96**, 826–834.

Bolger,A.M. *et al.* (2014) Trimmomatic: A flexible trimmer for illumina sequence data. *Bioinformatics*, **30**, 2114–2120.

Devesa,S.S. *et al.* (2005) International lung cancer trends by histologic type: Male: Female differences diminishing and adenocarcinoma rates rising. *International journal of cancer*, **117**, 294–299.

Durinck,S. *et al.* (2005) BioMart and bioconductor: A powerful link between biological databases and microarray data analysis. *Bioinformatics*, **21**, 3439–3440.

Li,Y. *et al.* (2015) RNA-seq analysis of lung adenocarcinomas reveals different gene expression profiles between smoking and nonsmoking patients. *Tumor Biology*, **36**, 8993–9003.

Patro,R. *et al.* (2014) Sailfish enables alignment-free isoform quantification from rna-seq reads using lightweight algorithms.

Peto,R. *et al.* (2000) Smoking, smoking cessation, and lung cancer in the uk since 1950: Combination of national statistics with two case-control studies. *Bmj*, **321**, 323–329.

Seo,J.-S. *et al.* (2012) The transcriptional landscape and mutational profile of lung adenocarcinoma. *Genome research*, **22**, 2109–2119.

Taguchi,A. *et al.* (2011) Lung cancer signatures in plasma based on proteome profiling of mouse tumor models. *Cancer cell*, **20**, 289–299.

Zang,E.A. and Wynder,E.L. (1996) Differences in lung cancer risk between men and women: Examination of the evidence. *JNCI: Journal of the National Cancer Institute*, **88**, 183–192.
