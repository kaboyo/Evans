**#Project 1**


**https://substack.com/app-link/post?publication_id=1508290&post_id=156174368&utm_source=post-email-title&utm_campaign=email-post-title&isFreemail=true&r=4bu9jt&token=eyJ1c2VyX2lkIjoyNjE3NTI1MzcsInBvc3RfaWQiOjE1NjE3NDM2OCwiaWF0IjoxNzM4MzYxMTI2LCJleHAiOjE3NDA5NTMxMjYsImlzcyI6InB1Yi0xNTA4MjkwIiwic3ViIjoicG9zdC1yZWFjdGlvbiJ9.KpHUBc2FS3GZkpQ3yemm634e9ThutDCDfRVyErkHESg**

#Uncovering Amiloride's Mechanisms of Action in Multiple Myeloma: A Bioinformatics Approach

This project investigates the mechanisms by which amiloride induces cell death in multiple myeloma cells through comprehensive bioinformatics analyses.


I'm excited to share an in-depth project walkthrough titled Uncovering Amiloride's Mechanisms of Action in Multiple Myeloma: A Bioinformatics Approach

If you're interested in learning how bioinformatics analysts and computational biologists contribute to the drug discovery and evaluation process, this project is for you.

This project investigates Amiloride's mechanisms of action in multiple myeloma cells through comprehensive bioinformatics analyses. By integrating gene co-expression network analysis, differential expression analysis, protein-protein interaction network analysis, gene ontology, and pathway analysis, the project identifies key regulatory hubs and disrupted pathways following Amiloride treatment.


**Github repository**
```
https://github.com/evanpeikon/Amilioride_Drug_MOA?utm_source=substack&utm_medium=email

```
```
🧬 Introduction
Background and Overview
Multiple myeloma is a form of blood cancer characterized by the uncontrolled proliferation of plasma cells in the bone marrow. This disease is associated with a variety of complications, including anemia, kidney dysfunction, and immune suppression. Despite advancements in treatment options, such as proteasome inhibitors, immunomodulatory drugs, and monoclonal antibodies, multiple myeloma remains largely incurable, and resistance to therapy often develops over time. As a result, there is a continued need for novel therapeutic approaches to improve patient outcomes.

Amiloride is a potassium-sparing diuretic commonly used to treat conditions like hypertension and heart failure. It works by inhibiting the epithelial sodium channel in the kidneys, preventing sodium reabsorption and promoting potassium retention. Interestingly, multiple studies (1, 2, 3) have suggested that Amiloride may have potential as a cancer therapeutic, including in multiple myeloma. It is believed that Amiloride exerts its antitumor effects by targeting multiple cellular pathways involved in tumor progression, including ion channel regulation, apoptosis, and cell migration. Additionally, Amiloride has shown synergistic effects when combined with other chemotherapeutic agents, such as dexamethasone and melphalan.

**Given Amiloride’s promising effects in preclinical studies and its established safety profile, there is growing interest in investigating its potential as a therapeutic agent for multiple myeloma. Thus, the purpose of this project is to explore Amiliorides mechanisms of action at the molecular level to determine how it can be effectively integrated into existing treatment regimens.**

By leveraging various bioinformatics approaches, including gene co-expression network, differential expression, protein-protein interaction network, and functional enrichment analysis, this project will provide a comprehensive view of how Amiloride impacts cellular processes at different levels of biological organization. The results of this project will enhance our understanding of Amiloride’s therapeutic potential and guide future research on its use in the treatment of multiple myeloma.

```
```
This repository is structured into four main sections:

Introduction: Provides an overview of the project, scientific background, the data sources used in this project, and the analyses employed.
Project Staging: Details project dependencies, data loading and preperation, exploratory data analysis (EDA), quality control, filtering, and normalization steps.
Analyses: Includes weighted gene co-expression network analysis, differential expression analysis, protein-protein interaction network analysis, and functional enrichment analysis. Summarizations of the results for each analysis are also provided here.
Results and Conclusion: Synthesizes the results from all analyses and discusses their broader biological context and implications.


Data Sources and Analyses
The data for this project is derived from a study titled "Amiloride, An Old Diuretic Drug, Is a Potential Therapeutic Agent for Multiple Myeloma," which is publicly available through the Gene Expression Omnibus (GEO) under accession number GSE95077. This study explores the effects of Amiloride on multiple myeloma cell lines and a xenograft mouse model, providing valuable insights into the drug's potential as a therapeutic agent.

```

```
The bulk RNA-sequencing data for the JJN3 cell lines used in this project are accessible through the following sample accession numbers: GSM2495770 , GSM2495771 , GSM2495772 , GSM2495767 , GSM2495768 , GSM2495769.

Note: The JJN3 lines, used in this analysis, are established from the bone marrow of a 57-year-old woman with plasma cell leukemia, an aggresive form of multiple myeloma.

In this study, Amiloride’s impact was assessed in multiple myeloma cell lines, alongside a xenograft mouse model, where human tumor cells were implanted into immunocompromised mice. The research aimed to evaluate the drug’s toxicity, specifically its ability to induce cell death (apoptosis) in tumor cells. Additionally, the study sought to understand the underlying mechanisms through a combination of RNA sequencing, quantitative PCR (qRT-PCR), immunoblotting, and immunofluorescence assays.

The investigators observed that Amiloride treatment led to apoptosis in a broad panel of multiple myeloma cell lines and in the xenograft mouse model. Furthermore, they discovered that Amiloride exhibited a synergistic effect when combined with other chemotherapeutic agents, including dexamethasone and melphalan. These findings suggest that Amiloride could enhance the efficacy of existing treatments, potentially improving therapeutic outcomes for multiple myeloma patients.

For this project, I will analyze the RNA-sequencing data from this study to gain deeper insights into Amiloride’s mechanisms of action. Specifically, the following analyses will be performed:

Weighted Gene Co-Expression Network Analysis: WGCNA will be used to identify coordinated patterns of gene expression across both the treatment and control groups. This analysis will help uncover underlying biological processes and regulatory networks that control gene expression. WGCNA is particularly useful for identifying upstream regulatory factors that may influence the expression of multiple genes, and it provides a global view of transcriptional organization, revealing how genes with spreimilar expression profiles are grouped together, often reflecting shared functions or regulatory mechanisms.

Differential Expression Analysis: This analysis will focus on identifying genes that show significant changes in expression between the treatment (amiloride) and control groups. Differential expression analysis will highlight genes that are directly affected by the drug, identifying those that are upregulated or downregulated as a result of amiloride treatment. By isolating these changes, we can pinpoint the most biologically relevant genes that may play key roles in amiloride’s effects on multiple myeloma cells.

Protein-Protein Interaction (PPI) Network Analysis: Building on the results of the differential expression analysis, PPI network analysis will shift the focus to the interactions between the proteins encoded by the differentially expressed genes. This analysis will help us understand how amiloride’s effects may be mediated at the protein level, revealing how these proteins interact with each other to form complexes that drive cellular processes. By mapping out these interactions, we can connect the gene expression changes to functional outcomes and identify critical protein networks that are altered by amiloride.

Functional Enrichment Analysis: To provide a broader biological context for the differentially expressed genes and proteins, functional enrichment analysis will categorize them into known biological functions and pathways. Gene Ontology (GO) analysis will assign genes to specific functional terms, such as biological processes, molecular functions, and cellular components. Additionally, pathway analysis (e.g., KEGG or Reactome) will map genes to known signaling and metabolic pathways. This step will abstract the results from the gene and protein levels to a higher-level understanding of the biological processes and molecular pathways that are impacted by amiloride treatment, linking gene expression changes and protein interactions to functional consequences within the cell.
By performing these analyses, we will examine Amiloride’s effects at multiple levels of biological organization—from gene expression changes to protein interactions and functional pathway alterations—allowing us to gain a comprehensive understanding of its mechanisms of action in multiple myeloma.
```

```
#🧬 Project Staging
This section is broken into three sub-sections:

Project Dependencies: Import the necessary libraries for downstream analyses.
Load, Inspect, and Prepare Data: Steps to load raw datasets, inspect data structure, and prepare it for downstream processing.
Quality Control, Filtering, and Normalization: Applying quality control measures and data transformations to generate high-quality, reliable data for down stream analysis.

```
```
SmallSeqFlow: A streamlined, notebook-based bulk RNA-seq analysis pipeline optimized for small-sample studies.

https://github.com/evanpeikon/SmallSeqFlow
I recently released SmallSeqFlow, a streamlined notebook-based pipeline for analyzing bulk RNA-seq data from small-sample studies. Based on community feedback, I've now expanded its capabilities with three new modules that enable deeper biological insights:

The first addition is a gene co-expression analysis module that reveals gene regulatory networks through correlation-based approaches. It calculates pairwise Spearman correlations between genes, applies significance thresholds to identify meaningful relationships, and provides tools for network visualization and hub gene identification. These features help researchers understand how treatments affect gene regulatory networks, even with limited samples.

The second module focuses on protein-protein interaction (PPI) networks, connecting transcriptional changes to protein-level interactions through the STRING database. It offers tools for network construction, community detection, and centrality analysis to identify key regulatory proteins. The module's visualization components make it easy to examine both global network structures and focused subnetworks of highly connected proteins, helping identify potential therapeutic targets.

The third extension adds functional enrichment analysis capabilities. This module integrates with multiple databases (GO, KEGG, and Reactome) to identify enriched biological processes, molecular functions, and pathways within differentially expressed genes or network modules.

You can access SmallSeqFlow with the link below. If you find this resource helpful, please give it a 🌟 or share it with someone who might benefit from it!



```
