# Matching AML RNA expression subtypes to genomic features and to AML cell lines in DepMap/CCLE

Sources:
--------
TCGA data: https://gdac.broadinstitute.org/
DepMap/CCLE data: https://depmap.org

TCGA AML subtype genomic feature mapping:
-----------------------------------------
- Description: Matches consensus cluster plus (CCP)-determined expression subtypes from AML RNA-seq data to cytogenetic features
  and mutations for samples. Produces an output HTML file with plotly interactive heatmaps aligning sample subtype assignments with these 
  features, along with reporting enrichment statistics. 
- Code: get_AML_feature_enrichments.py
- Example output: html/AML_exp_subtype_features_heatmaps.html

TCGA AML subtype matching to CCLE AML cell lines:
-------------------------------------------------
- Description: Matches TCGA AML RNA expression subtype markers to CCLE AML cell lines with three different metrics. Produces output HTML
  files with plotly interactive heatmaps indicating strength of association between cell line expression data and that of patient 
  expression subtype markers.
- Code: match_AML_expression_subtypes_to_cells.py
- Example output: html/DMP_protein_coding_all21Q3AML_RNAseq_CCP_6_clusters_markers_diff_1.0_pval_0.0001.html


