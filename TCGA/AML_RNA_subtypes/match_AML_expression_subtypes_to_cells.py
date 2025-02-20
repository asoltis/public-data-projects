import sys, os, _pickle
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics.pairwise import cosine_similarity as COS
sys.path.insert(0, '/Users/anthonysoltis/code/arsoltis_hunter_github/common_utils/')
from normalization import zscore_matrix
import plotly
import plotly.graph_objects as go
from plotly import offline

# CCLE cell line and transcriptomics data
ccle_dir = '' # Path python pickle objects containing sample info and RNA-seq processed data
CLD = _pickle.load(open(ccle_dir + 'DepMap_Public_21Q3_sample_info.pkl', 'rb'))
which_rna = 'DMP_protein_coding'
#which_rna = 'CCLE_RSEM'
if which_rna == 'DMP_protein_coding':
    RNAD = _pickle.load(open(ccle_dir + 'RNAseq/processed/' + 'DMP_21Q3_CCLE_expression_protein_coding.pkl', 'rb'))
    samples = RNAD['samples']
    genes = RNAD['genes_symbols']
    MAT = RNAD['matrix']
elif which_rna == 'CCLE_RSEM':
    RNAD = _pickle.load(open(ccle_dir + 'RNAseq/processed/' + 'CCLE_RNAseq_rsem_genes_tpm_20180929.pkl', 'rb'))
    MAT = RNAD['matrix']
    MAT = np.log2(MAT + 1) # Convert to logarithm
    samples = []
    for ccle_name in RNAD['samples']:
        cln = ccle_name.split('_')[0]
        if cln == 'D341MED': cln = 'D341'
        dmid = CLD['cellLineNameStrip2ID'][cln]
        samples.append(dmid)
    genes = RNAD['genes']

# Standardize RNA expression data according to tissue/cancer type
dmid_use = [] # list of cell line DepMap IDs to use
cell_lines0 = []
for dmid in CLD['depmapID2info']:
    cl_info = CLD['depmapID2info'][dmid]
    prim_dis = cl_info['primary_disease']
    lin_sub = cl_info['lineage_subtype']
    if prim_dis == 'Leukemia' and lin_sub == 'AML':
        dmid_use.append(dmid)
        cell_lines0.append(CLD['depmapID2info'][dmid]['stripped_cell_line_name'])
si_use = []
for si, samp in enumerate(samples):
    if samp in dmid_use: si_use.append(si)
MATRED = MAT[:, si_use] # reduce to select cell lines
samples = np.array(samples)[si_use]
genes = np.array(genes)
MATREDZ = zscore_matrix(MATRED) # z-score matrix

# TCGA expression clusters and marker genes
mdiff, mq = 1.0, 1e-3 # thresholds for marker genes
TCGA_markers_fn = 'TCGA/LAML/Firebrowse/mRNAseq/CCP/LAML-TB.subclassmarkers.txt'
markers = {}
for li, line in enumerate(open(TCGA_markers_fn)):

    if li in [0, 1]: continue
    gstr, p, diff, q, clust = line.strip().split('\t')
    gene = gstr.split('|')[0]
    if gene == '?': continue

    if abs(float(diff)) > mdiff and float(q) < mq:
        if clust not in markers: markers[clust] = {}
        markers[clust][gene] = float(diff)
for clust in sorted(markers):
    print(clust, len(markers[clust]))
print()

# Loop over cell lines of interest and compute match statistics
cell_lines0 = sorted(cell_lines0)
METRICS = ['Pearson', 'Spearman', 'NTP'] # metrics for comparisons
RES_MATS = {}
cell_lines, ci_use = [], []
for mi, METRIC in enumerate(METRICS):
    print(METRIC)
    RES_MATS[METRIC] = {}
    RES_MATS[METRIC]['metric'] = np.zeros((len(cell_lines0), len(markers)))
    RES_MATS[METRIC]['pvals'] = np.zeros((len(cell_lines0), len(markers))) * np.nan

    for clii, cl in enumerate(cell_lines0):
        try:
            cli = np.argwhere(samples == CLD['cellLineNameStrip2ID'][cl])[0][0]
        except:
            cli = None
        print(cl, cli)
        if cli == None: continue
        else:
            if mi == 0:
                ci_use.append(clii)
                cell_lines.append(cl)

        for clusti, clust in enumerate(sorted(markers)):
            X, Y = [], [] # lists for cell line and marker expression
            for gene in markers[clust]:
                try: gi = np.argwhere(genes == gene)[0][0]
                except: continue

                if METRIC in ['Pearson','Spearman']:
                    X.append(markers[clust][gene])
                elif METRIC == 'NTP':
                    if markers[clust][gene] > 0:
                        X.append(1)
                    else: X.append(-1)

                # cell line expression
                Y.append(MATREDZ[gi, cli])

            # Compute metrics
            if METRIC == 'Pearson':
                metric, pval = pearsonr(X, Y)
            elif METRIC == 'Spearman':
                metric, pval = spearmanr(X, Y)
            elif METRIC == 'NTP':
                metric = COS(np.array(X).reshape(-1,1).T, np.array(Y).reshape(-1,1).T)[0][0]
                pval = np.nan

            RES_MATS[METRIC]['metric'][clii, clusti] = metric
            RES_MATS[METRIC]['pvals'][clii, clusti] = pval
            #print(cl)
            #print(clust, metric, pval)

    print()
for METRIC in RES_MATS:
    RES_MATS[METRIC]['metric'] = RES_MATS[METRIC]['metric'][ci_use, :]
    RES_MATS[METRIC]['pvals'] = RES_MATS[METRIC]['pvals'][ci_use, :]
    print(RES_MATS[METRIC]['pvals'].shape, len(ci_use))

print('\t'.join(cell_lines))
for METRIC in METRICS:
    print(METRIC)
    print(RES_MATS[METRIC]['metric'])
    print()

# Create plotly plots
plotly_plots = []
plot_titles = []
collabels = ['C%s' % (x) for x in sorted(markers)]
cmap = 'RdBu_r' #'Inferno'
zmin, zmax = -0.5, 0.5
for METRIC in METRICS:
    mat = RES_MATS[METRIC]['metric']
    pmat = RES_MATS[METRIC]['pvals']

    text = np.array(['P-value: {}'.format('%0.3g'%(x)) for x in pmat[::-1,:].reshape(-1)]).reshape(pmat.shape)
    hmaply = go.Heatmap(z = mat[::-1, :], x = collabels, y = cell_lines[::-1],
                        colorscale = cmap, zmin = zmin, zmax = zmax,
                        hovertemplate = 'Cell line: %{y}<br>Cluster: %{x}</br>Metric: %{z}<br>%{text}</br><extra></extra>',
                        text = text)

    plotly_plots.append(hmaply)
    plot_titles.append('Cell line vs. TCGA cluster: %s' % (METRIC))

# Save plotly html
#plot_fn = 'results/%s_RNAseq_CCP_%d_clusters_markers_diff_%s_pval_%s.html' % (which_rna, len(markers), str(mdiff), str(mq))
plot_fn = 'results/%s_all21Q3AML_RNAseq_CCP_%d_clusters_markers_diff_%s_pval_%s.html' % (which_rna, len(markers), str(mdiff), str(mq))
plotlywidth = 500
plotlyheight = 1000 #500
dashboard = open(plot_fn, 'w')
dashboard.write("<html><head></head><body>" + "\n")
for pi, pplot in enumerate(plotly_plots):
    fig = go.Figure(data = pplot)
    layout = go.Layout(plot_bgcolor='#777777', xaxis = {'showgrid':False, 'tickangle':-90}, yaxis = {'showgrid':False},
                       title = plot_titles[pi], width = plotlywidth, height = plotlyheight,
                       font = {'family':'Arial'})
    fig['layout'].update(layout)
    inner_html = fig.to_html(include_plotlyjs='cdn').split('<body>')[1].split('</body>')[0]
    dashboard.write(inner_html)
dashboard.write("</body></html>" + "\n")
dashboard.close()


