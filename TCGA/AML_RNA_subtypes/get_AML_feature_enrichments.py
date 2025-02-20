import sys, os, _pickle, csv
import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import chi2_contingency as CHI2
from scipy.stats import fisher_exact as FET
import plotly
import plotly.graph_objects as go
from plotly import offline
from matplotlib import pyplot as plt

'''
Associate AML expression subtypes (from CCP assignments) with mutations and cytogenetic features.
Create HTML with plotly displays for enriched features with expression subtypes.
'''

# RNA clusters
clusters = {}
cluster_totals = {}
use_clust = 'optimal'
clust_vals = set()
if use_clust == 'optimal':
    clust_fn = 'LAML-TB.bestclus.txt' # optimal CCP assignments
    for line in open(clust_fn):
        if line.startswith('SampleName'): continue
        samp, clust, sil = line.strip().split('\t')
        clust = 'C' + clust
        clusters[samp] = clust
        clust_vals.add(clust)
        if clust not in cluster_totals: cluster_totals[clust] = 1
        else: cluster_totals[clust] += 1
clust_vals = sorted(list(clust_vals))
clust_samples = np.array(sorted(clusters))
clust_assign = np.array([clusters[x] for x in clust_samples])

# cytogenetics from cBioPortal sample data
sample_fn = 'TCGA/LAML/cBioPortal/laml_tcga_pub_clinical_data.tsv'
cytogen = {}
for line in open(sample_fn):
    if line.startswith('Study'): continue
    l = line.strip().split('\t')
    sid, cyt = l[2], l[9] # sample id and cytogenetic code
    cytogen[sid] = cyt
clist = []
for samp in cytogen:
    clist.append(cytogen[samp])
ccount = Counter(clist)
cyt_class = sorted(ccount)

# associations with cytogenetics
clust_i = {x:i for i,x in enumerate(clust_vals)}
cyt_j = {x:j for j,x in enumerate(cyt_class)}
clust_cyt_mat = np.zeros((len(clust_i), len(cyt_j)))
CYT_MAT_BIN = np.zeros((len(cyt_class), len(clust_samples)))
for samp in clusters:
    try:
        sc = clusters[samp]
        cyt = cytogen[samp]
        si = np.argwhere(clust_samples == samp)[0][0]
    except KeyError: continue
    i, j = clust_i[sc], cyt_j[cyt]
    clust_cyt_mat[i, j] += 1
    CYT_MAT_BIN[j, si] = 1
print('\t'.join(clust_vals))
for j in range(0, clust_cyt_mat.shape[1]):
    vals = ['%d'%(x) for x in clust_cyt_mat[:, j]]
    cyt = cyt_class[j]
    ol = '\t'.join(vals + [cyt])
    print(ol)

# Overall
chi2, pval, dof, expected = CHI2(clust_cyt_mat)
print('X^2: %0.2f (p-value %0.3g)' %(chi2, pval))

# Individual
cyt_stats = []
for clust in clust_i:
    i = clust_i[clust]
    otheri = [ii for ii in range(0, clust_cyt_mat.shape[0]) if ii != i]
    for feat in cyt_j:
        if feat == 'N.D.': continue
        j = cyt_j[feat]
        otherj = [jj for jj in range(0, clust_cyt_mat.shape[1]) if jj != j]
        PP = clust_cyt_mat[i, j]
        PN = sum(clust_cyt_mat[i, otherj])
        NP = sum(clust_cyt_mat[otheri, j])
        nnrc = clust_cyt_mat[otheri, :]
        nnrc = nnrc[:, otherj]
        NN = np.sum(nnrc)

        table = np.array([[PP, PN], [NP, NN]])
        odds, fpv = FET(table, alternative = 'greater')
        if odds > 1 and fpv < 0.05:
            print(clust, feat, odds, fpv, PP, PN, NP, NN, PP+PN+NP+NN)
            try:
                odds_str = '%0.3f'%(odds)
            except: odds_str = 'Inf'
            cyt_stats.append([clust, feat, odds_str, '%0.3g'%(fpv), '%d:%d:%d:%d'%(PP,PN,NP,NN)])
cyt_stats_df = pd.DataFrame(data = cyt_stats,
                            columns = ['Cluster','Cytogenic class','Odds ratio','p-value','ClustCyt:ClustNoCyt:OtherCyt:OtherNoCyt'])

### Mutations
maf_fn = 'mutations/LAML-TB.final_analysis_set.maf'
GOI = ['FLT3-ITD','FLT3-Other','FLT3-All','NPM1','DNMT3A','IDH2','IDH1','TET2','RUNX1','TP53','NRAS','PTPN11','ASXL1','CEBPA','WT1','KIT','KRAS']
GOIG = ['FLT3','NPM1','DNMT3A','IDH2','IDH1','TET2','RUNX1','TP53','NRAS','PTPN11','ASXL1','CEBPA','WT1','KIT','KRAS']
goi_j = {x:j for j, x in enumerate(GOI)}
goig_j = {x:j for j, x in enumerate(GOIG)}
vset = set()
MUT_MAT = np.zeros((len(clust_i), len(goi_j)))
MUT_MAT_BIN = np.zeros((len(GOIG), len(clust_samples)))
MUT_MAT_EFF = np.chararray((len(GOIG), len(clust_samples)), unicode = True, itemsize = 256)
for line in open(maf_fn):
    if line.startswith('Hugo'): continue
    l = line.strip().split('\t')
    gene, vclass, samp = l[0], l[8], l[15]
    samp = '-'.join(samp.split('-')[0:4])[0:-1] # remove last letter character

    if vclass in ['Silent', 'Intron', "3'UTR", "5'UTR", 'RNA', "5'Flank", 'IGR']: continue
    vset.add(vclass)

    if gene in ['FLT3']:
        try:
            sc = clusters[samp]
            i = clust_i[sc]
            si = np.argwhere(clust_samples == samp)[0][0]
        except KeyError: continue

        if vclass == 'In_Frame_Ins':
            j = goi_j['FLT3-ITD']
            MUT_MAT[i, j] += 1
        else:
            j = goi_j['FLT3-Other']
            MUT_MAT[i, j] += 1

        jj = goi_j['FLT3-All']
        MUT_MAT[i, jj] += 1

        gi = goig_j[gene]
        MUT_MAT_BIN[gi, si] = 1
        if MUT_MAT_EFF[gi, si] != '':
            MUT_MAT_EFF[gi, si] = MUT_MAT_EFF[gi, si] + ';' + vclass
        else:
            MUT_MAT_EFF[gi, si] = vclass

    elif gene in GOI:
        try:
            sc = clusters[samp]
            i = clust_i[sc]
            si = np.argwhere(clust_samples == samp)[0][0]
        except KeyError: continue
        j = goi_j[gene]
        #MUT_MAT[i, j] += 1

        gi = goig_j[gene]
        MUT_MAT_BIN[gi, si] = 1
        if MUT_MAT_EFF[gi, si] != '':
            MUT_MAT_EFF[gi, si] = MUT_MAT_EFF[gi, si] + ';' + vclass
        else:
            MUT_MAT_EFF[gi, si] = vclass
            MUT_MAT[i, j] += 1

print('\t'.join(clust_vals))
for j in range(0, MUT_MAT.shape[1]):
    vals = ['%d'%(x) for x in MUT_MAT[:, j]]
    gene = GOI[j]
    ol = '\t'.join(vals + [gene])
    print(ol)
# Individual
mut_stats = []
for clust in clust_i:
    i = clust_i[clust]
    otheri = [ii for ii in range(0, MUT_MAT.shape[0]) if ii != i]
    for gene in goi_j:
        j = goi_j[gene]
        otherj = [jj for jj in range(0, MUT_MAT.shape[1]) if jj != j]
        PP = MUT_MAT[i, j]
        PN = cluster_totals[clust] - PP
        NP = sum(MUT_MAT[otheri, j])
        NN = sum([cluster_totals[x] for x in cluster_totals if x != clust]) - NP

        table = np.array([[PP, PN], [NP, NN]])
        odds, fpv = FET(table, alternative = 'greater')
        if odds > 1 and fpv < 0.05:
            print(clust, gene, odds, fpv, PP, PN, NP, NN, PP+PN+NP+NN)
            try:
                odds_str = '%0.3f'%(odds)
            except: odds_str = 'Inf'
            mut_stats.append([clust, gene, odds_str, '%0.3g'%(fpv), '%d:%d:%d:%d'%(PP,PN,NP,NN)])
mut_stats_df = pd.DataFrame(data = mut_stats,
                            columns = ['Cluster','Mutated gene','Odds ratio','p-value','ClustMut:ClustNoMut:OtherMut:OtherNoMut'])
print()

##########
## Plot ##
##########
plotly_plots = []
plot_titles = []
plot_sizes = []

# Combine cytogenetic and mutation matrices
cyt_row_sums = np.sum(CYT_MAT_BIN, axis = 1)
cyt_sorti = np.argsort(cyt_row_sums)[::-1]
CYT_MAT_BIN = CYT_MAT_BIN[cyt_sorti, :]
cyt_class = np.array(cyt_class)[cyt_sorti]

mut_row_sums = np.sum(MUT_MAT_BIN, axis = 1)
mut_sorti = np.argsort(mut_row_sums)[::-1]
MUT_MAT_BIN = MUT_MAT_BIN[mut_sorti, :]
MUT_MAT_EFF = MUT_MAT_EFF[mut_sorti, :]
GOIG = np.array(GOIG)[mut_sorti]

COMB_MAT_BIN = np.concatenate((CYT_MAT_BIN, MUT_MAT_BIN), axis = 0)

# Sort by individual clusters
sort_inds = []
for clust in clust_vals:
    ci = [i for i,x in enumerate(clust_assign) if x == clust]
    comb_mat_bin_clust = COMB_MAT_BIN[:, ci]
    sorti = np.lexsort(comb_mat_bin_clust[::-1, :])[::-1]
    sort_inds += [ci[si] for si in sorti]
CYT_MAT_BIN = CYT_MAT_BIN[:, sort_inds]
MUT_MAT_BIN = MUT_MAT_BIN[:, sort_inds]
MUT_MAT_EFF = MUT_MAT_EFF[:, sort_inds]
clust_samples = clust_samples[sort_inds]
clust_assign = clust_assign[sort_inds]
CLUST_MAT_INT = np.array([int(x[1]) - 1 for x in clust_assign]).reshape((len(clust_assign), 1)).T
clust_assign = clust_assign.reshape((len(clust_assign), 1)).T

# Create mutation effects matrix with integers
MUT_MAT_EFF_INT = np.zeros(MUT_MAT_EFF.shape)
effect_conv = {x:i+1 for i,x in enumerate(sorted(list(vset)) + ['Multiple'])}
effect_conv['None'] = 0
effects = ['None'] + sorted(list(vset)) + ['Multiple']
for eff in effect_conv:
    inds = np.argwhere(MUT_MAT_EFF == eff)
    for ind in inds:
        i,j = ind
        MUT_MAT_EFF_INT[i, j] = effect_conv[eff]
for i in range(0, MUT_MAT_EFF.shape[0]):
    for j in range(0, MUT_MAT_EFF.shape[1]):
        if ';' in MUT_MAT_EFF[i,j]:
            MUT_MAT_EFF_INT[i, j] = effect_conv['Multiple']

## Cytogenetics
# clusters
cmap = plt.get_cmap('Set2').colors
cmapcolors = ['rgb(%d,%d,%d)' % (tuple([int(round(255*x)) for x in cmap[i]])) for i in range(0, len(cmap))]
colors = []
for ci in range(len(clust_vals)):
    color = cmapcolors[ci]
    colors.append([ci / len(clust_vals), color])
    colors.append([(ci + 1) / len(clust_vals), color])
tickvals = [np.mean([k, k+1]) for k in range(0, len(clust_vals))]
ticktext = [x + ' ' for x in clust_vals]
text = np.array(['Cluster: {}'.format('%s'%(x)) for x in clust_assign.reshape(-1)]).reshape(clust_assign.shape)
maxL = 0
for cc in cyt_class:
    if len(cc) > maxL: maxL = len(cc)
ylab = ' '.join(['']*(maxL - len('Cluster'))) + 'Cluster'
clhmaply = go.Heatmap(z = CLUST_MAT_INT, x = clust_samples, y = [ylab],
                      colorscale = colors, zmin = 0, zmax = len(clust_vals),
                      colorbar = dict(tickvals = tickvals, ticktext = ticktext, tickfont = {'family':'Courier'}),
                      hovertemplate = 'Sample: %{x}<br>%{text}</br><extra></extra>',
                      text = text)
plotly_plots.append(clhmaply)
plot_titles.append('Cytogenetics')
plot_sizes.append([1000, 250])

# cytogenetics
CMAT_CYT_STR = np.tile(clust_assign, (CYT_MAT_BIN.shape[0], 1))
text = np.array(['inClass: {}<br>Cluster: {}'.format('%d'%(x), '%s'%(y)) \
                 for (x,y) in zip(CYT_MAT_BIN[::-1, :].reshape(-1), CMAT_CYT_STR[::-1, :].reshape(-1))]).reshape(CYT_MAT_BIN.shape)
colors = [[0.0, 'rgb(255,255,255)'], [0.5, 'rgb(255,255,255)'], [0.5, 'rgb(0,0,0)'], [1.0, 'rgb(0,0,0)']]
tickvals = [0.25, 0.75]
ticktext = ['No','Yes']
chmaply = go.Heatmap(z = CYT_MAT_BIN[::-1, :], x = clust_samples, y = cyt_class[::-1],
                     colorscale = colors, zmin = 0,
                     colorbar = dict(tickvals = tickvals, ticktext = ticktext, tickfont = {'family':'Courier'}),
                     hovertemplate = 'Cytogenetic class: %{y}<br>Sample: %{x}</br>%{text}<extra></extra>',
                     text = text)
plotly_plots.append(chmaply)
plot_titles.append('')
plot_sizes.append([1000, 750])

# cytogenic associations
ctablely = go.Table(header = dict(values = ['<b>%s</b>'%(x) for x in list(cyt_stats_df.columns)],
                                  fill_color = 'lightskyblue', line_color = 'darkslategray', align = 'left', font = {'family':'Arial'}),
                    cells = dict(values = [cyt_stats_df['Cluster'], cyt_stats_df['Cytogenic class'], cyt_stats_df['Odds ratio'],
                                           cyt_stats_df['p-value'], cyt_stats_df['ClustCyt:ClustNoCyt:OtherCyt:OtherNoCyt']],
                                 fill_color = 'lavender', line_color = 'darkslategray', align = 'left', font = {'family':'Arial'}),
                    columnwidth = [1, 3, 1, 1, 3])
plotly_plots.append(ctablely)
plot_titles.append('Cytogenic features enriched in clusters (p < 0.05)')
plot_sizes.append([1000, 500])

## Mutations
# clusters
cmap = plt.get_cmap('Set2').colors
cmapcolors = ['rgb(%d,%d,%d)' % (tuple([int(round(255*x)) for x in cmap[i]])) for i in range(0, len(cmap))]
colors = []
for ci in range(len(clust_vals)):
    color = cmapcolors[ci]
    colors.append([ci / len(clust_vals), color])
    colors.append([(ci + 1) / len(clust_vals), color])
tickvals = [np.mean([k, k+1]) for k in range(0, len(clust_vals))]
#ticktext = clust_vals
maxXL = 0
for eff in effects:
    if len(eff) > maxXL: maxXL = len(eff)
ticktext = []
for cv in clust_vals:
    lab = cv + ' '.join(['']*(maxXL - len(cv)))
    ticktext.append(lab)
text = np.array(['Cluster: {}'.format('%s'%(x)) for x in clust_assign.reshape(-1)]).reshape(clust_assign.shape)
maxYL = 0
for g in GOIG:
    if len(g) > maxYL: maxYL = len(g)
ylab = ' '.join(['']*(maxYL - len('Cluster'))) + 'Cluster'
clhmaply = go.Heatmap(z = CLUST_MAT_INT, x = clust_samples, y = [ylab],
                      colorscale = colors, zmin = 0, zmax = len(clust_vals),
                      colorbar = dict(tickvals = tickvals, ticktext = ticktext, tickfont = {'family':'Courier'}),
                      hovertemplate = 'Sample: %{x}<br>%{text}</br><extra></extra>',
                      text = text)
plotly_plots.append(clhmaply)
plot_titles.append('Mutations')
plot_sizes.append([1000, 250])

# mutations
cmap = plt.get_cmap('tab20').colors
cmapcolors = ['rgb(255, 255, 255)'] + ['rgb(%d,%d,%d)' % (tuple([int(round(255*x)) for x in cmap[i]])) for i in range(0, len(cmap))]
colors = []
num_eff = len(effect_conv)
for ei in range(0, num_eff):
    color = cmapcolors[ei]
    colors.append([ei / num_eff, color])
    colors.append([(ei + 1) / num_eff, color])
tickvals = [np.mean([k, k+1]) for k in range(0, num_eff)]
ticktext = effects

CMAT_CYT_STR = np.tile(clust_assign, (MUT_MAT_EFF.shape[0], 1))
text = np.array(['Mutations: {}<br>Cluster: {}'.format('%s'%(x), '%s'%(y)) \
                 for (x,y) in zip(MUT_MAT_EFF[::-1, :].reshape(-1), CMAT_CYT_STR[::-1, :].reshape(-1))]).reshape(MUT_MAT_EFF.shape)
mhmaply = go.Heatmap(z = MUT_MAT_EFF_INT[::-1, :], x = clust_samples, y = GOIG[::-1],
                     colorscale = colors, zmin = 0, zmax = num_eff,
                     colorbar = dict(tickvals = tickvals, ticktext = ticktext, tickfont = {'family':'Courier'}),
                     hovertemplate = 'Gene: %{y}<br>Sample: %{x}</br>%{text}</br><extra></extra>',
                     text = text)
plotly_plots.append(mhmaply)
plot_titles.append('')
plot_sizes.append([1000, 750])

# mutation associations table
mtablely = go.Table(header = dict(values = ['<b>%s</b>'%(x) for x in list(mut_stats_df.columns)],
                                  fill_color = 'lightskyblue', line_color = 'darkslategray', align = 'left', font = {'family':'Arial'}),
                    cells = dict(values = [mut_stats_df['Cluster'], mut_stats_df['Mutated gene'], mut_stats_df['Odds ratio'],
                                           mut_stats_df['p-value'], mut_stats_df['ClustMut:ClustNoMut:OtherMut:OtherNoMut']],
                                 fill_color = 'lavender', line_color = 'darkslategray', align = 'left', font = {'family':'Arial'}),
                    columnwidth = [1, 3, 1, 1, 3])
plotly_plots.append(mtablely)
plot_titles.append('Gene mutations enriched in clusters (p < 0.05)')
plot_sizes.append([1000, 500])

## Create final plot(s)
plot_fn = 'AML_exp_subtype_features_heatmaps.html'
dashboard = open(plot_fn, 'w')
dashboard.write("<html><head></head><body>" + "\n")
for pi, pplot in enumerate(plotly_plots):
    fig = go.Figure(data = pplot)
    showxtick = False
    if pi % 3 == 1: showxtick = True
    layout = go.Layout(plot_bgcolor='#777777',
                       xaxis = {'showgrid':False, 'tickangle':-90, 'showticklabels':showxtick, 'tickfont':{'family':'Arial'}},
                       yaxis = {'showgrid':False, 'side':'left', 'tickfont':{'family':'Courier'}},
                       title = plot_titles[pi], width = plot_sizes[pi][0], height = plot_sizes[pi][1],
                       font = {'family':'Arial'})
    fig['layout'].update(layout)
    inner_html = fig.to_html(include_plotlyjs='cdn').split('<body>')[1].split('</body>')[0]
    dashboard.write(inner_html)
dashboard.write("</body></html>" + "\n")
dashboard.close()




