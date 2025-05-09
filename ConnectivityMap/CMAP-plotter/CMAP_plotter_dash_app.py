import sys, os, csv, glob, datetime, re, io
import pandas as pd
import numpy as np
from cmapPy.pandasGEXpress.parse import parse

import platform
from base64 import b64encode

from dash import Dash, dash_table, html, dcc
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import plotly
import plotly.graph_objects as go

# initialize app
external_stylesheets = [dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME]
app = Dash(__name__, external_stylesheets=external_stylesheets)

## Data files
# GCTx file
DIR = 'data/' # Main path to CMAP data files
level5_gctx = DIR + 'level5_beta_trt_cp_n720216x12328.gctx'

# Load meta data frames
cell_info = pd.read_csv(DIR + 'cellinfo_beta.txt', sep = '\t')
comp_info = pd.read_csv(DIR + 'compoundinfo_beta.txt', sep = '\t')
gene_info = pd.read_csv(DIR + 'geneinfo_beta.txt', sep = '\t')
sig_info = pd.read_csv(DIR + 'siginfo_beta.txt', sep = '\t')

COMPOUNDS = sorted(list(set(['%s:%s' % (a, b) for a, b in zip(list(comp_info['pert_id']), list(comp_info['cmap_name']))])))
gene_conv = {}
GENES = []
for (gid, gsym) in zip(list(gene_info['gene_id']), list(gene_info['gene_symbol'])):
    gene_conv[str(gid)] = gsym
    GENES.append(gsym)
landmark_gids = [str(x) for x in list(gene_info.loc[gene_info['feature_space']=='landmark']['gene_id'])]

## Components
title1 = html.H1('Connectivity Map expression data plotting tool',
                 style = {'font_family': 'Arial', 'color':'blue'})
button_plot = dbc.Button("Plot line graphs", id = "gene-plots-button", color = 'success', className = 'me-1',
                         style={"marginBottom": 10, 'marginTop':10})
button_hmap = dbc.Button("Create Heatmap", id = "heatmap-button", color = 'danger', className = 'me-1',
                         style = {'marginBottom': 10, 'marginTop':10})
buttons = [button_plot, button_hmap]

compound_select = html.Div([
    html.P(),
    html.Div([
    dcc.Dropdown(
        id = 'select-which-compound',
        options = [{'label':i, 'value': i} for i in COMPOUNDS],
        placeholder = 'Select compound')]),
])

inter_gene =  html.Div(id = 'interactivity-gene')
inter_gplot = html.Div(id = 'interactivity-plot')
inter_hmap =  html.Div(id = 'interactivity-hmap')
inters = [inter_gene, inter_gplot, inter_hmap]

## Set up app layout
app_items = [title1, compound_select] + buttons + inters
app.layout = dbc.Container(app_items)

## CALLBACKS
@app.callback(
    Output('interactivity-gene', 'children'),
    State("select-which-compound", 'value'),
    Input("gene-plots-button", 'n_clicks'),
    prevent_initial_call = True)
def get_gene_options(compound, n_clicks):

    if compound == None:
        return
    if n_clicks == None:
        return

    pert_id, pert_name = compound.split(':')
    cl_drug_df = sig_info.loc[(sig_info['pert_id'] == pert_id)][['sig_id','pert_idose','pert_itime','cell_mfc_name']]
    cell_lines = sorted(list(set(list(cl_drug_df['cell_mfc_name']))))

    elements = html.Div([
        html.P(),
        html.Div([
        dcc.Dropdown(id = 'select-which-gene', placeholder = 'Select gene',
                     options = [{'label':gene, 'value': gene} for gene in GENES]),
        dcc.Dropdown(id = 'select-which-cells', placeholder = 'Select cell lines (optional)',
                     options = [{'label':i, 'value':i} for i in cell_lines], multi = True)
        ],
        style = {'width': '100%', 'display': 'inline-block', 'columnCount':2}),
    ])

    return elements

# Resets
#@app.callback(
#    Output("interactivity-plot", 'children', allow_duplicate = True),
#    Input("select-which-gene", 'value'),
#    prevent_initial_call = True)
#def reset_plot(x):
#    return
#@app.callback(
#    Output("interactivity-plot", 'children', allow_duplicate = True),
#    Output('select-which-gene', 'value'),
#    Input("select-which-compound", 'value'),
#    prevent_initial_call = True)
#def reset_plot2(x):
#    return None, None

@app.callback(
    Output('interactivity-plot', 'children'),
    Input('select-which-compound', 'value'),
    Input('select-which-gene', 'value'),
    Input('select-which-cells', 'value'))
def plot_drug_profiles(compound, gene, cell_lines):

    if compound == None:
        return
    if gene == None:
        return

    print(compound, gene)

    gid = [x for x in gene_conv if gene_conv[x] == gene]
    pert_id, pert_name = compound.split(':')
    cl_drug_df = sig_info.loc[(sig_info['pert_id'] == pert_id)][['sig_id','pert_idose','pert_itime','cell_mfc_name']]
    if cell_lines != None:
        if len(cell_lines) > 0:
            cl_drug_df = cl_drug_df.loc[(cl_drug_df['cell_mfc_name'].isin(cell_lines))]
    cl_drug_ids = list(cl_drug_df['sig_id'])
    cell_lines = list(cl_drug_df['cell_mfc_name'])
    doses, times = list(cl_drug_df['pert_idose']), list(cl_drug_df['pert_itime'])

    new_cols = {}
    cn2cid = {}
    for (cl, sid, dose, time, sid) in zip(cell_lines, cl_drug_ids, doses, times, cl_drug_ids):
        nid = '%s:%s:%s:%s>%s' % (cl, pert_name, dose.replace(' ',''), time.replace(' ',''), sid)
        #print(sid, nid)
        new_cols[sid] = nid
    col_ordered0 = sorted([new_cols[x] for x in new_cols], key = lambda x:
                         (x.split(':')[0], float(x.split(':')[2].replace('uM',''))))
    col_ordered = [x.split('>')[-1] for x in col_ordered0]
    #print(col_ordered)

    # Parse data from GCTx
    gcto = parse(level5_gctx, cid = cl_drug_ids)#, rid = gid)
    gcto_df = gcto.data_df
    gcto_df = gcto_df[gcto_df.index.isin(gid)]
    gcto_df = gcto_df.rename(index = gene_conv)
    gcto_df = gcto_df.sort_index()
    gcto_df = gcto_df[col_ordered]
    gcto_df = gcto_df.rename(columns = new_cols)

    print(gcto_df)
    #gcto_dict0 = gcto_df.to_dict(orient = 'list')

    PLOT_DATA = {}
    durations = set()
    for col in gcto_df.columns:
        coll, sid = col.split('>')
        cl, dn, dose, time = coll.split(':')
        durations.add(time)
        if cl not in PLOT_DATA: PLOT_DATA[cl] = {}
        if time not in PLOT_DATA[cl]:
            PLOT_DATA[cl][time] = {}
            PLOT_DATA[cl][time]['X'] = []
            PLOT_DATA[cl][time]['Y'] = []
        PLOT_DATA[cl][time]['X'].append(float(dose.replace('uM','')))
        y = list(gcto_df[col])[0]
        #print(gcto_df[col])
        #print(list(gcto_df[col]))
        #print()
        PLOT_DATA[cl][time]['Y'].append(y)
    durations = sorted(list(durations))
    print(durations)
    t_figs = []
    for t in durations:
        fig = go.Figure()
        fig.update_layout(title = '%s expression change induced by %s after %s' % (gene, pert_name, t),
                          xaxis_title = 'Concentration (uM)', yaxis_title = '%s expression change' % (gene))
        t_figs.append(fig)
    for ti, t in enumerate(durations):
        for cl in PLOT_DATA:
            if t not in PLOT_DATA[cl]: continue
            X = PLOT_DATA[cl][t]['X']
            Y = PLOT_DATA[cl][t]['Y']
            tfmt = 'Cell line: %s<br>Concentration: {}</br>Expression change: {}' % (cl)
            text = [tfmt.format('%s' % (conc), '%0.4f' % (expr)) for (conc, expr) in zip(X, Y)]
            t_figs[ti].add_trace(go.Scatter(
                x = X, y = Y, text = text, hovertemplate = '%{text}<extra></extra>',
                name = cl, mode = 'lines+markers',
            ))

    OUT = []
    for fig in t_figs:
        fig.update_xaxes(type = "log")
        pld = html.Div([
            dcc.Graph(figure = fig),
        ])
        OUT.append(pld)

    return OUT

## Heatmap creation callbacks
@app.callback(
    Output('interactivity-hmap', "children"),
    Input('heatmap-button', 'n_clicks'),
    State("select-which-compound", 'value'),
    prevent_initial_call = True)
def create_hmap_options(n_clicks, compound):

    if compound == None:
        return

    pert_id, pert_name = compound.split(':')
    cl_drug_df = sig_info.loc[(sig_info['pert_id'] == pert_id)][['sig_id','pert_idose','pert_itime','cell_mfc_name']]
    cell_lines = sorted(list(set(list(cl_drug_df['cell_mfc_name']))))

    # Create elements
    elements = html.Div([
        html.P(),
        html.Div([
        dcc.Dropdown(
            id = 'hmap-which-cells',
            options = [{'label': i, 'value': i} for i in cell_lines],
            placeholder = 'Select cell lines (optional)', multi = True),
        dcc.Dropdown(
            id = 'hmap-which-genes',
            options = [{'label': g, 'value': g} for g in GENES],
            placeholder = 'Select genes (optional)', multi = True),
        ],
        style = {'width': '100%', 'display': 'inline-block', 'columnCount':2}),
        html.P(),
        dcc.Input(
            id = 'hmap-height-input', type = 'number', placeholder = 'Set heatmap height', debounce = True),
        dcc.Input(
            id = 'hmap-width-input', type = 'number', placeholder = 'Set heatmap width', debounce = True),
        html.Div([
        dcc.Slider(
            0, 2.5, 0.025,
            id = 'hmap-range-slider',
            value = 1.5,
            marks = {0:'0',0.25:'0.25',0.5:'0.5',0.75:'0.75',1:'1',1.25:'1.25',1.5:'1.5',1.75:'1.75',2:'2',2.25:'2.25',2.5:'2.5'},
            included = False, tooltip={"placement": "bottom", "always_visible": False})],
        style = {'display':'inline-block','width': '30%'}),
        html.P(),
        html.Div([
            dcc.RadioItems(['None','Cluster columns','Cluster rows', 'Cluster both'], 'None', id = 'hmap-cluster-radio'),
            dcc.RadioItems(['All genes','Landmark genes'], 'Landmark genes', id = 'hmap-genes-radio'),
            dcc.RadioItems(['All samples', 'Exemplars'], 'All samples', id = 'hmap-samples-radio'),
        ],
        style = {'display':'inline-block','width': '100%', 'height': '100%','columnCount':3}),
        html.Div(id = 'interactivity-hmap-plot'),
        html.P(),
    ])

    return elements

@app.callback(
    Output('interactivity-hmap-plot', "children"),
    Input('select-which-compound', 'value'),
    Input("hmap-which-cells", 'value'),
    Input("hmap-which-genes", 'value'),
    Input("hmap-range-slider", 'value'),
    Input("hmap-height-input", 'value'),
    Input("hmap-width-input", 'value'),
    Input("hmap-cluster-radio", 'value'),
    Input("hmap-genes-radio", 'value'),
    Input("hmap-samples-radio", 'value'))
def create_heatmap(compound, cell_lines, gene_list, dispr, plotH, plotW, cluster, which_genes, which_samps):

    if compound == None:
        return

    pert_id, pert_name = compound.split(':')
    cl_drug_df = sig_info.loc[(sig_info['pert_id'] == pert_id)][['sig_id','pert_idose','pert_itime','cell_mfc_name','is_exemplar_sig']]
    if cell_lines != None:
        if len(cell_lines) > 0:
            cl_drug_df = cl_drug_df.loc[(cl_drug_df['cell_mfc_name'].isin(cell_lines))]
    if which_samps == 'Exemplars':
        cl_drug_df = cl_drug_df.loc[(cl_drug_df['is_exemplar_sig'] == 1)]
    print(cl_drug_df)
    cl_drug_ids = list(cl_drug_df['sig_id'])
    cell_lines = list(cl_drug_df['cell_mfc_name'])
    doses, times = list(cl_drug_df['pert_idose']), list(cl_drug_df['pert_itime'])

    new_cols = {}
    cn2cid = {}
    for (cl, sid, dose, time, sid) in zip(cell_lines, cl_drug_ids, doses, times, cl_drug_ids):
        nid = '%s:%s:%s:%s>%s' % (cl, pert_name, dose.replace(' ',''), time.replace(' ',''), sid)
        new_cols[sid] = nid
    col_ordered0 = sorted([new_cols[x] for x in new_cols], key = lambda x:
                         (x.split(':')[0], float(x.split(':')[2].replace('uM',''))))
    col_ordered = [x.split('>')[-1] for x in col_ordered0]
    #print(col_ordered)

    # Parse data from GCTx
    gcto = parse(level5_gctx, cid = cl_drug_ids)
    gcto_df = gcto.data_df
    if which_genes == 'Landmark genes':
        gcto_df = gcto_df[gcto_df.index.isin(landmark_gids)]
    if gene_list != None:
        if len(gene_list) > 0:
            gids = [x for x in gene_conv if gene_conv[x] in gene_list]
            gcto_df = gcto_df[gcto_df.index.isin(gids)]
    gcto_df = gcto_df.rename(index = gene_conv)
    gcto_df = gcto_df.sort_index()
    gcto_df = gcto_df[col_ordered]
    gcto_df = gcto_df.rename(columns = new_cols)

    # Create matrix
    mat = gcto_df.to_numpy()
    genes = np.array(gcto_df.index.tolist())
    samples = np.array([new_cols[sid].split('>')[0] for sid in col_ordered])
    samp2idx = {s:i for i,s in enumerate(samples)}
    gene2idx = {g:i for i,g in enumerate(genes)}
    print(mat)
    print(mat.shape)
    print(genes[0:10])

    # Plot dimensions
    width, height = 1500, 1200
    if plotH != None:
        height = max(400, plotH)
    if plotW != None:
        width = max(400, plotW)

    xlabs = samples
    ylabs = genes

    # Check if clustering selected
    dengo = None
    if cluster.startswith('Cluster'):
        import plotly.figure_factory as ff
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import pdist

        if (cluster == 'Cluster columns' or cluster == 'Cluster both') and mat.shape[1] > 1:
            dengo = ff.create_dendrogram(np.nan_to_num(mat, nan = 0.0, posinf = 0.0, neginf = 0.0).T,
                                         orientation = 'bottom', labels = xlabs, #color_threshold = 1.25,
                                         distfun = lambda x: pdist(x, metric = 'euclidean'),
                                         linkagefun = lambda x: linkage(x, method = 'ward'))

            leaves = dengo['layout']['xaxis']['ticktext']
            print(leaves)
            cinds = [samp2idx[s] for s in leaves]
            mat = mat[:, cinds]
            print(mat.shape)
            #text = text[:, cinds]
            xlabs = leaves

            # Attempt to even widths
            max_lab_len, mlab = 0, None
            for lab in ylabs:
                if len(lab) > max_lab_len:
                    max_lab_len = len(lab)
                    mlab = lab
            ytickl = mlab
            dengo['layout']['yaxis']['tickvals'] = [0.0]
            dengo['layout']['yaxis']['ticktext'] = np.array([ytickl])

            dlayout = go.Layout(plot_bgcolor='#FFFFFF', xaxis = {'showgrid':False, 'tickangle':-90, 'ticks':""},
                                yaxis = {'showgrid':False, 'ticks':""},
                                width = width*0.965, height = 300,
                                font = {'family':'Arial'})
            dengo.update_layout(dlayout)
            dengo.update_yaxes(color = 'white') # make dummy label "disappear"
            dbuffer = io.StringIO()
            dengo.write_html(dbuffer, include_plotlyjs = 'cdn')

    if (cluster == 'Cluster rows' or cluster == 'Cluster both') and mat.shape[0] > 1:

        dengo_ = ff.create_dendrogram(np.nan_to_num(mat, nan = 0.0, posinf = 0.0, neginf = 0.0),
                                      orientation = 'bottom', labels = ylabs, #color_threshold = 1.25,
                                      distfun = lambda x: pdist(x, metric = 'euclidean'),
                                      linkagefun = lambda x: linkage(x, method = 'ward'))

        leaves = dengo_['layout']['xaxis']['ticktext']
        rinds = [gene2idx[g] for g in leaves]
        mat = mat[rinds, :]
        #text = text[rinds, :]
        ylabs = leaves

    # Create heatmap
    hmaply = go.Heatmap(z = mat[::-1, :], x = xlabs, y = ylabs[::-1],
                        colorscale='RdBu_r', zmin = -dispr, zmax = dispr,
                        hovertemplate = 'Value: %{z}<br>Sample: %{x}</br>Gene: %{y}<extra></extra>', #<br>%{text}</br><extra></extra>',
                        #text = text[::-1, :],
                        )
    fig = go.Figure(data = hmaply)

    # Adjust size, style
    title = '%s (%s)' % (pert_name, pert_id)
    layout = go.Layout(plot_bgcolor='#777777', xaxis = {'showgrid':False, 'tickangle':-90}, yaxis = {'showgrid':False},
                       title = title, width = width, height = height,
                       font = {'family':'Arial'})
    fig.update_layout(layout)
    buffer = io.StringIO()
    fig.write_html(buffer, include_plotlyjs = 'cdn')

    # Create elements object
    if dengo != None:
        elements = html.Div([
            dcc.Graph(figure = dengo),
            html.A(
                html.Button("Download as HTML"), id = "downloadDEN",
                href = "data:text/html;base64,"+b64encode(dbuffer.getvalue().encode()).decode(),
                download = pert_name + '_' + pert_id + '_dendrogram.html'),
            dcc.Graph(figure = fig),
            html.A(
                html.Button("Download as HTML"), id = "downloadHMAP",
                href = "data:text/html;base64,"+b64encode(buffer.getvalue().encode()).decode(),
                download = pert_name + '_' + pert_id + '_heatmap.html')
        ])
    else:
        elements = html.Div([
            dcc.Graph(figure = fig),
            html.A(
                html.Button("Download as HTML"), id = "downloadHMAP",
                href = "data:text/html;base64,"+b64encode(buffer.getvalue().encode()).decode(),
                download = pert_name + '_' + pert_id + '_heatmap.html')
        ])

    # return elements object
    return elements


### MAIN
if __name__ == '__main__':
    app.run(debug=False, dev_tools_silence_routes_logging = True, port = 8051)



