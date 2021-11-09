from flask import Flask, session, jsonify, request
from flask_session import Session
import redis
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, find
import statistics
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smt
from flask_cors import CORS, cross_origin

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_PERMANENT'] = False
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
server_session = Session(app)

app.secret_key = "pictctrtest"


@app.route('/')
def hello_world():
    return 'Hello World!'


@app.route('/load/<path:directory>/')
def load(directory=None):
    # given directory, load data
    session.clear()
    if directory:
        adata = sc.read("/" + directory)
        sc.tl.umap(adata)

        gene_matrix = csr_matrix(adata.X)
        nonzero = find(gene_matrix)

        gene_matrix = pd.DataFrame(
            list(zip(*nonzero)), columns=['cell_idx', 'gene_idx', 'expression'])

        genes = adata.var.reset_index().rename(columns={'index': 'gene_id'})
        gene_matrix = gene_matrix.merge(
            genes[['gene_id']], left_on='gene_idx', right_index=True)

        cells = adata.obs.reset_index().rename(columns={'index': 'cell_id'})

        gene_matrix = gene_matrix.merge(
            cells[['cell_id']], left_on='cell_idx', right_index=True)

        session['gene_matrix'] = gene_matrix

        # data = load_qc_data("/" + directory)
        # session['bins'] = get_bins_data(data)
        # gc_bias = get_gc_bias_data(data)
        # segs = get_segs_data(data)
        # qc = get_qc_data(data)

        # grouped_segs = segs.groupby(['id'])

        # qc_records = []
        # for record in qc.to_dict(orient="records"):
        #     clean_nans(record)
        #     cell_segs = []
        #     for seg_record in grouped_segs.get_group(record['id']).to_dict(orient='record'):
        #         clean_nans(seg_record)
        #         cell_segs.append(seg_record)

        #     record['segs'] = cell_segs
        #     qc_records.append(record)

        return 'Loaded'


@app.route('/test')
def test():
    data = session['gene_matrix']
    lasso_ids = ["Pre_P1_C1_P4_M11"]
    lasso_dict = { cell_id : True for cell_id in lasso_ids }

    included = []
    excluded = []
    all_data = data.to_dict(orient="record")

    for d in all_data:
        if  d["cell_id"] in lasso_dict:
            included.append(d)
        else:
            excluded.append(d)

    included_exp = [ exp["expression"] for exp in included]
    excluded_exp = [ exp["expression"] for exp in excluded]
    mean_included = statistics.mean(included_exp)
    mean_excluded = statistics.mean(excluded_exp)

    t = ttest_ind(mean_included,mean_excluded)

    p_val = t.pvalue
    response = jsonify({"p": p_val, "mean_included": mean_included, "mean_excluded": mean_excluded, "table_included":included })
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response

@app.route('/ttest/<cell_ids>')
def ttest(cell_ids):
    # data = session['genes']
    return []
