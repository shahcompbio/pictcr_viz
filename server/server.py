from flask import Flask, session, jsonify
from flask_session import Session
import redis
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, find

app = Flask(__name__)
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

    return jsonify(data[0:10].to_dict(orient="record"))


@app.route('/ttest/<cell_ids>')
def ttest(cell_ids):
    # data = session['genes']

    return []
