from flask import Flask, session, jsonify, request, make_response
from flask_session import Session
from datetime import timedelta
import redis
import sys
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, find
import statistics
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smt
from flask_cors import CORS, cross_origin
import anndata
import ast
import json
from flask_cors import CORS

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_PERMANENT'] = True
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_COOKIE_NAME'] = "permanent_cookie_2"
#app.config['SESSION_REDIS'] = redis.from_url('redis://pic-redis:6379')
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=False)
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=6)

app.secret_key = "test2"
server_session = Session(app)

COLUMNS = ['phenotype', 'clone_id']


adataSantosh = sc.read("./data/newfile.h5ad")
@app.route('/getSantoshData/<path:file>/')
def getSantoshData(file=None):
        if file:
            metadata = get_santosh_metadata(adataSantosh).to_dict(orient="records")
            filters = get_santosh_filter(adataSantosh)

            response = jsonify({"metadata":metadata,"filters": filters})
            response.headers.add('Access-Control-Allow-Credentials', 'true')
            response.headers.add("Access-Control-Allow-Origin", "http://localhost:3001")
            return response

def get_santosh_metadata(adata):
    umap = pd.DataFrame(adata.obsm['X_umap'])
    umap.columns = ['UMAP_1', 'UMAP_2']
    add_columns = ['response', 'patient', 'timepoint', 'treatment','cell_type']

    df = adata.obs[add_columns]
    df = df.reset_index()
    df = df.merge(umap, left_index=True, right_index=True)
    df = df.rename(columns={'index': 'cell_id',
                            'pgen': 'log10_probability', 'cell_type': 'subtype'})
    df = df.replace(to_replace="nan", value="None")
    return df

def get_santosh_filter(adata):
    columns = ['response', 'patient', 'timepoint', 'treatment','cell_type']
    records = []

    for column in columns:
        record = {
            "name": column if column != "phenotype" else "subtype",
            "values": list(adata.obs[column].unique())
        }
        records.append(record)

    return records

@app.route('/ttestSantosh/', methods=['POST'])
def ttestSantosh():
    data = json.loads(request.data)
    included = data["data"].split(",")

    adataSantosh.obs["included"] = adataSantosh.obs.index.isin(included)
    adataSantosh.obs['included'] = adataSantosh.obs.apply(lambda x: "included" if x.included else "excluded", axis=1)

    sc.tl.rank_genes_groups(adataSantosh,"included")

    genes = pd.DataFrame(adataSantosh.uns['rank_genes_groups']["names"])[['included']].rename(columns={'included': 'gene'})
    adjpvals = pd.DataFrame(adataSantosh.uns['rank_genes_groups']["pvals_adj"])[['included']].rename(columns={'included': 'p'})
    logfc = pd.DataFrame(adataSantosh.uns['rank_genes_groups']["logfoldchanges"])[['included']].rename(columns={'included': 'fc'})
    df = logfc.merge(genes, left_index=True, right_index=True).merge(adjpvals, left_index=True, right_index=True)
    df = df[df['p'] < 0.05]
    df = df[df['fc'] > 0.25]
    final = df.to_dict(orient='records')

    final = sorted(final, key=lambda x: x['p'])
    response = jsonify({"data":final})

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3001")
    return response


#adata = sc.read("../data/viz.h5ad")
adata = sc.read("./data/viz.h5ad")
@app.route('/getData/<path:file>/')
def getMetadata(file=None):
        if file:
            #adata = sc.read("/"+file)
            metadata = get_metadata(adata).to_dict(orient="records")
            filters = get_filter(adata)

            response = jsonify({"metadata":metadata,"filters": filters})
            response.headers.add('Access-Control-Allow-Credentials', 'true')
            #response.headers.add("Access-Control-Allow-Origin", "https://spectrum-staging.shahlab.mskcc.org")
            response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")

            return response

def get_filter(adata):
    columns = list(adata.uns['viz_columns']) + COLUMNS
    records = []

    for column in columns:
        record = {
            "name": column if column != "phenotype" else "subtype",
            "values": list(adata.obs[column].unique())
        }
        records.append(record)

    return records

def clean_record(record):
    floats = [field for field in record if isinstance(record[field], float)]
    for field in floats:
        if np.isnan(record[field]):
            record[field] = "None"

    return record

def get_metadata(adata):
    umap = pd.DataFrame(adata.obsm['X_umap'])
    umap.columns = ['UMAP_1', 'UMAP_2']
    add_columns = list(adata.uns['viz_columns'])
    df = adata.obs[COLUMNS + add_columns + ['pgen']]
    df = df.reset_index()
    df = df.merge(umap, left_index=True, right_index=True)
    df = df.rename(columns={'index': 'cell_id',
                            'pgen': 'log10_probability', 'phenotype': 'subtype'})
    df = df.replace(to_replace="nan", value="None")
    return df

@app.route('/ttest/', methods=['POST'])
def ttest():
    data = json.loads(request.data)
    included = data["data"].split(",")
    print(request.environ['HTTP_ORIGIN'])
    adata.obs["included"] = adata.obs.index.isin(included)
    adata.obs['included'] = adata.obs.apply(lambda x: "included" if x.included else "excluded", axis=1)

    sc.tl.rank_genes_groups(adata,"included")

    genes = pd.DataFrame(adata.uns['rank_genes_groups']["names"])[['included']].rename(columns={'included': 'gene'})
    adjpvals = pd.DataFrame(adata.uns['rank_genes_groups']["pvals_adj"])[['included']].rename(columns={'included': 'p'})
    logfc = pd.DataFrame(adata.uns['rank_genes_groups']["logfoldchanges"])[['included']].rename(columns={'included': 'fc'})
    df = logfc.merge(genes, left_index=True, right_index=True).merge(adjpvals, left_index=True, right_index=True)
    df = df[df['p'] < 0.05]
    df = df[df['fc'] > 0.25]
    final = df.to_dict(orient='records')

    final = sorted(final, key=lambda x: x['p'])
    response = jsonify({"data":final})


    response.headers.add('Access-Control-Allow-Credentials', 'true')
    #response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")
    return response

@app.route('/l/<path:directory>/')
def load(directory=None):
    response = jsonify({"data":False})
    print(request.environ['HTTP_ORIGIN'])
    if directory:
        adata = sc.read("/" + directory)
        session['gene_matrix'] = adata
        response = jsonify({"data":True})

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    #response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")

    return response

#if __name__ == "__main__":
#    app.run(debug=True, host='0.0.0.0')
#if __name__ == "__main__":
#    app.run(
#        host=os.environ.get("BACKEND_HOST", "172.0.0.1"),
#        port=5000,
#        debug=True,
#    )
