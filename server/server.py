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
import argparse
from dotenv import dotenv_values

config = dict(dotenv_values(".env"))

COLUMNS = ['phenotype', 'clone_id']

def open_data(filepath, timepoint, phenotype, clone_id):
    adata = sc.read(filepath)
    #assert clone_id in adata.obs.columns, 'Missing field in obs: ' + clone_id
    #assert timepoint in adata.obs.columns, 'Missing field in obs: ' + timepoint
    #assert phenotype in adata.obs.columns, 'Missing field in obs: ' + phenotype
    return adata

def configure_app(path, project_id, timepoint, phenotype, clone, app):
    adata = open_data(path, timepoint, phenotype, clone)
    metadata = get_metadata(adata).to_dict(orient="records")
    filters = get_filter(adata)

    app.config.update(
        project_id=project_id,
        adata=adata,
        metadata=metadata,
        filters=filters,
        timepoint=timepoint,
        phenotype=phenotype,
        clone=clone
    )

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

def create_app(config=config):
    print("creating application with env file")

    app = Flask(__name__)

    CORS(app, resources={r"/*": {"origins": "*"}})
    app.config['SESSION_TYPE'] = 'redis'
    app.config['SESSION_PERMANENT'] = True
    app.config['SESSION_USE_SIGNER'] = True
    app.config['SESSION_COOKIE_NAME'] = "permanent_cookie_2"
    app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
    app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=False)
    app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=6)
    app.secret_key = "test2"

    server_session = Session(app)

    configure_app(config["path"], config["project_id"], config["timepoint"], config["phenotype"], config["clone"],app)

    print("created app with config")
    print(config)

    return app

app = create_app()

@app.route('/getData/')
def getData():
    metadata = app.config['metadata']
    filters = app.config['filters']

    response = jsonify({"metadata":metadata,"filters": filters})
    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")

    return response


def clean_record(record):
    floats = [field for field in record if isinstance(record[field], float)]
    for field in floats:
        if np.isnan(record[field]):
            record[field] = "None"

    return record


@app.route('/ttest/', methods=['POST'])
def ttest():
    data = json.loads(request.data)
    included = data["data"].split(",")
    adata = app.config['adata']
    print(request.environ)
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
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")
    return response

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Output HTML with data file')
    parser.add_argument('dashboard_id', type=str)
    parser.add_argument('project_id', default="pic", help='project settings variant')
    parser.add_argument('path', type=str)
    parser.add_argument('--width', type=int, default=700, help='Pixel width of sankey plot')
    parser.add_argument('--height', type=int, default=600, help='Pixel height of sankey plot')
    parser.add_argument('--timepoint', type=str,action="store", default='tp', help='Column name for timepoint')
    parser.add_argument('--clone', type=str, action="store",default='trb', help='Column name for clone ID')
    parser.add_argument('--phenotype', type=str,action="store", default='phenotype', help='Column name for phenotype')

    args = parser.parse_args()
    print(args.timepoint)
    configure_app(args.path, args.project_id, args.timepoint, args.phenotype, args.clone)

    app.run()
