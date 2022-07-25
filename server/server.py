from flask import Flask, session, jsonify, request, make_response, stream_with_context, Response, g, session
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
from itertools import islice
from dotenv import dotenv_values

from multiprocessing.managers import BaseManager
import flask
from server.blueprint import blueprint

config = dict(dotenv_values(".env"))
port = "3006"
MAX_RESPONSE_LENGTH = 10000
#clone_id = 'IR_VDJ_1_junction_aa'
COLUMNS = ['IR_VDJ_1_junction_aa','cell_type']

def open_data(filepath, phenotype, clone_id):
    adata = sc.read(filepath)
    #assert clone_id in adata.obs.columns, 'Missing field in obs: ' + clone_id
    #assert phenotype in adata.obs.columns, 'Missing field in obs: ' + phenotype
    return adata

def get_filter(adata, range=[]):
    add_columns = list(adata.uns['viz_columns']) if 'viz_columns' in adata.uns else []
    columns = add_columns + COLUMNS
    records = []

    for column in columns:
        record = {
            "name": column if column != "phenotype" else "subtype",
            "values": list(adata.obs[column].unique().sort_values(ascending=False))[:20]
        }
        records.append(record)
    return records

def configure_app(path, project_id, phenotype, clone_id, app):
    print(clone_id)
    adata = open_data(path, phenotype, clone_id)

    #replace NAN
    adata.obs[clone_id] = adata.obs[clone_id].cat.add_categories('NA')
    adata.obs[clone_id].fillna('NA', inplace =True)

    #get different params
    olga = get_olga(adata, clone_id)
    metadata = get_metadata(adata).to_dict(orient="records")
    filters = get_filter(adata)
    stats = get_stats(adata, clone_id)

    app.config.update(
        project_id=project_id,
        adata=adata,
        metadata=metadata,
        filters=filters,
        stats=stats,
        phenotype=phenotype,
        clone_id=clone_id
    )

def get_olga(adata, clone_id):
    data = adata.obs[clone_id]
    df = pd.DataFrame(data).dropna().reset_index(drop=True)
    df.to_csv("./data/example_seqs.tsv", sep="\t")
    return adata


def get_metadata(adata):
    umap = pd.DataFrame(adata.obsm['X_umap'])
    umap.columns = ['UMAP_1', 'UMAP_2']
    add_columns = list(adata.uns['viz_columns']) if 'viz_columns' in adata.uns else []
    df = adata.obs[COLUMNS + add_columns]
    df = df.reset_index()
    df = df.merge(umap, left_index=True, right_index=True)
    df = df.rename(columns={'index': 'cell_id',
                            'pgen': 'log10_probability', 'phenotype': 'subtype'})
    df = df.replace(to_replace="nan", value="None")

    return df

def clean_record(record):
    floats = [field for field in record if isinstance(record[field], float)]
    for field in floats:
        if np.isnan(record[field]):
            record[field] = "None"

    return record

def get_stats(adata, clone_id):
    df = adata.obs[COLUMNS]
    count = df.groupby([clone_id]).size().sort_values(ascending=False).head(11)
    return {clone_id+"-stats": count.to_json()}

def create_app(config=config):
    print(config)
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

    app.register_blueprint(blueprint)
    print(app.url_map)
    server_session = Session(app)

    configure_app(config["path"], config["project_id"], config["phenotype"], config["clone_id"],app)

    print("created app with config")
    print(config)

    return app

app = create_app()
