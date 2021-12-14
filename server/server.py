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
from flask_cors import CORS

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

@app.route('/isLoaded/')
def isLoaded():
    response = jsonify({"data":False})
    if 'gene_matrix' in session:
        data = session['gene_matrix']

        response = jsonify({"data":True})
        #response = jsonify(data[0:100].to_dict(orient="record"))

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    #response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    response.headers.add('Access-Control-Allow-Origin', 'http://localhost:3000')
    return response


@app.route('/table/')
@cross_origin(supports_credentials=True)
def table():
    data = session['gene_matrix']
    response = jsonify(data[0:100].to_dict(orient="record"))

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    #response.headers.add('Access-Control-Allow-Origin', 'http://localhost:3000')

    return response

adata = sc.read("")
@app.route('/testing/<path:included>/')
def testing(included=None):
    if included:
        included = included.split(",")

        #adata = session['gene_matrix']

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
    else:
        response = jsonify({"data":"None"})
        print("none")

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    #response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")
    return response

@app.route('/l/<path:directory>/')
def load(directory=None):
    # given directory, load data
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
