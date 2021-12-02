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
app.config['SESSION_COOKIE_NAME'] = "permanent_cookie_1"
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=False)
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=5)

app.secret_key = "test"
server_session = Session(app)

@app.route('/isLoaded/')
def isLoaded():
    response = jsonify({"data":False})
    if 'gene_matrix' in session:
        data = session['gene_matrix']

        response = jsonify({"data":True})
        #response = jsonify(data[0:100].to_dict(orient="record"))

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    #response.headers.add('Access-Control-Allow-Origin', 'http://localhost:3000')
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

@app.route('/testing/<path:included>/')
def testing(included=None):
    if included:
        included = included.split(",")

        adata = session['gene_matrix']

        adata.obs["included"] = adata.obs.index.isin(included)
        adata.obs['included'] = adata.obs.apply(lambda x: "included" if x.included else "excluded", axis=1)

        sc.tl.rank_genes_groups(adata,"included")


        ## merge all records names into one large df

        genes = pd.DataFrame(adata.uns['rank_genes_groups']["names"])
        genes = pd.DataFrame(data=genes['excluded'].append(genes['included']), columns=['gene']).reset_index(drop=True)

        adjpvals = pd.DataFrame(adata.uns['rank_genes_groups']["pvals_adj"])
        adjpvals = pd.DataFrame(data=adjpvals['excluded'].append(adjpvals['included']), columns=['adjpvals']).reset_index(drop=True)

        logfc = pd.DataFrame(adata.uns['rank_genes_groups']["logfoldchanges"])
        logfc = pd.DataFrame(data=logfc['excluded'].append(logfc['included']), columns=['logfc']).reset_index(drop=True)

        df = logfc.merge(genes, left_index=True, right_index=True).merge(adjpvals, left_index=True, right_index=True)

        included_df = adata[adata.obs["included"]=="included"]
        included_total_exp = included_df.X.sum(axis=0)
        included_nonzero_counts = (included_df.X != 0).sum(0)
        included_means = included_total_exp / included_nonzero_counts

        excluded_df = adata[adata.obs["included"]=="excluded"]
        excluded_total_exp = excluded_df.X.sum(axis=0)
        excluded_nonzero_counts = (excluded_df.X != 0).sum(0)
        excluded_means = excluded_total_exp / excluded_nonzero_counts

        df['mean_included'] = df.apply(lambda x: included_means[0, included_df.var.index.tolist().index(x.gene)] , axis=1)
        df['mean_included'] = df['mean_included'].fillna(0)
        df['mean_excluded'] = df.apply(lambda x: excluded_means[0, excluded_df.var.index.tolist().index(x.gene)] , axis=1)
        df['mean_excluded'] = df['mean_excluded'].fillna(0)

        df['fc_str'] = df.apply(lambda x: "*+" if x['mean_included'] < x['mean_excluded'] else "*-", axis=1)
        df['logfc'] = df['logfc'].fillna(df['fc_str'])

        final_df = df[['gene','mean_included', 'mean_excluded','logfc', 'adjpvals']].rename(columns={'mean_included': 'includedMean', 'mean_excluded': 'excludedMean', 'logfc': 'fc', 'adjpvals': 'p'})
        final_df = final_df[final_df['p'] > 0]
        final_df = final_df[final_df['includedMean'] > 0]

        final = final_df.to_dict(orient='records')

        final = sorted(final, key=lambda x: x['p'])
        response = jsonify({"data":final})
    else:
        response = jsonify({"data":"None"})
        print("none")

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    #response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")
    return response

@app.route('/l/<path:directory>/')
def load(directory=None):
    # given directory, load data
    response = jsonify({"data":False})

    if directory:
        adata = sc.read("/" + directory)
        session['gene_matrix'] = adata
        response = jsonify({"data":True})

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin",request.environ['HTTP_ORIGIN'])
    #response.headers.add("Access-Control-Allow-Origin", "http://localhost:3000")

    return response
