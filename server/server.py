from flask import Flask, session, jsonify, request, make_response
from flask_session import Session
from datetime import timedelta
import redis
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, find
import statistics
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smt
from flask_cors import CORS, cross_origin

from flask_cors import CORS

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})
app.config['SESSION_TYPE'] = 'redis'
app.config['SESSION_PERMANENT'] = True
app.config['SESSION_USE_SIGNER'] = True
app.config['SESSION_COOKIE_NAME'] = "permanent_cookie"
app.config['SESSION_REDIS'] = redis.from_url('redis://localhost:6379')
app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=True)
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
    response.headers.add('Access-Control-Allow-Origin', 'http://localhost:3001')
    return response


@app.route('/table/')
def table():
    data = session['gene_matrix']
    response = jsonify(data[0:100].to_dict(orient="record"))

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add('Access-Control-Allow-Origin', 'http://localhost:3001')

    return response

@app.route('/testing/<path:included>/')
def testing(included=None):
    if included:
        included = included.split(",")
        adata = session['gene_matrix']
        gene_matrix = csr_matrix(adata.X)
        nonzero = find(gene_matrix)

        inc = []
        for d in adata.obs.index.tolist():
            if  d in included:
                inc.append("included")
            else:
                inc.append("excluded")
        adata.obs["included"] = inc

        sc.tl.rank_genes_groups(adata,"included")

        logfc = list(zip(*adata.uns['rank_genes_groups']["logfoldchanges"].tolist()))
        genes = list(zip(*adata.uns['rank_genes_groups']["names"].tolist()))
        adjpvals = list(zip(*adata.uns['rank_genes_groups']["pvals_adj"].tolist()))
        order = ["included","excluded"]


        timepoint = sc.tl.rank_genes_groups(adata,"included")

        final = []

        for cond,fcs,gene,pvalues in zip(order,logfc,genes,adjpvals):
            for fc, g, pval in zip(fcs,gene,pvalues):
                pre = adata[adata.obs["included"]=="included"]
                post = adata[adata.obs["included"]=="excluded"]
                premean = np.mean(pre.X[:,pre.var.index.tolist().index(g)])
                postmean = np.mean(post.X[:,post.var.index.tolist().index(g)])
                if str(fc) == "nan":
                    if premean < postmean:
                        fc = "*+"
                    if premean > postmean:
                        fc = "*-"
                if (pval != 0):
                    final.append({"gene":g, "includedMean":str(premean),"excludedMean":str(postmean),"fc":str(fc),"p":str(pval)})
        final = sorted(final, key=lambda x: float(x['p']))
        response = jsonify({"data":final})
    else:
        response = jsonify({"data":"None"})
        print("none")

    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3001")
    return response

@app.route('/load/<path:directory>/')
def load(directory=None):
    # given directory, load data
    response = jsonify({"data":False})
    if directory:
        adata = sc.read("/" + directory)
        session['gene_matrix'] = adata
        response = jsonify({"data":True})


    response.headers.add("Access-Control-Allow-Credentials", "true")
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:3001")

    return response
