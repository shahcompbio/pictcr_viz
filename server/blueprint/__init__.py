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
from flask import Blueprint
blueprint = Blueprint('blueprint', __name__)
from flask import current_app

config = dict(dotenv_values(".env"))
port = "3004"
MAX_RESPONSE_LENGTH = 10000
COLUMNS = ['IR_VDJ_1_junction_aa','cell_type']


def addHeaders(response):
    response.headers.add('Access-Control-Allow-Credentials', 'true')
    response.headers.add("Access-Control-Allow-Origin", "http://localhost:"+port)
    return response

@blueprint.route('/getFilterValuesByType/<type>-<int:rangestart>-<int:rangeend>')
def getFilters(type,rangestart, rangeend):
    response = ""
    filters = current_app.config['filters']
    for currfilter in filters:
        if currfilter["name"] == type:
            filters = currfilter["values"][rangestart:rangeend]
            print(filters)
            response = jsonify({"filters": filters})
    if len(filters) == 0:
        response = jsonify({"error": "filter type or range does not exist"})

    response = addHeaders(response)

    return response


@blueprint.route('/getStreamData/')
def getStreamData():
    adata = current_app.config["adata"]
    metadata = current_app.config['metadata']
    print(len(metadata))
    def generate():

        #d = list(adata.obs["IR_VDJ_1_junction_aa"].unique().sort_values(ascending=False))[:20]
        #b = list(adata.obs["IR_VDJ_1_junction_aa"].unique().sort_values(ascending=False))[21:100]


        for i in range(0, len(metadata), MAX_RESPONSE_LENGTH):
            yield json.dumps(metadata[i:i + MAX_RESPONSE_LENGTH])+"&/"

        #for element in newList:
        #    print("hello - again")
        #    yield json.dumps(list(element))
    #flask.Response.headers.add('Access-Control-Allow-Credentials', 'true')
    #flask.Response.headers.add("Access-Control-Allow-Origin", "http://localhost:"+port)
    headers = {"Access-Control-Allow-Credentials": "true", "Access-Control-Allow-Origin":"http://localhost:"+port}
    return Response(stream_with_context(generate()),headers=headers, status="200")

@blueprint.route('/getStats/')
def getStats():
    print(current_app)
    stats = current_app.config['stats']
    response = jsonify({"stats": stats})
    response = addHeaders(response)
    return response

@blueprint.route('/getFilterData/')
def getFilterData():
    filters = current_app.config['filters']
    response = jsonify({"filters": filters})
    response = addHeaders(response)
    return response

@blueprint.route('/getData/')
def getData():
    print(current_app)
    metadata = current_app.config['metadata']
    filters = current_app.config['filters']
    stats = current_app.config['stats']

    response = jsonify({"stats":stats, "metadata":metadata,"filters": filters})
    response = addHeaders(response)

    return response

@blueprint.route('/ttest/', methods=['POST'])
def ttest():
    data = json.loads(request.data)
    included = data["data"].split(",")
    adata = current_app.config['adata']
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


    response = addHeaders(response)
    return response
