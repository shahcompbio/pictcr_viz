from flask import Flask, session, jsonify, request, make_response, stream_with_context, Response, g, session
from flask_session import Session
from datetime import timedelta
import flask
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

from server.blueprint import blueprint


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
    print("url map:")
    print(app.url_map)
    
    configure_app(args.path, args.project_id, args.timepoint, args.phenotype, args.clone)

    app.run()
