from jinja2 import Template
import scanpy as sc
import numpy as np
import pandas as pd
import math
import os
import sys
import json

# import subprocess
# import pandas as pd
# import sys
# import os
# import json
# import shutil

PTHRESH = 0.01
FCTHRESH = 0.01
TOP_N = 20

def get_data(filepath):
    adata = sc.read(filepath)
    metadata = get_metadata(adata).to_dict(orient="records")
    degs = get_degs(adata).to_dict(orient="records")
    probabilities = get_probabilities(adata).to_dict(orient="records")

    return {
        "metadata": [clean_record(record) for record in metadata],
        "degs": [clean_record(record) for record in degs],
        "probabilities": [clean_record(record) for record in probabilities]
    }


def clean_record(record):
    floats = [field for field in record if isinstance(record[field], float)]
    for field in floats:
        if np.isnan(record[field]):
            record[field] = "None" 
    
    return record

def get_metadata(adata):
    umap = pd.DataFrame(adata.obsm['X_umap'])
    umap.columns = ['UMAP_1', 'UMAP_2']
    df = adata.obs[['IR_VDJ_1_cdr3', 'IR_VDJ_1_v_gene', 'subtype']]
    df = df.reset_index()
    df = df.merge(umap, left_index=True, right_index=True)
    df = df.rename(columns={'IR_VDJ_1_v_gene':'v_gene', 'IR_VDJ_1_cdr3': 'cdr3s_aa', 'index': 'cell_id'})
    df = df.replace(to_replace="nan", value="None")
    return df

def get_degs(adata):
    subtypes = adata.uns['rank_genes_groups']["names"].dtype.names

    genes = pd.DataFrame(adata.uns['rank_genes_groups']['names'].tolist(), columns=adata.uns['rank_genes_groups']['names'].dtype.names)
    adjpvals = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'].tolist(), columns=adata.uns['rank_genes_groups']['pvals_adj'].dtype.names)
    logfc = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'].tolist(), columns=adata.uns['rank_genes_groups']['logfoldchanges'].dtype.names)

    data = pd.DataFrame()
    for subtype in subtypes:
        df = pd.DataFrame()
        df['gene'] = genes[subtype]
        df['adj_pval'] = adjpvals[subtype]
        df['log_fc'] = logfc[subtype]
        df['subtype'] = subtype

        df = df[(df['adj_pval'] < PTHRESH) & (df['log_fc'] > FCTHRESH)]
        df = df.reset_index(drop=True)
        df = df.sort_values('log_fc', ascending=False)
        df = df[:TOP_N]

        data = pd.concat([data, df], ignore_index=True)

    return data

def get_probabilities(adata):
    df = adata.obs[["IR_VDJ_1_cdr3", "pgen", "subtype"]]  
    df = df[(df["IR_VDJ_1_cdr3"] != "None") & (df["pgen"] != "None")]
    df["log10_probability"] = [math.log10(float(prob)) for prob in df["pgen"].tolist()]

    df = df.rename(columns={'IR_VDJ_1_cdr3': 'cdr3s_aa', 'pgen': 'probability'})
    df = df.reset_index(drop=True)

    return df

# metadata = pd.read_csv(os.path.join(data_dir, "metadata.tsv"), sep="\t")
# metadata = metadata.to_dict('records')

# probabilities = pd.read_csv(os.path.join(data_dir, "probabilities.tsv"), sep="\t")
# probabilities = probabilities.to_dict("records")

# degs = pd.read_csv(os.path.join(data_dir, "degs.tsv"), sep="\t")
# degs = degs.to_dict("records")

# data = {
#     "metadata": metadata,
#     "probabilities": probabilities,
#     "degs": degs
# }


if __name__ == "__main__":
    filename = sys.argv[1]
    data = get_data(filename)
    data = json.dumps(data, indent=4)

    app_dir = os.path.dirname(os.path.abspath(__file__))
    # app_dir = os.path.abspath(os.path.join(app_dir, "../..", "build"))
    index_template=os.path.join(app_dir, "build", "index.html")
    template = Template(open(index_template,"r").read())


    html = template.render(data=data)
    output_html=os.path.join(app_dir,"build", "pictcr.html")
    output = open(output_html,"w")

    js_txt = open(os.path.join(app_dir, "build", "main.js"), 'r').read()
    css_txt = open(os.path.join(app_dir, "build", "main.css"), 'r').read()

    html = html.replace('<script src="./main.js"></script>', f"<script>{js_txt}</script>")
    html = html.replace('<link href="./main.css" rel="stylesheet">', f"<style>{css_txt}</style>")

    output.write(html)
    output.close()
# datalake_build=os.path.join(sys.argv[1], "build")
# shutil.copytree(build_folder, datalake_build)
