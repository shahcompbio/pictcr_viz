export const CONSTANTS = {
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  clonotypeParam: "IR_VDJ_1_junction_aa",
  cellIdParam: "cell_id",
  subtypeParam: "cell_type",
  logProbParam: "log10_probability",
};
export const PHENOTYPE_COLORS = [
  "#5E4FA2",
  "#3288BD",
  "#66C2A5",
  "#FEE08B",
  "#FDAE61",
  "#F46D43",
  "#D53E4F",
  "#c9cc76",
  "#9E0142",
  "#C6AEFF",
  "#BDD8FF",
  "#BDFFB2",
  "#FFC8AE",
  "#FF9FBB",
  "#b2dbd6",
  "#ffd470",
];
export const CLONOTYPE_COLORS = [
  "[0.369,0.31,0.635,1.0]",
  "[0.196,0.533,0.741,1.0]",
  "[0.4,0.761,0.647,1.0]",
  "[0.996,0.878,0.545,1.0]",
  "[0.957,0.427,0.263,1.0]",
  "[0.835,0.243,0.31,1.0]",
  "[0.788,0.8,0.463,1.0]",
  "[0.62,0.004,0.259,1.0]",
  "[0.776,0.682,1.0,1.0]",
  "[0.741,0.847,1.0,1.0]",
  "[0.741,1.0,0.698,1.0]",
  "[1.0,0.784,0.682,1.0]",
  "[1.0,0.624,0.733,1.0]",
  "[0.698,0.859,0.839,1.0]",
  "[1.0,0.831,0.439,1.0]",
];
/*export const CLONOTYPE_COLORS = [
  "#674172",
  "#098dde",
  "#fa832f",
  "#0e5702",
  "#c20c1e",
  "#911eb4",
  "#fc97bc",
  "#469990",
  "#b5762a",
  "#5aebed",
  "#8f8f3f",
  "#ed1a1a",
];*/
export const INFO = {
  SUBTYPEDOUGH: {
    title: "Phenotype Distribution",
    text: "",
  },
  CLONOTYPEDOUGH: {
    title: "Clone Distribution",
    text: "",
  },
  UMAP: {
    title: "Clone UMAP",
    text: "UMAP embedding generated from the subset of CD8/CD4 T cells identified with CellAssign. The top 10 clones by frequency are mapped to cells in the embedding. The radius slider allows the user to highlight areas of high density for a given sequence that can be referenced against the phenotype UMAP.",
  },
  SUBTYPEUMAP: {
    title: "Phenotype UMAP",
    text: "UMAP embedding generated from the subset of CD8/CD4 T-cells identified with CellAssign. Phenotypes are assigned from differential gene expression derived from Leiden clustering. The list of current phenotypes and known gene markers is provided below.",
  },
  HEATMAP: {
    title: "Phenotype to Clone Heatmap",
    text: "Frequency of the top clones is shown for each of the identified phenotype. The heatmap intensity captures the relative proportion of the given sequence with respect to total number of TCR+ sequences within each subtype (Frequency of Sequence in Subtype / Total TCR+ Cells in Subtype).",
  },
  HISTOGRAM: {
    title: "Generation Probabilities",
    text: "Distribution of log10 generation probabilities generated from OLGA for all TRB sequences. Description of the OLGA method can be found here: https://academic.oup.com/bioinformatics/article/35/17/2974/5292315/. The log10 probability of each top 10 TRB sequence can be referenced with respect to the entire distribution of the sample. Sequences with a lower probability of generation will be found in the left tail of the distribution.",
  },
  BARPLOT: {
    title: "Clone Expansion",
    text: "Clonal expansion in each subtype by plotting the fraction of cells that belong to sequences expanded to N cells. The fraction of cells with sequences that belong to a single cell (labeled 1 in legend) describes the relative proportion of unexpanded clonotypes. The fraction of cells with sequences that belong to >10 cells describes the relative proportion of highly expanded clonotypes within the subtype.",
  },
  TABLE: {
    title: "Differentially Expressed Genes",
    text: "Differentially expressed genes for each subtype (1 vs. All) using the Wilcoxon test. The top 50 statistically significant (p < 0.001) with a minimum log fold change of 0.25 and expression found in >= 50% of the subtype cells are included for each subtype. It is possible to sort DEGs by both p-value and log fold change.",
  },
  RANKED: {
    title: "Ranked Clone Frequency",
    text: "Ranked order of clonotype frequency for each subtype",
  },
  SANKEY: {
    title: "Sankey",
    text: "",
  },
};
