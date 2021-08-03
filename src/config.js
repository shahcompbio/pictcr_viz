export const CONSTANTS = {
  xParam: "UMAP_1",
  yParam: "UMAP_2",
  clonotypeParam: "clone_id",
  cellIdParam: "cell_id",
  subtypeParam: "subtype",
  logProbParam: "log10_probability",
};

export const CLONOTYPE_COLORS = [
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
];

export const INFO = {
  UMAP: {
    title: "Cdr3 Amino Acid Sequence UMAP",
    text: "UMAP embedding generated from the subset of CD8/CD4 T cells identified with CellAssign. The top 10 TRB sequences by frequency are mapped to cells in the embedding. The radius slider allows the user to highlight areas of high density for a given sequence that can be referenced against the subtype UMAP.",
  },
  SUBTYPEUMAP: {
    title: "Subtype UMAP",
    text: "UMAP embedding generated from the subset of CD8/CD4 T-cells identified with CellAssign. Subtypes are assigned from differential gene expression derived from Leiden clustering. The list of current subtypes and known gene markers is provided below.",
  },
  HEATMAP: {
    title: "Subtype to Sequence Heatmap",
    text: "Frequency of the top 10 TRB sequences is shown for each of the identified subtypes. The heatmap intensity captures the relative proportion of the given sequence with respect to total number of TCR+ sequences within each subtype (Frequency of Sequence in Subtype / Total TCR+ Cells in Subtype).",
  },
  HISTOGRAM: {
    title: "Generation Probabilities",
    text: "Distribution of log10 generation probabilities generated from OLGA for all TRB sequences. Description of the OLGA method can be found here: https://academic.oup.com/bioinformatics/article/35/17/2974/5292315/. The log10 probability of each top 10 TRB sequence can be referenced with respect to the entire distribution of the sample. Sequences with a lower probability of generation will be found in the left tail of the distribution.",
  },
  BARPLOT: {
    title: "Clonotype Expansion",
    text: "Clonal expansion in each subtype by plotting the fraction of cells that belong to sequences expanded to N cells. The fraction of cells with sequences that belong to a single cell (labeled 1 in legend) describes the relative proportion of unexpanded clonotypes. The fraction of cells with sequences that belong to >10 cells describes the relative proportion of highly expanded clonotypes within the subtype.",
  },
  TABLE: {
    title: "Differentially Expressed Genes",
    text: "Differentially expressed genes for each subtype (1 vs. All) using the Wilcoxon test. The top 50 statistically significant (p < 0.001) with a minimum log fold change of 0.25 and expression found in >= 50% of the subtype cells are included for each subtype. It is possible to sort DEGs by both p-value and log fold change.",
  },
  RANKED: {
    title: "Ranked Clonotype Frequency",
    text: "Ranked order of clonotype frequency for each subtype",
  },
};
