#!/usr/bin/env python3

import numpy as np
import argparse as ap
from sklearn import linear_model



'''Usage: ./01_initial_eRP_calc.py --genes_exp ~/taylorLab/VISION/data/RNAseq/scriptseq3.v3.kw2.IDlocexp.bed --TSS ~/mm10_genome/TSS_vM4_wgn_woh_fc.bed --TSS_intersect ~/taylorLab/VISION/replicating_regression/TSS_ccRE_w100_within200.bed --cCRE ~/taylorLab/VISION/data/ccRE/VISIONmusHem_ccREs_filterkw2.bed'''

parser = ap.ArgumentParser(description='Replicating Regression Step 1: Initial eRP calculation')
parser.add_argument('--genes_exp', action='store', nargs=1, type=str, required = True, help="file scriptseq3.v3.kw2.IDlocexp.bed which has the format of geneID_tab_chr_tab_start_tab_end_tab_cellTypei=TPMi;cellTypej=TPMj;")
parser.add_argument('--TSS', action='store', nargs=1, type=str, required=True, help='TSSs')
parser.add_argument('--TSS_intersect', action='store', nargs=1, type=str, required = True, help='file with TSS for genes and ccREs within 200bp')
parser.add_argument('--cCRE', action='store', nargs=1, type=str, required = True, help='VISIONmusHem_ccREs_filterkw2.bed file with candidate cis regulatory elements ID, index_set_label (accessibility), epigenetic_state_level_by_cellType')
#parser.add_argument('--w200', action='store', nargs=1, type=str, required = True, help='bedtools window -w 200 result of scriptseq3.v3.kw2.IDlocexp.bed and VISIONmusHem_ccREs_filterkw.bed')
#parser.add_argument('--w1kb', action='store', nargs=1, type=str, required = True, help='bedtools window -w 1000 result of scriptseq3.v3.kw2.IDlocexp.bed and VISIONmusHem_ccREs_filterkw.bed')
#parser.add_argument('--w1Mb', action='store', nargs=1, type=str, required = True, help='bedtools window -w 1000000 result of scriptseq3.v3.kw2.IDlocexp.bed and VISIONmusHem_ccREs_filterkw.bed but 1kb windows removed')
args=parser.parse_args()
genes_Y_file = args.genes_exp[0]
TSS = args.TSS[0]
TSS_intersect = args.TSS_intersect[0]
cCRE_file = args.cCRE[0]
#within_200 = args.w200[0]
#within_1kb = args.w1kb[0]
#within_1Mb = args.w1Mb[0] #excludes those that are within_1kb


cellTypesOI = ["Lsk", "Cmp", "Gmp", "Mep", "Cfue", "Eryad", "Cfum", "Imk", "Mon", "Neu", "G1e", "Er4"]
num_cellTypesOI = len(cellTypesOI)
cellTypesOI_index = {}
num_states = len(range(1,27))
#states_index = dict(zip((str(x) for x in range(1,27)), range(num_states)))

'''prepare a list of gene IDs, assign an index to correspond to a gene ID'''
genes = set()
with open(genes_Y_file) as f:
    for line in f:
        fields = line.strip('\r\n').split('\t')
        genes.add(fields[3])

genes=list(genes)
num_genes = len(genes)
gene_index = dict(zip(genes, range(num_genes)))


'''prepare a matrix mapping a gene ID to a vector of expression by cell type'''

Y_exp_matrix = np.zeros((num_cellTypesOI, num_genes), dtype = np.float64)
with open(genes_Y_file) as f:
    for j, line in enumerate(f):
        fields = line.strip('\r\n').split('\t')
        gene_ID = fields[3]
        expression_fields = fields[4].split(";")[:-1]
        for i, minifield in enumerate(expression_fields):
            cellType,expression = minifield.split("=")
            if j == 0:
                cellTypesOI_index[cellType] = i
            Y_exp_matrix[cellTypesOI_index[cellType], gene_index[gene_ID]] = np.log2(float(expression)+0.001)
'''prepare a list of ccREs, assign an index to correspond to these'''
ccREs = set()
with open (cCRE_file) as f:
    for line in f:
        fields = line.strip('\r\n').split('\t')
        ccREs.add(fields[0])

ccREs = list(ccREs)
ccRE_index = dict(zip(ccREs, range(len(ccREs))))

'''QUESTIONS CONCERNING "the proportions of each IDEAS state within the 200bp regions that cover the TSSs of genes"
1) is this any ideas call or IDEAS state of ccREs within 200 bp? Or actually intersecting (i.e. "cover"ing?) the TSS
2) proportion by count or by base pairs?
3) possible to have a zero vector if nothing intersects/within 200 bp?
ASSUME of ccREs that are within a 200 bp of TSS (100 on each side), count by base pairs within that 200, and yes'''


def get_overlap(range_1, range_2):
    range_1_len = range_1[1] - range_1[0]
    overlap_region = [max(range_1[0], range_2[0]), min(range_1[1], range_2[1])]
    len_overlap_region = overlap_region[1] - overlap_region[0]
    overlap = len_overlap_region/range_1_len
    return(overlap)


TSS_intersect_dict = {}
with open(TSS_intersect) as f:
    for line in open(TSS_intersect):
        fields = line.strip('\r\n').split('\t')
        geneID = fields[3]
        if geneID not in TSS_intersect_dict:
            TSS_intersect_dict[geneID] = []
        TSS_intersect_dict[geneID].append((fields[8], int(fields[1]), int(fields[5]), int(fields[6]))) #states, TSS_pos, ccRE_start, ccRE_stop

initial_props = np.zeros((num_cellTypesOI, num_genes, num_states),dtype=np.float64)
initial_regression_coeffs = np.ndarray((num_cellTypesOI, num_states))
initial_regression_intercept_resids = np.ndarray((num_cellTypesOI,2))
for key in gene_index:
    '''get the proportions of IDEAS states'''
    if key in TSS_intersect_dict:
        num_intersect = len(TSS_intersect_dict[key])
        for i in range(num_intersect):
            states, TSS_start, ccRE_start, ccRE_stop = TSS_intersect_dict[key][i]
            prop_overlap = get_overlap([TSS_start-100, TSS_start+100], [ccRE_start, ccRE_stop])
            states_minifields = states.split(";")[:-1]
            for k, element in enumerate(states_minifields):
                cellType, state = element.split("=")
                if state != "0" and cellType in cellTypesOI:
                    initial_props[cellTypesOI_index[cellType], gene_index[key], k] += prop_overlap
                    # if isinstance(initial_props[cellTypesOI_index[cellType], gene_index[key]], int):
                    #     state_array = np.zeros(num_states)
                    # else:
                    #     state_array = initial_props[cellTypesOI_index[cellType], gene_index[key]]
                    # state_array[states_index[state]] += prop_overlap
                    # initial_props[cellTypesOI_index[cellType], gene_index[key]] = state_array
    # else:
    #     for j in range(num_cellTypesOI):
    #         initial_props[j,gene_index[key]]= np.zeros(num_states)

'''linear regression'''

for cellType in cellTypesOI:
    indice = cellTypesOI_index[cellType]
    y = np.array(Y_exp_matrix[indice])
    #print(y.shape)
    #print(y.ndim)
    #print(y.size)

    X = np.array(initial_props[indice])
    #print(X.shape)
    #print(X.ndim)
    #print(X.size)

    fit_init_model = linear_model.LinearRegression().fit(X,y)

    '''metrics other than the residuals'''
    init_intercept = fit_init_model.intercept_
    init_coeffs = fit_init_model.coef_
    initial_regression_intercept_resids[indice,0] = init_intercept
    initial_regression_coeffs[indice] = init_coeffs

    '''get residuals using predict to figure out how far off it is'''
    predicted_y = fit_init_model.predict(X)
    residual = sum(y - predicted_y)
    initial_regression_intercept_resids[indice,1] = residual
