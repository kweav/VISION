#!/usr/bin/env python3

'''Use gene expression in 12 cell types to ...
1) score the ccREs for their regulatory potentials based on epigenetic states
2) map the ccREs to candidate genes
3) further select the most likely subset of ccREs for predicting gene expression

"investigated the effectiveness of the ccREs in explaining levels of gene expression....
Developed a modeling approach to evaluate how well the ccREs could account for levels of expression in the 12 cell types
for which the RNAseq measurements were determined in the same manner.... additional benefit of making predictions of target genes for each ccRE"'''

import numpy as np
import argparse as ap
from sklearn import linear_model


'''01_initial_eRP_calc'''
'''Usage: ./01_initial_eRP_calc.py --genes_exp ~/taylorLab/VISION/data/RNAseq/scriptseq3.v3.kw2.IDlocexp.bed --TSS ~/mm10_genome/TSS_vM4_wgn_woh_fc.bed --TSS_intersect ~/taylorLab/VISION/replicating_regression/TSS_ccRE_w100_within200.bed --cCRE ~/taylorLab/VISION/data/ccRE/VISIONmusHem_ccREs_filterkw2.bed'''

parser = ap.ArgumentParser(description='Replicating Regression Step 1: Initial eRP calculation')
parser.add_argument('--genes_exp', action='store', nargs=1, type=str, required = True, help="file scriptseq3.v3.kw2.IDlocexp.bed which has the format of geneID_tab_chr_tab_start_tab_end_tab_cellTypei=TPMi;cellTypej=TPMj;")
parser.add_argument('--TSS', action='store', nargs=1, type=str, required=True, help='file ~/mm10_genome/TSS_vM4_wgn_woh_fc.bed')
parser.add_argument('--TSS_intersect', action='store', nargs=1, type=str, required = True, help='file TSS_ccRE_ENSG_w100_within200.bed with TSS for genes and ccREs within 200bp')
parser.add_argument('--cCRE', action='store', nargs=1, type=str, required = True, help='VISIONmusHem_ccREs_wct_prop_filterkw2.bed file with candidate cis regulatory elements ID, index_set_label (accessibility), epigenetic_state_level_by_cellType')
parser.add_argument('--w1kb', action='store', nargs=1, type=str, required = True, help='file genes_ccRE_ID_1kb_window.bed bedtools window -w 1000 result of scriptseq3.v3.kw2.IDlocexp.bed and VISIONmusHem_ccREs_filterkw.bed')
parser.add_argument('--w1Mb', action='store', nargs=1, type=str, required = True, help='file genes_ccRE_ID_1Mb_diff_1kb.bed bedtools window -w 1000000 result of scriptseq3.v3.kw2.IDlocexp.bed and VISIONmusHem_ccREs_filterkw.bed but 1kb windows removed')
args=parser.parse_args()
genes_Y_file = args.genes_exp[0]
TSS = args.TSS[0]
TSS_intersect = args.TSS_intersect[0]
cCRE_file = args.cCRE[0]
within_1kb = args.w1kb[0]
within_1Mb = args.w1Mb[0] #excludes those that are within_1kb


cellTypesOI = ["Lsk", "Cmp", "Gmp", "Mep", "Cfue", "Eryad", "Cfum", "Imk", "Mon", "Neu", "G1e", "Er4"]
num_cellTypesOI = len(cellTypesOI)
cellTypesOI_index = {}
num_states = len(range(1,27))
#states_index = dict(zip((str(x) for x in range(1,27)), range(num_states)))

def encode_genomeloc_to_ID(chr, start, end):
    chrEncode = {'chr1':'A',
                'chr2':'B',
                'chr3':'C',
                'chr4':'D',
                'chr5':'E',
                'chr6':'F',
                'chr7':'G',
                'chr8':'H',
                'chr9':'I',
                'chr10':'J',
                'chr11':'K',
                'chr12':'L',
                'chr13':'M',
                'chr14':'N',
                'chr15':'O',
                'chr16':'P',
                'chr17':'Q',
                'chr18':'R',
                'chr19':'S',
                'chrX':'T',
                'chrY':'U'}
    ID = chrEncode[chr]
    h_start = hex(start)[2:]
    if len(h_start) < 7:
        for i in range(7 - len(h_start)):
            h_start = "."+h_start
    h_start = "h."+h_start
    ID+= h_start
    h_end = hex(end)[2:]
    if len(h_end) < 7:
        for i in range(7 - len(h_end)):
            h_end = "."+h_end
    h_end = "h."+h_end
    ID += h_end

    return (ID)

'''prepare a list of gene IDs, assign an index to correspond to a gene ID'''
genes = set()
with open(genes_Y_file) as f:
    for line in f:
        fields = line.strip('\r\n').split('\t')
        genes.add(fields[3])

genes=list(genes)
num_genes = len(genes)
gene_index = dict(zip(genes, range(num_genes)))
del genes


'''prepare a matrix mapping a gene ID to a vector of expression by cell type'''

Y_exp_matrix = np.zeros((num_cellTypesOI, num_genes), dtype = np.float64)
with open(genes_Y_file) as f:
    for j, line in enumerate(f):
        fields = line.strip('\r\n').split('\t')
        gene_ID = fields[3].strip()
        expression_fields = fields[4].split(";")[:-1]
        for i, minifield in enumerate(expression_fields):
            cellType,expression = minifield.split("=")
            if j == 0:
                cellTypesOI_index[cellType] = i
            Y_exp_matrix[cellTypesOI_index[cellType], gene_index[gene_ID]] = np.log2(float(expression)+0.001)

'''prepare a list of ccREs, assign an index to correspond to these
Also make an array that stores the proportion of the ccRE that is in each state in each cell type'''

num_ccRE = 205019
ccREs_state_prop = np.zeros((num_ccRE,num_cellTypesOI,num_states+1))
ccRE_index = {}
with open (cCRE_file) as f:
    for i, line in enumerate(f):
        fields = line.strip('\r\n').split('\t')
        #ccREs.add(fields[0])
        ccRE_index[fields[0]] = i
        subfields = fields[2].split(";")[:-1] #Lsk=(0_0.821)(9_0.179);Cmp=(7_1.0);Mep=(....)
        for field in subfields: #Lsk=(0_0.821)(9_0.179)
            cellType, statesRepresented = field.split("=") #Lsk (0_0.821)(9_0.179)
            if cellType in cellTypesOI:
                for stateProp in statesRepresented.split(")")[:-1]: #(0_0.821   (9_0.179
                    state, prop = stateProp.replace("(","").split("_") #0  0.821
                    state = int(state)
                    prop = float(prop)
                    ccREs_state_prop[i, cellTypesOI_index[cellType], state] = prop

'''QUESTIONS CONCERNING "the proportions of each IDEAS state within the 200bp regions that cover the TSSs of genes"
1) is this any ideas call or IDEAS state of ccREs within 200 bp? Or actually intersecting (i.e. "cover"ing?) the TSS
2) proportion by count or by base pairs?
3) possible to have a zero vector if nothing intersects/within 200 bp?
ASSUME of ccREs that are within a 200 bp of TSS (100 on each side), count by base pairs within that 200 (divide by min(200, len_of_ccRE)), and yes
NOTICE file --TSS_intersect ~/taylorLab/VISION/replicating_regression/TSS_ccRE_ENSG_w100_within200.bed includes whether a ccRE overlaps the same gene multiple times just at its different TSS's, and here I'm reflecting that by adding the prop overlap every time'''

def get_overlap(range_1, range_2): #range_1 is the one you want to divide by its length
    range_1_len = range_1[1] - range_1[0]
    range_2_len = range_2[1] - range_2[0]
    overlap_region = [max(range_1[0], range_2[0]), min(range_1[1], range_2[1])]
    len_overlap_region = overlap_region[1] - overlap_region[0]
    overlap = len_overlap_region/min(range_1_len, range_2_len)
    return (overlap)

TSS_intersect_dict = {}
with open(TSS_intersect) as f:
    for line in open(TSS_intersect):
        fields = line.strip('\r\n').split('\t')
        geneID = fields[3].strip()
        if geneID not in TSS_intersect_dict:
            TSS_intersect_dict[geneID] = []
        TSS_intersect_dict[geneID].append((fields[8], int(fields[1]), fields[4], int(fields[5]), int(fields[6]))) #states, TSS_pos, chr, ccRE_start, ccRE_stop

initial_props = np.zeros((num_cellTypesOI, num_genes, num_states),dtype=np.float64)
for key in gene_index:
    '''get the proportions of IDEAS states'''
    if key in TSS_intersect_dict:
        num_intersect = len(TSS_intersect_dict[key])
        for i in range(num_intersect):
            states, TSS_start, chr, ccRE_start, ccRE_stop = TSS_intersect_dict[key][i]
            ccRE_ID = encode_genomeloc_to_ID(chr, ccRE_start, ccRE_stop)
            states_minifields = states.split(";")[:-1]
            for k, element in enumerate(states_minifields):
                cellType, state = element.split("=")
                if state != "0" and cellType in cellTypesOI:
                    prop_overlap = get_overlap([TSS_start-100, TSS_start+100], [ccRE_start, ccRE_stop])
                    initial_props[cellTypesOI_index[cellType], gene_index[key], int(state)-1] += prop_overlap #state-1 because state 0 not in this array, therefore state 1 in index 0, etc....

initial_regression_coeffs = np.ndarray((num_cellTypesOI, num_states+1))
initial_regression_coeffs[:,0] = 0 #state 0 has coefficient 0

'''linear regression'''
for cellType in cellTypesOI:
    print(cellType)
    indice = cellTypesOI_index[cellType]
    #y = np.array(Y_exp_matrix[indice])
    #print(y)
    #print(y.shape)
    #print(y.ndim)
    #print(y.size)

    #X = np.array(initial_props[indice])
    #print(X)
    #print(X.shape)
    #print(X.ndim)
    #print(X.size)

    fit_init_model = linear_model.LinearRegression().fit(np.array(initial_props[indice]),np.array(Y_exp_matrix[indice]))

    '''metrics other than the residuals'''
    init_coeffs = fit_init_model.coef_ #these seem to be 26 zeros for each cell type
    print(init_coeffs)
    initial_regression_coeffs[indice,1:] = init_coeffs #set Beta coefficients for states 1-26; state 0 previously set to 0

    '''get residuals using predict to figure out how far off it is'''
    #predicted_y = fit_init_model.predict(X)
    #residual = sum(y - predicted_y)
    #print(residual)

'''want to find the initial eRP_0 score using initial_regression_coeffs and ccREs_state_prop to find weighted sum of beta coefficients == DOT PRODUCT
where weight of each state is the proportion of that state which occurred within the ccRE
initial_regression_coeffs has shape (num_cellTypesOI, num_states+1), so initial_regression_coeffs[celltypej,:] = [B0, B1, B2, ..., B26]
ccREs_state_prop has shape (num_ccRE, num_cellTypesOI, num_states+1), so ccREs_state_prop[ccREi, celltypej,:] = [P0, P1, P2, ..., P26]
result should have shape (num_ccRE, num_cellTypesOI), so result[ccREi, cellTypej] = P0*B0 + P1*B1 + P2*B2 + ... + P26*B26
'''

initial_eRP = np.sum(ccREs_state_prop*initial_regression_coeffs.reshape((1,num_cellTypesOI,num_states+1)),axis=2) #shape should be (num_ccRE, num_cellTypesOI)

'''QUESTION2: Still unsure what is meant by: "For every genome location that has at least one ccRE in some cell types, we further assigned the eRP​0 scores to the DNA segments at the same location in all other cell types.
This yields a 12-dimensional vector of eRP​0 s​cores for each location with at least one ccRE.
12-dimensional vector == vector with 12 elements?
Assuming that my method handles this....
"'''


'''02_Preselect_ccRE_gene_pairs'''
''' Find the correlations of each ccRE with gene expression of all genes within 1Mb '''
within_distance_array = np.zeros((num_ccRE, num_genes,3)) #3rd dimension is within1kb, within1Mb, within1Mbdiff1kb.
for line in open(within_1kb):
    fields = line.strip('\r\n').split('\t')
    geneID = fields[3]
    gIndex = gene_index[geneID]
    ccRE_ID = fields[5]
    cIndex = ccRE_index[ccRE_ID]
    within_distance_array[cIndex, gIndex, [0,1]] = 1

for line in open(within_1Mb): #REMINDER: this file has ccREs that are within 1Mb but NOT 1kb.
    fields = line.strip('\r\n').split('\t')
    geneID = fields[3]
    gIndex = gene_index[geneID]
    ccRE_ID = fields[5]
    cIndex = ccRE_index[ccRE_ID]
    within_distance_array[cIndex,gIndex,[2]]=1

correlations = np.zeros((num_ccRE, num_genes))
#tile the ccRE erp vector to match how many genes or the second dimension of the expression_vectors array??
expression_vectors = Y_exp_matrix[:,np.nonzero(within_distance_array[57,:,1])[0]]
print(expression_vectors)
print(expression_vectors.shape)
