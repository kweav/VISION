#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model
import os
from scipy import stats


def main():
    parser = generate_parser()
    args = parser.parse_args()
    setup_file_locs(args, args.where_run, args.other_path)
    setup_threads(args.threads)
    model = regress_gene_cre(args.cre_state_file, args.tss_state_file, args.exp_file)

    if args.chroms == 'all':
        chrom_list = []
        for i in range(1,20):
            chrom_list.append('chr{}'.format(i))
        chrom_list.append('chrX')
    else:
        chrom_list = args.chroms

    for chrom in chrom_list:
        model.subset_based_on_chrom(chrom)
        model.subset_based_on_group(args.exp_type, args.m_thresh, args.s_thresh)
        model.set_initial_betas()
        model.find_weighted_sum()
        model.find_initial_pairings(args.dist_gamma, args.lessone, args.cre_dist)
    #     model.refine_pairs(args.iterations)

def generate_parser():
    parser = ap.ArgumentParser(description = 'VISION regression of state and distance to assign CREs to genes based on ability to predict gene expression')
    parser.add_argument('--where', action='store', dest='where_run', type=str, default='mine', help='{comp, mine, marcc, other}; adds path to files; if "other" used provide the other path in the --otherpath argument')
    parser.add_argument('--threads', action='store', dest='threads', type=str, default="1")
    parser.add_argument('-c', '--cre_state_file', action='store', dest='cre_state_file', type=str, default='ccre_state_prop.npz', help='default stores the file name, path will be added based on --where argument')
    parser.add_argument('-t', '--tss_state_file', action='store', dest='tss_state_file', type=str, default='TSS_window_state_prop.npz', help='default stores the file name, path will be added based on ---where argument')
    parser.add_argument('-e', '--exp_file', action='store', dest='exp_file', type=str, default='TSS_expression.npz', help='default stores the file name, path will be added based on --where argument')
    parser.add_argument('-d', '--dist', action='store', type=int, dest='cre_dist', default='1500000')
    parser.add_argument('-g', '--group', action='store', dest='exp_type', type=int, default=0, help='{0, 1, 2, 3, 4} Group subselection. 0 is for all and is default')
    parser.add_argument('--chroms', action='store', nargs = '+', dest='chroms', type=str, default='all')
    parser.add_argument('-m', '--mean_thresh', action='store', dest='m_thresh', type=int, default=-4)
    parser.add_argument('-s', '--stdev_thresh', action='store', dest='s_thresh', type=int, default=2)
    parser.add_argument('--gamma', action='store', dest='dist_gamma', type=float, default=0.7)
    parser.add_argument('-i', '--iterations', action='store', dest='iterations', type=int, default=100)
    parser.add_argument('-l', '--lessone', action='store', dest='lessone', type=int, default=0, help='Cell type to leave out with 0 indexing. 0 is default')
    parser.add_argument('--otherpath', action='store', dest='other_path', default='NA', help='use to give other path if not one of the 3 previously specified. MUST be the same for ALL 3 input files')
    parser.add_argument('--correlation', action='store', dest='correlation', type=float, default=0.1)
    return parser

def setup_threads(threads):
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["OPENBLAS_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = threads
    os.environ["NUMEXPR_NUM_THREADS"] = threads

def setup_file_locs(args, where_run, other_path):
    argumentToAdd = {'mine': '/Users/kateweaver/taylorLab/VISION/regression_with_dist/ccRE_centric/inputs/',
                     'marcc': '/home-3/kweave23@jhu.edu/data/kweave23/regression_with_dist/inputs/',
                     'comp': '/home/kweave23/regression_with_dist/ccRE_centric/inputs',
                     'other': other_path}

    args.exp_file = argumentToAdd[where_run] + args.exp_file
    args.cre_state_file = argumentToAdd[where_run] + args.cre_state_file
    args.tss_state_file = argumentToAdd[where_run] + args.tss_state_file

class regress_gene_cre():
    def __init__(self, cre_state_file, tss_state_file, exp_file):
        '''
        self.exp_values_all:
        self.cellIndex:
        self.cell_to_index:
        self.TSS_chr_all:
        self.TSSs_all:
        self.cellN:
        self.TSS_window_props_all:
        self.TSS_window_chr_all:
        self.TSS_window_coord_all:
        self.stateN:
        '''
        self.exp_values_all, self.cellIndex, self.cell_to_index, self.TSS_chr_all, self.TSSs_all = self.load_expression(exp_file)
        self.cellN = self.cellIndex.shape[0]
        self.TSS_window_props_all, self.TSS_window_chr_all, self.TSS_window_coord_all = self.load_TSS_window_states(tss_state_file)
        self.cre_props_all, self.cre_chr_all, self.cre_coords_all = self.load_cre_states(cre_state_file)
        self.stateN  = self.cre_props_all.shape[2] #Note this value includes state 0 in the count although we're going to ignore the contribution of state 0

        self.chrSizes = {'chr1':195471971,
                         'chr2':182113224,
                         'chrX':171031299,
                         'chr3':160039680,
                         'chr4':156508116,
                         'chr5':151834684,
                         'chr6':149736546,
                         'chr7':145441459,
                         'chr10':130694993,
                         'chr8':129401213,
                         'chr14':124902244,
                         'chr9':124595110,
                         'chr11':122082543,
                         'chr13':120421639,
                         'chr12':120129022,
                         'chr15':104043685,
                         'chr16':98207768,
                         'chr17':94987271,
                         'chrY':91744698,
                         'chr18':90702639,
                         'chr19':61431566}

    def find_valid_cellTypes(self, cellIndexOI):
        valid = np.array([np.where(cellIndexOI == x) for x in self.cellIndex if np.isin(x, cellIndexOI)]).reshape(-1)
        return valid

    def load_expression(self, exp_file):
        npzfile = np.load(exp_file, allow_pickle = True)
        exp_values = npzfile['exp']
        exp_cellIndex = npzfile['cellIndex'] #index to celltype
        exp_cell_to_index = {} #celltype to index
        for i in range(exp_cellIndex.shape[0]):
            exp_cell_to_index[exp_cellIndex[i]] = i
        TSS_info = npzfile['TSS']
        TSS_chr = TSS_info[:,0]
        TSSs = np.around(TSS_info[:,1].astype(np.float32)).astype(np.int32)
        return exp_values, exp_cellIndex, exp_cell_to_index, TSS_chr, TSSs

    def load_TSS_window_states(self, tss_state_file):
        npzfile = np.load(tss_state_file, allow_pickle = True)
        TSS_window_props = npzfile['props'].astype(np.float32)
        TSS_cellIndex = npzfile['cellIndex']
        valid = self.find_valid_cellTypes(TSS_cellIndex)
        TSS_window_props_valid = TSS_window_props[:,valid,:]
        TSS_window_index = npzfile['ccREIndex']
        TSS_window_chr = TSS_window_index[:,0]
        TSS_window_coord = TSS_window_index[:,1:].astype(np.int32)
        return TSS_window_props_valid, TSS_window_chr, TSS_window_coord

    def load_cre_states(self, cre_state_file):
        npzfile = np.load(cre_state_file, allow_pickle = True)
        cre_props = npzfile['props'].astype(np.float32)
        cre_cellIndex = npzfile['cellIndex']
        valid = self.find_valid_cellTypes(cre_cellIndex)
        cre_props_valid = cre_props[:,valid,:]
        cre_index = npzfile['ccREIndex']
        cre_chr = cre_index[:,0]
        cre_coords = cre_index[:,1:].astype(np.int32)
        return cre_props_valid, cre_chr, cre_coords

    def subset_based_on_chrom(self, chrom):
        self.chrom = chrom
        where_chr = np.where(self.TSS_chr_all == self.chrom)[0]
        self.exp_values = self.exp_values_all[where_chr]
        self.TSSs = self.TSSs_all[where_chr]

        self.TSS_window_props = self.TSS_window_props_all[where_chr]
        self.TSS_window_coord = self.TSS_window_coord_all[where_chr]

        where_chr = np.where(self.cre_chr_all == self.chrom)[0]
        self.cre_props = self.cre_props_all[where_chr]
        self.cre_coords = self.cre_coords_all[where_chr]

    def subset_based_on_group(self, group, m_thresh, s_thresh):
        self.m_thresh = m_thresh
        self.s_thresh = s_thresh
        self.group = group
        if self.group != 0:
            m = np.mean(self.exp_values, axis=1)
            s = np.std(self.exp_values, axis=1, ddof=1)
            group_1 = np.where((m <= self.m_thresh) & ( s <= self.s_thresh))[0]
            group_2 = np.where((m <= self.m_thresh) & (s > self.s_thresh))[0]
            group_3 = np.where((m > self.m_thresh) & (s > self.s_thresh))[0]
            group_4 = np.where((m > self.m_thresh) & (s <= self.s_thresh))[0]
            group_dict = {1: group_1,
                          2: group_2,
                          3: group_3,
                          4: group_4}
            self.exp_values = self.exp_values[group_dict[self.group]]
            self.TSSs = self.TSSs[group_dict[self.group]]
            self.TSS_window_props = self.TSS_window_props[group_dict[self.group]]
            self.TSS_window_coord = self.TSS_window_coord[group_dict[self.group]]
        self.tssN = self.TSSs.shape[0]

    def linear_regression(self, X, Y, intercept=False):
        '''reshape arrays to match function input requirements
        remove contribution of state 0 from X; will set this state's contribution to zero afterwards'''
        X = X.reshape(-1, self.stateN)[:,1:]
        Y = Y.reshape(-1,)
        fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X,Y)
        model_coeffs = fit_model.coef_
        r_squared = fit_model.score(X,Y)
        return {'coeffs': model_coeffs, 'rsquare': r_squared}

    def set_initial_betas(self):
        '''do initial regression of proportion of states within two-sided 75kbp window against expression to find initial coefficients for all states except 0
        set contribution of state 0 to be 0'''
        initial_coeffs = self.linear_regression(self.TSS_window_props, self.exp_values)['coeffs']
        self.initial_coeffs = np.hstack(([0], initial_coeffs)).reshape((1,-1))

    def find_weighted_sum(self):
        '''find weighted sum of ccRE props with initial_coeffs
        given a 3 dimensional array ccREn * cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                1 dimensional array stateN with c_i equal to the initial coefficient for state_i
        return a 2 dimensional array ccREn * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''

        self.cre_weighted_sum = np.sum(self.cre_props * self.initial_coeffs, axis=2)

    def compute_adj_distance(self, starts, stops, TSS):
        locsTA = (starts + stops) // 2
        distance = np.abs(locsTA- TSS) #using abs because don't care if upstream or downstream
        adjusted_distance = distance**self.dist_gamma
        adjusted_distance[np.where(adjusted_distance == 0)[0]] = 1 #anywhere that distance is 0, set to identity.
        return adjusted_distance

    def adjust_by_distance(self, to_adjust, TSS, starts, stops):
        adj_dist = self.compute_adj_distance(starts, stops, TSS)
        Y = np.ones(to_adjust.shape[0]) #adjustment array
        Y /= adj_dist
        adjusted = np.multiply(to_adjust, Y.reshape(-1,1))
        return adjusted

    def find_initial_pairings(self, dist_gamma, lessone, cre_dist):
        '''for each TSS:
            1) find CREs within distance of interest
            2) adjust the CREs weighted sum of proporitions by their distance from the TSS
            3) for each CRE within distance of interest
                3.1) find correlation of expression with this adjsuted weighted sum
                3.2) decide whether to initally pair or not '''
        self.dist_gamma = dist_gamma
        self.lessone = lessone
        writeCorrelations = open('correlations_wrn_{}.txt'.format(self.chrom), 'a')
        writeMetrics = open('correlation_wrn_metrix_{}.txt'.format(self.chrom), 'a')
        for i in range(self.tssN): #can I collapse this from a for loop to just fancy numpy?
            TSS = self.TSSs[i]
            #find CREs within distance of interest using containment
            windowMin = max(0, TSS - cre_dist)
            windowMax = min(TSS + cre_dist, self.chrSizes[self.chrom])
            CREs_within = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
            #writeNum.write('{}\t{}\t'.format(i, np.sum(CREs_within)))
            CREs_within_weighted_sum = self.cre_weighted_sum[CREs_within]
            CREs_within_starts = self.cre_coords[CREs_within, 0]
            CREs_within_stops = self.cre_coords[CREs_within, 1]
            #adjust their weighted_sums by their distance
            CREs_within_adjusted = self.adjust_by_distance(CREs_within_weighted_sum, TSS, CREs_within_starts, CREs_within_stops)
            where_row_not_zero = np.sum(CREs_within_adjusted, axis=1) != 0
            #for each CRE within distance of interest and as long as its state isn't 0 across all 12 cell types find correlation of expression with this adjusted weighted sum
            #corr_matrix, pvalues = stats.spearmanr(self.exp_values[i].reshape((1,-1)), CREs_within_adjusted, axis=1)
            corr_matrix, pvalues = stats.spearmanr(self.exp_values[i].reshape((1,-1)), CREs_within_adjusted[where_row_not_zero], axis=1)
            if isinstance(corr_matrix, np.ndarray):
                exp_cre_corr = corr_matrix[0,1:]
                np.savetxt(writeCorrelations, exp_cre_corr)
            else:
                exp_cre_corr = corr_matrix
                writeCorrelations.write('{}\n'.format(exp_cre_corr))
            avg = np.mean(exp_cre_corr)
            stdev = np.std(exp_cre_corr, ddof=1)
            median = np.median(exp_cre_corr)
            q25 = np.quantile(exp_cre_corr, 0.25)
            q5 = np.quantile(exp_cre_corr, 0.5)
            q75 = np.quantile(exp_cre_corr, 0.75)
            min_val = np.amin(exp_cre_corr)
            max_val = np.amax(exp_cre_corr)
            num_nan = np.sum(np.isnan(exp_cre_corr))
            num_inf = np.sum(np.isinf(exp_cre_corr))
            writeMetrics.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(avg, stdev, median, q25, q5, q75, min_val, max_val, num_nan, num_inf))
        writeCorrelations.close()
        writeMetrics.close()





    def refine_pairs(self, iterations):
        return 0


main()
