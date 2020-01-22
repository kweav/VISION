#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model, metrics
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
        model.build_pairing_array(args.iterations)
        model.build_metric_arrays()
        model.set_initial_betas()
        model.find_initial_weighted_sum()
        model.drive_pairing(args.dist_gamma, args.lessone, args.cre_dist, args.correlation)

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
    parser.add_argument('--correlation', action='store', dest='correlation', type=float, default=0.4)
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
                     'comp': '/home/kweave23/regression_with_dist/ccRE_centric/inputs/',
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
        TSSs = np.around(TSS_info[:,1].astype(np.float64)).astype(np.int32)
        return exp_values, exp_cellIndex, exp_cell_to_index, TSS_chr, TSSs

    def load_TSS_window_states(self, tss_state_file):
        npzfile = np.load(tss_state_file, allow_pickle = True)
        TSS_window_props = npzfile['props'].astype(np.float64)
        TSS_cellIndex = npzfile['cellIndex']
        valid = self.find_valid_cellTypes(TSS_cellIndex)
        TSS_window_props_valid = TSS_window_props[:,valid,:]
        TSS_window_index = npzfile['ccREIndex']
        TSS_window_chr = TSS_window_index[:,0]
        TSS_window_coord = TSS_window_index[:,1:].astype(np.int32)
        return TSS_window_props_valid, TSS_window_chr, TSS_window_coord

    def load_cre_states(self, cre_state_file):
        npzfile = np.load(cre_state_file, allow_pickle = True)
        cre_props = npzfile['props'].astype(np.float64)
        cre_cellIndex = npzfile['cellIndex']
        valid = self.find_valid_cellTypes(cre_cellIndex)
        cre_props_valid = cre_props[:,valid,:]
        #filter out any CREs which are fully state 0 in all 12 cell types
        mask = np.array(np.hstack(([0], np.tile([1], cre_props_valid.shape[2]-1))))
        cre_props_adjusted = np.sum(cre_props_valid * mask, axis=2)
        where_row_not_zero = np.sum(cre_props_adjusted, axis=1) != 0
        cre_props_valid = cre_props_valid[where_row_not_zero]
        cre_index = npzfile['ccREIndex']
        cre_chr = cre_index[where_row_not_zero,0]
        cre_coords = cre_index[where_row_not_zero,1:].astype(np.int32)
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
        self.creM = self.cre_props.shape[0]
        self.creIndex_range = np.arange(self.creM)

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
        return {'coeffs': model_coeffs, 'rsquare': r_squared, 'fit_model': fit_model}

    def build_pairing_array(self, iterations):
        self.iter = iterations
        self.pairing_array = np.full((self.tssN, self.creM, self.iter+1), -1, dtype=np.int32) #3D array of tssN * creM * iter+1 to store pairings for initial pairing and refinement iterations

    def build_metric_arrays(self):
        self.mse_array = np.empty((self.tssN, self.iter+1), dtype=np.float32) #2D array of tssN * iter+1 to store MSE for initial pairings and refined pairings

    def find_yhat(self, linear_model, X):
        X = X.reshape(-1, self.stateN)[:,1:]
        return linear_model.predict(X)

    def find_MSE(self, i, iter_val): #iter_val should be -1 if calling for initial pairing.
        cre_inclusions = self.creIndex_range[self.pairing_array[i,:,iter_val+1] == 1]
        X = self.cre_props[cre_inclusions, :, :][:,self.lessone_range,:]
        #adjust by distance
        cre_starts = self.cre_coords[cre_inclusions, 0]
        cre_ends = self.cre_coords[cre_inclusions, 1]
        X = self.adjust_by_distance(X, self.TSSs[i], cre_starts, cre_ends)
        #sum across cres for states
        X = np.sum(X, axis=0)
        Y = self.exp_values[i, self.lessone_range]
        lin_reg = self.linear_regression(X,Y)
        yhat = self.find_yhat(lin_reg['fit_model'], X)
        self.mse_array[i, iter_val+1] = metrics.mean_squared_error(Y, yhat)

    def set_initial_betas(self):
        '''do initial regression of proportion of states within two-sided 75kbp window against expression to find initial coefficients for all states except 0
        set contribution of state 0 to be 0'''
        initial_coeffs = self.linear_regression(self.TSS_window_props, self.exp_values)['coeffs']
        self.initial_coeffs = np.hstack(([0], initial_coeffs)).reshape((1,-1))

    def find_initial_weighted_sum(self):
        '''find weighted sum of ccRE props with initial_coeffs
        given a 3 dimensional array ccREn * cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                1 dimensional array stateN with c_i equal to the initial coefficient for state_i
        return a 2 dimensional array ccREn * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''

        self.cre_weighted_sum = np.sum(self.cre_props * self.initial_coeffs, axis=2)

    def find_gen_weighted_sum(self, subset_props, coeffs):
        '''find weighted sum of general ccRE props with general coeffs
        given a 3 dimensional array ccREn * cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                1 dimensional array stateN with c_i equal to the coefficient for state_i
        return a 2 dimensional array ccREn * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''
        weighted_sum = np.sum(subset_props * coeffs, axis=2)
        return weighted_sum

    def compute_adj_distance(self, starts, stops, TSS):
        locsTA = (starts + stops) // 2
        distance = np.abs(locsTA- TSS) #using abs because don't care if upstream or downstream
        adjusted_distance = distance**self.dist_gamma
        adjusted_distance[np.where(adjusted_distance == 0)[0]] = 1 #anywhere that distance is 0, set to identity.
        return adjusted_distance

    def compute_distance(self, starts, stops, TSS):
        midLocs = (starts + stops) //2
        distance = midLocs - TSS
        return distance

    def adjust_by_distance(self, to_adjust, TSS, starts, stops):
        adj_dist = self.compute_adj_distance(starts, stops, TSS)
        Y = np.ones(to_adjust.shape[0]) #adjustment array
        Y /= adj_dist
        if to_adjust.ndim == 2:
            adjusted = np.multiply(to_adjust, Y.reshape(-1,1))
        elif to_adjust.ndim == 3:
            adjusted = np.multiply(to_adjust, Y.reshape(-1,1,1))
        return adjusted

    def find_initial_pairings(self, cre_dist, correlation, i):
        '''for given TSS:
            1) find CREs within distance of interest
            2) adjust the CREs weighted sum of proporitions by their distance from the TSS
            3) for each CRE within distance of interest
                3.1) find correlation of expression with this adjusted weighted sum
                      -> using spearmanr so no distribution assumptions are made and monotonicity is consistent with linearity which we hold this relationship to be linear
                3.2) decide whether to initally pair or not '''

        TSS = self.TSSs[i]
        #find CREs within distance of interest using containment
        windowMin = max(0, TSS - cre_dist)
        windowMax = min(TSS + cre_dist, self.chrSizes[self.chrom])
        CREs_within = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
        if np.sum(CREs_within) == 0: #no pairing possible
            return
        CREs_within_weighted_sum = self.cre_weighted_sum[CREs_within]
        CREs_within_starts = self.cre_coords[CREs_within, 0]
        CREs_within_stops = self.cre_coords[CREs_within, 1]

        #adjust their weighted_sums by their distance
        CREs_within_adjusted = self.adjust_by_distance(CREs_within_weighted_sum, TSS, CREs_within_starts, CREs_within_stops)

        #for each CRE within distance of interest (and as long as its state isn't 0 across all 12 cell types)
        #find correlation of expression with this adjusted weighted sum
        corr_matrix, pvalues = stats.spearmanr(self.exp_values[i].reshape((1,-1)), CREs_within_adjusted, axis=1)
        if isinstance(corr_matrix, np.ndarray):
            exp_cre_corr = corr_matrix[0,1:] #take slice of correlations between that TSS and the CREs within the window
        else:
            exp_cre_corr = corr_matrix
        subset_no_nan = ~np.isnan(exp_cre_corr) #find loc of NaN values so correlations can be compared to threshold without these
        if np.sum(subset_no_nan) == 0: #no pairing possible
            return
        #compare abs value of correlation of expression with threshold value
        corr_passes = np.abs(exp_cre_corr[subset_no_nan]) >= correlation
        if np.sum(corr_passes) == 0: #no pairing possible
            return

        #subset self.creIndex_range by the 3 boolean masks to mark which CREs are initially paired in self.pairing_array
        self.pairing_array[i, self.creIndex_range[CREs_within][subset_no_nan][corr_passes], 0] = 1

    def get_dim_of_pairing(self, tss_i, iter_val): #iter_val of -1 for initial_pairing
        return np.sum(self.pairing_array[tss_i, :, iter_val+1] == 1) #will return number of CREs paired in that slice

    def refine_pairs(self, tss_i, iter_val):
        return 0

    def utilize_IC(self):
        return 0

    def bootstrap(self):
        return 0



    def drive_pairing(self, dist_gamma, lessone, cre_dist, correlation):
        '''for each TSS - Model Selection:
            1) Find its initial pairings
                1.1) If no pairings - notate and move on
                1.2) If pairings - move to 2
            2) Iteratively refine its pairings
            3) utlize IC on refined pairing models to find a model that hopefully isn't overfit
        '''
        #instead of looping through all TSS to find the initial pairing and then looping through them again to refine the pairings, will loop through a single time.
        self.dist_gamma = dist_gamma
        self.lessone = lessone
        self.lessone_range = np.r_[np.arange(self.lessone),
                                    np.arange(self.lessone+1, self.cellN)]
        #writeNoPairings = open('noPairings_possible_{}_{}.txt'.format(correlation, self.chrom), 'a')
        for i in range(self.tssN):
            self.find_initial_pairings(cre_dist, correlation, i)
            if self.get_dim_of_pairing(i, -1) == 0: #notate no pairings
                #writeNoPairings.write('{}\tTSS: {}\n'.format(self.chrom, TSS))
                continue
            else: #From this point forward, do not include lessone cell type in ANYTHING
                #find MSE for initial pairing
                self.find_MSE(i, -1)
                #iteratively refine
                for j in range(self.iter):
                    self.refine_pairs(i, j)

        #writeNoPairings.close()
main()
