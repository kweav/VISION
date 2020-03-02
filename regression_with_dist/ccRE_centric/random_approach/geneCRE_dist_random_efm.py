#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model, metrics
import os
from scipy import stats
from itertools import chain, combinations
import datetime


def main():
    parser = generate_parser()
    args = parser.parse_args()
    setup_file_locs(args, args.where_run, args.other_path)
    setup_threads(args.threads)
    model = regress_gene_cre(args.cre_state_file, args.tss_state_file, args.exp_file, args.output_pair_file, args.output_beta_file)

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
        model.find_initial_weighted_sum()
        model.build_pairing_array()
        model.build_beta_array()
        model.drive_pairing(args.dist_gamma, args.cre_dist, args.correlation, args.iterations, args.psubset)

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
    parser.add_argument('--p', '--psubset', action='store', dest='psubset', type=int, default=200)
    parser.add_argument('-i', '--iterations', action='store', dest='iterations', type=int, default=1000)
    #parser.add_argument('-l', '--lessone', action='store', dest='lessone', type=int, default=0, help='Cell type to leave out with 0 indexing. 0 is default')
    parser.add_argument('--otherpath', action='store', dest='other_path', default='NA', help='use to give other path if not one of the 3 previously specified. MUST be the same for ALL 3 input files')
    parser.add_argument('--correlation', action='store', dest='correlation', type=float, default=0.4)
    parser.add_argument('--output_pair', action='store', dest='output_pair_file', type=str, default='dist_regression_randAp_pairs_{}.txt')
    parser.add_argument('--output_beta', action='store', dest='output_beta_file', type=str, default='dist_regression_randAp_betas_{}.txt')
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
                     'comp': '/project/vision/Target_Genes/random_with_dist/inputs/',
                     'other': other_path}

    args.exp_file = argumentToAdd[where_run] + args.exp_file
    args.cre_state_file = argumentToAdd[where_run] + args.cre_state_file
    args.tss_state_file = argumentToAdd[where_run] + args.tss_state_file

class regress_gene_cre():
    def __init__(self, cre_state_file, tss_state_file, exp_file, output_pair_file, output_beta_file):
        '''
        self.output_pair_file:
        self.output_beta_file:
        self.exp_values_all:
        self.cellIndex:
        self.cell_to_index:
        self.TSS_chr_all:
        self.TSSs_all:
        self.cellN: REPORT
        self.TSS_window_props_all:
        self.TSS_window_chr_all:
        self.stateN: REPORT
        self.chrSizes:
        self.chrom:
        self.exp_values:
        self.TSSs:
        self.TSS_window_props:
        self.cre_props:
        self.cre_coords:
        self.creM: REPORT
        self.creIndex_range:
        self.m_thresh: REPORT
        self.s_thresh: REPORT
        self.group: REPORT
        self.tssN: REPORT
        self.initial_coeffs:
        self.cre_weighted_sum:
        self.pairing_array:
        self.beta_array:
        self.passing_corrs:
        self.dist_gamma: REPORT
        self.iter: REPORT
        self.psubset: REPORT
        self.correlation: REPORT
        self.cre_dist: REPORT
        self.lessone:
        self.lessone_range:
        '''
        self.output_pair_file = output_pair_file
        self.output_beta_file = output_beta_file
        self.exp_values_all, self.cellIndex, self.cell_to_index, self.TSS_chr_all, self.TSSs_all = self.load_expression(exp_file)
        self.cellN = self.cellIndex.shape[0]
        self.TSS_window_props_all, self.TSS_window_chr_all = self.load_TSS_window_states(tss_state_file)
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
        TSS_window_props = npzfile['props'].astype(np.float32)
        TSS_cellIndex = npzfile['cellIndex']
        valid = self.find_valid_cellTypes(TSS_cellIndex)
        TSS_window_props_valid = TSS_window_props[:,valid,:]
        TSS_window_index = npzfile['ccREIndex']
        TSS_window_chr = TSS_window_index[:,0]
        return TSS_window_props_valid, TSS_window_chr

    def load_cre_states(self, cre_state_file):
        npzfile = np.load(cre_state_file, allow_pickle = True)
        cre_props = npzfile['props'].astype(np.float32)
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

    def build_pairing_array(self):
        self.pairing_array = np.zeros((self.creM, self.cellN, 3), dtype=np.bool) #3D array of creM * cellN * 3. Last dimension is for initial, best, and current

    def build_beta_array(self):
        self.beta_array = np.zeros((self.cellN+1, self.stateN), dtype=np.float32) #store betas for initial pairing 0 in first dim and then self.lessone+1 after that


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

    def find_initial_pairings(self, i):
        '''for given TSS:
            1) find CREs within distance of interest
            2) adjust the CREs weighted sum of proporitions by their distance from the TSS
            3) for each CRE within distance of interest
                3.1) find correlation of expression with this adjusted weighted sum
                      -> using spearmanr so no distribution assumptions are made and monotonicity is consistent with linearity which we hold this relationship to be linear
                3.2) decide whether to initally pair or not '''

        TSS = self.TSSs[i]
        #find CREs within distance of interest using containment
        windowMin = max(0, TSS - self.cre_dist)
        windowMax = min(TSS + self.cre_dist, self.chrSizes[self.chrom])
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
        corr_passes = np.abs(exp_cre_corr[subset_no_nan]) >= self.correlation
        if np.sum(corr_passes) == 0: #no pairing possible
            return

        #subset self.creIndex_range by the 3 boolean masks to mark which CREs are initially paired in self.pairing_array
        self.pairing_array[self.creIndex_range[CREs_within][subset_no_nan][corr_passes], :, 0] = 1
        self.passing_corrs = exp_cre_corr[subset_no_nan][corr_passes]

    def find_yhat(self, linear_model, X):
        X = X.reshape(-1, self.stateN)[:,1:]
        return linear_model.predict(X)

    def find_MSE_and_betas(self, tss_i, last_dim): #if calling for inital pairing last_dim = 0, if for best last_dim = 1, if for current last_dim=2
        cre_inclusions = self.creIndex_range[self.pairing_array[:,self.lessone, last_dim] == 1]
        X = self.cre_props[cre_inclusions, :, :][:,self.lessone_range,:]
        #adjust by distance
        cre_starts = self.cre_coords[cre_inclusions, 0]
        cre_ends = self.cre_coords[cre_inclusions, 1]
        X = self.adjust_by_distance(X, self.TSSs[tss_i], cre_starts, cre_ends)
        #sum across cres for states
        X = np.sum(X, axis=0)
        Y = self.exp_values[tss_i, self.lessone_range]
        lin_reg = self.linear_regression(X,Y)
        yhat = self.find_yhat(lin_reg['fit_model'], X)
        return metrics.mean_squared_error(Y, yhat), np.hstack(([0], lin_reg['coeffs'])).reshape((1,-1))

    def get_dim_of_pairing(self, lessone, last_dim): #if calling for initial pairing last_dim = 0, if calling for best so far last_dim = 1, if for current last_dim =2
        return np.sum(self.pairing_array[:, lessone, last_dim] == 1) #will return number of CREs paired in that slice

    def utilize_IC(self, MSE, num_examples, pairing_dim):
        BIC = (-2 * MSE) + (np.log(num_examples) * self.get_dim_of_pairing(self.lessone, pairing_dim))
        return BIC

    def refine_pairs(self, tss_i, inclusions_to_test, given_bestIC):
        bestIC = given_bestIC
        for n in range(inclusions_to_test.shape[0]):
            self.pairing_array[inclusions_to_test[n], self.lessone, 2] = 1 #set current pairing
            current_MSE, current_betas = self.find_MSE_and_betas(tss_i, 2)
            current_IC = self.utilize_IC(current_MSE, self.cellN-1, 2)
            if current_IC < bestIC:
                bestIC = current_IC
                self.pairing_array[:, self.lessone, 1] = np.copy(self.pairing_array[:, self.lessone, 2])
                self.beta_array[self.lessone+1, :] = current_betas
        return bestIC

    def refine_pairs_iter(self, tss_i, generated_powerset, initialIC):
        npz_dict = {}
        prob = np.full(generated_powerset.shape[0], 1/generated_powerset.shape[0], dtype=np.float32)
        # prob = np.full(generated_powerset.shape[0], 1, dtype=np.float32)
        # prob /= np.square(np.abs([len(s) for s in generated_powerset] - np.mean([len(s) for s in generated_powerset])))
        # prob /= np.sum(prob)
        for iter_val in range(self.iter):
            npz_dict["iter{}".format(iter_val)] = np.random.choice(generated_powerset, self.psubset, p=prob, replace=False) #allows for replacement across iterations
            if iter_val == 0:
                bestIC = self.refine_pairs(tss_i, npz_dict["iter{}".format(iter_val)], initialIC)
            else:
                bestIC = self.refine_pairs(tss_i, npz_dict["iter{}".format(iter_val)], bestIC)
            if iter_val % 50 == 0:
                print("Finished iteration {} for {} TSS {}. BestIC currently {}".format(iter_val, self.chrom, self.TSSs[tss_i], bestIC), flush=True)

            # if np.sum(prob != 0) > 0:
            #     random_ints = np.random.choice(np.arange(generated_powerset.shape[0]), self.psubset, p=prob, replace=False)
            #     npz_dict['iter{}'.format(iter_val)] = generated_powerset[random_ints]
            #     '''remove from prob and generated_powerset for next iteration - set their probabilities to 0!'''
            #     prob[random_ints] = 0
            #     prob /= np.sum(prob) #reweight probabilities to sum to one
            #     if iter_val == 0:
            #         bestIC = self.refine_pairs(tss_i, npz_dict["iter{}".format(iter_val)], initialIC)
            #     else:
            #         bestIC = self.refine_pairs(tss_i, npz_dict["iter{}".format(iter_val)], bestIC)
            #     if iter_val % 50 == 0:
            #         print("Finished iteration {} for {} TSS {}".format(iter_val, self.chrom, self.TSSs[tss_i]), flush=True)
            # else:
            #     for i in range(iter_val, self.iter):
            #         npz_dict['iter{}'.format(i)] = np.array([])
            #     break

        npz_file = open('{}_{}_iteration_subsets.npz', 'wb')
        np.savez(npz_file, **npz_dict)
        npz_file.close()
        return bestIC

    def powerset(self, iterable): #avoids the empty set and full set
        "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        return np.array(list(chain.from_iterable(combinations(s, r) for r in range(1, len(s)))))

    def set_up_refinement(self, tss_i, initialIC):
        cre_inclusions = self.creIndex_range[self.pairing_array[:,self.lessone, 0] == 1]
        powerset_of_inclusions = self.powerset(cre_inclusions)
        if powerset_of_inclusions.shape[0] > self.psubset:
            bestIC = self.refine_pairs_iter(tss_i, powerset_of_inclusions, initialIC)
        else:
            bestIC = self.refine_pairs(tss_i, powerset_of_inclusions, initialIC)
        return bestIC

    def find_gen_weighted_sum(self, subset_props, coeffs):
        '''find weighted sum of general ccRE props with general coeffs
        given a 3 dimensional array ccREn * cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                1 dimensional array stateN with c_i equal to the coefficient for state_i
        return a 2 dimensional array ccREn * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl
        given 2 1 dimensional arrays, will return a single scalar value'''
        weighted_sum = np.sum(subset_props * coeffs, axis=subset_props.ndim-1)
        return weighted_sum

    def report_metadata(self):
        if self.group == 0:
            metadata = '#{}\ttssN: {}\tcreM: {}\tcellN: {}\tstateN: {}\tGroup: {}\tcre_dist: {}\tcorrelation: {}\tdist_gamma: {}\titer: {}\tpsubset: {}\n'.format(datetime.datetime.now(), self.tssN, self.creM, self.cellN, self.stateN, self.group, self.cre_dist, self.correlation, self.dist_gamma, self.iter, self.psubset)
        else:
            metadata = '#{}\ttssN: {}\tcreM: {}\tcellN: {}\tstateN: {}\tGroup: {}\tm_thresh: {}\ts_thresh: {}\tcre_dist: {}\tcorrelation: {}\tdist_gamma: {}\titer: {}\tpsubset: {}\n'.format(datetime.datetime.now(), self.tssN, self.creM, self.cellN, self.stateN, self.group, self.m_thresh, self.s_thresh, self.cre_dist, self.correlation, self.dist_gamma, self.iter, self.psubset)
        filePair = open(self.output_pair_file.format(self.chrom), 'a')
        filePair.write(metadata)
        filePair.close()
        fileBeta = open(self.output_beta_file.format(self.chrom), 'a')
        fileBeta.write(metadata)
        fileBeta.close()

    def report_pairing(self, tss_i, bestIC, pairing_flag):
        TSS = self.TSSs[tss_i]
        toWriteTo = open(self.output_pair_file.format(self.chrom), 'a')
        if pairing_flag == False:
            toWriteTo.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.chrom, TSS, '.', '.', '.', '.', '.', '.', '.'))
        else:
            cre_initial_inclusions = self.creIndex_range[self.pairing_array[:,0, 0] == 1]
            cre_initial_starts = self.cre_coords[cre_initial_inclusions, 0]
            cre_initial_stops = self.cre_coords[cre_initial_inclusions, 1]
            numCells_paired = np.sum(self.pairing_array[cre_initial_inclusions,:,1] == 1, axis=1)
            for cre_i in cre_initial_inclusions.shape[0]:
                whichCells_named = ",".join(self.cellIndex[np.where(self.pairing_array[cre_initial_inclusions[cre_i], :, 1] == 1)[0]])
                indexOCs = [self.cell_to_index[cellType] for cellType in whichCells_named.split(',')]
                weighted_sums = [self.find_gen_weighted_sum(self.cre_props[cre_initial_inclusions[cre_i], indexOC, :], self.beta_array[indexOC+1, :]) for indexOC in indexOCs]
                toWriteTo.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(self.chrom, TSS, cre_initial_starts[cre_i], cre_initial_ends[cre_i], numCells_paired[cre_i], whichCells_named, bestIC, self.passing_corrs[cre_i], ",".join(weighted_sums)))
        toWriteTo.close()

    def report_betas(self, i):
        toWriteTo = open(self.output_beta_file.format(self.chrom), 'a')
        TSS = self.TSSs[i]
        for k in range(-1, self.cellN):
            if k == -1:
                ct = 'initial'
            else:
                ct = self.cellIndex[k]
            betas = ",".join(self.beta_array[k, :])
            toWriteTo.write('{}\t{}\t{}\t{}\n'.format(self.chrom, TSS, ct, betas))
        toWriteTo.close()

    def drive_pairing(self, dist_gamma, cre_dist, correlation, iterations, psubset):
        '''for each TSS - Model Selection:
            1) Find its initial pairings
                1.1) If no pairings - notate and move on
                1.2) If pairings - move to 2
            2) Iteratively refine its pairings
            3) utlize IC on refined pairing models to find a model that hopefully isn't overfit
        '''
        #instead of looping through all TSS to find the initial pairing and then looping through them again to refine the pairings, will loop through a single time.
        self.dist_gamma = dist_gamma
        self.iter = iterations
        self.psubset = psubset
        self.correlation=correlation
        self.cre_dist = cre_dist
        self.report_metadata()
        for i in range(self.tssN):
            self.find_initial_pairings(i)
            if self.get_dim_of_pairing(0, 0) == 0: #notate no pairings
                pairing_flag = False
                bestIC = 'NA'
            else: #From this point forward, do not include lessone cell type in ANYTHING
                pairing_flag = True
                for k in range(self.cellN):
                    self.lessone = k
                    self.lessone_range = np.r_[np.arange(self.lessone),
                                                np.arange(self.lessone+1, self.cellN)]
                    #find MSE and IC for initial pairing
                    initial_MSE, initial_betas = self.find_MSE_and_betas(i, 0)
                    self.beta_array[0, :] = initial_betas
                    initialIC = self.utilize_IC(initial_MSE, self.cellN-1, 0)
                    #iteratively refine
                    bestIC = self.set_up_refinement(i, initialIC)
            self.report_pairing(i, bestIC, pairing_flag)
            self.build_pairing_array() #reset the array
            self.report_betas(i)
            self.build_beta_array() #reset the array

main()
