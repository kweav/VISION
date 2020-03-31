#!/usr/bin/env python3

'''Neighbor Joining procedure to cluster each TSS with up to n cCREs
Unlike classical neighbor joining, only looking to join CREs with a TSS or a TSS-CRE cluster
D_ij = d_ij - (r_i + r_j); r_i = 1/(|L| - 2) * Sigma_{k in L, not = i,j} d_ik
Will not include state 0 in anything'''

import numpy as np
import argparse as ap
import os
from sklearn import linear_model
import itertools

def main():
    parser = generate_parser()
    args = parser.parse_args()
    setup_file_locs(args, args.where_run, args.other_path)
    setup_threads(args.threads)
    model = regress_sampler(args.train_cre, args.train_tss, args.train_exp, args.cre_dist)

    if args.chroms == 'all':
        chrom_list = []
        for i in range(1, 20):
            chrom_list.append('chr{}'.format(i))
        chrom_list.append('chrX')
    else:
        chrom_list = args.chroms

    for chrom in chrom_list:
        model.subset_based_on_chrom(chrom)
        model.run_clustering(args.n_thresh, args.dist_thresh)

def generate_parser():
    parser = ap.ArgumentParser(description='looking at distances for clustering')

    parser.add_argument('--where', action='store', dest='where_run', type=str, default='mine', help='{comp, mine, marcc, other}; adds path to files; if "other" used provide the other path in the --otherpath argument')
    parser.add_argument('--threads', action='store', dest='threads', type=str, default="1")
    parser.add_argument('--otherpath', action='store', dest='other_path', type=str, default='NA', help='use to give other path if not one of the 3 previously specified. MUST be the same for ALL input files')
    parser.add_argument('-d', '--dist', action='store', type=int, dest='cre_dist', default='1000000')
    parser.add_argument('--train_cre', action='store', dest='train_cre', type=str, default='train_cre_state_prop.npz')
    parser.add_argument('--test_cre', action='store', dest='test_cre', type=str, default='test_cre_state_prop.npz')
    parser.add_argument('--train_tss', action='store', dest='train_tss', type=str, default='trainTSS_window_state_prop.npz')
    parser.add_argument('--test_tss', action='store', dest='test_tss', type=str, default='testTSS_window_state_prop.npz')
    parser.add_argument('--train_exp', action='store', dest='train_exp', type=str, default='trainTPM.npz')
    parser.add_argument('--test_exp', action='store', dest='test_exp', type=str, default='testTPM.npz')
    parser.add_argument('--chroms', action='store', nargs='+', dest='chroms', type=str, default=['chr1'], help='use all if you want to run all')
    parser.add_argument('--n_thresh', action='store', dest='n_thresh', type=int, default=30)
    parser.add_argument('--dist_thresh', action='store', dest='dist_thresh', type=float, default=10000)
    return parser

def setup_threads(threads):
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["OPENBLAS_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = threads
    os.environ["NUMEXPR_NUM_THREADS"] = threads

def setup_file_locs(args, where_run, other_path):
    argumentToAdd = {'mine': '/Users/kateweaver/taylorLab/VISION/regress_dist_gibbs/inputs/',
                     'marcc': 'NA',
                     'comp': '/project/vision/Target_Genes/npz_inputs/',
                     'other': other_path}

    args.train_cre = argumentToAdd[where_run] + args.train_cre
    args.test_cre = argumentToAdd[where_run] + args.test_cre
    args.train_exp = argumentToAdd[where_run] + args.train_exp
    args.test_exp = argumentToAdd[where_run] + args.test_exp
    args.train_tss = argumentToAdd[where_run] + args.train_tss
    args.test_tss = argumentToAdd[where_run] + args.test_tss

class regress_sampler():
    def __init__(self, train_cre, train_tss, train_exp, cre_dist):
        self.cre_dist = cre_dist
        self.exp_values_all, self.cellIndex, self.cell_to_index, self.TSS_chr_all, self.TSSs_all = self.load_expression(train_exp)
        self.cellN = self.cellIndex.shape[0]
        self.TSS_window_props_all, self.TSS_window_chr_all = self.load_TSS_window_states(train_tss)
        self.cre_props_all, self.cre_chr_all, self.cre_coords_all = self.load_cre_states(train_cre)
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
        return TSS_window_props_valid, TSS_window_chr

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

    def linear_fit(self, X, Y, intercept=False):
        X = X.reshape(-1, self.stateN)[:,1:]
        Y = Y.reshape(-1,)
        fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X, Y)
        return fit_model.coef_

    def subset_based_on_chrom(self, chrom):
        self.chrom = chrom
        where_chr = np.where(self.TSS_chr_all == self.chrom)[0]
        self.exp_values = self.exp_values_all[where_chr]
        self.TSSs = self.TSSs_all[where_chr]
        self.tssN = self.TSSs.shape[0]
        self.TSS_window_props = self.TSS_window_props_all[where_chr]

        where_chr = np.where(self.cre_chr_all == self.chrom)[0]
        self.cre_props = self.cre_props_all[where_chr]
        self.cre_coords = self.cre_coords_all[where_chr]
        self.creM = self.cre_props.shape[0]
        self.creIndex_range = np.arange(self.creM)

        #self.beta_e = np.hstack(([0], self.linear_fit(self.TSS_window_props, self.exp_values))).reshape((1,-1))
        self.beta_e = self.linear_fit(self.TSS_window_props, self.exp_values).reshape((1, -1)) #doesn't include hstack of 0 for state 0

    def cut_quiescent_state(self, toCut):
        '''if given a 2 dimensional array (cellN * stateN): return a 2 dimensional array (cellN * stateN-1) without the quiescent state
        if given a 3 dimensional array (elementN * cellN * stateN): return a 3 dimensional array (elementN * cellN * stateN-1) without the quiescent state'''
        if toCut.ndim == 2:
            return toCut[:,1:]
        elif toCut.ndim == 3:
            return toCut[:, :, 1:]
        else:
            return toCut[1:]

    def find_cres_within(self, i):
        '''#find CREs within distance of interest using containment
        start location of CRE is less than or equal to window end
        AND end locatiion of CRE is greater than or equal to window beginning'''
        TSS = self.TSSs[i]
        windowMin = max(0, TSS - self.cre_dist)
        windowMax = min(TSS + self.cre_dist, self.chrSizes[self.chrom])
        CREs_within = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
        return CREs_within

    def find_weighted_sum(self, E):
        '''find weighted sum of ccRE props with coeffs
        given a 2 dimensional array 1 ccRE  cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                1 dimensional array stateN with c_i equal to the initial coefficient for state_i
        return a 1 dimensional array cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''
        cre_weighted_sum = np.sum(E * self.beta_e, axis=1)
        return cre_weighted_sum

    def magnitude(self, C):
        magnitude = np.sqrt(np.sum(np.square(C)))
        if magnitude == 0:
            magnitude += 0.1
        return magnitude

    def dot_product(self, a, b):
        dot = np.dot(a,b)
        return dot

    def distance_E_E(self, E1, E2):
        '''Euclidean distance of the difference between two matrices'''
        C = self.cut_quiescent_state(E1) - self.cut_quiescent_state(E2)
        distance = self.magnitude(C)
        return distance

    def distance_T_E(self, tss_i, E):
        '''Find weighted sum/dot product'''
        weighted_sum = self.find_weighted_sum(self.cut_quiescent_state(E))
        C = self.exp_values[tss_i] - weighted_sum
        '''Euclidean distance of the difference between the dot product (of the E matrix and betas vector) and TPM'''
        DE_distance = self.magnitude(C)
        '''Euclidean distance of the difference between the two matrices'''
        Z = self.cut_quiescent_state(self.TSS_window_props[tss_i]) - self.cut_quiescent_state(E)
        TE_distance = self.magnitude(Z)
        return DE_distance, TE_distance

    def angle_T_E(self, tss_i, E):
        '''Find weighted sum/dot product'''
        weighted_sum = self.find_weighted_sum(self.cut_quiescent_state(E))
        dot = self.dot_product(self.exp_values[tss_i], weighted_sum)
        mag_i = self.magnitude(self.exp_values[tss_i])
        mag_w = self.magnitude(weighted_sum)
        y = dot/(mag_i * mag_w)
        if y <= 1 and y >= -1:
            angle = np.arccos(y)
        else:
            angle = np.nan
        return angle

    def get_r_i(self, index_j):
        bool_all_but_j = np.full((self.distances.shape[1]), True)
        bool_all_but_j[index_j] = False
        summed = np.sum(self.distances[0, bool_all_but_j])
        return summed/(bool_all_but_j.shape[0] - 1)

    def get_r_j(self, index_j):
        summand = np.hstack((self.distances[1:index_j+1, index_j].reshape(-1,), self.distances[index_j+1, index_j+1:].reshape(-1,)))
        return np.sum(summand)/summand.shape[0]

    def get_D_ij(self, index_j):
        sum = self.distances[0, index_j] - self.get_r_j(index_j) - self.get_r_i(index_j)
        return sum

    def find_initial_little_d(self, tss_i):
        within_bool = self.find_cres_within(tss_i)
        withinN = np.sum(within_bool)
        print('initial withinN ', withinN, flush=True)
        self.distances = np.zeros((withinN+1, withinN), dtype=np.float32)
        print('initial little d ', self.distances.shape, flush=True)

        cre_j_to_index = {} #these indices are zero-based; for rows, must add one; for columns good as is
        index_to_cre_j = {} #these indices are zero-based; for rows, must add one; for columns good as is

        '''find the distances between the TSS and each CRE...'''
        for enum_j, cre_j in enumerate(np.arange(self.creM)[within_bool]):
            DEdistance, TEdistance = self.distance_T_E(tss_i, self.cre_props[cre_j])
            adjusted_distance = DEdistance + TEdistance
            '''only storing in upper right of matrix'''
            self.distances[0, enum_j] = adjusted_distance
            cre_j_to_index[cre_j] = enum_j
            index_to_cre_j[enum_j] = cre_j

        '''...as well as then the distances between each CRE and each CRE'''
        for cre_a,cre_b in itertools.combinations(np.arange(self.creM)[within_bool], 2):
            distance = self.distance_E_E(self.cre_props[cre_a], self.cre_props[cre_b])
            '''only storing in upper right of matrix'''
            enum_a, enum_b = cre_j_to_index[cre_a], cre_j_to_index[cre_b]
            self.distances[min(enum_a, enum_b)+1, max(enum_a, enum_b)] = distance

        return within_bool, withinN, cre_j_to_index, index_to_cre_j

    def find_initial_capital_D(self, dist_thresh, within_bool, index_dict):
        min_D_ij = dist_thresh
        index_min = -1
        for cre_j in np.arange(self.creM)[within_bool]:
            capital_dist = self.get_D_ij(index_dict[cre_j])
            if capital_dist < min_D_ij:
                min_D_ij = capital_dist
                index_min = index_dict[cre_j]
        return index_min

    def update_little_d(self, withinN, current_n, index_min):
        print('withinN update ', withinN, flush=True)
        new_distances = np.zeros((withinN+1-current_n, withinN-current_n))
        '''copy information for not clustered distances'''
        new_distances[1:index_min+1, 0:index_min] = self.distances[1:index_min+1, 0:index_min]
        new_distances[index_min+1:, 0:index_min] = self.distances[index_min+2:, 0:index_min] #might already be zeros and therefore an unnecessary step????
        new_distances[1:index_min+1, index_min:] = self.distances[1:index_min+1, index_min+1:]
        new_distances[index_min+1:, index_min:] = self.distances[index_min+2:, index_min+1:]
        '''update information for TSS/CRE cluster'''
        d_im = np.hstack((self.distances[0,0:index_min], self.distances[0,index_min+1:]))
        d_jm = np.hstack((self.distances[1:index_min+1, index_min].reshape(-1,), self.distances[index_min+1, index_min+1:].reshape(-1,)))
        d_ij = self.distances[0, index_min]
        updated_cluster = 1/2*((d_im + d_jm) - d_ij)
        new_distances[0, :] = updated_cluster
        self.distances = new_distances
        print('update little d ', self.distances.shape, flush=True)

    def update_index_dict(self, tss_i, within_bool, clustered):
        clustered_bool = clustered[tss_i, within_bool] == 0
        print(np.sum(clustered_bool), flush=True)
        print(np.sum(clustered[tss_i, within_bool]), flush=True)
        cre_j_to_index = {} #these indices are zero-based; for rows, must add one; for columns good as is
        index_to_cre_j = {} #these indices are zero-based; for rows, must add one; for columns good as is
        for enum_j, cre_j in enumerate(np.arange(self.creM)[within_bool][clustered_bool]):
            print(cre_j, enum_j, flush=True)
            cre_j_to_index[cre_j] = enum_j
            index_to_cre_j[enum_j] = cre_j
        return clustered_bool, cre_j_to_index, index_to_cre_j

    def find_next_capital_D(self, dist_thresh, within_bool, index_dict, clustered_bool):
        min_D_ij = dist_thresh
        index_min = -1
        for cre_j in np.arange(self.creM)[within_bool][clustered_bool]:
            capital_dist = self.get_D_ij(index_dict[cre_j])
            if capital_dist < min_D_ij:
                min_D_ij = capital_dist
                index_min = index_dict[cre_j]
        return index_min

    def run_clustering(self, n_thresh, dist_thresh, clustering_available=False):
        clustered = np.zeros((self.tssN, self.creM), dtype=np.int32)
        for tss_i in range(self.tssN):
            current_n = 0
            '''find little d distances'''
            within_bool, withinN, cre_j_to_index, index_to_cre_j = self.find_initial_little_d(tss_i)
            if withinN > 0:
                clustering_available = True
                '''Find close and far capital D distances'''
                index_min = self.find_initial_capital_D(dist_thresh, within_bool, cre_j_to_index)
                if index_min != -1:
                    print('got here', index_min, type(index_min), flush=True)
                    print(type(clustered), flush=True)
                    print(clustered.dtype, flush=True)
                    print(clustered.shape, flush=True)
                    print(clustered[tss_i].shape, flush=True)
                    print(clustered[tss_i, within_bool].shape, flush=True)
                    print('before set ', clustered[tss_i, within_bool][index_min], flush=True)
                    clustered[tss_i, within_bool][index_min] =
                    print('after set ', clustered[tss_i, within_bool][index_min], flush=True)
                    current_n += 1
                    quit()
                else:
                    clustering_available = False
                while (current_n < n_thresh) and clustering_available:
                    self.update_little_d(withinN, current_n, index_min)
                    clustered_bool, updated_cre_j_to_index, updated_index_to_cre_j = self.update_index_dict(tss_i, within_bool, clustered)
                    index_min = self.find_next_capital_D(dist_thresh, within_bool, updated_cre_j_to_index, clustered_bool)
                    if index_min != -1:
                        clustered[tss_i, within_bool][clustered_bool][index_min] = 1
                        current_n += 1
                    else:
                        clustering_available = False
                print(tss_i, np.sum(clustered[tss_i]), flush=True)
                quit()
            else:
                continue

main()
