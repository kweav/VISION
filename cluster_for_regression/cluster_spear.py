#!/usr/bin/env python3

import numpy as np
import argparse as ap
import os
from sklearn import linear_model
import itertools
from scipy import stats

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
        model.compute_spearmanr()
        model.save_distances()

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
    parser.add_argument('--chroms', action='store', nargs='+', dest='chroms', type=str, default='chr1', help='use all if you want to run all')
    parser.add_argument('--which_with', action='store', dest='which_with', type=str, default='within', help='{TT, EE, within, not}')
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

        self.beta_e = np.hstack(([0], self.linear_fit(self.TSS_window_props, self.exp_values))).reshape((1,-1))

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
        given a 2 dimensional array 1 ccRE  cellN * stateN with p_ij equal to the proporition of the 1 ccRE for cellType_i in state_j
                1 dimensional array stateN with c_i equal to the initial coefficient for state_i
        return a 1 dimensional array cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl
        given a 3 dimensional array k ccRE * cellN * stateN with p_kij equal to the proportion of ccRE_k for cellType_i in state_j
                1 dimensional array stateN with c_i equal to the initial coefficient for state_i
        return a 1 dimensional array cellN with w_ij equal to the sum of l=0 to 26 of c_l*p_ijl'''
        cre_weighted_sum = np.sum(E * self.beta_e, axis=(E.ndim-1))
        return cre_weighted_sum

    def compute_spearmanr(self):
        self.corr_matrix, pvalues = stats.spearmanr(self.exp_values, self.find_weighted_sum(self.cre_props), axis=1)

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
        C = E1 - E2
        distance = self.magnitude(C)
        return distance

    def distance_T_T(self, tss_i, tss_j):
        '''Euclidean distance of the difference between two vectors'''
        C = self.exp_values[tss_i] - self.exp_values[tss_j]
        T_distance = self.magnitude(C)
        '''Euclidean distance of the difference between two matrices'''
        Z = self.TSS_window_props[tss_i] - self.TSS_window_props[tss_j]
        P_distance = self.magnitude(Z)
        return P_distance, T_distance

    def angle_T_T(self, tss_i, tss_j):
        '''use dot product and magnitude of vector to find angel b/w two vectors'''
        dot = self.dot_product(self.exp_values[tss_i], self.exp_values[tss_j])
        mag_i = self.magnitude(self.exp_values[tss_i])
        mag_j = self.magnitude(self.exp_values[tss_j])
        y = dot/(mag_i * mag_j)
        if y <= 1 and y >= -1:
            angle = np.arccos(y)
        else:
            angle = np.nan
        return angle

    def distance_T_E(self, tss_i, E):
        '''Find weighted sum/dot product'''
        weighted_sum = self.find_weighted_sum(E)
        C = self.exp_values[tss_i] - weighted_sum
        '''Euclidean distance of the difference between the dot product (of the E matrix and betas vector) and TPM'''
        DE_distance = self.magnitude(C)
        '''Euclidean distance of the difference between the two matrices'''
        Z = self.TSS_window_props[tss_i] - E
        TE_distance = self.magnitude(Z)
        return DE_distance, TE_distance

    def angle_T_E(self, tss_i, E):
        '''Find weighted sum/dot product'''
        weighted_sum = self.find_weighted_sum(E)
        dot = self.dot_product(self.exp_values[tss_i], weighted_sum)
        mag_i = self.magnitude(self.exp_values[tss_i])
        mag_w = self.magnitude(weighted_sum)
        y = dot/(mag_i * mag_w)
        if y <= 1 and y >= -1:
            angle = np.arccos(y)
        else:
            angle = np.nan
        return angle

    def report_T_T_values(self, Pdistance, Tdistance, angle):
        for arg, argT in zip([Pdistance, Tdistance, angle], ['Pdistance', 'Tdistance', 'angle']):
            toWriteTo = 'outputTT_{}_{}.txt'.format(argT, self.chrom)
            np.savetxt(toWriteTo, arg)

    def report_T_E_W(self, DEdistance, TEdistance, angle):
        for arg, argT in zip([DEdistance, TEdistance, angle], ['DEdistance', 'TEdistance', 'angle']):
            toWriteTo = 'outputTEW_{}_{}.txt'.format(argT, self.chrom)
            np.savetxt(toWriteTo, arg)

    def report_T_E_NW(self, DEdistance, TEdistance, angle):
        for arg, argT in zip([DEdistance, TEdistance, angle], ['DEdistance', 'TEdistance', 'angle']):
            toWriteTo = 'outputTENW_{}_{}.txt'.format(argT, self.chrom)
            np.savetxt(toWriteTo, arg)

    def report_E_E_values(self, distance):
        toWriteTo = 'outputEE_distance_{}.txt'.format(self.chrom)
        np.savetxt(toWriteTo, distance)

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

    #def save_distances(self, which_with):
    def save_distances(self):
        to_save = np.full((self.tssN, self.creM, 3), np.nan) #3 dim 0 index is little d distance, 1 index is Cap D distance assuming only two clustered, 2 index is spearmanr
        for i in range(self.tssN):
            CREs_within = self.find_cres_within(i)
            withinN = np.sum(CREs_within)
            self.distances = np.zeros((withinN+1, withinN), dtype=np.float32) #save all distances so can compute cap D
            cre_j_to_index = {} #these indices are zero-based; for rows, must add one; for columns good as is
            index_to_cre_j = {} #these indices are zero-based; for rows, must add one; for columns good as is
            '''find the distances between the TSS and each CRE...'''
            for enum_j, cre_j in enumerate(np.arange(self.creM)[CREs_within]):
                DEdistance, TEdistance = self.distance_T_E(i, self.cre_props[cre_j])
                adjusted_distance = DEdistance + TEdistance
                to_save[i, cre_j, 0] = adjusted_distance
                to_save[i, cre_j, 2] = self.corr_matrix[i, self.tssN+cre_j]
                '''only storing in upper right of matrix'''
                self.distances[0, enum_j] = adjusted_distance
                cre_j_to_index[cre_j] = enum_j
                index_to_cre_j[enum_j] = cre_j

            '''...as well as then the distances between each CRE and each CRE'''
            for cre_a,cre_b in itertools.combinations(np.arange(self.creM)[CREs_within], 2):
                distance = self.distance_E_E(self.cre_props[cre_a], self.cre_props[cre_b])
                '''only storing in upper right of matrix'''
                enum_a, enum_b = cre_j_to_index[cre_a], cre_j_to_index[cre_b]
                self.distances[min(enum_a, enum_b)+1, max(enum_a, enum_b)] = distance

            for cre_j in np.arange(self.creM)[CREs_within]:
                capital_dist = self.get_D_ij(cre_j_to_index[cre_j])
                to_save[i, cre_j, 1] = capital_dist

        npz_file = open('{}_corrs_and_dists.npz'.format(self.chrom), 'wb')
        np.savez(npz_file, corrs_and_dists=to_save)
        npz_file.close()

        # if which_with == "TT":
        #     TT_Pdistances = []
        #     TT_Tdistances = []
        #     TT_angles = []
        #     for a, b in itertools.combinations(range(self.tssN), 2):
        #         Pdistance, Tdistance = self.distance_T_T(a, b)
        #         TT_Pdistances.append(Pdistance)
        #         TT_Tdistances.append(Tdistance)
        #         angle = self.angle_T_T(a, b)
        #         TT_angles.append(angle)
        #     TT_Pdistances = np.array(TT_Pdistances).reshape(-1,)
        #     TT_Tdistances = np.array(TT_Tdistances).reshape(-1,)
        #     TT_angles = np.array(TT_angles).reshape(-1,)
        #     self.report_T_T_values(TT_Pdistances, TT_Tdistances, TT_angles)
        #
        # if which_with == 'EE':
        #     distances = []
        #     for a,b in itertools.combinations(range(self.creM), 2):
        #         distance = self.distance_E_E(self.cre_props[a], self.cre_props[b])
        #         distances.append(distance)
        #     distances = np.array(distances).reshape(-1,)
        #     self.report_E_E_values(distances)
        #
        # if which_with == 'within':
        #     DE_dist_within = []
        #     TE_dist_within = []
        #     angle_within = []
        #     for i in range(self.tssN):
        #         CREs_within = self.find_cres_within(i)
        #         for j in np.arange(self.creM)[CREs_within]:
        #             DEdistance, TEdistance = self.distance_T_E(i, self.cre_props[j])
        #             DE_dist_within.append(DEdistance)
        #             TE_dist_within.append(TEdistance)
        #             angle = self.angle_T_E(i, self.cre_props[j])
        #             angle_within.append(angle)
        #     DE_dist_within = np.array(DE_dist_within).reshape(-1,)
        #     TE_dist_within = np.array(TE_dist_within).reshape(-1,)
        #     angle_within = np.array(angle_within).reshape(-1,)
        #     self.report_T_E_W(DE_dist_within, TE_dist_within, angle_within)
        #
        # if which_with == 'not':
        #     DE_dist_nw = []
        #     TE_dist_nw = []
        #     angle_nw = []
        #     for i in range(self.tssN):
        #         CREs_within = self.find_cres_within(i)
        #         CREs_notwithin = ~CREs_within
        #         for j in np.arange(self.creM)[CREs_notwithin]:
        #             DEdistance, TEdistance = self.distance_T_E(i, self.cre_props[j])
        #             DE_dist_nw.append(DEdistance)
        #             TE_dist_nw.append(TEdistance)
        #             angle = self.angle_T_E(i, self.cre_props[j])
        #             angle_nw.append(angle)
        #     DE_dist_nw = np.array(DE_dist_nw).reshape(-1,)
        #     TE_dist_nw = np.array(TE_dist_nw).reshape(-1,)
        #     angle_nw = np.array(angle_nw).reshape(-1,)
        #     self.report_T_E_NW(DE_dist_nw, TE_dist_nw, angle_nw)

main()
