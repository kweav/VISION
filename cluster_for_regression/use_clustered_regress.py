#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model
import os

def main():
    parser = generate_parser()
    args = parser.parse_args()
    setup_file_locs(args, args.where_run, args.other_path)
    setup_threads(args.threads)
    model = regress_sampler(args.train_cre, args.train_tss, args.train_exp)

    if args.chroms == 'all':
        chrom_list = []
        for i in range(1, 20):
            chrom_list.append('chr{}'.format(i))
        chrom_list.append('chrX')
    else:
        chrom_list = args.chroms

    for chrom in chrom_list:
        model.subset_based_on_chrom(chrom, args.train_clustered)
        model.run_fit()


def generate_parser():
    parser = ap.ArgumentParser(description='VISION Gibbs Sampler of parameters for regression of state and distance to assign CREs to genes based on the ability to predict gene expression')

    parser.add_argument('--where', action='store', dest='where_run', type=str, default='mine', help='{comp, mine, marcc, other}; adds path to files; if "other" used provide the other path in the --otherpath argument')
    parser.add_argument('--threads', action='store', dest='threads', type=str, default="1")
    parser.add_argument('--otherpath', action='store', dest='other_path', type=str, default='NA', help='use to give other path if not one of the 3 previously specified. MUST be the same for ALL input files')
    parser.add_argument('--train_cre', action='store', dest='train_cre', type=str, default='train_cre_state_prop.npz')
    parser.add_argument('--test_cre', action='store', dest='test_cre', type=str, default='test_cre_state_prop.npz')
    parser.add_argument('--train_tss', action='store', dest='train_tss', type=str, default='trainTSS_window_state_prop.npz')
    parser.add_argument('--test_tss', action='store', dest='test_tss', type=str, default='testTSS_window_state_prop.npz')
    parser.add_argument('--train_exp', action='store', dest='train_exp', type=str, default='trainTPM.npz')
    parser.add_argument('--test_exp', action='store', dest='test_exp', type=str, default='testTPM.npz')
    parser.add_argument('--train_clustered', action='store', dest='train_clustered', type=str, default='{}_with_dists_20_-5.0_nj_clustered.npz')
    parser.add_argument('--chroms', action='store', nargs='+', dest='chroms', type=str, default='chr1', help='use all if you want to run all')

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
    args.train_clustered = argumentToAdd[where_run] + args.train_clustered

class regress_sampler():
    def __init__(self, train_cre, train_tss, train_exp):
        self.exp_values_all, self.cellIndex, self.cell_to_index, self.TSS_chr_all, self.TSSs_all = self.load_expression(train_exp)
        self.cellN = self.cellIndex.shape[0]
        self.TSS_window_props_all, self.TSS_window_chr_all = self.load_TSS_window_states(train_tss)
        self.cre_props_all, self.cre_chr_all, self.cre_coords_all = self.load_cre_states(train_cre)
        self.stateN  = self.cre_props_all.shape[2] #Note this value includes state 0 in the count although we're going to ignore the contribution of state 0

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

    def load_clustered(self, clustered_file):
        npz_file = np.load(clustered_file.format(self.chrom))
        clustered = npz_file['clustered']
        return clustered

    def subset_based_on_chrom(self, chrom, clustered_file):
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

        self.clustered = self.load_clustered(clustered_file)

    def compute_adj_distance(self, starts, stops, TSS):
        locsTA = (starts + stops) // 2
        distance = np.abs(locsTA- TSS) #using abs because don't care if upstream or downstream
        adjusted_distance = distance**self.gamma
        adjusted_distance[np.where(adjusted_distance == 0)[0]] = 1 #anywhere that distance is 0, set to identity.
        return adjusted_distance

    def adjust_by_distance(self, to_adjust, TSS, starts, stops):
        adj_dist = self.compute_adj_distance(starts, stops, TSS)
        Y = np.ones(to_adjust.shape[0]) #adjustment array
        Y /= adj_dist
        if to_adjust.ndim == 2:
            adjusted = np.multiply(to_adjust, Y.reshape(-1,1))
        elif to_adjust.ndim == 3:
            adjusted = np.multiply(to_adjust, Y.reshape(-1,1,1))
        return adjusted

    def take_logcosh_loss(self, true, yhat):
        '''rewrote to implement without keras and tf
        keras documentation says: logcosh = log((exp(x) + exp(-x))/2) where x is yhat - true'''
        x = yhat - true
        loss = np.sum(np.log((np.exp(x)+np.exp(-x))/2))
        return loss

    def indicator_function(self, tss_i):
        return self.clustered[tss_i]

    def get_beta_e(self):
        beta_e = np.hstack(([0], self.stacked_beta[self.stateN-1:])).reshape((1, -1))
        return beta_e

    def get_beta_p(self):
        beta_p = np.hstack(([0], self.stacked_beta[0:self.stateN-1])).reshape((1,-1))
        return beta_p

    def justP(self, tss_i):
        y_hat = np.sum(self.TSS_window_props[tss_i] * self.get_beta_p(), axis=1)
        if np.sum(np.isnan(y_hat)) > 0:
            print('justP has NaNs', flush=True)
            quit()
        return y_hat

    def justPair(self, tss_i, indicator_boolean):
        subset_weighted = np.sum(self.cre_props * indicator_boolean.reshape((-1,1,1))*self.get_beta_e(), axis=2)
        subset_adjusted = self.adjust_by_distance(subset_weighted, self.TSSs[tss_i], self.cre_coords[:, 0], self.cre_coords[:,1])
        sum_all_E = np.sum(subset_adjusted, axis=0)
        if np.sum(np.isnan(sum_all_E)) > 0:
            print('justPair has NaNs', flush=True)
            quit()
        return sum_all_E

    def regression_equation(self, tss_i):
        PairingFlag = True
        indicator_boolean = self.indicator_function(tss_i)
        if np.sum(indicator_boolean) == 0:
            PairingFlag = False
            return self.justP(tss_i), PairingFlag
        yhat = self.justP(tss_i) + self.justPair(tss_i, indicator_boolean)

    def update_yhats(self):
        yhats = np.zeros((self.tssN, self.cellN), dtype=np.float32)
        for i in range(self.tssN):
            yhats[i], self.flags[i] = self.regression_equation(i)
        self.logcosh_expression = self.take_logcosh_loss(self.exp_values, yhats)

    def stack_X_data(self):
        stacked_XE_data = np.zeros((self.tssN, self.cellN, self.stateN-1))
        for i in range(self.tssN):
            indicator_boolean = self.indicator_function(i)
            if np.sum(indicator_boolean) == 0:
                continue
            stacked_XE_data[i] = np.sum(self.adjust_by_distance(self.cre_props * indicator_boolean.reshape((-1,1,1)), self.TSSs[i], self.cre_coords[:, 0], self.cre_coords[:, 1]), axis=0)[:,1:]
        stacked_data = np.hstack((self.TSS_window_props.reshape((-1, self.stateN))[:,1:], stacked_XE_data.reshape((-1,self.stateN-1))))
        return stacked_data

    def linear_fit(self, X, Y, intercept=False):
        Y = Y.reshape(-1,)
        fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X, Y)
        return fit_model.coef_

    def report_iteration_hyperparameters(self, stacked_beta):
        toWriteTo_beta = open('output_beta_hyperparameters_{}_g{}_clustering_20_-5.txt'.format(self.chrom, self.gamma), 'a')
        np.savetxt(toWriteTo_beta, stacked_beta)
        toWriteTo_beta.close()

    def report_metrics(self, logcosh_sum, numNP):
        toWriteTo = open('output_metrics_clustering_20_-5.txt', 'a')
        toWriteTo.write('chrom:\t{}\tgamma:\t{}\tLoss:\t{}\tnumNP:\t{}\tnumNPRatio:\t{}\n'.format(self.chrom, self.gamma, logcosh_sum, numNP, numNP/self.tssN))
        toWriteTo.close()

    def run_fit(self):
        gamma_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5]
        for gamma_val in gamma_values:
            self.gamma = gamma_val
            self.flags = np.zeros(self.tssN, dtype=np.bool)
            stacked_data = self.stack_X_data()
            self.stacked_beta = self.linear_fit(stacked_data, self.exp_values)
            self.update_yhats()
            self.report_metrics(self.logcosh_expression, np.sum(~self.flags))
            self.report_iteration_hyperparameters(self.stacked_beta)


main()
