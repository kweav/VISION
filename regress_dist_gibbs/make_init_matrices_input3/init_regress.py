#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model, metrics
import os
from scipy import stats, linalg
import logging
import datetime


'''outstanding questions/todo:
1) Optimize - where is the time being spent?
One iteration appears to take ~0.0399 minutes
2) How are NaN's in the correlation matrix being propagated?
- Comparing a NaN in greater than or equal returns a False which is what I want, but the warnings are still passed
- Should I suppress/catch the warnings?
3) Since I'm using a normal for k,
- abs value of whatever is returned
- genereate from uniform or beta(5,2) if sampled value > 1
4) Should I update correlations at each iteration or use the initial ones?
- Spearmanr is the longest part of the initial setup ~185sec, so take that as a last resort'''

def main():
    #logging.info('began main(): ' + str(datetime.datetime.now()))
    parser = generate_parser()
    #logging.info('generate_parser() complete: ' + str(datetime.datetime.now()))
    args = parser.parse_args()
    #logging.info('parse_args() complete: ' + str(datetime.datetime.now()))
    setup_file_locs(args, args.where_run, args.other_path)
    #logging.info('setup_file_locs() complete: ' + str(datetime.datetime.now()))
    setup_threads(args.threads)
    #logging.info('setup_threads() complete: ' + str(datetime.datetime.now()))
    model = regress_sampler(args.train_tss, args.train_exp)
    #logging.info('regress_sampler model initiated: ' + str(datetime.datetime.now()))
    return model.find_initial_coeffs()


def generate_parser():
    parser = ap.ArgumentParser(description ='initial regression to get beta mu')
    parser.add_argument('--where', action='store', dest='where_run', type=str, default='mine', help='{comp, mine, marcc, other}; adds path to files; if "other" used provide the other path in the --otherpath argument')
    parser.add_argument('--threads', action='store', dest='threads', type=str, default="1")
    parser.add_argument('--otherpath', action='store', dest='other_path', type=str, default='NA', help='use to give other path if not one of the 3 previously specified. MUST be the same for ALL input files')
    parser.add_argument('--train_tss', action='store', dest='train_tss', type=str, default='trainTSS_window_state_prop.npz')
    parser.add_argument('--test_tss', action='store', dest='test_tss', type=str, default='testTSS_window_state_prop.npz')
    parser.add_argument('--train_exp', action='store', dest='train_exp', type=str, default='trainTPM.npz')
    parser.add_argument('--test_exp', action='store', dest='test_exp', type=str, default='testTPM.npz')
    parser.add_argument('--chroms', action='store', nargs='+', dest='chroms', type=str, default='all')

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
                     'comp': '/project/vision/Target_Genes/simple_sampler/npz_inputs/',
                     'other': other_path}

    args.train_exp = argumentToAdd[where_run] + args.train_exp
    args.test_exp = argumentToAdd[where_run] + args.test_exp
    args.train_tss = argumentToAdd[where_run] + args.train_tss
    args.test_tss = argumentToAdd[where_run] + args.test_tss



class regress_sampler():
    def __init__(self, train_tss, train_exp):
        self.exp_values, self.cellIndex, self.cell_to_index, self.TSS_chr, self.TSSs = self.load_expression(train_exp)
        self.cellN = self.cellIndex.shape[0]
        self.TSS_window_props, self.TSS_window_chr = self.load_TSS_window_states(train_tss)
        self.stateN  = self.TSS_window_props.shape[2] #Note this value includes state 0 in the count although we're going to ignore the contribution of state 0
        print(self.stateN, flush=True)

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
        #logging.info('begin find_valid_cellTypes(): ' + str(datetime.datetime.now()))
        valid = np.array([np.where(cellIndexOI == x) for x in self.cellIndex if np.isin(x, cellIndexOI)]).reshape(-1)
        #logging.info('end find_valid_cellTypes(): ' + str(datetime.datetime.now()))
        return valid

    def load_expression(self, exp_file):
        #logging.info('begin load_expression(): ' + str(datetime.datetime.now()))
        npzfile = np.load(exp_file, allow_pickle = True)
        exp_values = npzfile['exp']
        exp_cellIndex = npzfile['cellIndex'] #index to celltype
        exp_cell_to_index = {} #celltype to index
        for i in range(exp_cellIndex.shape[0]):
            exp_cell_to_index[exp_cellIndex[i]] = i
        TSS_info = npzfile['TSS']
        TSS_chr = TSS_info[:,0]
        TSSs = np.around(TSS_info[:,1].astype(np.float64)).astype(np.int32)
        #logging.info('end load_expression(): ' + str(datetime.datetime.now()))
        return exp_values, exp_cellIndex, exp_cell_to_index, TSS_chr, TSSs

    def load_TSS_window_states(self, tss_state_file):
        #logging.info('begin load_TSS_window_states(): ' + str(datetime.datetime.now()))
        npzfile = np.load(tss_state_file, allow_pickle = True)
        TSS_window_props = npzfile['props'].astype(np.float64)
        TSS_cellIndex = npzfile['cellIndex']
        valid = self.find_valid_cellTypes(TSS_cellIndex)
        TSS_window_props_valid = TSS_window_props[:,valid,:]
        TSS_window_index = npzfile['ccREIndex']
        TSS_window_chr = TSS_window_index[:,0]
        #TSS_window_coord = TSS_window_index[:,1:].astype(np.int32)
        #logging.info('end load_TSS_window_states(): ' + str(datetime.datetime.now()))
        return TSS_window_props_valid, TSS_window_chr#, TSS_window_coord


    def linear_regression(self, X, Y, intercept=False):
        '''reshape arrays to match function input requirements
        remove contribution of state 0 from X; will set this state's contribution to zero afterwards'''
        #logging.info('begin linear_regression(): ' + str(datetime.datetime.now()))
        X = X.reshape(-1, self.stateN)[:,1:]
        Y = Y.reshape(-1,)
        fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X,Y)
        model_coeffs = fit_model.coef_
        r_squared = fit_model.score(X,Y)
        #logging.info('end linear_regression(): ' + str(datetime.datetime.now()))
        return {'coeffs': model_coeffs, 'rsquare': r_squared, 'fit_model': fit_model}

    def find_initial_coeffs(self):
        '''do initial regression of proportion of states within two-sided 75kbp window against expression to find initial coefficients for all states except 0
        set contribution of state 0 to be 0'''
        initial_coeffs = self.linear_regression(self.TSS_window_props, self.exp_values)['coeffs']
        return (initial_coeffs)

if __name__ == "__main__":
    main()
