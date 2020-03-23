#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model
import os
from scipy import stats, linalg

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
        model.subset_based_on_chrom(chrom)
        model.run_sampler(args.init_beta, args.init_gamma, args.init_k, args.init_sigma_sqr, args.iters, args.burn_in, args.cre_dist)

def generate_parser():
    parser = ap.ArgumentParser(description='VISION Gibbs Sampler of parameters for regression of state and distance to assign CREs to genes based on the ability to predict gene expression')

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
    parser.add_argument('--burn_in', action='store', dest='burn_in', type=int, default=500)
    parser.add_argument('--iterations', action='store', dest='iters', type=int, default=10**5)
    parser.add_argument('--init_beta', action='store', dest='init_beta', type=str, default='init_beta.npy')
    parser.add_argument('--init_gamma', action='store', dest='init_gamma', type=float, default=0.7)
    parser.add_argument('--init_k', action='store', dest='init_k', type=float, default=0.2)
    parser.add_argument('--init_sigma_sqr', action='store', dest='init_sigma_sqr', type=float, default=1.0)
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
    args.init_beta = argumentToAdd[where_run] + args.init_beta

class regress_sampler():
    def __init__(self, train_cre, train_tss, train_exp):
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

        def find_cres_within(self, i):
            '''#find CREs within distance of interest using containment
            start location of CRE is less than or equal to window end
            AND end locatiion of CRE is greater than or equal to window beginning'''
            TSS = self.TSSs[i]
            windowMin = max(0, TSS - self.cre_dist)
            windowMax = min(TSS + self.cre_dist, self.chrSizes[self.chrom])
            CREs_within = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
            return CREs_within

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

        def find_weighted_sum(self, beta_e, CREs_within_bool):
            '''find weighted sum of ccRE props with coeffs
            given a 3 dimensional array ccREn * cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                    1 dimensional array stateN with c_i equal to the initial coefficient for state_i
            return a 2 dimensional array ccREn * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''
            cre_weighted_sum = np.sum(self.cre_props[CREs_within_bool] * beta_e, axis=2)
            return cre_weighted_sum

        def take_inv(self, x):
            '''approximate inverses with pinv if matrix singular'''
            try:
                inv = linalg.inv(x)
            except linalg.LinAlgError:
                inv = linalg.pinv(x)
            return inv

        def take_logcosh_loss(self, true, yhat):
            '''rewrote to implement without keras and tf
            keras documentation says: logcosh = log((exp(x) + exp(-x))/2) where x is yhat - true'''
            x = yhat - true
            loss = np.sum(np.log((np.exp(x)+np.exp(-x))/2))
            return loss

        def posterior_k(self):
            '''prior of beta(3,1)
            use inverse transform method to generate from the approximate posterior
            use logcosh loss rather than S^2/RSS/MSE as to be less sensitive to outliers'''
            generated_u = stats.uniform.rvs()
            sampled = np.sqrt(generated_u/3 * np.exp(1/(2*self.sigma_sqr)*self.logcosh_expression))
            return sampled

        def posterior_gamma(self):
            '''prior of normal(0.7, 0.075)
            use inverse transform method to generate from the approximate posterior
            use logcosh loss rather than S^2/RSS/MSE as to be less senstive to outliers'''
            generated_u = stats.uniform.rvs()
            sampled = 0.7 + np.sqrt((-6.67)*(np.log(generated_u/2*np.sqrt(0.6*np.pi))+(1/(2*self.sigma_sqr)*self.logcosh_expression)))
            return sampled

        def posterior_sigma_sqr(self):
            '''prior of gamma(1,1)
            this setup generates the precision which is why we return 1/precision
            use inverse transform method to generate from the approximate posterior
            use logcosh loss rather than S^2/RSS/MSE as to be less sensitive to outliers'''
            generated_u = stats.uniform.rvs()
            precision = -np.log(generated_u)/(0.5*self.logcosh_expression + 1)
            return 1/precision

        def compute_spearmanr(self, tss_i, cre_weighted_sum):
            corr_matrix, pvalues = stats.spearmanr(self.exp_values[tss_i].reshape((1,-1)), cre_weighted_sum, axis=1)
            if isinstance(corr_matrix, np.ndarray):
                corr = corr_matrix[0, 1:]
            else:
                corr = corr_matrix
            return corr

        def indicator_function(self, tss_i, CREs_within_bool):
            weighted_sum = self.find_weighted_sum(self.get_beta_e(), CREs_within_bool)
            spearman = self.compute_spearmanr(tss_i, weighted_sum)
            indicator_boolean = np.abs(spearman) >= self.k
            return indicator_boolean

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

        def justPair(self, tss_i, CREs_within, indicator_boolean):
            subset_weighted = np.sum(self.cre_props[CREs_within] * indicator_boolean.reshape((-1,1,1))*self.get_beta_e(), axis=2)
            subset_adjusted = self.adjust_by_distance(subset_weighted, self.TSSs[tss_i], self.cre_coords[CREs_within, 0], self.cre_coords[CREs_within, 1])
            sum_all_E = np.sum(subset_adjusted, axis=0)
            if np.sum(np.isnan(sum_all_E)) > 0:
                print('justPair has NaNs', flush=True)
                quit()
            return sum_all_E

        def regression_equation(self, tss_i):
            PairingFlag = True
            CREs_within = self.find_cres_within(tss_i)
            if np.sum(CREs_within) == 0:
                PairingFlag = False
                return self.justP(tss_i), PairingFlag
            indicator_boolean = self.indicator_function(tss_i, CREs_within)
            if np.sum(indicator_boolean) == 0:
                PairingFlag = False
                return self.justP(tss_i), PairingFlag
            yhat = self.justP(tss_i) + self.justPair(tss_i, CREs_within, indicator_boolean)
            return yhat, PairingFlag

        def update_yhats(self):
            yhats = np.zeros((self.tssN, self.cellN), dtype=np.float32)
            for i in range(self.tssN):
                yhats[i], self.flags[i] = self.regression_equation(i)
            self.logcosh_expression = self.take_logcosh_loss(self.exp_values, yhats)

        def stack_X_data(self):
            stacked_XE_data = np.zeros((self.tssN, self.cellN, self.stateN-1))
            for i in range(self.tssN):
                CREs_within_bool = self.find_cres_within(tss_i)
                if np.sum(CREs_within_bool) == 0:
                    continue
                indicator_boolean = self.indicator_function(tss_i, CREs_within_bool)
                if np.sum(indicator_boolean) == 0:
                    continue
                stacked_XE_data[i] = np.sum(self.adjust_by_distance(self.cre_props[CREs_within_bool] * indicator_boolean.reshape((-1,1,1)), self.TSSs[i], self.cre_coords[CREs_within_bool, 0], self.cre_coords[CREs_within_bool, 1]), axis=0)[:,1:]
            stacked_data = np.hstack((self.TSS_window_props.reshape((-1, self.stateN))[:,1:], stacked_XE_data.reshape((-1,self.stateN-1))))
            return stacked_data

        def linear_fit(self, X, Y, intercept=False):
            Y = Y.reshape(-1,)
            fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X, Y)
            return fit_model.coef_

        def update_parameters(self):
            self.sigma_sqr = self.posterior_sigma_sqr()
            self.update_yhats()
            self.k = self.posterior_k()
            self.update_yhats()
            self.gamma = self.posterior_gamma()
            stacked_data = self.stack_X_data()
            self.stacked_beta = self.linear_fit(stacked_data, self.exp_values)

        def report_iteration_hyperparameters(self, iteration, param_dict):
            toWriteTo_scalar = open('output_scalar_hyperparameters.txt', 'a')
            toWriteTo_scalar.write('Iteration:\t{}\tsigma_sqr:\t{}\tk:\t{}\tgamma:\t{}\n'.format(iteration, param_dict['sigma_sqr'], param_dict['k'], param_dict['gamma']))
            toWriteTo_scalar.close()

            toWriteTo_beta = open('output_beta_hyperparameters.txt', 'a')
            toWriteTo_beta.write('Iteration:\t{}\n'.format(iteration))
            np.savetxt(toWriteTo_beta, param_dict['stacked_beta'], fmt='%.5f')
            toWriteTo_beta.close()

        def report_metrics(self, iteration, logcosh_sum, numNP):
            toWriteTo = open('output_metrics.txt', 'a')
            toWriteTo.write('Iteration:\t{}\tLoss:\t{}\tnumNP:\t{}\tnumNPRatio:\t{}\n'.format(iteration, logcosh_sum, numNP, numNP/self.tssN))
            toWriteTo.close()

        def run_sampler(self, init_beta, init_gamma, init_k, init_sigma_sqr, iters, burn_in, cre_dist):
            self.cre_dist = cre_dist
            self.stacked_beta, self.gamma, self.k, self.sigma_sqr = np.load(init_beta), init_gamma, init_k, init_sigma_sqr
            argmin = {'stacked_beta': np.copy(self.stacked_beta),
                      'gamma': self.gamma,
                      'k': self.k,
                      'sigma_sqr': self.sigma_sqr}

            self.flags = np.zeros(self.tssN, dtype=np.bool)

            self.update_yhats()
            minLoss = self.logcosh_expression
            minNP = np.sum(~self.flags)
            self.report_metrics(-1, minLoss, minNP)
            self.report_iteration_hyperparameters(-1, argmin)

            for iteration in range(iters):
                self.update_parameters()
                self.update_yhats()
                if (iteration > burn_in) and (self.logcosh_expression < minLoss):
                    argmin['stacked_beta'] = np.copy(self.stacked_beta)
                    argmin['gamma'] = self.gamma
                    argmin['k'] = self.k
                    argmin['sigma_sqr'] = self.sigma_sqr
                    minLoss = self.logcosh_expression
                    minNP = np.sum(~self.flags)
                    self.report_metrics(iteration, minLoss, minNP)
                    self.report_iteration_hyperparameters(iteration, argmin)
                else:
                    iter_dict = {'stacked_beta': np.copy(self.stacked_beta),
                                 'gamma': self.gamma,
                                 'k': self.k,
                                 'sigma_sqr': self.sigma_sqr}
                    self.report_metrics(iteration, self.logcosh_expression, np.sum(~self.flags))
                    self.report_iteration_hyperparameters(iteration, iter_dict)

            self.report_metrics('argmin', minLoss, minNP)
            self.report_iteration_hyperparameters('argmin', argmin)

main()
