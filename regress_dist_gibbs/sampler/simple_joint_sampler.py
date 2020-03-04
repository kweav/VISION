#!/usr/bin/env python3

import argparse as ap
import numpy as np
from sklearn import linear_model, metrics
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
        for i in range(1,20):
            chrom_list.append('chr{}'.format(i))
        chrom_list.append('chrX')
    else:
        chrom_list = args.chroms

    for chrom in chrom_list:
        model.subset_based_on_chrom(chrom)
        model.find_initial_weighted_sum()
        model.compute_spearmanr()
        model.set_up_prior_info(args.sigma_gamma_alpha, args.sigma_gamma_beta, args.gamma_norm_mu, args.gamma_norm_var, args.k_norm_mu, args.k_norm_var, args.Sigma_invwishart_v_0, args.Sigma_invwishart_S_0, args.theta_MVN_Lambda_0, args.theta_MVN_mu_0, args.bsfs)
        model.run_sampler(args.init_beta, args.init_theta, args.init_Sigma, args.init_gamma, args.init_k, args.init_sigma_sqr, args.iters, args.burn_in)

def generate_parser():
    parser = ap.ArgumentParser(description ='VISION Gibbs Sampler of paramters for regression of state and distance to assign CREs to genes based on ability to predict gene expression')
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
    parser.add_argument('--burn_in', action='store', dest='burn_in', type=int, default=250)
    parser.add_argument('--iterations', action='store', dest='iters', type=int, default=10**4)
    parser.add_argument('--init_beta', action='store', dest='init_beta', type=str, default='init_beta.npy')
    parser.add_argument('--init_theta', action='store', dest='init_theta', type=str, default='init_theta.npy')
    parser.add_argument('--init_Sigma', action='store', dest='init_Sigma', type=str, default='init_Sigma.npy')
    parser.add_argument('--init_gamma', action='store', dest='init_gamma', type=float, default=0.7)
    parser.add_argument('--init_k', action='store', dest='init_k', type=float, default=0.4)
    parser.add_argument('--init_sigma_sqr', action='store', dest='init_sigma_sqr', type=float, default=1.0)
    parser.add_argument('--sigma_gamma_alpha', action='store', dest='sigma_gamma_alpha', type=float, default=1.0)
    parser.add_argument('--sigma_gamma_beta', action='store', dest='sigma_gamma_beta', type=float, default=1.0)
    parser.add_argument('--gamma_norm_mu', action='store', dest='gamma_norm_mu', type=float, default=0.7)
    parser.add_argument('--gamma_norm_var', action='store', dest='gamma_norm_var', type=float, default=0.075)
    parser.add_argument('--k_norm_mu', action='store', dest='k_norm_mu', type=float, default=0.5)
    parser.add_argument('--k_norm_var', action='store', dest='k_norm_var', type=float, default=0.175)
    parser.add_argument('--Sigma_invwishart_v_0', action='store', dest='Sigma_invwishart_v_0', type=int, default=54)
    parser.add_argument('--Sigma_invwishart_S_0', action='store', dest='Sigma_invwishart_S_0', type=str, default='Sigma_invwishart_S_0.npy')
    parser.add_argument('--theta_MVN_Lambda_0', action='store', dest='theta_MVN_Lambda_0', type=str, default='theta_MVN_Lambda_0.npy')
    parser.add_argument('--theta_MVN_mu_0', action='store', dest='theta_MVN_mu_0', type=str, default='theta_MVN_mu_0.npy')
    parser.add_argument('--beta_shape_forSamp', action='store', nargs='+', dest='bsfs', type=int, default=52)
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

    args.train_cre = argumentToAdd[where_run] + args.train_cre
    args.test_cre = argumentToAdd[where_run] + args.test_cre
    args.train_exp = argumentToAdd[where_run] + args.train_exp
    args.test_exp = argumentToAdd[where_run] + args.test_exp
    args.train_tss = argumentToAdd[where_run] + args.train_tss
    args.test_tss = argumentToAdd[where_run] + args.test_tss

class regress_sampler():
    def __init__(self, train_cre, train_exp, train_tss):
        self.exp_values_all, self.cellIndex, self.cell_to_index, self.TSS_chr_all, self.TSSs_all = self.load_expression(train_exp)
        self.cellN = self.cellIndex.shape[0]
        self.TSS_window_props_all, self.TSS_window_chr_all, self.TSS_window_coord_all = self.load_TSS_window_states(train_tss)
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
        self.tssN = self.TSSs.shape[0]
        self.TSS_window_props = self.TSS_window_props_all[where_chr]
        self.TSS_window_coord = self.TSS_window_coord_all[where_chr]

        where_chr = np.where(self.cre_chr_all == self.chrom)[0]
        self.cre_props = self.cre_props_all[where_chr]
        self.cre_coords = self.cre_coords_all[where_chr]
        self.creM = self.cre_props.shape[0]
        self.creIndex_range = np.arange(self.creM)

    def linear_regression(self, X, Y, intercept=False):
        '''reshape arrays to match function input requirements
        remove contribution of state 0 from X; will set this state's contribution to zero afterwards'''
        X = X.reshape(-1, self.stateN)[:,1:]
        Y = Y.reshape(-1,)
        fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X,Y)
        model_coeffs = fit_model.coef_
        r_squared = fit_model.score(X,Y)
        return {'coeffs': model_coeffs, 'rsquare': r_squared, 'fit_model': fit_model}

    def find_initial_weighted_sum(self):
        '''do initial regression of proportion of states within two-sided 75kbp window against expression to find initial coefficients for all states except 0
        set contribution of state 0 to be 0'''
        initial_coeffs = self.linear_regression(self.TSS_window_props, self.exp_values)['coeffs']
        initial_coeffs_f = np.hstack(([0], initial_coeffs)).reshape((1,-1))
        '''find weighted sum of ccRE props with initial_coeffs
        given a 3 dimensional array ccREn * cellN * stateN with p_ijk equal to the proporition of ccRE_i for cellType_j in state_k
                1 dimensional array stateN with c_i equal to the initial coefficient for state_i
        return a 2 dimensional array ccREn * cellN with w_ij equal to the sum of l = 0 to 26 of c_l*p_ijl'''
        self.cre_weighted_sum = np.sum(self.cre_props * initial_coeffs_f, axis=2)

    def compute_spearmanr(self): #note we are not adjusting for distance when looking at this correlation
        self.corr_matrix, pvalues = stats.spearmanr(self.exp_values, self.cre_weighted_sum, axis=1)

    def find_cres_within(self, i):
        TSS = self.TSSs[i]
        #find CREs within distance of interest using containment
        windowMin = max(0, TSS - self.cre_dist)
        windowMax = min(TSS + self.cre_dist, self.chrSizes[self.chrom])
        CREs_within = (self.cre_coords[:,1]>=windowMin) & (self.cre_coords[:,0]<=windowMax)
        return CREs_within

    def indicator_function(self, tss_i, CREs_within_bool):
        return np.abs(self.corr_matrix[tss_i, self.tssN:][CREs_within_bool]) >= self.k #take slice of correlations between that TSS and the CREs within the window

    def no_pair_justP(self, tss_i):
        beta_p = np.hstack(([0],self.stacked_beta[0:self.stateN-1])).reshape((1,-1))
        y_hat = np.sum(self.TSS_window_props[tss_i] * beta_p, axis=2)
        return y_hat

    def compute_adj_distance(self, starts, stops, TSS):
        locsTA = (starts + stops) // 2
        distance = np.abs(locsTA- TSS) #using abs because don't care if upstream or downstream
        adjusted_distance = distance**self.gamma
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

    def pair_noP(self, tss_i, CREs_within, indicator_boolean, firstTimeCalled=False, withinFirstSet=False):
        beta_e = np.hstack(([0],self.stacked_beta[self.stateN-1:])).reshape((1,-1))
        #weighted_sum including indicator
        if firstTimeCalled:
            self.build_X_e = np.zeros((self.tssN, self.cellN, self.stateN-1))
        if withinFirstSet:
            CREs_within_starts = self.cre_coords[CREs_within, 0]
            CREs_within_stops = self.cre_coords[CREs_within, 1]
            self.build_X_e[tss_i] = np.sum(self.adjust_by_distance(self.cre_props[CREs_within] * inidicator_boolean, self.TSSs[tss_i], CREs_within_starts, CREs_within_stops), axis=2)[:,1:]
        subset_weighted = np.sum(self.cre_props[CREs_within] * inidicator_boolean * beta_e, axis=2)
        #adjust by distance
        CREs_within_starts = self.cre_coords[CREs_within, 0]
        CREs_within_stops = self.cre_coords[CREs_within, 1]
        subset_adjusted = self.adjust_by_distance(subset_weighted, self.TSSs[tss_i], CREs_within_starts, CREs_within_stops)
        return np.sum(subset_adjusted)

    def regression_equation(self, tss_i, firstTimeCalled=False, withinFirstSet=False):
        PairingFlag = True
        CREs_within = self.find_cres_within(tss_i)
        if np.sum(CREs_within) == 0: #no pairing possible
            PairingFlag = False
            return self.no_pair_justP(tss_i), PairingFlag
        inidicator_boolean = self.indicator_function(tss_i, CREs_within)
        if np.sum(indicator_boolean) == 0: #no pairing happened
            PairingFlag = False
            return self.no_pair_justP(tss_i), PairingFlag
        yhat = self.no_pair_justP(tss_i) + self.pair_noP(tss_i, CREs_within, indicator_boolean, firstTimeCalled, withinFirstSet)
        print(yhat.shape)
        return yhat, PairingFlag

    def find_MSE(self, i, yhat):
        return metrics.mean_squared_error(self.exp_values[i], yhat)

    def run_regression_equation(self, initialTime=False):
        #find yhat, note if pairing was possible, and find MSE
        totalMSE = 0
        flags = np.zeros(self.tssN, dtype=np.bool)
        yhats = np.zeros((self.tssN, self.cellN), dtype=np.float32)
        for i in range(tssN):
            if initialTime and i == 0:
                yhat, PairingFlag = self.regression_equation(i, firstTimeCalled=True, withinFirstSet=True)
            elif initialTime and i != 0:
                yhat, PairingFlag = self.regression_equation(i, withinFirstSet=True)
            else:
                yhat, PairingFlag = self.regression_equation(i)
            flags[i] = PairingFlag
            yhats[i] = yhat
            totalMSE += self.find_MSE(i, yhat)
        return yhats, totalMSE, np.sum(flags)

    def update_yhats(self):
        for i in range(tssN):
            self.yhats[i], PairingFlag = self.regression_equation(i)

    # def sample_MNV(self, mu, Sigma):
    #     lower_cholesky = linalg.cholesky(Sigma, lower=True)
    #     X = scipy.norm.rvs(size=self.stacked_beta_n)
    #     Z = mu + lower_cholesky*X
    #     return Z

    def posterior_gamma(self, sigma_sqr, mu_data):
        norm_mu = (((self.sigma_norm_mu/self.sigma_norm_var)+(np.sum(mu_data)/sigma_sqr))/((1/self.gamma_norm_var) + (self.tssN*self.cellN/sigma_sqr)))
        norm_var = 1/((1/self.sigma_norm_var) + (self.tssN*self.cellN/sigma_sqr))
        sampled = stats.norm.rvs(loc=norm_mu, scale=np.sqrt(norm_var))
        return sampled

    def posterior_k(self, sigma_sqr, mu_data):
        norm_mu = (((self.k_norm_mu/self.k_norm_var)+(np.sum(mu_data)/sigma_sqr))/((1/self.k_norm_var) + (self.tssN*self.cellN/sigma_sqr)))
        norm_var = 1/((1/self.k_norm_var) + (self.tssN*self.cellN/sigma_sqr))
        sampled = stats.norm.rvs(loc=norm_mu, scale=np.sqrt(norm_var))
        return sampled

    def posterior_theta(self, Beta, Sigma):
        '''process taken from A First Course in Bayesian Statistical Methods - Peter D. Hoff ISBN: 978-0-387-92299-7
           Chapter 7 The multivariate normal model 7.2 A semiconjugate prior distribution for the mean and 7.4 Gibbs sampling of the mean and covariance'''
        #to Sample Theta^(s+1)
        lambda_0_inv = linalg.inv(self.theta_MVN_Lambda_0)
        Sigma_inv = linalg.inv(Sigma)
        #a - compute mu_n and Lambda_n from beta and Sigma^(s)
        #note: self.stacked_beta_n is a constant
        #      lambda_0_inv is a matrix
        #      Sigma_inv is a matrix
        #      Beta is a vector
        #      self.theta_MVN_mu_0 is a vector
        #      Lambda_n is a matrix
        Lambda_n = linalg.inv(lambda_0_inv + self.stacked_beta_n*Sigma_inv)
        mu_n = Lambda_n.dot(lambda_0_inv.dot(self.theta_MVN_mu_0) + self.stacked_beta_n*Sigma_inv.dot(Beta))
        #above here the Sigma_inv.dot(Beta) should be the dot product of the Sigma_inv and the vector of variable_specific averages. Since we assume that Beta is the same for all cell types, this average is just Beta.
        #b - sample Theta^(s+1) ~ MVN(mu_n, Lambda_n)
        sampled = stats.multivariate_normal.rvs(mean = mu_n, cov=Lambda_n)
        return sampled

    def posterior_Sigma(self, Beta, theta):
        '''process taken from A First Course in Bayesian Statistical Methods - Peter D. Hoff ISBN: 978-0-387-92299-7
           Chapter 7 The multivariate normal model 7.3 The inverse-Wishart distribution and 7.4 Gibbs sampling of the mean and covariance
           Note that within this process I am assuming that the true covariance is only loosely centered around our initial covariance input by saying Sigma_0 = S_0 and v_0 equals p+2'''
        #to Sample Sigma^(s+1)
        #a - compute S_n from Beta and Theta^(s+1)
        S_theta = np.sum(np.square(Beta-theta))
        S_n = self.Sigma_invwishart_S_0 + S_theta
        #b - sample Sigma^(s+1) ~ inverse-Wishart(v_0 +n, S_n^-1)
        df = self.Sigma_invwishart_v_0 + self.stacked_beta_n
        scale = linalg.inv(S_n)
        sampled = stats.invwishart.rvs(df=df, scale=scale)
        return sampled

    def posterior_beta(self, sigma_sqr, stacked_X_data, Sigma, theta, mu_data):
        '''process taken from A First Course in Bayesian Statistical Methods - Peter D. Hoff ISBN: 978-0-387-92299-7
           Chapter 9 Linear Regression 9.2 Bayesian estimation for a regression model 9.2.1 A semiconugate prior distribution
           However, the process interchanges the most recent iteration of theta (from the MVN prior) and Sigma (from the inv-wishart) rather than beta_0 and Sigma_0 '''
        Sigma_inv = linalg.inv(Sigma)
        variance = linalg.inv(Sigma_inv + (stacked_X_data.T).dot(stacked_X_data)/sigma_sqr)
        expectation_inner_dot = Sigma_inv.dot(theta) + stacked_X_data.T.dot(mu_data.reshape(-1,))/sigma_sqr
        expectation = variance.dot(expectation_inner_dot)
        sampled = stats.multivariate_normal(mean=expectation, cov=variance)
        return sampled

    def posterior_sigma_sqr(self, u):
        gamma_alpha = 1 + (self.tssN*self.cellN/2)
        gamma_beta = 1 + ((1/2) * np.sum(self.exp_values-u))
        precision = stats.gamma.rvs(x, a=gamma_alpha, scale=1/gamma_beta)
        return 1/precision

    def set_up_prior_info(self, sigma_gamma_alpha, sigma_gamma_beta, gamma_norm_mu, gamma_norm_var, k_norm_mu, k_norm_var, Sigma_invwishart_v_0, Sigma_invwishart_S_0, theta_MVN_Lambda_0, theta_MVN_mu_0, bsfs):
        self.sigma_gamma_alpha, self.sigma_gamma_beta = sigma_gamma_alpha, sigma_gamma_beta,
        self.gamma_norm_mu, self.gamma_norm_var = gamma_norm_mu, gamma_norm_var
        self.k_norm_mu, self.k_norm_var = k_norm_mu, k_norm_var
        self.Sigma_invwishart_v_0, self.Sigma_invwishart_S_0 = Sigma_invwishart_v_0, np.load(Sigma_invwishart_S_0)
        self.theta_MVN_Lambda_0, self.theta_MVN_mu_0 = np.load(theta_MVN_Lambda_0), np.load(theta_MVN_mu_0)
        self.stacked_beta_n = bsfs

    def report_iteration_hyperparameters(self, iteration):
        toWriteTo_scalar = open('output_scalar_hyperparameters.txt', 'a')
        toWriteTo_scalar.write('Iteration:\t{}\tsigma_sqr:\t{}\tk:\t{}\tgamma:\t{}\n'.format(iteration, self.sigma_sqr, self.k, self.gamma))
        toWriteTo_scalar.close()

        toWriteTo_Sigma = open('output_Sigma_hyperparameters.txt', 'ab')
        toWriteTo_Sigma.write('Iteration:\t{}\n'.format(iteration))
        np.savetxt(toWriteTo_Sigma, self.Sigma, fmt='%.5f')
        toWriteTo_Sigma.close()

        toWriteTo_theta = open('output_theta_hyperparameters.txt', 'ab')
        toWriteTo_theta.write('Iteration:\t{}\n'.format(iteration))
        np.savetxt(toWriteTo_theta, self.theta, fmt='%.5f')
        toWriteTo_theta.close()

        toWriteTo_beta = open('output_beta_hyperparameters.txt', 'ab')
        toWriteTo_beta.write('Iteration:\t{}\n'.format(iteration))
        np.savetxt(toWriteTo_beta, self.stacked_beta, fmt='%.5f')
        toWriteTo_beta.close()

    def report_metrics(self, iteration, MSE_sum, numNP, minNotPaired):
        toWriteTo = open('output_metrics.txt', 'a')
        toWriteTo.write('Iteration:\t{}\tMSE_sum:\t{}\tnotPaired:\t{}\tnotPairedRatio:\t{}\n'.format(iteration, MSE_sum, numNP, (minNotPaired - numNP)/self.tssN ))


    def get_stacked_X_data(self):
        return np.hstack((self.TSS_window_props.reshape((-1, stateN))[:,1:], self.build_X_e.reshape((-1, stateN-1))))

    def update_parameters(self):
        self.sigma_sqr = self.posterior_sigma_sqr(self.yhats)
        self.update_yhats()
        self.k = self.posterior_k(self.sigma_sqr, self.yhats)
        self.update_yhats()
        self.gamma = self.posterior_gamma(self.sigma_sqr, self.yhats)
        self.update_yhats()
        self.theta = self.posterior_theta(self.stacked_beta, self.Sigma)
        #don't update_yhats because theta doesn't directly affect yhat
        self.Sigma = self.posterior_Sigma(self.stacked_beta,  self.theta)
        #don't update yhats because Sigma doesn't directly affect yhat
        self.stacked_beta = self.posterior_beta(self.sigma_sqr, self.get_stacked_X_data(), self.Sigma, self.theta, self.yhats)
        #self.update_yhats() unnecessary as no more parameters to update

    def run_sampler(self, init_beta, init_theta, init_Sigma, init_gamma, init_k, init_sigma_sqr, iters, burn_in):
        self.stacked_beta, self.theta, self.Sigma, self.gamma, self.k, self.sigma_sqr = np.load(init_beta), np.load(init_theta), np.load(init_Sigma), init_gamma, init_k, init_sigma_sqr
        self.yhats, minMSE, minNotPaired = self.run_regression_equation(initialTime=True)
        argmin = (self.stacked_beta, self.theta, self.Sigma, self.gamma, self.k, self.sigma)
        for iteration in range(iters):
            # update hyperparameters
            self.update_parameters()
            #write hyperparameters for plotting later
            self.report_iteration_hyperparameters(iteration)
            yhats, MSE_sum, numNP = self.run_regression_equation()
            #update argmin if appropriate
            if (iteration > burn_in) and (MSE_sum < minMSE):
                argmin = (self.stacked_beta, self.theta, self.Sigma, self.gamma, self.k, self.sigma)
                minMSE = MSE_sum
                minNotPaired = numNP
            #write MSE, numNP, and ratioNP for plotting later
            self.report_metrics(iteration, MSE_sum, numNP, minNotPaired)
        #report argmin
        self.report_iteration_hyperparameters('argmin')

main()
