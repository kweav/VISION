#!/usr/bin/env python3

import argparse as ap
import os
import numpy as np
from numpy.linalg import qr, inv, det
import pickle
from sklearn import linear_model

def main():
    parser = generate_parser()
    args = parser.parse_args()
    setup_file_locs(args, args.comp)
    setup_threads(args.threads)
    model = regress_gene_cre(args.statefile, args.exp_file, args.cre_file, args.binsize)

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
        model.set_initial_betas(args.initdist)
        model.find_initial_pairings(args.lessone, args.tss_dist, args.distal_dist, args.correlation, args.gamma)
        #model.refine_pairs()


def generate_parser():
    parser = ap.ArgumentParser(description = 'VISION regression to predict gene expression through IDEAS states and cCREs incorporating distance')
    #parser.add_argument('-s','--statepref', dest='statefile', action='store', type=str, default='/Users/kateweaver/taylorLab/VISION/replicating_regression/their_inputs/pknorm_2_16lim_ref1mo_0424_lesshet.state', help='default is /home/kweave23/VISION_regression/their_stuff/pknorm_2_16lim_ref1mo_0424_lesshet.state; either this argument or -f')
    parser.add_argument('--where', action='store', dest='comp', type=str, default='mine', help='{comp, mine}; adds path to files')
    parser.add_argument('-f','--state_by_chr_file', action='store', dest='statefile', type=str, default='state_all_and_pos_all_by_chr.pickle', help='default is state_all_and_pos_all_by_chr.pickle; either this argument or -s; pickling assumed 200bp bins')
    parser.add_argument('-r','--exp_file', action='store', dest='exp_file', type=str, default='rnaTPM.txt', help=' default is rnaTPM.txt')
    parser.add_argument('-c','--cre_file', action='store', dest='cre_file', type=str, default='vision_cres.txt', help='default is vision_cres.txt')
    parser.add_argument('--chroms', action='store', nargs='+', dest='chroms', type=str, default='all')
    parser.add_argument('-o','--output_file_name', action='store', dest='output_name', type=str, default='vision_rna_tss2k_ccreunit', help='default is vision_rna_tss2k_ccreunit')
    parser.add_argument('-g', '--group', dest='exp_type', type=int, action='store', default=0, help='{0, 1, 2, 3, 4} Group subselection. 0 is for all and is default')
    parser.add_argument('-e', '--lessone', dest='lessone', type=int, action='store', default=0, help='Cell type to leave out with 0 indexing. 0 is default')
    parser.add_argument('-i', '--iterations', dest='iters', type=int, action='store', default=100, help='Refinement iterations')
    parser.add_argument('-t', '--threads', dest='threads', type=str, action='store', default="1", help='default is 1; max number of threads to use')
    parser.add_argument('-z', '--initialization_dist', dest='initdist', type=int, action='store', default=1000, help='Beta initialization distance cutoff')
    parser.add_argument('-b', '--binsize', dest='binsize', type=int, action='store', default=200)
    parser.add_argument('-p', '--tss_dist', dest='tss_dist', type=int, action='store', default=1000, help='tss_window cutoff')
    parser.add_argument('-d', '--distal_dist', dest='distal_dist', type=int, action='store', default=1000000, help='Distal distance cutoff')
    parser.add_argument('--gamma', dest='gamma', type=float, action='store', default=0.7, help='gamma for distance function. Suggested 0.7 or 1 for the Extrusion or Fractal Globule models respectively')
    parser.add_argument('--tss', dest='tss', action='store_true', help='Force inclusion of TSS in cCRE set')
    parser.add_argument('--correlation', dest='correlation', action='store', type=float, default=0.2, help='Initial pairing correlation cutoff')
    parser.add_argument('--mean_thresh', dest='m_thresh', action='store', type=int, default=-4)
    parser.add_argument('--stdev_thresh', dest='s_thresh', action='store', type=int, default=2)
    parser.add_argument('-l', '--log', dest="log", action='store_true', help='Log2-transform predictors')
    #parser.add_argument('-k', '--skip', dest='which_skip', action='store', type=str, default='pm', help='options in the following set: {pm, px, pmd, pxd, pmpx} where pm is promoter and px is proximal and d is distal')
    parser.add_argument('-m', '--multi_refinement', dest='multirefine', action='store_true', help='Change multiple cCRE inclusion status per gene per round')
    return parser

def setup_threads(threads):
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["OPENBLAS_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = threads
    os.environ["NUMEXPR_NUM_THREADS"] = threads

def setup_file_locs(args, comp):
    if comp=='mine':
        args.exp_file = '/Users/kateweaver/taylorLab/VISION/replicating_regression/their_inputs/' + args.exp_file
        args.cre_file = '/Users/kateweaver/taylorLab/VISION/replicating_regression/their_inputs/' + args.cre_file
        args.statefile = '/Users/kateweaver/taylorLab/VISION/regression_with_dist/' + args.statefile
    elif comp=='comp':
        args.exp_file = '/home/kweave23/VISION_regression/their_stuff/' + args.exp_file
        args.cre_file = '/home/kweave23/VISION_regression/their_stuff/' + args.cre_file
        args.statefile = '/home/kweave23/VISION_regression/' + args.statefile


class regress_gene_cre():
    '''Instance variables ordered by time of definition
    self.binsize: binsize set in __init__
    self.rna_all: dictionary of TPM values by chromosome set in __init__
    self.rna_names: set in __init__
    self.tss_all: set in __init__
    self.rna_info_all: set in __init__
    self.cellN: number of cell types set in __init__
    self.state_all: set in __init__
    self.pos_all: set in __init__
    self.cre_all: set in __init__
    self.rna: set in subset_based_on_chrom first, refined in subset_based_on_group
    self.state: set in subset_based_on_chrom
    self.stateN: set in subset_based_on_chrom
    self.cre: set in subset_based_on_chrom
    self.creN: set in subset_based_on_chrom
    self.pos: set in subset_based_on_chrom
    self.tss: location of TSS's in bins (not coordinates) set in subset_based_on_chrom first, refined in subset_based_on_group
    self.m_thresh: mean threshold for RNA expression set in subset_based_on_group
    self.s_thresh: standard deviation threshold for RNA expression set in subset_based_on_group
    self.group: RNA expression group set in subset_based_on_group
    self.tssN: number of TSS's set in subset_based_on_group
    self.initbins: half of symmetric initial window size in bins set in set_initial_betas
    self.largest_bin: largest_possible_bin (based on available genome state information) set in set_initial_betas
    self.initial_coeffs: coefficients based on initial regression of states in window around TSS and expression. Note state 0's Beta set to 0
    self.init_betas_acBinsCT: Betas across bins of the genome for every cell type based on self.initial_coeffs set in set_initial_betas
    self.norm_rna: normalized_rna set in set_inital_betas
    self.norm_init_betas: normalized init_betas_acBinsCT set in set_initial_betas
    self.gamma: gamma for distance function set in find_initial_pairings
    self.maxdist: distance for full windows around TSS in bins set in find_initial_pairings
    self.tssdist: distance for tss windows arounds TSS in bins set in find_initial_pairings
    self.lessone: leave-one-out cell type index set in find_initial_pairings
    self.lessone_range: np range exluding leave-one-out cell type index set in find_initial_pairings
    self.tss_windows: dictionary whose key is the tss_n from i for i in range(self.tssN).
    Set in find_initial_pairings: Values include {"tss window": tss_bin_window,
                                                  "full window": full_bin_window,
                                                  "F_booleans": {'mx': full_max,
                                                                 'mn': full_min},
                                                  "T_booleans": {'mx': tss_max,
                                                                'mn': tss_min}}
    self.pair: initial pairings set in find_initial_pairings
    '''
    def __init__(self, statefile, exp_file, cre_file, binsize):
        self.binsize=binsize
        self.rna_all, self.rna_names, self.tss_all, self.rna_info_all = self.load_RNA(exp_file)
        chroms_possible = list(self.rna_all)
        self.cellN = self.rna_all[chroms_possible[0]].shape[1]
        if "pickle" in statefile:
            self.state_all, self.pos_all = self.load_state_from_pickle(statefile)
        else:
            self.state_all, self.pos_all = self.load_state(statefile)
        self.cre_all = self.load_CREs(cre_file)

    def load_RNA(self, exp_file): #rnaTPM.txt
        rna = []
        chromosomes = []
        tss = []
        rna_info = []
        with open(exp_file) as f:
            rna_names = f.readline().split()[4:]
            for line in f:
                fields = line.strip('\r\n').split()
                rna.append(fields[4:])
                chromosomes.append(fields[0])
                tss_val = round(float(fields[1]))
                tss.append(tss_val)
                rna_info.append((fields[0], int(tss_val), fields[2], fields[3]))


        rna = np.array(rna, dtype=np.float64)
        chromosomes = np.array(chromosomes)
        tss = np.array(tss, dtype=np.int32)

        rna_info = np.array(rna_info, dtype=np.dtype([('chr', 'U5'), ('rTSS', np.int32), ('gene_cat', 'U25'), ('strand', 'U1')]))
        tss = (tss // self.binsize).astype(np.int32).reshape(-1,1) # convert TSS coordinates into 200 bp bins

        '''subselect chromosomes to build dictionaries'''
        chroms = np.unique(chromosomes)
        rna_all = {}
        tss_all = {}
        rna_info_all = {}
        for chrom in chroms:
            where = np.where(chromosomes == chrom)[0]
            rna_all[chrom] = rna[where]
            tss_all[chrom] = tss[where]
            rna_info_all[chrom] = rna_info[where]

        return (rna_all, rna_names, tss_all, rna_info_all)


    def load_state(self, statepref): #pknorm_2_16lim_ref1mo_0424_lesshet.state
        file_to_load = "%s.state" %statepref
        state = []
        with open(file_to_load) as f:
            header = f.readline().split()[4:]
            for line in f:
                state.append(line.strip('\r\n').split())

        state = np.array(state, dtype=np.object)

        #state data is a file with columns for binID, chr, posStart, posEnd, celltype1, ..., cellTypeN, posClass
        scr = state[:,1] #state chromosome names
        pos = state[:,2].astype(np.int32)//self.binsize #state bins (coordinate --> 200bp bin space)
        state = state[:,4:].astype(np.int32) #remove state position data, so just states left
        valid = [header.index(x) for x in self.rna_names if x in header]
        #CFU_E_ad   CFUMK   CMP ERY_fl  GMP MK_imm_ad   LSK_BM  MEP MONO_BM NEU ER4 G1E
        #CFU_E_ad   CFUMK   CMP ERY_fl  GMP MK_imm_ad   LSK_BM  MEP MONO_BM NEU ER4 G1E
        state = state[:,valid] #select only cell types RNA data; in essence subselecting and reordering the state cell types to match the rna cell types

        '''subselect chromosomes to build dictionaries'''
        chroms = np.unique(scr)
        state_all = {}
        pos_all = {}
        for chrom in chroms:
            t = np.where(scr == chrom)[0] #extract only relevant chromosomes genes/states
            pos_all[chrom] = pos[t] #extract only relevatn chromosome state bins
            state_all[chrom] = state[t,:] #extract only relevant chromosomes states

        return (state_all, pos_all)

    def load_state_from_pickle(self, state_by_chr_file):
        '''load pickled object with dictionary of state bin position data and states with each of those dictionaries having sub dictionaries with chromosome keys
            arguments:
                - state_by_chr_file: pickled file'''
        with open(state_by_chr_file, 'rb') as picklefile:
            pickled = pickle.load(picklefile)
            state_all = pickled['state_all']
            pos_all = pickled['pos_all']
        return(state_all, pos_all)

    def load_CREs(self, cre_file): #vision_cres.txt
        '''load CREs into dictionary structure where dictionary key is chromosome and value is numpy array object with chr and bin coordinates
            arguments:
                - cre_file: vision_cres.txt '''
        tcre = []
        with open(cre_file) as f:
            for line in f:
                fields = line.strip('\r\n').split()
                tcre.append((fields[0], (int(fields[1]), int(fields[2]))))

        tcre = np.array(tcre, dtype=np.dtype([('chr', 'U5'), ('coord', np.int32, (2,))]))
        tcre['coord'][:,0] = tcre['coord'][:,0] // self.binsize
        tcre['coord'][:,1] = (tcre['coord'][:,1]+199) // self.binsize +1

        '''subselect chromosomes to build dictionary'''
        tcre_all = {}
        chroms = np.unique(tcre['chr'])
        for chrom in chroms:
            where = np.where(tcre['chr']==chrom)[0]
            tcre_all[chrom] = tcre['coord'][where, :]
        return (tcre_all)

    def subset_based_on_chrom(self, chrom):
        '''subselect all data for specific chromosome
            arguments:
                - chrom: '''
        self.rna = self.rna_all[chrom]
        self.state = self.state_all[chrom]
        self.stateN = np.amax(self.state)+1 #Determine number of states; adding 1 so this is including state 0 though want to ignore state 0
        self.cre = self.cre_all[chrom]
        self.creN = self.cre.shape[0]
        self.pos = self.pos_all[chrom]
        self.tss = self.tss_all[chrom]

    def subset_based_on_group(self, group, m_thresh, s_thresh):
        '''subselect expression data for specific average and stdev
            arguments:
                - group:
                - m_thresh:
                - s_thresh: '''
        self.m_thresh = m_thresh
        self.s_thresh = s_thresh
        self.group = group
        if self.group != 0:
            m = np.mean(self.rna, axis=1)
            s = np.std(self.rna, axis=1, ddof=1)
            group_1 = np.where((m <= self.m_thresh) & ( s <= self.s_thresh))[0]
            group_2 = np.where((m <= self.m_thresh) & (s > self.s_thresh))[0]
            group_3 = np.where((m > self.m_thresh) & (s > self.s_thresh))[0]
            group_4 = np.where((m > self.m_thresh) & (s <= self.s_thresh))[0]
            group_dict = {1: group_1,
                          2: group_2,
                          3: group_3,
                          4: group_4}

            self.rna = self.rna[group_dict[self.group]]
            self.tss = self.tss[group_dict[self.group]]
        self.tssN = self.tss.shape[0]

    def linear_regression(self, X, Y, intercept=False):
        '''reshape arrays to match function; remove contribution of 0 state from X; will set this state's contribution to zero afterwards'''
        X = X.reshape(-1, self.stateN)[:,1:]
        Y = Y.reshape(-1,)
        fit_model = linear_model.LinearRegression(fit_intercept=intercept).fit(X, Y)
        model_coeffs = fit_model.coef_
        r_squared = fit_model.score(X, Y)
        return {'coeffs': model_coeffs, 'rsquare': r_squared}

    def find_TSS_bin(self, max_bool, min_bool):

        return TSS_bin

    def adjust_by_distance(self, to_adjust, tss_n, FoT):
        '''
        to_adjust: X array to be adjusted
        tss_n: which i for i in range(self.tssN) because key of self.tss_windows based on this
        FoT: 'F' or 'T' where 'F' is for a full_bin_window based on self.maxdist
                              'T' is for a tss_window based on self.tssdist
        '''
        tss_bin = self.find_TSS_bin(self.tss_windows[tss_n]['{}_booleans'.format(FoT)]['mx'], self.tss_windows[tss_n]['{}_booleans'.format(FoT)]['mn'])
        Y = np.zeros(to_adjust.shape[0]) #adjustment
        return adjusted

    # def get_overlap(self, range_1, range_2):
    #     range_1_len = range_1[1] - range_1[0]
    #     range_2_len = range_2[1] - range_2[0]
    #     overlap_region = [max(range_1[0], range_2[0]), min(range_1[1], range_2[1])]
    #     len_overlap_region = overlap_region[1] - overlap_region[0]
    #     overlap = len_overlap_region/min(range_1_len, range_2_len)
    #     return (overlap)

    def set_initial_betas(self, initdist):
        self.initbins = initdist//self.binsize #put initial distance into bins
        self.largest_bin = max(np.amax(self.tss), np.amax(self.pos)) + 1 #Determine largest number of bins
        allbin_state_array = np.zeros((self.largest_bin, self.cellN), dtype=np.int32)
        allbin_state_array[self.pos, :] = self.state #fill in valid & known states. Everything else is 0 (most likely the most prevalent)

        tss_initial_props = np.zeros((self.tssN, self.cellN, self.stateN))
        for i in range(self.tssN):
            bin_window = np.arange(max(0, self.tss[i]-self.initbins), min(self.tss[i]+self.initbins+1, self.largest_bin))
            state_counts = allbin_state_array[bin_window,:] #all cell types
            single_initial_props = np.apply_along_axis(lambda x: np.bincount(x, minlength=self.stateN), axis=0, arr=state_counts).astype(np.float32).T
            single_initial_props /= np.sum(single_initial_props, axis=1, keepdims=True) #turn into proportions
            tss_initial_props[i] = single_initial_props

        initial_coeffs = self.linear_regression(tss_initial_props, self.rna)['coeffs']
        '''set contribution of state 0 to be 0'''
        self.initial_coeffs = np.hstack(([0], initial_coeffs))
        self.init_betas_acBinsCT = self.initial_coeffs[allbin_state_array]
        self.norm_rna = ((self.rna - np.mean(self.rna, axis=1, keepdims=True))
                    / ((np.std(self.rna, axis=1, ddof=1, keepdims=True) + 1e-5) * (self.cellN  - 1) ** 0.5))
        self.norm_init_betas = ((self.init_betas_acBinsCT - np.mean(self.init_betas_acBinsCT, axis=1, keepdims=True))
                    / ((np.std(self.init_betas_acBinsCT, axis=1, keepdims=True, ddof=1) + 1e-5) * (self.cellN -1) ** 0.5))

    def find_initial_pairings(self, lessone, tss_dist, distal_dist, correlation, gamma):
        self.gamma = gamma #will be used when adjusting based on distance function
        self.maxdist = distal_dist//self.binsize #put max distance into bins
        self.tssdist = tss_dist//self.binsize #put tss distance into bins
        self.lessone = lessone
        self.lessone_range = np.r_[np.arange(self.lessone),
                                    np.arange(self.lessone+1, self.cellN)]

        cre_loc = np.zeros((self.largest_bin, self.creN))
        for i in range(self.creN):
            cre_loc[self.cre[i,0]:min(self.cre[i,1]+1, self.largest_bin),i] = 1

        cre_loc_by_bin = np.sum(cre_loc, axis=1).reshape(-1, 1)

        self.tss_windows = {} #key is tss index, value is range of windows around that TSS
        pair = [] #want to append TSS, CRE, TSSidx, predicted_contribution
        for i in range(self.tssN):
            '''want to track if tss is ceiling(n/2) for a window length n'''
            full_max, tss_max, full_min, tss_min = False, False, False, False
            if max(0, self.tss[i] - self.tssdist) == 0: #if subtraction of tssdist lowers below 0, surely maxdist will too; if tssdist doesn't, maxdist still could
                full_max = True
                tss_max = True
            elif max(0, self.tss[i]- self.maxdist) == 0:
                full_max = True
            if min(self.tss[i] + self.tssdist + 1, self.largest_bin) == self.largest_bin: #if addition of tssdist goes over max, surely maxdist will to. if tssdist doesn't, maxdist still could
                full_min = True
                tss_min = True
            elif min(0, self.tss[i] + self.maxdist +1, self.largest_bin) == self.largest_bin:
                full_min = True

            full_bin_window = np.arange(max(0, self.tss[i] - self.maxdist), min(self.tss[i]+self.maxdist+1, self.largest_bin))
            tss_bin_window = np.arange(max(0, self.tss[i] - self.tssdist), min(self.tss[i]+self.tssdist+1, self.largest_bin))
            self.tss_windows[i] = {"tss window": tss_bin_window,
                                   "full window": full_bin_window,
                                   "F_booleans": {'mx': full_max,
                                                'mn': full_min},
                                   "T_booleans": {'mx': tss_max,
                                                'mn': tss_min}}

            '''these are based just on state'''
            predicted_contribution = np.dot(self.norm_rna[i:i+1,:], self.norm_init_betas[full_bin_window,:].T).ravel(order='F')
            position_correlation = np.dot(self.norm_rna[i:i+1, self.lessone_range],
                                          self.norm_init_betas[full_bin_window,:][:,self.lessone_range].T).ravel(order='F')
            locsOI = np.where((position_correlation >= correlation) & (cre_loc_by_bin >= 1))[0]
            #these locsOI are bins. Want to translate to CREs
            if locsOI.shape[0] > 0:
                for j, bin in enumerate(locsOI):
                    if j ==  0:
                        cres = np.where(cre_loc[bin])[0]
                    else:
                        cres = np.hstack((cres, np.where(cre_loc[bin])[0]))
                valid_cres = np.unique(cres) #unique CREidx's

            '''Let's add state and loc in contribution/correlation'''
            adjusted_state_window_array = self.adjust(self.norm_init_betas[full_bin_window,:], i, 'F')























main()
