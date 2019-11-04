#!/usr/bin/env python3


import sys
import numpy as np
from numpy.linalg import qr, inv
import pickle


class regress_gene_cre():
    def __init__(self, statepref, exp_file, cre_file, state_by_chr_file, m_thresh, s_thresh, e=None, lessone=0, tssdist=5, maxdist=int(1e6), cut=0.2, B=100, fixEffect=False):
        """ statepref: string for file with columns for binID, chr, posStart, posEnd, celltype1, ..., cellTypeN, posClass; pknorm_2_16lim_ref1mo_0424_lesshet
           exp_file: string for file with chr tss gene_category strand celltype1, ..., celltypej; rnaTPM.txt
           cre_file: string for file with chr, start, stop; vision_cres.txt
           state_by_chr_file: string for npz file with preprocessed state and pos arrays
           m_thresh: int, threshold for mean of expression (gene group subselection)
           s_thresh: int, threshold for stdev of expression (gene group subselection)
           e: float array of beta coefficient, dim is # states; must be None initially
           lessone: int, leave-one-out cellType value {0,1,2,3,4,5,6,7,8,9,10,11}; default 11
           tssdist: int, number of bins around TSS for proximal CREs where bin is 200 bp; default 5
           maxdist: int, distance in bp for distal CRE window; default 1e6
           cut: float, correlation cutoff for initial CRE filtering; default 0.2
           B: int, number of refinement iterations; default 100
           fixEffect: boolean, ?; default False
        """
        self.m_thresh = m_thresh
        self.s_thresh = s_thresh
        self.e = e
        self.lessone = lessone
        self.tssdist = tssdist
        self.maxdist = maxdist
        self.cut = cut
        self.B = B
        self.fixEffect = fixEffect

        '''read in CREs -> self.tcre (subselected by self.chr)'''
        self.tcre_all = self.load_CREs(cre_file)

        '''read in the RNA exp info -> self.rna (subselected by self.chr)
                        cell type header info -> self.rna_names
                        Transcription start sites -> self.tss
            set
            set number of cell types -> self.cellN
            '''
        self.rna_all, self.rna_names, self.tss_all = self.load_RNA(exp_file)
        self.chroms = list(self.rna_all)
        self.cellN = self.rna_all[self.chroms[0]].shape[1]

        '''read in the genomic state info -> self.state (subselected by self.chr)
                        bin information -> self.pos (subselected by self.chr)
            '''
        #self.state, self.pos = self.load_state(statepref)
        self.state_all, self.pos_all = self.load_state_from_pickle(state_by_chr_file) #sooo much faster

        print('Done Initializing')

    def load_RNA(self, exp_file): #rnaTPM.txt
        rna = []
        chromosomes = []
        tss = []
        with open(exp_file) as f:
            rna_names = f.readline().split()[4:]
            for line in f:
                fields = line.strip('\r\n').split()
                rna.append(fields[4:])
                chromosomes.append(fields[0])
                tss.append(round(float(fields[1])))


        rna = np.array(rna, dtype=np.float64)
        chromosomes = np.array(chromosomes)
        tss = np.array(tss, dtype=np.int32)
        tss = (tss // 200).astype(np.int32).reshape(-1,1) # convert TSS coordinates into 200 bp bins

        '''subselect chromosomes'''
        chroms = np.unique(chromosomes)
        rna_all = {}
        tss_all = {}
        for chrom in chroms:
            where = np.where(chromosomes == chrom)[0]
            rna_all[chrom] = rna[where]
            tss_all[chrom] = tss[where]

        return (rna_all, rna_names, tss_all)


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
        pos = state[:,2].astype(np.int32)//200 #state bins (coordinate --> 200bp bin space)
        state = state[:,4:].astype(np.int32) #remove state position data, so just states left
        valid = [header.index(x) for x in self.rna_names if x in header]
        #CFU_E_ad   CFUMK   CMP ERY_fl  GMP MK_imm_ad   LSK_BM  MEP MONO_BM NEU ER4 G1E
        #CFU_E_ad   CFUMK   CMP ERY_fl  GMP MK_imm_ad   LSK_BM  MEP MONO_BM NEU ER4 G1E
        state = state[:,valid] #select only cell types RNA data; in essence subselecting and reordering the state cell types to match the rna cell types

        '''subselect chromosomes'''
        chroms = np.unique(scr)
        state_all = {}
        pos_all = {}
        for chrom in chroms:
            t = np.where(scr == chrom)[0] #extract only relevant chromosomes genes/states
            pos_all[chrom] = pos[t] #extract only relevatn chromosome state bins
            state_all[chrom] = state[t,:] #extract only relevant chromosomes states

        return (state_all, pos_all)

    def load_state_from_pickle(self, state_by_chr_file):
        with open(state_by_chr_file, 'rb') as picklefile:
            pickled = pickle.load(picklefile)
            state_all = pickled['state_all']
            pos_all = pickled['pos_all']
        return(state_all, pos_all)

    def load_CREs(self, cre_file): #vision_cres.txt
        tcre = []
        with open(cre_file) as f:
            for line in f:
                fields = line.strip('\r\n').split()
                tcre.append((fields[0], (int(fields[1]), int(fields[2]))))

        tcre = np.array(tcre, dtype=np.dtype([('chr', 'U5'), ('coord', np.int32, (2,))]))
        tcre['coord'][:,0] = tcre['coord'][:,0] // 200
        tcre['coord'][:,1] = (tcre['coord'][:,1]+199) // 200 +1

        '''subselect chromosomes'''
        tcre_all = {}
        chroms = np.unique(tcre['chr'])
        print(chroms)
        for chrom in chroms:
            where = np.where(tcre['chr']==chrom)[0]
            tcre_all[chrom] = tcre['coord'][where, :]
        return (tcre_all)

    def lm(self, x,y):
        """x: 2dim rows is samples; columns is features
        y: 1dim, samples"""

        print("I got to lm")

        Q,R = qr(x) #Use QR decomposition to solve linear regression
        coeff = inv(R).dot(Q.T).dot(y.reshape(-1,1))
        coeff[np.where(np.isnan(coeff))] = 0
        yhat = np.sum(x * coeff.T, axis=1)
        R2 = np.corrcoef(y,yhat)[0,1]**2
        n,k = x.shape
        #R2adj = R2 - k/ (n-k-1) * (1-R2)
        R2adj = 1 - ((1-R2)*(n-1)/(n-k-1))
        print("I finished lm")
        return {'coeff': coeff, 'R2': R2, 'R2adj': R2adj}

    def genereg(self, rna, state, pair, prior=None):
        """
        rna: float array of RNA values, #TSS by #cell_types
        state: float array of state masks, #positions by #states*#cell_types; mid position of CREs are mean state values across CRE position range
        pair: numpy object of (TSS, CRE positions, TSS index, predicted CRE position contributions)
        prior: ?
        """

        some_stateN = state.shape[1] // self.cellN
        some_l = state.shape[0]
        y = []
        x0 = []
        x = []
        ut = np.unique(pair['TSSidx']) #Find unique TSS indices -- all the TSS that have at least one ccRE pair to it

        self.utN = ut.shape[0]
        self.ut = ut

        for i in range(self.utN): #for each unique TSS
            tt = np.where(pair['TSSidx'] == ut[i])[0] #Find entries in pair for current TSS idx
            if tt.shape[0] == 0:
                continue
            y.append(rna[ut[i],:]) #add cell type expression data to y for current TSS
            ttt = np.arange(max(0, pair['TSS'][tt[0]] - self.tssdist),
                            min(some_l, pair['TSS'][tt[0]] + self.tssdist +1)) #get position range window around TSS
            x0.append(np.mean(state[ttt,:], axis=0)) #Find state mean across TSS window for each cell type by state
            if prior is not None:
                x.append(np.mean(prior[pair['TSS'][tt], pair['CRE'][tt]].reshape(-1,1)
                                * state[pair['CRE'][tt],:], axis=0))
            else:
                x.append(np.mean(state[pair['CRE'][tt],:], axis=0)) #Get CRE mean states for TSS
        y = np.array(y, dtype=np.float32)
        x = np.array(x, dtype=np.float32)
        x0 = np.array(x0, dtype=np.float32)
        xx = []
        xx0 = []
        for i in range(self.cellN):
            xx += list(x[:,(i*self.stateN):((i+1)*self.stateN)])
            xx0 += list(x0[:,(i*self.stateN):((i+1)*self.stateN)])
            # Reshape x and x0 arrays to 2D [cellN * # TSSs, stateN]
            # So, for each cell type, there is a 2D matrix of TSSs by TSS(CRE) mean state values
        xx = np.array(xx, dtype=np.float32)
        xx0 = np.array(xx0, dtype=np.float32)

        rt = {
            'y': y.ravel(order='F'), #float array of size #TSS * #celltypes
            'z': ((y-np.mean(y, axis=1, keepdims=True))
                /(np.std(y, axis=1, keepdims=True, ddof =1) + 1e-3)).ravel(order='F'),
            'x': xx,
            'x0': xx0
             }

        np.savetxt('my_rt_x.txt', rt['x'])
        np.savetxt('my_rt_x0.txt', rt['x0'])
        np.savetxt('my_rt_z.txt', rt['z'])
        np.savetxt('my_rt_y.txt', rt['y'])

        return (rt)



    def refineList(self, rt, state, pair, sel, itern):
                        #rt, ss, self.pair
        """
        rt: dictionary with
            y: float array, RNA expression across cell types, dim (cellN * #TSSs)
            x: float array, mean states across all CREs for a given TSS, dim (cellN * #TSSs, stateN)
            x0: float array, mean states across TSSdist window for a given TSS, dim (cellN* #TSSs, stateN)
        state: float array, state masks, dim (# positions, stateN * cellN); mid position of CREs are mean state values across CRE position range
        pair: numpy object of TSS-CRE pairs (TSS, CRE positions, TSS index, predicted CRE position contributions)
        sel: int array, currently selected CREs
        lessone: int, leave one out cell type value
        itern: int, iteration number
        """

        y = rt['y']
        x = rt['x']
        x0 = rt['x0']

        ut = np.unique(pair['TSSidx'])
        n = ut.shape[0]
        cellN = y.shape[0] / n
        stateN = x.shape[1]

        tx = np.hstack((np.log2(x+0.001), np.log2(x0 + 0.001))) #n by stateN *2
        #print("tx shape: ", tx.shape)
        #****#*********#*************#***********#
        e = np.zeros(((self.stateN * 2) + 1, self.cellN -1), dtype=np.float32) # the new betas in the refined list for x and x0s effect on transcription across cell cell_types
        #****#*********#*************#***********#
        #print("e shape: ", e.shape)
        f = np.zeros((self.utN, self.cellN-1), dtype=np.float32) # the new expression predictions
        #print("f shape: ", f.shape)
        index = 0
        for i in range(self.cellN):
            if i == self.lessone:
                continue
            mask = np.r_[np.arange(self.utN * min(i, self.lessone)),
                        np.arange(self.utN * (min(i,self.lessone)+1), n*max(i,self.lessone)),
                        np.arange(self.utN * (max(i, self.lessone) +1), y.shape[0])]
            #print("Mask: ", mask.shape)
            # create a mask for all TSSs not in i or lessone cell type; leave out lessone cell type and cell type i
            r = self.lm(np.hstack((np.ones((tx.shape[0],1), dtype=np.float32),tx))[mask, :], y[mask])
            #do linear regression with limited set of cell type data
            te = r['coeff'].reshape(-1)
            #****#*********#*************#***********#
            te[np.isnan(te)==True] = 0
            #****#*********#*************#***********#
            #print("te: ", te.shape)
            e[:,index] = te #uses index becasue it's not sure which i
            # test1 = np.ones((self.utN,1))
            # test2 = tx[(i*self.utN):((i+1)*self.utN),:]
            # print("1: ",test1.shape)
            # print("2: ", test2.shape)
            # test3 = np.hstack((test1, test2))
            # print("3: ", test3.shape)
            # test4 = np.dot(test3, te)
            # print("4: ", test4.shape)
            f[:,index] = np.dot(np.hstack((np.ones((self.utN, 1), dtype=np.float32),
                                            tx[(i*self.utN):((i+1)*self.utN),:])), te)
            #create expression prediction for current cell type; for each TSS predict the left out cell type i expression
            index += 1

        #****#*********#*************#***********#
        mm = np.sum((np.r_[y[:(self.lessone * self.utN)], y[((self.lessone + 1) * self.utN):]].reshape(self.utN, self.cellN-1, order='F')
                    - f) ** 2, axis=1) #Find mean squared error for total set
        #****#*********#*************#***********#
        nx = np.copy(x)

        for i in range(self.utN): #For each TSS
            t = np.where(pair['TSSidx'] == ut[i])[0] #list of ccREs for a given TSS - mask used with sel
            if t.shape[0] == 0:
                continue
            me = np.zeros(t.shape[0], dtype=np.float32)
            for j in range(t.shape[0]):
                ttt = ((x[np.arange(self.cellN) * self.utN  + i, :] * np.sum(sel[t])
                        - (2 * sel[t[j]] - 1) * state[pair['CRE'][t[j]], :].reshape(self.cellN, self.stateN, order='C'))
                        / (np.sum(sel[t]) - (2 * sel[t[j]] - 1) + 1e-10))
                ttt[np.where(ttt < 0)[0]] = 0
                #np.savetxt('my_ttt_298.txt', ttt)
                f = np.dot(np.hstack((np.ones((self.cellN,1), dtype=np.float32),
                                      np.log2(ttt + 0.001),
                                      np.log2(x0[np.arange(self.cellN)* self.utN + i, :] + 0.001))), e)
                me[j] = np.sum((y[np.r_[np.arange(self.lessone), np.arange(self.lessone + 1, self.cellN)]
                                * self.utN + i] - f) ** 2)
            j = np.where(me < mm[i])[0]
            if j.shape[0] == 0:
                continue
            tp = np.exp((mm[i] - me[j] - np.amax(mm[i] - me[j])) / 2)
            tp[np.where(np.isnan(tp))] = 0
            tp += 1e-10 / j.shape[0]
            if j.shape[0] > np.round(1000 / (itern+1)) + 1:
                j = j[np.random.choice(j.shape[0], int(np.round(100 / (itern+1))) + 1, p=tp)]
            ttt = ((x[np.arange(self.cellN) * self.utN + i, :] * np.sum(sel[t])
                    - np.sum(state[pair['CRE'][t[j]], :].reshape(1, -1)
                            * (2 * sel[t[j]] - 1). reshape(-1, 1),
                            axis=0).reshape(self.utN, self.cellN, order='F').T)
                    / (np.sum(sel[t]) - np.sum(2 * sel[t[j]] -1) + 1e-10))
            ttt[np.where(ttt < 0)] = 0
            f = np.dot(np.hstack((np.ones((ttt.shape[0], 1), dtype=np.float32),
                                  np.log2(ttt + 0.001),
                                  np.log2(x0[np.arange(self.cellN) * self.utN + i, :] + 0.001))), e)
            f = np.diag(f[:,0])
            mm[i] = np.sum((y[np.r_[np.arange(self.lessone),
                                    np.arange(self.lessone + 1, self.cellN)] * self.utN + i] -f) ** 2)
            nx[np.arange(self.cellN) * self.utN + i, :] = ttt
            sel[t[j]] = 1 - sel[t[j]]
            
        rt = {'x': nx,
              'sel': sel}
        return rt

    def runRefine(self, rt, ss, pair):
        """
        rt: dict of (y: RNA expression for each TSS, z: normalized RNA, x: mean CRE state, x0: mean TSS state)
        cellN: number of cell types
        ss: float array of state masks, #positions by #states *#cell_types; mid position of CREs are mean state values across CRE position range
        pair: numpy object of TSS-CRE pairs (TSS, CRE positions, TSS index, predicted CRE position contributions)
        lessone: int, cell type to skip in leave-one-out
        """
        r = {'x':np.copy(rt['x'])}
        print('x shape: ', rt['x'].shape)
        #k = rt['y'].shape[0]/ self.cellN k==self.utN

        #a00_x = rt['x0'][(self.lessone * self.utN):((self.lessone + 1)*self.utN),:]
        #print(a00_x.shape)
        #np.savetxt('my_a00_x_minus4_2.txt', a00_x, fmt = '%.5f')

        a00 = self.lm(np.log2(rt['x0'][(self.lessone * self.utN):((self.lessone + 1)*self.utN),:]+0.001),
                        rt['y'][(self.lessone*self.utN):((self.lessone+1)*self.utN)])['R2adj']
        ma = self.lm(np.hstack((np.log2(r['x'][(self.lessone*self.utN):((self.lessone+1)*self.utN),:] + 0.001),
                                np.log2(rt['x0'][(self.lessone*self.utN):((self.lessone+1)*self.utN),:] + 0.001))),
                    rt['y'][(self.lessone*self.utN):((self.lessone + 1)*self.utN)])['R2adj']
        ms = np.ones(self.pair.shape[0], dtype=np.int32)
        sel = np.ones(self.pair.shape[0], dtype=np.int32)
        a0 = ma
        n0 = 1.0
        a = np.zeros(self.B, dtype=np.float32)
        n = np.zeros(self.B, dtype=np.float32)
        print("\r%s\rRound:%i n:%f a:%f" % (' ' * 80, 0, n0, a0), end='', file=sys.stderr)

        for i in range(self.B):
            r = self.refineList(rt, ss, pair, sel, i)
            if np.sum(sel != r['sel']) == 0: #no change in the selected list
                break #No more changes, end early
            sel = np.copy(r['sel']) #sel equal to current/returned
            a[i] = self.lm(np.hstack((np.log2(r['x'][(self.lessone*self.utN):((self.lessone+1)*self.utN),:]+0.001),
                                      np.log2(rt['x0'][(self.lessone*self.utN):((self.lessone+1)*self.utN),:]+0.001))),
                            rt['y'][(self.lessone*self.utN):((self.lessone+1)*self.utN)])['R2adj'] #track the R2 performance for each run
            n[i] = np.mean(sel) #Percentage of retained CRMs
            print("\r%s\rRound:%i n:%f a:%f" % (' ' * 80, i + 1, n, a), end='', file=sys.stderr)
            if a[i] > ma: #retain best performing set; if you have achieved a better performance, update the best performing
                ma = a[i]
                ms = np.copy(sel)

        r = {'n': n[:(i + 1)],
             'n0': n0,
             'a': a[:(i + 1)],
             'a0': a0,
             'a00': a00,
             'ma': ma,
             'msel': ms,
             'sel': sel}

        return r

    def run(self, chrom):
        #****#*********#*************#***********#
        self.e = None
        #****#*********#*************#***********#

        '''subselect all data for specific chromosome '''
        self.rna = self.rna_all[chrom]
        self.state = self.state_all[chrom]
        self.tcre = self.tcre_all[chrom]
        self.pos = self.pos_all[chrom]
        self.tss = self.tss_all[chrom]

        '''subselect expression data for specific average and stdev '''
        m = np.mean(self.rna, axis=1)
        s = np.std(self.rna, axis=1, ddof=1)
        tt = np.where((m>self.m_thresh) & (s> self.s_thresh))[0]

        self.rna = self.rna[tt]
        self.tss = self.tss[tt]
        self.tssN = self.tss.shape[0]

        '''now set up pairs '''
        G = np.amax(self.state) + 1 #Determine number of states
        self.stateN = G
        t = np.bincount(self.state.ravel(order='F'), minlength=self.stateN) #Count occurence of each state
        k = np.where(t == np.amax(t))[0][0] #Determine most prevelant state
        l = max(np.amax(self.tss), np.amax(self.pos)) + 1 #Determine largest number (because of adding 1) of bins
        self.l = l

        ss = np.full((l, self.cellN), k, dtype=np.int32) #Create array of state bins (including skipped regions and fill with most prevelant/ common state)
        ss[self.pos, :] = self.state #Fill in valid & known states

        cre = np.full(self.l, -1, dtype=np.int32) #Create a cre array of size equal to # bins
        for i in range(self.tcre.shape[0]): #indice of the ccre is the row in tcre
            cre[self.tcre[i,0]:min(self.tcre[i,1],self.l)] = i #mark positions and indices of CREs, in this vector, if the crre overlaps that genome bin, set the element equal to the index (upstream will be overwritten by downstream with this approach if overlap)

        if self.e is None:
            sss = np.zeros((self.tssN * self.cellN, self.stateN), dtype=np.float32) #stacked two dimensional array, blocked by cell type (consider moving to 3dim array?)
            for i in range(self.tssN):
                ttt = np.arange(max(0, self.tss[i] - self.tssdist), min(self.tss[i]+ self.tssdist+1, self.l)) #Find bin window around TSS (1000bp) looking 5 bins down and 5 bins up and the one it's actually in, so 11 total
                temp = ss[ttt,:] #Find state counts for each cell type in window; temp is an array of bin x celltype; should be shape(11, 12)
                aaa = np.empty((temp.shape[1], self.stateN), dtype=np.float32)
                for j in range(self.cellN):
                    aaa[j,:] = np.bincount(temp[:,j], minlength=self.stateN) #finding the fraction of each state in that window around the TSS for each cell type
                aaa /= np.sum(aaa, axis=1, keepdims=True) #convert to state proportions
                sss[np.arange(self.cellN)*self.tssN + i,:] = aaa #sss has list of TSS state proportions, grouped by cell type
            self.e = self.lm(sss, self.rna.ravel(order='F'))['coeff'][:,0] #Use QR decomposition to solve linear regression; e is a vector of 27 coefficients
        sse = self.e[ss] #get state coefficients for each bin and cell type; returns an array of exact same shape as ss, but sse_i,j takes on the values of e[ss_i,j]
        tr = ((self.rna - np.mean(self.rna, axis=1, keepdims=True))
                / ((np.std(self.rna, axis=1, keepdims=True, ddof =1) + 1e-5) *(self.cellN -1) ** 0.5 ))

        # rna is unnoramlized array, tr is normalized array
        #sse is unnormalized beta coefficient, ts is normalized beta coeff per genomic bin across cell types
        ts = ((sse - np.mean(sse, axis=1, keepdims=True))
                / ((np.std(sse, axis=1, keepdims=True, ddof =1) + 1e-5) * (self.cellN -1) ** 0.5 ))#.astype(np.int32)

        pair = []
        for i in range(self.tssN):
            a = np.arange(max(0, self.tss[i] - self.maxdist // 200), min(self.l, self.tss[i] + self.maxdist // 200 + 1)) #put maxdist into bins and get bin range of window maxdist around TSS
            #****#*********#*************#***********#
            tss_local = np.where(a==self.tss[i])[0]
            #****#*********#*************#***********#
            #np.savetxt('my_a_{}.txt'.format(i), a)
            rr = np.dot(tr[i:(i+1), :], ts[a,:].T).ravel(order='F') #predict state contribution to rna
            #np.savetxt('my_rr_{}.txt'.format(i), rr)
            trr = np.dot(tr[i:(i+1), np.r_[np.arange(self.lessone), #predict position correlation to RNA leaving out one cell type specified by lessone
                                           np.arange(self.lessone +1, self.cellN)]],
                        (ts[a,:][:,np.r_[np.arange(self.lessone),
                                    np.arange(self.lessone +1, self.cellN)]]).T).ravel(order='F')

            t = np.where((trr >= self.cut) & (cre[a] >= 0))[0]

            if t.shape[0] > 0:
                nt = []
                for j in np.unique(cre[a[t]]): #a[t] is global position, cre[a[t]] is the cre index at that position. Since a ccRE can overlap multiple bins, take unique
                    tt = t[np.where(cre[a[t]]==j)] #Identify each CRE position(s) of valid state in window; tt is local position of valid CRE
                    nt.append(tt[np.where(np.around(trr[tt], decimals=5)==np.amax(np.around(trr[tt], decimals=5)))[0]][0]) #determine CRE position (if it spans multiple bins) with highest contribution/correlation and append to nt

                    #****#*********#*************#***********#
                    ntt = tt[np.where(np.around(trr[tt], decimals=5)==np.amax(np.around(trr[tt], decimals=5)))[0]]
                    if tss_local in ntt:
                        nt.append(tss_local)
                    else:
                        nt.append(ntt[0])
                    #nt.append(ntt[0])
                    #****#*********#*************#***********#

                t = np.array(nt) #refining t to be the single position that has maximum correlation for each ccRE
            #t = np.unique(np.r_[np.where(a == self.tss[i])[0], t]) #appending the TSS bin if it's not there; make position list only unique entries, including TSS bin

            #****#*********#*************#***********#
            t = np.unique(t)
            ttss = np.where(a == self.tss[i])[0]
            if not np.isin(ttss, t)[0]:
                t = np.r_[ttss,t]
            else:
                where = np.where(t == ttss)[0][0]
                t = np.r_[ttss, t[:where], t[(where+1):]]
                #t = np.delete(t, np.where(t==ttss)[0])
                #t = np.r_[ttss,t]
            #****#*********#*************#***********#

            #np.savetxt('my_t_{}.txt'.format(i), t)

            for j in range(t.shape[0]):
                pair.append((self.tss[i], a[t[j]], i, rr[t[j]], cre[a[t[j]]])) #TSS, global cre bin position, tss index, predicted correlation without leave one out
                            #TSS          CRE   TSSidx  rr      CREidx


        pair = np.array(pair, dtype=np.dtype([('TSS', np.int32), ('CRE', np.int32), #convert pair into a numpy array
                                                ('TSSidx', np.int32), ('rr', np.float32),
                                                ('CREidx', np.int32)]))

        #np.savetxt('my_pair_TSS.txt', pair['TSS'])
        #np.savetxt('my_pair_CRE.txt', pair['CRE'])
        #np.savetxt('my_pair_TSSidx.txt', pair['TSSidx'])
        #np.savetxt('my_pair_rr.txt', pair['rr'])
        #np.savetxt('my_pair_CREidx.txt', pair['CREidx'])
        print("PAIR SHAPE: ", pair.shape)

        self.pair = pair
        #output = open('tmp_%s_py.txt' % chrom, 'a')
        #self.pair['TSS'] += 1
        #self.pair['CRE'] += 1
        #self.pair['TSSidx'] += 1
        #self.pair['CREidx'] += 1
        #for line in self.pair:
        #   print( "%s %i %i %i %0.15f %i" % (chrom, line[0], line[1], line[2], line[3], line[4]), file=output)
        #output.close()
        #return

        kk = np.zeros(self.stateN, dtype=np.int32) #create an array of size #ofstates
        #print("kk: ", kk.shape)
        kk[k] = 1 #mark most prevelant state
        #ss = np.tile(np.repeat(kk.reshape(1, -1), repeats=self.l, axis=0), self.cellN) #make a state mask array; horizontal stacking; 2d array
        ss = np.tile(np.repeat(kk.reshape(1, -1), repeats=self.l, axis=0), (1, self.cellN)).astype(np.float32)

        for i in range(self.state.shape[0]): #for each global position
            tt = np.zeros(self.cellN*self.stateN, dtype=np.int32)

            tt[np.arange(self.cellN) * self.stateN + self.state[i,:]] = 1 #one hot encoding all the states for all the cell types
            ss[self.pos[i], :] = tt #update state mask with valid position states for each cell type
        # if np.product(self.tcre.shape) > 2: #if more than 1 tcre
        #     oss = np.copy(ss) #storing your mask array before you mess with it
        #     ss = ss.astype(np.float32)
        #     for i in range(self.tcre.shape[0]): #for each CRE
        #         ss[(self.tcre[i,0] + self.tcre[i,1]) //2, :] = np.mean(oss[self.tcre[i,0]:self.tcre[i,1],:], axis=0).astype(np.float32) #one-hot components and float components because ss is all genomic positions instead of all ccREs
                #Find mean state values for CREs across all bins they span for each cell type


        if self.fixEffect:
            ss = 2 ** sse

        rt = self.genereg(self.rna, ss, self.pair, prior=None)
        r = self.runRefine(rt, ss, self.pair)

        rt['nx'] = r['x']
        rt['sel'] = r['sel']
        rt['msel'] = r['msel']
        rt['n'] = r['n']
        rt['n0'] = r['n0']
        rt['a'] = r['a']
        rt['a0'] = r['a0']
        rt['a00'] = r['a00']
        rt['ma'] = r['ma']
        rt['pair'] = r['pair']
        print("\r%s\r" % (' ' * 80), end='', file=sys.stderr)
        return rt

statepref = '/home/kweave23/VISION_regression/their_stuff/pknorm_2_16lim_ref1mo_0424_lesshet'
exp_file = '/home/kweave23/VISION_regression/their_stuff/rnaTPM.txt'
cre_file = '/home/kweave23/VISION_regression/their_stuff/vision_cres.txt'
#atacsig_file = '/home/kweave23/VISION_regression/their_stuff/vision_cres.mat.atacsig.txt'
state_by_chr_file = '/home/kweave23/VISION_regression/state_all_and_pos_all_by_chr.pickle'

test1 = regress_gene_cre(statepref, exp_file, cre_file, state_by_chr_file, -4, 2)


for chrom in ['chr1', 'chr11', 'chr19']:
    test1.run(chrom)
