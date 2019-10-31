#!/usr/bin/env python


import sys
import numpy as np
from numpy.linalg import qr, inv

class regress_gene_cre():
    def __init__(self, statepref, exp_file, cre_file, m_thresh, s_thresh, e=None, lessone=0, tssdist=5, maxdist=int(1e6), cut=0.2, B=100, fixEffect=False):
        """ statepref: string for file with columns for binID, chr, posStart, posEnd, celltype1, ..., cellTypeN, posClass; pknorm_2_16lim_ref1mo_0424_lesshet
           exp_file: string for file with chr tss gene_category strand celltype1, ..., celltypej; rnaTPM.txt
           cre_file: string for file with chr, start, stop; vision_cres.txt
           chr: string; chrJ where J is the chromosome number/letter
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

        '''read in the RNA exp info -> self.rna (subselected by self.chr)
                        cell type header info -> self.rna_names
                        Transcription start sites -> self.tss
            set number of cell types -> self.cellN
            set number of TSSs -> self.tssN'''
        self.tcre_all = self.load_CREs(cre_file)
        self.rna_all, self.rna_names, self.tss_all = self.load_RNA(exp_file)
        self.chroms = list(self.rna_all)
        #print("rna: ",self.rna.shape)
        #print("rna_names: ",len(self.rna_names))
        #print("tss: ", self.tss.shape)
        self.cellN = self.rna_all[self.chroms[0]].shape[1]
        #self.tssN = self.tss.shape[0]

        '''read in the genomic state info -> self.state (subselected by self.chr)
                        bin information -> self.pos (subselected by self.chr)
            '''
        self.state_all, self.pos_all = self.load_state(statepref)
        #self.state, self.pos = self.load_state_from_npz(state_by_chr_file) #sooo much faster
        #print("state: ",self.state.shape)
        #print("pos: ",self.pos.shape)
        #self.stateN = self.state.shape[1] // cellN
        #self.l = self.state.shape[0]

        '''read in CREs -> self.tcre (subselected by self.chr)'''
        #print("tcre: ", self.tcre.shape)
        print('Done initializing')

    def load_RNA(self, exp_file): #rnaTPM.txt
        rna = []
        chromosomes = []
        tss = []
        with open(exp_file) as f:
            rna_names = f.readline().split()[4:]
            for line in f:
                fields = line.strip('\r\n').split()
                rna.append([float(x) for x in fields[4:]])
                chromosomes.append(fields[0])
                tss.append(round(float(fields[1])))


        rna = np.array(rna, dtype=np.float64)
        '''subselect only the correct chromosomes'''
        chromosomes = np.array(chromosomes)
        tss = np.array(tss, dtype=np.int32)
        tss = (tss // 200).astype(np.int32).reshape(-1, 1)
        chroms = np.unique(chromosomes)
        rna_all = {}
        tss_all = {}
        for chrom in chroms:
            where = np.where(chromosomes == chrom)[0]
            rna_all[chrom] = rna[where]
            tss_all[chrom] = tss[where]

        '''subselect the gene expression type'''
        #m = np.mean(rna, axis=1)
        #s = np.std(rna, axis=1, ddof =1)
        #tt = np.where((m > self.m_thresh) & (s> self.s_thresh))[0]

        #rna = rna[tt]
        #tss = tss[tt]
        #tss = tss.reshape((-1,1))
        #tss = (tss // 200).astype(np.int32) # convert TSS coordinates into 200 bp bins
        return (rna_all, rna_names, tss_all)


    def load_state(self, statepref): #pknorm_2_16lim_ref1mo_0424_lesshet.state
        file_to_load = "%s" %statepref
        state = []
        scr = []
        pos = []
        with open(file_to_load) as f:
            header = f.readline().split()
            for line in f:
                line = line.strip('\r\n').split()
                scr.append(line[1])
                pos.append(int(line[2]))
                state.append([int(x) for x in line[4:]])

        state = np.array(state, dtype=np.int32)
        scr = np.array(scr, dtype='U10')
        pos = np.array(pos, dtype=np.int32)

        #state data is a file with columns for binID, chr, posStart, posEnd, celltype1, ..., cellTypeN, posClass
        #scr = state[:,1] #state chromosome names
        pos = pos // 200 #state bins (coordinate --> 200bp bin space)
        #state = state[:,4:].astype(np.int32) #remove state position data, so just states left

        ########################
        header = header[4:]
        ########################

        valid = [header.index(x) for x in self.rna_names if x in header]
        #CFU_E_ad   CFUMK   CMP ERY_fl  GMP MK_imm_ad   LSK_BM  MEP MONO_BM NEU ER4 G1E
        #CFU_E_ad   CFUMK   CMP ERY_fl  GMP MK_imm_ad   LSK_BM  MEP MONO_BM NEU ER4 G1E
        state = state[:,valid] #select only cell types RNA data; in essence subselecting and reordering the state cell types to match the rna cell types
        chroms = np.unique(scr)
        state_all = {}
        pos_all = {}
        for chrom in chroms:
            t = np.where(scr == chrom)[0] #extract only relevant chromosomes genes/states;
            pos_all[chrom] = pos[t] #extract only relevant chromosome state bins
            state_all[chrom] = state[t,:] #extract only relevant chromosomes states
        return (state_all, pos_all)

    def load_state_from_npz(self, state_by_chr_file):
        with np.load("{}".format(state_by_chr_file)) as npzfile:
            state = npzfile['{}_state'.format(self.chr)]
            pos = npzfile['{}_pos'.format(self.chr)]
        return(state, pos)

    def load_CREs(self, cre_file): #vision_cres.txt
        tcre = []
        with open(cre_file) as f:
            for line in f:
                line = line.strip('\r\n').split()
                tcre.append((line[0], (int(line[1]), int(line[2]))))

        tcre = np.array(tcre, dtype=np.dtype([('chr', 'U5'), ('coord', np.int32, (2,))]))
        tcre['coord'][:, 0] = tcre['coord'][:, 0] // 200
        tcre['coord'][:, 1] = (tcre['coord'][:, 1] + 199) // 200 + 1
        tcre_all = {}
        chroms = np.unique(tcre['chr'])
        for chrom in chroms:
            where = np.where(tcre['chr']==chrom)[0]
            tcre_all[chrom] = tcre['coord'][where, :] #keep only coordinates for relevant chromosome
        #tcre[:,0] = tcre[:,0] // 200 #convert start coord to bin
        #tcre[:,1] = (tcre[:,1] + 199) // 200 #convert end coord to bin, add 199 because of open interval; equiv to (end-1)//bin_size +1

        #########################
        #tcre[:,1] = (tcre[:,1] + 199) // 200 + 1 #convert end coord to bin, add 199 because of open interval; equiv to (end-1)//bin_size +1
        #########################

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
        #print("some_stateN: ", some_stateN)
        #print("self.stateN: ", self.stateN)
        some_l = state.shape[0]
        #print("some_l: ", some_l)
        #print("self.l: ", self.l)
        y = []
        x0 = []
        x = []
        ut = np.unique(pair['TSSidx']) #Find unique TSS indices -- all the TSS that have at least one ccRE pair to it
        self.utN = ut.shape[0]
        #print("utN: ", self.utN)
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
        #print("y: ", y.shape)
        x = np.array(x, dtype=np.float32)
        #print('x: ', x.shape)
        x0 = np.array(x0, dtype=np.float32)
        #print('x0: ', x0.shape)
        xx = []
        xx0 = []
        for i in range(self.cellN):
            xx += list(x[:,(i*self.stateN):((i+1)*self.stateN)])
            xx0 += list(x0[:,(i*self.stateN):((i+1)*self.stateN)])
            # Reshape x and x0 arrays to 2D [cellN * # TSSs, stateN]
            # So, for each cell type, there is a 2D matrix of TSSs by TSS(CRE) mean state values
        xx = np.array(xx, dtype=np.float32)
        #print("xx: ", xx.shape)
        xx0 = np.array(xx0, dtype=np.float32)
        #print("xx0: ", xx0.shape)

        rt = {
            'y': y.ravel(order='F'), #float array of size #TSS * #celltypes
            'z': ((y-np.mean(y, axis=1, keepdims=True))
                /(np.std(y, axis=1, keepdims=True, ddof =1) + 1e-3)).ravel(order='F'),
            'x': xx,
            'x0': xx0
             }

        #print("rt_y: ", rt['y'].shape)
        #print("rt_xx: ", rt['x'].shape)
        #print("rt_xx0: ", rt['x0'].shape)
        #print("z: ", rt['z'].shape)
        return (rt)



    def refineList(self, rt, state, pair, sel, intern):
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
        intern: int, iteration number
        """

        y = rt['y']
        x = rt['x']
        x0 = rt['x0']

        ut = np.unique(pair['TSSidx'])
        n = ut.shape[0]
        cellN = y.shape[0] / n
        stateN = x.shape[1]

        #print("some_n: ", n)
        #print("some_cellN: ", cellN)
        #print("some_stateN: ", stateN)

        #quit()

        tx = np.hstack((np.log2(x+0.001), np.log2(x0 + 0.001))) #n by stateN *2
        e = np.zeros((self.stateN * 2, self.cellN -1), dtype=np.float32) # the new betas in the refined list for x and x0s effect on transcription across cell cell_types
        f = np.zeros((self.utN, self.cellN-1), dtype=np.float32) # the new expression predictions

        index = 0
        for i in range(self.cellN):
            if i == self.lessone:
                continue
            mask = np.r_[np.arange(self.utN * min(i, self.lessone)),
                        np.arange(self.utN * (min(i,self.lessone)+1), n*max(i,self.lessone)),
                        np.arange(self.utN * (max(i, self.lessone) +1), y.shape[0])]
            # create a mask for all TSSs not in i or lessone cell type; leave out lessone cell type and cell type i
            r = self.lm(tx[mask, :], y[mask])
            #do linear regression with limited set of cell type data
            te = r['coeff']
            e[:,index] = te #uses index becasue it's not sure which i
            f[:,index] = np.dot(np.hstack((np.ones((self.utN, 1), dtype=np.flloat32),
                                            tx[(i*n):((i+1)*n),:])), te)[:,0]
            #create expression prediction for current cell type; for each TSS predict the left out cell type i expression
            index += 1

    def runRefine(self, rt, ss, pair):
        """
        rt: dict of (y: RNA expression for each TSS, z: normalized RNA, x: mean CRE state, x0: mean TSS state)
        cellN: number of cell types
        ss: float array of state masks, #positions by #states *#cell_types; mid position of CREs are mean state values across CRE position range
        pair: numpy object of TSS-CRE pairs (TSS, CRE positions, TSS index, predicted CRE position contributions)
        lessone: int, cell type to skip in leave-one-out
        """
        #print("PAIR SHAPE: ", pair.shape)
        #sys.exit()
        r = {'x':np.copy(rt['x'])}
        #k = rt['y'].shape[0]/ self.cellN k==self.utN
        #print(self.lessone*self.utN)
        #print((self.lessone+1)*self.utN)


        #print(rt['x0'].shape)
        #a00_x = rt['x0'][(self.lessone * self.utN):((self.lessone + 1)*self.utN),:]
        #print(a00_x.shape)

        #print(rt['y'].shape)
        #print(rt['y'][(self.lessone*self.utN):((self.lessone+1)*self.utN)].shape)

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

        for i in range(self.B):
            r = self.refineList(rt, ss, pair, sel, i)

    def run(self, chrom):
        self.rna = self.rna_all[chrom]
        self.state = self.state_all[chrom]
        self.tcre = self.tcre_all[chrom]
        self.pos = self.pos_all[chrom]
        self.tss = self.tss_all[chrom]

        m = np.mean(self.rna, axis=1)
        s = np.std(self.rna, axis=1, ddof =1)
        tt = np.where((m > self.m_thresh) & (s> self.s_thresh))[0]

        self.rna = self.rna[tt]
        self.tss = self.tss[tt]
        self.tssN = self.tss.shape[0]

        G = np.amax(self.state) + 1 #Determine number of states
        self.stateN = G
        t = np.bincount(self.state.ravel(order='F'), minlength=self.stateN) #Count occurence of each state
        #print("t_F: ", t)
        k = np.where(t == np.amax(t))[0][0] #Determine most prevelant state
        #print("k: ", k)
        l = max(np.amax(self.tss), np.amax(self.pos)) + 1 #Determine largest number (because of adding 1) of bins
        #print("l: ", l)
        self.l = l

        ss = np.full((l, self.cellN), k, dtype=np.int32) #Create array of state bins (including skipped regions and fill with most prevelant/ common state)
        ss[self.pos, :] = self.state #Fill in valid & known states
        #where_nonzero = np.where(ss != 0)
        #print(where_nonzero)
        #print(where_nonzero[0].shape)
        #print(where_nonzero[1].shape)
        #print("TCRE: ", self.tcre.shape)
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
        #print("tr: ", tr.shape)
        # rna is unnoramlized array, tr is normalized array
        #sse is unnormalized beta coefficient, ts is normalized beta coeff per genomic bin across cell types
        ts = ((sse - np.mean(sse, axis=1, keepdims=True))
                / ((np.std(sse, axis=1, keepdims=True, ddof =1) + 1e-5) * (self.cellN -1) ** 0.5 ))#.astype(np.int32)
        #print("ts: ", ts.shape)
        pair = []
        for i in range(self.tssN):
            a = np.arange(max(0, self.tss[i] - self.maxdist // 200), min(self.l, self.tss[i] + self.maxdist // 200 + 1)) #put maxdist into bins and get bin range of window maxdist around TSS

            rr = np.dot(tr[i:(i+1), :], ts[a,:].T).ravel(order='F') #predict state contribution to rna
            trr = np.dot(tr[i:(i+1), np.r_[np.arange(self.lessone), #predict position correlation to RNA leaving out one cell type specified by lessone
                                           np.arange(self.lessone +1, self.cellN)]],
                        (ts[a,:][:,np.r_[np.arange(self.lessone),
                                    np.arange(self.lessone +1, self.cellN)]]).T).ravel(order='F')
            #print("1st: ", tr[i:(i+1), np.r_[np.arange(self.lessone),
            #                               np.arange(self.lessone +1, self.cellN)]].shape)
            #print("2nd: ", ((ts[a,:][:,np.r_[np.arange(self.lessone),
            #            np.arange(self.lessone +1, self.cellN)]]).T).shape)

            #print("trr: ", trr.shape)
            #print("a: ", a.shape)
            t = np.where((trr >= self.cut) & (cre[a] >= 0))[0]
            #print("t: ", t.shape)
            #Find predicted position correlations that are above cutoff and contain a CRE within distal window around TSS
            #print("a_thing: ", a[t].shape)
            #print("cre_thing: ", cre[a[t]].shape)
            #print("unique_thing: ", np.unique(cre[a[t]]).shape)

            if t.shape[0] > 0:
                nt = []
                for j in np.unique(cre[a[t]]): #a[t] is global position, cre[a[t]] is the cre index at that position. Since a ccRE can overlap multiple bins, take unique
                    tt = t[np.where(cre[a[t]]==j)] #Identify each CRE position(s) of valid state in window; tt is local position of valid CRE
                    #print("tt: ",tt.shape)
                    #print("where thing: ",np.where(trr[tt]==np.amax(trr[tt]))[0].shape)
                    #if j == 60 and i == 1:
                    #    print trr[tt], a[tt]
                    nt.append(tt[np.where(trr[tt]==np.amax(trr[tt]))[0][0]]) #determine CRE position (if it spans multiple bins) with highest contribution/correlation and append to nt
                #print("tt shape: ", tt.shape)
                t = np.array(nt) #refining t to be the single position that has maximum correlation for each ccRE
            t = np.unique(t)
            ttss = np.where(a == self.tss[i])[0]
            if not np.isin(ttss, t)[0]:
                t = np.r_[ttss, t]
            else:
                where = np.where(t == ttss)[0][0]
                t = np.r_[ttss, t[:where], t[(where + 1):]]
            #if i == 299:
            #    c = cre[438785]
            #    j=np.where(cre[a] == c)[0]
            #    print("Cut: %s  C: %i  A: %s  TRR: %s. T: %s" % (j+1, c+1, a[j]+1, trr[j], t+1))

            #t = np.unique(np.r_[np.where(a == self.tss[i])[0], t]) #appending the TSS bin if it's not there; make position list only unique entries, including TSS bin
            #print("t shape: ",t.shape)
            for j in range(t.shape[0]):
                pair.append((self.tss[i], a[t[j]], i, rr[t[j]], cre[a[t[j]]])) #TSS, global cre bin position, tss index, predicted correlation without leave one out

        #print("a: ", a.shape)
        pair = np.array(pair, dtype=np.dtype([('TSS', np.int32), ('CRE', np.int32), #convert pair into a numpy array
                                                ('TSSidx', np.int32), ('rr', np.float32),
                                                ('CREidx', np.int32)]))
        np.savetxt('my_pair_TSS.txt', pair['TSS'])
        np.savetxt('my_pair_CRE.txt', pair['CRE'])
        np.savetxt('my_pair_TSSidx.txt', pair['TSSidx'])
        np.savetxt('my_pair_rr.txt', pair['rr'])
        self.pair = pair
        print('pair: ', self.pair.shape[0])
        #output = open('tmp_%s_py.txt' % chrom, 'a')
        #self.pair['TSS'] += 1
        #self.pair['CRE'] += 1
        #self.pair['TSSidx'] += 1
        #self.pair['CREidx'] += 1
        #for line in self.pair:
        #    print >> output, "%s %i %i %i %0.15f %i" % (chrom, line[0], line[1], line[2], line[3], line[4])
        #output.close()
        #return
        kk = np.zeros(self.stateN, dtype=np.int32) #create an array of size #ofstates
        #print("kk: ", kk.shape)
        kk[k] = 1 #mark most prevelant state
        #ss = np.tile(np.repeat(kk.reshape(1, -1), repeats=self.l, axis=0), self.cellN) #make a state mask array; horizontal stacking; 2d array
        ss = np.tile(np.repeat(kk.reshape(1, -1), repeats=self.l, axis=0), (1, self.cellN)).astype(np.float32)
        #print("ss: ", ss.shape)
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

        #print("ss: ", ss.shape)
        rt = self.genereg(self.rna, ss, self.pair, prior=None)
        r = self.runRefine(rt, ss, self.pair)

    def drive(self):
        """
        """

#statepref = '/Users/msauria/projects/VISION/Gene_Prediction/Data/pknorm_2_16lim_ref1mo_0424_lesshet.state'
#exp_file = '/Users/msauria/projects/VISION/Gene_Prediction/Data/rnaTPM.txt'
#cre_file = '/Users/msauria/projects/VISION/Gene_Prediction/Data/vision_cres.txt'
#atacsig_file = '/home/kweave23/VISION_regression/their_stuff/vision_cres.mat.atacsig.txt'
#state_by_chr_file = '/home/kweave23/VISION_regression/state_and_pos_by_chr.npz'
chrom = sys.argv[1]
#for i in range(1, 20):
#    chroms.append('chr%i' % i)
#chroms.append('chrX')

statepref = '/home/kweave23/VISION_regression/their_stuff/pknorm_2_16lim_ref1mo_0424_lesshet.state'
exp_file = '/home/kweave23/VISION_regression/their_stuff/rnaTPM.txt'
cre_file = '/home/kweave23/VISION_regression/their_stuff/vision_cres.txt'
test1 = regress_gene_cre(statepref, exp_file, cre_file, -4, 2)
test1.run(chrom)
