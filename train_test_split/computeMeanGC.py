#!/usr/bin/env python3

import fasta
#from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
chromosomes = ['1', '2', '3', '4', '5','6', '7', '8', '9', '10', '11','12', '13', '14', '15', '16','17','18', '19', 'X', 'Y']

def computeGC(sequence):
    Gs = sequence.casefold().count('G'.casefold())
    Cs = sequence.casefold().count('C'.casefold())
    if Gs==0 and Cs==0:
        GCcontent = 0
    else:
        length = len(sequence)
        GCcontent = (Gs+Cs)/length
    return GCcontent

'''sliding windows''' #1 Mbp windows; slide by 500bp

for chromosome in chromosomes:
    file = '/Users/kateweaver/mm10_genome/chr{}.fa'.format(chromosome)
    #file = '/home-3/kweave23@jhu.edu/work/users/kweave23/mm10_genome/chr{}.fa'.format(chromosome)
    reader = fasta.FASTAReader(open(file))
    for ident, sequence in reader:
        window = 0
        slides = 0
        seqLen = len(sequence)
        gcList = []
        gcMeans = []
        starts = []
        for i in range(0, seqLen-1000000, 500):
            gc = computeGC(sequence[i:i+1000000])
            gcList.append(gc)
            starts.append(i)
            window += gc
            slides += 1
            if (i+1000000)%10000000 == 0 and i != 0:
                gcMean = window/slides
                gcMeans.append(gcMean)
                print(ident, gcMean, i, slides)
                window = 0
                slides = 0
        #How to handle final segment of chromosome
        gc=computeGC(sequence[list(range(0,seqLen-10000000,500))[-1]:seqLen])
        window += gc
        slides += 1
        gcMean = window/slides
        gcMeans.append(gcMean)
        print(ident, gcMean, slides, 'final')


    fig, ax = plt.subplots()
    plt.plot(starts, gcList)
    plt.scatter(np.arange(len(gcMeans))*10000000, gcMeans)
    plt.suptitle('1 Mbp Windows (sliding 500 bp) chr {}'.format(chromosome))
    ax.set_ylabel('GC content')
    ax.set_xlabel('Nucleotide Start Position')
    fig.savefig('slidingWindow_{}.png'.format(chromosome))
    plt.close(fig)
