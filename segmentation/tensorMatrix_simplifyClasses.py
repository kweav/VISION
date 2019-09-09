#!/usr/bin/env python3

import sys
import numpy as np
import argparse as ap
import subprocess #use check_output() since bedtools sends to stdout


"""
Usage: ./tensorMatrix_simplifyClasses.py --outfile preprocessed_matrices_sc3.npz --IDEAScalls /project/vision/Data/IDEAS/segmentation_kw/*getfa.bed --RNAseq /project/vision/Data/RNAseq/scriptseq3.v3.kw2_sorted.tab  --ATACseq /project/vision/Data/ccRE/VISIONmusHem_ccREs_wct_filterkw.bed
"""
"""allows for -h for this file and for specific files to be given for annotation type"""
parser = ap.ArgumentParser(description='Preprocessing of data for machine learnings')
parser.add_argument('--IDEAScalls', action='store', nargs='+', type=str, required = True, help='List of files with IDEAS calls')
parser.add_argument('--RNAseq', action='store', nargs='+', type=str, required = True, help="File with RNAseq data") #
parser.add_argument('--ATACseq', action='store', nargs='+', type=str, required = True, help="File with ATACseq data")
parser.add_argument('--outfile', action='store', nargs=1, type=str, required = True, help="Name of file to save Annotated/matched matrices")
args=parser.parse_args()

fileList1 = args.IDEAScalls #IDEAScalls files
file2 = args.RNAseq[0] #RNAseq files
file3 = args.ATACseq[0] #ATACseq files
outfile = args.outfile[0]

print("First print statement")
"""annotating data"""
passedSequences=0
passedSequencesGenomeLength=0
numTotalGenomeLength = 0
numTotal = 0
cell_types =[]
loc={} #dictionary key: cellType, chrom, start, end; value: [ideasLabel, sequence]
# newLabel={17:0, 1:0, 25:0, 14:0, 8:0,
#           5:1, 4:1,
#           20:2, 3:2, 16:2, 22:2, 2:2,
#           9:3, 7:3 ,26:3, 13:3,
#           11:4, 12:4, 24:4, 15:4, 10:4, 23:4, 19:4, 21:4, 18:4, 6:4}
# newLabel2={17:0, 1:0, 25:0, 14:0, 8:0,
#           5:1, 4:1,
#           20:2, 3:2, 16:2, 22:2,
#           9:3, 7:3 ,26:3, 13:3,
#           11:4, 24:4, 15:4, 10:4, 21:4, 18:4,
#           12:5, 6:5, 23:5, 19:5,
#           2:6}
newLabel3={17:0, 1:0, 25:0, 14:0, 8:0,
          5:1, 4:1,
          20:2, 3:2, 16:2, 22:2,
          9:3, 7:3 ,26:3, 13:3,
          11:4, 24:4, 15:4, 10:4, 21:4, 18:4, 12:4, 6:4, 23:4, 19:4,
          2:5}


for file in fileList1: #IDEAScalls ##filter out if label == 0
    lhs,rhs=file.split("Seg") #file names are in the format ideasVisionV20p8Seg$CELLTYPEgetfa.bed -> ideasVisionV20p8 $CELLTYPEgetfa.bed
    cellTypeI, rhs2 = rhs.split("getfa") # -> $CELLTYPE .bed
    cell_types.append(cellTypeI)
    for line in open(file):
        fields=line.strip("\r\n").split("\t")
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        ideasLabel = int(fields[3])
        numTotalGenomeLength += (end-start)
        numTotal += 1
        if ideasLabel == 0:
            passedSequencesGenomeLength += (end-start)
            passedSequences += 1
            pass
        else:
            ideasLabel = newLabel3[ideasLabel]
            #ideasLabel = newLabel2[ideasLabel]
            if end - start <= 1000:
                loc[cellTypeI,chrom,start,end] = [ideasLabel, fields[9]] #ideaslabel=fields[3], sequence=fields[9]

print("percent passed, genomic length: ",passedSequencesGenomeLength/numTotalGenomeLength*100)
print("percent passed, count: ", passedSequences/numTotal*100)

print("IDEAS sequences added to dictionary")

"""annotating with RNAseq - file ex: chr    geneStart    geneEnd  Lsk=0.0;Cmp=72.0;Mep=0.0;G1e=180.0;Er4=0.0;Cfue=0.0;Eryad=0.0;Cfum=0.0;Imk=0.0;Gmp=0.0;Mon=0.0;Neu=0.0;"""
for cell_type in cell_types:
    bedtools_closest_out = subprocess.check_output("bedtools closest -k 3 -a /project/vision/Data/IDEAS/segmentation_kw/ideasVisionV20p8Seg{}getfa.bed -b {} -t last".format(cell_type, file2), shell=True).decode("utf-8").splitlines()
    for line in bedtools_closest_out: #3 lines correspond to the same location in the IDEAS file
        fields=line.strip("\r\n;").split("\t")
        chrom, startL, endL, startf, endf = fields[0], int(fields[1]), int(fields[2]), int(fields[11]), int(fields[12])
        ideasLabel = int(fields[3])
        if ideasLabel == 0:
            pass
        elif endL-startL > 1000:
            pass
        else:
            cellTypeIndex = fields[13].split(";")
            tpm = dict([x.split('=') for x in cellTypeIndex])[cell_type]
            loc[cell_type,chrom,startL,endL].append(float(tpm))

print("RNAseq values added to dictionary")


#cellTypesOI = ['Lsk', 'Cmp', 'Mep','G1e','Er4','Cfue','Eryad','Cfum','Imk','Gmp','Mon','Neu']

"""annotating with ATACseq - file ex: chr    geneStart    geneEnd    Lsk=0;Hpc7=1;Cmp=0;Mep=0;G1e=0;Er4=1;Cfue=1;Eryad=0;Eryfl=0;Cfum=1;Imk=1;Gmp=1;Mono=1;Neu=1;Nk=0;B=0;Tcd4=1;Tcd8=1;"""
ATACseqContainment = 0.5 #minimum containment or sequence overlap for the annotation to be added to the loc dictionary list for that genome location
for cell_type in cell_types:
    bedtools_out = subprocess.check_output("bedtools intersect -loj -a /project/vision/Data/IDEAS/segmentation_kw/ideasVisionV20p8Seg{}getfa.bed -b {}".format(cell_type, file3), shell=True).decode("utf-8").splitlines()
    for line in bedtools_out:
        fields=line.strip("\r\n;").split("\t")
        chrom, startL, endL, startf, endf = fields[0], int(fields[1]), int(fields[2]), int(fields[11]), int(fields[12])
        ideasLabel = int(fields[3])
        if ideasLabel == 0:
            pass
        elif endL - startL > 1000:
            pass
        elif fields[10] == ".":
            loc[cell_type, chrom, startL, endL].append(0) #no containment/or overlap between IDEAS and ATAC
        else:
            cellTypeIndex = fields[13].split(";")[:-1]
            aonab = dict([x.split('=') for x in cellTypeIndex])[cell_type]
            containment = (min(endL, endf) - max(startL, startf))/(endL-startL)
            if containment >= ATACseqContainment:
                loc[cell_type,chrom,startL,endL].append(int(aonab))
            else:
                loc[cell_type,chrom,startL,endL].append(0) #below specified containment level

print("ATACseq values added to dictionary")

"""separating annotated data into matched arrays"""
cellTypeIndex = []
labels = []
sequences = []
RNA_seq = []
ATAC_seq = []

cellType2Index = {
    "Cfue":0,
    "Cfum":1,
    "Cmp":2,
    "Er4":3,
    "Eryad":4,
    "G1e":5,
    "Gmp":6,
    "Imk":7,
    "Lsk":8,
    "Mep":9,
    "Mon":10,
    "Neu":11,
    "B":12,
    "Cd4":13,
    "Cd8":14,
    "Clp":15,
    "Eryfl":16,
    "Hpc7":17,
    "Mk":18,
    "Nk":19,
}

sequenceDict={
    "a":0,
    "c":1,
    "g":2,
    "t":3,
    'n':[0,1,2,3],
}

#loc - key:cellType,chrom,start,end
#loc - value=[label, sequence, RNAseq1, RNAseq2, RNAseq3, ATACseq1, ATACseq2,...,ATACseqi] where i is unknown
for (cellType, chrom, start, end) in loc:
    cellTypeIndex.append(cellType2Index.setdefault(cellType, 20))
    labels.append(loc[cellType, chrom, start, end][0])
    RNA_seq.append(loc[cellType,chrom,start,end][2:5]) #Append a list of 3
    if 1 in loc[cellType,chrom,start,end][5:]: #necessity for only 1 for the whole thing to be 1
        ATAC_seq.append(1)
    else:
        ATAC_seq.append(0)
    sequence = loc[cellType,chrom,start,end][1]
    sequenceTensorial = np.zeros((4,1000)) #2D array of shape (4,10000)

    for i,n in enumerate(sequence.lower()):
        sequenceTensorial[sequenceDict[n],i]=1


    sequenceTensorial[:,len(sequence)+1:1000] = np.nan #for any column values that are blank because the sequence isn't 10000 nucleotides, NaN

    sequences.append(sequenceTensorial) #append the 2D array to the list

print("Lists made")

cellTypeIndex = np.array(cellTypeIndex, dtype=np.intp)
labels = np.array(labels, dtype=np.intp)
sequences = np.array(sequences) #make the list of 2D arrays a 3D array
RNA_seq = np.array(RNA_seq, dtype=np.float_)
ATAC_seq = np.array(ATAC_seq, dtype=np.bool_)

print("Converted to arrays")

f = open(outfile, 'wb')
np.savez(f, cellTypeIndex = cellTypeIndex, labels = labels, sequences=sequences, RNA_seq = RNA_seq, ATAC_seq = ATAC_seq)
f.close()

print("File should be saved")
