import numpy as np
import math
import random
from DCutils.utils import findIndels,genotypeProbIns,genotypeProbSnv,genotypeProbDel
import array


class ReadGroup:

    def __init__(self,seqs,bc1,bc2,chrom):
        seqNameDict = {}
        self.seqs = seqs
        NMs = [seq.get_tag('NM') for seq in self.seqs]
        self.meanNM = sum(NMs)/len(NMs)
        ASs = [seq.get_tag('AS') for seq in self.seqs]
        self.meanAS = sum(ASs)/len(ASs)
        XSs = [seq.get_tag('XS') for seq in self.seqs]
        self.meanXS = sum(XSs)/len(XSs)
        self.readLengths = [seq.query_length for seq in self.seqs]
        self.readLength = max(self.readLengths)
        self.count = len(self.seqs)
        self.bc1 = bc1
        self.bc2 = bc2
        self.seqs = seqs
        cigars = [seq.query_sequence for seq in seqs]
        pivotIndex = cigars.index(max(set(cigars),key=cigars.count))
        self.pivotIndex = pivotIndex
        pivotSeq = seqs[pivotIndex]
        self.pivotSeq = pivotSeq
        #self.pivotName = pivotSeq.query_name
        self.cigar = pivotSeq.cigarstring
        pivotCigar = pivotSeq.cigartuples
        self._cigartuples = pivotCigar
        self.start = pivotSeq.reference_start
        self.chrom = chrom
        if self._cigartuples[-1][0] == 4:
            self.effLength = self.readLength - self._cigartuples[-1][1]
        else:
            self.effLength = self.readLength

    def _padSequence(self,sequence):
        return sequence + 'N'* (self.readLength-len(sequence))

    def _padQuality(self,qualities):
        #print(qualities,qualities[10])
        return qualities + array.array('B',[0] * (self.readLength-len(qualities)))
    
    def _padCigarArray(self,cigarArray):
        return cigarArray + ['S'] * (self.readLength-len(cigarArray))

    def call(self,fasta=False,germlineVcf=False,somaticMutRate = 3*10E-7,ampErr = 10E-5, pcut = 0.01):
        def cigarTuplesToArray(cigartuples):
            array = []
            cigarLists = [list(t) for t in cigartuples]
            while True:
                if cigarLists[0][1] == 0:
                    cigarLists.pop(0)
                    if len(cigarLists) == 0:
                        break
                convertList = ['M','I','D','N','S','H']
                if cigarLists[0][0] in [0,1,3,4,5]:
                    array += [convertList[cigarLists[0][0]]]
                    cigarLists[0][1] -= 1
                else:
                    array += [str(cigarLists[0][1]) + 'D']
                    cigarLists[0][1] = 0
            return array

        strands = ['F' if seq.query_name.split('_')[1].split('+')[0] == self.bc1 else 'R' for seq in self.seqs]   
        chrom = self.chrom

        if 'I' in self.cigar or 'D' in self.cigar:
            indels = [findIndels(seq) for seq in self.seqs]
            for indel in indels[self.pivotIndex].keys():
                indelPrior = [somaticMutRate/50,1-somaticMutRate/50]
                gnomadFlag = "0"
                sfr = 0 # supporting forward read
                ufr = 0 # unsupporting forward read
                srr = 0 # supporting reverse read
                urr = 0 # unsupporting reverse read
                """
                call insertions
                """
                if indel[-1] == "I":
                    FinsBQ = []
                    RinsBQ = []
                    currentPos = int(indel[:-1]) # 1-based position
                    insSeq = indels[self.pivotIndex][indel][0]
                    ref = fasta[self.chrom][currentPos-1].upper()
                    altNow = ref+insSeq
                    readPos = indels[self.pivotIndex][indel][2]
                    readPos3p = self.effLength - readPos
                    for nn,seq in enumerate(self.seqs):
                        hasIndel = 0
                        if indels[nn].get(indel):
                            hasIndel = 1
                            if indels[nn][indel][0] == insSeq:
                                if strands[nn] == 'F':
                                    sfr += 1
                                    FinsBQ.append(sum(indels[nn][indel][1]))
                                else:
                                    srr += 1
                                    RinsBQ.append(sum(indels[nn][indel][1]))   
                        elif hasIndel == 0:
                            if strands[nn] == 'F':
                                ufr += 1
                            else:
                                urr += 1
                    if chrom in germlineVcf.header.contigs:
                        for rec in germlineVcf.fetch(chrom,currentPos-1,currentPos):
                            gnomadFlag = "1"
                            #print("rec",rec.ref,"ref",ref,'I')
                            for nn,alt in enumerate(rec.alts):
                                if altNow == alt:
                                    indelPrior[0] = rec.info['AF'][nn]
                                    indelPrior[1] = 1 - indelPrior[0]
                    insFPosterior = genotypeProbIns(sfr,ufr,FinsBQ,indelPrior,ampErr/50)
                    insRPosterior = genotypeProbIns(srr,urr,RinsBQ,indelPrior,ampErr/50) 
                    if (insFPosterior[1]<= pcut and insRPosterior[1] <= 0.1) or (insRPosterior[1] <= pcut and insFPosterior[1] <= 0.1):
                        yield [chrom]+[str(currentPos)]+[ref,altNow]+ \
                        ["0","0","0","0",str(insFPosterior[0])] + ["0","0","0","0",str(sfr)]+ \
                        ["0","0","0","0",str(insRPosterior[0])] + ["0","0","0","0",str(srr)]+ \
                            [str(self.count),gnomadFlag,str(self.start),str(readPos),self.bc1,self.bc2,str(self.meanNM),str(self.meanAS),str(self.meanXS),str(readPos3p)]       

                """
                call deletions
                """
                if indel[-1] == "D":
                    currentPos = int(indel[:-1]) # 1-based position
                    delLen = indels[self.pivotIndex][indel][0]
                    ref = str(fasta[self.chrom][currentPos-1:currentPos+delLen].seq).upper()
                    altNow = ref[0]
                    readPos = indels[self.pivotIndex][indel][1]
                    readPos3p = self.effLength - readPos
                    for nn,seq in enumerate(self.seqs):
                        if indels[nn].get(indel):
                            if indels[nn][indel][0] == delLen:
                                if strands[nn] == 'F':
                                    sfr += 1
                                else:
                                    srr += 1 
                        else:
                            if strands[nn] == 'F':
                                ufr += 1
                            else:
                                urr += 1
                    if chrom in germlineVcf.header.contigs:
                        for rec in germlineVcf.fetch(chrom,currentPos-1,currentPos):
                            gnomadFlag = "1"
                            #print("rec",rec.ref,"ref",ref,'D')
                            if rec.ref == ref:
                                for nn,alt in enumerate(rec.alts):
                                    if altNow == alt:
                                        indelPrior[0] = rec.info['AF'][nn]
                                        indelPrior[1] = 1 - indelPrior[0]
                    delFPosterior = genotypeProbDel(sfr,ufr,delLen,indelPrior,ampErr/50)
                    delRPosterior = genotypeProbDel(srr,urr,delLen,indelPrior,ampErr/50) 
                    if (delFPosterior[1]<= pcut and delRPosterior[1] <= 0.1) or (delRPosterior[1] <= pcut and delFPosterior[1] <= 0.1):
                        yield [chrom]+[str(currentPos)]+[ref,altNow]+ \
                        ["0","0","0","0",str(delFPosterior[0])] + ["0","0","0","0",str(sfr)]+ \
                        ["0","0","0","0",str(delRPosterior[0])] + ["0","0","0","0",str(srr)]+ \
                            [str(self.count),gnomadFlag,str(self.start),str(readPos),self.bc1,self.bc2,str(self.meanNM),str(self.meanAS),str(self.meanXS),str(readPos3p)]               
 

        else:
            cPointers = [0]* self.count
            sPointers = [0]* self.count
            sequences = [self._padSequence(s.query_sequence) for s in self.seqs]
            quals = [self._padQuality(s.query_qualities) for s in self.seqs]
            cigars = [self._padCigarArray(cigarTuplesToArray(s.cigartuples)) for s in self.seqs]
            #sequences = [s.query_sequence.upper() for s in self.seqs]
            #quals = [s.query_qualities for s in self.seqs]
            #cigars = [cigarTuplesToArray(s.cigartuples) for s in self.seqs]
            #print(sequences,quals,cigars)
            completed = 0
            currentPos = self.start
            contigLen = len(fasta[self.chrom])
            while completed < self.count and currentPos < contigLen:
                observedFBases = []
                observedFQuals = []
                observedRBases = []
                observedRQuals = []
                ref = fasta[self.chrom][currentPos].upper()
                #print(ref)
                Dnums = []
                for nn in range(len(cPointers)):
                    if cPointers[nn] == len(cigars[nn]): continue
                    while cigars[nn][cPointers[nn]] == 'S' or cigars[nn][cPointers[nn]] == 'I':
                        cPointers[nn] += 1    
                        sPointers[nn] += 1  
                        if cPointers[nn] == len(cigars[nn]): 
                            completed += 1
                            break
                    if cPointers[nn] == len(cigars[nn]): continue
                        #if cigars[nn][cPointers[nn]] == 'I':
                            #while cigars[nn][cPointers[nn]] == 'I':
                                #sPointers[nn] += 1
                                #cPointers[nn] += 1
                    if cigars[nn][cPointers[nn]] == 'M':
                        if sequences[nn][sPointers[nn]] in ['A','T','C','G']:
                            if strands[nn] == 'F':
                                observedFBases.append(sequences[nn][sPointers[nn]])
                                observedFQuals.append(quals[nn][sPointers[nn]])
                            else:
                                observedRBases.append(sequences[nn][sPointers[nn]])
                                observedRQuals.append(quals[nn][sPointers[nn]])             
                        sPointers[nn] += 1
                        cPointers[nn] += 1
                        if cPointers[nn] == len(cigars[nn]): completed += 1

                    elif cigars[nn][cPointers[nn]][-1] == 'D':
                        delNum = int(cigars[nn][cPointers[nn]][:-1]) - 1
                        if delNum == 0: 
                            cPointers[nn] += 1
                        else:
                            cigars[nn][cPointers[nn]] = str(delNum) + 'D'
                            Dnums.append(delNum)
                currentPos += 1
                
                if len(Dnums) + completed >= self.count and Dnums:
                    DSkip = min(Dnums)
                    for nn in range(len(cPointers)):
                        if cPointers[nn] < len(cigars[nn]):
                            if cigars[nn][cPointers[nn]][-1] == 'D':
                                delNum = int(cigars[nn][cPointers[nn]][:-1]) - DSkip
                                if delNum == 0: 
                                    cPointers[nn] += 1
                                else:
                                    cigars[nn][cPointers[nn]] = str(delNum) + 'D'
                    currentPos += DSkip   
            

                readPos = sPointers[0]
                readPos3p = self.effLength - readPos
                if ref not in ['A','T','C','G']: continue
                if any([b != ref for b in observedFBases]) and any([b != ref for b in observedRBases]):
                    baseProb = [0,0,0,0]
                    gnomadFlag = "0"
                    if chrom in germlineVcf.header.contigs:
                        for rec in germlineVcf.fetch(chrom,currentPos-1,currentPos):
                            #print(rec.ref,ref)
                            gnomadFlag = "1"
                            if len(rec.ref) != 1: continue
                            AF = rec.info['AF']
                            for nn,alt in enumerate(rec.alts):
                                if len(alt) != 1: continue
                                if alt == 'A': baseProb[0] = AF[nn]
                                elif alt == 'T': baseProb[1] = AF[nn]
                                elif alt == 'C': baseProb[2] = AF[nn]
                                elif alt == 'G': baseProb[3] = AF[nn]
                    for nn,Base in enumerate(['A','T','C','G']):
                        if baseProb[nn] <= somaticMutRate and Base != ref:
                            baseProb[nn] = somaticMutRate
                    refIndex = ['A','T','C','G'].index(ref)
                    baseProb[refIndex] = 1 - sum(baseProb)
                    if observedFBases:
                        baseFPosterior = genotypeProbSnv(ObservedBases=observedFBases,quals=observedFQuals,baseProb=baseProb,ampErr=ampErr)
                    else:
                        baseFPosterior = baseProb
                    if observedRBases:
                        baseRPosterior = genotypeProbSnv(ObservedBases=observedRBases,quals=observedRQuals,baseProb=baseProb,ampErr=ampErr)
                    else:
                        baseRPosterior = baseProb
                    if (baseFPosterior[refIndex] <= pcut and baseRPosterior[refIndex] <= 0.1) or (baseRPosterior[refIndex] <= pcut and baseFPosterior[refIndex] <= 0.1):
                        yield [chrom]+[str(currentPos)]+[ref,"*"]+ \
                        ["%7f" % p for p in baseFPosterior] + ['0.0000000'] + [str(observedFBases.count(b)) for b in ['A','T','C','G']]+ ['0'] +\
                        ["%7f" % p for p in baseRPosterior] + ['0.0000000'] + [str(observedRBases.count(b)) for b in ['A','T','C','G']]+ ['0'] +\
                            [str(self.count),gnomadFlag,str(self.start),str(readPos),self.bc1,self.bc2,str(self.meanNM),str(self.meanAS),str(self.meanXS),str(readPos3p)]

