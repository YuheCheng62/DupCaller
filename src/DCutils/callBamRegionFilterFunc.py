import argparse 
from DCutils.ReadGroup import ReadGroup
from Bio import SeqIO
from pysam import AlignmentFile as BAM
from pysam import VariantFile as VCF
from pysam import index as indexBam
#from pysam import FastaFile as FASTA
import time
from DCutils.filterRawCall import filterOutput
from DCutils.utils import createVcfStrings

"""
def calculateMismatches(seqs,referenceSequence):
    mismatch = 0
    for seq in seqs:
        sequence = seq.query_sequence
        for ii in range(len(sequence)):
            if sequence[ii] != referenceSequence[ii]:
                mismatch += 1
                break
    return mismatch
"""


def bamIterateMultipleRegion(bamObject,regions):
    for region in regions:
        for rec in bamObject.fetch(*region):
            if len(region)  >= 2:
                if rec.reference_start < region[1]:
                    continue
            yield rec


def callBam(params,processNo,chunkSize):

    bam = params['tumorBam']
    regions = params['regions']
    germline = params['germline']
    reference = params['reference']
    minMapq = params["mapq"]
    mutRate = params["mutRate"]
    perror = params["perror"]
    pcut = params["pcutoff"]
    nn = processNo
    output = 'tmp/' + params["output"] + '_' + str(nn)
    
    print("Process" + str(processNo)+": Initiated")
    #print(regionStrings)
    #regionStringsOne = ' '.join(regionStrings)
    starttime = time.time()
    #starttime = time.time()
    tumorBam = BAM(bam,'rb')
    currentReadSet = []
    currentBc = ''
    currentStart = -1
    currentReadDict = {}
    #sameBcFlag = 1
    chromPast = 'startChrom'
    fasta = reference
    #fasta = FASTA(reference)
    #fasta = SeqIO.to_dict(SeqIO.parse(reference,'fasta'))
    germline = VCF(germline,'rb')
    dedupBam = BAM(output+'_dedupped.bam','wb',header=tumorBam.header)
    duplexBam = BAM(output+'_duplex.bam','wb',header=tumorBam.header)
    with open(output+'_muttable.tsv','w') as mutTable:
        mutTable.write('')
    outTable = open(output+'_muttable.tsv','a')
    recCount = 0
    #checkPoints = [nn * int(chunkSize/20) for nn in range(1,21)]
    #currentCheckPoint = 1
    currentCheckPoint = 100000

    passRead = [0] * 31
    passLength = [0] * 31
    for rec in bamIterateMultipleRegion(tumorBam,regions):
        #print(recCount)
        recCount += 1
        if recCount == currentCheckPoint:
            print("Process"+str(processNo)+": processed "+str(recCount)+" reads in "+\
                str((time.time()-starttime)/60)+" minutes" )
            currentCheckPoint += 100000
        if rec.mapping_quality <= minMapq or \
        rec.is_supplementary or \
        rec.is_secondary or \
        rec.is_duplicate or \
        not rec.is_proper_pair or \
        rec.is_qcfail: 
            continue
        #if rec.get_tag('AS') - rec.get_tag('XS') <= 50: continue
        if rec.cigartuples[0][0] == 4: continue
        if rec.cigarstring.count('I') + rec.cigarstring.count('D') >= 2: continue
        start = rec.reference_start
        bc = rec.query_name.split('_')[1]
        bcsplit = bc.split('+')
        bc1 = bcsplit[0]
        bc2 = bcsplit[1]
        #if currentBc =='': currentBc = bc
        if currentStart == -1: currentStart = start 
        if chromPast == 'startChrom':
            chromPast = tumorBam.get_reference_name(rec.reference_id)
        else:
            chromPast = chrom
        chrom = tumorBam.get_reference_name(rec.reference_id)
        #print(start,currentStart)
        if start == currentStart:
            #print(currentReadDict,1)
            if currentReadDict.get(bc1+'+'+bc2) != None:
                currentReadDict[bc1+'+'+bc2]["seqs"].append(rec)
                currentReadDict[bc1+'+'+bc2]["F1R2"] += 1
            elif currentReadDict.get(bc2+'+'+bc1) != None:
                currentReadDict[bc2+'+'+bc1]["seqs"].append(rec)
                currentReadDict[bc2+'+'+bc1]["F2R1"] += 1
            else:
                currentReadDict.update({bc:{"seqs":[rec],"F1R2":1,"F2R1":0}})

        else:
            for key in currentReadDict.keys():
                readSet = currentReadDict[key]["seqs"]
                setBc = key.split('+')
                setBc1 = setBc[0]
                setBc2 = setBc[1]
                F2R1 = currentReadDict[key]["F2R1"]
                F1R2 = currentReadDict[key]["F1R2"]
                if setBc1 != setBc2 and F2R1 > 0 and F1R2+F2R1 > 1:
                    NMs = [seq.get_tag('NM') for seq in readSet]
                    meanNM = sum(NMs)/len(NMs)
                    dsi = F1R2 + F2R1
                    if dsi <= 29:
                        passRead[dsi] += 1
                        passLength[dsi] += readSet[0].reference_length
                    else:
                        passRead[30] += 1
                        passLength[30] += readSet[0].reference_length

                    if meanNM >= 0.5:
                        RG = ReadGroup(readSet,setBc1,setBc2,chrom=chromPast)
                        caller = RG.call(fasta=fasta,germlineVcf=germline, somaticMutRate=mutRate, ampErr=perror, pcut=pcut)
                        for mut in caller:
                            #print(mut)
                            outTable.write('\t'.join(mut)+'\n')
                        dedupBam.write(RG.pivotSeq) 
                        duplexBam.write(RG.pivotSeq)
                    else:
                        dedupBam.write(readSet[0]) 
                        duplexBam.write(readSet[0])
                else:  dedupBam.write(readSet[0])         
            currentReadDict = {bc:{"seqs":[rec],"F1R2":1,"F2R1":0}}
            currentStart = start

    for key in currentReadDict.keys():
        readSet = currentReadDict[key]["seqs"]
        setBc = key.split('+')
        setBc1 = setBc[0]
        setBc2 = setBc[1]
        F2R1 = currentReadDict[key]["F2R1"]
        F1R2 = currentReadDict[key]["F1R2"]
        if setBc1 != setBc2 and F2R1 > 0 and F1R2+F2R1 > 1:
            referenceSeq = fasta[chromPast][readSet[0].reference_start:readSet[0].reference_start + readSet[0].query_length]
            NMs = [seq.get_tag('NM') for seq in readSet]
            meanNM = sum(NMs)/len(NMs)
            if dsi <= 29:
                passRead[dsi] += 1
                passLength[dsi] += readSet[0].reference_length
            else:
                passRead[30] += 1
                passLength[30] += readSet[0].reference_length
    
            if meanNM >= 0.5:
                RG = ReadGroup(readSet,setBc1,setBc2,chrom=chromPast)
                caller = RG.call(fasta=fasta,germlineVcf = germline, somaticMutRate = mutRate, ampErr = perror, pcut = pcut)
                for mut in caller:
                    outTable.write('\t'.join(mut)+'\n')
                dedupBam.write(RG.pivotSeq) 
            else:
                dedupBam.write(readSet[0]) 
        else:  dedupBam.write(readSet[0])  

    dedupBam.close()
    indexBam(output+'_dedupped.bam')
    print("Process" + str(processNo)+": Completed raw calling for "+str(recCount)+" reads,filtering results......")
    
    muts = filterOutput(params,nn)

    print("Process" + str(processNo)+": Completed")
    return passRead,passLength,muts
