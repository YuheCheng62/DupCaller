import argparse 
from DCutils.ReadGroup import ReadGroup
from pysam import AlignmentFile as BAM
from pysam import VariantFile as VCF
from pysam import TabixFile as BED
import pandas as pd
import numpy as np
from multiprocessing import Pool
from DCutils.utils import extractDepthSnv,extractDepthIndel
from DCutils.utils import createVcfStrings
import pysam


"""
def extractStrand(bam,chrom,pos,ref,alt):
    fRead = 0
    rRead = 0
    for pileupcolumn in bam.pileup(chrom,pos-1,pos,ignore_orphans=False):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if pileupread.alignment.is_reverse:
                    rRead += 1
                else:
                    fRead += 1
            break  
    return rRead,fRead
    
    for pileupcolumn in bam.pileup(chrom,pos-1,pos,ignore_orphans=False):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                #if pileupread.
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:
                    allele.append(pileupread.alignment.query_sequence[pileupread.query_position])
    alt_count = allele.count(alt)
    ref_count = allele.count(ref)
    alts_count = len(allele) - ref_count
    depth = len(allele)
    return alt_count,alts_count,ref_count,depth
"""


def filterOutput(params,nn):
    inputPrefix = params['output'] + '_' + str(nn)
    tBam = 'tmp/' + params['output'] + '_' + str(nn) + '_dedupped.bam'
    nBam = params['normalBam'] 
    noise = params['noise']
    germline = params['germline'] 
    trim5 = params['trim5']
    trim3 = params['trim3']
    trim5DBS = params['trim5DBS']
    trim3DBS = params['trim3DBS']
    minNdepth = params['minNdepth']
    maxAltCount = params['maxAltCount']
    readLength = params['readLength']
 
    muts = []
    print("Process "+str(nn)+": Initiated filtering jobs")
    normalBam = BAM(nBam,'rb')
    tumorBam = BAM(tBam,'rb')    
    mutDf = pd.read_csv("tmp/"+inputPrefix+"_muttable.tsv",sep='\t',\
        names = ["Chromosome","Position","Ref","Alt","AFP","TFP",'CFP',"GFP", "IDFP",\
            "AFC","TFC","CFC","GFC","IDFC", \
            "ARP","TRP",'CRP',"GRP","IDRP", \
            "ARC","TRC","CRC","GRC","IDRC", \
                "Count","gnomad","readStart",\
                "readPosition","barcode1","barcode2",\
                "meanNM","meanAS","meanXS","readPosition3p"],dtype={'Chromosome':str})
    germlineVcf =VCF(germline,'rb')
    mutDf = mutDf.drop_duplicates(subset=['Chromosome','Position'],ignore_index=True)
    if len(mutDf.index) == 0:
        return []
    """
    with open("tmp/tmp.bed",'w') as tmpBed:
        for row in mutDf.index:
            chrom = mutDf["Chromosome"]
            pos = mutDf["Position"]
            tmpBed.write(chrom,str(pos),str(pos))

    tPileupLines = check_output(["samtools","mpileup","-Q","0",tbam,"-R","tmp/tmp.bed"],stderr=subprocess.DEVNULL).decode('ascii').split('\n')[:-1]
    tDepthDict = {}
    for line in pileupLines:
        infos = pileupLine.strip('\n').split('\t')
        chrom = infos[0]
        pos = infos[1]
        depth = int(infos[3])
        alleles = infos[4].upper()
        altCount = alleles.count(alt)
        refCount = alleles.count(ref)
        depthDict.update({chrom+":"+"pos":[altCount,refCount,depth]})
    return altCount,refCount,depth
    """
    currentStart = mutDf.loc[0,"readStart"]
    currentBarcode = mutDf.loc[0,"barcode1"]+ mutDf.loc[0,"barcode2"]
    currentCount = 0
    currentMuts = []
    currentRecNum = 0
    for row in mutDf.index:
        currentRecNum += 1
        chrom = str(mutDf.loc[row,'Chromosome'])
        #print(mutDf.loc[row,:])
        dc = mutDf.loc[row,'Count']
        pos = int(mutDf.loc[row,'Position'])
        if currentRecNum%100000 == 0: 
            print("Process {} current progress {}:{}".format(nn,chrom,pos))
        meanNM = mutDf.loc[row,'meanNM']
        meanAS = mutDf.loc[row,'meanAS']
        meanXS = mutDf.loc[row,'meanXS']
        #if meanNM > 2: continue
        if meanAS - meanXS < 50: continue
        readPos = int(mutDf.loc[row,'readPosition'])
        readPos3p = int(mutDf.loc[row,'readPosition3p'])
        ref = mutDf.loc[row,'Ref']
        alt = mutDf.loc[row,'Alt']
        start = str(mutDf.loc[row,'readStart'])
        bc1 = mutDf.loc[row,'barcode1']
        bc2 = mutDf.loc[row,'barcode2']
        dupCount = mutDf.loc[row,'Count']
        #print(alt)
        if bc1+bc2 != currentBarcode or start != currentStart:
            currentBarcode = bc1 + bc2
            currentStart = start
            passFlag = 1
            if len(currentMuts) > 1:
                positions = [mut["pos"] for mut in currentMuts]
                for nnn in range(1,len(positions)):
                    if positions[nnn] - positions[nnn-1] != 1:
                        passFlag = 0
                        break
                    elif currentMuts[0]["samples"][0][6] < trim3DBS or \
                        currentMuts[0]["samples"][0][7] < trim5DBS:
                        passFlag = 0 
                        break
            if passFlag == 1:
                muts += currentMuts
            currentMuts = []
        #print(readPos,trim5,trim3)
        if readPos <= trim5 or readPos3p <= trim3: continue 
        #print("pass")
        maskFlag = 0
        if noise:
            mask = BED(noise)
            for rec in mask.fetch(chrom,pos-1,pos+1,parser=pysam.asBed()):
                maskFlag = 1
                break
        if maskFlag == 1: continue
        if alt == "*":
            probsF = mutDf.loc[row,['AFP','TFP','CFP','GFP']].apply(lambda x: float(x))
            probsR = mutDf.loc[row,['ARP','TRP','CRP','GRP']].apply(lambda x: float(x))
            altF = probsF.idxmax()[0]
            altR = probsR.idxmax()[0]
            #if meanNM > 3: continue
            if altF != altR: continue
            else: alt = altF

            if alt == ref: continue

            gnomadFilter = 0
            if mutDf.loc[row,'gnomad'] == 1:
                for rec in germlineVcf.fetch(chrom,pos-1,pos):
                    if alt in rec.alts:
                        if rec.info['AF'][rec.alts.index(alt)] >= 0.001:
                            gnomadFilter = 1
            if gnomadFilter == 1:
                continue
            #fRead,rRead = extractStrand(tumorBam,chrom,int(pos),ref,alt)
            #if fRead == 0 or rRead == 0: continue
            #print(normalBam,chrom,pos)
            na,nr,ndp = extractDepthSnv(nBam,chrom,int(pos),ref,alt)
            ta,tr,tdp = extractDepthSnv(tBam,chrom,int(pos),ref,alt)
            if na <= maxAltCount and ndp > minNdepth:
                mut = {"chrom":chrom,"pos":pos,"ref":ref,"alt":alt,"infos":{},"formats":['AC',"RC","DP","NM","AS","XS","RP5","RP3","DC"],"samples":[[ta,tr,tdp,meanNM,meanAS,meanXS,readPos,readPos3p,dupCount],[na,nr,ndp,0,0,0,0,1]],"start":start,"barcode":bc1+bc2}
                currentMuts.append(mut)
            
            """
        elif len(alt) == 1:
            gnomadFilter = 0
            if mutDf.loc[row,'gnomad'] == 1:
                for rec in germlineVcf.fetch(chrom,pos-1,pos):
                    if ref == rec.ref:
                        if rec.info['AF'][0] >= 0.01:
                            gnomadFilter = 1
            if gnomadFilter == 1:
                continue
            na,nr,ndp = extractDepth(normalBam,chrom,int(pos),ref,alt)
            if na == 0 and ndp > minNdepth:
                outTable.write('\t'.join(['Exnano','HepG2','.','GRCh38','SNP',chrom,str(pos),str(pos+len(ref)-1),ref,alt,'SOMATIC',str(na),str(nr),str(readPos)])+'\n')
            """
        else:
           # gnomadFilter = mutDf.loc[row,'gnomad']
            gnomadFilter = 0
            if mutDf.loc[row,'gnomad'] == 1:
                for rec in germlineVcf.fetch(chrom,pos-1,pos):
                    for altNow in rec.alts:
                        if (len(altNow) - len(rec.ref) == len(alt) - len(ref)) and rec.info['AF'][rec.alts.index(altNow)] >= 0.001:
                            gnomadFilter = 1
            if gnomadFilter == 1:
                continue
            if min(meanNM - len(ref),meanNM - len(alt)) > 2:
                continue
            if readPos3p <=15:
                continue
            na,nr,ndp = extractDepthIndel(normalBam,chrom,int(pos),ref,alt)
            ta,tr,tdp = extractDepthIndel(tumorBam,chrom,int(pos),ref,alt)
            #print(na,nr,ndp,ta,tr,tdp)
            if na == 0 and ndp > minNdepth:
                mut = {"chrom":chrom,"pos":pos,"ref":ref,"alt":alt,"infos":{},"formats":['AC',"RC","DP","NM","AS","XS","RP5","RP3","DC"],"samples":[[ta,tr,tdp,meanNM,meanAS,meanXS,readPos,readPos3p,dupCount],[na,nr,ndp,0,0,0,0,0,1]]}
                currentMuts.append(mut)
                #outTable.write('\t'.join(['Exnano','HepG2','.','GRCh38','SNP',chrom,str(pos),str(pos),ref,alt,'SOMATIC',str(na),str(nr),str(readPos)])+'\n')  
    passFlag = 1
    if len(currentMuts) > 1:
        positions = [mut["pos"] for mut in currentMuts]
        for nn in range(1,len(positions)):
            if positions[nn] - positions[nn-1] != 1:
                passFlag = 0
                break
    if passFlag == 1:
        muts += currentMuts
    currentMuts = []
    return muts
    print("Process "+str(nn)+": finished filtering jobs")          
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Call any mutations from a bam file')
    parser.add_argument('-i','--prefix',type=str,help='mutation table generated from callBam')
    parser.add_argument('-b','--tumorBam',type=str,help='deduplicated normal bam')   
    parser.add_argument('-n','--normalBam',type=str,help='deduplicated normal bam')
    parser.add_argument('-g','--germline',type=str,help='indexed germline af vcf')
    #parser.add_argument('-f','--reference',type=str,help='reference fasta file')
    parser.add_argument('-m','--noise',type=str,help='noise mask')
    parser.add_argument('-o','--output',type=str,help='prefix of the output files')
    parser.add_argument('-p','--threads',type=int,help='prefix of the output files')
    parser.add_argument('-t5','--trim5',type=int,help='ignore mutation if it is less than n bps from 5 prime',default=0)
    parser.add_argument('-t3','--trim3',type=int,help='ignore mutation if it is less than n bps from 3 prime',default=0)
    parser.add_argument('-t52','--trim5DBS',type=int,help='ignore mutation if it is less than n bps from 5 prime',default=30)
    parser.add_argument('-t32','--trim3DBS',type=int,help='ignore mutation if it is less than n bps from 3 prime',default=30)
    parser.add_argument('-d','--minNdepth',type=int,help='minumum coverage in normal for called variants',default=10)
    parser.add_argument('-r','--readLength',type=int,help='trimmed read length')
    parser.add_argument('-nad','--maxAltCount',type=int,default=0)
    #parser.add_argument('-re','--repeatBed',type=int,default=1)
    args = parser.parse_args()
    #tBam = BAM(args.tumorBam,'rb') 
    #contigs = tBam.references
    #print(contigs)
    #chromDict = {contig:tBam.get_reference_length(contig) for contig in contigs}
    #print(chromDict)
    #mutDf = mutDf.drop_duplicates(subset=['Chromosome','Position'],ignore_index=True)
    #somaticVcf = VCF(output+'.somatic.vcf','w')
    pool = Pool()
    callArguments = [(args.prefix+'_'+str(nn),args.tumorBam,args.normalBam,args.noise,args.germline,\
        args.trim5,args.trim3,args.trim5DBS,args.trim3DBS,args.minNdepth,args.maxAltCount,args.readLength,nn) for nn in range(args.threads)]
    muts = pool.starmap(filterOutput,callArguments)
    mutsAll = sum(muts,[])
    pool.close()
    pool.terminate()
    pool.join()
    """
    mutTableNames = ["tmp/"+args.prefix+'_'+str(nn)+"_muttable_filtered.tsv" for nn in range(args.threads)]
    with open(args.prefix+"_muttable_filtered.tsv",'w+') as merged:
        merged.write("\t".join(["Project","Sample","ID","Genome","mut_type",\
        "chrom","pos_start","pos_end","ref","alt","Type","nalt","nref","readPosition"])+'\n')
        for table in mutTableNames:
            with open(table,'r') as myTable:
                for line in myTable.readlines():
                    merged.write(line)
    """
    tBam = BAM(args.tumorBam,'rb')
    contigs = tBam.references
    #print(contigs)
    chromDict = {contig:tBam.get_reference_length(contig) for contig in contigs}
    infoDict = {}
    formatDict = {"AC":[1,"Integer","Count of alt allele"],"RC":[1,"Integer","Count of ref allele"],"DP":[1,"Integer","Depth at the location"],\
    "NM":[1,"Float","mean NM of alt read group"],"AS":[1,"Float","mean AS of the read group"],"XS":[1,"Float","mean XS of the read group"],\
    "RP5":[1,"Integer","read position"],"RP3":[1,"Integer","distance from 3p"],"DC":[1,"Integer","Number of reads in the duplex group"]}
    filterDict = {"PASS":"All filter Passed"}
    vcfLines = createVcfStrings(chromDict,infoDict,formatDict,filterDict,mutsAll)
    with open(args.output,'w') as vcf:
        vcf.write(vcfLines)