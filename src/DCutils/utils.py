
import random
import math
from subprocess import check_output
import subprocess

def findIndels(seq):
    refPos = seq.reference_start
    readPos = 0
    indels = {}
    for cigar in seq.cigartuples:
        if cigar[0] == 0:
            refPos += cigar[1]
            readPos += cigar[1]
        if cigar[0] in [3,4]:
            readPos += cigar[1]
        if cigar[0] == 1:
            pos = refPos
            sequence = seq.query_sequence[readPos:readPos + cigar[1]]
            quals = seq.query_qualities[readPos:readPos + cigar[1]]
            indels.update({str(pos)+'I':[sequence,quals,readPos]})
            readPos += cigar[1]
        if cigar[0] == 2:
            pos = refPos
            indels.update({str(pos)+'D':[cigar[1],readPos]})
            refPos += cigar[1]  
    return indels      
            

def genotypeProbSnv(ObservedBases,quals,baseProb,ampErr):
    PO = 0
    while PO == 0:
        PONs = [0,0,0,0]
        for mm,base in enumerate(['A','T','C','G']):
            POB = 1
            for nn,ObservedBase in enumerate(ObservedBases):
                # Calculate P( Observation | base)
                seqErr = 10**(-quals[nn]/10)
                if base != ObservedBase:
                    Pnow = seqErr + ampErr - 2 * seqErr * ampErr
                else:
                    Pnow = 1 - (seqErr + ampErr - 2 * seqErr * ampErr)
                POB *= Pnow
            PONs[mm] = POB * baseProb[mm]
        PO = sum(PONs)
        if PO == 0:
            sampleIndex = random.sample(range(len(ObservedBases)),math.ceil(len(ObservedBases)/2))
            ObservedBases = [ObservedBases[ii] for ii in sampleIndex]
            quals = [quals[ii] for ii in sampleIndex]
            if len(sampleIndex) == 1: return [0.25,0.25,0.25,0.25]
    return [pon/PO for pon in PONs]

def genotypeProbIns(sr,ur,insBQ,indelPrior,ampErr,skipRate = 10**-(3.5)):
    PO = 0
    while PO == 0:
        PONs = [0,0]
        POBAlt = 1
        POBRef = 1
        for nn,bq in enumerate(insBQ):
            # Calculate P( Observation | base)
            seqErr = 10**(-bq/10)
            PnowAlt = 1 - (seqErr + ampErr - 2 * seqErr * ampErr)
            PnowRef = seqErr + ampErr - 2 * seqErr * ampErr
            POBAlt *= PnowAlt
            POBRef *= PnowRef
        for nn in range(ur):
            seqErr = skipRate
            PnowAlt = seqErr + ampErr - 2 * seqErr * ampErr
            PnowRef = 1 - (seqErr + ampErr - 2 * seqErr * ampErr)
            POBAlt *= PnowAlt
            POBRef *= PnowRef     
        PONs[0] = POBAlt * indelPrior[0]
        PONs[1] = POBRef * indelPrior[1]
        PO = sum(PONs)
        if PO == 0:
            sampleIndex = random.sample(range(sr+ur),math.ceil((sr+ur)/2))
            boolReads = [(nn < sr) for nn in range(sr+ur)]
            boolReadsSub = [boolReads[ii] for ii in sampleIndex]
            newsr = boolReadsSub.count(True)
            newur = len(sampleIndex) - newsr
            sr = newsr
            ur = newur
            #quals = [quals[ii] for ii in sampleIndex]
            if len(sampleIndex) == 1: return [0.5,0.5]
    return [pon/PO for pon in PONs]

def genotypeProbDel(sr,ur,delLen,indelPrior,ampErr,seqErrpb = 35):
    PO = 0
    while PO == 0:
        PONs = [0,0]
        POBAlt = 1
        POBRef = 1
        seqErr = 10**(-seqErrpb*delLen/10)
        for nn in range(sr):
            # Calculate P( Observation | base)
            PnowAlt = 1 - (seqErr + ampErr - 2 * seqErr * ampErr)
            PnowRef = seqErr + ampErr - 2 * seqErr * ampErr
            POBAlt *= PnowAlt
            POBRef *= PnowRef
        for nn in range(ur):
            #seqErr = skipRate
            PnowAlt = seqErr + ampErr - 2 * seqErr * ampErr
            PnowRef = 1 - (seqErr + ampErr - 2 * seqErr * ampErr)
            POBAlt *= PnowAlt
            POBRef *= PnowRef     
        PONs[0] = POBAlt * indelPrior[0]
        PONs[1] = POBRef * indelPrior[1]
        PO = sum(PONs)
        if PO == 0:
            sampleIndex = random.sample(range(sr+ur),math.ceil((sr+ur)/2))
            boolReads = [(nn < sr) for nn in range(sr+ur)]
            boolReadsSub = [boolReads[ii] for ii in sampleIndex]
            newsr = boolReadsSub.count(True)
            newur = len(sampleIndex) - newsr
            sr = newsr
            ur = newur
            #quals = [quals[ii] for ii in sampleIndex]
            if len(sampleIndex) == 1: return [0.5,0.5]
    return [pon/PO for pon in PONs]
"""
def extractDepthSnv(bam,chrom,pos,ref,alt):
    allele = []
    for pileupcolumn in bam.pileup(chrom,pos-1,pos,ignore_orphans=False,min_base_quality=0):
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                #if pileupread.
                if not pileupread.is_del and not pileupread.is_refskip and not pileupread.alignment.is_secondary and not pileupread.alignment.is_supplementary:
                    allele.append(pileupread.alignment.query_sequence[pileupread.query_position])
            break
    alt_count = allele.count(alt)
    ref_count = allele.count(ref)
    alts_count = len(allele) - ref_count
    depth = len(allele)
    return alt_count,alts_count,ref_count,depth
"""

def extractDepthSnv(bam,chrom,pos,ref,alt):
    pileupLine = check_output(["samtools","mpileup","-Q","1",bam,"-r",chrom+":"+str(pos)+"-"+str(pos)],stderr=subprocess.DEVNULL).decode('ascii')
    infos = pileupLine.strip('\n').split('\t')
    if len(infos) <= 1: return 0,0,0
    depth = int(infos[3])
    alleles = infos[4].upper()
    altCount = alleles.count(alt)
    refCount = alleles.count(ref)
    return altCount,refCount,depth

def extractDepthIndel(bam,chrom,pos,ref,alt):
    if len(ref) == 1:
        indelPos = str(pos)+'I'
    else:
        indelPos = str(pos)+'D'
    sr = 0
    ur = 0
    for read in bam.fetch(chrom,pos-1,pos):
        matchFlag = 0
        if read.mapping_quality <= 30: continue
        if 'I' in read.cigarstring or 'D' in read.cigarstring:
            indels = findIndels(read)
            if indels.get(indelPos):
                matchFlag = 1
        if matchFlag == 1:
            sr += 1
        else: ur += 1
    return sr,ur,sr+ur



def createVcfStrings(chromDict,infoDict,formatDict,filterDict,recs):
    lines = ["##fileformat=VCFv4.2"]
    for filterr in filterDict.keys():
        lines.append("##FILTER=<ID={},Description=\"{}\">".format(filterr,filterDict[filterr]))
    for info in infoDict.keys():
        lines.append("##INFO=<ID={},Number={},Type={},Description=\"{}\">".format(info,infoDict[info][0],infoDict[info][1],infoDict[info][2]))
    for form in formatDict.keys():
        lines.append("##FORMAT=<ID={},Number={},Type={},Description=\"{}\">".format(form,formatDict[form][0],formatDict[form][1],formatDict[form][2]))
    for chrom in chromDict.keys():
        lines.append("##contig=<ID={},length={}>".format(chrom,chromDict[chrom]))
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL")
    for rec in recs:
        chrom = rec['chrom']
        pos = rec['pos']
        alt = rec['alt']
        ref = rec['ref']
        infos = "."
        formats = rec['formats']
        samples = rec["samples"]
        lineEntries = [chrom,str(pos),'.',ref,alt,".","PASS",".",':'.join(formats),":".join([str(s) for s in samples[0]]),":".join([str(s) for s in samples[1]])]
        lines.append('\t'.join(lineEntries))
    return '\n'.join(lines)+"\n"




