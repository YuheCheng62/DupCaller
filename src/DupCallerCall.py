from DCutils.splitBamRegions import splitBamRegions
from DCutils.callBamRegionFilterFunc import callBam
from DCutils.utils import createVcfStrings

import argparse
from Bio import SeqIO
import os
from multiprocessing import Pool
from pysam import AlignmentFile as BAM
import pysam
import time


if __name__ == "__main__":

    """
    Parse Arguments
    """
    parser = argparse.ArgumentParser(description='Call any mutations from a bam file')
    parser.add_argument('-b','--bam',type=str,help='bam file')
    parser.add_argument('-g','--germline',type=str,help='indexed germline af vcf')
    parser.add_argument('-f','--reference',type=str,help='reference fasta file')
    parser.add_argument('-o','--output',type=str,help='prefix of the output files')
    parser.add_argument('-r','--regions',nargs='+',type=str,help='contigs to consider for variant calling',default = ['chr'+str(_) for _ in range(1,23,1)] + ['chrX','chrY'])
    parser.add_argument('-p','--threads',type=int,help='prefix of the output files',default = 1)
    parser.add_argument('-pe','--perror',type=float,help='estimated polymerase error rate',default = 10E-5)
    parser.add_argument('-mr','--mutRate',type=float,help='estimated somatic mutation rate per base',default = 3*10E-7)
    parser.add_argument('-pc','--pcut',type=float,help='probability cut-off for making a variant call',default = 0.1)
    parser.add_argument('-mq','--mapq',type=float,help='minumum mapq for an alignment to be considered',default = 30)


    parser.add_argument('-n','--normalBam',type=str,help='deduplicated normal bam')
    parser.add_argument('-m','--noise',type=str,help='noise mask')
    parser.add_argument('-t5','--trim5',type=int,help='ignore mutation if it is less than n bps from 5 prime',default=8)
    parser.add_argument('-t3','--trim3',type=int,help='ignore mutation if it is less than n bps from 3 prime',default=8)
    parser.add_argument('-t52','--trim5DBS',type=int,help='ignore mutation if it is less than n bps from 5 prime',default=30)
    parser.add_argument('-t32','--trim3DBS',type=int,help='ignore mutation if it is less than n bps from 3 prime',default=30)
    parser.add_argument('-d','--minNdepth',type=int,help='minumum coverage in normal for called variants',default=10)
    parser.add_argument('-rl','--readLength',type=int,help='trimmed read length')
    parser.add_argument('-nad','--maxAltCount',type=int,default=0)

    parser.add_argument('-od','--outputDuplex',type=bool,default=False,help='Output bam file for all duplex reads')
    parser.add_argument('-odd','--outputDedup',type=bool,default=False,help='Output bam file for de-duplicated reads. (Warning: may take very long time for large volume data)')
    
    args = parser.parse_args()

    """
    Store Parameters
    """
    params = {"tumorBam": args.bam,\
            "normalBam": args.normalBam,\
            "germline": args.germline,\
            "reference": args.reference,\
            "output": args.output,\
            "regions": args.regions,\
            "threads": args.threads,\
            "perror": args.perror,\
            "mutRate": args.mutRate,\
            "pcutoff": args.pcut,\
            "mapq": args.mapq,\
            "noise": args.noise,\
            "trim5": args.trim5,\
            "trim3": args.trim3,\
            "trim5DBS": args.trim5DBS,\
            "trim3DBS": args.trim3DBS,\
            "minNdepth": args.minNdepth,\
            "readLength": args.readLength,\
            "maxAltCount": args.maxAltCount,\
             }
    """
    Initialze run
    """
    print("..............Loading reference genome.....................")
    fasta = SeqIO.to_dict(SeqIO.parse(args.reference,'fasta'))
    startTime = time.time()
    if not os.path.exists('tmp'):
        os.mkdir('tmp')
    bamObject = BAM(args.bam,'rb')
    
    """
    Execulte variant calling
    """
    if args.threads == 1:
        """
        Single-thread execution
        """
        print(".........Starting variant calling..............")
        #contigs = [(r.strip('\n'),) for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        paramsNow = params
        paramsNow['reference'] = fasta
        #paramsNow['regions'] = contigs
        results = callBam(paramsNow,0,1000000)
        passRead = readStat[0]
        passLength = readStat[1]
        muts  = readStat[2]
        with open(args.output + '_samplestats.txt','w') as stats:
            stats.write("minRBreads\t"+"\t".join([str(_) for _ in range(31)])+'\n')
            stats.write("DuplexReadsNo\t"+"\t".join([str(_) for _ in passRead])+'\n')
            stats.write("EffLength\t"+"\t".join(str(_) for _ in passLength)+'\n')
        mutsAll = muts
        if args.outputDuplex:
            os.rename("tmp/"+args.output+"_0_duplex.bam",args.output+"_duplex.bam")
        if args.outputDedup:
            os.rename("tmp/"+args.output+"_0_dedupped.bam",args.output+"_dedupped.bam")
    else:
        """
        Multi-thread execution
        """
        #contigs = [r.strip('\n') for r in open(args.regions,'r').readlines()] # Only process contigs in region file
        contigs = args.regions
        contigLengths = [bamObject.get_reference_length(contig) for contig in contigs]
        print("...........Spliting genomic regions for parallel execution................")
        cutSites,chunkSize = splitBamRegions(args.bam,args.threads,contigs) #Split the whole genome for parallel execution
        regionSequence = []
        currentContigIndex = 0

        """
        Determine regions for each process
        """

        for nn,site in enumerate(cutSites[1:]):
            pSite = cutSites[nn]
            if site[0] == pSite[0]:
                regionSequence.append((contigs[site[0]],pSite[1],site[1]))
            else:
                if pSite[1] != 0:
                    regionSequence.append((contigs[pSite[0]],pSite[1]))
                else:
                    regionSequence.append((contigs[pSite[0]],))
                for ii in range(pSite[0]+1,site[0]):
                    regionSequence.append((contigs[ii],))
                regionSequence.append((contigs[site[0]],0,site[1]))
        regionSequence.append((contigs[site[0]],site[1]))
        for ii in range(site[0]+1,len(contigs)):
            regionSequence.append((contigs[ii],))
        print("............Completed region splitting in " + str((time.time()-startTime)/60) + " minutes,loading reference genome..................")

        """
        Start variant calling 
        """

        callArguments = []
        startTime2 = time.time()
        print(".........Starting variant calling..............")
        pool = Pool()
        for nn in range(args.threads):
            regions = []
            while len(regionSequence) != 0:
                if len(regionSequence[0]) != 3:
                    regions.append(regionSequence.pop(0))
                else:
                    regions.append(regionSequence.pop(0))
                    break
            print(regions)
            chroms = [region[0] for region in regions]
            fastaNow = {chrom:fasta[chrom] for chrom in chroms} # Takes partial fasta to reduce memory usage
            paramsNow = params
            paramsNow['reference'] = fastaNow
            paramsNow['regions'] = regions
            callArgument = (paramsNow.copy(),nn,chunkSize)
            callArguments.append(callArgument)
            regions = []
        results = pool.starmap(callBam,callArguments)# each result return three list: number of duplex reads, effective lengths, list of mutations
        readStats = [_[0:2] for _  in results]
        muts = [_[2] for _  in results]
        print("..............Completed bam calling in " + str((time.time()-startTime2)/60) + " minutes,merging results................." )
        pool.close()
        pool.terminate()
        pool.join()
        

        """
        Write sample stats
        """        
        passRead = [0] * 31
        passLength = [0] * 31
        for nn in range(args.threads):
            for binn in range(31):
                passRead[binn] += readStats[nn][0][binn]
                passLength[binn] += readStats[nn][1][binn]
        
        with open(args.output + '_samplestats.txt','w') as stats:
            stats.write("minRBreads\t"+"\t".join([str(_) for _ in range(31)])+'\n')
            stats.write("DuplexReadsNo\t"+"\t".join([str(_) for _ in passRead])+'\n')
            stats.write("EffLength\t"+"\t".join(str(_) for _ in passLength)+'\n')

        
        if args.outputDuplex:
            bamNames = ["tmp/"+args.output+'_'+str(nn)+"_dedupped.bam" for nn in range(args.threads)]
            pysam.merge("-@",str(args.threads),"-f",args.output+'_dedupped.bam',*bamNames)
            pysam.index("-@",str(args.threads),args.output+'_dedupped.bam')
        
        if args.outputDedup:
            bamNames = ["tmp/"+args.output+'_'+str(nn)+"_duplex.bam" for nn in range(args.threads)]
            pysam.merge("-@",str(args.threads),"-f",args.output+'_duplex.bam',*bamNames)
            pysam.index("-@",str(args.threads),args.output+'_duplex.bam')
    
        mutsAll = sum(muts,[])
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
    
    """
    write VCF
    """
    tBam = BAM(args.bam,'rb')
    contigs = tBam.references
    #print(contigs)
    chromDict = {contig:tBam.get_reference_length(contig) for contig in contigs}
    infoDict = {}
    formatDict = {"AC":[1,"Integer","Count of alt allele"],"RC":[1,"Integer","Count of ref allele"],"DP":[1,"Integer","Depth at the location"],\
    "NM":[1,"Float","mean NM of alt read group"],"AS":[1,"Float","mean AS of the read group"],"XS":[1,"Float","mean XS of the read group"],\
    "RP5":[1,"Integer","read position"],"RP3":[1,"Integer","distance from 3p"],"DC":[1,"Integer","Number of reads in the duplex group"]}
    filterDict = {"PASS":"All filter Passed"}
    vcfLines = createVcfStrings(chromDict,infoDict,formatDict,filterDict,mutsAll)
    with open(args.output + '.vcf','w') as vcf:
        vcf.write(vcfLines)

    print("..............Completed variant calling " + str((time.time()-startTime)/60) + " minutes..............." )