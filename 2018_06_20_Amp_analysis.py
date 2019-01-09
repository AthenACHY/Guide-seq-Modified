import os, sys, argparse
import HTSeq
import subprocess
import pyfaidx

inbam="/data/AChu_data/2018_06_19_Vivian_amp/Lib01_on_target.bam"
wkdir="/data/AChu_data/2018_06_19_Vivian_amp/"
samplename="Lib01"
def create_sort_sam(inbam, wkdir, samplename):
    insam=wkdir + samplename + ".sam"
### remove supplementary alignment ###
    cmd=["samtools", "view", "-F", "2048", inbam]
    subprocess.call(cmd, stdout=open(insam, "w"))
    cmd=["sort", "-k1" , insam]
    sortsam=wkdir + samplename + "_sort.sam"
    subprocess.call(cmd, stdout=open(sortsam, "w"))
    cmd=["rm", insam]
    subprocess.call(cmd)
    return sortsam

def match_base_filter(aln):
    """ check if "M" >= 30 bp """
    m=0
    for c in aln.cigar:
        if c.type=="M":
            m+=c.size
    if m >=30:
        return True
    else:
        return False

def check_AS_filter(aln):
    """check AS filter"""
    AS=[f for f in aln.optional_fields if f[0]=="AS"]
    if AS !=[]:
        if AS[0][1] >=50:
            return True
        else:
            return False
    else:
        return False

def check_MD_filter(aln, region_start, region_end):
    """check MD filter"""
    SNP_genomic_coord=[]
    SNP_bool=[]
    MD=[f for f in aln.optional_fields if f[0]=="MD"]
    if MD !=[]:
        SNPs=re.findall("[ATGC]|\d+", MD[0][1])
        matches_stretches=[[c.ref_iv, c.size] for c in aln.cigar if c.type=="M"]       
        qcoord=0
        for s in SNPs:
            try:
### intergers mean matches, all to the coord ### 
                qcoord+= int(s)
            except:
### encounter a SNP ###
                qcoord +=1
### now identify which stretch of matches the SNP is in ###
                match_location=0
                for m in matches_stretches:
                    match_location +=m[1]
                    if match_location >=qcoord:
### Qcoord not adjusted properly ###
                        adjust_qcoord=qcoord - match_location +m[1] 
                        SNP_location=m[0].start + adjust_qcoord
                        SNP_genomic_coord.append([s, SNP_location])
                        if SNP_location >= region_start and SNP_location <=region_end:
                            SNP_bool.append(True)
                            break
                        else:
                            SNP_bool.append(False)
    if SNP_bool.count(True) >0:
        return SNP_genomic_coord, True
    else:
        return SNP_genomic_coord, False 
                    



def soft_clip_filter(ODNS, aln):
    """check if "S" part matches with ODNS"""
    """ cigar string coordinate follows the orientation of the alignment (+) """
    """ but aln.read follow the read orientation, can be (+) or (-) """
    s=[]
    for c in aln.cigar:
        if c.type=="S":
            start=c.query_from
            end=c.query_to
            try:
                s.append(aln.read_as_aligned.seq[start:end+1])
            except:
                s.append(aln.read_as_aligned.seq[start:end+1])
    odn_tag=[t in p for t in ODNS for p in s]
    if odn_tag.count(True) > 0:
        return True
    else:
        return False

def check_ODN_presence_insertion(aln):
    """ check if the insertion contain any ODNS middle part sequence """
    inserted_seq=[]
    ODNS_middle=["GACAACTCAA", "GTCATATGTTAAT"]
    for c in aln.cigar:
        if c.type=="I" and c.size >20:
            start=c.query_from
            end=c.query_to
            inserted_seq0=aln.read_as_aligned.seq[start:end+1]
            inserted_seq.extend([O for O in ODNS_middle if O in inserted_seq0])
    if inserted_seq !=[]:              
        return True
    else:
        return False





amplicon_indels=HTSeq.GenomicArrayOfSets("auto", stranded=False )
amplicon_cleavage=HTSeq.GenomicArrayOfSets("auto", stranded=False )
amplicon_SNP_spectrum=HTSeq.GenomicArrayOfSets("auto", stranded=False )

ODNplus="ACATATGACAACTCAATTAAAC"
ODNminus="TTGAGTTGTCATATGTTAATAACGGTA"
ODNplus_RC="GTTTAATTGAGTTGTCATATGT"
ODNminus_RC="TACCGTTATTAACATATGACAACTCAA"
ODNS=["ACATATGACAACTCAATTAAAC", "TTGAGTTGTCATATGTTAATAACGGTA", "GTTTAATTGAGTTGTCATATGT", "TACCGTTATTAACATATGACAACTCAA"]
### readingin sam line by line or should read in bundle? ###

def report_indel_from_amp(sortsam, ODNS, region_start, region_end):
    Total_read=0
    ODNS_inserted=0
    INDEL_read=0
    substitution_read=0
    for aln in HTSeq.SAM_Reader(sortsam):
### one alignment at a time, alignment independent of read name ###
### only assign one type of mutation for one read according to priority ODNs> INDEL> substitution ###
        if aln.aligned and match_base_filter(aln) and check_AS_filter(aln):
### pass AS and # matches filter ###
            Total_read += 1
### count number of reads being edited ###
            if soft_clip_filter(ODNS, aln):
### soft-clipped part contain ODNs == ODN insertion ###
                ODNS_inserted +=1
            else:
                cigar_property=[(c.type, c.size) for c in aln.cigar]
                if [p for p in cigar_property if p[0] == "I" and p[1] > 20] !=[]:
                    if check_ODN_presence_insertion(aln):
### count ODNS insertion marked as I ###
                        ODNS_inserted +=1
                    else:
                        INDEL_read +=1
### these reads were either modified or not modified, so record all the indel mutations ###   
                elif [p[0] for p in cigar_property if p[1] < 20 and p[0] =="D" or p[0]=="I"] :
### distinguish insertion from ODNS insertion ###
                    indel_region=[[c.ref_iv.start, c.ref_iv.end] for c in aln.cigar if c.type == "D" or c.type == "I" and c.size < 20]
                    valid_indel=[r for r in indel_region if r[0] >=region_start and r[1] <= region_end]
                    if len(valid_indel) >0: 
### count indel over substitutions ###
                        INDEL_read +=1             
            #### Cannot sort out substitution ###
                else:
                    if "M" in [p[0] for p in cigar_property]:   
### these are reads having everything as matches, so check MD to see if there are any mutations ###
                        SNP_genomic_coord, MD_bool= check_MD_filter(aln, region_start, region_end)
                        print SNP_genomic_coord
### check if subtitution within pcr region???###
                        if MD_bool ==True and [s[0] for s in SNP_genomic_coord if s[1] != 43737486] !=[]:
                            substitution_read +=1
    return Total_read, ODNS_inserted, INDEL_read, substitution_read


import numpy as np
wkdir="/data/AChu_data/2018_06_19_Vivian_amp/"
region_start=43737419
region_end=43737497
outdata=[]
for n in range(1, 9):
    samplename="Lib0" +str(n)
    inbam=wkdir + samplename + "_on_target.bam"
    sortsam=create_sort_sam(inbam, wkdir, samplename)
    Total_read, ODNS_inserted, INDEL_read, substitution_read=report_indel_from_amp(sortsam, ODNS, region_start, region_end)
    editing_eff=np.divide(INDEL_read+ODNS_inserted+substitution_read, Total_read, dtype="f4")
    ODNs_insert_eff=np.divide(ODNS_inserted, INDEL_read+ODNS_inserted+substitution_read, dtype="f4")
    outdata.append([samplename, str(Total_read), str(ODNS_inserted), str(INDEL_read), str(substitution_read), str(editing_eff), str(ODNs_insert_eff)])
              
### GUIDE1 ###
region_start=73160892
region_end=73160962
samples=["09", "10", "11", "12", "13", "14", "15", "16"]
for n in samples:
    samplename="Lib" +n
    inbam=wkdir + samplename + "_on_target.bam"
    sortsam=create_sort_sam(inbam, wkdir, samplename)
    Total_read, ODNS_inserted, INDEL_read, substitution_read=report_indel_from_amp(sortsam, ODNS, region_start, region_end)
    editing_eff=np.divide(INDEL_read+ODNS_inserted+substitution_read, Total_read, dtype="f4")
    ODNs_insert_eff=np.divide(ODNS_inserted, INDEL_read+ODNS_inserted+substitution_read, dtype="f4")
    outdata.append([samplename, str(Total_read), str(ODNS_inserted), str(INDEL_read), str(substitution_read), str(editing_eff), str(ODNs_insert_eff)])

### GUIDE_3 ###
region_start=73161080
region_end=73161150
samples=["17", "18", "19", "20", "21", "22", "23", "24"]
for n in samples:
    samplename="Lib" +n
    inbam=wkdir + samplename + "_on_target.bam"
    sortsam=create_sort_sam(inbam, wkdir, samplename)
    Total_read, ODNS_inserted, INDEL_read, substitution_read=report_indel_from_amp(sortsam, ODNS, region_start, region_end)
    editing_eff=np.divide(INDEL_read+ODNS_inserted+substitution_read, Total_read, dtype="f4")
    ODNs_insert_eff=np.divide(ODNS_inserted, INDEL_read+ODNS_inserted+substitution_read, dtype="f4")
    outdata.append([samplename, str(Total_read), str(ODNS_inserted), str(INDEL_read), str(substitution_read), str(editing_eff), str(ODNs_insert_eff)])


outfile="/data/AChu_data/2018_06_19_Vivian_amp/out_data.txt"

o=open(outfile, "w")
o.write("samplename\tTotal_read\tODNS_inserted_reads\tINDEL_read\tsubstitution_read\tediting_efficiency\tODNs_insert_efficiency") 
for g in outdata:
    o.write("\t".join(g)+"\n")

o.close()



          
