###read in on-target.bam ###
###set range as 6:43737420 - 43737490 ###
###at a 20 bp window, characterize the type of indels ###



### read in region of interest, create cigar per base for each read ###
def read_in_sortsam(sortsam):
    M_array = HTSeq.GenomicArray(  ["6" ], stranded=False, typecode="i" )
    I_array = HTSeq.GenomicArray( ["6" ], stranded=False, typecode="i" )
    D_array = HTSeq.GenomicArray( ["6" ], stranded=False, typecode="i" )
    bundles=[]
    Total_read=0
    indel_read=0
    for aln in HTSeq.SAM_Reader(sortsam):
        """ only take alignment map right at the location """
        if aln.iv.start < 43737440 and aln.iv.end > 43737380:
            Total_read+=1
            for c in aln.cigar:
                if c.ref_iv.start==c.ref_iv.end:
                    aln_end=c.ref_iv.end+1
                    aln_start=c.ref_iv.start  
                else:
                    aln_end=c.ref_iv.end
                    aln_start=c.ref_iv.start  
                if c.type != "S" and c.type != "H" and aln_start < aln_end:
                    if c.type== "I":
                        I_array[HTSeq.GenomicInterval( c.ref_iv.chrom, aln_start, aln_end) ] += 1
                    elif c.type== "M":
                        M_array[HTSeq.GenomicInterval( c.ref_iv.chrom, aln_start, aln_end) ] += 1
                    elif c.type== "D":
                        D_array[HTSeq.GenomicInterval( c.ref_iv.chrom, aln_start, aln_end) ] += 1
            if [c.type for c in aln.cigar if c.type=="D" or c.type=="I"] != []:
                bundles.append(aln)
                indel_read +=1
    return Total_read, indel_read, bundles, M_array, I_array, D_array

    
### clean up the alignment, only pick the most relevant type of indel for the read ###    
### set hot region at 43737367 ###
def determine_key_indel(cigars):
    relevant_indel=[c for c in cigars if c.type=="D" or c.type=="I"]
    ### if multiple indels, calculate mid-point of the indel ###
    #relevant_indel.sort(key=lambda x: abs(np.divide(x.ref_iv.start+x.ref_iv.end, 2) - 43737367)) 
    #relevant_indel_by_location=relevant_indel[0]
    relevant_indel.sort(key=lambda x: x.size, reverse=True)
    relevant_indel_by_size=relevant_indel[0]
    return relevant_indel_by_size

def characteriz_indel(bundles):
    insert_bin={i:[] for i in range(1,21)}
    delete_bin={i:[] for i in range(1,21)}
    for b in bundles:
        cigars=b.cigar
        key_indel=determine_key_indel(cigars)
        size=key_indel.size
        if size > 20:
            size=20
        if key_indel.type == "I":
            insert_bin[size].append(key_indel)
        elif key_indel.type == "D":
            delete_bin[size].append(key_indel)
    return insert_bin, delete_bin

def count_per_base_Guide(I_array, D_array, M_array):
    I_list=[]
    M_list=[]
    D_list=[]  
    for j in range(43737444, 43737491):
        read_iv = HTSeq.GenomicInterval( "6", j, j+1 )
        I_list.append([[iv.start, val] for iv, val in I_array[ read_iv ].steps()][0])
        M_list.append([[iv.start, val] for iv, val in M_array[ read_iv ].steps()][0])
        D_list.append([[iv.start, val] for iv, val in D_array[ read_iv ].steps()][0])
    total_list=[[m[0], m[1]+i[1]+d[1]] for m, i, d in zip(I_list, M_list, D_list)]
    final_list=[[str(t[0]), str(t[1]), str(i[1]), str(d[1])] for t, i, d in zip(total_list, I_list, D_list)]
    return final_list

def output_indel_profile_perbase(samplename, final_list, out_profile):
    header="samplename\tchrom\tpos\ttotal_cov\tinsertion\tdeletion\n"
    [i.insert(0, samplename) for i in final_list]
    [i.insert(1, "6") for i in final_list]
    output=["\t".join(j) for j in final_list]
    output="\n".join(output)
    o=open(out_profile, "w")
    o.write(header)
    o.write(output)
    o.close()


def output_indel_details(samplename, Total_read, indel_read, insert_bin, delete_bin, outfile):
    header="samplename\tTotal_read\tIndel_read\tIndel_type\tIndel_size\treads\n"
    output=[]
    Total_read=str(Total_read)
    indel_read=str(indel_read)
    for i in range(1, 21):
        row=[samplename, Total_read, indel_read, "I", str(i), str(len(insert_bin[i]))]
        output.append(row)
    for i in range(1, 21):
        row=[samplename, Total_read, indel_read, "D", str(i), str(len(delete_bin[i]))]
        output.append(row)
    output=["\t".join(j) for j in output]
    output="\n".join(output)
    o=open(outfile, "w")
    o.write(header)
    o.write(output)
    o.close()


        
if __name__ == "__main__":    
    import os, sys, argparse
    import HTSeq
    import subprocess
    import pyfaidx
    import numpy as np
    import re
        

    sortsam=sys.argv[1]
    samplename=sys.argv[2]
    outfile=sys.argv[3]
    out_profile=sys.argv[3] + ".out_profile"
    Total_read, indel_read, bundles, M_array, I_array, D_array=read_in_sortsam(sortsam)
    insert_bin, delete_bin=characteriz_indel(bundles)
    output_indel_details(samplename, Total_read, indel_read, insert_bin, delete_bin, outfile)
    final_list=count_per_base_Guide(I_array, D_array, M_array)
    output_indel_profile_perbase(samplename, final_list, out_profile)
