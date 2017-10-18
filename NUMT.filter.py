#!/usr/bin/python

import argparse
import distutils.util
import re

def main():
    usage = "usage: %prog [options] \nheteroplasmy calling rewrite by python\n\
            Author: Ruoyu Zhang \n\
            last update Oct 18 2017"
    
    parser = argparse.ArgumentParser(description = usage)
    
    parser.add_argument("--sam",type=str,required=True,help="input .sam file")
    parser.add_argument("--out",type=str,required=True,help="output .sam file")
    parser.add_argument("--max_mis",type=float,default=0.05,help="mapping mismathes rate cutoff [0.01]")
    parser.add_argument("--chr",type=str,default='MT',help="chromosome name [MT]")
    
    opts = parser.parse_args()
    
    #begin:
    sam=open(opts.sam,"a+")
    out=open(opts.out,"w")
    of2=open(opts.out,"w")
    
    [last, read1, read2, isread2, discard] = ['','','','no','no']
    
    for i, line in enumerate(sam):
        line.strip()
        if line[0] == "@":
            if line.find('PG') != -1 or line.find(opts.chr) != -1:
                out.write(line)
            continue
        
        info = line.split()
        #print info
        [name,flag,c,pos,ms,match,pair,pair_p,dis,seq,qua] = info[0:11]
        tags = info[11:]
        
        if name == '' or name != last:
            [read1, last, isread2, discard] = [line, name, 'no', 'no']
        else:
            [read2, last, isread2] = [line, name, 'yes']
            if discard == 'yes':
                continue
        
        # discard non mitochondrial reads
        if c == "*" or c != opts.chr:
            discard = 'yes'
            
        # discard reads with mismatch > cutoff
        if line.find('NM:i') != -1:
            mismatch = re.search(r'NM:i:(\d+)',line).group(1)
            if int(mismatch) > len(seq) * opts.max_mis:
                discard = 'yes'
        else:
            discard = 'yes'
                    
        if discard != 'yes' and isread2 == 'yes':
            out.write(read1 + read2 )


if __name__ == "__main__":
    main()  