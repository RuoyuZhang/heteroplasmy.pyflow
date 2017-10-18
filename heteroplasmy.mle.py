#!/usr/bin/python

from optparse import OptionParser
#import sys
#import os
import re
import math

def frequency_likelihood(f,error_ref,error_mut):
    if (f>=0 and f<=1):
        part1 = 1
        for j in range(len(error_mut)):
            QQ = ord(error_mut[j])-33
            error_rate= 10**(-QQ/10.0)
            part1=part1*((1-f)*error_rate + f*(1-error_rate))
#        print part1
        part2 = 1
        for j in range(len(error_ref)):
            QQ = ord(error_ref[j])-33
            error_rate= 10**(-QQ/10.0)
            part2=part2*((1-f)*(1-error_rate) + f*error_rate)
#        print part2
        lh = part1 * part2
        return lh
    else:
        return 0.0


def main():
    usage = "usage: %prog [options] \nheteroplasmy calling rewrite by python\n\
            Author: Ruoyu Zhang \n\
            April 9th 2015\n\
            last modified April 25th 2017\n\
            this is only used for produce all raw counts, use heteroplsmy.filter.py to extract heteroplasy sites"
    
    parser = OptionParser(usage=usage)
    
    parser.add_option("-i",help="input mpileup file")
    parser.add_option("-o",help="output file")
    parser.add_option("-Q",type=int,default=20,help="quality cutoff[20]")
    #parser.add_option("--rmrate",type=float,default=0.7,help="Ratio of bases remained after quality filter[0.7]")
    
    (opts,args)=parser.parse_args()
    
    
    Q=opts.Q
    
    mpfile=open(opts.i,"a+")
    hfile=open(opts.o,"w")
    change={'A':'a','T':'t','G':'g','C':'c','a':'A','t':'T','g':'G','c':'C'}
    
    
    for i, line in enumerate(mpfile):
        line.rstrip()
    #    print "yes"
        if len(line.split())<6:
            continue
        else:
            [ch,site,ref,depth,bases,qual]=line.split()
        
        if ref == "N" or ref == "n":
            continue
        bases=re.sub('\^.|\$','',bases)
    
        indel=''
        indelseq=''
        allindel=dict()
    
        while re.search('(\+|-)(\d+)',bases):
            m=re.search('(\+|-)(\d+)',bases)
            indel=m.group(1)
            num=m.group(2)
            if indel=="+":
                capindel=re.search('(\+)(%s)([ATGCNatgcn]{%s})' %(num,num),bases).group()
                indelnum=bases.count(capindel)
                allindel[capindel]=indelnum
                bases=bases.replace(capindel,'')
    #            bases=re.sub('(\+)(%s)([ATGCNatgcn]{%s})' %(num,num),'',bases,count=1)
            elif indel=="-":
                capindel=re.search('(-)(%s)([ATGCNatgcn]{%s})' %(num,num),bases).group()
                indelnum=bases.count(capindel)
                allindel[capindel]=indelnum
                bases=bases.replace(capindel,'')
    
        bd={'A':0,'T':0,'G':0,'C':0,'a':0,'t':0,'g':0,'c':0}
        #base quality
        bdq={'A':'','T':'','G':'','C':'','a':'','t':'','g':'','c':''}
    
        for i in range(len(bases)):
            q=ord(qual[i])-33
            if (q >= Q and bases[i] != "n" and bases[i] != "N"):
                if re.search('[ATGCatgc]',bases[i]):
                    bd[bases[i]]=bd[bases[i]] + 1
                    bdq[bases[i]]=bdq[bases[i]] + qual[i]
                elif re.search('\.',bases[i]):
                    bd[ref]=bd[ref]+1
                    bdq[ref]=bdq[ref]+qual[i]
                elif re.search(',',bases[i]):
                    bd[change[ref]]=bd[change[ref]]+1
                    bdq[change[ref]]=bdq[change[ref]]+qual[i]
    
        dep=0
        for ba in ['A','T','G','C']:
            bd[ba+change[ba]]=bd[ba]+bd[change[ba]]
            bdq[ba+change[ba]]=bdq[ba]+bdq[change[ba]]
            dep=dep+bd[ba+change[ba]]
        
        if dep==0:
            continue
    
        hfile.write("%-10s %-10s %-7s %-10s " %(ch,site,ref,dep))
    
        for ba in ['A','T','G','C']:
            if dep==0:
                freq=0
            else:
                freq = float(bd[ba+change[ba]])/dep
            
            hfile.write("%s:%s:%s:%s:%-7.5f " %(ba,bd[ba+change[ba]],bd[ba],bd[change[ba]],freq))
    
    #calculate mle
        dict1={'Aa':bd['Aa'],'Tt':bd['Tt'],'Gg':bd['Gg'],'Cc':bd['Cc']}
        dict_sort= sorted(dict1.iteritems(), key=lambda d:d[1], reverse = True)
        #print dict
        
        if dep == 0:
            continue
        
        mle=float(dict_sort[1][1])/dep
        hfile.write("%-7.5f " %mle)
    
        L1 = frequency_likelihood(f=mle,error_ref=bdq[dict_sort[0][0]],error_mut=bdq[dict_sort[1][0]])
        L0 = max(frequency_likelihood(f=0.0,error_ref=bdq[dict_sort[0][0]],error_mut=bdq[dict_sort[1][0]]),frequency_likelihood(f=1.0,error_ref=bdq[dict_sort[0][0]],error_mut=bdq[dict_sort[1][0]]))
        if L0 == 0:
            LLR=10
        elif L1 == 0:
            LLR=0
            
        else:
            LLR = -math.log10(float((L0)/(L1)))
    
        hfile.write("%-7.5f " %LLR)
    
    #    print allindel 
        hfile.write(depth+" ")
        for key1 in allindel.keys():
           # for key2 in allindel[key1]:
            hfile.write(key1+":"+str(allindel[key1])+" ")
        
        #if float(dep)/(float(depth)+0.000001) < float(opts.rmrate):
        #    print >>hfile,"not_reliable",
    
        hfile.write('\n')


if __name__ == "__main__":
    main()  


