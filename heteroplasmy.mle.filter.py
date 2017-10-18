#!/usr/bin/python

from optparse import OptionParser
from scipy import stats
#import sys
#import os
#import re
import math


def probe_mode(data,low,f):
    result=dict()
    for key in ['A','T','G','C']:
        ndep=data[key][2]
        frequency=ndep/float(data['cover'])
        if ndep >= low and frequency >= f:
            result[key]=frequency
        
    return result


def loose_mode(data,depth,low,f):
    result=dict()
    for key in ['A','T','G','C']:
        if data[key][0] >= math.floor(low[0]*depth) and data[key][1] >= math.floor(low[1]*depth):
            ndep=data[key][2] 
            frequency=ndep/float(data['cover'])
            if frequency >= f:
                result[key]=frequency

    return result


def chi_mode(data,depth,low,alpha,f):
    result=dict()
    plus=data['A'][0]+data['T'][0]+data['G'][0]+data['C'][0]
    minus=data['A'][1]+data['T'][1]+data['G'][1]+data['C'][1]
    for key in ['A','T','G','C']:
        if data[key][0] >= low[0]*depth and data[key][1] >=low[1]*depth:
            ndep=data[key][2]
            frequency=ndep/float(data['cover'])
            if frequency >= f:
            #add chi square test:
                if frequency > 0.5:
                    result[key]=frequency
                else:
                    a=data[key][0]
                    b=data[key][1]
                    c=plus-data[key][0]
                    d=minus-data[key][1]

                    least=sorted([a,b,c,d])[0]
                    table=[[a,b],[c,d]]

                    if least < 5:
                        pvalue=stats.fisher_exact(table)[1]
                    else:
                        pvalue=stats.chi2_contingency(table)[1]

                    if pvalue > alpha:
                        result[key]=frequency
    return result


def main():
    usage = "usage: %prog [options] \nheteroplasmy calling rewrite by python\n\
            Author: Ruoyu Zhang \n\
            last update April 25th 2017 \n\
            This is used for pick up heteroplasmy and mutation site from raw data produced by heteroplasmy.mle.py"
    
    parser = OptionParser(usage=usage)
    
    parser.add_option("-i",help="input het.raw file")
    parser.add_option("-o",help="output file prefix, output o.heteroplasmy and o.homoplasmy")
    parser.add_option("-f",type=float,default=0.01,help="heteroplasmy cutoff[0.01]")
    parser.add_option("--loose",default="0.003,0.003",help="loose mode, defaul on [0.003,0.003]")
    parser.add_option("--chi",type=float,default=0.0,help="chi-square test mode default off, set to 0 to disable [0.0]")
    parser.add_option("--probe",type=float,default=0,help="probe mode, defaul off, set cutoff [0]")
    parser.add_option("-d",type=int,default=1000,help="site depth cutoff[1000]")
    parser.add_option("--mle",type=float,default=0,help="mle value, set to 0 to disable [0]")
    
    (opts,args)=parser.parse_args()
    
    lcut=opts.loose.split(',')
    lcut=map(float,lcut)
    chi=opts.chi
    freq=opts.f
    probe=opts.probe
    d=opts.d
    mlecut=opts.mle
    
    #begin:
    raw=open(opts.i,"a+")
    of1=open(opts.o+".heteroplasmy","w")
    of2=open(opts.o+".homoplasmy","w")
    
    for i, line in enumerate(raw):
        line.rstrip()
        
        if line.find("not_reliable")!= -1:
            continue
    
        [ch,pos,ref,cover,an,tn,gn,cn,maf,mle]=line.split()[0:10]
        mle=float(mle)
        maf=float(maf)
        cover=int(cover)
    
        if int(cover)<d:
            continue
    
        content={'ref':ref,'cover':int(cover)}
        for tide in [an,tn,gn,cn]:
            tidelist=tide.split(':')
            content[tidelist[0]]=[int(tidelist[2]),int(tidelist[3]),int(tidelist[1])]
        
        if probe != 0:
            chi = 0
    
        if chi > 0:
            result=chi_mode(data=content,depth=cover,low=lcut,alpha=chi,f=freq)
        elif probe != 0:
            result=probe_mode(data=content,low=probe,f=freq)
        else:
            result=loose_mode(data=content,depth=cover,low=lcut,f=freq)
        
            
        if len(result.keys())>=2 and mle >= mlecut:
            mtype='H'
        elif len(result.keys())==1:
            if result.keys()[0] == ref:
                #print pos,result.keys()[0],ref
                continue
            else:
    #            print result.keys()[0],ref
                mtype='M'
        else:
            continue
        
        if mtype=='M':
            of2.write("%-10s %-10s %-7s %-10s %-15s %-15s %-15s %-15s %-7.5f %-7.5f" %(ch,pos,ref,cover,an,tn,gn,cn,maf,mle))
            of2.write('\t'+ mtype)
            for key in result.keys():
                of2.write('\t'+key+' '+str(result[key]))
    
            of2.write('\n')
    
        elif mtype=='H':
            of1.write("%-10s %-10s %-7s %-10s %-15s %-15s %-15s %-15s %-7.5f %-7.5f" %(ch,pos,ref,cover,an,tn,gn,cn,maf,mle))
            of1.write('\t'+ mtype)
            for key in result.keys():
                of1.write('\t'+key+' '+str(result[key]))
    
            of1.write('\n')


if __name__ == "__main__":
    main()  