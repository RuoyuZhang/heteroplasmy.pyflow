#!/usr/bin/python

import argparse
import distutils.util
import sys

def main():
    usage = "usage: %prog [options] \nheteroplasmy calling rewrite by python\n\
            Author: Ruoyu Zhang \n\
            last update Oct 18 2017"
    
    parser = argparse.ArgumentParser(description = usage)
    
    parser.add_argument("--pyflow",type=str,required=False,default='/home/fs01/rz253/project/general/heteroplasmy.pyflow/heteroplasmy.pyflow.py',help="pyflow script dir")
    parser.add_argument("--sample_list",type=str,required=True,help='sample list file')
    parser.add_argument("--out_dir",type=str,required=True,help='output dir')
    parser.add_argument("--config",type=str,required=False,default='/home/fs01/rz253/project/general/heteroplasmy.pyflow/config',help='config file dir')
    parser.add_argument("--het_config",type=str,required=False,default='/home/fs01/rz253/project/general/heteroplasmy.pyflow/het_config',help='output dir')
    parser.add_argument("--stages",type=str,required=False,default='all_rm',help='output dir')
    parser.add_argument("--ref",type=str,required=False,default='/home/fs01/rz253/project/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa',help='reference genome')
    parser.add_argument("--ref_mt",type=str,required=False,default='/home/fs01/rz253/project/reference/GRCh38/chrM.fa',help='mtDNA reference')
    parser.add_argument('--isContinue', type=distutils.util.strtobool, default='False')
    parser.add_argument('--isDryRun', type=distutils.util.strtobool, default='False')
    parser.add_argument('--run', type=int, default=1,help='run task number in parallel')
    parser.add_argument('--script', type=str, default='./pyflow.sh',help='output shell script file')
    
    opts = parser.parse_args()
    
    #begin:
    sample=open(opts.sample_list,"a+")
    out=open(opts.script,"w")
    
    for i, line in enumerate(sample):
        line.strip()
        [name, fq1, fq2] = line.split()
        sh = '%s %s --sample %s --fastq %s %s --out_dir %s/%s --config_file %s --het_config %s --stages %s --ref %s --ref_mt %s --isContinue %s --isDryRun %s' \
                 %(sys.executable,opts.pyflow,name,fq1,fq2,opts.out_dir,name,opts.config,opts.het_config,opts.stages,opts.ref,opts.ref_mt,opts.isContinue,opts.isDryRun)
        
        if i % opts.run != 0:
            sh = sh + '&'
        out.write(sh + '\n')
        

if __name__ == "__main__":
    main()  