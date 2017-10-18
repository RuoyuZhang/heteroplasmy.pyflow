#!/usr/bin/env python

import os
import sys
sys.path.append("/home/fs01/rz253/.local/lib/python2.7/site-packages/pyflow/pyflow/src")
from pyflow import WorkflowRunner
import argparse
import distutils.util

def Config(config_file):
    config = dict()
    with open(config_file) as fh:
        for i,line in enumerate(fh):
            line = line.strip()
            [name,content] = line.split('=')
            config[name] = content
    
    return config


class Heteroplasmy(WorkflowRunner):

    def __init__(
            self,
            sample,
            fastq,
            out_dir,
            ref,
            ref_mt,
            config_file,
            het_config,
            stages,
    ):

        for i in ('sample','fastq','out_dir','ref','ref_mt','config_file','het_config','stages'):
            if i not in locals():
                raise Exception('--%s is required for %s' % (i, self.__class__.__name__))
            setattr(self, i, locals()[i])

        self.config = Config(self.config_file)
        self.het_config = Config(self.het_config)
        
        if not os.path.exists(self.out_dir + '/stat'):
            os.makedirs(self.out_dir + '/stat')
              
        self.prefix = '%s/%s' % (self.out_dir, self.sample)
            
        # trimmatic
        if len(self.fastq) == 2:
            self.trim_command = '%s PE -threads 3 %s %s %s.trimmed.R1.fastq %s.trimmed.R1.un.fastq %s.trimmed.R2.fastq %s.trimmed.R2.un.fastq ILLUMINACLIP:/home/fs01/rz253/bin/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>%s.trim.log' \
                                %(self.config['trim'], self.fastq[0], self.fastq[1], self.prefix, self.prefix, self.prefix,self.prefix,self.prefix)
        else:
            self.trim_command = '%s SE -threads 3 %s %s.trimmed.fastq ILLUMINACLIP:/home/fs01/rz253/bin/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>%s.trim.log' \
                                %(self.config['trim'], self.fastq[0], self.prefix, self.prefix)
        # bwa
        if len(self.fastq) == 2:
            self.bwa_command = '%s mem -t 3 %s %s.trimmed.R1.fastq %s.trimmed.R2.fastq > %s.bwa.sam 2>%s.bwa.log' \
                               % (self.config['bwa'], self.ref, self.prefix, self.prefix, self.prefix, self.prefix)
        else:
            self.bwa_command = '%s mem -t 3 %s %s.trimmed.fastq > %s.bwa.sam 2>%s.bwa.log' \
                               % (self.config['bwa'], self.ref, self.prefix, self.prefix, self.prefix)
            
        # NUMTs filter
        self.numt_command = '%s %s --sam %s.bwa.sam --out %s.bwa.filter.sam --max_mis %s --chr %s' \
                            %(sys.executable, self.config['numt'], self.prefix, self.prefix, self.het_config['map_max_mis'], self.het_config['chr'])
        
        # sam_dedup
        self.sam_dedup_cm1 = '%s SamFormatConverter INPUT=%s.bwa.filter.sam OUTPUT=%s.bwa.bam 2>%s.log' %(self.config['picard'],self.prefix,self.prefix,self.prefix)
        self.sam_dedup_cm2 = '%s SortSam INPUT=%s.bwa.bam OUTPUT=%s.bwa.sort.bam SORT_ORDER=coordinate 2>>%s.log'\
                        % (self.config['picard'],self.prefix,self.prefix,self.prefix)
        self.sam_dedup_cm3 = '%s AddOrReplaceReadGroups INPUT=%s.bwa.sort.bam OUTPUT=%s.sort.gp.bam SORT_ORDER=coordinate RGID=%s.gp RGLB=%s.gb RGPL=illumina RGSM=%s.gp RGPU=barcode 2>>%s.log'\
                           % (self.config['picard'],self.prefix,self.prefix,self.sample,self.sample,self.sample,self.prefix)
        self.sam_dedup_cm4 = '%s MarkDuplicates REMOVE_DUPLICATES=true INPUT=%s.sort.gp.bam OUTPUT=%s.sort.gp.rmdup.bam M=%s.duplicate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=False 2>>%s.log'\
                             % (self.config['picard'],self.prefix,self.prefix,self.prefix,self.prefix)
        self.sam_dedup_cm5 = '%s BuildBamIndex INPUT=%s.sort.gp.rmdup.bam 2>>%s.log' %(self.config['picard'],self.prefix,self.prefix) 
        
        self.sam_dedup_cm6 = '%s CollectAlignmentSummaryMetrics R=%s I=%s.sort.gp.bam O=%s/stat/alignment_metrics.txt' % (self.config['picard'],self.ref_mt,self.prefix,self.out_dir)
        self.sam_dedup_cm7 = '%s CollectInsertSizeMetrics INPUT=%s.sort.gp.rmdup.bam OUTPUT=%s/stat/insert_metrics.txt HISTOGRAM_FILE=%s/stat/insert_size_histogram.pdf'\
                             % (self.config['picard'],self.prefix,self.out_dir,self.out_dir)
        self.sam_dedup_cm8 = '%s depth -a %s.sort.gp.rmdup.bam > %s/stat/depth_out.txt' % (self.config['samtools'],self.prefix,self.out_dir)
           
        # indel realigner
        self.realign_cm1 = '%s -T RealignerTargetCreator -I %s.sort.gp.rmdup.bam -R %s -o %s.indel.list 2>>%s.log' \
                           % (self.config['gatk'],self.prefix,self.ref_mt,self.prefix,self.prefix)
        self.realign_cm2 = '%s -T IndelRealigner -l INFO -I %s.sort.gp.rmdup.bam -R %s -targetIntervals %s.indel.list -o %s.sort.gp.rmdup.realign.bam 2>>%s.log'\
                           % (self.config['gatk'],self.prefix,self.ref_mt,self.prefix,self.prefix,self.prefix)
           
        # variants calling
        self.vcf_cm1 = '%s -T HaplotypeCaller -R %s -I %s.sort.gp.rmdup.realign.bam -o %s.raw_variants.vcf 2>>%s.log'\
                       % (self.config['gatk'],self.ref_mt,self.prefix,self.prefix,self.prefix)
        self.vcf_cm2 = '%s -T SelectVariants -R %s -V %s.raw_variants.vcf -selectType SNP -o %s.raw_snps.vcf 2>>%s.log'\
                         % (self.config['gatk'],self.ref_mt,self.prefix,self.prefix,self.prefix)
        self.vcf_cm3 = '%s -T SelectVariants -R %s -V %s.raw_variants.vcf -selectType INDEL -o %s.raw_indels.vcf 2>>%s.log'\
                         % (self.config['gatk'],self.ref_mt,self.prefix,self.prefix,self.prefix)
        self.vcf_cm4 = "%s -T VariantFiltration -R %s -V %s.raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName \"basic_snp_filter\" -o %s.filtered_snps_final.vcf 2>>%s.log" \
                         % (self.config['gatk'],self.ref_mt,self.prefix,self.prefix,self.prefix)
        self.vcf_cm5 = "%s -T VariantFiltration -R %s -V %s.raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName \"basic_indel_filter\" -o %s.filtered_indels_final.vcf 2>>%s.log" \
                         % (self.config['gatk'],self.ref_mt,self.prefix,self.prefix,self.prefix)
        
        # mpileup
        self.mp_command = '%s mpileup -d 1000000 -f %s %s.sort.gp.rmdup.realign.bam > %s.mp 2>>%s.log'\
                          %(self.config['samtools'],self.ref_mt,self.prefix,self.prefix,self.prefix)
        
        # heteroplasmy
        self.het_raw_cm = '%s %s -i %s.mp -o %s.raw' %(sys.executable, self.config['het_raw'], self.prefix, self.prefix)
        self.het_filter_cm = '%s %s --loose %s --chi %s -d %s --mle %s -i %s.raw -o %s' \
                             %(sys.executable, self.config['het_filter'], self.het_config['loose'], self.het_config['chi'], self.het_config['d'], self.het_config['mle'], self.prefix, self.prefix)
  
        # rm
        self.rm_command = 'rm %s.trimmed.R1.fastq %s.trimmed.R1.un.fastq %s.trimmed.R2.fastq %s.trimmed.R2.un.fastq %s.bwa.sam %s.bwa.filter.sam %s.bwa.bam %s.bwa.sort.bam %s.sort.gp.bam %s.sort.gp.rmdup.bam %s.sort.gp.rmdup.bai '\
                          %(self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,self.prefix,)
                  
    def workflow(self):

        trim_task = None
        if 'all' in self.stages or 'trim' in self.stages:
            trim_task = self.addTask('trim', self.trim_command, isForceLocal = True, nCores=1)
        
        bwa_task = None    
        if 'all' in self.stages or 'bwa' in self.stages:
            bwa_task = self.addTask('bwa', self.bwa_command, dependencies=trim_task)
            
        numt_task = None
        if 'all' in self.stages or 'NUMT' in self.stages:
            numt_task = self.addTask('numt_filter', self.numt_command, dependencies=bwa_task)
        
        sam_dedup_task = {'sam_convert':None, 'sam_sort':None, 'sam_add_RG':None, 'sam_dedup':None, 'sam_index':None, 'align_stat':None, 'insert_stat':None, 'depth_stat':None}
        if 'all' in self.stages or 'sam_dedup' in self.stages:
            sam_dedup_task['sam_convert'] = self.addTask('sam_convert', self.sam_dedup_cm1, dependencies=bwa_task)
            sam_dedup_task['sam_sort'] = self.addTask('sam_sort', self.sam_dedup_cm2, dependencies=sam_dedup_task['sam_convert'])
            sam_dedup_task['sam_add_RG'] = self.addTask('sam_add_RG', self.sam_dedup_cm3, dependencies=sam_dedup_task['sam_sort'])
            sam_dedup_task['sam_dedup'] = self.addTask('sam_dedup', self.sam_dedup_cm4, dependencies=sam_dedup_task['sam_add_RG'])
            sam_dedup_task['sam_index'] = self.addTask('sam_index', self.sam_dedup_cm5, dependencies=sam_dedup_task['sam_dedup'])
            sam_dedup_task['align_stat'] = self.addTask('align_stat', self.sam_dedup_cm6, dependencies=sam_dedup_task['sam_add_RG'])
            sam_dedup_task['insert_stat'] = self.addTask('insert_stat', self.sam_dedup_cm7, dependencies=sam_dedup_task['sam_dedup'])
            sam_dedup_task['depth_stat'] = self.addTask('depth_stat', self.sam_dedup_cm8, dependencies=sam_dedup_task['sam_dedup'])
        
        realign_task = {'indel_list':None,'realign':None}
        if 'all' in self.stages or 'realign' in self.stages:
            realign_task['indel_list'] = self.addTask('indel_list', self.realign_cm1, dependencies=sam_dedup_task['sam_index'])
            realign_task['realign'] = self.addTask('realign', self.realign_cm2, dependencies=realign_task['indel_list'])
                
        vcf_task = {'vcf_vcf':None,'vcf_take_snp':None,'vcf_take_indel':None,'vcf_filter_snp':None,'vcf_filter_indel':None}
        if 'all' in self.stages or 'vcf' in self.stages:
            vcf_task['vcf_vcf'] = self.addTask('vcf_vcf', self.vcf_cm1, dependencies=realign_task['realign'])
            vcf_task['vcf_take_snp'] = self.addTask('vcf_take_snp', self.vcf_cm2, dependencies=vcf_task['vcf_vcf'])
            vcf_task['vcf_take_indel'] = self.addTask('vcf_take_indel', self.vcf_cm3, dependencies=vcf_task['vcf_vcf'])
            vcf_task['vcf_filter_snp'] = self.addTask('vcf_filter_snp', self.vcf_cm4, dependencies=vcf_task['vcf_take_snp'])
            vcf_task['vcf_filter_indel'] = self.addTask('vcf_filter_indel', self.vcf_cm5, dependencies=vcf_task['vcf_take_indel'])
        
        mp_task = None
        if 'all' in self.stages or 'mpileup' in self.stages:
            mp_task = self.addTask('mpileup', self.mp_command, dependencies=realign_task['realign'])
            
        het_task = {'het_raw':None, 'het_filter':None}
        if 'all' in self.stages or 'het' in self.stages:
            het_task['het_raw'] = self.addTask('het_raw', self.het_raw_cm, dependencies=mp_task)
            het_task['het_filter'] = self.addTask('het_filter', self.het_filter_cm, dependencies=het_task['het_raw'])
            
        rm_task = None
        if 'all_rm' in self.stages or 'rm' in self.stages:
            dep = None
            if het_task['het_filter'] is not None:
                dep = [het_task['het_filter'],sam_dedup_task['align_stat'],sam_dedup_task['insert_stat'],sam_dedup_task['depth_stat']]
            rm_task = self.addTask('rm', self.rm_command, dependencies=dep)
            #rm_task = self.addTask('rm',self.rm_command)
            
if __name__ == "__main__":
    usage = "usage: %prog [opstatistics of rawtions] \n\
            heteroplasmy pyflow\n\
            Oct 18 2017 \n\
            "

    #parser = OptionParser(usage=usage)
    parser = argparse.ArgumentParser(description = usage)
    
    parser.add_argument('--sample', type=str, required=True, help="the sample name")
    parser.add_argument('--fastq', type=str, required=True, nargs="+", help="the fastq file[s]")
    parser.add_argument('--out_dir', type=str, required=True, help="output directory")
    parser.add_argument('--ref', type=str, required=True, help="reference genome for mapping")
    parser.add_argument('--ref_mt', type=str, required=True, help="mtDNA reference genome")
    parser.add_argument('--config_file', type=str, required=True, help="config")
    parser.add_argument('--het_config', type=str, required=True, help="config file for heteroplasmy calling")
    parser.add_argument('--stages', type=str, required=True, nargs="+", help="stages to run")
    parser.add_argument('--isContinue', type=distutils.util.strtobool, default='False')
    parser.add_argument('--isDryRun', type=distutils.util.strtobool, default='False')
            
    opts = parser.parse_args()
    
    print opts
        
    heteroplasmy_wf = Heteroplasmy(
        opts.sample,
        opts.fastq,
        opts.out_dir,
        opts.ref,
        opts.ref_mt,
        opts.config_file,
        opts.het_config,
        opts.stages,
        
    )
    
    sys.exit(heteroplasmy_wf.run(mode='local',dataDirRoot=opts.out_dir,isContinue=opts.isContinue,isDryRun=opts.isDryRun))