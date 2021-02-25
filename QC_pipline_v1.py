#!/usr/bin/env python
import os
import sys
import argparse
import ConfigParser

def getConfig(section,key):
	config = ConfigParser.ConfigParser()
	path=os.path.abspath('./softwares.config')
	config.read(path)
	return config.get(section,key)

def parse_short_read_dir(inputs,outs,PE_SE='PE',fastp_p=''):
	input_path=os.path.abspath(inputs)+'/'
	out_path=os.path.abspath(outs)+'/'
	lst=os.listdir(input_path)
	parameters=fastp_p.strip()
	if PE_SE == 'PE':
		if lst[0].endswith('.fq.gz'):
			samples = [i.replace('_1.fq.gz','').replace('_2.fq.gz','') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,input2=input_path+sample+'_1.fq.gz',input_path+sample+'_2.fq.gz'
				output1,output2=out_path+sample+'_1.clean.fq.gz',out_path+sample+'_2.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
		elif lst[0].endswith('.fq'):
			samples = [i.replace('_1.fq','').replace('_2.fq','') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,input2=input_path+sample+'_1.fq',input_path+sample+'_2.fq'
				output1,output2=out_path+sample+'_1.clean.fq.gz',out_path+sample+'_2.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
		elif lst[0].endwith('.fastq.gz'):
			samples = [i.replace('_1.fastq.gz','').replace('_2.fastq.gz','') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,input2=input_path+sample+'_1.fastq.gz',input_path+sample+'_2.fastq.gz'
				output1,output2=out_path+sample+'_1.clean.fq.gz',out_path+sample+'_2.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
		elif lst[0].endwith('.fastq'):
			samples = [i.replace('_1.fastq','').replace('_2.fastq','') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,input2=input_path+sample+'_1.fastq',input_path+sample+'_2.fastq'
				output1,output2=out_path+sample+'_1.clean.fq.gz',out_path+sample+'_2.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
	elif sys.argv[3] == 'SE':
		if lst[0].endswith('.fq.gz'):
			samples = [i.replace('.fq.gz') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,output1 = input_path+sample+'.fq.gz',out_path+sample+'.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters)
		elif lst[0].endswith('.fq'):
			samples = [i.replace('.fq') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,output1 = input_path+sample+'.fq',out_path+sample+'.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters)
		elif lst[0].endwith('.fastq.gz'):
			samples = [i.replace('.fastq.gz') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,output1 = input_path+sample+'.fastq.gz',out_path+sample+'.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters)
	 	elif lst[0].endwith('.fastq'):
			samples = [i.replace('.fastq') for i in lst]
			samples = sorted(samples)
			for sample in samples:
				input1,output1 = input_path+sample+'.fastq',out_path+sample+'.clean.fq.gz'
				html,json=out_path+sample+'.html',out_path+sample+'.json'
				fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters)

def parse_pacbio_read_dir(inputs,outs,seq_type='Sequel',omic_type='DNA',mode='CCS',illumina='',ccs_p='',lima_p='',lordec_p='',isoseq3=';--verbose;'):
	input_path=os.path.abspath(inputs)+'/'
	out_path=os.path.abspath(outs)+'/'
	lst=os.listdir(input_path)
	isoseq3_refine,isoseq3_cluster,isoseq3_polish=isoseq3.split(';')
	if seq_type=='Sequel':
		samples = [i if i.endswith('.subreads.bam') for i in lst] # or .subreadset.xml
		if omic_type=='DNA':
			for sample in samples:
				input1=input_path+sample
				fw=open(sample+'.dna.ccs.sh','a')
				fw.write('#!/usr/bin/bash\n{lima} {input1} barcoded_primers.fasta {sample}.demux.subreads.bam {lima_p}\n
for i in {sample}*.demux.primer_5p--primer_3p.bam;do {ccs} $i $i.ccs.bam {ccp_p}\n
for i in {sample}*.demux.primer_5p--primer_3p.bam.ccs.bam ;do {bam2fasta} $i $i.pacbio.fasta'.format(lima=lima,input1=input1,sample=sample,lima_p=lima_p,ccs=ccs,ccs_p=ccs_p,bam2fasta=bam2fasta))
				if mode=='CCS':continue
				elif mode=='CLR':
					fw.write('for i in {sample}*.pacbio.fasta;do {} -2 {} -i $i -o $i.corrected.fasta {}'.format(lordec_correct,illumina,sample,sample,lordec_correct_p)

		elif omic_type=='RNA':
			for sample in samples:
				input1=input_path+sample
				fw=open(sample+'.rna.ccs.sh','a')
				fw.write('#!/usr/bin/bash\n{ccs} {input1} {sample}.ccs.bam {ccs_p}\n
{lima} {sample}.ccs.bam primers.fasta {sample}.demux.ccs.bam {lima_p}\n
{isoseq3} refine {sample}.demux.primer_5p--primer_3p.bam primers.fasta {sample}.flnc.bam {isoseq3_p_refine}\n
{isoseq3} cluster {sample}.flnc.bam {sample}.unpolished.bam {isoseq3_p_cluster}\n
{isoseq3} polish {sample}.unpolished.bam {input1} {sample}.polished.bam {isoseq3_p_polish}\n
{bam2fasta} {sample}.polished.bam {sample}.pacbio.fasta
'.format(ccs=ccs,input1=input1,sample=sample,ccs_p=ccs_p,lima=lima,lima_p=lima_p,isoseq3=isoseq3,isoseq3_p_refine=isoseq3_p_refine,isoseq3_p_cluster=isoseq3_p_cluster,isoseq3_p_polish=isoseq3_p_polish,bam2fasta=bam2fasta))
				if mode=='CCS':continue  ## not corrected
				elif mode=='CLR':        ## corrected by illumina
					fw.write('{} -2 {} -i {}.pacbio.fasta -o {}.pacbio.corrected.fasta {}'.format(lordec_correct,illumina,sample,sample,lordec_correct_p)
	elif seq_type=='RS':
		samples = [i.replace('.1.bax.h5').replace('.2.bax.h5').replace('.3.bax.h5') if i.endswith('.bax.h5') for i in lst]  ## one bas.h5 file and three bax.h5 files in each sample run
		samples = sorted(samples)
		if omic_type=='DNA':
			for sample in samples:
				input1,input2,input3=input_path+sample+'.1.bax.h5',input_path+sample+'.2.bax.h5',input_path+sample+'.3.bax.h5'
				fh=open(sample+'.dna.ccs.sh','a')
				fw.write('#!/usr/bin/bash\n{} {} {} {}\n'.format(bax2bam,input1,input2,input3))    ## bax2bam m.1.bax.h5 m.2.bax.h5 m.3.bax.h5 --> m.subreads.bam	
				fw.write('#!/usr/bin/bash\n{lima} {input1} barcoded_primers.fasta {sample}.demux.subreads.bam {lima_p}\n
for i in {sample}*.demux.primer_5p--primer_3p.bam;do {ccs} $i $i.ccs.bam {ccp_p}\n
for i in {sample}*.demux.primer_5p--primer_3p.bam.ccs.bam ;do {bam2fasta} $i $i.pacbio.fasta'.format(lima=lima,input1=input1,sample=sample,lima_p=lima_p,ccs=ccs,ccs_p=ccs_p,bam2fasta=bam2fasta))
				if mode=='CCS':continue
				elif mode=='CLR':
					fw.write('for i in {sample}*.pacbio.fasta;do {} -2 {} -i $i -o $i.corrected.fasta {}'.format(lordec_correct,illumina,sample,sample,lordec_correct_p)
		elif omic_type=='RNA':
			for sample in samples:
				input1,input2,input3=input_path+sample+'.1.bax.h5',input_path+sample+'.2.bax.h5',input_path+sample+'.3.bax.h5'
				fh=open(sample+'.dna.ccs.sh','a')
				fw.write('#!/usr/bin/bash\n{ccs} {input1} {sample}.ccs.bam {ccs_p}\n
{lima} {sample}.ccs.bam primers.fasta {sample}.demux.ccs.bam {lima_p}\n
{isoseq3} refine {sample}.demux.primer_5p--primer_3p.bam primers.fasta {sample}.flnc.bam {isoseq3_p_refine}\n
{isoseq3} cluster {sample}.flnc.bam {sample}.unpolished.bam {isoseq3_p_cluster}\n
{isoseq3} polish {sample}.unpolished.bam {input1} {sample}.polished.bam {isoseq3_p_polish}\n
{bam2fasta} {sample}.polished.bam {sample}.pacbio.fasta
'.format(ccs=ccs,input1=input1,sample=sample,ccs_p=ccs_p,lima=lima,lima_p=lima_p,isoseq3=isoseq3,isoseq3_p_refine=isoseq3_p_refine,isoseq3_p_cluster=isoseq3_p_cluster,isoseq3_p_polish=isoseq3_p_polish,bam2fasta=bam2fasta))
				if mode=='CCS':continue  ## not corrected
				elif mode=='CLR':        ## corrected by illumina
					fw.write('{} -2 {} -i {}.pacbio.fasta -o {}.pacbio.corrected.fasta {}'.format(lordec_correct,illumina,sample,sample,lordec_correct_p)

def parse_nanopore_read_dir():
	


if __name__ == '__main__':
#### Parse arguments ####
	parser = argparse.ArgumentParser(description='QC pipline v1.0')
	parser.add_argument('inputs',type=str,
                    help='The input directiory(rawdata),the suffixes of short reads in the directory must be _1/2.fq or _1/2.fq.gz or _1/2.fastq or _1/2.fastq.gz;the suffixes of pacbio long reads must be .subreads.bam(Sequel platform) or .1.bax.h5,.2.bax.h5,.3.bax.h5(RS/RSII platform);the suffiexes of nanopore long reads must be .fastq')
	parser.add_argument('outputs',type=str,
                    help='The output directiory(cleandata and report)')

	parser.add_argument('-c','--choice',type=str,default='illumina',choices=['illumina', 'pacbio', 'nanopore'],
                    help='Choose a sequencing machine,illumina for short reads,pacbio and nanopore for long reads,default is illumina')
	parser.add_argument('-s','--seq_type',type=str,default='PE',
                    help='Paired end(PE) or singe end(SE) short reads in illumina sequencing platform,default is PE')
	parser.add_argument('-omic','--omic_type',type=str,default='DNA',
                    help='The omics type:DNA(WGS) or RNA(transcipt)')
	parser.add_argument('-fastp_p','--fastp_parameters',type=str,default='',
                    help='The parameters for fastq softwares,e.g: -fastp_p "-w 1 --html fastp.html"')
	parser.add_argument('-ccs_p','--ccs_parameters',type=str,default='',
                    help='The parameters for SMRTlinks ccs tools')
	parser.add_argument('-lima_p','--lima_parameters',type=str,default='',
                    help='The parameters for SMRTlinks lima tools')
	parser.add_argument('-isoseq3_p','--isoseq3_parameters',type=str,default='',
                    help='The parameters for SMRTlinks isoseq3 tools')
	parser.add_argument('-lordec_p','--lordec_correct_parameters',type=str,default='',
                    help='The parameters for lordec_correct tools')
	parser.add_argument('-nanoplot_p','--nanoplot_parameters',type=str,default='',
                    help='The parameters for Nanoplot software')
	args = parser.parse_args()
	if args.choice=='pacbio' or args.choice =='nanopore' and args.seq_type:
		print 'Long reads have no PE/SE type,please choose illumina'
		exit()
##### Parse software config ####
	fastp=getConfig('QC','fastp')
	ccs=getConfig('QC','ccs')
	lima=getConfig('QC','lima')
	bax2bam=getConfig('QC','bax2bam')
	bam2fastq=getConfig('QC','bam2fastq')
	isoseq3=getConfig('QC','isoseq3')
	lordec_correct=getConfig('QC','lordec_correct')
	lordec_trim=getConfig('QC','lordec_trim')
	nanoplot=getConfig('QC','nanoplot')
#### Main ####
	if args.choice=='illumina':
		if args.omic_type=='DNA':
			parse_short_read_dir(args.inputs,args.outputs,PE_SE='PE',fastp_p='')
		elif args.omic_type=='RNA': ## remove rRNA first
			input_path=os.path.abspath(args.inputs)+'/'
			lst=os.listdir(input_path)
			
			
	elif args.choice=='pacbio':
		parse_pacbio_read_dir()

	elif args.choice=='nanopore':
		parse_nanopore_read_dir()

