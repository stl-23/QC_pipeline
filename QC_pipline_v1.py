#!/usr/bin/env python
import os
import sys
import argparse
import ConfigParser


def getConfig(section, key):
    config = ConfigParser.ConfigParser()
    path = os.path.abspath('./softwares.config')
    config.read(path)
    return config.get(section, key)


def remove_rRNA(sample, rRNA_out, rRNA_db, PE_SE='PE', hisat_p='', *samples):
    outname = rRNA_out + sample + '.gz'
    out = ''
    if PE_SE == 'PE':
        input1, input2 = samples[0], samples[1]
        out = 'hisat2 -x {} -1 {} -2 {} -S {}.sam --un-conc-gz {} {}'.format(rRNA_db, input1, input2, sample, outname,
                                                                             hisat_p)
    elif PE_SE == 'SE':
        input1 = samples[0]
        out = 'hisat2 -x {} -U {} -S {}.sam --un-gz {} {}'.format(rRNA_db, input1, sample, outname, hisat_p)
    return out


def parse_short_read_dir(inputs, outs, PE_SE='PE'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    # parameters=fastp_p.strip()
    samples = []
    input1 = []
    input2 = []
    output1 = []
    output2 = []
    html = []
    json = []
    if PE_SE == 'PE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('_1.fq.gz', '').replace('_2.fq.gz', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                input2.append(input_path + sample + '_2.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
        elif lst[0].endswith('.fq'):
            samples = [i.replace('_1.fq', '').replace('_2.fq', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                input2.append(input_path + sample + '_2.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
        elif lst[0].endwith('.fastq.gz'):
            samples = [i.replace('_1.fastq.gz', '').replace('_2.fastq.gz', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                input2.append(input_path + sample + '_2.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
        elif lst[0].endwith('.fastq'):
            samples = [i.replace('_1.fastq', '').replace('_2.fastq', '') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                input2.append(input_path + sample + '_2.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -I {} -o {} -O {} {} --html {} --json {}'.format(fastp,input1,input2,output1,output2,html,json,parameters))
    elif sys.argv[3] == 'SE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('.fq.gz') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters))
        elif lst[0].endswith('.fq'):
            samples = [i.replace('.fq') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters))
        elif lst[0].endwith('.fastq.gz'):
            samples = [i.replace('.fastq.gz') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters))
        elif lst[0].endwith('.fastq'):
            samples = [i.replace('.fastq') for i in lst]
            samples = sorted(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        #	fh=open(sample+'.sh','w').write('{} -i {} -o {} {} --html {} --json {}'.format(fastp,input1,output1,html,json,parameters))
    return samples, input1, input2, output1, output2, html, json

#def barcode_bam(raw_bam, barcode_file,lima):   ### https://www.biostars.org/p/349068/
    #fw = open('demultplexing.sh', 'w')
#    out = '{lima} {raw_bam} {barcode} --split-bam'.format(lima=lima,raw_bam=raw_bam,barcode=barcode_file)
#    return out

def parse_pacbio_read_dir(inputs, outs, seq_type='Sequel', omic_type='DNA', mode='CCS', correction='', samtools,Cal_subreads,lima, ccs,
                          lordec, isoseq3, bax2bam, bam2fasta, ccs_p='',lima_p='', lordec_p='', isoseq3_p=';--verbose;'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    isoseq3_refine, isoseq3_cluster, isoseq3_polish = isoseq3_p.split(';')
    if seq_type == 'Sequel':
        samples = [i for i in lst if i.endswith('.subreads.bam')]  # or .subreadset.xml
        if omic_type == 'DNA':
            for sample in samples:
                input1 = input_path + sample
                output1 = out_path + sample
                fw = open(sample + '.dna.ccs.sh', 'a')
                fw.write(f'#!/usr/bin/bash\n{samtools} stats {input1}.subreads.bam > {output1}.stat.tmp\n\
                python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list\n\
                Rscript {read_distributon_hist} {output1}.stat.list {output1}\n\
                {ccs} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}\n\
                {bam2fasta} {output1}.ccs.bam {output1}.pacbio.fasta'.format(ccs=ccs,
                                                                                    samtools=samtools,
                                                                                    read_distributon_hist=read_distributon_hist,
                                                                                    Cal_subreads=Cal_subreads,
                                                                                    input1=input1,
                                                                                    output1=output1,
                                                                                    ccs_p=ccs_p,
                                                                                    bam2fasta=bam2fasta))
                if mode == 'CCS':
                    continue
                elif mode == 'CLR':
                    fw.write('{} -2 {} -i {}.pacbio.fasta -o {}.corrected.fasta {}'.format(lordec,output1,output1,correction,lordec_p))

        elif omic_type == 'RNA':
            for sample in samples:
                input1 = input_path + sample
                output1 = out_path + sample
                fw = open(sample + '.rna.ccs.sh', 'a')
                fw.write('#!/usr/bin/bash\n{samtools} stats {input1}.subreads.bam > {output1}.stat.tmp\n\
python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list\n\
Rscript {read_distributon_hist} {output1}.stat.list {output1}\n\
{ccs} {input1} {output1}.ccs.bam {ccs_p}\n\
{lima} {output1}.ccs.bam primers.fasta {output1}.demux.ccs.bam {lima_p}\n\
{isoseq3} refine {output1}.demux.primer_5p--primer_3p.bam primers.fasta {output1}.flnc.bam {isoseq3_p_refine}\n\
{isoseq3} cluster {output1}.flnc.bam {output1}.unpolished.bam {isoseq3_p_cluster}\n\
{isoseq3} polish {output1}.unpolished.bam {input1} {output1}.polished.bam {isoseq3_p_polish}\n\
{bam2fasta} {output1}.polished.bam {output1}.pacbio.fasta\
'.format(ccs=ccs, input1=input1, output1=output1, ccs_p=ccs_p, lima=lima, lima_p=lima_p, isoseq3=isoseq3,
         read_distributon_hist=read_distributon_hist,samtools=samtools,Cal_subreads=Cal_subreads,
         isoseq3_p_refine=isoseq3_p_refine, isoseq3_p_cluster=isoseq3_p_cluster, isoseq3_p_polish=isoseq3_p_polish,
         bam2fasta=bam2fasta))
                if mode == 'CCS':
                    continue  ## not corrected
                elif mode == 'CLR':  ## corrected by illumina
                    fw.write('{} -2 {} -i {}.pacbio.fasta -o {}.corrected.fasta {}'.format(lordec,output1,output1,correction,lordec_p))

    elif seq_type == 'RS':
        samples = [i.replace('.1.bax.h5').replace('.2.bax.h5').replace('.3.bax.h5') for i in lst if
                   i.endswith('.bax.h5')]  ## one bas.h5 file and three bax.h5 files in each sample run
        samples = sorted(samples)
        if omic_type == 'DNA':
            for sample in samples:
                input1, input2, input3 = input_path + sample + '.1.bax.h5', input_path + sample + '.2.bax.h5', input_path + sample + '.3.bax.h5'
                output1 = out_path + sample
                fw = open(sample + '.dna.ccs.sh', 'a')
                fw.write('#!/usr/bin/bash\n{} {} {} {} && mv {}.subreads.bam {}\n'.format(bax2bam, input1, input2,
                                                                 input3, sample, out_path))  ## bax2bam m.1.bax.h5 m.2.bax.h5 m.3.bax.h5 --> m.subreads.bam
                fw.write('{samtools} stats {output1}.subreads.bam > {output1}.stat.tmp\n\
                python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list\n\
                Rscript {read_distributon_hist} {output1}.stat.list {output1}\n\
                {ccs} {output1}.subreads.bam {output1}.ccs.bam {ccs_p}\n\
                {bam2fasta} {output1}.ccs.bam {output1}.pacbio.fasta'.format(ccs=ccs,
                                                                             samtools=samtools,
                                                                             read_distributon_hist=read_distributon_hist,
                                                                             Cal_subreads=Cal_subreads,
                                                                             output1=output1,
                                                                             ccs_p=ccs_p,
                                                                             bam2fasta=bam2fasta))
                if mode == 'CCS':
                    continue
                elif mode == 'CLR':
                    fw.write('{} -2 {} -i {}.pacbio.fasta -o {}.corrected.fasta {}'.format(lordec,output1,output1,correction,lordec_p))

        elif omic_type == 'RNA':
            for sample in samples:
                input1, input2, input3 = input_path + sample + '.1.bax.h5', input_path + sample + '.2.bax.h5', input_path + sample + '.3.bax.h5'
                output1 = out_path + sample
                fw = open(sample + '.rna.ccs.sh', 'a')
                fw.write('#!/usr/bin/bash\n{} {} {} {} && mv {}.subreads.bam {}\n'.format(bax2bam, input1, input2,
                                                                                          input3, sample, out_path))
                fw.write('{samtools} stats {output1}.subreads.bam > {output1}.stat.tmp\n\
python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list\n\
Rscript {read_distributon_hist} {output1}.stat.list {output1}\n\
{ccs} {output1}.subreads.bam {output1}.ccs.bam {ccs_p}\n\
{lima} {output1}.ccs.bam primers.fasta {output1}.demux.ccs.bam {lima_p}\n\
{isoseq3} refine {output1}.demux.primer_5p--primer_3p.bam primers.fasta {output1}.flnc.bam {isoseq3_p_refine}\n\
{isoseq3} cluster {output1}.flnc.bam {output1}.unpolished.bam {isoseq3_p_cluster}\n\
{isoseq3} polish {output1}.unpolished.bam {input1} {output1}.polished.bam {isoseq3_p_polish}\n\
{bam2fasta} {output1}.polished.bam {output1}.pacbio.fasta\
'.format(ccs=ccs, input1=input1, output1=output1, ccs_p=ccs_p, lima=lima, lima_p=lima_p, isoseq3=isoseq3,
         read_distributon_hist=read_distributon_hist, samtools=samtools, Cal_subreads=Cal_subreads,
         isoseq3_p_refine=isoseq3_p_refine, isoseq3_p_cluster=isoseq3_p_cluster, isoseq3_p_polish=isoseq3_p_polish,
         bam2fasta=bam2fasta))
                if mode == 'CCS':
                    continue  ## not corrected
                elif mode == 'CLR':  ## corrected by illumina
                    fw.write('{} -2 {} -i {}.pacbio.fasta -o {}.corrected.fasta {}'.format(lordec,output1,output1,correction,lordec_p))


def parse_nanopore_read_dir(inputs, outputs, nanoplot_p='--plots hex dot pauvre kde'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    samples = [i.replace('.fastq.gz') for i in lst if i.endswith('.fastq.gz')]
    samples = sorted(samples)
    for sample in samples:
        input_sample = input_path + sample
        fw = open(sample + '.nanoplot.sh', 'a')
        fw.write('{} --fastq {}.fastq.gz -o {} {}'.format(nanoplot, input_sample, out_path, nanoplot_p))


if __name__ == '__main__':
    #### Parse arguments ####
    parser = argparse.ArgumentParser(description='QC pipline v1.0')
    parser.add_argument('inputs', type=str,
                        help='The input directiory(rawdata),the suffixes of short reads in the directory must be _1/2.fq or _1/2.fq.gz or _1/2.fastq or _1/2.fastq.gz;the suffixes of pacbio long reads must be .subreads.bam(Sequel platform) or .1.bax.h5,.2.bax.h5,.3.bax.h5(RS/RSII platform);the suffiexes of nanopore long reads must be .fastq.gz')
    parser.add_argument('outputs', type=str,
                        help='The output directiory(cleandata and report)')

    parser.add_argument('-c', '--choice', type=str, default='illumina', choices=['illumina', 'pacbio', 'nanopore'],
                        help='Choose a sequencing machine,illumina for short reads,pacbio and nanopore for long reads,default is illumina')
    parser.add_argument('-cr', '--correction', type=str, default='',
                        help='The illumina reads for pacbio long reads correction(used in CLR mode,the self-correction is enough for CCS mode in normal condition),FASTA/Q file(s),e.g. -cr "reads1.fa,reads2.fq,reads3.fq.gz" or -cr "read.fa"')
    parser.add_argument('-s', '--seq_type', type=str, default='PE',
                        help='Paired end(PE) or singe end(SE) short reads in illumina sequencing platform,default is PE')
    parser.add_argument('-omic', '--omic_type', type=str, default='DNA',
                        help='The omics type:DNA(genome) or RNA(transciptome),default is DNA')
    parser.add_argument('-fastp_p', '--fastp_parameters', type=str,
                        help='The parameters for fastp softwares,e.g: -fastp_p "-w 1", the parameters:--html and --json are defaultly set and named by the samples names to avoid overlapping in this pipline,please do not set again')
    parser.add_argument('-ccs_p', '--ccs_parameters', type=str, default='',
                        help='The parameters for SMRTlinks ccs tools,defalut parameters in DNA,"-noPolish --minPasses 1" in RNA')
    parser.add_argument('-lima_p', '--lima_parameters', type=str, default='',
                        help='The parameters for SMRTlinks lima tools,default parameters in DNA,"--isoseq --no-pbi" in RNA')
    parser.add_argument('-isoseq3_p', '--isoseq3_parameters', type=str, default=';--verbose;',
                        help='The parameters for SMRTlinks isoseq3 tools(refine/cluster/polish),seperated by semicolon,e.g. -isoseq3_p "-j 1;--verbose;-r 0.99",default is ";--verbose;"')
    parser.add_argument('-lordec_p', '--lordec_correct_parameters', type=str, default='',
                        help='The parameters for lordec_correct tools,"-k 19 -s 3" in DNA,default parameters in RNA')
    parser.add_argument('-nanoplot_p', '--nanoplot_parameters', type=str, default='',
                        help='The parameters for Nanoplot software')
    args = parser.parse_args()
    input_dir = args.inputs
    output_dir = args.outputs
    platform_choice = args.choice
    correct_dir = args.correction
    seq_type = args.seq_type
    omic_type = args.omic_type
    fastp_parameters = args.fastp_parameters
    ccs_parameters = args.ccs_parameters
    lima_parameters = args.lima_parameters
    isoseq3_parameters = args.isoseq3_parameters
    lordec_correct_parameters = args.lordec_correct_parameters
    nanoplot_parameters = args.nanoplot_parameters
    if seq_type and platform_choice == 'pacbio' or platform_choice == 'nanopore':
        print('Long reads have no PE/SE type,please choose illumina')
        exit()
    ##### Parse software config ####
    fastp = getConfig('QC', 'fastp')
    ccs = getConfig('QC', 'ccs')
    lima = getConfig('QC', 'lima')
    bax2bam = getConfig('QC', 'bax2bam')
    bam2fastq = getConfig('QC', 'bam2fastq')
    isoseq3 = getConfig('QC', 'isoseq3')
    lordec_correct = getConfig('QC', 'lordec_correct')
    lordec_trim = getConfig('QC', 'lordec_trim')
    nanoplot = getConfig('QC', 'nanoplot')
    rRNA_data = getConfig('Database', 'rRNA')
    cal_subreads = getConfig('QCscripts', 'Cal_subreads')
    #### Main ####
    if platform_choice == 'illumina':
        samples, input1, input2, output1, output2, html, json = parse_short_read_dir(input_dir, output_dir,
                                                                                     PE_SE=seq_type)
        if omic_type == 'DNA':
            if seq_type == 'PE':
                for sample, index in enumerate(samples):
                    fw = open(sample + '.sh', 'w').write(
                        '{} -i {} -I {} -o {} -O {} --html {} --json {} {}'.format(fastp, input1[index], input2[index],
                                                                                   output1[index], output2[index],
                                                                                   html[index], json[html],
                                                                                   fastp_parameters))
            elif seq_type == 'SE':
                for sample, index in enumerate(samples):
                    fw = open(sample + '.sh', 'w').write(
                        '{} -i {} -o {} --html {} --json {} {}'.format(fastp, input1[index], output1[index],
                                                                       html[index], json[index], fastp_parameters))
        elif omic_type == 'RNA':  ## remove rRNA first
            os.system('mkdir -p {}/rRNA_dir/'.format(output_dir))
            tmp_dir = output_dir + '/rRNA_dir'
            if seq_type == 'PE':
                for sample, index in enumerate(samples):
                    out_file1 = remove_rRNA(sample, tmp_dir, rRNA_data, seq_type, '', input1[index], input2[index])
                    out_file2 = '{} -i {} -I {} -o {} -O {} --html {} --json {} {}'.format(fastp, input1[index],
                                                                                           input2[index],
                                                                                           output1[index],
                                                                                           outpt2[index], html[index],
                                                                                           json[html], fastp_parameters)
                    fw = open(sample + '.sh', 'w').write(out_file1 + '\n' + out_file2)
            elif seq_type == 'SE':
                for sample, index in enumerate(samples):
                    out_file1 = remove_rRNA(sample, tmp_dir, rRNA_data, seq_type, '', input1[index], input2[index])
                    out_file2 = '{} -i {} -o {} --html {} --json {} {}'.format(fastp, input1[index], output1[index],
                                                                               html[index], json[index],
                                                                               fastp_parameters)
                    fw = open(sample + '.sh', 'w').write(out_file1 + '\n' + out_file2)
    elif platform_choice == 'pacbio':
        if omic_type == 'DNA':
            if isoseq3_parameters == '':
                if lordec_correct_parameters == '':
                    lordec_correct_parameters = '-k 19 -s 3'
                parse_pacbio_read_dir(input_dir, output_dir, seq_type, omic_type, mode, correction, ccs_parameters,
                                      lima_parameters, lordec_correct_parameters, isoseq3_parameters)
            else:
                print('iso-seq3 parameters set in RNA analysis not DNA')
                exit()
        elif omic_type == 'RNA':
            if ccs_parameters == '':
                ccs_parameters = '-noPolish --minPasses 1'
            if lima_parameters == '':
                lima_parameters = '--isoseq --no-pbi'
            parse_pacbio_read_dir(input_dir, output_dir, seq_type, omic_type, mode, correction, ccs_parameters,
                                  lima_parameters, lordec_correct_parameters, isoseq3_parameters)

    elif platform_choice == 'nanopore':
        parse_nanopore_read_dir(input_dir, output_dir, nanoplot_parameters)

