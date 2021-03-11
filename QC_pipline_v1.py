#!/usr/bin/env python
import os
import argparse
import configparser
import time
#import QC.Cal_subreads as Cal_subreads



def getConfig(section, key):
    config = configparser.ConfigParser()
    path = os.path.abspath(os.path.dirname(__file__))+'/softwares.config'
    config.read(path)
    return config.get(section, key)


def remove_rRNA(sample, output_dir, rRNA_db, seq_type='PE', hisat_p='', *samples):
    rRNA_out = output_dir + 'rm_rRNA_dir/'
    out = ''
    if seq_type == 'PE':
        input1, input2 = samples[0], samples[1]
        out = """mkdir -p {rRNA_out}
        hisat2 -x {rRNA_db} -1 {input1} -2 {input2} -S {rRNA_out}/{sample}.sam --un-conc-gz {rRNA_out} {hisat_p}
        mv {rRNA_out}/un-conc-mate.1 {rRNA_out}/{sample}.1.rm.fq.gz 
        mv {rRNA_out}/un-conc-mate.2 {rRNA_out}/{sample}.2.rm.fq.gz
        """.format(rRNA_db=rRNA_db, input1=input1, input2=input2, sample=sample, rRNA_out=rRNA_out, hisat_p=hisat_p)
    elif seq_type == 'SE':
        input1 = samples[0]
        out = """mkdir -p {rRNA_out}
        hisat2 -x {rRNA_db} -U {input1} -S {rRNA_out}/{sample}.sam --un-gz {rRNA_out}/{sample}.rm.fq.gz {hisat_p}
        """.format(rRNA_db=rRNA_db, input1=input1, sample=sample, rRNA_out=rRNA_out, hisat_p=hisat_p)
    return out


def parse_short_read_dir(inputs, outs, seq_type='PE'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    samples = []
    input1 = []
    input2 = []
    output1 = []
    output2 = []
    html = []
    json = []
    seq_suffix = ['.fq', '.fq.gz', '.fastq', '.fastq.gz']
    lst = [file for file in lst for suffix in seq_suffix if file.endswith(suffix)]
    if not lst:
        print("No such file or directory or wrong file format,must be .fq/.fq.gz/.fastq/.fastq.gz")
        exit()
    if seq_type == 'PE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('_1.fq.gz', '').replace('_2.fq.gz', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                input2.append(input_path + sample + '_2.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        elif lst[0].endswith('.fq'):
            samples = [i.replace('_1.fq', '').replace('_2.fq', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                input2.append(input_path + sample + '_2.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        elif lst[0].endswith('.fastq.gz'):
            samples = [i.replace('_1.fastq.gz', '').replace('_2.fastq.gz', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                input2.append(input_path + sample + '_2.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        elif lst[0].endswith('.fastq'):
            samples = [i.replace('_1.fastq', '').replace('_2.fastq', '') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                input2.append(input_path + sample + '_2.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                output2.append(out_path + sample + '_2.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
    elif seq_type == 'SE':
        if lst[0].endswith('.fq.gz'):
            samples = [i.replace('.fq.gz','') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        elif lst[0].endswith('.fq'):
            samples = [i.replace('.fq','') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        elif lst[0].endswith('.fastq.gz'):
            samples = [i.replace('.fastq.gz','') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq.gz')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
        elif lst[0].endswith('.fastq'):
            samples = [i.replace('.fastq','') for i in lst]
            samples = set(samples)
            for sample in samples:
                input1.append(input_path + sample + '_1.fastq')
                output1.append(out_path + sample + '_1.clean.fq.gz')
                html.append(out_path + sample + '.html')
                json.append(out_path + sample + '.json')
    return samples, input1, input2, output1, output2, html, json

#def barcode_bam(raw_bam, barcode_file,lima):   ### https://www.biostars.org/p/349068/
#    fw = open('demultplexing.sh', 'w')
#    out = '{lima} {raw_bam} {barcode} --split-bam'.format(lima=lima,raw_bam=raw_bam,barcode=barcode_file)
#    return out

def parse_pacbio_read_dir(inputs, outs, mt_type='Sequel'):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    input1=[]
    input2=[]
    input3=[]
    output1=[]
    samples=[]
    seq_suffix = ['.subreads.bam','1.bax.h5','.2.bax.h5','.3.bax.h5']
    lst = [file for file in lst for suffix in seq_suffix if file.endswith(suffix)]
    if not lst:
        print("No such file or directory or wrong file format,must be .subreads.bam or .1/2/3.bax.h5")
        exit()
    if mt_type == 'Sequel':
        samples = [i.replace('.subreads.bam','') for i in lst if i.endswith('.subreads.bam')]  # or .subreadset.xml
        samples = set(samples)
        for sample in samples:
            input1.append(input_path + sample)
            output1.append(out_path + sample)
    elif mt_type == 'RS':
        samples = [i.replace('.1.bax.h5','').replace('.2.bax.h5','').replace('.3.bax.h5','') for i in lst if i.endswith('.bax.h5')]  ## one bas.h5 file and three bax.h5 files in each sample run
        samples = set(samples)
        for sample in samples:
            input1.append(input_path + sample + '.1.bax.h5')
            input2.append(input_path + sample + '.2.bax.h5')
            input3.append(input_path + sample + '.3.bax.h5')
            output1.append(out_path + sample)

    return samples,input1,input2,input3,output1

def parse_nanopore_read_dir(inputs, outs):
    input_path = os.path.abspath(inputs) + '/'
    out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    samples = [i.replace('.fastq.gz','') for i in lst if i.endswith('.fastq.gz')]
    samples = set(samples)
    input_samples = []
    for sample in samples:
        input_samples.append(input_path + sample + '.fastq.gz')
    return samples, input_samples


if __name__ == '__main__':
    #### Parse arguments ####
    examplelog="""EXAMPLES:
    python QC_pipline_v1.py /root/my_data/example_inputs/ /root/my_data/example_outputs/ -c illumina -omic DNA -s PE -fastp_p "-w 3"
    python QC_pipline_v1.py /root/my_data/example_inputs/ /root/my_data/example_outputs/ -c pacbio -omic RNA -mt Sequel -ccs_p "-noPolish --minPasses 1" -lima_p "--isoseq --no-pbi" -isoseq3_p ";--verbose;" -lordec_p "-m 2G"
    python QC_pipline_v1.py /root/my_data/example_inputs/ /root/my_data/example_outputs/ -c nanopore -omic DNA -nanoplot_p "--plots hex dot pauvre kde"
    """
    parser = argparse.ArgumentParser(description='QC pipline v1.0', epilog=examplelog, add_help=False)
    general = parser.add_argument_group(title='General options')
    general.add_argument('-inputs', type=str,
                         help='The input directory(rawdata),the suffixes of short reads in the directory must be _1/2.fq or _1/2.fq.gz or _1/2.fastq or _1/2.fastq.gz;the suffixes of pacbio long reads must be .subreads.bam(Sequel platform) or .1.bax.h5,.2.bax.h5,.3.bax.h5(RS/RSII platform);the suffiexes of nanopore long reads must be .fastq.gz')
    general.add_argument('-outputs', type=str,
                        help='The output directory')
    general.add_argument('-h', '--help', action="help",
                         help="show the help and exit")
    general.add_argument('-c', '--choice', type=str, default='illumina', choices=['illumina', 'pacbio', 'nanopore'],
                        help='Choose a sequencing machine,illumina for short reads,pacbio and nanopore for long reads,default is illumina')
    general.add_argument('-omic', '--omic_type', type=str, default='DNA',
                        help='The omics type:DNA(genome) or RNA(transciptome),default is DNA')
    illumina = parser.add_argument_group(title='Short reads options (illumina)')
    illumina.add_argument('-s', '--seq_type', type=str, default='PE',
                        help='Paired end(PE) or singe end(SE) short reads in illumina sequencing platform,default is PE')
    illumina.add_argument('-fastp_p', '--fastp_parameters', type=str,default='',
                        help='The parameters for fastp softwares,e.g: -fastp_p "-w 1", the parameters:--html and --json are defaultly set and named by the samples names to avoid overlapping in this pipline,please do not set again')
    pacbio = parser.add_argument_group(title='Long reads options (Pacbio)')
    pacbio.add_argument('-mt', '--machine_type', type=str,default='Sequel',
                        help='The sequencing platform:Sequel/RS,default is Sequel')
    pacbio.add_argument('-cr', '--correction', type=str, default='no', choices=['no', 'yes'],
                        help='Choose illumina short reads for pacbio long reads correction or not (used in CLR mode,the self-correction is enough for CCS mode in normal condition),FASTA/Q file(s),default is false')
    pacbio.add_argument('-ccs_p', '--ccs_parameters', type=str, default='',
                        help='The parameters for SMRTlinks ccs tools,defalut parameters in DNA,"-noPolish --minPasses 1" in RNA')
    pacbio.add_argument('-lima_p', '--lima_parameters', type=str, default='',
                        help='The parameters for SMRTlinks lima tools,default parameters in DNA,"--isoseq --no-pbi" in RNA')
    pacbio.add_argument('-isoseq3_p', '--isoseq3_parameters', type=str, default='',
                        help='The parameters for SMRTlinks isoseq3 tools(refine/cluster/polish),seperated by semicolon,e.g. -isoseq3_p "-j 1;--verbose;-r 0.99",default is ";--verbose;"')
    pacbio.add_argument('-lordec_p', '--lordec_correct_parameters', type=str, default='',
                        help='The parameters for lordec_correct tools,"-k 19 -s 3" in DNA,default parameters in RNA')
    nanopore = parser.add_argument_group(title='Long reads options (Nanopore)')
    nanopore.add_argument('-nanoplot_p', '--nanoplot_parameters', type=str, default='--plots hex dot',
                        help='The parameters for Nanoplot software')
    nanopore.add_argument('-nanofilt_p', '--nanofilt_parameters', type=str, default=' -q 7 -l 1000 --headcrop 50 --tailcrop 50',
                          help='The parameters for NanoFilt software')
    args = parser.parse_args()
    input_dir = args.inputs
    output_dir = args.outputs
    platform_choice = args.choice
    correct = args.correction
    seq_type = args.seq_type
    omic_type = args.omic_type
    mt_type = args.machine_type
    fastp_parameters = args.fastp_parameters
    ccs_parameters = args.ccs_parameters
    lima_parameters = args.lima_parameters
    isoseq3_parameters = args.isoseq3_parameters
    lordec_correct_parameters = args.lordec_correct_parameters
    nanoplot_parameters = args.nanoplot_parameters
    nanofilt_parameters = args.nanofilt_parameters

    seq_suffix = ['.ct.fa','.ct.fq','.ct.fa.gz','.ct.fq.gz','.ct.fasta','.ct.fasta.gz','.ct.fastq','.ct.fastq.gz']
    correct_data_lst = [os.path.abspath(input_dir)+'/'+ file for file in os.listdir(os.path.abspath(input_dir)) for suffix in seq_suffix if file.endswith(suffix)]

    #if seq_type and platform_choice == 'pacbio' or platform_choice == 'nanopore':
    #    print(seq_type)
    #    print('Long reads have no PE/SE type,please choose illumina')
    #    exit()
    #if platform_choice == 'pacbio' and omic_type == 'DNA' and isoseq3_parameters:
    #    print('Iso-seq3 parameters must be set in RNA analysis not DNA')
    #    exit()

    ##### Parse software and database config ####
    fastp = getConfig('QC', 'fastp').strip("'")
    ccs = getConfig('QC', 'ccs').strip("'")
    lima = getConfig('QC', 'lima').strip("'")
    bax2bam = getConfig('QC', 'bax2bam').strip("'")
    bam2fasta = getConfig('QC', 'bam2fasta').strip("'")
    isoseq3 = getConfig('QC', 'isoseq3').strip("'")
    lordec_correct = getConfig('QC', 'lordec_correct').strip("'")
    lordec_trim = getConfig('QC', 'lordec_trim').strip("'")
    nanoplot = getConfig('QC', 'nanoplot').strip("'")
    nanofilt = getConfig('QC','nanofilt').strip("'")
    samtools = getConfig('QC', 'samtools').strip("'")
    #rRNA_data = getConfig('Database', 'rRNA')
    #cal_subreads = getConfig('QCscripts', 'Cal_subreads')
    #read_distributon_hist = getConfig('Report','read_distributon_hist')
    #scripts_dir = os.path.abspath(os.path.dirname(__file__))+'/tmp_' + time.strftime("%Y-%m-%d_%H:%M:%S", time.localtime())
    #os.system('mkdir -p {}'.format(scripts_dir))
    rRNA_data = '/obs/gene-container-qc-test/database/rRNA/rRNAs'
    cal_subreads = os.path.abspath(os.path.dirname(__file__))+'/Cal_subreads.py'
    read_distributon_hist = os.path.abspath(os.path.dirname(__file__))+'/read_distributon_hist.R'
    #### Main ####
    if platform_choice == 'illumina':
        samples, input1, input2, output1, output2, html, json = parse_short_read_dir(input_dir, output_dir,
                                                                                     seq_type)
        if omic_type == 'DNA':
            if seq_type == 'PE':
                for index, sample in enumerate(samples):
                    fw = open(sample + '.sh', 'w').write(
                        '{} -i {} -I {} -o {} -O {} --html {} --json {} {}'.format(fastp, input1[index], input2[index],
                                                                                   output1[index], output2[index],
                                                                                   html[index], json[index],
                                                                                   fastp_parameters))
            elif seq_type == 'SE':
                for index, sample in enumerate(samples):
                    fw = open(sample + '.sh', 'w').write(
                        '{} -i {} -o {} --html {} --json {} {}'.format(fastp, input1[index], output1[index],
                                                                       html[index], json[index], fastp_parameters))
        elif omic_type == 'RNA':  ## remove rRNA first
            rm_rRNA_seq1 = []
            rm_rRNA_seq2 = []
            if seq_type == 'PE':
                for index, sample in enumerate(samples):
                    #os.system('mkdir -p {}/rm_rRNA_dir/{}'.format(output_dir,sample))
                    #tmp_dir = output_dir + '/rm_rRNA_dir/'+sample
                    rm_rRNA_seq1.append(output_dir + '/rm_rRNA_dir/'+sample+'/'+sample+'.1.rm.fq.gz')
                    rm_rRNA_seq2.append(output_dir + '/rm_rRNA_dir/'+sample+'/'+sample+'.2.rm.fq.gz')
                    out_file1 = remove_rRNA(sample, output_dir, rRNA_data, seq_type, '', input1[index], input2[index])
                    out_file2 = '{} -i {} -I {} -o {} -O {} --html {} --json {} {}'.format(fastp, rm_rRNA_seq1[index],
                                                                                           rm_rRNA_seq2[index],
                                                                                           output1[index],
                                                                                           output2[index], html[index],
                                                                                           json[index], fastp_parameters)
                    fw = open(sample + '.sh', 'w').write(out_file1 + '\n' + out_file2)
            elif seq_type == 'SE':
                for index, sample in enumerate(samples):
                    #os.system('mkdir -p {}/rm_rRNA_dir/{}'.format(output_dir, sample))
                    #tmp_dir = output_dir + '/rm_rRNA_dir/' + sample
                    rm_rRNA_seq1.append(output_dir + '/rm_rRNA_dir/' + sample + '/' + sample + '.rm.fq.gz')
                    out_file1 = remove_rRNA(sample, output_dir, rRNA_data, seq_type, '', input1[index])
                    out_file2 = '{} -i {} -o {} --html {} --json {} {}'.format(fastp, rm_rRNA_seq1[index], output1[index],
                                                                               html[index], json[index],
                                                                               fastp_parameters)
                    fw = open(sample + '.sh', 'w').write(out_file1 + '\n' + out_file2)
    elif platform_choice == 'pacbio':
        if omic_type == 'DNA':
            if lordec_correct_parameters == '':
                lordec_correct_parameters = "-k 19 -s 3"
            else:
                lordec_correct_parameters = lordec_correct_parameters
            if mt_type == 'Sequel':
                samples, input1, input2, input3, output1 = parse_pacbio_read_dir(input_dir, output_dir, mt_type)
                for index, sample in enumerate(samples):
                    out_file1 = """#!/usr/bin/bash
{samtools} stats {input1}.subreads.bam > {output1}.stat.tmp
python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list
Rscript {read_distributon_hist} {output1}.stat.list {output1}
{ccs} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}
{bam2fasta} {output1}.ccs.bam {output1}.pacbio.fasta
""".format(samtools=samtools,input1=input1[index],Cal_subreads=cal_subreads,output1=output1[index],
            read_distributon_hist=read_distributon_hist,ccs=ccs,ccs_p=ccs_parameters,
            bam2fasta=bam2fasta)
                    if correct == 'no':
                        out_file2 = ''
                    elif correct == 'yes':
                        if correct_data_lst:
                                correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input1[index])]
                                correct_data = ','.join(correct_data_lst_tmp)
                                out_file2 = """{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim,correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters)
                        else:
                            print('correct data not accessible!')
                            exit()
                    fw = open(sample+'.Sequel_DNA.sh', 'w').write(out_file1+'\n'+out_file2)
            elif mt_type == 'RS':
                samples, input1, input2, input3, output1 = parse_pacbio_read_dir(input_dir, output_dir, mt_type)
                for index, sample in enumerate(samples):
                    input_1, input_2, input_3 = input_dir + sample + '.1.bax.h5', input_dir + sample + '.2.bax.h5', input_dir + sample + '.3.bax.h5'
                    out_file1 = """#!/usr/bin/bash
{bax2bam} {input_1} {input_2} {input_3} && mv {sample}.subreads.bam {output_dir}
{samtools} stats {output1}.subreads.bam > {output1}.stat.tmp
python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list
Rscript {read_distributon_hist} {output1}.stat.list {output1}
{ccs} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}
{bam2fasta} {output1}.ccs.bam {output1}.pacbio.fasta
""".format(bax2bam=bax2bam, input_1=input_1, input_2=input_2, input_3=input_3, sample=sample, output_dir=output_dir,
            samtools = samtools, input1 = input1[index], Cal_subreads = cal_subreads, output1 = output1[index],
            read_distributon_hist = read_distributon_hist, ccs = ccs, ccs_p = ccs_parameters,
            bam2fasta=bam2fasta)
                    if correct == 'no':
                        out_file2 = ''
                    elif correct == 'yes':
                        if correct_data_lst:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input_dir+sample)]
                            correct_data = ','.join(correct_data_lst_tmp)
                            out_file2 = """{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim, correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters)
                        else:
                            print('correct data not accessible!')
                            exit()
                    fw = open(sample+'.RS_DNA.sh', 'w').write(out_file1 + '\n' + out_file2)
        elif omic_type == 'RNA':
            isoseq3_refine = ''
            isoseq3_cluster = ''
            isoseq3_polish = ''
            if ccs_parameters == '':
                ccs_parameters = "-noPolish --minPasses 1"
            else:
                ccs_parameters = ccs_parameters
            if lima_parameters == '':
                lima_parameters = "--isoseq --no-pbi"
            else:
                lima_parameters = lima_parameters
            if isoseq3_parameters == '':
                isoseq3_parameters = ";--verbose;"
            else:
                isoseq3_parameters = isoseq3_parameters
            isoseq3_refine, isoseq3_cluster, isoseq3_polish = isoseq3_parameters.strip().split(';')

            if mt_type == 'Sequel':
                samples, input1, input2, input3, output1 = parse_pacbio_read_dir(input_dir, output_dir, mt_type)
                for index, sample in enumerate(samples):
                    out_file1 = """#!/usr/bin/bash
{samtools} stats {input1}.subreads.bam > {output1}.stat.tmp
python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list
Rscript {read_distributon_hist} {output1}.stat.list {output1}
{ccs} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}
{lima} {output1}.ccs.bam primers.fasta {output1}.demux.ccs.bam {lima_p}
{isoseq3} refine {output1}.demux.primer_5p--primer_3p.bam primers.fasta {output1}.flnc.bam {isoseq3_refine}
{isoseq3} cluster {output1}.flnc.bam {output1}.unpolished.bam {isoseq3_cluster}
{isoseq3} polish {output1}.unpolished.bam {input1} {output1}.polished.bam {isoseq3_polish}
{bam2fasta} {output1}.polished.bam {output1}.pacbio.fasta
""".format(samtools=samtools,input1=input1[index],Cal_subreads=cal_subreads,output1=output1[index],
            read_distributon_hist=read_distributon_hist,ccs=ccs, ccs_p=ccs_parameters, lima=lima,
            lima_p=lima_parameters, isoseq3=isoseq3, isoseq3_refine=isoseq3_refine,
            isoseq3_cluster=isoseq3_cluster,isoseq3_polish=isoseq3_polish, bam2fasta=bam2fasta)
                    if correct == 'no':
                        out_file2 = ''
                    elif correct == 'yes':
                        if correct_data_lst:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input1[index])]
                            correct_data = ','.join(correct_data_lst_tmp)
                            out_file2 = """{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim,correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters)
                        else:
                            print('correct data not accessible!')
                            exit()
                    fw = open(sample+'.Sequel_RNA.sh', 'w').write(out_file1 + '\n' + out_file2)
            elif mt_type == 'RS':
                samples, input1, input2, input3, output1 = parse_pacbio_read_dir(input_dir, output_dir, mt_type)
                for index, sample in enumerate(samples):
                    input_1, input_2, input_3 = input_dir + sample + '.1.bax.h5', input_dir + sample + '.2.bax.h5', input_dir + sample + '.3.bax.h5'
                    out_file1 = """#!/usr/bin/bash
{bax2bam} {input_1} {input_2} {input_3} && mv {sample}.subreads.bam {output_dir}
{samtools} stats {output1}.subreads.bam > {output1}.stat.tmp
python {Cal_subreads} {output1}.stat.tmp {output1} > {output1}.stat.list
Rscript {read_distributon_hist} {output1}.stat.list {output1}
{ccs} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}
{lima} {output1}.ccs.bam primers.fasta {output1}.demux.ccs.bam {lima_p}
{isoseq3} refine {output1}.demux.primer_5p--primer_3p.bam primers.fasta {output1}.flnc.bam {isoseq3_refine}
{isoseq3} cluster {output1}.flnc.bam {output1}.unpolished.bam {isoseq3_cluster}
{isoseq3} polish {output1}.unpolished.bam {input1} {output1}.polished.bam {isoseq3_polish}
{bam2fasta} {output1}.polished.bam {output1}.pacbio.fasta
""".format(bax2bam=bax2bam, input_1=input_1, input_2=input_2, input_3=input_3, sample=sample, output_dir=output_dir,
            samtools = samtools, input1 = input1[index], Cal_subreads = cal_subreads, output1 = output1[index],
            read_distributon_hist = read_distributon_hist, ccs = ccs, ccs_p = ccs_parameters, lima=lima,
            lima_p=lima_parameters, isoseq3=isoseq3, isoseq3_refine=isoseq3_refine,
            isoseq3_cluster=isoseq3_cluster, isoseq3_polish=isoseq3_polish, bam2fasta=bam2fasta)
                    if correct == 'no':
                        out_file2 = ''
                    elif correct == 'yes':
                        if correct_data_lst:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input_dir+sample)]
                            print(correct_data_lst_tmp)
                            correct_data = ','.join(correct_data_lst_tmp)
                            out_file2 = """{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim,correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters)
                        else:
                            print('correct data not accessible!')
                            exit()
                    fw = open(sample+'.RS_RNA.sh', 'w').write(out_file1 + '\n' + out_file2)

    elif platform_choice == 'nanopore':
        samples, input_samples = parse_nanopore_read_dir(input_dir, output_dir)
        Before_filter = output_dir+'Before_filter'
        After_filter = output_dir+'After_filter'
        if omic_type == 'DNA' or omic_type == 'RNA':
            for index, sample in enumerate(samples):
                outfile1 = """mkdir -p {output_dir}/{BF} {output_dir}/{AF}
{nanoplot} --fastq {input_sample} -o {BF} {nanoplot_p}
gunzip -c {input_sample} | {nanofilt} {nanofilt_p} |gzip > {AF}/{sample}.filtered.fastq.gz
{nanoplot} --fastq {AF}/{sample}.filtered.fastq.gz -o {AF} {nanoplot_p}
""".format(nanoplot=nanoplot, input_sample=input_samples[index], BF=Before_filter, AF=After_filter,output_dir=output_dir,
            sample=sample, nanofilt=nanofilt, nanoplot_p=nanoplot_parameters, nanofilt_p=nanofilt_parameters)
                fw = open(sample+'.nanopore_qc.sh', 'w').write(outfile1)