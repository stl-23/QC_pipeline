#!/usr/bin/env python
import os
import re
import argparse
import configparser
import subprocess


def getConfig(section, key):
    config = configparser.ConfigParser()
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'softwares.config')
    config.read(path)
    return config.get(section, key)


def remove_rRNA(sample, output_dir, rRNA_db, seq_type='PE', hisat_p=None, *samples):
    rRNA_out = os.path.join(output_dir, 'rm_rRNA_dir', sample)
    cmd = ''
    if seq_type == 'PE':
        input1, input2 = samples[0], samples[1]
        cmd = """mkdir -p {rRNA_out}
{hisat2} -x {rRNA_db} -1 {input1} -2 {input2} -S {rRNA_out}/{sample}.sam --un-conc-gz {rRNA_out} {hisat_p}
mv {rRNA_out}/un-conc-mate.1 {rRNA_out}/{sample}.1.rm.fq.gz 
mv {rRNA_out}/un-conc-mate.2 {rRNA_out}/{sample}.2.rm.fq.gz
rm {rRNA_out}/{sample}.sam
""".format(hisat2=hisat2,rRNA_db=rRNA_db, input1=input1, input2=input2, sample=sample,
                   rRNA_out=rRNA_out, hisat_p=hisat_p)
    elif seq_type == 'SE':
        input1 = samples[0]
        cmd = """mkdir -p {rRNA_out}
{hisat2} -x {rRNA_db} -U {input1} -S {rRNA_out}/{sample}.sam --un-gz {rRNA_out}/{sample}.rm.fq.gz {hisat_p}
rm {rRNA_out}/{sample}.sam
""".format(hisat2=hisat2,rRNA_db=rRNA_db, input1=input1, sample=sample, rRNA_out=rRNA_out, hisat_p=hisat_p)
    return cmd


def parse_short_read_dir(inputs, outs, seq_type='PE'):
    input_path = os.path.abspath(inputs)
    out_path = os.path.abspath(outs)
    lst = os.listdir(input_path)
    dic = {}
    seq_suffix = [seq + zz for seq in ['.fq', '.fastq'] for zz in ['', '.gz', '.bz2']]
    lst = [file for file in lst for s in seq_suffix if file.endswith(s)]
    if lst:
        samples = [ i.split('_1')[0] for i in lst if '_1' in i]
        suffix = [i.split('_1')[1] for i in lst if '_1' in i]
    else:
        raise IOError('No such file or directory or wrong file suffix (e.q. sample_1.fq.gz/_1.fastq.gz/_1.fq/_1.fastq)')
    if samples:
        if seq_type == 'PE':
            for index, sample in enumerate(samples):
                dic[sample] = [os.path.join(input_path,sample+'_1'+suffix[index]),
                               os.path.join(input_path,sample+'_2'+suffix[index]),
                               os.path.join(out_path,sample+'_1'+suffix[index]),
                               os.path.join(out_path,sample+'_2'+suffix[index]),
                               os.path.join(out_path,sample+'.html'),
                               os.path.join(out_path,sample+'.json')]
        elif seq_type == 'SE':
            for index, sample in enumerate(samples):
                dic[sample] = [os.path.join(input_path,sample+'_1'+suffix[index]),
                               os.path.join(out_path,sample+'_1'+suffix[index]),
                               os.path.join(out_path,sample+'.html'),
                               os.path.join(out_path,sample+'.json')]
    else:
        raise IOError('Wrong file name (e.q. sample_1.fq.gz/_1.fastq.gz/_1.fq/_1.fastq)')

    return dic
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
        raise Exception("No such file or directory or wrong file name,must be .subreads.bam or .1/2/3.bax.h5")

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
    #out_path = os.path.abspath(outs) + '/'
    lst = os.listdir(input_path)
    samples = [i.replace('.fastq.gz','') for i in lst if i.endswith('.fastq.gz')]
    samples = set(samples)
    input_samples = []
    for sample in samples:
        input_samples.append(input_path + sample + '.fastq.gz')
    return samples, input_samples

def shell_cmd(cmds,cmd_files,script):
    ## if script was set True, commands will be written into shell scripts, you can submit or run locally;
    ## else the commands will run automatically.
    for index, cmd in enumerate(cmds):
        if script == 'True':
            with open(cmd_files[index],'w') as fh:
                fh.write(cmd+'\n')
        else:
            p = subprocess.check_call(cmd, shell=True)
            if p:
                raise Exception("%s has not been done,something went wrong"%cmd_files[index])
            else:
                print("%s completed"%cmd_files[index])

 #       if len(cmds) <= bin_size:
 #           bin_num = 1
 #           bin_size = len(cmds)
 #       else:
 #           bin_num = (len(cmds) - len(cmds)%bin_size) / bin_size
 #       x = 0
 #       while x < bin_num:
 #           bin_cmd = cmds[bin_size*x:bin_size*(x+1)]

 #           x += 1
 #       if len(cmds)%bin_size > 0:
 #           last_bin_cmd = cmds[-len(cmds)%bin_size:]
 #           for cmd in last_bin_cmd:
 #               subprocess.check_call(cmd, shell=True)


if __name__ == '__main__':
    #### Parse arguments ####
    examplelog="""EXAMPLES:
    python QC_pipline_v1.py -inputs /root/my_data/example_inputs/ -outputs /root/my_data/example_outputs/ -c illumina -omic DNA -s PE -fastp_p "-w 3" \
    python QC_pipline_v1.py -inputs /root/my_data/example_inputs/ -outputs /root/my_data/example_outputs/ -c pacbio -omic RNA -mt Sequel -ccs_p "--skip-polish --minPasses 1" -lima_p "--isoseq --no-pbi" -isoseq3_p ";--verbose;" -lordec_p "-m 2G" \
    python QC_pipline_v1.py -inputs /root/my_data/example_inputs/ -outputs /root/my_data/example_outputs/ -c nanopore -omic DNA -nanoplot_p "--plots hex dot pauvre kde"
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
    general.add_argument('--script',type=str,default='True',choices=['True','False'],
                         help="Write commands into shell files,do not run in local")

    illumina = parser.add_argument_group(title='Short reads options (illumina)')
    illumina.add_argument('-s', '--seq_type', type=str, default='PE',
                        help='Paired end(PE) or singe end(SE) short reads in illumina sequencing platform,default is PE')
    illumina.add_argument('-fastp_p', '--fastp_parameters', type=str,default='-q 20 -u 50',
                        help='The parameters for fastp softwares,e.g: -fastp_p "-g -x -5 -3", the parameters:--html and --json are defaultly set and named by the samples names to avoid overlapping in this pipline,please do not set again')
    pacbio = parser.add_argument_group(title='Long reads options (Pacbio)')
    pacbio.add_argument('-mt', '--machine_type', type=str,default='Sequel',
                        help='The sequencing platform:Sequel/RS,default is Sequel')
    pacbio.add_argument('-cr', '--correction', type=str, default='no', choices=['no', 'yes'],
                        help='Use illumina data to correct the pacbio long reads or not, default is no;if yes, please upload illumina reads with the pacbio long reads in the same directory and have a format of {same_name_as_pacbio_subreads}.xx.fq.gz (used in CLR mode,the self-correction is enough for CCS mode in normal condition)')
    pacbio.add_argument('-ccs_p', '--ccs_parameters', type=str, default='',
                        help='The parameters for SMRTlinks ccs tools,defalut parameters in DNA,"--skip-polish --minPasses 1" in RNA')
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
    script = args.script


    ##### Parse software and database config ####
    fastp = getConfig('QC','fastp').strip("'")
    smrtlink9 = getConfig('QC','smrtlink9').strip("'")
    smrtlink7 = getConfig('QC','smrtlink7').strip("'")
    ccs9 = os.path.join(smrtlink9,'ccs')
    ccs7 = os.path.join(smrtlink7,'ccs')
    lima = os.path.join(smrtlink9,'lima')
    bax2bam = os.path.join(smrtlink7,'bax2bam')
    bam2fasta = os.path.join(smrtlink9,'bam2fasta')
    isoseq3 = os.path.join(smrtlink9,'isoseq3')
    samtools = os.path.join(smrtlink9,'samtools')
    pbindex = os.path.join(smrtlink9,'pbindex')

    lordec_correct = getConfig('QC','lordec_correct').strip("'")
    lordec_trim = getConfig('QC','lordec_trim').strip("'")
    nanoplot = getConfig('QC','nanoplot').strip("'")
    nanofilt = getConfig('QC','nanofilt').strip("'")
    hisat2 = getConfig('QC','hisat2').strip("'")

    rRNA_data = os.path.join(os.path.abspath(os.path.dirname(os.path.dirname(__file__))),'database','rRNA','rRNAs')
    primers = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'primers.fasta')
    seq_suffix = ['.ct.fa','.ct.fq','.ct.fa.gz','.ct.fq.gz','.ct.fasta','.ct.fasta.gz','.ct.fastq','.ct.fastq.gz']
    correct_data_lst = [os.path.abspath(input_dir)+'/'+ file for file in os.listdir(os.path.abspath(input_dir))
                        for suffix in seq_suffix if file.endswith(suffix)]
    cal_subreads = os.path.join(os.path.abspath(os.path.dirname(__file__)),'Cal_subreads.py')
    read_distribution_hist = os.path.join(os.path.abspath(os.path.dirname(__file__)),'read_distribution_hist.R')

    #### Main ####
    cmds = []
    cmd_files = []
    if platform_choice == 'illumina':
        dic = parse_short_read_dir(input_dir, output_dir, seq_type)
        if omic_type == 'DNA':
            if seq_type == 'PE':
                for sample, path in dic.items():
                    cmds.append('{} -i {} -I {} -o {} -O {} --html {} --json {} {}'.format(fastp, path[0], path[1],
                                                                                   path[2], path[3],
                                                                                   path[4], path[5],
                                                                                   fastp_parameters))
                    cmd_files.append(sample + '.sh')
                shell_cmd(cmds, cmd_files, script)

            elif seq_type == 'SE':
                for sample, path in dic.items():
                    cmds.append('{} -i {} -o {} --html {} --json {} {}'.format(fastp, path[0], path[1],
                                                                       path[2], path[3], fastp_parameters))
                    cmd_files.append(sample + '.sh')
                shell_cmd(cmds, cmd_files, script)

        elif omic_type == 'RNA':  ## remove rRNA first
            rm_rRNA_seq1 = {}
            rm_rRNA_seq2 = {}
            rm_RNA_cmds = []
            run_cmds = []
            if seq_type == 'PE':
                for sample, path in dic.items():
                    rm_rRNA_seq1[sample] = os.path.join(output_dir, 'rm_rRNA_dir', sample, sample+'.1.rm.fq.gz')
                    rm_rRNA_seq2[sample] = os.path.join(output_dir, 'rm_rRNA_dir', sample, sample+'.2.rm.fq.gz')
                    rm_RNA_cmds.append(remove_rRNA(sample, output_dir, rRNA_data, seq_type, '', path[0], path[1]))
                    run_cmds.append('{} -i {} -I {} -o {} -O {} --html {} --json {} {}'.format(fastp, rm_rRNA_seq1[sample],
                                                                                           rm_rRNA_seq2[sample],
                                                                                           path[2],path[3],
                                                                                           path[4], path[5],
                                                                                           fastp_parameters))

                    cmd_files.append(sample + '.sh')
                cmds = map(lambda x,y:x+'\n'+y,rm_RNA_cmds,run_cmds)
                shell_cmd(cmds, cmd_files, script)

            elif seq_type == 'SE':
                for sample, path in dic.items():

                    rm_rRNA_seq1[sample] = os.path.join(output_dir, 'rm_rRNA_dir',sample, sample+'.rm.fq.gz')
                    rm_RNA_cmds.append(remove_rRNA(sample, output_dir, rRNA_data, seq_type, '', path[0]))
                    run_cmds.append('{} -i {} -o {} --html {} --json {} {}'.format(fastp, rm_rRNA_seq1[sample], path[1],
                                                                               path[2], path[3],
                                                                               fastp_parameters))
                    cmd_files.append(sample + '.sh')
                cmds = map(lambda x,y:x+'\n'+y,rm_RNA_cmds,run_cmds)
                shell_cmd(cmds, cmd_files, script)
    elif platform_choice == 'pacbio':
        index_cmd = []
        static_cmd = []
        correct_cmd = []
        samples, input1, input2, input3, output1 = parse_pacbio_read_dir(input_dir, output_dir, mt_type)

        if omic_type == 'DNA':
            if lordec_correct_parameters == '':
                lordec_correct_parameters = "-k 19 -s 3"
            else:
                lordec_correct_parameters = lordec_correct_parameters

            if mt_type == 'Sequel':
                for index, sample in enumerate(samples):
                    if not re.search('--report-file', ccs_parameters):
                        ccs_parameters += ' --report-file %s.ccs.report.txt' % output1[index]
                    else:
                        ccs_parameters = ccs_parameters
                    static_cmd.append("""#!/usr/bin/bash
{{
{samtools} stats {input1}.subreads.bam > {output1}.subreads.stat.tmp
python {Cal_subreads} {output1}.subreads.stat.tmp {output1}.subreads_length.list > {output1}.subreads.stat.list
Rscript {read_distribution_hist} {output1}.subreads_length.list {output1}.subreads
}}&
{{
{ccs9} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}
{samtools} stats {output1}.ccs.bam > {output1}.ccs.stat.tmp
python {Cal_subreads} {output1}.ccs.stat.tmp {output1}.ccs_length.list > {output1}.ccs.stat.list
Rscript {read_distribution_hist} {output1}.ccs_length.list {output1}.ccs
{bam2fasta} {output1}.ccs.bam -o {output1}.pacbio
}}&
rm *.stat.tmp
""".format(samtools=samtools,input1=input1[index],Cal_subreads=cal_subreads,output1=output1[index],
            read_distribution_hist=read_distribution_hist,ccs9=ccs9,ccs_p=ccs_parameters,
            bam2fasta=bam2fasta))

                    if correct == 'yes':
                        try:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input1[index])]
                            correct_data = ','.join(correct_data_lst_tmp)
                            correct_cmd.append("""{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta.gz -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim,correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters))
                        except IOError as e:
                            print(e)
                    cmd_files.append(sample + '.sh')
                if correct_cmd:
                    cmds = map(lambda x,y:x+'\n'+y,static_cmd,correct_cmd)
                else:
                    cmds = static_cmd
                shell_cmd(cmds, cmd_files,script)

            elif mt_type == 'RS':
                for index, sample in enumerate(samples):
                    input_1, input_2, input_3 = input_dir + sample + '.1.bax.h5', input_dir + sample + '.2.bax.h5', \
                                                input_dir + sample + '.3.bax.h5'
                    if not re.search('--reportFile', ccs_parameters):
                        ccs_parameters += ' --reportFile %s.ccs.report.txt' % output1[index]
                    else:
                        ccs_parameters = ccs_parameters
                    static_cmd.append("""#!/usr/bin/bash
{bax2bam} {input_1} {input_2} {input_3} -o {output1}
{{
{samtools} stats {output1}.subreads.bam > {output1}.subreads.stat.tmp
python {Cal_subreads} {output1}.subreads.stat.tmp {output1}.subreads_length.list > {output1}.subreads.stat.list
Rscript {read_distribution_hist} {output1}.subreads_length.list {output1}.subreads
}}&
{{
{ccs7} {output1}.subreads.bam {output1}.ccs.bam {ccs_p}
{samtools} stats {output1}.ccs.bam > {output1}.ccs.stat.tmp
python {Cal_subreads} {output1}.ccs.stat.tmp {output1}.ccs_length.list > {output1}.ccs.stat.list
Rscript {read_distribution_hist} {output1}.ccs_length.list {output1}.ccs
{bam2fasta} {output1}.ccs.bam -o {output1}.pacbio
}}&
rm *.stat.tmp
""".format(bax2bam=bax2bam, input_1=input_1, input_2=input_2, input_3=input_3, sample=sample, output_dir=output_dir,
            samtools = samtools, input1 = input1[index], Cal_subreads = cal_subreads, output1 = output1[index],
            read_distribution_hist = read_distribution_hist, ccs7 = ccs7, ccs_p = ccs_parameters,
            bam2fasta=bam2fasta))

                    if correct == 'yes':
                        try:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input_dir+sample)]
                            correct_data = ','.join(correct_data_lst_tmp)
                            correct_cmd.append("""{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta.gz -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim, correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters))
                        except IOError as e:
                            print(e)
                    cmd_files.append(sample + '.sh')
                if correct_cmd:
                    cmds = map(lambda x,y:x+'\n'+y,static_cmd,correct_cmd)
                else:
                    cmds = static_cmd
                shell_cmd(cmds, cmd_files,script)

        elif omic_type == 'RNA':
            isoseq3_refine = ''
            isoseq3_cluster = ''
            isoseq3_polish = ''

            if lima_parameters == '':
                #lima_parameters = "--isoseq --no-pbi" ## smrt-link v7.0
                lima_parameters = "--isoseq"  ## smrt-link v9.0
            else:
                lima_parameters = lima_parameters
            if isoseq3_parameters == '':
                isoseq3_parameters = ";--verbose;"
            else:
                isoseq3_parameters = isoseq3_parameters
            isoseq3_refine, isoseq3_cluster, isoseq3_polish = isoseq3_parameters.strip().split(';')

            if mt_type == 'Sequel':
                if not re.search('--skip-polish', ccs_parameters):
                    ccs_parameters += " --skip-polish"  ## smrt-link v9.0
                if not re.search('--min-passes', ccs_parameters):
                    ccs_parameters += " --min-passes 1"  ## smrt-link v9.0

                for index, sample in enumerate(samples):
                    if not re.search('--report-file', ccs_parameters):
                        ccs_parameters += ' --report-file %s.ccs.report.txt' % output1[index]
                    else:
                        ccs_parameters = ccs_parameters

                    static_cmd.append("""#!/usr/bin/bash
{{
{samtools} stats {input1}.subreads.bam > {output1}.subreads.stat.tmp
python {Cal_subreads} {output1}.subreads.stat.tmp {output1}.subreads_length.list > {output1}.subreads.stat.list
Rscript {read_distribution_hist} {output1}.subreads_length.list {output1}.subreads
}}&
{{
{ccs9} {input1}.subreads.bam {output1}.ccs.bam {ccs_p}
{samtools} stats {output1}.ccs.bam > {output1}.ccs.stat.tmp
python {Cal_subreads} {output1}.ccs.stat.tmp {output1}.ccs_length.list > {output1}.ccs.stat.list
Rscript {read_distribution_hist} {output1}.ccs_length.list {output1}.ccs

{lima} {output1}.ccs.bam {primers} {output1}.demux.ccs.bam {lima_p}
{isoseq3} refine {output1}.demux.ccs.primer_5p--primer_3p.bam {primers} {output1}.flnc.bam {isoseq3_refine}
{isoseq3} cluster {output1}.flnc.bam {output1}.unpolished.bam {isoseq3_cluster}
{isoseq3} polish {output1}.unpolished.bam {input1}.subreads.bam {output1}.polished.bam {isoseq3_polish}
{bam2fasta} {output1}.polished.bam -o {output1}.pacbio
}}&
rm *.stat.tmp
""".format(samtools=samtools,input1=input1[index],Cal_subreads=cal_subreads,output1=output1[index],
            read_distribution_hist=read_distribution_hist,ccs9=ccs9, ccs_p=ccs_parameters, lima=lima,
            lima_p=lima_parameters, isoseq3=isoseq3, isoseq3_refine=isoseq3_refine,
            isoseq3_cluster=isoseq3_cluster,isoseq3_polish=isoseq3_polish, bam2fasta=bam2fasta,primers=primers))

                    if not os.path.exists('%s.subreads.bam.pbi' %input1[index]):
                        index_cmd.append(f"""{pbindex} {input1}.subreads.bam""".format(pbindex=pbindex,
                                                                                       input1=input1[index]))
                    if correct == 'yes':
                        try:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input1[index])]
                            correct_data = ','.join(correct_data_lst_tmp)
                            correct_cmd.append("""{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta.gz -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim,correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters))
                        except IOError as e:
                            print(e)
                    cmd_files.append(sample + '.sh')
                if correct_cmd:
                    cmds = map(lambda x,y,z:x+'\n'+y+'\n'+z,index_cmd,static_cmd,correct_cmd)
                else:
                    cmds = map(lambda x,y:x+'\n'+y,index_cmd,static_cmd)
                shell_cmd(cmds, cmd_files, script)

            elif mt_type == 'RS':
                if not re.search('--noPolish', ccs_parameters):
                    ccs_parameters += " --noPolish"  ## smrt-link v7.0
                if not re.search('--minPasses', ccs_parameters):
                    ccs_parameters += " --minPasses 1"  ## smrt-link v7.0
                for index, sample in enumerate(samples):
                    input_1, input_2, input_3 = input_dir + sample + '.1.bax.h5', input_dir + sample + '.2.bax.h5', \
                                                input_dir + sample + '.3.bax.h5'
                    if not re.search('--reportFile', ccs_parameters):
                        ccs_parameters += ' --reportFile %s.ccs.report.txt' % output1[index]
                    else:
                        ccs_parameters = ccs_parameters
                    static_cmd.append("""#!/usr/bin/bash
{bax2bam} {input_1} {input_2} {input_3} -o {output1}
{{
{samtools} stats {output1}.subreads.bam > {output1}.subreads.stat.tmp
{pbindex} {output1}.subreads.bam
python {Cal_subreads} {output1}.subreads.stat.tmp {output1}.subreads_length.list > {output1}.subreads.stat.list
Rscript {read_distribution_hist} {output1}.subreads_length.list {output1}.subreads
}}&
{{
{ccs7} {output1}.subreads.bam {output1}.ccs.bam {ccs_p}
{samtools} stats {output1}.ccs.bam > {output1}.ccs.stat.tmp
python {Cal_subreads} {output1}.ccs.stat.tmp {output1}.ccs_length.list > {output1}.ccs.stat.list
Rscript {read_distribution_hist} {output1}.ccs_length.list {output1}.ccs

{lima} {output1}.ccs.bam {primers} {output1}.demux.ccs.bam {lima_p}
{isoseq3} refine {output1}.demux.ccs.primer_5p--primer_3p.bam {primers} {output1}.flnc.bam {isoseq3_refine}
{isoseq3} cluster {output1}.flnc.bam {output1}.unpolished.bam {isoseq3_cluster}
{isoseq3} polish {output1}.unpolished.bam {output1}.subreads.bam  {output1}.polished.bam {isoseq3_polish}

{bam2fasta} {output1}.polished.bam -o {output1}.pacbio
}}&
rm *.stat.tmp
""".format(bax2bam=bax2bam, input_1=input_1, input_2=input_2, input_3=input_3, sample=sample, output_dir=output_dir,
            samtools = samtools, input1 = input1[index], Cal_subreads = cal_subreads, output1 = output1[index],
            read_distribution_hist = read_distribution_hist, ccs7 = ccs7, ccs_p = ccs_parameters, lima=lima,
            lima_p=lima_parameters, isoseq3=isoseq3, isoseq3_refine=isoseq3_refine,pbindex=pbindex,
            isoseq3_cluster=isoseq3_cluster, isoseq3_polish=isoseq3_polish, bam2fasta=bam2fasta,primers=primers))
                    if correct == 'yes':
                        try:
                            correct_data_lst_tmp = [file for file in correct_data_lst if file.startswith(input_dir+sample)]
                            correct_data = ','.join(correct_data_lst_tmp)
                            correct_cmd.append("""{lordec_correct} -2 {correct_data} -i {output1}.pacbio.fasta.gz -o {output1}.corrected.fasta {lordec_p}
{lordec_trim} -i {output1}.corrected.fasta -o {output1}.corrected.trimed.fasta
""".format(lordec_correct=lordec_correct, lordec_trim=lordec_trim,correct_data=correct_data,
            output1=output1[index], lordec_p=lordec_correct_parameters))
                        except IOError as e:
                            print(e)
                    cmd_files.append(sample + '.sh')
                if correct_cmd:
                    cmds = map(lambda x,y:x+'\n'+y,static_cmd,correct_cmd)
                else:
                    cmds = static_cmd
                shell_cmd(cmds, cmd_files, script)

    elif platform_choice == 'nanopore':
        samples, input_samples = parse_nanopore_read_dir(input_dir, output_dir)
        Before_filter = os.path.join(output_dir,'Before_filter')
        After_filter = os.path.join(output_dir,'After_filter')
        if omic_type == 'DNA' or omic_type == 'RNA':
            for index, sample in enumerate(samples):
                outfile1 = """mkdir -p {BF}/{sample} {AF}/{sample}
{nanoplot} --fastq {input_sample} -o {BF}/{sample} {nanoplot_p}
gunzip -c {input_sample} | {nanofilt} {nanofilt_p} |gzip > {AF}/{sample}/{sample}.filtered.fastq.gz
{nanoplot} --fastq {AF}/{sample}/{sample}.filtered.fastq.gz -o {AF}/{sample} {nanoplot_p}
""".format(nanoplot=nanoplot, input_sample=input_samples[index], BF=Before_filter, AF=After_filter,output_dir=output_dir,
            sample=sample, nanofilt=nanofilt, nanoplot_p=nanoplot_parameters, nanofilt_p=nanofilt_parameters)
                with open(sample + '.sh', 'w') as fh:
                    fh.write(outfile1)

