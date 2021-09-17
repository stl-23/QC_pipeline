# QC_pipelinev1.0
The quality control pipeline for sequencing reads(DNA and RNA library) of illumina,Pacbio (RS/Sequel),Nanopore platform

# Dependencies download and usage
1. Pull our docker image from Docker Hub 
```
docker pull stl23/qc:v1.7
```
2. Run the pipeline
```
docker run -v "YOUR_INPUT_DIR":/input\
-v "YOUR_OUTPUT_DIR":/output \
-v “YOUR_WORK_DIR”:/work \
stl23/qc:v1.7 bash -c 'cd /work/ && python3 /script/qc/QC_pipeline_v1.py \
-inputs /input \
-outputs /output \
-c illumina \
-omic DNA \
-s PE \
-fastp_p "-q 20 -u 50" \
--script False'
```
注意：有时候系统的安全模块selinux会把权限禁掉，导致出现无法读取目录或文件“cannot open directory xx: Permission denied”，可以在运行时加 --privileged=true，即
```
docker run -v "YOUR_INPUT_DIR":/input\
-v "YOUR_OUTPUT_DIR":/output \
-v “YOUR_WORK_DIR”:/work \
--privileged=true \
stl23/qc:v1.7 bash -c 'cd /work/ && python3 /script/qc/QC_pipeline_v1.py \
-inputs /input \
-outputs /output \
-c illumina \
-omic DNA \
-s PE \
-fastp_p "-q 20 -u 50" \
--script False'
```
# Parameters
```
usage: QC_pipeline_v1.py [-inputs INPUTS] [-outputs OUTPUTS] [-h] [-c {illumina,pacbio,nanopore}] [-omic OMIC_TYPE] [-s SEQ_TYPE] [-fastp_p FASTP_PARAMETERS] [-mt MACHINE_TYPE]
                        [-cr {no,yes}] [-ccs_p CCS_PARAMETERS] [-lima_p LIMA_PARAMETERS] [-isoseq3_p ISOSEQ3_PARAMETERS] [-lordec_p LORDEC_CORRECT_PARAMETERS]
                        [-nanoplot_p NANOPLOT_PARAMETERS] [-nanofilt_p NANOFILT_PARAMETERS]

General options:
  -inputs INPUTS        The input directory(rawdata),the suffixes of short reads in the directory must be _1/2.fq or _1/2.fq.gz or                               _1/2.fastq or _1/2.fastq.gz;the suffixes of pacbio
                        long reads must be .subreads.bam(Sequel platform) or .1.bax.h5,.2.bax.h5,.3.bax.h5(RS/RSII platform);the                                 suffiexes of nanopore long reads must be .fastq.gz
  -outputs OUTPUTS      The output directory
  -h, --help            show the help and exit
  -c {illumina,pacbio,nanopore}, --choice {illumina,pacbio,nanopore}
                        Choose a sequencing machine,illumina for short reads,pacbio and nanopore for long reads,default is illumina
  -omic OMIC_TYPE, --omic_type OMIC_TYPE
                        The omics type:DNA(genome) or RNA(transciptome),default is DNA
  --script {True,False}
                        Write commands into shell files,do not run in local
Short reads options (illumina):
  -s SEQ_TYPE, --seq_type SEQ_TYPE
                        Paired end(PE) or singe end(SE) short reads in illumina sequencing platform,default is PE
  -fastp_p FASTP_PARAMETERS, --fastp_parameters FASTP_PARAMETERS
                        The parameters for fastp softwares,e.g: -fastp_p "-q 20 -u 50", the parameters:--html and --json are defaultly set and named by the samples names to avoid
                        overlapping in this pipline,please do not set again

Long reads options (Pacbio):
  -mt MACHINE_TYPE, --machine_type MACHINE_TYPE
                        The sequencing platform:Sequel/RS,default is Sequel
  -cr {no,yes}, --correction {no,yes}
                        Use illumina data to correct the pacbio long reads or not, default is no;if yes, please upload illumina reads with the pacbio long reads in the same directory and have a format of {same_name_as_pacbio_subreads}.xx.fq.gz (used in CLR mode,the self-correction is enough for CCS mode in normal condition)
  -ccs_p CCS_PARAMETERS, --ccs_parameters CCS_PARAMETERS
                        The parameters for SMRTlinks ccs tools,defalut parameters in DNA,"--skip-polish --minPasses 1" in RNA（smrtlink9)
  -lima_p LIMA_PARAMETERS, --lima_parameters LIMA_PARAMETERS
                        The parameters for SMRTlinks lima tools,default parameters in DNA,"--isoseq" in RNA(smrtlink9)
  -isoseq3_p ISOSEQ3_PARAMETERS, --isoseq3_parameters ISOSEQ3_PARAMETERS
                        The parameters for SMRTlinks isoseq3 tools(refine/cluster/polish),seperated by semicolon,e.g. -isoseq3_p "-j 1;--verbose;-r 0.99",default is ";--verbose;"
  -lordec_p LORDEC_CORRECT_PARAMETERS, --lordec_correct_parameters LORDEC_CORRECT_PARAMETERS
                        The parameters for lordec_correct tools,"-k 19 -s 3" in DNA,default parameters in RNA

Long reads options (Nanopore):
  -nanoplot_p NANOPLOT_PARAMETERS, --nanoplot_parameters NANOPLOT_PARAMETERS
                        The parameters for Nanoplot software
  -nanofilt_p NANOFILT_PARAMETERS, --nanofilt_parameters NANOFILT_PARAMETERS
                        The parameters for NanoFilt software

```
# Notes
- `-ccs_p` to date, the SMRT-Link sortware has different versions,which have some changes in the parameters list. We used version 7.0.0 to deal with data from RS platform (ccs parameters: --noPolish --minPasses); version 9.0 to deal with data from Sequel platform (ccs paramters: --skip-polish --min-passes).
- `-lima_p` same as above,set with different versions.

# EXAMPLES: 
```
python3 QC_pipeline_v1.py -inputs /root/my_data/example_inputs/ -outputs /root/my_data/example_outputs/ 
-c illumina -omic DNA -s PE -fastp_p "-q 20 -u 50" 

python3 QC_pipeline_v1.py -inputs /root/my_data/example_inputs/ -outputs /root/my_data/example_outputs/
-c pacbio -omic RNA -mt Sequel -ccs_p "--skip-polish --minPasses 1" -lima_p "--isoseq" -isoseq3_p ";--verbose;" -lordec_p "-m 2G" 

python3 QC_pipeline_v1.py -inputs /root/my_data/example_inputs/ -outputs /root/my_data/example_outputs/ 
-c nanopore -omic DNA -nanoplot_p "--plots hex dot " -nanofilt_p "-q 7 -l 1000 --headcrop 50 --tailcrop 50"
```
