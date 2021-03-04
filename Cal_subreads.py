#!/usr/bin/env python
import os
import sys

def cal(file):
    read_total_num = 0
    read_total_length = 0
    read_distribution = {}
    read_mean = 0
    read_max = 0
    N50 = ''
    read_lst=[]
    with open(file) as fh:
        for lines in fh:
            line = lines.strip().split('\t')
            if lines.startswith('SN'):
                if line[0]=='raw total sequences:':
                    read_total_num = int(line[1])
                elif line[0]=='total length:':
                    read_total_length = int(line[1])
                elif line[0] == 'maximum length:':
                    read_max = int(line[1])
            elif lines.startswith('RL'):
                read_distribution[line[2]] = int(line[3])
        read_mean = round(float(read_total_length)/read_total_num,2)

        for k,v in read_distribution.items():
            if v == 1:
                read_lst.append(k)
            elif v > 1:
                while v > 0:
                    read_lst.append(k)
                    v-=1
            read_lst = sorted(read_lst, reverse=True)
        sum_tmp=0
        lst_tmp=[]
        for i in range(len(read_lst)):
            if sum_tmp < read_total_length/2.0:
                sum_tmp+=read_lst[i]
                lst_tmp.append(i)
        N50 = read_lst[lst_tmp.index(lst_tmp[-1])]
    return read_total_num,read_total_length,read_mean,read_max,N50,read_lst


if __name__ == '__main__':
    #samtools = getConfig('QC', 'samtools')
    #os.system('{samtools} stats {subreads_bam} > {prefix}.stat.tmp'.format(samtools=samtools,subreads_bam=bam,prefix=prefix))
    input_file = sys.argv[1]  ## file that generates by samtools stats
    read_total_num,read_total_length,read_mean,read_max,N50,read_lst = cal(input_file)
    print('total subreads number: ':str(read_total_num))
    print('total_subreads length: ':str(read_total_length))
    print('average subreads length: ': str(read_mean))
    print('longest subreads length: ': str(read_max))
    print('subreads N50: ':str(N50))
    out_prefix=sys.argv[2]
    outfile_name = out_prefix+'subread_length.list'
    fw = open(outfile_name,'w')
    fw.write('\n'.join(read_lst))
    fw.close()
