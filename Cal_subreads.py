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
                if line[1]=='raw total sequences:':
                    read_total_num = int(line[2])
                elif line[1]=='total length:':
                    read_total_length = int(line[2])
                elif line[1] == 'maximum length:':
                    read_max = int(line[2])
            elif lines.startswith('RL'):
                read_distribution[line[1]] = int(line[2])
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
                sum_tmp+=int(read_lst[i])
                lst_tmp.append(i)
        N50 = read_lst[lst_tmp.index(lst_tmp[-1])]
    return read_total_num,read_total_length,read_mean,read_max,N50,read_lst


if __name__ == '__main__':
    #samtools = getConfig('QC', 'samtools')
    #os.system('{samtools} stats {subreads_bam} > {prefix}.stat.tmp'.format(samtools=samtools,subreads_bam=bam,prefix=prefix))
    input_file = sys.argv[1]  ## file that generates by samtools stats
    read_total_num,read_total_length,read_mean,read_max,N50,read_lst = cal(input_file)
    print('total subreads number: '+str(read_total_num))
    print('total subreads length: '+str(read_total_length))
    print('average subreads length: '+str(read_mean))
    print('longest subreads length: '+str(read_max))
    print('subreads N50: '+str(N50))
    out_name=sys.argv[2]
    fw = open(out_name,'w')
    fw.write('\n'.join(read_lst))
    fw.close()
