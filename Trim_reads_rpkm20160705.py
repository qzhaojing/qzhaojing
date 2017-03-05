# -*- coding: utf-8 -*-
"""
Spyder Editor
用于对前面生成的rpkm结果卡值
可以针对rpkm，reads数量的卡值
This is a temporary script file.


"""

def trim_rpkm(fileFN,rpkm=0,reads=0):
    with open(fileFN) as f:
        alist=f.readlines()
    result_list=[]
    for line in alist:
        # line=alist[0]
        line=line.split()
        line[1],line[2]=int(line[1]),float(line[2])
        if reads==0 and line[2]>=rpkm:
            line=[str(i) for i in line]
            result_list.append('\t'.join(line)+'\n')
        elif rpkm==0 and line[1]>=reads:
            line=[str(i) for i in line]
            result_list.append('\t'.join(line)+'\n')
        elif rpkm!=0 and reads!=0 and line[1]>=reads and line[2]>=rpkm:
            line=[str(i) for i in line]
            result_list.append('\t'.join(line)+'\n')
    result_name=os.path.splitext(fileFN)
    result_name=result_name[0]+'_rpkm'+str(rpkm)+'_reads'+str(reads)+result_name[1]
    with open(result_name,'w') as f:
        f.writelines(result_list)
        
        
        
        
import os,glob
rpkm=1
reads=10
filedir='\\\\fs\\DATA\\ZhaoJing\\rice_tissue201604\\total_50bp_20160906\\fanse3_IRGSP20160805\\2.Qualification\\rpkm\\un_union'    
filelist=glob.glob(filedir+'\\'+'*.normaltxt')

for fileFN in filelist:
    trim_rpkm(fileFN,rpkm,reads)