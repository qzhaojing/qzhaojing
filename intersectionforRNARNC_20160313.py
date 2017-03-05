# -*- coding: utf-8 -*-
"""
Created on Sat Mar 12 22:11:00 2016
给出两个文件RNA RNC 算出rpkm的两个file，求两者的交集,对交集的基因画散点图，存储图片
横轴纵轴名称是 文件名取20个字符。算R2#计算的是R2，不是R，

基因名            read数    rpkm
Os01t0851400	265	9.32851920899
Os03t0666500	110	5.23929433967
Os06t0109400	434	8.54697744037


@author: ZhaoJing
"""
import glob,math,os
import matplotlib
import matplotlib.pyplot as plt 
from matplotlib.pyplot import plot,savefig 
import pylab as pl
from sklearn.linear_model import LinearRegression
def union(file1,file2):
    list1={i.split()[0]:i.split()[1:] for i in open(file1).readlines()}
    list2={i.split()[0]:i.split()[1:] for i in open(file2).readlines()}
    unionset1=list(set(list1.keys())&set(list2.keys()))   #效率好高的，秒级
    finallist=[[x]+list1[x]+list2[x]+[float(list2[x][1])/float(list1[x][1])] for x in unionset1] #[name,reads1,rpkm1,reads2,rpkm2,TR]

    def linefit(x , y):
        N = float(len(x))
        sx,sy,sxx,syy,sxy=0,0,0,0,0
        for i in range(0,int(N)):
            sx  += x[i]
            sy  += y[i]
            sxx += x[i]*x[i]
            syy += y[i]*y[i]
            sxy += x[i]*y[i]
        a = (sy*sx/N -sxy)/( sx*sx/N -sxx)
        b = (sy - a*sx)/N
        r = abs(sy*sx/N-sxy)/math.sqrt((sxx-sx*sx/N)*(syy-sy*sy/N))
        return a,b,r   
        
    plotlist1=[math.log10(float(x[2])) for x in  finallist]    #画图
    plotlist2=[math.log10(float(x[4])) for x in  finallist]
    a,b,r=linefit(plotlist1 , plotlist2)
    r2=round(r**2,3)              #计算的是R2，不是R
    plottxt='R2='+str(r2)    
    pl.plot(plotlist1,plotlist2,'k.',markersize=1) 
    pl.title(os.path.basename(file1).split('.')[0]+" vs. "+os.path.basename(file2).split('.')[0])# give plot a title
    pl.xlabel(os.path.basename(file1).split('.')[0])# make axis labels
    pl.ylabel(os.path.basename(file2).split('.')[0])
    pl.xlim(-2, 6.0)# set axis limits
    pl.ylim(-2, 6.0)

    pl.text(1,5.5, plottxt ,color='blue',ha='center')#    pl.show()# show the plot on the screen
    savefig(os.path.join(filedir, os.path.basename(file1).split('.')[0]+"_vs_"+os.path.basename(file2).split('.')[0]+'.eps'),dpi=300)    
    savefig(os.path.join(filedir, os.path.basename(file1).split('.')[0]+'_vs_'+os.path.basename(file2).split('.')[0]+'.png'),dpi=300) 
    pl.close('all')
#    unionset=[x for x in list1.keys() if x in list2.keys()]   #约1min
#    finallist1=[[x]+list1[x]+list2[x] for x in list1.keys() if x in list2.keys()]  #最简练，但是速度上不是很快，不如上面两行的


#做一个文件夹里的所有文件用这个
filedir='\\\\fs\\data\\ZhaoJing\\rice_tissue201604\\total_50bp_20160906\\fanse3_IRGSP20160805\\2.Qualification\\rpkm\\un_union'   
filelist=glob.glob(filedir+'\\'+'*.normaltxt')
flag=1
for fil in filelist:
    if flag==1:
        file1=fil
        flag=2
    elif flag==2:
        file2=fil
        flag=1
        union(file1,file2)
        print os.path.basename(file1),'  and  ', os.path.basename(file2),'\n'
        
##单独做两个文件的用下面这个-----------------------------------------------------------------------------
#filedir='\\\\fs\\data\\ZhaoJing\\rice_tissue201604\\total_50bp_20160906\\fanse3_IRGSP20160805\\rpkm\\union\\zhangstyle_not' 
#file1=filedir+'\\9311_HY_RNA.fanse3_10.Reducetxt'
#file2=filedir+'\\9311_HY_RNC.fanse3_10.Reducetxt'
#union(file1,file2)
##1和7,2和9这样子，RNA和本身的RNC对应
#filedir='d:\\Zhaojing\\Translatome_mapping\\2.Mappingdata\\R1_R12DATA\\fanse3data\\fanse3rpkm'
#filelist=glob.glob(filedir+'\\'+'*.uniontxt')
#for i in range(6):
#    file1=filelist[i]
#    file2=filelist[i+6]
#    union(file1,file2)
#    print os.path.basename(file1),'  and  ', os.path.basename(file2),'\n'


