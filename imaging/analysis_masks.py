# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 09:59:47 2020

@author: johns

This program will generate figures that compare the segmentation of each brain by 
slice in the z-direction. It will plot the slices of the segmentation as images.

"""

import os

import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np

def graph_inator(binmask,size):
    ''' 
    This function seperates a binary mask numpy array into a
    vector for x and y. Pass in size as a tuple.
    '''
    xvec=[]
    yvec=[]
    for i in range(size[0]):
        for j in range(size[1]):
            if binmask[i,j]==1:
                xvec.append(i)
                yvec.append(j)

    return xvec,yvec


def switch_inator(array):
    '''Pass in a 3-dimensional array, and the 2nd and 3rd dimensions will
    be swapped and the old 3rd dimension flipped. '''
    
    out1=np.transpose(array,(0,2,1))
    out=np.flip(out1,axis=1)
   
    return out


def brainplotter(maskdatlst,slc,
                 clst=['crimson','orange','gold','green','blue','purple'],
                 switchlst=[False,False,False,False,False,False]):
    '''
    This function takes in a list of 6 binary masks and graphs a particular slice 'slc'. This is built for 
    the Holland mouse brain research project. The switchlst keyword allows users to flip the axes 
    in a way useful for the project.
    '''
    
    gdat=[]
    for x in range(6):
        if switchlst[x]:
            gdat.append(switch_inator(maskdatlst[x]))
        else:
            gdat.append(maskdatlst[x])
    
    plt.figure()
    
    for icnt in range(6):
        axi=plt.subplot(2, 3, icnt+1)
        (xv,yv)=graph_inator(gdat[icnt][:,:,slc],(180,90))
        axi.plot(xv,yv,'.',color=clst[icnt])
        axi.axis([0,200,0,150])

    return

if __name__ == '__main__':

    brain_lst=['F01','F02','F03','M01','M02','M03']
    avgs=[]

    #l_add is the length of the final addition vector, and offs is the slice distance from the center that 
    # should be used to compute the average. 

    l_add=90
    offs=20

    datlist04=[]
    datlist08=[]
    datlist12=[]
    datlist16=[]
    datlist20=[]

    for name in brain_lst:
        xdata04=nib.load('all masks/Wk04'+name+'_mask.nii.gz').get_fdata()
        xdata08=nib.load('all masks/Wk08'+name+'_mask.nii.gz').get_fdata()
        xdata12=nib.load('all masks/Wk12'+name+'_mask.nii.gz').get_fdata()
        xdata16=nib.load('all masks/Wk16'+name+'_mask.nii.gz').get_fdata()    
        xdata20=nib.load('all masks/Wk20'+name+'_mask.nii.gz').get_fdata()
        
        datlist04.append(xdata04)
        datlist08.append(xdata08)
        datlist12.append(xdata12)
        datlist16.append(xdata16)
        datlist20.append(xdata20)
        
        if name!='F01':
            addlist=sum(sum(switch_inator(xdata20)))
        else:
            addlist=sum(sum(xdata20))
        
        i1=int(l_add/2-offs)
        i2=int(l_add/2+offs)
        
        xavg=sum(addlist[i1:i2])/(i2-i1)
        
        avgs.append(xavg)

    myslc=70

    brainplotter(datlist04,myslc)
                 #clst=['firebrick','maroon','darkred','hotpink','deeppink','lightpink'])
    plt.savefig('mask comparison/Week04Compare.jpg')

    brainplotter(datlist08,myslc)
                 #clst=['skyblue','powderblue','cadetblue','lime','lawngreen','mediumseagreen'])
    plt.savefig('mask comparison/Week08Compare.jpg')

    brainplotter(datlist12,myslc)
                 #clst=['dodgerblue','royalblue','navy','forestgreen','seagreen','darkgreen'])
    plt.savefig('mask comparison/Week12Compare.jpg')

    brainplotter(datlist16,myslc,
                 #clst=['chocolate','sandybrown','sienna','teal','darkturquoise','lightseagreen'],
                switchlst=[True,False,False,False,False,False])
    plt.savefig('mask comparison/Week16Compare.jpg')

    brainplotter(datlist20,myslc,
                 #clst=['mediumslateblue','darkviolet','rebeccapurple','crimson','coral','orangered'],
                 switchlst=[False,True,True,True,True,True])
    plt.savefig('mask comparison/Week20Compare.jpg')


