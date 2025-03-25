#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import sys,os
import os.path
from os import path
import subprocess as sp
from os import listdir
from os.path import isfile, join
import pickle,gzip
import matplotlib.pyplot as plt
import gc


def sem2d_read_fault(model_name,fault_name):

    # length of the tag at the begining and end of a binary record
    # in number of single precision words (4*bytes)
    LENTAG = 2; # gfortran older versions
    LENTAG = 1;

    # assumes header file name is FltXX_sem2d.hdr
    if not os.path.isdir(model_name):
        print("Wrong path to the model directory...")
        quit()
    headfile_exist = os.path.isfile(model_name+"/"+fault_name+"_sem2d.hdr")
    initfile_exist = os.path.isfile(model_name+"/"+fault_name+"_init_sem2d.tab")
    datafile_exist = os.path.isfile(model_name+"/"+fault_name+"_sem2d.dat")
    if (not headfile_exist):
        print("Miss head file in this directory...")
        quit()
    elif (not initfile_exist):
        print("Miss init file in this directory...")
        exit()
    elif (not datafile_exist):
        print("Miss fault data files in this directory...")
        exit()

    data = {}
    f = open(model_name+"/"+fault_name+"_sem2d.hdr")
    lines = f.readlines()
    data['nx'] = int(lines[1].split()[0])
    ndat       = int(lines[1].split()[1])
    data['nt'] = int(lines[1].split()[2])
    data['dt'] = float(lines[1].split()[3])
    xyz = []
    for line in lines[4::]:
        xyz.append(line.split())
    xyz = np.asarray(xyz).astype(float)
    data['x'] = xyz[:,0]
    data['z'] = xyz[:,1]
    # Read initial fault data
    f = open(model_name+"/"+fault_name+"_init_sem2d.tab")
    lines = f.readlines()
    xyz = []
    #for line in lines[4::]:
    #    xyz.append(line.split())
    for i in range(len(lines)):
        if (i % 2) == 0:
            xyz.append(lines[i].split())
        else:
            xyz[-1].append(lines[i].split()[0])
    xyz = np.asarray(xyz).astype(float)
    data['st0'] = xyz[:,0]
    data['sn0'] = xyz[:,1]
    data['mu0'] = xyz[:,2]
    data['theta0'] = xyz[:,3]

    # Read fault data in a big matrix
    f   = open(model_name+"/"+fault_name+"_sem2d.dat", "rb")
    dt  = np.dtype((np.float32, data['nx']+2*LENTAG))
    raw = np.fromfile(f, dtype=dt)
    raw = np.reshape(raw[:,LENTAG:LENTAG+data['nx']],(int(raw.shape[0]/ndat),ndat, data['nx']));
    # Reformat each field [nx,nt]
    data['d']  = raw[:,0,:]
    data['v']  = raw[:,1,:]
    data['st'] = raw[:,2,:]
    data['sn'] = raw[:,3,:]
    data['mu'] = raw[:,4,:]
    data['theta'] = raw[:,5,:]
    if (ndat == 5+4):
        data['d1t'] = raw[:,5,:]
        data['d2t'] = raw[:,6,:]
        data['v1t'] = raw[:,7,:]
        data['v2t'] = raw[:,8,:]
    elif (ndat == 5+4*2):
        data['d1t'] = raw[:,5,:]
        data['d1n'] = raw[:,6,:]
        data['d2t'] = raw[:,7,:]
        data['d2n'] = raw[:,8,:]
        data['v1t'] = raw[:,9,:]
        data['v1n'] = raw[:,10,:]
        data['v2t'] = raw[:,11,:]
        data['v2n'] = raw[:,12,:]
    return data


def sem2d_read_STF(model_name):
    
    # length of the tag at the begining and end of a binary record
    # in number of single precision words (4*bytes)
    LENTAG = 2; # gfortran older versions
    LENTAG = 1;
    
    # assumes header file name is FltXX_sem2d.hdr
    if not os.path.isdir(model_name):
        print("Wrong path to the model directory...")
        quit()
    initfile_exist = os.path.isfile(model_name+"/SourcesTime_sem2d.tab")
    if (not initfile_exist):
        print("Miss init file in this directory...")
        exit()
    
    data = {}
    # Read initial fault data
    f = open(model_name+"/SourcesTime_sem2d.tab")
    #f = open("SourcesTime_sem2d.tab")
    lines = f.readlines()
    xyz = []
    for line in lines[4::]:
        xyz.append(line.split())
    xyz = np.asarray(xyz).astype(np.float)
    data['t'] = xyz[:,0]
    data['stf'] = xyz[:,1]
    return data