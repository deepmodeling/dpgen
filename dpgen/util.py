#!/usr/bin/env python
# coding: utf-8

from dpgen import dlog

"""
some common utilities for generator, auto_test and data
"""

# constants define
MaxLength=70

def sepline(ch='-',sp='-',screen=False):
    r'''
    seperate the output by '-'
    '''
    if screen:
       print(ch.center(MaxLength,sp))
    else:
       dlog.info(ch.center(MaxLength,sp))

def box_center(ch='',fill=' ',sp="|"):
    r'''
    put the string at the center of |  |
    '''
    strs=ch.center(Len,fill)
    dlog.info(sp+strs[1:len(strs)-1:]+sp)
