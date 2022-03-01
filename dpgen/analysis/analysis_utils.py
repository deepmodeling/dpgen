#!/usr/bin/env python3

import os, sys

###########################################################
def select_logs(root_path, fname='OUTCAR'):
    """
    Find all fname file in root_path and its all sub-dirs
    This code is written by dfz

    Adopted at 20211217
    """
    logs = []
    for r, d, fs in os.walk(root_path):
        for f in fs:
            if os.path.islink(f):
                continue
            if f == fname:
                logs.append(os.path.join(r, f))
    return sorted(logs)
###########################################################


###########################################################
def process_bar(I, Max, Process=' Reading File:'):
    percent = I * 100.0 / Max
    e = int(I * 20 / Max)
    y = 20 - e
    process_bar = '{:s}'.format(Process) + '[' + '>' * e + '-' * y + ']' + '{:.2f}'.format(percent) + '%' +  '\r'
    sys.stdout.write(process_bar)
    sys.stdout.flush()
###########################################################
