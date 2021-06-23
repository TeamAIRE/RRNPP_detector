#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RRNPP_detector: a tool to detect RRNPP-Type quorum sensing systems in gram-positive bacteria and bacteriophages genomes
@author: Bernard Charles
Team AIRE (Adaptation, Integration, Reticulation, Evolution): http://www.evol-net.fr/
"""

import os
import pandas
import shutil
import subprocess
    

def signalp(faa, out_dir, rrnpp_detector_dir, current_dir):
    try:
        signalp_dir = os.path.dirname(shutil.which('signalp'))
    except:
        print('warning: the signalp binary is not found in your $PATH. RRNPP_detector will try to run built-in signalp version 5.0b Linux x86_64')
        signalp_dir = os.path.join(rrnpp_detector_dir, 'signalp-5.0b', 'bin')
    # signalp must be run from the directory where the binary is stored
    os.chdir(signalp_dir)
    signalp_args = ['./signalp', '-org', 'gram+', '-format', 'short', '-fasta', faa, '-prefix', os.path.join(out_dir, 'propeptides')]
    subprocess.run(signalp_args, check = True)
    os.chdir(current_dir)
    df = pandas.read_csv(os.path.join(out_dir, 'propeptides_summary.signalp5'), comment='#', sep='\t', header=None)
    df.columns = ['protein_id', 'Prediction',	 'SP(Sec/SPI)', 'TAT(Tat/SPI)', 'LIPO(Sec/SPII)', 'OTHER', 'CS Position']
    return df 

