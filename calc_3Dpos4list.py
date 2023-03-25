#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 12:25:20 2023

@author: arest
"""

import sys, os, re, copy, shutil,io,math
import argparse
if 'PIPE_PYTHONSCRIPTS' in os.environ:
    sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
    sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
from pdastro import pdastroclass
from calc_3Dpos import calc3Dposclass

class calc3Dpos4listclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        
        self.EventInfo = calc3Dposclass()
        
        self.racol='RA'
        self.deccol='Dec'
    
    def loadCoordinates(self,filename,checkcolnames=True,racol=None,deccol=None):
        if racol is not None: self.racol=racol
        if deccol is not None: self.deccol=deccol
        self.load(filename)

        # get rid of leading # in first column name, old format!
        if re.search('^\#',self.t.columns[0]) is not None:
            self.t = self.t.rename(columns={self.t.columns[0]:re.sub('^\#','',self.t.columns[0])})

        if checkcolnames:
            if not (self.racol in self.t.columns):
                raise RuntimeError(f'RA column {self.racol} does not exist in coordinate table (columns: {list(self.t.columns)})')
            if not (self.deccol in self.t.columns):
                raise RuntimeError(f'Dec column {self.deccol} does not exist in coordinate table (columns: {list(self.t.columns)})')

        # convert ra,dec into decimal degrees
        self.assert_radec_cols_decimal_degrees(racol=self.racol,deccol=self.deccol)
        
    def calc3Dpos4list(self,eventname,year,indices=None):
        ixs=self.getindices(indices)
        for ix in ixs:
            results=self.EventInfo.calc_3Dpos(eventname,self.t.loc[ix,self.racol],self.t.loc[ix,self.deccol],year,verbose=0)
            for key in ['angsep_deg','PA','alpha_deg','theta_deg','rho_ly','x_ly','y_ly','z_ly','r_ly']:
                self.t.loc[ix,key]=results[key]
        for col in ['angsep_deg','PA','alpha_deg','theta_deg','rho_ly','x_ly','y_ly','z_ly','r_ly']:
            self.default_formatters[col]='{:.2f}'.format
        for col in [self.t.loc[ix,self.racol],self.t.loc[ix,self.deccol]]:
            self.default_formatters[col]='{:.6f}'.format

    def saveresults(self,outfilename,outdir=None):
        if outfilename is None:
            raise RuntimeError('output filename is None, cannot save results!')
        if outfilename.lower() == 'auto':
            outfilename = re.sub('\.txt$','',self.filename)
            outfilename += '.xyz.txt'
        if outdir is not None:
            outfilename = os.path.abspath(f'{outdir}/{os.path.basename(outfilename)}')
        print(f'Saving results to {outfilename}')
        self.write(outfilename)
        
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coordfilename', help="Filename of the list of coordinates")
    parser.add_argument('eventname', help="Name of the event, e.g. Tycho. This name needs to be in the Event info table ")
    parser.add_argument('year', type=float, help="Year of light echo observation")
    parser.add_argument('--event_info_filename', default='EventInfo.txt', help='Event info table filename (default=%(default)s)')
    parser.add_argument('-r','--racol', default='RA', help='RA colmun name of coordinate list (default=%(default)s)')
    parser.add_argument('-d','--deccol', default='Dec', help='Dec column name of coordinate list (default=%(default)s)')
    parser.add_argument('-o','--outfilename', default=None, help='output filename. If "auto", then the output filename is the input filename, with .txt removed, and .xyz.txt added (default=%(default)s)')
    parser.add_argument('--outdir', default=None, help='if not None, then the outputfilename is stripped of the path, and saved in outdir (default=%(default)s)')
    args = parser.parse_args()
    
    
    calc3D = calc3Dpos4listclass()
    calc3D.EventInfo.loadEventInfo(filename=args.event_info_filename)
    calc3D.EventInfo.write()
    calc3D.loadCoordinates(args.coordfilename,racol=args.racol,deccol=args.deccol)
    calc3D.calc3Dpos4list(args.eventname,args.year)
    #calc3D.calc_3Dpos(args.eventname,args.RA,args.Dec,args.year)
    calc3D.write()
    if args.outfilename is not None:
        calc3D.saveresults(args.outfilename,outdir=args.outdir)