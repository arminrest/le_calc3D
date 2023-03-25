#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 19:00:54 2023

@author: arest
"""

import sys, os, re, copy, shutil,io,math
import argparse
if 'PIPE_PYTHONSCRIPTS' in os.environ:
    sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/LESPEC')
    sys.path.append(os.environ['PIPE_PYTHONSCRIPTS']+'/tools')
from pdastro import pdastroclass
#from tools import skydist_degree,calcPA
from lightechoprocs import z_ly_atD,rho_ly_zt,pc2ly,deg2rad,rad2deg

### Converts sexigesimal notation to decimal degrees or decimal hours (if option 'ra=True' invoked)
def sex2deg(sexigecimal, ra=False):
    #print(sexigecimal)
    ### 2005/12/02 - AR: make sure it is in sexagesimal format!
    # is it a string? if not check if it is None
    if not (type(sexigecimal) is str): #types.StringType):
        if type(sexigecimal) == None:
            raise RuntimeError("ERROR: sex2deg cannot handle 'None'!")
        return sexigecimal
    # Does it have ':' or ' '? If not, it must be a float in string format, just convert it and return
    if re.search('\:|\s',sexigecimal) == None:
        return(float(sexigecimal))

    s1, s2, s3 = list(map(float, re.split('[DHMSdhms:\s]', sexigecimal.strip())[0:3]))
    # Get the sign
    if re.search('-', sexigecimal):
        sign = -1
    else:
        sign = 1

    deg = abs(s1) + s2 / 60. + s3 / 3600.

    deg *= sign

    if ra:
        deg *= 15.

    return deg

### Converts decimal degrees or hours (ra=True) to sexigesimal notation
###  [-+]DD:MM:SS.ss
### the format is fixed at two decimals of precision for the decimal seconds.
### No -/+ if ra=True
def deg2sex(degrees, ra=False, outputformatRA='%02d:%02d:%06.3f',outputformatDEC='%1s%02d:%02d:%05.2f'):
    if type(degrees) is str: # types.StringType:
        degrees=float(degrees)
    if ra:
        # a.k.a. indeg and outhours
        if degrees < 0.0:
            while degrees<0:degrees+=360.
        if degrees > 360.0:
            while degrees>360.0:degrees-=360.
        degrees /= 15.

    if degrees < 0:
        sign = '-'
    else:
        sign = '+'

    degrees = abs(degrees)

    d1  = (degrees - (degrees % 1))
    rem = (degrees % 1) * 60
    d2  = rem - (rem % 1)
    srem = (rem % 1) * 60
    d3 = srem

    if ra:
      return outputformatRA % (d1, d2, d3)
    else:
      return outputformatDEC % (sign, d1, d2, d3)

# Returns the passed in RA in decimal degrees
# input RA can be in 'HH:MM:SS.ss', 'HH MM SS.ss' or in decimal degrees
def RaInDeg(Ra):
    if type(Ra)==str: #types.StringType:
        if re.search('[DHMShms: ]',Ra.strip()):
            return(sex2deg(Ra,ra=True))
        return(float(Ra))
    else:
        return(Ra)

# Returns the passed in Dec in decimal degrees
# input Dec can be in 'DD:MM:SS.ss', 'DD MM SS.ss' or in decimal degrees
def DecInDeg(Dec):
    if type(Dec)==str:#types.StringType:
        if re.search('[DHMShms: ]',Dec.strip()):
            return(sex2deg(Dec,ra=False))
        return(float(Dec))
    else:
        return(Dec)

# true angular distance between two position in the sky
#http://www2.sjsu.edu/faculty/watkins/sphere.htm
#cos(A) = sin(phi2)sin(phi1)+ cos(phi2)cos(phi1)cos(theta2 - theta1)
#longitude=theta=RA
#latitude=phi=DEC
#cos(A) = sin(DEC2)sin(DEC1)+ cos(DEC2)cos(DEC1)cos(RA2 - RA1)
#works also at RA=0 etc:
#skydist_degree(ra1,dec1,ra2,dec2) = skydist_degree(ra1+360,dec1+360,ra2-360,dec2-360) etc!
def skydist_degree(ra1,dec1,ra2,dec2):
    deg2rad=0.0174532925199  # math.pi / 180.
    #convert all values to radians
    ra1 =RaInDeg(ra1)*deg2rad
    ra2 =RaInDeg(ra2)*deg2rad
    dec1=DecInDeg(dec1)*deg2rad
    dec2=DecInDeg(dec2)*deg2rad
    if ((ra1==ra2) and (dec1==dec2)): return 0
    return math.acos(math.sin(dec2)*math.sin(dec1)+ math.cos(dec2)*math.cos(dec1)*math.cos(ra2 - ra1))*180.0/math.pi

def calcPA(ra1,dec1,ra2,dec2):
    deg2rad=0.0174532925199      # math.pi / 180.
    rad2deg=57.2957795131      # 180.0 / math.pi
    ra1=RaInDeg(ra1)
    ra2=RaInDeg(ra2)
    dec1=DecInDeg(dec1)
    dec2=DecInDeg(dec2)
    dec1rad=dec1*deg2rad
    dec2rad=dec2*deg2rad

    if ra1-ra2>180:
        ra2=ra2+360.0
    if ra1-ra2<-180:
        ra2=ra2-360.0

    dRa = ra2*deg2rad - ra1*deg2rad

    pa_deg = rad2deg*math.atan2(math.sin(dRa),math.cos(dec1rad)*math.sin(dec2rad)/math.cos(dec2rad)-math.sin(dec1rad)*math.cos(dRa));
    return pa_deg



class calc3Dposclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        self.info_colnames = {}
        self.info_colnames['RA']='RA'
        self.info_colnames['Dec']='Dec'
        self.info_colnames['Distance_pc']='Distance_pc'
        self.info_colnames['Eventname']='Name'
        self.info_colnames['yr_expl']='year_explosion'
        
    def loadEventInfo(self,filename=None,checkcolnames=True):
        if filename is None:
            filename = 'EventInfo.txt'
        self.load(filename)
        if checkcolnames:
            for col in self.info_colnames:
                if not (self.info_colnames[col] in self.t.columns):
                    raise RuntimeError(f'Column {self.info_colnames[col]} does not exist in Info table (columns: {list(self.t.columns)})')
        self.t[self.info_colnames['yr_expl']]=self.t[self.info_colnames['yr_expl']].astype(float)
        self.t[self.info_colnames['Distance_pc']]=self.t[self.info_colnames['Distance_pc']].astype(float)

    def calc_3Dpos(self,eventname,RA,Dec,Year,verbose=1):
        # get the index for the row in the Info Table for event
        ixs = self.ix_equal(self.info_colnames['Eventname'],eventname)
        if len(ixs)==0:
            raise RuntimeError(f'Could not find event "{eventname}" in column {self.info_colnames["Eventname"]} of info table {self.filename}')
        elif len(ixs)==1:
            ix=ixs[0]
        else:
            raise RuntimeError(f'More than 1 entry for event "{eventname}" in column {self.info_colnames["Eventname"]} of info table {self.filename}')
        
        # delta time in years
        delta_t_years = Year-self.t.loc[ix,self.info_colnames['yr_expl']]
        if verbose: print(f'delta_t_years = {delta_t_years} years = {Year}-{self.t.loc[ix,self.info_colnames["yr_expl"]]}')
        
        Distance_lty = self.t.loc[ix,self.info_colnames['Distance_pc']]*pc2ly
        if verbose: print(f'Distance_lty = {Distance_lty:0.1f} lty = {self.t.loc[ix,self.info_colnames["Distance_pc"]]:0.1f} pc * {pc2ly} lty/pc')
        
        PA_deg = calcPA(self.t.loc[ix,self.info_colnames['RA']],self.t.loc[ix,self.info_colnames['Dec']],
                        RA,Dec)
        if PA_deg<0.0:PA_deg+=360.0
        if verbose: print(f'PA = {PA_deg:.02f} degree')

        angsep_deg = skydist_degree(self.t.loc[ix,self.info_colnames['RA']],self.t.loc[ix,self.info_colnames['Dec']],
                                    RA,Dec)
        if verbose: print(f'angsep_deg = {angsep_deg:.04f} degree')

        z_ly = z_ly_atD(angsep_deg,delta_t_years,Distance_lty)
        #rho_ly = rho_ly_aDz(angsep_deg,Distance_lty,z_ly)
        rho_ly = rho_ly_zt(z_ly,delta_t_years)
        x_ly = -math.sin(PA_deg*deg2rad) * rho_ly
        y_ly = math.cos(PA_deg*deg2rad) * rho_ly
        r_ly = math.sqrt(rho_ly*rho_ly+z_ly*z_ly)
        if verbose: 
            print(f'x,y,z={x_ly:.1f},{y_ly:.1f},{z_ly:.1f} ly')
            print(f'rho={rho_ly:.1f}, r={r_ly:.1f} ly')

        alpha_deg = math.atan(z_ly/rho_ly)*rad2deg
        theta_deg = math.atan(rho_ly/z_ly)*rad2deg
        if verbose: 
            print(f'alpha_deg={alpha_deg:.2f} degree')
            print(f'theta_deg={theta_deg:.2f} degree')
        
        
        results3Dinfo = {'RA':RA,
                  'Dec':Dec,
                  'delta_t_yr':delta_t_years,
                  'angsep_deg':angsep_deg,
                  'PA':PA_deg,
                  'alpha_deg':alpha_deg,
                  'theta_deg':theta_deg,
                  'x_ly':x_ly,
                  'y_ly':y_ly,
                  'z_ly':z_ly,
                  'rho_ly':rho_ly,
                  'r_ly':r_ly
                  }
        return(results3Dinfo)
        

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('eventname', help="Name of the event, e.g. Tycho. This name needs to be in the Event info table ")
    parser.add_argument('RA', help="RA of light echo")
    parser.add_argument('Dec', help="Dec of light echo")
    parser.add_argument('year', type=float, help="Year of light echo observation")
    parser.add_argument('--event_info_filename', default='EventInfo.txt', help='Event info table filename (default=%(default)s)')
    args = parser.parse_args()
    
    
    calc3D = calc3Dposclass()
    calc3D.loadEventInfo(filename=args.event_info_filename)
    calc3D.write()
    #print(calc3D.t.dtypes)
    calc3D.calc_3Dpos(args.eventname,args.RA,args.Dec,args.year)