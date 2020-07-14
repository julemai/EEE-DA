#!/usr/bin/env python
from __future__ import print_function

# Copyright 2016-2019 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../lib')

import argparse
import numpy        as np      # to perform numerics
import datetime     as datetime
import xarray   as xr

import netcdf4    as     nc4                 # in lib/
from   date2dec   import date2dec            # in lib/
from   dec2date   import dec2date            # in lib/
from   autostring import astr                # in lib/

def get_temp_prec(input_file,time_step,all_time_steps,date_time_step,period_time_step,grid_cells):

    time_var_name = 't'  # name of the time variable
    time_dim_name = 't'  # name of the time dimension

    prec_var_name = 'pTot'
    temp_var_name = 'tMoy'  # 'tMax' #
    
    ncin  = nc4.NcDataset(input_file, "r")

    # check if given time_step is valid
    ntime = len(ncin.dimensions[time_dim_name])      # number of time steps in original file

    if (time_step > ntime-1 or  -time_step > ntime ):
        raise ValueError("The time step you want to extract is not valid!")
    else:
        if all_time_steps:
            time_step_start = 0
            time_step_end   = ntime
        else:
            if (time_step > -1):
                time_step_start = time_step
                time_step_end   = time_step+1
            else:
                time_step_start = ntime-1
                time_step_end   = ntime

    if ( date_time_step != '' ):
        # if date is given, find the index of time matching this date
        time_unit = ncin.variables[time_var_name].units
        
        if ('0000-00-00' in time_unit or 'Matlab' in time_unit): 
            # Matlab time stamp
            time_in_matlab = ncin.variables[time_var_name][:]
            date_in_matlab = date2dec(eng=date_time_step,calendar='proleptic_gregorian')+367
            # print("file:                      ",input_file)
            # print("dates in file:             ",dec2date(time_in_matlab[0]-367,calendar='proleptic_gregorian',eng=True)," to ",
            #                                     dec2date(time_in_matlab[-1]-367,calendar='proleptic_gregorian',eng=True))
            # print("date specified:            ",date_time_step)
        
            if (date_in_matlab < time_in_matlab[0] or date_in_matlab > time_in_matlab[-1]):
                raise ValueError("Date you have given is not available in the file specified with option -i.")
            time_step = list(time_in_matlab).index(date_in_matlab)
            print("time step to be extracted: ",time_step)
        else:
            # Normal time
            time_in_file   = ncin.variables[time_var_name][:]                               # e.g. [0.0, 1.0, 2.0, ....]
            ref_date       = ' '.join(time_unit.split(' ')[2:])                             # e.g. u'1953-01-02 00:00:00'
            ref_date_jul   = date2dec(eng=ref_date,calendar='proleptic_gregorian')          # e.g. 712953.0
            
            date_in_julian = date2dec(eng=date_time_step,calendar='proleptic_gregorian')    # e.g. 734868.0 for 2013-01-01
            # print("file:                      ",input_file)
            # print("dates in file:             ",dec2date(time_in_file[0] +ref_date_jul,calendar='proleptic_gregorian',eng=True)," to ",
            #                                     dec2date(time_in_file[-1]+ref_date_jul,calendar='proleptic_gregorian',eng=True))
            # print("date specified:            ",date_time_step)
        
            if (date_in_julian < time_in_file[0]+ref_date_jul or date_in_julian > time_in_file[-1]+ref_date_jul):
                raise ValueError("Date you have given is not available in the file specified with option -i.")
            time_step = list(time_in_file+ref_date_jul).index(date_in_julian)
            
        print("time step to be extracted: ",time_step)       
        time_step_start = time_step
        time_step_end   = time_step+1

    elif ( period_time_step != '' ):
        # if period given find both indexes for start and end date
        time_unit = ncin.variables[time_var_name].units
        
        if ('0000-00-00' in time_unit or 'Matlab' in time_unit): 
            # Matlab time stamp
            time_in_matlab = ncin.variables[time_var_name][:]

            # (1) start day
            date_time_step_start = period_time_step.split('_')[0]
            date_in_matlab = date2dec(eng=date_time_step_start,calendar='proleptic_gregorian')+367
            # print("file:                      ",input_file)
            # print("dates in file:             ",dec2date(time_in_matlab[0]-367,calendar='proleptic_gregorian',eng=True)," to ",
            #                                     dec2date(time_in_matlab[-1]-367,calendar='proleptic_gregorian',eng=True))
            # print("date specified:            ",date_time_step_start)  
  
            if (date_in_matlab < time_in_matlab[0] or date_in_matlab > time_in_matlab[-1]):
                raise ValueError("Date you have given is not available in the file specified with option -i.")    
            time_step_start = list(time_in_matlab).index(date_in_matlab)

            # (2) end day
            date_time_step_end = period_time_step.split('_')[1]
            date_in_matlab = date2dec(eng=date_time_step_end,calendar='proleptic_gregorian')+367
            # print("file:                      ",input_file)
            # print("dates in file:             ",dec2date(time_in_matlab[0]-367,calendar='proleptic_gregorian',eng=True)," to ",
            #                                     dec2date(time_in_matlab[-1]-367,calendar='proleptic_gregorian',eng=True))
            # print("date specified:            ",date_time_step_end)    
            if (date_in_matlab < time_in_matlab[0] or date_in_matlab > time_in_matlab[-1]):
                raise ValueError("Date you have given is not available in the file specified with option -i.")    
            time_step_end = list(time_in_matlab).index(date_in_matlab)+1
            
        else:
            # Normal time
            time_in_file   = ncin.variables[time_var_name][:]                               # e.g. [0.0, 1.0, 2.0, ....]
            ref_date       = ' '.join(time_unit.split(' ')[2:])                             # e.g. u'1953-01-02 00:00:00'
            ref_date_jul   = date2dec(eng=ref_date,calendar='proleptic_gregorian')          # e.g. 712953.0

            # (1) start day
            date_time_step_start = period_time_step.split('_')[0]
            date_in_julian = date2dec(eng=date_time_step_start,calendar='proleptic_gregorian')    # e.g. 734868.0 for 2013-01-01
            # print("file:                      ",input_file)
            # print("dates in file:             ",dec2date(time_in_file[0] +ref_date_jul,calendar='proleptic_gregorian',eng=True)," to ",
            #                                     dec2date(time_in_file[-1]+ref_date_jul,calendar='proleptic_gregorian',eng=True))
            # print("date specified:            ",date_time_step_start)
        
            if (date_in_julian < time_in_file[0]+ref_date_jul or date_in_julian > time_in_file[-1]+ref_date_jul):
                raise ValueError("Date you have given is not available in the file specified with option -i.")
            time_step_start = list(time_in_file+ref_date_jul).index(date_in_julian)

            # (2) end day
            date_time_step_end = period_time_step.split('_')[1]
            date_in_julian = date2dec(eng=date_time_step_end,calendar='proleptic_gregorian')    # e.g. 734868.0 for 2013-01-01
            # print("file:                      ",input_file)
            # print("dates in file:             ",dec2date(time_in_file[0] +ref_date_jul,calendar='proleptic_gregorian',eng=True)," to ",
            #                                     dec2date(time_in_file[-1]+ref_date_jul,calendar='proleptic_gregorian',eng=True))
            # print("date specified:            ",date_time_step_end)
        
            if (date_in_julian < time_in_file[0]+ref_date_jul or date_in_julian > time_in_file[-1]+ref_date_jul):
                raise ValueError("Date you have given is not available in the file specified with option -i.")
            time_step_end = list(time_in_file+ref_date_jul).index(date_in_julian)+1

        # print extraction range (-1 because last time step won't be extracted by python)
        # print("time steps to be extracted: ",time_step_start, " to ", time_step_end-1)

    # make grid cell string an array; start now from 0
    grid_cells = np.array(grid_cells.split('[')[1].split(']')[0].split(','),dtype=np.int)-1

    # identify which axis is time; other one will be basins
    if ( len(ncin.variables[prec_var_name].dimensions) != 2 ):
        raise ValueError('Precipitation variable has more than 2 dimensions (should be only time and basins)')
    else:
        idx_t = np.where(np.array(ncin.variables[prec_var_name].dimensions) == time_dim_name)[0][0]

        if idx_t == 0:
            # -4 because we might need prec of 4 days before as well
            prec = ncin.variables[prec_var_name][max(time_step_start-4,0):time_step_end-1,grid_cells]
            # average over grid cells
            prec = np.mean(prec,axis=1)
        elif idx_t == 1:
            # -4 because we might need prec of 4 days before as well
            prec = ncin.variables[prec_var_name][grid_cells,max(time_step_start-4,0):time_step_end-1]
            # average over grid cells
            prec = np.mean(prec,axis=0)
        else:
            raise ValueError('Precipitation: Dimension that contains time needs to be 0 or 1.')
    if ( len(ncin.variables[temp_var_name].dimensions) != 2 ):
        raise ValueError('Temperature variable has more than 2 dimensions (should be only time and basins)')
    else:
        idx_t = np.where(np.array(ncin.variables[temp_var_name].dimensions) == time_dim_name)[0][0]

        if idx_t == 0:
            temp = ncin.variables[temp_var_name][time_step_start:time_step_end,grid_cells]
            # average over grid cells
            temp = np.mean(temp,axis=1)
        elif idx_t == 1:
            temp = ncin.variables[temp_var_name][grid_cells,time_step_start:time_step_end]
            # average over grid cells
            temp = np.mean(temp,axis=0)
        else:
            raise ValueError('Temperature: Dimension that contains time needs to be 0 or 1.')


    return [time_step_start,time_step_end,prec,temp]

    


def categorize_meteo_decision_tree_expert(input_file='', 
                                     time_step=-1, date_time_step='',
                                     period_time_step='', all_time_steps=False,
                                     grid_cells='[1]',
                                     parameters=False,parameter_ranges=False,
                                     return_all_para_info=False,
                                     return_time=False,
                                     return_all_recipe_info=False,
                                     recipe=False):

    # categorizer used for paper
    # recipe 1: expert knowledge ==> recipe #5 from categorize_meteo_decision_tree_1 (performed best; see EGU 2019 presentation)
    # recipe 2: identified using EEE

    # initial values need to be chosen such that this adjustment is the one-operator, i.e., X = initial * X
    ranges = {
            #          min ,   max , ini
            'v1': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of niveauEauSol
            'v2': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of niveauEauNappe
            'v3': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of G (previously stockNeige)
            'v4': [  -5.0 ,  5.0 , 0.0 ],   #           intercept in linear adjustment of eTg (previously indexMurissementNeige)
            'v5': [ -10.0 , 10.0 , 0.0 ],   # removed from code (previously indexTempNeige)
            'v6': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of precipitation (a in a + b x prec)
            'v7': [   0.0 ,  2.0 , 1.0 ],   # not used; slope     in linear adjustment of precipitation (b in a + b x prec)
            'v8': [  -5.0 ,  5.0 , 0.0 ],   #           intercept in linear adjustment of temperature   (a in a + b x prec)
            'v9': [   0.0 ,  2.0 , 1.0 ]    # not used; slope     in linear adjustment of temperature   (b in a + b x prec)
            }

    recipes = {
        #            C1               C2              C3                          C4                          C5                          C6                          C7
        # recipe 1: expert knowledge ==> recipe #1 from categorize_meteo_decision_tree_1 (performed worst; see EGU 2019 presentation)
        # recipe 5: expert knowledge ==> recipe #5 from categorize_meteo_decision_tree_1 (performed best; see EGU 2019 presentation)
        'recipe_1':  [    [ 'v1' ],
                          [ 'v1' ],
                          [ 'v3' ],
                          [ 'v3' ],
                          [ 'v2' ],
                          [ 'v2' ],
                          [ 'v2' ]  ],
        'recipe_2':  [    [ 'v2' ],
                          [ 'v2' ],
                          [ 'v4' ],
                          [ 'v4' ],
                          [ 'v1' ],
                          [ 'v1' ],
                          [ 'v1' ]  ],
        'recipe_3':  [    [ 'v1', 'v2' ],
                          [ 'v1', 'v2' ],
                          [ 'v3', 'v4' ],
                          [ 'v3', 'v4' ],
                          [ 'v1', 'v2' ],
                          [ 'v1', 'v2' ],
                          [ 'v1', 'v2' ]  ],
		'recipe_5':  [    [ 'v2' ],
                          [ 'v1' ],
                          [ 'v3', 'v4', 'v6', 'v8' ],
                          [ 'v3', 'v4', 'v6', 'v8' ],
                          [ 'v1', 'v2', 'v6', 'v8' ],
                          [ 'v1', 'v2', 'v6', 'v8' ],
                          [ 'v1', 'v2', 'v6', 'v8' ]  ]
        }

    if return_all_para_info: # return only information on all parameters

        parameters       = [ pp         for pp in ranges ]
        parameter_ranges = [ ranges[pp] for pp in ranges ]
        return [parameters, parameter_ranges]

    elif return_all_recipe_info: # return only recipe information

        rrecipes      = [ rr          for rr in recipes ]
        rrecipes_info = [ recipes[rr] for rr in recipes ]
        return [rrecipes, rrecipes_info]

    else:   # categorize and return parameter info of parameters active

        if not(recipe):
            raise ValueError("Recipe must be given, e.g. 'recipe_1'")

        # -----------------------------------------------
        # get the meteo from file
        # -----------------------------------------------
        [time_step_start,time_step_end,prec,temp] = get_temp_prec(input_file,time_step,all_time_steps,date_time_step,period_time_step,grid_cells)

        # -------------
        # This is how is should be...
        # -------------
        ttime = xr.open_dataset(input_file)['t'][time_step_start:time_step_end].data
        ttime = np.array( [ datetime.datetime(int(str(tt)[0:4]), int(str(tt)[5:7]), int(str(tt)[8:10]), int(str(tt)[11:13]), int(str(tt)[14:16])) for tt in ttime ])

        ttime_months = [ tt.month for tt in ttime ]

        # initialize categories
        icategories       = [ -9999 for tt in range(time_step_end-time_step_start)]
        iparameters       = [ -9999 for tt in range(time_step_end-time_step_start)]
        iparameter_ranges = [ [] for tt in range(time_step_end-time_step_start)]

        for tt in range(time_step_end-time_step_start):
            
            if temp[tt] <= -10.0:
                icategories[tt] = 1
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 1 # ['v1', 'v2']  # 'v2' is new
            elif temp[tt] <= -5.0:
                icategories[tt] = 2
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 2 # ['v1', 'v2']
            elif temp[tt] <= 5.0 and ttime_months[tt] >= 3 and ttime_months[tt] <= 5:
                icategories[tt] = 3
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 3 # ['v1', 'v4', 'v6', 'v8']
            elif ttime_months[tt] >= 3 and ttime_months[tt] <= 5:
                icategories[tt] = 4
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 4 # ['v3', 'v4', 'v6', 'v8']
            elif ttime_months[tt] >= 6 and ttime_months[tt] <= 11:
                # 4 because of offset of precipitation values read in
                if np.sum(prec[tt:tt+4]) > 20.0:
                    icategories[tt] = 5
                    iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 5 # ['v1', 'v2']  # 'v2' is new
                else:
                    icategories[tt] = 6
                    iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 6 # ['v1', 'v2', 'v6', 'v8']  # 'v2' is new
            else:
                icategories[tt] = 7
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  # Category 7 # ['v2', 'v6']
                
            iparameter_ranges[tt] = [ ranges[pp] for pp in iparameters[tt] ]

        out = []
        out += [icategories]
    
        if parameters:       out += [iparameters]
        if parameter_ranges: out += [iparameter_ranges]
        if return_time:      out += [ttime]

        if len(out) == 1:
            return np.array(out[0])
        else:
            return out

def categorize_meteo_decision_tree_EEE(input_file='', 
                                     time_step=-1, date_time_step='',
                                     period_time_step='', all_time_steps=False,
                                     grid_cells='[1]',
                                     parameters=False,parameter_ranges=False,
                                     return_all_para_info=False,
                                     return_time=False,
                                     return_all_recipe_info=False,
                                     recipe=False):

    # categorizer used for paper
    # recipe 1: identified using EEE

    # initial values need to be chosen such that this adjustment is the one-operator, i.e., X = initial * X
    ranges = {
            #          min ,   max , ini
            'v1': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of niveauEauSol
            'v2': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of niveauEauNappe
            'v3': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of G (previously stockNeige)
            'v4': [  -5.0 ,  5.0 , 0.0 ],   #           intercept in linear adjustment of eTg (previously indexMurissementNeige)
            'v5': [ -10.0 , 10.0 , 0.0 ],   # removed from code (previously indexTempNeige)
            'v6': [ -25.0 , 25.0 , 0.0 ],   #           intercept in linear adjustment of precipitation (a in a + b x prec)
            'v7': [   0.0 ,  2.0 , 1.0 ],   # not used; slope     in linear adjustment of precipitation (b in a + b x prec)
            'v8': [  -5.0 ,  5.0 , 0.0 ],   #           intercept in linear adjustment of temperature   (a in a + b x prec)
            'v9': [   0.0 ,  2.0 , 1.0 ]    # not used; slope     in linear adjustment of temperature   (b in a + b x prec)
            }

    recipes = {
        # recipe 1: identified using EEE
        'recipe_1':  [      ['v1','v2'],                         # C1:  DOY in [350, 70)
                            ['v1','v2','v8'],                    # C2:  DOY in [ 70, 83)
                            ['v1','v2','v8','v4'],               # C3:  DOY in [ 83, 95)
                            ['v1','v2','v4','v6','v8'],          # C4:  DOY in [ 95,116)
                            ['v1','v2','v3','v4','v6','v8'],     # C5:  DOY in [116,143)
                            ['v1','v2','v3','v6','v8'],          # C6:  DOY in [143,155)
                            ['v1','v2','v3','v6'],               # C7:  DOY in [155,166)
                            ['v1','v2','v6'],                    # C8:  DOY in [166,260)
                            ['v1','v2','v6','v8'],               # C9:  DOY in [260,275)
                            ['v1','v2','v3','v6','v8'],          # C10: DOY in [275,319)
                            ['v1','v2','v3','v8'],               # C11: DOY in [319,325)
                            ['v1','v2','v8'],                    # C12: DOY in [325,350)
                     ],
        # recipe 2: recipe 1 identified using EEE but only most informative parameters v1 and v6 (maybe try also to add v3 and then v8)
        'recipe_2':  [      ['v1'],                              # C1:  DOY in [350, 70)
                            ['v1'],                              # C2:  DOY in [ 70, 83)
                            ['v1'],                              # C3:  DOY in [ 83, 95)
                            ['v1','v6'],                         # C4:  DOY in [ 95,116)
                            ['v1','v6'],                         # C5:  DOY in [116,143)
                            ['v1','v6'],                         # C6:  DOY in [143,155)
                            ['v1','v6'],                         # C7:  DOY in [155,166)
                            ['v1','v6'],                         # C8:  DOY in [166,260)
                            ['v1','v6'],                         # C9:  DOY in [260,275)
                            ['v1','v6'],                         # C10: DOY in [275,319)
                            ['v1'],                              # C11: DOY in [319,325)
                            ['v1'],                              # C12: DOY in [325,350)
                     ],
        # recipe 3: recipe 1 identified using EEE but only most informative parameters of state variables (no forcings) --> v1 and v3
        'recipe_3':  [      ['v1'],                              # C1:  DOY in [350, 70)
                            ['v1'],                              # C2:  DOY in [ 70, 83)
                            ['v1'],                              # C3:  DOY in [ 83, 95)
                            ['v1'],                              # C4:  DOY in [ 95,116)
                            ['v1','v3'],                         # C5:  DOY in [116,143)
                            ['v1','v3'],                         # C6:  DOY in [143,155)
                            ['v1','v3'],                         # C7:  DOY in [155,166)
                            ['v1'],                              # C8:  DOY in [166,260)
                            ['v1'],                              # C9:  DOY in [260,275)
                            ['v1','v3'],                         # C10: DOY in [275,319)
                            ['v1','v3'],                         # C11: DOY in [319,325)
                            ['v1'],                              # C12: DOY in [325,350) 
                     ],
        # recipe 4: recipe 1 identified using EEE but without forcings adjustments --> v1, v2, v3, v4
        'recipe_4':  [      ['v1','v2'],                         # C1:  DOY in [350, 70)
                            ['v1','v2'],                         # C2:  DOY in [ 70, 83)
                            ['v1','v2','v4'],                    # C3:  DOY in [ 83, 95)
                            ['v1','v2','v4'],                    # C4:  DOY in [ 95,116)
                            ['v1','v2','v3','v4'],               # C5:  DOY in [116,143)
                            ['v1','v2','v3'],                    # C6:  DOY in [143,155)
                            ['v1','v2','v3'],                    # C7:  DOY in [155,166)
                            ['v1','v2'],                         # C8:  DOY in [166,260)
                            ['v1','v2'],                         # C9:  DOY in [260,275)
                            ['v1','v2','v3'],                    # C10: DOY in [275,319)
                            ['v1','v2','v3'],                    # C11: DOY in [319,325)
                            ['v1','v2'],                         # C12: DOY in [325,350)
                     ]
        }

    if return_all_para_info: # return only information on all parameters

        parameters       = [ pp         for pp in ranges ]
        parameter_ranges = [ ranges[pp] for pp in ranges ]
        return [parameters, parameter_ranges]

    elif return_all_recipe_info: # return only recipe information

        rrecipes      = [ rr          for rr in recipes ]
        rrecipes_info = [ recipes[rr] for rr in recipes ]
        return [rrecipes, rrecipes_info]

    else:   # categorize and return parameter info of parameters active

        if not(recipe):
            raise ValueError("Recipe must be given, e.g. 'recipe_1'")

        # -----------------------------------------------
        # get the meteo from file
        # -----------------------------------------------
        [time_step_start,time_step_end,prec,temp] = get_temp_prec(input_file,time_step,all_time_steps,date_time_step,period_time_step,grid_cells)

        # -------------
        # This is how is should be...
        # -------------
        ttime = xr.open_dataset(input_file)['t'][time_step_start:time_step_end].data
        ttime = np.array( [ datetime.datetime(int(str(tt)[0:4]), int(str(tt)[5:7]), int(str(tt)[8:10]), int(str(tt)[11:13]), int(str(tt)[14:16])) for tt in ttime ])

        ttime_months = [ tt.month for tt in ttime ]

        # initialize categories
        icategories       = [ -9999 for tt in range(time_step_end-time_step_start)]
        iparameters       = [ -9999 for tt in range(time_step_end-time_step_start)]
        iparameter_ranges = [ [] for tt in range(time_step_end-time_step_start)]

        for tt in range(time_step_end-time_step_start):

            doy_tt = ttime[tt].timetuple().tm_yday
            
            if (doy_tt >= 350 or doy_tt < 70):
                icategories[tt] = 1
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  
            elif (doy_tt >= 70 and doy_tt < 83):
                icategories[tt] = 2
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  
            elif (doy_tt >= 83 and doy_tt < 95):
                icategories[tt] = 3
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  
            elif (doy_tt >= 95 and doy_tt < 116):
                icategories[tt] = 4
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  
            elif (doy_tt >= 116 and doy_tt < 143):
                icategories[tt] = 5
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 143 and doy_tt < 155):
                icategories[tt] = 6
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 155 and doy_tt < 166):
                icategories[tt] = 7
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 166 and doy_tt < 260):
                icategories[tt] = 8
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 260 and doy_tt < 275):
                icategories[tt] = 9
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 275 and doy_tt < 319):
                icategories[tt] = 10
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 319 and doy_tt < 325):
                icategories[tt] = 11
                iparameters[tt] = recipes[recipe][icategories[tt]-1]
            elif (doy_tt >= 325 and doy_tt < 350):
                icategories[tt] = 12
                iparameters[tt] = recipes[recipe][icategories[tt]-1]  
            else:
                print('doy  = ',doy_tt)
                print('temp = ',temp[tt])
                raise ValueError('Category not known!')
                
            iparameter_ranges[tt] = [ ranges[pp] for pp in iparameters[tt] ]

        out = []
        out += [icategories]
    
        if parameters:       out += [iparameters]
        if parameter_ranges: out += [iparameter_ranges]
        if return_time:      out += [ttime]

        if len(out) == 1:
            return np.array(out[0])
        else:
            return out
        

if __name__ == '__main__':
    
    verbose          = False
    input_file       = ''
    time_step        = -1
    date_time_step   = ''
    period_time_step = ''
    all_time_steps   = False
    grid_cells       = '[1]'
    recipe           = False
    parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                      description='''Extracts last time step from input file and writes it to output file. The variables and values are unchanged in the output file.''')
    parser.add_argument('-i', '--input_file', action='store',
                        default=input_file, dest='input_file', metavar='input_file', nargs=1,
                        help='Name of input file (standard is a CEQUEAU file with groups and variables).')
    parser.add_argument('-t', '--time_step', action='store',
                        default=time_step, dest='time_step', metavar='time_step',
                        help='Time step which should be extracted. The default is the last time step (-1). Be aware that the numbering starts with 0.')
    parser.add_argument('-d', '--date_time_step', action='store',
                        default=date_time_step, dest='date_time_step', metavar='date_time_step',
                        help='Date of the time step which should be extracted (format: YYYY-MM-DD).')
    parser.add_argument('-p', '--period_time_step', action='store',
                        default=period_time_step, dest='period_time_step', metavar='period_time_step',
                        help='Period of time steps which should be extracted (format: YYYY-MM-DD_YYYY-MM-DD).')
    parser.add_argument('-a', '--all_time_steps', action='store_true',
                        default=all_time_steps, dest='all_time_steps', 
                        help='If -a used, all time steps will be categorized.')
    parser.add_argument('-g', '--grid_cells', action='store',
                        default=grid_cells, dest='grid_cells', metavar='grid_cells',
                        help='Grid cells that should be considered. List with comma as separator. E.g., [20,21,22,23]. Numbering starts with 1.')
    parser.add_argument('-r', '--recipe', action='store',
                        default=recipe, dest='recipe', metavar='recipe',
                        help='Recipe to use for decision variables. E.g. "recipe_1"')

    # example:
    #     all timesteps for ID=2 = PERIB = Peribonka
    #         run categorize_meteo_decision_tree.py -i ../../data/meteo.nc -g [6,7] -a -r recipe_1
    #     certain timestep for ID=2 = PERIB = Peribonka
    #         run categorize_meteo_decision_tree.py -i ../../data/meteo.nc -g [6,7] -t 2 -r recipe_1
    #
    # python categorize_meteo_decision_tree.py -i ../../data/meteo.nc -g [1,2,3,4,5] -t -1 -a -r recipe_1
    # python categorize_meteo_decision_tree.py -i ../../data/meteo.nc -g [6,7]       -t -1 -a -r recipe_1
    # python categorize_meteo_decision_tree.py -i ../../data/meteo.nc -g [8,9,10]    -t -1 -a -r recipe_1
    # python categorize_meteo_decision_tree.py -i ../../data/meteo.nc -g [11]        -t -1 -a -r recipe_1

    args             = parser.parse_args()
    input_file       = args.input_file[0]
    grid_cells       = args.grid_cells
    all_time_steps   = args.all_time_steps
    time_step        = np.int(args.time_step)
    date_time_step   = args.date_time_step
    period_time_step = args.period_time_step
    recipe           = args.recipe


    categories = categorize_meteo_decision_tree_EEE(input_file=input_file, 
                                               time_step=time_step, date_time_step=date_time_step,
                                               period_time_step=period_time_step, all_time_steps=all_time_steps,
                                               grid_cells=grid_cells, recipe=recipe)

    if all_time_steps:
        count = { cc: np.shape(np.where(categories==cc)[0])[0] for cc in np.unique(categories) }
        percentage = { cc: np.shape(np.where(categories==cc)[0])[0]/np.float(len(categories))*100.0 for cc in np.unique(categories) }
        print("count:      ",count)
        print("percentage: ",percentage)

        print("")
        for cc in np.unique(categories):
            print(str(count[cc])+" & ("+astr(percentage[cc],prec=0)+"\%) & ")

    else:
        print("")
        print('category: ',categories, '   (numbers start with 1)')




