#!/usr/bin/env python
from __future__ import print_function

# Copyright 2019-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
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
#    to dump results to NetCDF
#    -----------------------------
#    python figure_8.py -o "/Users/j6mai/Documents/GitHub/CRD-DA/scripts/forward_run/open-loop_performance/meteoCorrected/" -i "[forecast_performance_EEE.nc, forecast_performance_expert.nc]" -b "[Passes Dangereuses, Peribonka, Lac Manouane, Montagnes Blanches]" -p figure_8.pdf -t pdf   --only_summer
#
#    to read results from NetCDF
#    -----------------------------
#    python figure_8.py -i "[../data/forecast/forecast_performance_EEE.nc, ../data/forecast/forecast_performance_expert.nc]" -b "[Passes Dangereuses, Peribonka, Lac Manouane, Montagnes Blanches]" -p figure_8.pdf -t pdf   --only_summer

import argparse

dolog            = False
dowhite          = False
verbose          = False
openloop_folder  = ['']
basin_name       = None
plotname         = ''
outtype          = ''
usetex           = False
serif            = False
only_summer      = False
infiles          = ''

parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Plots results of exeriment for Engage Plus (calibration results and adjusted optimal parameters).''')
parser.add_argument('-o', '--openloop_folder', action='store',
                    default=openloop_folder, dest='openloop_folder', metavar='openloop_folder', nargs=1,
                    help='Name of input folder containing all open-loop performance measures (expected to have subfolders named with basins short names, e.g. ["ASAM", "MBLANC"].')
parser.add_argument('-i', '--infiles', action='store',
                    default=infiles, dest='infiles', metavar='infiles', nargs=1,
                    help='Name of input files, e.g. "[initial_state_adjustments_EEE.nc, initial_state_adjustments_expert.nc]".')
parser.add_argument('-b', '--basin_name', action='store',
                    default=basin_name, dest='basin_name', metavar='basin_name',
                    help='Name of the basin, e.g. ["Ashuapmushuan Amont","Montagnes Blanches"].')
parser.add_argument('-p', '--plotname', action='store',
                    default=plotname, dest='plotname', metavar='plotname',
                    help='Name of plot output file for types pdf, html or d3, '
                    'and name basis for type png (default: '+__file__[0:__file__.rfind(".")]+').')
parser.add_argument('-s', '--serif', action='store_true', default=serif, dest="serif",
                    help="Use serif font; default sans serif.")
parser.add_argument('-t', '--type', action='store',
                    default=outtype, dest='outtype', metavar='outtype',
                    help='Output type is pdf, png, html, or d3 (default: open screen windows).')
parser.add_argument('-u', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf, png and html.")
parser.add_argument('-c', '--only_summer', action='store_true',
                    default=only_summer, dest="only_summer",
                    help='Only summer data (July 1 to Nov 30) are used. Default: False')

args             = parser.parse_args()
openloop_folder  = args.openloop_folder[0]
basin_name       = args.basin_name
plotname         = args.plotname
outtype          = args.outtype
serif            = args.serif
usetex           = args.usetex
only_summer      = args.only_summer
infiles          = args.infiles[0]

del parser, args

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/../scripts')
sys.path.append(dir_path+'/../scripts/lib')

# -----------------------
# Load additional packages
# -----------------------
import numpy as np
import datetime
import collections
import glob
import time
import netCDF4 as nc4
import color                                                      # in lib/
from abc2plot      import abc2plot                                # in lib/
from position      import position                                # in lib/
from str2tex       import str2tex                                 # in lib/
from autostring    import astr                                    # in lib/

rename_recipes_EEE    = { 'recipe_1': 'recipe_1',    # all states/inputs
                          'recipe_2': None,          # not used
                          'recipe_3': 'recipe_3',    # most sensitive states
                          'recipe_4': 'recipe_2'     # all states             (rename actual recipe_4 to recipe_2)
                        }
rename_recipes_expert = { 'recipe_1': 'recipe_3',    # one state at a time    (rename actual recipe_1 to recipe_3)
                          'recipe_2': None,          # not used
                          'recipe_3': 'recipe_2',    # all states             (rename actual recipe_3 to recipe_2)
                          'recipe_5': 'recipe_1'     # all states/inputs      (rename actual recipe_5 to recipe_1)
                        }
rename_vars           = { 'v_1': 'v_1', 'v1': 'v1',              
                          'v_2': 'v_2', 'v2': 'v2',
                          'v_3': 'v_3', 'v3': 'v3',
                          'v_4': 'v_4', 'v4': 'v4',
                          'v_6': 'v_5', 'v6': 'v5',  # (rename actual v_6 to v_5)
                          'v_8': 'v_6', 'v8': 'v6'   # (rename actual v_8 to v_6)
                        }

basins = collections.OrderedDict()
basins['Passes Dangereuses'] = ['PD',     1, [1,2,3,4,5]]
basins['Peribonka']          = ['PERIB',  2, [6,7]]
basins['Lac Manouane']       = ['LM',     3, [8,9,10]]
basins['Montagnes Blanches'] = ['MBLANC', 4, [11]]

if (basin_name == None):
    basin_name = [ bb for bb in basins ]
else:
    # convert string inputs into lists
    basin_name = basin_name.split('[')[1].split(']')[0].split(',')
    basin_name = [ ii.strip() for ii in basin_name ]

# convert string inputs into lists
infiles = infiles.split('[')[1].split(']')[0].split(',')
infiles = [ ii.strip() for ii in infiles ]
netcdf_file_EEE          = infiles[0]
netcdf_file_expert       = infiles[1]

# reference date for all time variables in NetCDF files
refdate = datetime.datetime(1950,1,1,0,0)

# which decision tree is used
from categorize_meteo_decision_tree import categorize_meteo_decision_tree_expert as categorizer_expert
from categorize_meteo_decision_tree import categorize_meteo_decision_tree_EEE as categorizer_EEE


print('-----------------------------------------------------------------')
print('READ: categorize_meteo_decision_tree_EEE')

input_folder_EEE = "/Users/j6mai/Documents/GitHub/CRD-DA/scripts/engage_plus_runs/categorize_meteo_decision_tree_EEE"
recipe = '[recipe_1, recipe_2, recipe_3, recipe_4]'
recipe = '[recipe_1, recipe_4, recipe_3]'

# convert string inputs into lists
recipe_EEE = recipe.split('[')[1].split(']')[0].split(',')
recipe_EEE = [ ii.strip() for ii in recipe_EEE ]
nrecipe_EEE = len(recipe_EEE)
    
nbasins = len(basin_name)


# overwrite read_from_netcdf if file does not exist
if not(os.path.exists(netcdf_file_EEE)):
    read_from_netcdf = False
else:
    read_from_netcdf = True
    
    
obj_before_EEE          = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
obj_after_EEE           = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
para_before_EEE         = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
para_after_EEE          = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
para_active_EEE         = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
category_EEE            = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
simtime_EEE             = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
objective_EEE           = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   
objective_open_loop_EEE = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_EEE) ]   

para_info_EEE   = categorizer_EEE(return_all_para_info=True)

recipe_info_tmp_EEE  = categorizer_EEE(return_all_recipe_info=True)
recipe_info_EEE = {}
for ir,rr in enumerate(recipe_info_tmp_EEE[0]):
    recipe_info_EEE[rr] = recipe_info_tmp_EEE[1][ir]

if not(read_from_netcdf): 

    # write to netcdf
    opt_results_EEE  = nc4.Dataset(netcdf_file_EEE,  "w", format="NETCDF4")

    for iirecipe, irecipe in enumerate(recipe_EEE):

        print('--------------------------------------')
        print('Read recipe: ',irecipe)

        grp_recipe = opt_results_EEE.createGroup(irecipe)

        for iibasin, ibasin in enumerate(basin_name):

            # create group
            grp_basin = opt_results_EEE.createGroup(irecipe+'/'+ibasin)

            # path = os.path.normpath(os.path.join(dir_path,'../engage_plus_runs',pre_string+categorizer.__name__,basins[ibasin][0]))
            path_forecast = os.path.normpath(os.path.join(input_folder_EEE,irecipe,basins[ibasin][0]))
            path_openloop = os.path.normpath(os.path.join(openloop_folder,basins[ibasin][0]))

            # files and number of files
            adjust_token = np.sort(glob.glob(os.path.join(path_forecast,'*','adjustments.token')))
            nfiles = len(adjust_token)
            objectives_token = np.sort(glob.glob(os.path.join(path_forecast,'*','simulation','objectives.out')))
            if (nfiles > len(objectives_token)): # should be the same but sometimes objective was already written but calibration not finished and hence adjustements file is not existing
                raise ValueError('Number of adjustment file does not match number of files containing objective function values.')

            # save all the simulation start times T
            isimtime_EEE = np.array([ datetime.datetime(np.int(aa.split('/')[-2].split('-')[0]),
                                                    np.int(aa.split('/')[-2].split('-')[1]),
                                                    np.int(aa.split('/')[-2].split('-')[2]),0,0) for aa in adjust_token ])

            # read in one file to see number of parameters
            f = open(adjust_token[0], 'r')
            content = f.readlines()
            f.close()
            idx_v = np.where([ cc.startswith('v') for cc in content ])[0]
            npara = np.int(np.shape(idx_v)[0]/2)

            # read one file to see number of objectives
            f = open(objectives_token[0], 'r')
            content = f.readlines()
            f.close()
            nobj = np.shape(content)[0]
            objective_info = collections.OrderedDict()
            for cc in content:
                ccc = cc.split(',')
                ccc[2] = ccc[2].rstrip()   # remove trailing newline \n
                if 'week(s' in ccc[2]:
                    # track number of weeks to normalize objectives per week
                    nweeks = np.int(ccc[2].strip().split('week(s')[0].strip().split(' ')[-1])
                else:
                    nweeks = 1

                # if string too long add newline
                if len(ccc[2]) > 50:
                    tmp = ccc[2].split()
                    split_idx = np.sum(np.cumsum([ len(tt) for tt in tmp ]) < 35)
                    ccc[2] = ' '.join(tmp[0:split_idx])+'\n'+' '.join(tmp[split_idx:])
                    
                objective_info[ccc[0]] = [ ccc[2], nweeks ]

            # allocate result vectors and matrices
            iobj_before_EEE          = np.ones(nfiles) * -9999.0
            iobj_after_EEE           = np.ones(nfiles) * -9999.0
            ipara_before_EEE         = np.ones([nfiles,npara]) * -9999.0
            ipara_after_EEE          = np.ones([nfiles,npara]) * -9999.0
            ipara_active_EEE         = np.ones([nfiles,npara],dtype=bool)
            icalibrated_EEE          = np.ones(nfiles,dtype=bool)
            icategory_EEE            = np.ones(nfiles,dtype=int) * -9999
            iobjective_EEE           = np.ones([nfiles,nobj],dtype=float) * -9999.0
            iobjective_open_loop_EEE = np.ones([nfiles,nobj],dtype=float) * -9999.0

            # create dimensions
            dim_npara    = grp_basin.createDimension("npara" , npara)
            dim_ntime    = grp_basin.createDimension("ntime" , nfiles)
            dim_nobj     = grp_basin.createDimension("nobj"  , nobj)

            # create variables
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/obj_before_EEE'          , "f4",("ntime"),         zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/obj_after_EEE'           , "f4",("ntime"),         zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/para_before_EEE'         , "f4",("ntime","npara"), zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/para_after_EEE'          , "f4",("ntime","npara"), zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/para_active_EEE'         , "i4",("ntime","npara"), zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/calibrated_EEE'          , "i4",("ntime"),         zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/category_EEE'            , "i4",("ntime"),         zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/objective_EEE'           , "f4",("ntime","nobj"),  zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/objective_open_loop_EEE' , "f4",("ntime","nobj"),  zlib=True)
            grp_var = opt_results_EEE.createVariable(irecipe+'/'+ibasin+'/simtime_EEE'             , "i4",("ntime"),         zlib=True)
            grp_var.setncattr("long_name",     "time")
            grp_var.setncattr("units",         "days since "+str(refdate))
            grp_var.setncattr("calendar",      "gregorian")
            grp_var.setncattr("standard_name", "time")

            for ia, aa in enumerate(adjust_token):

                f = open(aa, 'r')
                content_adjust = f.readlines()
                f.close()

                calib_token = glob.glob(os.path.join(os.path.dirname(aa),'calib.token'))
                if len(calib_token) == 0:
                    icalibrated_EEE[ia] = False
                else:
                    icalibrated_EEE[ia] = True
                    f = open(calib_token[0], 'r')
                    content_calib = f.readlines()
                    f.close()
                    dict_calib = {}
                    for cc in content_calib:
                        ccc = cc.strip().split(',')
                        dict_calib[ccc[0]] = ','.join(ccc[1:])
                    icategory_EEE[ia] = np.int(dict_calib['category'])

                idx_obj_before_EEE = np.where([ 'obj_before' in cc for cc in content_adjust ])[0][0]
                idx_obj_after_EEE  = np.where([ 'obj_after'  in cc for cc in content_adjust ])[0][0]

                iobj_before_EEE[ia] = np.float(content_adjust[idx_obj_before_EEE].strip().split(',')[1])
                iobj_after_EEE[ia]  = np.float(content_adjust[idx_obj_after_EEE].strip().split(',')[1])

                for pp in range(npara):
                    idx_v = np.where([ cc.startswith('v'+str(pp+1)+',') for cc in content_adjust ])[0]

                    ipara_before_EEE[ia,pp] = np.float(content_adjust[idx_v[0]].strip().split(',')[1])
                    ipara_after_EEE[ia,pp]  = np.float(content_adjust[idx_v[1]].strip().split(',')[1])
                    ipara_active_EEE[ia,pp] = (content_adjust[idx_v[0]].strip().split(',')[2]=='True')

            # performance measures of forecast
            for io, oo in enumerate(objectives_token):

                if io < nfiles: # make sure that calibration was finished and adjustements.token exists
                    f = open(oo, 'r')
                    content_objectives = f.readlines()
                    f.close()
                    
                    iobjective_EEE[io,0:nobj] = np.array([ np.float(cc.split(',')[1]) for cc in content_objectives ])

            # performance measures of open-loop simulation
            for io, oo in enumerate(objectives_token):

                if io < nfiles: # make sure that calibration was finished and adjustements.token exists
                    current_day = oo.split('/')[-3]
                    oo_openloop = os.path.normpath(os.path.join(path_openloop,current_day,'objectives.out'))
                    f = open(oo_openloop, 'r')
                    content_objectives = f.readlines()
                    f.close()
                    
                    iobjective_open_loop_EEE[io,0:nobj] = np.array([ np.float(cc.split(',')[1]) for cc in content_objectives ])

            # apply mask on runs that were not calibrated
            iobj_before_EEE  = np.ma.array(iobj_before_EEE,  mask=~icalibrated_EEE)
            iobj_after_EEE   = np.ma.array(iobj_after_EEE,   mask=~icalibrated_EEE)
            #
            # only parameters that are calibrated (and active) 
            ipara_before_EEE = np.ma.array(ipara_before_EEE, mask=~(np.transpose(np.array([ icalibrated_EEE for ii in range(npara) ])) & ipara_active_EEE))
            ipara_after_EEE  = np.ma.array(ipara_after_EEE,  mask=~(np.transpose(np.array([ icalibrated_EEE for ii in range(npara) ])) & ipara_active_EEE))
            #
            # only mask out runs where no calibration happened
            ipara_active_EEE = np.ma.array(ipara_active_EEE, mask=~(np.transpose(np.array([ icalibrated_EEE for ii in range(npara) ]))))
            #
            icategory_EEE    = np.ma.array(icategory_EEE,    mask=~icalibrated_EEE)
            isimtime_EEE_nomask = copy.deepcopy(isimtime_EEE)
            isimtime_EEE     = np.ma.array(isimtime_EEE,     mask=~icalibrated_EEE)

            # some summary stats printed to screen
            print('   Basin: ',ibasin)
            print('   ---------------------------')
            print('   summer period (Jul-Nov)')
            print('   ---------------------------')
            tt_idx=((np.array([dd.month if dd != -9999 else -9999 for dd in isimtime_EEE.data ]) >= 7) & (np.array([dd.month if dd != -9999 else 9999 for dd in isimtime_EEE.data ]) <= 11))
            print('      number days                  calibrated: ',sum(icalibrated_EEE[tt_idx]),   '   --> initial Q(t) was more than 5% off')
            print('      number days     successfully calibrated: ', np.sum(iobjective_EEE[(icalibrated_EEE&tt_idx),2]<0.05),' (',100.0*np.sum(iobjective_EEE[(icalibrated_EEE&tt_idx),2]<0.05)/np.sum((icalibrated_EEE&tt_idx)),'%)   --> final Q(t) is less than 5% off')
            print('      number days with optimal initial states: ', np.sum(iobjective_EEE[tt_idx,2]<0.05),' (',100.0*np.sum(iobjective_EEE[tt_idx,2]<0.05)/len(icalibrated_EEE[tt_idx]),'%)   --> final Q(t) is less than 5% off')
            print('      number days              not calibrated: ', sum(~icalibrated_EEE[tt_idx]), '   --> initial Q(t) was less than 5% off')

            obj_before_EEE[iirecipe][iibasin]          = np.ma.round(iobj_before_EEE,6)  
            obj_after_EEE[iirecipe][iibasin]           = np.ma.round(iobj_after_EEE,6)   
            para_before_EEE[iirecipe][iibasin]         = np.ma.round(ipara_before_EEE,6) 
            para_after_EEE[iirecipe][iibasin]          = np.ma.round(ipara_after_EEE,6)  
            para_active_EEE[iirecipe][iibasin]         = ipara_active_EEE 
            category_EEE[iirecipe][iibasin]            = icategory_EEE
            simtime_EEE[iirecipe][iibasin]             = isimtime_EEE
            objective_EEE[iirecipe][iibasin]           = np.ma.round(iobjective_EEE,6)
            objective_open_loop_EEE[iirecipe][iibasin] = np.ma.round(iobjective_open_loop_EEE,6)

            # write float variables
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['obj_before_EEE'         ][:]    = np.ma.round(iobj_before_EEE,6)         
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['obj_after_EEE'          ][:]    = np.ma.round(iobj_after_EEE,6)          
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['para_before_EEE'        ][:]    = np.ma.round(ipara_before_EEE,6)        
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['para_after_EEE'         ][:]    = np.ma.round(ipara_after_EEE,6)         
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['para_active_EEE'        ][:]    = ipara_active_EEE        
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['calibrated_EEE'         ][:]    = icalibrated_EEE         
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['category_EEE'           ][:]    = icategory_EEE          
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['objective_EEE'          ][:]    = np.ma.round(iobjective_EEE,6)          
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['objective_open_loop_EEE'][:]    = np.ma.round(iobjective_open_loop_EEE,6)
            opt_results_EEE.groups[irecipe].groups[ibasin].variables['simtime_EEE'            ][:]    = np.array([(isimtime_EEE_nomask[ii]-refdate).days for ii in np.arange(nfiles)])           

            print("      shape: simtime:             ",np.shape(isimtime_EEE))
            print("      shape: objective:           ",np.shape(iobjective_EEE))
            print("      shape: objective_open_loop: ",np.shape(iobjective_open_loop_EEE))

    # close file
    opt_results_EEE.close()    # close results file

    print('Wrote optimization results for EEE-based recipes to: ',netcdf_file_EEE)
else:
    # read from netcdf
    opt_results_EEE  = nc4.Dataset(netcdf_file_EEE,  "r")

    for iirecipe, irecipe in enumerate(recipe_EEE):

        print('--------------------------------------')
        print('Read recipe: ',irecipe)

        for iibasin, ibasin in enumerate(basin_name):

            iobj_before_EEE           = opt_results_EEE.groups[irecipe].groups[ibasin].variables['obj_before_EEE'         ][:]
            iobj_after_EEE            = opt_results_EEE.groups[irecipe].groups[ibasin].variables['obj_after_EEE'          ][:]
            ipara_before_EEE          = opt_results_EEE.groups[irecipe].groups[ibasin].variables['para_before_EEE'        ][:]
            ipara_after_EEE           = opt_results_EEE.groups[irecipe].groups[ibasin].variables['para_after_EEE'         ][:]
            ipara_active_EEE          = opt_results_EEE.groups[irecipe].groups[ibasin].variables['para_active_EEE'        ][:]
            icalibrated_EEE           = opt_results_EEE.groups[irecipe].groups[ibasin].variables['calibrated_EEE'         ][:]
            icategory_EEE             = opt_results_EEE.groups[irecipe].groups[ibasin].variables['category_EEE'           ][:]
            iobjective_EEE            = opt_results_EEE.groups[irecipe].groups[ibasin].variables['objective_EEE'          ][:]
            iobjective_open_loop_EEE  = opt_results_EEE.groups[irecipe].groups[ibasin].variables['objective_open_loop_EEE'][:]
            isimtime_EEE              = opt_results_EEE.groups[irecipe].groups[ibasin].variables['simtime_EEE'            ][:]

            isimtime_EEE = [ refdate+datetime.timedelta(days=np.int(ii)) if ii != -9999 else -9999 for ii in isimtime_EEE ]    
            isimtime_EEE = np.ma.array(isimtime_EEE, mask=icategory_EEE.mask)

            icalibrated_EEE  = np.array(icalibrated_EEE,dtype=bool)
            ipara_active_EEE = np.array(ipara_active_EEE,dtype=bool)

            obj_before_EEE[iirecipe][iibasin]          = np.ma.round(iobj_before_EEE,6)
            obj_after_EEE[iirecipe][iibasin]           = np.ma.round(iobj_after_EEE,6)   
            para_before_EEE[iirecipe][iibasin]         = np.ma.round(ipara_before_EEE,6) 
            para_after_EEE[iirecipe][iibasin]          = np.ma.round(ipara_after_EEE,6)  
            para_active_EEE[iirecipe][iibasin]         = ipara_active_EEE 
            category_EEE[iirecipe][iibasin]            = icategory_EEE
            objective_EEE[iirecipe][iibasin]           = np.ma.round(iobjective_EEE,6)
            objective_open_loop_EEE[iirecipe][iibasin] = np.ma.round(iobjective_open_loop_EEE,6)
            simtime_EEE[iirecipe][iibasin]             = isimtime_EEE

            # some summary stats printed to screen
            print('   Basin: ',ibasin)
            print('   ---------------------------')
            print('   summer period (Jul-Nov)')
            print('   ---------------------------')
            tt_idx=((np.array([dd.month if dd != -9999 else -9999 for dd in isimtime_EEE.data]) >= 7) & (np.array([dd.month if dd != -9999 else 9999 for dd in isimtime_EEE.data ]) <= 11))
            print('      number days                  calibrated: ',sum(icalibrated_EEE[tt_idx]),   '   --> initial Q(t) was more than 5% off')
            print('      number days     successfully calibrated: ', np.sum(iobjective_EEE[(icalibrated_EEE&tt_idx),2]<0.05),' (',100.0*np.sum(iobjective_EEE[(icalibrated_EEE&tt_idx),2]<0.05)/np.sum((icalibrated_EEE&tt_idx)),'%)   --> final Q(t) is less than 5% off')
            print('      number days with optimal initial states: ', np.sum(iobjective_EEE[tt_idx,2]<0.05),' (',100.0*np.sum(iobjective_EEE[tt_idx,2]<0.05)/len(icalibrated_EEE[tt_idx]),'%)   --> final Q(t) is less than 5% off')
            print('      number days              not calibrated: ', sum(~icalibrated_EEE[tt_idx]), '   --> initial Q(t) was less than 5% off')
            print("      shape: simtime:             ",np.shape(isimtime_EEE))
            print("      shape: objective:           ",np.shape(iobjective_EEE))
            print("      shape: objective_open_loop: ",np.shape(iobjective_open_loop_EEE))

    # close file
    opt_results_EEE.close()    # close results file

    objective_info = collections.OrderedDict()
    objective_info['waeT-7_0']      = ['abs. of weighted summed errors of days t-7 ... t', 1]
    objective_info['absT0']         = ['absolute error of first time step', 1]
    objective_info['relT0']         = ['relative absolute error of first time step', 1]
    objective_info['VE_T1_3']       = ['volume errors of days t+1 ... t+3', 1]
    objective_info['VE_T1_7']       = ['volume errors of days t+1 ... t+7', 1]
    objective_info['VE_14orWinter'] = ['volume errors for 14 days (Jul-Nov) or\nfor remaining winter (Dec-Jun)', 1]
    objective_info['VE_W01']        = ['volume errors of first  1 week(s)', 1]
    objective_info['VE_W02']        = ['volume errors of first  2 week(s)', 2]
    objective_info['VE_W04']        = ['volume errors of first  4 week(s)', 4]
    objective_info['VE_W08']        = ['volume errors of first  8 week(s)', 8]
    objective_info['VE_W12']        = ['volume errors of first 12 week(s)', 12]
    objective_info['VE_W16']        = ['volume errors of first 16 week(s)', 16]
    objective_info['VE_W20']        = ['volume errors of first 20 week(s)', 20]
    objective_info['VE_W26']        = ['volume errors of first 26 week(s)', 26]
    objective_info['VE_W52']        = ['volume errors of first 52 week(s)', 52]

ncategory_EEE = np.max([ np.max( [ np.max(category_EEE[ii][jj]) for jj in range(nbasins) ]) for ii in range(len(category_EEE)) ])

print('-----------------------------------------------------------------')
print('READ: categorize_meteo_decision_tree_expert')

input_folder_expert = "/Users/j6mai/Documents/GitHub/CRD-DA/scripts/engage_plus_runs/categorize_meteo_decision_tree_expert"
recipe = '[recipe_1, recipe_2, recipe_3, recipe_5]'
recipe = '[recipe_5, recipe_3, recipe_1]'

# convert string inputs into lists
recipe_expert = recipe.split('[')[1].split(']')[0].split(',')
recipe_expert = [ ii.strip() for ii in recipe_expert ]
nrecipe_expert = len(recipe_expert)
    
nbasins = len(basin_name)



# overwrite read_from_netcdf if file does not exist
if not(os.path.exists(netcdf_file_expert)):
    read_from_netcdf = False

obj_before_expert          = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
obj_after_expert           = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
para_before_expert         = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
para_after_expert          = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
para_active_expert         = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
category_expert            = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
simtime_expert             = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
objective_expert           = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   
objective_open_loop_expert = [ [ [] for ibasin in range(nbasins) ] for irecipe in range(nrecipe_expert) ]   

para_info_expert   = categorizer_expert(return_all_para_info=True)

recipe_info_tmp_expert  = categorizer_expert(return_all_recipe_info=True)
recipe_info_expert = {}
for ir,rr in enumerate(recipe_info_tmp_expert[0]):
    recipe_info_expert[rr] = recipe_info_tmp_expert[1][ir]

if not(read_from_netcdf):

    # write to netcdf
    opt_results_expert  = nc4.Dataset(netcdf_file_expert,  "w", format="NETCDF4")

    for iirecipe, irecipe in enumerate(recipe_expert):

        print('--------------------------------------')
        print('Read recipe: ',irecipe)

        grp_recipe = opt_results_expert.createGroup(irecipe)

        for iibasin, ibasin in enumerate(basin_name):

            # create group
            grp_basin = opt_results_expert.createGroup(irecipe+'/'+ibasin)

            # path = os.path.normpath(os.path.join(dir_path,'../engage_plus_runs',pre_string+categorizer.__name__,basins[ibasin][0]))
            path_forecast = os.path.normpath(os.path.join(input_folder_expert,irecipe,basins[ibasin][0]))
            path_openloop = os.path.normpath(os.path.join(openloop_folder,basins[ibasin][0]))

            # files and number of files
            adjust_token = np.sort(glob.glob(os.path.join(path_forecast,'*','adjustments.token')))
            nfiles = len(adjust_token)
            objectives_token = np.sort(glob.glob(os.path.join(path_forecast,'*','simulation','objectives.out')))
            if (nfiles > len(objectives_token)): # should be the same but sometimes objective was already written but calibration not finished and hence adjustements file is not existing
                raise ValueError('Number of adjustment file does not match number of files containing objective function values.')

            # save all the simulation start times T
            isimtime_expert = np.array([ datetime.datetime(np.int(aa.split('/')[-2].split('-')[0]),
                                                    np.int(aa.split('/')[-2].split('-')[1]),
                                                    np.int(aa.split('/')[-2].split('-')[2]),0,0) for aa in adjust_token ])

            # read in one file to see number of parameters
            f = open(adjust_token[0], 'r')
            content = f.readlines()
            f.close()
            idx_v = np.where([ cc.startswith('v') for cc in content ])[0]
            npara = np.int(np.shape(idx_v)[0]/2)

            # read one file to see number of objectives
            f = open(objectives_token[0], 'r')
            content = f.readlines()
            f.close()
            nobj = np.shape(content)[0]
            objective_info = collections.OrderedDict()
            for cc in content:
                ccc = cc.split(',')
                ccc[2] = ccc[2].rstrip()   # remove trailing newline \n
                if 'week(s' in ccc[2]:
                    # track number of weeks to normalize objectives per week
                    nweeks = np.int(ccc[2].strip().split('week(s')[0].strip().split(' ')[-1])
                else:
                    nweeks = 1

                # if string too long add newline
                if len(ccc[2]) > 50:
                    tmp = ccc[2].split()
                    split_idx = np.sum(np.cumsum([ len(tt) for tt in tmp ]) < 35)
                    ccc[2] = ' '.join(tmp[0:split_idx])+'\n'+' '.join(tmp[split_idx:])
                    
                objective_info[ccc[0]] = [ ccc[2], nweeks ]

            # allocate result vectors and matrices
            iobj_before_expert          = np.ones(nfiles) * -9999.0
            iobj_after_expert           = np.ones(nfiles) * -9999.0
            ipara_before_expert         = np.ones([nfiles,npara]) * -9999.0
            ipara_after_expert          = np.ones([nfiles,npara]) * -9999.0
            ipara_active_expert         = np.ones([nfiles,npara],dtype=bool)
            icalibrated_expert          = np.ones(nfiles,dtype=bool)
            icategory_expert            = np.ones(nfiles,dtype=int) * -9999
            iobjective_expert           = np.ones([nfiles,nobj],dtype=float) * -9999.0
            iobjective_open_loop_expert = np.ones([nfiles,nobj],dtype=float) * -9999.0

            # create dimensions
            dim_npara    = grp_basin.createDimension("npara" , npara)
            dim_ntime    = grp_basin.createDimension("ntime" , nfiles)
            dim_nobj     = grp_basin.createDimension("nobj"  , nobj)

            # create variables
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/obj_before_expert'          , "f4",("ntime"),         zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/obj_after_expert'           , "f4",("ntime"),         zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/para_before_expert'         , "f4",("ntime","npara"), zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/para_after_expert'          , "f4",("ntime","npara"), zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/para_active_expert'         , "i4",("ntime","npara"), zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/calibrated_expert'          , "i4",("ntime"),         zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/category_expert'            , "i4",("ntime"),         zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/objective_expert'           , "f4",("ntime","nobj"),  zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/objective_open_loop_expert' , "f4",("ntime","nobj"),  zlib=True)
            grp_var = opt_results_expert.createVariable(irecipe+'/'+ibasin+'/simtime_expert'             , "i4",("ntime"),         zlib=True)
            grp_var.setncattr("long_name",     "time")
            grp_var.setncattr("units",         "days since "+str(refdate))
            grp_var.setncattr("calendar",      "gregorian")
            grp_var.setncattr("standard_name", "time")

            for ia, aa in enumerate(adjust_token):

                f = open(aa, 'r')
                content_adjust = f.readlines()
                f.close()

                calib_token = glob.glob(os.path.join(os.path.dirname(aa),'calib.token'))
                if len(calib_token) == 0:
                    icalibrated_expert[ia] = False
                else:
                    icalibrated_expert[ia] = True
                    f = open(calib_token[0], 'r')
                    content_calib = f.readlines()
                    f.close()
                    dict_calib = {}
                    for cc in content_calib:
                        ccc = cc.strip().split(',')
                        dict_calib[ccc[0]] = ','.join(ccc[1:])
                    icategory_expert[ia] = np.int(dict_calib['category'])

                idx_obj_before_expert = np.where([ 'obj_before' in cc for cc in content_adjust ])[0][0]
                idx_obj_after_expert  = np.where([ 'obj_after'  in cc for cc in content_adjust ])[0][0]

                iobj_before_expert[ia] = np.float(content_adjust[idx_obj_before_expert].strip().split(',')[1])
                iobj_after_expert[ia]  = np.float(content_adjust[idx_obj_after_expert].strip().split(',')[1])

                for pp in range(npara):
                    idx_v = np.where([ cc.startswith('v'+str(pp+1)+',') for cc in content_adjust ])[0]

                    ipara_before_expert[ia,pp] = np.float(content_adjust[idx_v[0]].strip().split(',')[1])
                    ipara_after_expert[ia,pp]  = np.float(content_adjust[idx_v[1]].strip().split(',')[1])
                    ipara_active_expert[ia,pp] = (content_adjust[idx_v[0]].strip().split(',')[2]=='True')

            # performance measures of forecast
            for io, oo in enumerate(objectives_token):

                if io < nfiles: # make sure that calibration was finished and adjustements.token exists
                    f = open(oo, 'r')
                    content_objectives = f.readlines()
                    f.close()
                    
                    iobjective_expert[io,0:nobj] = np.array([ np.float(cc.split(',')[1]) for cc in content_objectives ])

            # performance measures of open-loop simulation
            for io, oo in enumerate(objectives_token):

                if io < nfiles: # make sure that calibration was finished and adjustements.token exists
                    current_day = oo.split('/')[-3]
                    oo_openloop = os.path.normpath(os.path.join(path_openloop,current_day,'objectives.out'))
                    f = open(oo_openloop, 'r')
                    content_objectives = f.readlines()
                    f.close()
                    
                    iobjective_open_loop_expert[io,0:nobj] = np.array([ np.float(cc.split(',')[1]) for cc in content_objectives ])

            # apply mask on runs that were not calibrated
            iobj_before_expert  = np.ma.array(iobj_before_expert,  mask=~icalibrated_expert)
            iobj_after_expert   = np.ma.array(iobj_after_expert,   mask=~icalibrated_expert)
            #
            # only parameters that are calibrated (and active) 
            ipara_before_expert = np.ma.array(ipara_before_expert, mask=~(np.transpose(np.array([ icalibrated_expert for ii in range(npara) ])) & ipara_active_expert))
            ipara_after_expert  = np.ma.array(ipara_after_expert,  mask=~(np.transpose(np.array([ icalibrated_expert for ii in range(npara) ])) & ipara_active_expert))
            #
            # only mask out runs where no calibration happened
            ipara_active_expert = np.ma.array(ipara_active_expert, mask=~(np.transpose(np.array([ icalibrated_expert for ii in range(npara) ]))))
            #
            icategory_expert    = np.ma.array(icategory_expert,    mask=~icalibrated_expert)
            isimtime_expert_nomask = copy.deepcopy(isimtime_expert)
            isimtime_expert     = np.ma.array(isimtime_expert,     mask=~icalibrated_expert)

            # some summary stats printed to screen
            print('   Basin: ',ibasin)
            print('   ---------------------------')
            print('   summer period (Jul-Nov)')
            print('   ---------------------------')
            tt_idx=((np.array([dd.month if dd != -9999 else -9999 for dd in isimtime_expert.data ]) >= 7) & (np.array([dd.month if dd != -9999 else 9999 for dd in isimtime_expert.data ]) <= 11))
            print('      number days                  calibrated: ',sum(icalibrated_expert[tt_idx]),   '   --> initial Q(t) was more than 5% off')
            print('      number days     successfully calibrated: ', np.sum(iobjective_expert[(icalibrated_expert&tt_idx),2]<0.05),' (',100.0*np.sum(iobjective_expert[(icalibrated_expert&tt_idx),2]<0.05)/np.sum((icalibrated_expert&tt_idx)),'%)   --> final Q(t) is less than 5% off')
            print('      number days with optimal initial states: ', np.sum(iobjective_expert[tt_idx,2]<0.05),' (',100.0*np.sum(iobjective_expert[tt_idx,2]<0.05)/len(icalibrated_expert[tt_idx]),'%)   --> final Q(t) is less than 5% off')
            print('      number days              not calibrated: ', sum(~icalibrated_expert[tt_idx]), '   --> initial Q(t) was less than 5% off')

            obj_before_expert[iirecipe][iibasin]          = np.ma.round(iobj_before_expert,6)
            obj_after_expert[iirecipe][iibasin]           = np.ma.round(iobj_after_expert,6)   
            para_before_expert[iirecipe][iibasin]         = np.ma.round(ipara_before_expert,6) 
            para_after_expert[iirecipe][iibasin]          = np.ma.round(ipara_after_expert,6)  
            para_active_expert[iirecipe][iibasin]         = ipara_active_expert 
            category_expert[iirecipe][iibasin]            = icategory_expert
            simtime_expert[iirecipe][iibasin]             = isimtime_expert
            objective_expert[iirecipe][iibasin]           = np.ma.round(iobjective_expert,6)
            objective_open_loop_expert[iirecipe][iibasin] = np.ma.round(iobjective_open_loop_expert,6)

            # write float variables
            opt_results_expert.groups[irecipe].groups[ibasin].variables['obj_before_expert'         ][:]    = np.ma.round(iobj_before_expert,6)         
            opt_results_expert.groups[irecipe].groups[ibasin].variables['obj_after_expert'          ][:]    = np.ma.round(iobj_after_expert,6)          
            opt_results_expert.groups[irecipe].groups[ibasin].variables['para_before_expert'        ][:]    = np.ma.round(ipara_before_expert,6)        
            opt_results_expert.groups[irecipe].groups[ibasin].variables['para_after_expert'         ][:]    = np.ma.round(ipara_after_expert,6)         
            opt_results_expert.groups[irecipe].groups[ibasin].variables['para_active_expert'        ][:]    = ipara_active_expert        
            opt_results_expert.groups[irecipe].groups[ibasin].variables['calibrated_expert'         ][:]    = icalibrated_expert         
            opt_results_expert.groups[irecipe].groups[ibasin].variables['category_expert'           ][:]    = icategory_expert          
            opt_results_expert.groups[irecipe].groups[ibasin].variables['objective_expert'          ][:]    = np.ma.round(iobjective_expert,6)          
            opt_results_expert.groups[irecipe].groups[ibasin].variables['objective_open_loop_expert'][:]    = np.ma.round(iobjective_open_loop_expert,6) 
            opt_results_expert.groups[irecipe].groups[ibasin].variables['simtime_expert'            ][:]    = np.array([(isimtime_expert_nomask[ii]-refdate).days for ii in np.arange(nfiles)])    

            print("      shape: simtime:             ",np.shape(isimtime_expert))
            print("      shape: objective:           ",np.shape(iobjective_expert))
            print("      shape: objective_open_loop: ",np.shape(iobjective_open_loop_expert))

    # close file
    opt_results_expert.close()    # close results file

    print('Wrote optimization results for expert-based recipes to: ',netcdf_file_expert)
else:
    # read from netcdf
    opt_results_expert  = nc4.Dataset(netcdf_file_expert,  "r")

    for iirecipe, irecipe in enumerate(recipe_expert):

        print('--------------------------------------')
        print('Read recipe: ',irecipe)

        for iibasin, ibasin in enumerate(basin_name):

            iobj_before_expert           = opt_results_expert.groups[irecipe].groups[ibasin].variables['obj_before_expert'         ][:]
            iobj_after_expert            = opt_results_expert.groups[irecipe].groups[ibasin].variables['obj_after_expert'          ][:]
            ipara_before_expert          = opt_results_expert.groups[irecipe].groups[ibasin].variables['para_before_expert'        ][:]
            ipara_after_expert           = opt_results_expert.groups[irecipe].groups[ibasin].variables['para_after_expert'         ][:]
            ipara_active_expert          = opt_results_expert.groups[irecipe].groups[ibasin].variables['para_active_expert'        ][:]
            icalibrated_expert           = opt_results_expert.groups[irecipe].groups[ibasin].variables['calibrated_expert'         ][:]
            icategory_expert             = opt_results_expert.groups[irecipe].groups[ibasin].variables['category_expert'           ][:]
            iobjective_expert            = opt_results_expert.groups[irecipe].groups[ibasin].variables['objective_expert'          ][:]
            iobjective_open_loop_expert  = opt_results_expert.groups[irecipe].groups[ibasin].variables['objective_open_loop_expert'][:]
            isimtime_expert              = opt_results_expert.groups[irecipe].groups[ibasin].variables['simtime_expert'            ][:]

            isimtime_expert = [ refdate+datetime.timedelta(days=np.int(ii)) if ii != -9999 else -9999 for ii in isimtime_expert ]    
            isimtime_expert = np.ma.array(isimtime_expert, mask=icategory_expert.mask)

            icalibrated_expert  = np.array(icalibrated_expert,dtype=bool)
            ipara_active_expert = np.array(ipara_active_expert,dtype=bool)

            obj_before_expert[iirecipe][iibasin]          = np.ma.round(iobj_before_expert,6)
            obj_after_expert[iirecipe][iibasin]           = np.ma.round(iobj_after_expert,6)   
            para_before_expert[iirecipe][iibasin]         = np.ma.round(ipara_before_expert,6) 
            para_after_expert[iirecipe][iibasin]          = np.ma.round(ipara_after_expert,6)  
            para_active_expert[iirecipe][iibasin]         = ipara_active_expert 
            category_expert[iirecipe][iibasin]            = icategory_expert
            objective_expert[iirecipe][iibasin]           = np.ma.round(iobjective_expert,6)
            objective_open_loop_expert[iirecipe][iibasin] = np.ma.round(iobjective_open_loop_expert,6)
            simtime_expert[iirecipe][iibasin]             = isimtime_expert

            # some summary stats printed to screen
            print('   Basin: ',ibasin)
            print('   ---------------------------')
            print('   summer period (Jul-Nov)')
            print('   ---------------------------')
            tt_idx=((np.array([dd.month if dd != -9999 else -9999 for dd in isimtime_expert.data ]) >= 7) & (np.array([dd.month if dd != -9999 else 9999 for dd in isimtime_expert.data ]) <= 11))
            print('      number days                  calibrated: ',sum(icalibrated_expert[tt_idx]),   '   --> initial Q(t) was more than 5% off')
            print('      number days     successfully calibrated: ', np.sum(iobjective_expert[(icalibrated_expert&tt_idx),2]<0.05),' (',100.0*np.sum(iobjective_expert[(icalibrated_expert&tt_idx),2]<0.05)/np.sum((icalibrated_expert&tt_idx)),'%)   --> final Q(t) is less than 5% off')
            print('      number days with optimal initial states: ', np.sum(iobjective_expert[tt_idx,2]<0.05),' (',100.0*np.sum(iobjective_expert[tt_idx,2]<0.05)/len(icalibrated_expert[tt_idx]),'%)   --> final Q(t) is less than 5% off')
            print('      number days              not calibrated: ', sum(~icalibrated_expert[tt_idx]), '   --> initial Q(t) was less than 5% off')
            print("      shape: simtime:             ",np.shape(isimtime_expert))
            print("      shape: objective:           ",np.shape(iobjective_expert))
            print("      shape: objective_open_loop: ",np.shape(iobjective_open_loop_expert))

    # close file
    opt_results_expert.close()    # close results file

    objective_info = collections.OrderedDict()
    objective_info['waeT-7_0']      = ['abs. of weighted summed errors of days t-7 ... t', 1]
    objective_info['absT0']         = ['absolute error of first time step', 1]
    objective_info['relT0']         = ['relative absolute error of first time step', 1]
    objective_info['VE_T1_3']       = ['volume errors of days t+1 ... t+3', 1]
    objective_info['VE_T1_7']       = ['volume errors of days t+1 ... t+7', 1]
    objective_info['VE_14orWinter'] = ['volume errors for 14 days (Jul-Nov) or\nfor remaining winter (Dec-Jun)', 1]
    objective_info['VE_W01']        = ['volume errors of first  1 week(s)', 1]
    objective_info['VE_W02']        = ['volume errors of first  2 week(s)', 2]
    objective_info['VE_W04']        = ['volume errors of first  4 week(s)', 4]
    objective_info['VE_W08']        = ['volume errors of first  8 week(s)', 8]
    objective_info['VE_W12']        = ['volume errors of first 12 week(s)', 12]
    objective_info['VE_W16']        = ['volume errors of first 16 week(s)', 16]
    objective_info['VE_W20']        = ['volume errors of first 20 week(s)', 20]
    objective_info['VE_W26']        = ['volume errors of first 26 week(s)', 26]
    objective_info['VE_W52']        = ['volume errors of first 52 week(s)', 52]
        
# ncategory = np.max([ np.max(category[ii]) for ii in range(len(category)) ])
ncategory_expert = np.max([ np.max( [ np.max(category_expert[ii][jj]) for jj in range(nbasins) ]) for ii in range(len(category_expert)) ])









# setup plot environment

outtype = outtype.lower()
outtypes = ['', 'pdf', 'png', 'html', 'd3']
if outtype not in outtypes:
    print('\nError: output type must be in ', outtypes)
    import sys
    sys.exit()

t1 = time.time()

if (outtype == 'd3'):
    try:
        import mpld3
    except:
        print("No mpld3 found. Use output type html instead of d3.")
        outtype = 'html'

if dowhite:
    fgcolor = 'white'
    bgcolor = 'black'
else:
    fgcolor = 'black'
    bgcolor = 'white'


# ------------------------------------------
# (E) Plot 
# ------------------------------------------

# -------------------------------------------------------------------------
# Customize plots
# -------------------------------------------------------------------------

nrow        = 5           # # of rows of subplots per figure
ncol        = 2           # # of columns of subplots per figure
hspace      = 0.1        # x-space between subplots
vspace      = 0.015       # y-space between subplots
left        = 0.125       # right space on page
right       = 0.9         # right space on page
bottom      = 0.25        # right space on page
# bottom      = 0.55        # only 5 catchments
top         = 0.9         # right space on page
textsize    = 9          # standard text size
dxabc       = 0.985       # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 1.02        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth      = 1.0         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 1.0         # axis line width
msize       = 2.5         # marker size
mwidth      = 1.0         # marker edge width
# ufz blue (0.0, 0.34509803921568627, 0.611764705882353)
# ufz dark blue (0.0, 0.24313725490196078, 0.43137254901960786)
# ufz red (0.8313725490196079, 0.17647058823529413, 0.07058823529411765)
mcol1       = color.get_brewer('rdylbu4').colors[1]   # observation          --> orange
mcol2       = color.get_brewer('rdylbu4').colors[3]   # optimal simulation   --> light blue
mcol3       = color.get_brewer('rdylbu4').colors[2]   # initial simulation   --> dark blue
mcol        = color.colours(['blue','red','darkblue','gray','yellow', 'green', 'darkgreen', 'darkblue','darkgrey','black'])
lcol1       = color.colours('gray')
lcol2       = color.colours('black')
lcol3       = '0.0'

textbox_x  = 0.5
textbox_y  = 1.1

# Legend
llxbbox     = 1.0         # x-anchor legend bounding box
llybbox     = 0.95        # y-anchor legend bounding box
llrspace    = 0.2         # spacing between rows in legend
llcspace    = 0.05        # spacing between columns in legend
llhtextpad  = 0.4         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

# PNG
dpi         = 300
transparent = False
bbox_inches = 'tight'
pad_inches  = 0

import matplotlib as mpl
if (outtype == 'pdf'):
    mpl.use('PDF') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    # Customize: http://matplotlib.sourceforge.net/users/customizing.html
    mpl.rc('ps', papersize='a4', usedistiller='xpdf') # ps2pdf
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------
outtype_ends = ['', '.pdf', '_', '.html', '.html']
if plotname == '':
    plotfile = __file__[0:__file__.rfind(".")] + outtype_ends[outtypes.index(outtype)]
else:
    plotfile = plotname
if outtype == '':
    print('    Plot X')
else:
    print('    Plot ', plotfile)

if (outtype == 'pdf'):
    pdf_pages = PdfPages(plotfile)
# figsize = mpl.rcParams['figure.figsize']

if outtype in ['html','d3']:
    print('    Write html file ', plotfile)
    fhtml = open(plotfile,'w')
    print('<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">', file=fhtml)
    print("<html>", file=fhtml)
    print("<head>", file=fhtml)
    print("<title>"+__file__[0:__file__.rfind(".")]+"</title>", file=fhtml)
    print("</head>", file=fhtml)
    print("<body>", file=fhtml)


alphamod    = 0.7

ctextsize = 'xx-small'    # textsize of basins and parameters
cbartsize = textsize #'small'   # textsize at color bar (default)
htextsize = 'large'  # textsize Hidden/Standard box
ftextsize = 'large'  # textsize of flux labeling
    
# Brewer sequential
cols_para = color.get_brewer('YlOrRd9', rgb=True)   # need to be at least 9 colors
cols_para = color.get_brewer('Paired9', rgb=True)   # need to be at least 9 colors
cols_basin = color.get_brewer('YlOrRd4', rgb=True)   # need to be at least 4 colors
cols_basin = color.get_brewer('Paired4', rgb=True)   # need to be at least 4 colors
cols = color.get_brewer('RdBu8', rgb=True)
cols_category = color.get_brewer('RdYlBu6', rgb=True)    # need to be at least 6 colors
cols_dense = color.get_brewer('BlueWhiteOrangeRed', rgb=True) 
#cols = color.get_brewer('Spectral8', rgb=True)
# if not dolog:
#     # remove two colors starting at the end
#     cols.pop(2)
#     cols.pop(0)

# add same alpha as for categories
def add_alpha(col, alpha):
    icol = list(col)
    icol.append(alpha)
    return tuple(icol)
cols_category = [ add_alpha(c, alphamod) for c in cols_category ]
# first color is gray
if dowhite:
    colors = [(0.3,0.3,0.3,1.)]
else:
    colors = [(0.7,0.7,0.7,1.)]
colors.extend(cols_category)

# color map
cmap       = mpl.colors.ListedColormap(colors[0:max(ncategory_EEE,ncategory_expert)+1]) # all the categories needed plus gray
cmap_dense = mpl.colors.ListedColormap(cols_dense) 

def metric_name(objective_short_name):

    if (objective_short_name == 'waeT-7_0'):
        return 'WAE'
    elif (objective_short_name == 'absT0'):
        return '$M_0^{abs}$'
        #return '$|Q^{obs}(0)-Q^{sim}(0)|$'
    elif (objective_short_name == 'relT0'):
        return '$M_0^{rel}$'
        #return '$|Q^{obs}(0)-Q^{sim}(0)|/Q^{obs}(0)$'
    elif (objective_short_name == 'VE_T1_3'):
        return '$M_1$'
    elif (objective_short_name == 'VE_T1_7'):
        return '$M_2$'
    elif (objective_short_name == 'VE_14orWinter'):
        return '$M_3$'
    elif (objective_short_name == 'VE_W01'):
        return '$M_4^{(1)}$'
    elif (objective_short_name == 'VE_W02'):
        return '$M_4^{(2)}$'
    elif (objective_short_name == 'VE_W04'):
        return '$M_4^{(4)}$'
    elif (objective_short_name == 'VE_W08'):
        return '$M_4^{(8)}$'
    elif (objective_short_name == 'VE_W12'):
        return '$M_4^{(12)}$'
    elif (objective_short_name == 'VE_W16'):
        return '$M_4^{(16)}$'
    elif (objective_short_name == 'VE_W20'):
        return '$M_4^{(20)}$'
    elif (objective_short_name == 'VE_W26'):
        return '$M_4^{(26)}$'
    elif (objective_short_name == 'VE_W52'):
        return '$M_4^{(52)}$'
    else:
        return '???'




# ---------------------------
# boxplots of improvements of objective function per basin and EEE_recipe_1, expert_recipe_1, expert_recipe_5
# ---------------------------

ifig = 0

ifig += 1
iplot = 0
fig = plt.figure(ifig)


# -----------------------------------------
# nominal improvements
# -----------------------------------------
iplot += 1

sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

# ylabel
ylabel = 'Nominal Improvement [$m^3 s^{-1}$]'
plt.setp(sub, ylabel=ylabel)

for iirecipe_EEE, irecipe_EEE in enumerate(recipe_EEE):

    # EEE_recipes
    for iibasin, ibasin in enumerate(basin_name):

        dxx = 1./(nrecipe_expert+nrecipe_EEE+1)

        if only_summer:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) >= 7) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) <= 11) )
        else:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) >= 1) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) <= 12) )

        color_method = colors[iirecipe_EEE+1]
        sub.boxplot((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx] ),
                        showfliers=False,
                        positions=[iibasin-dxx*(nrecipe_EEE-iirecipe_EEE)+dxx/2.],
                        whis=[10,90],
                        widths=dxx*0.8,
                        patch_artist=True,
                        boxprops=dict(facecolor=color_method, color=color_method),
                        capprops=dict(color=color_method),
                        whiskerprops=dict(color=color_method),
                        flierprops=dict(color=color_method, markeredgecolor=color_method),
                        medianprops=dict(color=color_method),
                        )
        median = np.median((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx]))
        p25 = np.percentile((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx]),25)
        p75 = np.percentile((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx]),75)
        sub.text(iibasin-dxx*(nrecipe_EEE-iirecipe_EEE)+dxx/2.,  (p25+p75)/2., str2tex(astr(median,prec=1),usetex=usetex), #transform=sub.transAxes,
                            rotation=90, fontsize=textsize-4,
                            horizontalalignment='center', verticalalignment='center')

for iirecipe_expert, irecipe_expert in enumerate(recipe_expert):

    # expert_recipes
    for iibasin, ibasin in enumerate(basin_name):

        dxx = 1./(nrecipe_expert+nrecipe_EEE+1)

        if only_summer:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) >= 7) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) <= 11) )
        else:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) >= 1) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) <= 12) )

        color_method = colors[nrecipe_EEE+iirecipe_expert+1]
        sub.boxplot((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin]) - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])),
                        showfliers=False,
                        positions=[iibasin+dxx*(iirecipe_expert+1)-dxx/2.],
                        whis=[10,90],
                        widths=dxx*0.8,
                        patch_artist=True,
                        boxprops=dict(facecolor=color_method, color=color_method),
                        capprops=dict(color=color_method),
                        whiskerprops=dict(color=color_method),
                        flierprops=dict(color=color_method, markeredgecolor=color_method),
                        medianprops=dict(color=color_method),
                        )
        median = np.median((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx]))
        p25 = np.percentile((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx]),25)
        p75 = np.percentile((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx]),75)
        sub.text(iibasin+dxx*(iirecipe_expert+1)-dxx/2., (p25+p75)/2., str2tex(astr(median,prec=1),usetex=usetex), #transform=sub.transAxes,
                            rotation=90, fontsize=textsize-4,
                            horizontalalignment='center', verticalalignment='center')

# Get artists and labels for legend and chose which ones to display
handles, labels = sub.get_legend_handles_labels()
display         = [0]

# Create custom artists
boxes = [plt.Line2D((0,2),(0,0), color=colors[ii+1], marker='s', linestyle='') for ii in range(nrecipe_EEE+nrecipe_expert)]
box_labels = [str2tex('EEE '+rename_recipes_EEE[ii].replace('_',' '),usetex=usetex) for ii in recipe_EEE] + [str2tex('Expert '+rename_recipes_expert[ii].replace('_',' '),usetex=usetex) for ii in recipe_expert]

# Create legend from custom artist/label lists
sub.legend([handle for i,handle in enumerate(handles) if i in display]+boxes, 
           [label  for i,label  in enumerate(labels)  if i in display]+box_labels, 
           frameon=frameon, ncol=2,
           fontsize=textsize-2,
           labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
           loc='upper center', bbox_to_anchor=(1.15,-0.15), scatterpoints=1, numpoints=1)

xlabel = ''
xticks = sub.get_xticks()
sub.set_xticks(range(len(basin_name)))
xticks = sub.get_xticks()
xticklabels = [ str2tex(basins[basin_name[ii]][0], usetex=usetex) for ii,itick in enumerate(xticks) ]
plt.setp(sub, xlabel=xlabel, xticklabels=xticklabels)

# set range of plot
plt.setp(sub, xlim=[-0.5,nbasins-0.4], ylim=[-20,180])

abc2plot(sub,dxabc,dyabc,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',horizontalalignment='right', verticalalignment='bottom',zorder=400)


# -----------------------------------------
# relative improvements
# -----------------------------------------
iplot += 1

sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

# ylabel
ylabel = 'Relative Improvement [%]'
plt.setp(sub, ylabel=ylabel)

for iirecipe_EEE, irecipe_EEE in enumerate(recipe_EEE):

    # EEE_recipes
    for iibasin, ibasin in enumerate(basin_name):

        dxx = 1./(nrecipe_expert+nrecipe_EEE+1)

        if only_summer:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) >= 7) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) <= 11) )
        else:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) >= 1) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_EEE[iirecipe_EEE][iibasin]) ]) <= 12) )

        color_method = colors[iirecipe_EEE+1]
        sub.boxplot((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx])/np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx]*100.,
                        showfliers=False,
                        positions=[iibasin-dxx*(nrecipe_EEE-iirecipe_EEE)+dxx/2.],
                        whis=[10,90],
                        widths=dxx*0.8,
                        patch_artist=True,
                        boxprops=dict(facecolor=color_method, color=color_method),
                        capprops=dict(color=color_method),
                        whiskerprops=dict(color=color_method),
                        flierprops=dict(color=color_method, markeredgecolor=color_method),
                        medianprops=dict(color=color_method),
                        )
        median = np.median((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx])/np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx]*100.)
        p25 = np.percentile((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx])/np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx]*100.,25)
        p75 = np.percentile((np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx] - np.ma.compressed(obj_after_EEE[iirecipe_EEE][iibasin])[tt_idx])/np.ma.compressed(obj_before_EEE[iirecipe_EEE][iibasin])[tt_idx]*100.,75)
        sub.text(iibasin-dxx*(nrecipe_EEE-iirecipe_EEE)+dxx/2., (p25+p75)/2., str2tex(astr(median,prec=1),usetex=usetex), #transform=sub.transAxes,
                            rotation=90, fontsize=textsize-4,
                            horizontalalignment='center', verticalalignment='center')

for iirecipe_expert, irecipe_expert in enumerate(recipe_expert):

    # expert_recipes
    for iibasin, ibasin in enumerate(basin_name):

        dxx = 1./(nrecipe_expert+nrecipe_EEE+1)

        if only_summer:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) >= 7) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) <= 11) )
        else:
            tt_idx=np.where( (np.array([dd.month if dd != -9999 else -9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) >= 1) &
                             (np.array([dd.month if dd != -9999 else  9999 for dd in np.ma.compressed(simtime_expert[iirecipe_expert][iibasin]) ]) <= 12) )

        color_method = colors[nrecipe_EEE+iirecipe_expert+1]
        sub.boxplot((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx])/np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx]*100.,
                        showfliers=False,
                        positions=[iibasin+dxx*(iirecipe_expert+1)-dxx/2.],
                        whis=[10,90],
                        widths=dxx*0.8,
                        patch_artist=True,
                        boxprops=dict(facecolor=color_method, color=color_method),
                        capprops=dict(color=color_method),
                        whiskerprops=dict(color=color_method),
                        flierprops=dict(color=color_method, markeredgecolor=color_method),
                        medianprops=dict(color=color_method),
                        )
        median = np.median((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx])/np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx]*100.)
        p25 = np.percentile((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx])/np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx]*100.,25)
        p75 = np.percentile((np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx] - np.ma.compressed(obj_after_expert[iirecipe_expert][iibasin])[tt_idx])/np.ma.compressed(obj_before_expert[iirecipe_expert][iibasin])[tt_idx]*100.,75)
        sub.text(iibasin+dxx*(iirecipe_expert+1)-dxx/2., (p25+p75)/2., str2tex(astr(median,prec=1),usetex=usetex), #transform=sub.transAxes,
                            rotation=90, fontsize=textsize-4,
                            horizontalalignment='center', verticalalignment='center')

# # Get artists and labels for legend and chose which ones to display
# handles, labels = sub.get_legend_handles_labels()
# display         = [0]

# # Create custom artists
# boxes = [plt.Line2D((0,2),(0,0), color=colors[ii+1], marker='s', linestyle='') for ii in range(nrecipe_EEE+nrecipe_expert)]
# box_labels = [str2tex('EEE '+ii.replace('_',' '),usetex=usetex) for ii in recipe_EEE] + [str2tex('Expert '+ii.replace('_',' '),usetex=usetex) for ii in recipe_expert]

# # Create legend from custom artist/label lists
# sub.legend([handle for i,handle in enumerate(handles) if i in display]+boxes, 
#            [label  for i,label  in enumerate(labels)  if i in display]+box_labels, 
#            frameon=frameon, ncol=4,
#            fontsize=textsize-3,
#            labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
#            loc='upper center', bbox_to_anchor=(1.15,-0.3), scatterpoints=1, numpoints=1)

xlabel = ''
xticks = sub.get_xticks()
sub.set_xticks(range(len(basin_name)))
xticks = sub.get_xticks()
xticklabels = [ str2tex(basins[basin_name[ii]][0], usetex=usetex) for ii,itick in enumerate(xticks) ]
plt.setp(sub, xlabel=xlabel, xticklabels=xticklabels)

# set range of plot
plt.setp(sub, xlim=[-0.5,nbasins-0.4], ylim=[-2,110])

abc2plot(sub,dxabc,dyabc,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',horizontalalignment='right', verticalalignment='bottom',zorder=400)

# # -----------------------------------------
# # count on how often one method is better than another
# # -----------------------------------------
# for iibasin, ibasin in enumerate(basin_name):

    # iplot += 1

    # sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

    # line1 = sub.plot( obj_after_EEE[0][iibasin] )
    # color_method = colors[1]
    # plt.setp(line1, linestyle='-', color=color_method, linewidth=lwidth,
    #              marker='None', markeredgecolor=color_method, markerfacecolor=color_method,
    #              markersize=msize, markeredgewidth=mwidth)
    
    # line2 = sub.plot( obj_after_expert[0][iibasin] )
    # color_method = colors[2]
    # plt.setp(line2, linestyle='-', color=color_method, linewidth=lwidth,
    #              marker='None', markeredgecolor=color_method, markerfacecolor=color_method,
    #              markersize=msize, markeredgewidth=mwidth)
    
    # line3 = sub.plot( obj_after_expert[1][iibasin] )
    # color_method = colors[3]
    # plt.setp(line3, linestyle='-', color=color_method, linewidth=lwidth,
    #              marker='None', markeredgecolor=color_method, markerfacecolor=color_method,
    #              markersize=msize, markeredgewidth=mwidth)

    # print('---------------------------------------------------------------------')
    # print('BASIN: ',ibasin)
    # #                        EEE   expert      expert 
    # #    > 1% Better than          recipe 1    recipe 2
    # # EEE                     -    XX %        XX %
    # # expert recipe 1        XX %    -         XX %
    # # expert recipe 2        XX %  XX %         -
    # perform_matrix = np.array([[-9.9,
    #                             np.ma.sum(obj_after_EEE[0][iibasin].data    < 0.99*obj_after_expert[0][iibasin].data)*100.0/np.ma.count(obj_after_EEE[0][iibasin].data   ),
    #                             np.ma.sum(obj_after_EEE[0][iibasin].data    < 0.99*obj_after_expert[1][iibasin].data)*100.0/np.ma.count(obj_after_EEE[0][iibasin].data   ) ],
    #                            [np.ma.sum(obj_after_expert[0][iibasin].data < 0.99*obj_after_EEE[0][iibasin].data   )*100.0/np.ma.count(obj_after_expert[0][iibasin].data),
    #                             -9.9,
    #                             np.ma.sum(obj_after_expert[0][iibasin].data < 0.99*obj_after_expert[1][iibasin].data)*100.0/np.ma.count(obj_after_expert[0][iibasin].data) ],
    #                            [np.ma.sum(obj_after_expert[1][iibasin].data < 0.99*obj_after_EEE[0][iibasin].data   )*100.0/np.ma.count(obj_after_expert[1][iibasin].data),
    #                             np.ma.sum(obj_after_expert[1][iibasin].data < 0.99*obj_after_expert[0][iibasin].data)*100.0/np.ma.count(obj_after_expert[1][iibasin].data),
    #                             -9.9]])
    # print('EEE better than Expert recipe 1:             ',perform_matrix[0,1])
    # print('EEE better than Expert recipe 2:             ',perform_matrix[0,2])
    # print('Expert recipe 1 better than EEE:             ',perform_matrix[1,0])
    # print('Expert recipe 1 better than Expert recipe 2: ',perform_matrix[1,2])
    # print('Expert recipe 2 better than EEE:             ',perform_matrix[2,0])
    # print('Expert recipe 2 better than Expert recipe 1: ',perform_matrix[2,1])

    # print(' \\begin{tabular}{llrrr}                                                                              ')
    # print('     \hline                                                                                          ')
    # print('        & \multirow{2}{* }{Better than}     & EEE      & Expert   & Expert   \\\\                      ')
    # print('        &                                   & recipe   & recipe 1 & recipe 2 \\\\ \hline               ')
    # print('     \multirow{3}{*}{'+basins[ibasin][0]+'}  & EEE recipe      & -        & $\mathbf{'+astr(perform_matrix[0,1],prec=1)+'\%}$ & $'+astr(perform_matrix[0,2],prec=1)+'\%$ \\\\             ')
    # print('                          & Expert recipe 1 & $'+astr(perform_matrix[1,0],prec=1)+'\%$ & -        & $'+astr(perform_matrix[1,2],prec=1)+'\%$ \\\\                      ')
    # print('                          & Expert recipe 2 & $\mathbf{'+astr(perform_matrix[2,0],prec=1)+'\%}$ & $\mathbf{'+astr(perform_matrix[2,1],prec=1)+'\%}$ & -        \\ \hline    ')
    # print(' \end{tabular}   
    #                                                                                   ')

    
    


if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# ------------------------------
# Finish plot
# ------------------------------
if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()





