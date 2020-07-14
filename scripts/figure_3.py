#!/usr/bin/env python
from __future__ import print_function

# Copyright 2019-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
#    This file is part of GitHub "EEE-DA" (https://github.com/julemai/EEE-DA) 
#    providing data and scripts to reproduce all figures of the publication:
#
#       J. Mai, R. Arsenault, B.A. Tolson, M. Latraverse, and K. Demeester (2020).
#       Application of Parameter Screening To Derive Optimal Initial State 
#       Adjustments for Streamflow Forecasting.
#       Water Resources Research, ??, ???-???.
#       https://doi.org/10.1002/???.
#
#    The EEE-DA codes are under MIT license.
#
#    Permission is hereby granted, free of charge, to any person obtaining a copy
#    of this software and associated documentation files (the "Software"), to deal
#    in the Software without restriction, including without limitation the rights
#    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#    copies of the Software, and to permit persons to whom the Software is
#    furnished to do so, subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.
#
# Usage
# python figure_3.py -b "[Peribonka, Montagnes Blanches]" -m ../data/open_loop/resultats_1953-01-01_original.nc -n ../data/open_loop/resultats_1953-01-01_corrected.nc -t 2012-01-01_2013-12-31 -i ../data/observations/ObservedData.nc -o figure_3.pdf

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path+'/lib')

# -----------------------
# Load additional packages
# -----------------------
import numpy as np
import glob  as glob
import collections
import netCDF4
import datetime
from matplotlib.dates import DateFormatter
import matplotlib.dates as mdates

import color                                                      # in lib/
import netcdf4     as nc4                                         # in lib/
from abc2plot      import abc2plot                                # in lib/
from fread         import fread                                   # in lib/
from position      import position                                # in lib/
from str2tex       import str2tex                                 # in lib/
from autostring    import astr                                    # in lib/               
from errormeasures import nse, kge, rmse, mse, mae, pear2, bias   # in lib/
from date2dec      import date2dec                                # in lib/
from dec2date      import dec2date                                # in lib/

import argparse
import textwrap                # nicer formatting of help text in parser
import numpy as np             # to perform numerics
import shutil                  # file operations

simulation_period      = ''
results_file_original  = ''
results_file_corrected = ''
inflow_file            = ''
output_file            = ''
basin_name             = None
only_summer            = False

parser       = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                  description='''Plots the simulated and observed inflows for a ceratin period of time.''')
parser.add_argument('-t', '--simulation_period', action='store',
                    default=simulation_period, dest='simulation_period', metavar='simulation_period', nargs=1,
                    help='Simulation period in format YYYY-MM-DD_YYYY-MM-DD.')
parser.add_argument('-m', '--results_file_original', action='store',
                    default=results_file_original, dest='results_file_original', metavar='results_file_original', nargs=1,
                    help='Name of CEQUEAU results file with original meteorology (with groups and variables).')
parser.add_argument('-n', '--results_file_corrected', action='store',
                    default=results_file_corrected, dest='results_file_corrected', metavar='results_file_corrected', nargs=1,
                    help='Name of CEQUEAU results file with corrected meteorology (with groups and variables).')
parser.add_argument('-b', '--basin_name', action='store',
                    default=basin_name, dest='basin_name', metavar='basin_name',
                    help='Name of the basin, e.g. ["Ashuapmushuan Amont","Montagnes Blanches"].')
parser.add_argument('-i', '--inflow_file', action='store',
                    default=inflow_file, dest='inflow_file', metavar='inflow_file', nargs=1,
                    help='Name of file containing observed inflows (must have same number of time steps as results).')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file',nargs=1,
                    help='Name of plot.')
parser.add_argument('-c', '--only_summer', action='store_true',
                    default=only_summer, dest="only_summer",
                    help='Only summer data (July 1 to Nov 30) are used. Default: False')


# example:
#     python figure_10.py -b "[Passes Dangereuses, Peribonka, Lac Manouane, Montagnes Blanches]" -m /Users/j6mai/Documents/GitHub/CRD-DA/scripts/forward_run/resultats_1953-01-01_original.nc -n /Users/j6mai/Documents/GitHub/CRD-DA/scripts/forward_run/resultats_1953-01-01_corrected.nc -t 1953-01-01_2017-04-10 -i /Users/j6mai/Documents/GitHub/CRD-DA/scripts/forward_run/Qobs.nc -o figure_10.pdf

args                       = parser.parse_args()
simulation_period          = args.simulation_period[0]
results_file_original      = args.results_file_original[0]
results_file_corrected     = args.results_file_corrected[0]
inflow_file                = args.inflow_file[0]
output_file                = args.output_file[0]
basin_name                 = args.basin_name
only_summer                = args.only_summer

del parser, args


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


    

start_day = simulation_period.split('_')[0]
start_day = datetime.datetime(np.int(start_day[0:4]),np.int(start_day[5:7]),np.int(start_day[8:10]),0,0)
end_day   = simulation_period.split('_')[1]
end_day   = datetime.datetime(np.int(end_day[0:4]),np.int(end_day[5:7]),np.int(end_day[8:10]),0,0)

ntime = (end_day-start_day).days+1


nse_original        = []
kge_original        = []
mse_original        = []
rmse_original       = []
mae_original        = []
bias_original       = []
pear2_original      = []
diff_T_0_original   = []
diff_T_1_original   = []
diff_T_2_original   = []

nse_corrected       = []
kge_corrected       = []
mse_corrected       = []
rmse_corrected      = []
mae_corrected       = []
bias_corrected      = []
pear2_corrected     = []
diff_T_0_corrected  = []
diff_T_1_corrected  = []
diff_T_2_corrected  = []

ttime_obs_all             = []
ttime_sim_all             = []
inflows_obs_original_all  = []
inflows_obs_corrected_all = []
inflows_sim_original_all  = []
inflows_sim_corrected_all = []
    

for basin in basin_name:

    basin_id = basins[basin][1]

    # ------------------------------------------
    # (A) Read simulated inflows
    # ------------------------------------------

    # open file
    nc_results_original  = nc4.NcDataset(results_file_original,  "r")
    nc_results_corrected = nc4.NcDataset(results_file_corrected, "r")

    # read inflows
    var    = "debitExutoire"
    ibasin = basin_id - 1
    group  = "etatsCP"
    ilag   = 0
    inflows_sim_original  = nc_results_original.groups[group].variables[var][:]          # [time,basin,lag]
    inflows_sim_original  = inflows_sim_original[:,ibasin,0]                             # lag = 0
    inflows_sim_corrected = nc_results_corrected.groups[group].variables[var][:]         # [time,basin,lag]
    inflows_sim_corrected = inflows_sim_corrected[:,ibasin,0]                            # lag = 0

    # read time
    ttime_original  = nc_results_original.variables[u't']
    ttime_original  = netCDF4.num2date(ttime_original[:],ttime_original.units)
    ttime_corrected = nc_results_corrected.variables[u't']
    ttime_corrected = netCDF4.num2date(ttime_corrected[:],ttime_corrected.units)
    

    if np.any(ttime_original != ttime_corrected):
        raise ValueError("Times in two simulation files is not the same. But needed for plot and getting Qobs. :(")
    else:
        ttime_sim = ttime_original

    # close file
    nc_results_original.close()    # close results file
    nc_results_corrected.close()   # close results file

    # cut the period from the data
    time_step_start       = list(ttime_sim).index(start_day)
    time_step_end         = list(ttime_sim).index(  end_day)

    ttime_sim             = ttime_sim[time_step_start:time_step_end+1]
    inflows_sim_original  = inflows_sim_original[time_step_start:time_step_end+1]
    inflows_sim_corrected = inflows_sim_corrected[time_step_start:time_step_end+1]

    # ------------------------------------------
    # (B) Read observed inflows
    # ------------------------------------------

    # open file
    nc_observ = nc4.NcDataset(inflow_file, "r")

    # read inflows
    ibasin = basin_id - 1
    inflows_obs = nc_observ.variables['Qobs'][:][ibasin,:]  # Qobs(idCP,t)

    # read time
    ttime_obs = nc_observ.variables[u't']
    ttime_obs = netCDF4.num2date(ttime_obs[:],ttime_obs.units)

    # close file
    nc_observ.close()   # close results file

    # cut the period from the data
    time_step_start = list(ttime_obs).index(start_day)
    time_step_end   = list(ttime_obs).index(  end_day)
    ttime_obs       = ttime_obs[time_step_start:time_step_end+1]
    inflows_obs     = inflows_obs[time_step_start:time_step_end+1]

    tt_idx_nonan_original  = np.where(~np.isnan(inflows_obs) & ~np.isnan(inflows_sim_original))[0]
    tt_idx_nonan_corrected = np.where(~np.isnan(inflows_obs) & ~np.isnan(inflows_sim_corrected))[0]

    inflows_obs_original   = np.ma.array( inflows_obs,           mask=(np.isnan(inflows_obs) | np.isnan(inflows_sim_original)) )
    inflows_sim_original   = np.ma.array( inflows_sim_original,  mask=(np.isnan(inflows_obs) | np.isnan(inflows_sim_original)) )
    inflows_obs_corrected  = np.ma.array( inflows_obs,           mask=(np.isnan(inflows_obs) | np.isnan(inflows_sim_corrected)) )
    inflows_sim_corrected  = np.ma.array( inflows_sim_corrected, mask=(np.isnan(inflows_obs) | np.isnan(inflows_sim_corrected)) )


    ttime_sim_all              += [ ttime_sim ]
    ttime_obs_all              += [ ttime_obs ]
    inflows_obs_original_all   += [ inflows_obs_original  ]
    inflows_sim_original_all   += [ inflows_sim_original  ]
    inflows_obs_corrected_all  += [ inflows_obs_corrected ]
    inflows_sim_corrected_all  += [ inflows_sim_corrected ]


    

    # ------------------------------------------
    # (C) Calculate objectives 
    # ------------------------------------------
    print("   ")
    print("   Basin: ",basin)
    print("   --------------------------------------------------------------------")
    print("   number of simulated inflows (original  meteo): ",np.shape(inflows_sim_original))
    print("   number of simulated inflows (corrected meteo): ",np.shape(inflows_sim_corrected))
    print("   number of observed  inflows (original  meteo): ",np.shape(inflows_obs_original))
    print("   number of observed  inflows (corrected meteo): ",np.shape(inflows_obs_corrected))
    print("   NSE (original meteo)                         = ",nse(inflows_obs_original[:],inflows_sim_original[:]))
    print("   NSE (corrected meteo)                        = ",nse(inflows_obs_corrected[:],inflows_sim_corrected[:]))
    print("   ")

    # weighted absolute error (WAE):
    #       weights differences of last time points higher than early time points:
    #         sum(  w_i * abs( Q_obs(i) - Q_sim(i) )  ) --> Min
    #         w_i = 2i / (N (N+1)) where N is the number of time points
    #                                    i = 1, ..., N is the index of the time point
    #         --> sum(w_i) = 1
    wi       = np.array([ (2.*ii)/((ntime+1.)*ntime) for ii in range(1,ntime+1) ])
    abs_diff_original = np.abs(inflows_obs[:]-inflows_sim_original[:])
    wae_original      = np.sum(wi*abs_diff_original)

    nse_original       += [ nse(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())    ]
    kge_original       += [ kge(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())    ]
    mse_original       += [ mse(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())    ]
    rmse_original      += [ rmse(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())   ]
    mae_original       += [ mae(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())    ]
    bias_original      += [ bias(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())   ]
    pear2_original     += [ pear2(inflows_obs_original[:].compressed(),inflows_sim_original[:].compressed())  ]
    diff_T_0_original  += [ np.abs(inflows_obs_original[-1]-inflows_sim_original[-1])                         ]
    diff_T_1_original  += [ np.abs(inflows_obs_original[-2]-inflows_sim_original[-2])                         ]
    diff_T_2_original  += [ np.abs(inflows_obs_original[-3]-inflows_sim_original[-3])                         ]

    # weighted absolute error (WAE):
    #       weights differences of last time points higher than early time points:
    #         sum(  w_i * abs( Q_obs(i) - Q_sim(i) )  ) --> Min
    #         w_i = 2i / (N (N+1)) where N is the number of time points
    #                                    i = 1, ..., N is the index of the time point
    #         --> sum(w_i) = 1
    wi       = np.array([ (2.*ii)/((ntime+1.)*ntime) for ii in range(1,ntime+1) ])
    abs_diff_corrected = np.abs(inflows_obs[:]-inflows_sim_corrected[:])
    wae_corrected      = np.sum(wi*abs_diff_corrected)

    nse_corrected       += [ nse(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())    ]
    kge_corrected       += [ kge(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())    ]
    mse_corrected       += [ mse(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())    ]
    rmse_corrected      += [ rmse(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())   ]
    mae_corrected       += [ mae(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())    ]
    bias_corrected      += [ bias(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())   ]
    pear2_corrected     += [ pear2(inflows_obs_corrected[:].compressed(),inflows_sim_corrected[:].compressed())  ]
    diff_T_0_corrected  += [ np.abs(inflows_obs_corrected[-1]-inflows_sim_corrected[-1])                         ]
    diff_T_1_corrected  += [ np.abs(inflows_obs_corrected[-2]-inflows_sim_corrected[-2])                         ]
    diff_T_2_corrected  += [ np.abs(inflows_obs_corrected[-3]-inflows_sim_corrected[-3])                         ]

    
# ------------------------------------------
# (D) Plot 
# ------------------------------------------

outtype  = 'pdf'
usetex   = False
serif    = False
dowhite  = False
alphamod = 0.7

# Brewer sequential
cols_para = color.get_brewer('YlOrRd9', rgb=True)   # need to be at least 9 colors
cols_para = color.get_brewer('Paired9', rgb=True)   # need to be at least 9 colors
cols_basin = color.get_brewer('YlOrRd4', rgb=True)   # need to be at least 4 colors
cols_basin = color.get_brewer('Paired4', rgb=True)   # need to be at least 4 colors
cols = color.get_brewer('RdBu8', rgb=True)
cols_category = color.get_brewer('RdYlBu8', rgb=True)    # need to be at least 8 colors
cols_dense = color.get_brewer('BlueWhiteOrangeRed', rgb=True) 

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

# -------------------------------------------------------------------------
# Customize plots
# -------------------------------------------------------------------------

nrow        = 5           # # of rows of subplots per figure
ncol        = 1           # # of columns of subplots per figure
hspace      = 0.10        # x-space between subplots
vspace      = 0.06       # y-space between subplots
textsize    = 12          # standard text size
dxabc       = 0.03        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.88        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth      = 1.0         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 1.0         # axis line width
msize       = 2.5         # marker size
mwidth      = 1.0         # marker edge width
mcol1       = color.colours('gray')      # primary marker colour
mcol2       = colors[2]                  # secondary
mcol3       = colors[8]                  # third
mcol        = color.colours(['red','blue','orange','gray','yellow', 'green', 'darkgreen', 'darkblue','darkgrey','black'])
lcol1       = color.colours('gray')
lcol2       = '0.0'
lcol3       = '0.0'

textbox_x  = 0.5
textbox_y  = 1.1

# Legend
llxbbox     = 1.0         # x-anchor legend bounding box
llybbox     = 1.0        # y-anchor legend bounding box
llrspace    = 0.4          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
llhtextpad  = 0.4         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

# PNG
dpi         = 300
transparent = False
bbox_inches = 'tight'
pad_inches  = 0

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
if (outtype == 'pdf'):
    print('Plot PDF ', output_file)
    pdf_pages = PdfPages(output_file)
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']
ifig = 0

ifig += 1
iplot = 0
fig = plt.figure(ifig)


for ibasin,basin in enumerate(basin_name):
    
    iplot += 1

    # discharge
    xlab   = '' #str2tex('Date',usetex=usetex)
    ylab   = str2tex('Discharge $Q$ [m$^3$ s$^{-1}$]',usetex=usetex)
    sub    = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))
      
    # plot observed inflows
    mark1 = sub.plot(ttime_obs_all[ibasin], inflows_obs_corrected_all[ibasin])
    plt.setp(mark1, linestyle='None', color=mcol1, linewidth=lwidth*0.0,
             alpha=0.5,
             marker='o', markeredgecolor=mcol1, markerfacecolor='None',
             markersize=msize, markeredgewidth=mwidth*0.5,
             label=str2tex('$Q_\mathrm{obs}$',usetex=usetex))

    # plot simulated inflows (original meteo)
    mark1 = sub.plot(ttime_sim_all[ibasin], inflows_sim_original_all[ibasin])
    plt.setp(mark1, linestyle='-', color=mcol2, linewidth=lwidth*0.8,
             marker='None', markeredgecolor=mcol2, markerfacecolor='None',
             markersize=msize, markeredgewidth=mwidth,
             label=str2tex('$Q_\mathrm{sim}^\mathrm{original}$',usetex=usetex))

    # plot simulated inflows (corrected meteo)
    mark1 = sub.plot(ttime_sim_all[ibasin], inflows_sim_corrected_all[ibasin])
    plt.setp(mark1, linestyle='-', color=mcol3, linewidth=lwidth*0.8,
             marker='None', markeredgecolor=mcol3, markerfacecolor='None',
             markersize=msize, markeredgewidth=mwidth,
             label=str2tex('$Q_\mathrm{sim}^\mathrm{corrected}$',usetex=usetex))

    plt.setp(sub, xlabel=xlab) # axis labels
    plt.setp(sub, ylabel=ylab)
    sub.grid(False)

    sub.text(0.0, 1.14, str2tex("Original: "), transform=sub.transAxes,
                 color=mcol2,
                 rotation=0, fontsize='small',
                 horizontalalignment='left', verticalalignment='center')
    sub.text(0.17, 1.14, str2tex("NSE = "+astr(nse_original[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                 color=mcol2,
                 rotation=0, fontsize='small',
                 horizontalalignment='left', verticalalignment='center')
    sub.text(0.34, 1.14, str2tex("KGE = "+astr(kge_original[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                 color=mcol2,
                 rotation=0, fontsize='small',
                 horizontalalignment='left', verticalalignment='center')
    sub.text(0.515, 1.14, str2tex("RMSE = "+astr(rmse_original[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                 color=mcol2,
                 rotation=0, fontsize='small',
                 horizontalalignment='left', verticalalignment='center')
    sub.text(0.7075, 1.14, str2tex("BIAS = "+astr(bias_original[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                 color=mcol2,
                 rotation=0, fontsize='small',
                 horizontalalignment='left', verticalalignment='center')
    sub.text(0.9, 1.14, str2tex("$r^2$ = "+astr(pear2_original[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                 color=mcol2,
                 rotation=0, fontsize='small',
                 horizontalalignment='left', verticalalignment='center')


    sub.text(0.0, 1.05, str2tex("Corrected: "), transform=sub.transAxes,
                     color=mcol3,
                     rotation=0, fontsize='small',
                     horizontalalignment='left', verticalalignment='center')
    sub.text(0.17, 1.05, str2tex("NSE = "+astr(nse_corrected[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                     color=mcol3,
                     rotation=0, fontsize='small',
                     horizontalalignment='left', verticalalignment='center')
    sub.text(0.34, 1.05, str2tex("KGE = "+astr(kge_corrected[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                     color=mcol3,             
                     rotation=0, fontsize='small',
                     horizontalalignment='left', verticalalignment='center')
    sub.text(0.515, 1.05, str2tex("RMSE = "+astr(rmse_corrected[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                     color=mcol3,
                     rotation=0, fontsize='small',
                     horizontalalignment='left', verticalalignment='center')
    sub.text(0.7075, 1.05, str2tex("BIAS = "+astr(bias_corrected[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                     color=mcol3,
                     rotation=0, fontsize='small',
                     horizontalalignment='left', verticalalignment='center')
    sub.text(0.9, 1.05, str2tex("$r^2$ = "+astr(pear2_corrected[ibasin],prec=2),usetex=usetex), transform=sub.transAxes,
                     color=mcol3,
                     rotation=0, fontsize='small',
                     horizontalalignment='left', verticalalignment='center')


    # basin name
    sub.text(1.02, 0.5, str2tex(basin,usetex=usetex), transform=sub.transAxes,
                     color='black',
                     rotation=90, #fontsize='small',
                     horizontalalignment='left', verticalalignment='center')

    ll = sub.legend(frameon=frameon, ncol=1,
                    labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                    loc='upper right', bbox_to_anchor=(llxbbox,llybbox), scatterpoints=1, numpoints=1,
                    fontsize = 'small')

    # Define the date format
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    years_fmt = mdates.DateFormatter('%Y')

    sub.xaxis.set_major_locator(years)
    sub.xaxis.set_major_formatter(years_fmt)
    sub.xaxis.set_minor_locator(months)
    date_form = DateFormatter("%Y")
    sub.xaxis.set_major_formatter(date_form)

    abc2plot(sub,dxabc,dyabc,iplot,bold=True,usetex=usetex,mathrm=True, large=True, parenthesis='none',
                 horizontalalignment='left',
                 verticalalignment='top',zorder=400)

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


