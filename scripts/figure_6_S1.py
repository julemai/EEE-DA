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
#    python figure_6_S1.py  -i ../data/eee/ -o ../data/eee/eee_results.nc -p figure_6_S1_ -t png -u

import argparse

dolog            = False
dowhite          = False
input_folder     = ''
plotname         = ''
outtype          = ''
usetex           = False
serif            = False
netcdf_file      = 'eee_results.nc'
read_from_netcdf = True

parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Plots results of exeriment for Engage Plus (calibration results and adjusted optimal parameters).''')
parser.add_argument('-i', '--input_folder', action='store',
                    default=input_folder, dest='input_folder', metavar='input_folder', nargs=1,
                    help='Name of input folder containing all experiment results (expected to have subfolders named with basins short names, e.g. ["ASAM", "MBLANC"].')
parser.add_argument('-o', '--netcdf_file', action='store',
                    default=netcdf_file, dest='netcdf_file', metavar='netcdf_file', 
                    help='File where EEE results are either dumped to or read from (default: <input_folder>/eee_results.nc).')
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

args             = parser.parse_args()
input_folder     = args.input_folder[0]
plotname         = args.plotname
outtype          = args.outtype
serif            = args.serif
usetex           = args.usetex
netcdf_file      = input_folder+'/eee/eee_results.nc'
netcdf_file      = args.netcdf_file

basin_name = None

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
import glob
import collections
import time
from sklearn.cluster import DBSCAN
import color                                                      # in lib/
import netcdf4     as nc4                                         # in lib/
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
#       Short         Long Name            ID   Upstream ID    Grid cells
basins['PD']      = ['Passes Dangereuses', 1,   None,          [1,2,3,4,5]]
basins['PERIB']   = ['Peribonka',          2,   1,             [6,7]]
basins['LM']      = ['Lac Manouane',       3,   None,          [8,9,10]]
basins['MBLANC']  = ['Montagnes Blanches', 4,   2,             [11]]

if (basin_name == None):
    basin_short_name   = [ bb for bb in basins ]
    basin_long_name    = [ basins[bb][0] for bb in basins ]
    basin_folder_names = [ str(basins[bb][1])+'_'+str(bb) for bb in basins ]
nbasins = len(basin_short_name)

# overwrite read_from_netcdf if file does not exist
if not(os.path.exists(netcdf_file)):
    read_from_netcdf = False

if not(read_from_netcdf):
    # -----------------------------------
    # Read EEE results
    # -----------------------------------
    dates = [ [] for ibasin in range(nbasins) ]
    for ibasin,ibasin_folder_name in enumerate(basin_folder_names):
        ff = glob.glob(input_folder+'/'+ibasin_folder_name+'/*/parameters.dat.final')
        dates[ibasin] = [ iff.split('/')[-2] for iff in ff ]
        
    parameters_inf_ninf    = [ [ [] for idate in dates[ibasin] ] for ibasin in range(nbasins) ] 
    parameters_eee_vals    = [ [ [] for idate in dates[ibasin] ] for ibasin in range(nbasins) ]   
    parameters_eee_counter = [ [ [] for idate in dates[ibasin] ] for ibasin in range(nbasins) ]

    for ibasin in range(nbasins):
        for idate in range(len(dates[ibasin])):

            ff = open(input_folder+'/'+basin_folder_names[ibasin]+'/'+dates[ibasin][idate]+'/'+'parameters.dat.final', 'r')
            content = ff.readlines()
            ff.close()
            inf = []
            for ic in content:
                if not( ic.startswith('#') ):
                    inf.append(np.int(ic.strip().split(' ')[-1]))

            parameters_inf_ninf[ibasin][idate] = inf
            npara = len(inf)

            iter_folders = glob.glob(input_folder+'/'+ibasin_folder_name+'/'+dates[ibasin][idate]+'/iter*')
            eee_count = np.array([ 0.0 for ii in range(npara) ])
            eee_vals  = np.array([ 0.0 for ii in range(npara) ])
            for iter_folder in iter_folders:
                ff = open(iter_folder+'/'+'eee_results.dat', 'r')
                content = ff.readlines()
                ff.close()

                ipar = 0
                for ic in content:
                    if not( ic.startswith('#') ):
                        ic_count = np.float(ic.strip().split()[3])
                        ic_val   = np.float(ic.strip().split()[2])

                        eee_count[ipar] += ic_count
                        eee_vals[ipar]  += ic_val*ic_count
                        ipar += 1

            eee_vals2 = np.where(eee_count>0.0, eee_vals/eee_count, 0.0)
            parameters_eee_vals[ibasin][idate]    = eee_vals2
            parameters_eee_counter[ibasin][idate] = eee_count

    ntime = np.shape(parameters_inf_ninf)[1]
    
    # write to netcdf
    eee_results  = nc4.NcDataset(netcdf_file,  "w", format="NETCDF4")

    # create dimensions
    dim_npara    = eee_results.createDimension("npara"     , npara)
    dim_nbasins  = eee_results.createDimension("nbasins"   , nbasins)
    dim_ntime    = eee_results.createDimension("ntime"     , ntime)

    # create variables
    nc4var = eee_results.createVariable("basin_short_name"         , str,  ("nbasins"),                   zlib=True)
    nc4var = eee_results.createVariable("basin_long_name"          , str,  ("nbasins"),                   zlib=True)
    nc4var = eee_results.createVariable("basin_folder_names"       , str,  ("nbasins"),                   zlib=True)
    nc4var = eee_results.createVariable("dates"                    , str,  ("ntime"),          zlib=True)
    nc4var = eee_results.createVariable("parameters"               , str,  ("npara"),                     zlib=True)
    nc4var = eee_results.createVariable("parameters_inf_ninf"      , int,  ("nbasins", "ntime", "npara"), zlib=True)
    nc4var = eee_results.createVariable("parameters_eee_vals"      , "f4", ("nbasins", "ntime", "npara"), zlib=True)
    nc4var = eee_results.createVariable("parameters_eee_counter"   , "f4", ("nbasins", "ntime", "npara"), zlib=True)

    # set some attributes
    eee_results.variables["basin_short_name"].setncattr(      "long_name", "Basin short name")
    eee_results.variables["basin_long_name"].setncattr(       "long_name", "Basin long name")
    eee_results.variables["basin_folder_names"].setncattr(    "long_name", "Name of folder containing results")
    eee_results.variables["dates"].setncattr(                 "long_name", "Start date in format YYYY-MM-DD of EEE simulation runs")
    eee_results.variables["parameters"].setncattr(            "long_name", "Description of parameters")
    eee_results.variables["parameters_inf_ninf"].setncattr(   "long_name", "(-1): parameter not analysed, (0): parameter informative, (1): parameter noninformative")
    eee_results.variables["parameters_eee_vals"].setncattr(   "long_name", "EEE sensitivity estimate")
    eee_results.variables["parameters_eee_counter"].setncattr("long_name", "Number of runs the EEE estimate is based on")

    # write info to variables; string variables need to be written sequentially
    for jj, val in enumerate(basin_short_name):    eee_results.variables["basin_short_name"][jj]   = val
    for jj, val in enumerate(basin_long_name):     eee_results.variables["basin_long_name"][jj]    = val
    for jj, val in enumerate(basin_folder_names):  eee_results.variables["basin_folder_names"][jj] = val
    for jj, val in enumerate(parameters):          eee_results.variables["parameters"][jj]         = val
    for ii in range(nbasins):
        for jj in range(ntime):
            val = dates[ii][jj]
            eee_results.variables["dates"][jj] = val

    

    # write float variables
    eee_results.variables['parameters_inf_ninf'][:]    = parameters_inf_ninf
    eee_results.variables['parameters_eee_vals'][:]    = parameters_eee_vals
    eee_results.variables['parameters_eee_counter'][:] = parameters_eee_counter

    eee_results.setncattr('references',  'Application of Parameter Screening To Derive Optimal Initial State Adjustments for Streamflow Forecasting (v1.0)')
    eee_results.setncattr('License',     'MIT licence')
    eee_results.setncattr('product',     'EEE-DA')
    eee_results.setncattr('title',       'Application of Parameter Screening To Derive Optimal Initial State Adjustments for Streamflow Forecasting')
    eee_results.setncattr('institution', 'University of Waterloo, Waterloo, ON, Canada')
    eee_results.setncattr('Conventions', 'CF-1.6')

    # close file
    eee_results.close()    # close results file

    print('Wrote EEE results to: ',netcdf_file)
else:
    # read from netcdf
    eee_results  = nc4.NcDataset(netcdf_file,  "r")
    
    npara                  = eee_results.dimensions['npara'].size               
    nbasins                = eee_results.dimensions['nbasins'].size                            
    basin_short_name       = eee_results.variables['basin_short_name'][:]      
    basin_long_name        = eee_results.variables['basin_long_name'][:]       
    basin_folder_names     = eee_results.variables['basin_folder_names'][:]    
    dates                  = eee_results.variables['dates'][:]                 
    parameters_inf_ninf    = eee_results.variables['parameters_inf_ninf'][:]   
    parameters_eee_vals    = eee_results.variables['parameters_eee_vals'][:]   
    parameters_eee_counter = eee_results.variables['parameters_eee_counter'][:]

    # close file
    eee_results.close()    # close results file

# parameter unused: True  --> Parameter never used (EE = -1)
#                   False --> Parameter at least once used in analysis (EE either 0 or 1)
para_unused = np.arange(npara)[ np.array([ np.all(np.array(parameters_inf_ninf[0])[:,ipara] == -1) for ipara in range(np.shape(parameters_inf_ninf[0])[1]) ])]
para_used   = np.arange(npara)[~np.array([ np.all(np.array(parameters_inf_ninf[0])[:,ipara] == -1) for ipara in range(np.shape(parameters_inf_ninf[0])[1]) ])]

# grep results only for parameters that are used
parameters_inf_ninf    = np.array(parameters_inf_ninf)[:,:,para_used]
parameters_eee_vals    = np.array(parameters_eee_vals)[:,:,para_used]
parameters_eee_counter = np.array(parameters_eee_counter)[:,:,para_used]
    
# -----------------------------------
# Read meteorological data (Kenjynized meteo)
# -----------------------------------
from categorize_meteo_decision_tree import get_temp_prec

input_file = "../data/observations/meteoCorrected.nc"
time_step  = -1
all_time_steps = True
date_time_step = ''

temp_basin = {}
prec_basin = {}
for ibasin in range(nbasins):
    
    period_time_step = dates[0]+'_'+dates[-1]
    grid_cells = str(basins[basin_short_name[ibasin]][3])
    # prec is 4 longer than period because we need precip of last 4 days
    [time_step_start,time_step_end,prec,temp] = get_temp_prec(input_file,time_step,all_time_steps,date_time_step,period_time_step,grid_cells)

    temp_basin[basin_short_name[ibasin]] = temp
    prec_basin[basin_short_name[ibasin]] = np.array([ np.sum(prec[ii:ii+4]) for ii in range(len(temp)) ])

def my_weighted_perc(data,perc,weights=None):
    if (weights is None):
        return np.nanpercentile(data,perc)
    else:
        d=data[(~np.isnan(data))&(~np.isnan(weights))]
        ix=np.argsort(d)
        d=d[ix]
        wei=weights[ix]
        wei_cum=100.*np.cumsum(wei*1./np.sum(wei))
        return np.interp(perc,wei_cum,d)
        
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

nrow        = npara       # # of rows of subplots per figure
ncol        = nbasins     # # of columns of subplots per figure
hspace      = 0.02        # x-space between subplots
vspace      = 0.015       # y-space between subplots
left        = 0.125       # right space on page
right       = 0.9         # right space on page
bottom      = 0.25        # right space on page
# bottom      = 0.55        # only 5 catchments
top         = 0.9         # right space on page
textsize    = 10          # standard text size
dxabc       = 0.99        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.90        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

lwidth      = 1.0         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 1.0         # axis line width
msize       = 2.5         # marker size
mwidth      = 1.0         # marker edge width
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
llrspace    = -0.02       # spacing between rows in legend
llcspace    = 0.05         # spacing between columns in legend
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
                
ifig = 0

# ----------------------------------------------------------------
# color-coded temp/precipitation according to parameters information content
# ----------------------------------------------------------------
ifig += 1
iplot = 0
fig = plt.figure(ifig)

col_ignor   = color.colours('lightgray')
col_inf     = color.get_brewer('RdYlBu7', rgb=True)[0] # dark red
col_ninf    = color.get_brewer('RdYlBu7', rgb=True)[2] # yellow
markers = ['o','o','o','o']

min_prec = np.floor(np.min([ np.percentile(prec_basin[ibasin],1) for ibasin in basin_short_name ])-10)
max_prec = np.ceil( np.max([ np.percentile(prec_basin[ibasin],99) for ibasin in basin_short_name ])+25)
min_temp = np.floor(np.min([ np.min(temp_basin[ibasin]) for ibasin in basin_short_name ]))
max_temp = np.ceil( np.max([ np.max(temp_basin[ibasin]) for ibasin in basin_short_name ]))

for ipara in range(len(para_used)):
    for iibasin, ibasin in enumerate(basin_short_name):

        iplot += 1

        sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

        # plot markers where parameter was ignored
        idx_ignor = np.where(np.array(parameters_inf_ninf[iibasin])[:,ipara]==-1)[0]
        idx_ninf  = np.where(np.array(parameters_inf_ninf[iibasin])[:,ipara]==1)[0]
        idx_inf   = np.where(np.array(parameters_inf_ninf[iibasin])[:,ipara]==0)[0]
        if len(idx_ignor) > 0:
            sub.scatter(prec_basin[ibasin][idx_ignor],temp_basin[ibasin][idx_ignor],marker=markers[iibasin],alpha=0.7,linewidths=0.0,s=msize,color=col_ignor)
        if len(idx_ninf) > 0:
            eee_vals_ninf = np.array(parameters_eee_vals[iibasin])[idx_ninf,ipara] 
            eee_vals_max_ninf = np.max(np.array(parameters_eee_vals[iibasin]))
            # eee_vals_max_ninf = np.max(eee_vals_ninf)           # <<<<<<<<<<<<<<<<<<<<< RETHINK!!! Should be max of all values in this basin not per parameter...
            eee_vals_med_ninf = np.median(eee_vals_ninf)
            eee_weights_ninf  = eee_vals_ninf/eee_vals_max_ninf
            # sub.scatter(prec_basin[ibasin][idx_ninf], temp_basin[ibasin][idx_ninf],
            #                 #markeredgecolor=col_ninf, #markerfacecolor=col_ninf, 
            #                 marker=markers[iibasin],alpha=0.7,linewidths=0.0,s=np.where(eee_weights_ninf>0.01,msize*12*eee_weights_ninf,0.01*12*msize),color=col_ninf)
            for ii in range(np.shape(prec_basin[ibasin][idx_ninf])[0]):
                sub.plot(prec_basin[ibasin][idx_ninf][ii],  temp_basin[ibasin][idx_ninf][ii],
                         linestyle='None',
                         marker='o',
                         markersize=np.sqrt(np.where(eee_weights_ninf>0.01,msize*12*eee_weights_ninf,0.01*12*msize)[ii]/np.pi),
                         markeredgewidth=0.1,
                         markeredgecolor=col_ninf, markerfacecolor=col_ninf, 
                         alpha=0.7,linewidth=0.0)
        if len(idx_inf) > 0:
            eee_vals_inf = np.array(parameters_eee_vals[iibasin])[idx_inf,ipara] 
            eee_vals_max_inf = np.max(np.array(parameters_eee_vals[iibasin]))
            # eee_vals_max_inf = np.max(eee_vals_inf)           # <<<<<<<<<<<<<<<<<<<<< RETHINK!!! Should be max of all values in this basin not per parameter...
            eee_vals_med_inf = np.median(eee_vals_inf)
            eee_weights_inf  = eee_vals_inf/eee_vals_max_inf
            # sub.scatter(prec_basin[ibasin][idx_inf],  temp_basin[ibasin][idx_inf],
            #                 #markeredgecolor=col_inf, #markerfacecolor=col_inf, 
            #                 marker=markers[iibasin],alpha=0.7,linewidths=0.0,s=np.where(eee_weights_inf>0.01,msize*12*eee_weights_inf,0.01*12*msize),color=col_inf)
            for ii in range(np.shape(prec_basin[ibasin][idx_inf])[0]):
                sub.plot(prec_basin[ibasin][idx_inf][ii],  temp_basin[ibasin][idx_inf][ii],
                         linestyle='None',
                         marker='o',
                         markersize=np.sqrt(np.where(eee_weights_inf>0.01,msize*12*eee_weights_inf,0.01*12*msize)[ii]/np.pi),
                         markeredgewidth=0.1,
                         markeredgecolor=col_inf, markerfacecolor=col_inf, 
                         alpha=0.7,linewidth=0.0)

        if len(idx_inf) > 0:

            # percentiles
            perc_inf_temp = [ np.percentile(temp_basin[ibasin][idx_inf],1), np.percentile(temp_basin[ibasin][idx_inf],99) ]
            # weighted percentiles
            perc_inf_temp = [ my_weighted_perc(temp_basin[ibasin][idx_inf],1 ,weights=eee_weights_inf),
                              my_weighted_perc(temp_basin[ibasin][idx_inf],99,weights=eee_weights_inf) ]
            perc_inf_prec = [ my_weighted_perc(prec_basin[ibasin][idx_inf],1 ,weights=eee_weights_inf),
                              my_weighted_perc(prec_basin[ibasin][idx_inf],99,weights=eee_weights_inf) ]

            sub.plot([min_prec,max_prec],[perc_inf_temp[0],perc_inf_temp[0]],color=col_inf, linewidth=0.3)
            sub.plot([min_prec,max_prec],[perc_inf_temp[1],perc_inf_temp[1]],color=col_inf, linewidth=0.3)
            sub.fill_between([min_prec,max_prec],[perc_inf_temp[0],perc_inf_temp[0]],[perc_inf_temp[1],perc_inf_temp[1]],color=col_inf,alpha=0.05)

            sub.plot([perc_inf_prec[0],perc_inf_prec[0]],[min_temp,max_temp],color=col_inf, linewidth=0.3)
            sub.plot([perc_inf_prec[1],perc_inf_prec[1]],[min_temp,max_temp],color=col_inf, linewidth=0.3)
            sub.fill_between([perc_inf_prec[0],perc_inf_prec[1]],[min_temp,min_temp],[max_temp,max_temp],color=col_inf,alpha=0.05)

            sub.text(max_prec-3,perc_inf_temp[0],str2tex('$p^T_{1}='+astr(perc_inf_temp[0],prec=1)+'$',usetex=usetex),fontsize='xx-small',va='top',ha='right')
            sub.text(max_prec-3,perc_inf_temp[1],str2tex('$p^T_{99}='+astr(perc_inf_temp[1],prec=1)+'$',usetex=usetex),fontsize='xx-small',va='bottom',ha='right')

            sub.text(perc_inf_prec[0],min_temp,str2tex('$p^P_{1}='+astr(perc_inf_prec[0],prec=1)+'$',usetex=usetex),rotation=90,fontsize='xx-small',va='bottom',ha='right')
            sub.text(perc_inf_prec[1],min_temp,str2tex('$p^P_{99}='+astr(perc_inf_prec[1],prec=1)+'$',usetex=usetex),rotation=90,fontsize='xx-small',va='bottom',ha='left')

            # average sensitivity of informative variables
            sub.text(0.5,1.01,str2tex('$\overline{EEE}_{inf}='+astr(eee_vals_med_inf,prec=1)+'$',usetex=usetex),color=col_inf, fontsize='xx-small',va='bottom',ha='center',transform=sub.transAxes)

        if iplot%ncol != 1:                # first col
            yticks = sub.get_yticks()
            yticklabels = [ "" for ii,itick in enumerate(yticks) ]
            plt.setp(sub, yticklabels=yticklabels)
        else:
            plt.setp(sub, ylabel=str2tex('Temp. [$^\circ$C]',usetex=usetex))

        # last col
        if iplot%ncol == 0:
            var_str = 'v'+str(para_used[ipara]+1)            # v1, v2, v3, v4, v6, v8
            var_str = rename_vars[var_str]                   # v1, v2, v3, v4, v5, v6
            var_str = '$'+var_str.replace('v','v_{')+'}$'    # $v_{1}$, $v_{2}$, $v_{3}$, $v_{4}$, $v_{5}$, $v_{6}$, 
            sub.text(1.06, 0.5, str2tex(var_str,usetex=usetex), fontsize=cbartsize, va='center', ha='left',rotation=90,transform=sub.transAxes)
        # first row
        if (iplot-1)//ncol == 0:
            sub.text(0.5, 1.15, str2tex(basin_long_name[iibasin],usetex=usetex), fontsize=cbartsize, va='bottom', ha='center',rotation=0,transform=sub.transAxes)

        # first plot
        if iplot == 1:
            # Get artists and labels for legend and chose which ones to display
            handles, labels = sub.get_legend_handles_labels()
            display         = [0]

            # Create custom artists
            box2legend_1 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,50)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_ninf, marker='o', markeredgecolor=col_ninf, markerfacecolor=col_ninf, linewidth=0.0, linestyle='')     # non-informative, small EEE
            box2legend_2 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,99)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_ninf, marker='o', markeredgecolor=col_ninf, markerfacecolor=col_ninf, linewidth=0.0, linestyle='')     # non-informative, large EEE
            box2legend_3 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,50)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_inf,  marker='o', markeredgecolor=col_inf,  markerfacecolor=col_inf,  linewidth=0.0, linestyle='')     #     informative, small EEE
            box2legend_4 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,99)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_inf,  marker='o', markeredgecolor=col_inf,  markerfacecolor=col_inf,  linewidth=0.0, linestyle='' )     #     informative, large EEE

            # Create legend from custom artist/label lists
            sub.legend([handle for i,handle in enumerate(handles) if i in display]+[box2legend_1,box2legend_2,box2legend_3,box2legend_4],
                       [label  for i,label  in enumerate(labels)  if i in display]+[str2tex('non-inf., $EEE = '+astr(np.percentile(eee_vals_inf,50), prec=1)+'$', usetex=usetex),
                                                                                    str2tex('non-inf., $EEE = '+astr(np.percentile(eee_vals_inf,99), prec=1)+'$', usetex=usetex),
                                                                                    str2tex(    'inf., $EEE = '+astr(np.percentile(eee_vals_inf,50), prec=1)+'$', usetex=usetex),
                                                                                    str2tex(    'inf., $EEE = '+astr(np.percentile(eee_vals_inf,99), prec=1)+'$', usetex=usetex)],
                       frameon=frameon, ncol=2,
                       fontsize=textsize, #-3,
                       labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                       loc='lower center', bbox_to_anchor=(2.0,1.3), scatterpoints=1, numpoints=1)
            

        if (iplot-1)//ncol != 5: #(nrow-1):    # last row
            xticks = sub.get_xticks()
            xticklabels = [ "" for ii,itick in enumerate(xticks) ]
            plt.setp(sub, xticklabels=xticklabels)
        else:
            plt.setp(sub, xlabel=str2tex('Accum. precip. over\n last 4 days [mm]',usetex=usetex))

        sub.set_xlim([min_prec,max_prec*1.1])
        sub.set_ylim([min_temp*1.1,max_temp*1.1])
    

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = plotname+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# ----------------------------------------------------------------
# color-coded temp/day-of-year according to parameters information content
# ----------------------------------------------------------------
ifig += 1
iplot = 0
fig = plt.figure(ifig)

col_ignor   = color.colours('lightgray')
col_inf     = color.get_brewer('RdYlBu7', rgb=True)[0] # dark red
col_ninf    = color.get_brewer('RdYlBu7', rgb=True)[2] # yellow
markers = ['o','o','o','o']

min_day  = 0
max_day  = 365
min_temp = np.floor(np.min([ np.min(temp_basin[ibasin]) for ibasin in basin_short_name ]))-13
max_temp = np.ceil( np.max([ np.max(temp_basin[ibasin]) for ibasin in basin_short_name ]))

for ipara in range(len(para_used)):
    for iibasin, ibasin in enumerate(basin_short_name):
        print("")

        iplot += 1

        sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

        doy = np.array([ ((
                      datetime.datetime(int(dates[ii][0:4]),int(dates[ii][5:7]),int(dates[ii][8:10]))-
                      datetime.datetime(int(dates[0][0:4]), int(dates[0][5:7]), int(dates[0][8:10]))).days+1)%365
                    for ii in range(len(dates)) ])
        doy = np.where(doy==0,365,doy)

        # plot markers where parameter was ignored
        idx_ignor = np.where(np.array(parameters_inf_ninf[iibasin])[:,ipara]==-1)[0]
        idx_ninf  = np.where(np.array(parameters_inf_ninf[iibasin])[:,ipara]==1)[0]
        idx_inf   = np.where(np.array(parameters_inf_ninf[iibasin])[:,ipara]==0)[0]

        print("basin: "+ibasin+"  para v"+str(para_used[ipara]+1)) 

        if len(idx_ignor) > 0:
            sub.scatter(doy[idx_ignor],temp_basin[ibasin][idx_ignor],marker=markers[iibasin],alpha=0.7,linewidths=0.0,s=msize,color=col_ignor)
        if len(idx_ninf) > 0:
            eee_vals_ninf = np.array(parameters_eee_vals[iibasin])[idx_ninf,ipara]
            eee_vals_max_ninf = np.max(np.array(parameters_eee_vals[iibasin]))
            eee_vals_med_ninf = np.median(eee_vals_ninf)
            eee_weights_ninf  = eee_vals_ninf/eee_vals_max_ninf
            for ii in range(np.shape(doy[idx_ninf])[0]):
                sub.plot(doy[idx_ninf][ii],  temp_basin[ibasin][idx_ninf][ii],
                         linestyle='None',
                         marker='o',
                         markersize=np.sqrt(np.where(eee_weights_ninf>0.01,msize*12*eee_weights_ninf,0.01*12*msize)[ii]/np.pi),
                         markeredgewidth=0.1,
                         markeredgecolor=col_ninf, markerfacecolor=col_ninf, 
                         alpha=0.7,linewidth=0.0)
        if len(idx_inf) > 0:
            eee_vals_inf = np.array(parameters_eee_vals[iibasin])[idx_inf,ipara]
            eee_vals_max_inf = np.max(np.array(parameters_eee_vals[iibasin]))
            eee_vals_med_inf = np.median(eee_vals_inf)
            eee_weights_inf  = eee_vals_inf/eee_vals_max_inf
            for ii in range(np.shape(doy[idx_inf])[0]):
                sub.plot(doy[idx_inf][ii],  temp_basin[ibasin][idx_inf][ii],
                         linestyle='None',
                         marker='o',
                         markersize=np.sqrt(np.where(eee_weights_inf>0.01,msize*12*eee_weights_inf,0.01*12*msize)[ii]/np.pi),
                         markeredgewidth=0.1,
                         markeredgecolor=col_inf, markerfacecolor=col_inf, 
                         alpha=0.7,linewidth=0.0)

        if len(idx_inf) > 0:
            # percentiles
            perc_inf_temp = [ np.percentile(temp_basin[ibasin][idx_inf],1), np.percentile(temp_basin[ibasin][idx_inf],99) ]
            # weighted percentiles
            perc_inf_temp = [ my_weighted_perc(temp_basin[ibasin][idx_inf],1 ,weights=eee_weights_inf),
                              my_weighted_perc(temp_basin[ibasin][idx_inf],99,weights=eee_weights_inf) ]
            perc_inf_prec = [ my_weighted_perc(prec_basin[ibasin][idx_inf],1 ,weights=eee_weights_inf),
                              my_weighted_perc(prec_basin[ibasin][idx_inf],99,weights=eee_weights_inf) ]
            perc_inf_doy  = [ my_weighted_perc(doy[idx_inf],1 ,weights=eee_weights_inf),
                              my_weighted_perc(doy[idx_inf],99,weights=eee_weights_inf) ]

            sub.plot([min_day,max_day],[perc_inf_temp[0],perc_inf_temp[0]],color=col_inf, linewidth=0.3)
            sub.plot([min_day,max_day],[perc_inf_temp[1],perc_inf_temp[1]],color=col_inf, linewidth=0.3)
            sub.fill_between([min_day,max_day],[perc_inf_temp[0],perc_inf_temp[0]],[perc_inf_temp[1],perc_inf_temp[1]],color=col_inf,alpha=0.05)

            sub.text(max_day-3,perc_inf_temp[0],str2tex('$p^T_{1}='+astr(perc_inf_temp[0],prec=1)+'$',usetex=usetex),fontsize='xx-small',va='top',ha='right')
            sub.text(max_day-3,perc_inf_temp[1],str2tex('$p^T_{99}='+astr(perc_inf_temp[1],prec=1)+'$',usetex=usetex),fontsize='xx-small',va='bottom',ha='right')

            if (para_used[ipara] == 2) or (para_used[ipara] == 3) or (para_used[ipara] == 5) or (para_used[ipara] == 7): 
                clustering=DBSCAN(eps=20, min_samples=2).fit(np.transpose(np.array([doy[idx_inf],temp_basin[ibasin][idx_inf]])),
                                                                 sample_weight=eee_weights_inf*2)
                cluster_labels = clustering.labels_
                iclus_plot = 0
                for iclus in np.unique(cluster_labels):
                    idx_clus=np.where(cluster_labels==iclus)
                    # print("  ------------------------------------------")
                    # print("  cluster label: ",iclus," DOYs: ",doy[idx_inf][idx_clus])
                    # print("  mean EEE of cluster: ",np.median(eee_weights_inf[idx_clus]))
                    # print("  max  EEE of cluster: ",np.max(eee_weights_inf[idx_clus]))
                    # print("  p95  EEE of cluster: ",np.percentile(eee_weights_inf[idx_clus],95))
                    
                    if np.percentile(eee_weights_inf[idx_clus],95) > 0.11: 
                        iclus_plot += 1
                        perc_inf_doy = [ my_weighted_perc(doy[idx_inf][idx_clus],1 ,weights=eee_weights_inf[idx_clus]),
                                         my_weighted_perc(doy[idx_inf][idx_clus],99 ,weights=eee_weights_inf[idx_clus]) ]
                        # print("      >>>>>>>>>> Plotted: DOY = ",perc_inf_doy)
                        sub.plot([perc_inf_doy[0],perc_inf_doy[0]],[min_temp,max_temp],color=col_inf, linewidth=0.3)
                        sub.plot([perc_inf_doy[1],perc_inf_doy[1]],[min_temp,max_temp],color=col_inf, linewidth=0.3)
                        sub.fill_between([perc_inf_doy[0],perc_inf_doy[1]],[min_temp,min_temp],[max_temp,max_temp],color=col_inf,alpha=0.1)
            
                        sub.text(perc_inf_doy[0],min_temp,str2tex('$p^{C_'+astr(iclus_plot)+'}_{1}='+astr(perc_inf_doy[0],prec=0)+'$',usetex=usetex),rotation=90,fontsize='xx-small',va='bottom',ha='right')
                        sub.text(perc_inf_doy[1],min_temp,str2tex('$p^{C_'+astr(iclus_plot)+'}_{99}='+astr(perc_inf_doy[1],prec=0)+'$',usetex=usetex),rotation=90,fontsize='xx-small',va='bottom',ha='left')
                
            # average sensitivity of informative variables
            sub.text(0.5,1.01,str2tex('$\overline{EEE}_{inf}='+astr(eee_vals_med_inf,prec=1)+'$',usetex=usetex),color=col_inf, fontsize='xx-small',va='bottom',ha='center',transform=sub.transAxes)

        if iplot%ncol != 1:                # first col
            yticks = sub.get_yticks()
            yticklabels = [ "" for ii,itick in enumerate(yticks) ]
            plt.setp(sub, yticklabels=yticklabels)
        else:
            plt.setp(sub, ylabel=str2tex('Temp. [$^\circ$C]',usetex=usetex))

        # last col
        if iplot%ncol == 0:
            var_str = 'v'+str(para_used[ipara]+1)            # v1, v2, v3, v4, v6, v8
            var_str = rename_vars[var_str]                   # v1, v2, v3, v4, v5, v6
            var_str = '$'+var_str.replace('v','v_{')+'}$'    # $v_{1}$, $v_{2}$, $v_{3}$, $v_{4}$, $v_{5}$, $v_{6}$, 
            sub.text(1.06, 0.5, str2tex(var_str,usetex=usetex), fontsize=cbartsize, va='center', ha='left',rotation=90,transform=sub.transAxes)
        # first row
        if (iplot-1)//ncol == 0:
            sub.text(0.5, 1.15, str2tex(basin_long_name[iibasin],usetex=usetex), fontsize=cbartsize, va='bottom', ha='center',rotation=0,transform=sub.transAxes)

        # first plot
        if iplot == 1:
            # Get artists and labels for legend and chose which ones to display
            handles, labels = sub.get_legend_handles_labels()
            display         = [0]

            # Create custom artists
            box2legend_1 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,50)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_ninf, marker='o', markeredgecolor=col_ninf, markerfacecolor=col_ninf, linewidth=0.0, linestyle='')     # non-informative, small EEE
            box2legend_2 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,99)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_ninf, marker='o', markeredgecolor=col_ninf, markerfacecolor=col_ninf, linewidth=0.0, linestyle='')     # non-informative, large EEE
            box2legend_3 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,50)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_inf,  marker='o', markeredgecolor=col_inf,  markerfacecolor=col_inf,  linewidth=0.0, linestyle='')     #     informative, small EEE
            box2legend_4 = plt.Line2D((0,2),(0,0),markersize=np.sqrt(msize*12*np.percentile(eee_weights_inf,99)/np.pi),
                                          #alpha=0.7,
                                          markeredgewidth=0.1,
                                          color=col_inf,  marker='o', markeredgecolor=col_inf,  markerfacecolor=col_inf,  linewidth=0.0, linestyle='' )     #     informative, large EEE

            # Create legend from custom artist/label lists
            sub.legend([handle for i,handle in enumerate(handles) if i in display]+[box2legend_1,box2legend_2,box2legend_3,box2legend_4],
                       [label  for i,label  in enumerate(labels)  if i in display]+[str2tex('non-inf., $EEE = '+astr(np.percentile(eee_vals_inf,50), prec=1)+'$', usetex=usetex),
                                                                                    str2tex('non-inf., $EEE = '+astr(np.percentile(eee_vals_inf,99), prec=1)+'$', usetex=usetex),
                                                                                    str2tex(    'inf., $EEE = '+astr(np.percentile(eee_vals_inf,50), prec=1)+'$', usetex=usetex),
                                                                                    str2tex(    'inf., $EEE = '+astr(np.percentile(eee_vals_inf,99), prec=1)+'$', usetex=usetex)],
                       frameon=frameon, ncol=2,
                       fontsize=textsize, #-3,
                       labelspacing=llrspace, handletextpad=llhtextpad, handlelength=llhlength,
                       loc='lower center', bbox_to_anchor=(2.0,1.3), scatterpoints=1, numpoints=1)
            

        if (iplot-1)//ncol != 5: #(nrow-1):    # last row
            xticks = sub.get_xticks()
            xticklabels = [ "" for ii,itick in enumerate(xticks) ]
            plt.setp(sub, xticklabels=xticklabels)
        else:
            plt.setp(sub, xlabel=str2tex('DOY',usetex=usetex))

        sub.set_xlim([min_day,max_day])
        sub.set_ylim([min_temp*1.1,max_temp*1.1])
    

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = plotname+"{0:04d}".format(ifig)+".png"
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


    

