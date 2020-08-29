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
#       Water Resources Research, 56, e2020WR027960.
#       https://doi.org/10.1029/2020WR027960
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
# python figure_5.py -i ../data/eee/ -o ../data/eee/eee_results.nc -p figure_5.pdf -t pdf -u

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

parameters = ["v_1: level of soil reservoir [mm]",
              "v_2: level of groundwater reservoir [mm]",
              "v_3: depth of snowpack [mm]",
              "v_4: energy contained in snowpack [deg C]",
              "not analysed: variable removed from code (was indexTempNeige)",
              "v_5: precipitation [mm] (intercept a)",
              "not analysed: precipitation [mm] (slope b=1)",
              "v_6: temperature [deg C] (intercept a)",
              "not analysed: temperature [deg C] (slope b=1)"]

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
para_unused = np.arange(npara)[ np.array([ np.all(np.array(parameters_inf_ninf[0])[:,ipara] == -1) for ipara in range(np.shape(parameters_inf_ninf[0])[1]) ])][::-1]
para_used   = np.arange(npara)[~np.array([ np.all(np.array(parameters_inf_ninf[0])[:,ipara] == -1) for ipara in range(np.shape(parameters_inf_ninf[0])[1]) ])][::-1]

# grep results only for parameters that are used
parameters_inf_ninf    = np.array(parameters_inf_ninf)[:,:,para_used]
parameters_eee_vals    = np.array(parameters_eee_vals)[:,:,para_used]
parameters_eee_counter = np.array(parameters_eee_counter)[:,:,para_used]


# ------------------------------------------
# setup plot environment
# ------------------------------------------
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

nrow        = 9           # # of rows of subplots per figure
ncol        = 1           # # of columns of subplots per figure
hspace      = 0.02        # x-space between subplots
vspace      = 0.015       # y-space between subplots
left        = 0.125       # right space on page
right       = 0.9         # right space on page
bottom      = 0.25        # right space on page
# bottom      = 0.55        # only 5 catchments
top         = 0.9         # right space on page
textsize    = 10          # standard text size
dxabc       = 0.99       # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
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
dpi         = 600
transparent = False
bbox_inches = 'tight'
pad_inches  = 0

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42   # use Type 42 (a.k.a. TrueType) fonts for PostScript and PDF files
mpl.rcParams['ps.fonttype'] = 42    # use Type 42 (a.k.a. TrueType) fonts for PostScript and PDF files

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
    mpl.rc('text.latex') #, unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        #mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex') #, unicode=True)
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
# informative/non-informative parameters over time for each basin
# ----------------------------------------------------------------
ifig += 1
iplot = 0
fig = plt.figure(ifig)
    
for iibasin, ibasin in enumerate(basin_short_name):

    iplot += 1

    sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace))

    cvals = [-1,0,1]
    cvals = [0,1]
    cmap = mpl.colors.ListedColormap([color.get_brewer('RdYlBu7', rgb=True)[0],color.get_brewer('RdYlBu7', rgb=True)[2]])
    norm = mpl.colors.BoundaryNorm(cvals, cmap.N)
    mesh = sub.pcolormesh(np.transpose(parameters_inf_ninf[iibasin]), cmap=cmap, vmin=0, vmax=1, linewidths=0.0, edgecolor='None') 

    if iplot != nbasins:
        xlabel = ''
    else:
        xlabel = str2tex('Simulation Period Start Day $T-7$',usetex=usetex)
        
    xticks = sub.get_xticks()
    if iplot != nbasins:
        xticklabels = [ "" for ii,itick in enumerate(xticks) ]
    else:
        xticklabels = [ str2tex(dates[min(len(dates)-1,int(itick))], usetex=usetex) for ii,itick in enumerate(xticks) ]
    plt.setp(sub, xlabel=xlabel, xticklabels=xticklabels)

    # ylabel
    ylabel = str2tex('Variables',usetex=usetex)
    npara = len(parameters_inf_ninf[iibasin][0])
    sub.set_yticks([ ii+0.5 for ii in range(len(para_used)) ])
    yticks = sub.get_yticks()
    yticklabels = [ str2tex('$'+rename_vars['v_'+str(para_used[int(itick)]+1)]+'$', usetex=usetex) for ii,itick in enumerate(yticks) ]
    plt.setp(sub, ylabel=ylabel, yticklabels=yticklabels)

    # basin name
    sub.text(1.03, 0.5, str2tex('\n'.join(basin_long_name[iibasin].split()),usetex=usetex), fontsize=cbartsize, va='center', ha='center',rotation=90,transform=sub.transAxes)
    
    if (iplot == 1):
        # Colour bar
        pos  = position(nrow, ncol, iplot, hspace=hspace, vspace=vspace,
                            left=left, right=right, bottom=bottom, top=top)
        shrink = 0.6
        pos[0] += 0.5*(1.-shrink)*pos[2] # shrink
        pos[2]  = shrink*pos[2]
        pos[1] += 1.30*pos[3]            # shift up
        pos[3] *= 0.125                  # squeeze
        csub = fig.add_axes(pos, frameon=False)
        cbar = plt.colorbar(mesh, cax=csub, orientation='horizontal')
        
        cticks = cvals[:-1] # remove last tick
        cticks = cvals      # dont remove last tick
        cticknames = [ str2tex('not\ used',usetex=usetex), str2tex('informative',usetex=usetex), str2tex('noninformative',usetex=usetex)]
        cticknames = [ str2tex('informative',usetex=usetex), str2tex('noninformative',usetex=usetex)]
        cticknames = str2tex(cticknames, usetex=usetex)
        cbar.set_ticks([-2./3.,0.0,2./3.]) # draw ticks in the middle
        cbar.set_ticks([0.25,0.75]) # draw ticks in the middle
        cbar.set_ticklabels(cticknames)
        cbar.ax.tick_params(direction='out', length=3, width=1, colors=fgcolor, labelsize=cbartsize, pad=-1.0,
                                bottom=False, top=True, labelbottom=False, labeltop=True)
        clab = str2tex('Information Content\n of Model Variable', usetex=usetex)
        cbar.ax.yaxis.set_label_text(clab, fontsize=cbartsize, va='center', ha='right', rotation=0)
        cbar.ax.yaxis.set_label_coords(-0.03, 0.5, transform=csub.transAxes)

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


    

