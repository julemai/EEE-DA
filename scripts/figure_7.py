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
#    python figure_7.py  -i ../data/observations/meteo.nc -a -p figure_7.pdf -t pdf -c "[categorize_meteo_decision_tree_expert, categorize_meteo_decision_tree_EEE]" -r "[[recipe_5,recipe_3,recipe_1],[recipe_1,recipe_4,recipe_3]]"

import argparse

dolog            = False
dowhite          = False
verbose          = False
input_file       = ''
all_time_steps   = False
period_time_step = ''
grid_cells       = None
basin_name       = None
recipes          = None # categorization should be same for all recipes
plotname         = ''
outtype          = ''
usetex           = False
serif            = False
categorizers     = None

parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Extracts last time step from input file and writes it to output file. The variables and values are unchanged in the output file.''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', nargs=1,
                    help='Name of input file (standard is a CEQUEAU file with groups and variables).')
parser.add_argument('-a', '--all_time_steps', action='store_true',
                    default=all_time_steps, dest='all_time_steps', 
                    help='If -a used, all time steps will be categorized.')
parser.add_argument('-e', '--period_time_step', action='store',
                        default=period_time_step, dest='period_time_step', metavar='period_time_step',
                        help='Period of time steps which should be extracted (format: YYYY-MM-DD_YYYY-MM-DD).')
parser.add_argument('-g', '--grid_cells', action='store',
                    default=grid_cells, dest='grid_cells', metavar='grid_cells',
                    help='Grid cells that should be considered. List with comma as separator. E.g., [[20,21,22,23],[48]]. Numbering starts with 1.')
parser.add_argument('-b', '--basin_name', action='store',
                    default=basin_name, dest='basin_name', metavar='basin_name',
                    help='Name of the basin (should match grid cells given), e.g. ["Ashuapmushuan Amont","Montagnes Blanches"].')
parser.add_argument('-r', '--recipes', action='store',
                        default=recipes, dest='recipes', metavar='recipes',
                        help='Recipe to use for decision variables. Must be specified for each categorizer. E.g. "[[recipe_1,recipe_2]]"')
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
parser.add_argument('-c', '--categorizers', action='store',
                    default=categorizers, dest='categorizers', metavar='categorizers',
                    help='List of categorizers of dates according to their meteorological conditions. E.g. "[categorize_meteo_decision_tree_1]"')

args             = parser.parse_args()
input_file       = args.input_file[0]
grid_cells       = args.grid_cells
basin_name       = args.basin_name
recipes          = args.recipes
all_time_steps   = args.all_time_steps
period_time_step = args.period_time_step
plotname         = args.plotname
outtype          = args.outtype
serif            = args.serif
usetex           = args.usetex
categorizers     = args.categorizers

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
import color                                                      # in lib/
from position      import position                                # in lib/
from str2tex       import str2tex                                 # in lib/
from autostring    import astr                                    # in lib/
from abc2plot      import abc2plot                                # in lib/

if (input_file is None):
    print('\nError: Input file must be given.\n')
    import sys
    sys.exit()

outtype = outtype.lower()
outtypes = ['', 'pdf', 'png', 'html', 'd3']
if outtype not in outtypes:
    print('\nError: output type must be in ', outtypes)
    import sys
    sys.exit()

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

if (grid_cells == None) or (basin_name == None):
    basin_name = [ bb for bb in basins ]
    grid_cells = [ basins[bb][2] for bb in basins ]
else:
    # convert string inputs into lists
    basin_name = basin_name.split('[')[1].split(']')[0].split(',')
    strs = grid_cells.replace('[','').split('],')
    grid_cells = [map(int, s.replace(']','').split(',')) for s in strs]

nbasins = len(basin_name)

if recipes is None:
    raise ValueError("Recipe must be given, e.g. 'recipe_1'")

recipes = recipes.replace('[','').split('],')
recipes = [ s.replace(']','').split(',') for s in recipes]

if categorizers is None:
    raise ValueError("Name of routine used to categorize days based on their meteorology must be given, e.g. 'categorize_meteo_decision_tree_1'. Should be a function defined in 'categorize_meteo_decision_tree.py'.")

categorizers = categorizers.split('[')[1].split(']')[0].split(',')
categorizers = [ii.strip() for ii in categorizers]

if len(recipes) != len(categorizers):
    raise ValueError('Recipes need to be specified for each categorizer!')

all_cat_ttime       = [ [] for ii in categorizers ]
all_cat_categories  = [ [] for ii in categorizers ]
all_cat_recipe_info = [ [] for ii in categorizers ]
for icategorizer,categorizer in enumerate(categorizers):
    
    # which decision tree is used
    if categorizer == 'categorize_meteo_decision_tree_expert':
        from categorize_meteo_decision_tree import categorize_meteo_decision_tree_expert as categorizer
    elif categorizer == 'categorize_meteo_decision_tree_EEE':
        from categorize_meteo_decision_tree import categorize_meteo_decision_tree_EEE as categorizer
    else:
        raise ValueError('This categorization method is not implemented yet!')

    all_cat_recipe_info[icategorizer] = categorizer(input_file=input_file,return_all_recipe_info=True)

    all_ttime      = [ [] for ii in range(nbasins) ]
    all_categories = [ [] for ii in range(nbasins) ]
    for ibasin in range(nbasins):
        # get categories
        [ categories,ttime ] = categorizer(input_file=input_file,
                                           period_time_step=period_time_step,
                                           all_time_steps=all_time_steps,
                                           grid_cells=str(grid_cells[ibasin]),
                                           return_time=True,
                                           recipe=recipes[icategorizer][0])   # just a random recipe --> categorizing of days independent of that

        # append NODATA at beginning if first day is not Jan 1
        if ( not(ttime[0].month == 1 and ttime[0].day == 1) ):

            n_missing_days = (ttime[0] - datetime.datetime(ttime[0].year, 1, 1, 0, 0)).days
            missing_days   = np.array([ datetime.datetime(ttime[0].year, 1, 1, 0, 0) + datetime.timedelta(days=ii) for ii in range(n_missing_days) ])
            missing_data   = np.array([ -9999                                                                      for ii in range(n_missing_days) ])
            ttime          = np.append(missing_days, ttime)
            categories     = np.append(missing_data, categories)
            
        # append NODATA at end       if last  day is not Dec 31
        if ( not(ttime[-1].month == 12 and ttime[-1].day == 31) ):

            n_missing_days = (datetime.datetime(ttime[-1].year, 12, 31, 0, 0) - ttime[-1]).days
            missing_days   = np.array([ ttime[-1] + datetime.timedelta(days=ii+1) for ii in range(n_missing_days) ])
            missing_data   = np.array([ -9999                                     for ii in range(n_missing_days) ])
            ttime          = np.append(ttime,      missing_days)
            categories     = np.append(categories, missing_data)


        # remove leap days
        idx_leap   = np.where( (np.array([ tt.month for tt in ttime ]) != 2) | (np.array([ tt.day for tt in ttime ]) != 29) )[0]
        ttime      = ttime[idx_leap]
        categories = categories[idx_leap]

        if (np.shape(ttime)[0] % 365):
            raise ValueError('It seems that not every year has 365 days?!')

        if (np.shape(ttime)[0]/365 != np.shape(categories)[0]/365):
            raise ValueError('Time and categories have different number of data.')
        else:
            nyears = np.int(np.shape(categories)[0]/365)

        ttime      = np.reshape(ttime,     [nyears,365])
        categories = np.reshape(categories,[nyears,365])

        all_ttime[ibasin] = ttime
        all_categories[ibasin] = categories

    all_cat_ttime[icategorizer]      = all_ttime
    all_cat_categories[icategorizer] = all_categories


import numpy as np
from collections import OrderedDict
import time
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

nrow        = 8           # # of rows of subplots per figure
ncol        = 1           # # of columns of subplots per figure
hspace      = 0.10        # x-space between subplots
vspace      = 0.010       # y-space between subplots
left        = 0.125       # right space on page
right       = 0.9         # right space on page
bottom      = 0.25        # right space on page
# bottom      = 0.55        # only 5 catchments
top         = 0.9         # right space on page
textsize    = 10          # standard text size
dxabc       = 0.98        # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
dyabc       = 0.85        # % of (max-min) shift up from lower x-axis for a,b,c,... labels

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
lcol2       = '0.0'
lcol3       = '0.0'

textbox_x  = 0.5
textbox_y  = 1.1

# Legend
llxbbox     = 1.0         # x-anchor legend bounding box
llybbox     = 0.95        # y-anchor legend bounding box
llrspace    = 0.          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
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
cbartsize = 'small'       #'small'   # textsize at color bar (default)
htextsize = 'large'  # textsize Hidden/Standard box
ftextsize = 'large'  # textsize of flux labeling
    
# Brewer sequential
cols = color.get_brewer('YlOrRd9', rgb=True)
cols = color.get_brewer('Paired8', rgb=True)
cols = color.get_brewer('RdBu8', rgb=True)
cols = color.get_brewer('RdYlBu8', rgb=True)
cols = color.get_brewer('BlueYellowRed', rgb=True)

# add same alpha as for categories
def add_alpha(col, alpha):
    icol = list(col)
    icol.append(alpha)
    return tuple(icol)
cols = [ add_alpha(c, alphamod) for c in cols ]
# first color is gray
if dowhite:
    colors = [(0.3,0.3,0.3,1.)]
else:
    colors = [(0.7,0.7,0.7,1.)]
colors.extend(cols)

# color map
cmap = mpl.colors.ListedColormap(colors)
    
ifig = 0

# ------------------------------------------
# only first basin's categorization
# ------------------------------------------
ifig += 1
iplot = 0
fig = plt.figure(ifig)

vvspace  = 0.10       # y-space between subplots
nnrow     = 4
for icategorizer,categorizer_name in enumerate(categorizers):
    for ibasin in range(1):
        
        iplot += 1

        pos  = position(nnrow, ncol, iplot, hspace=hspace, vspace=vvspace,
                            left=left, right=right, bottom=bottom, top=top)
        sub    = fig.add_axes(pos)

        cvals = np.append(np.unique(all_cat_categories[icategorizer][ibasin]),np.max(all_cat_categories[icategorizer][ibasin])+1)  
        norm = mpl.colors.BoundaryNorm(cvals, cmap.N)
        mesh = sub.pcolormesh(all_cat_categories[icategorizer][ibasin], cmap=cmap, norm=norm, linewidths=0.0, edgecolor='grey')

        # average percentage of days in category over this basin and all days
        percent_in_cat = np.array([ np.sum(all_cat_categories[icategorizer][ibasin]==cval)*1.0 for cval in cvals ])/np.prod(np.shape(all_cat_categories[icategorizer][ibasin]))*100.

        if (iplot == len(categorizers)):
            xlabel = str2tex('Day of Year',usetex=usetex)
            xticks = np.arange(0,365,50)
            xticklabels = [ str2tex(str(ii), usetex=usetex) for ii in xticks ]
        else:
            xlabel = str2tex('',usetex=usetex)
            xticks = np.arange(0,365,50)
            xticklabels = [ str2tex('', usetex=usetex) for ii in xticks ]
            
        yticks      = np.arange((all_cat_ttime[icategorizer][ibasin][0,0].year-all_cat_ttime[icategorizer][ibasin][0,0].year%10+10)-all_cat_ttime[icategorizer][ibasin][0,0].year,nyears,10)
        yticklabels = [ str2tex(str(ii+all_cat_ttime[icategorizer][ibasin][0,0].year), usetex=usetex) for ii in yticks ]
        ylabel = str2tex('Year',usetex=usetex)

        plt.setp(sub, xlabel=xlabel, xticks=xticks, xticklabels=xticklabels)
        plt.setp(sub, ylabel=ylabel, yticks=yticks, yticklabels=yticklabels)

        # numbering of subplots
        abc2plot(sub, dxabc, dyabc, iplot, lower=False, bold=True, va='top', ha='right', usetex=usetex, mathrm=True, parenthesis='none')      

        # plot variables that are changed for different recipes
        recipe_names = all_cat_recipe_info[icategorizer][0]
        recipe_vars  = all_cat_recipe_info[icategorizer][1]

        # determine on how much to shift color bar up
        cticks = cvals      
        dy = -0.5
        for irecipe,recipe_cat in enumerate(recipes[icategorizer]):

            # skip if this recipe is renames with None
            if ( categorizer_name == 'categorize_meteo_decision_tree_expert'):
                if rename_recipes_expert[recipe_cat] is None:
                    continue
                else:
                    idx = recipe_names.index(recipe_cat)
                    vars_recipe = recipe_vars[idx]
                    nvars = [len(vars_recipe[iitick]) for iitick,itick in enumerate(cticks[1:-1]) ]

                    # set dy for next recipe
                    dy -= 0.6
                    if any(np.array(nvars) > 4):
                        dy -= 0.6
            elif ( categorizer_name == 'categorize_meteo_decision_tree_EEE'):
                if rename_recipes_EEE[recipe_cat] is None:
                    continue
                else:
                    idx = recipe_names.index(recipe_cat)
                    vars_recipe = recipe_vars[idx]
                    nvars = [len(vars_recipe[iitick]) for iitick,itick in enumerate(cticks[1:-1]) ]

                    # set dy for next recipe
                    dy -= 0.6
                    if any(np.array(nvars) > 4):
                        dy -= 0.6
            else:
                raise ValueError('Categorizer not known!')

            

        # Colour bar
        pos  = position(nnrow, ncol, iplot, hspace=hspace, vspace=vvspace,
                            left=left, right=right, bottom=bottom, top=top)
        shrink = 1.0 #0.8
        # [left, bottom, width, height]
        pos[0] += 0.5*(1.-shrink)*pos[2] # shrink
        pos[2]  = shrink*pos[2]
        pos[1] = pos[1]+pos[3]-0.015*dy 
        pos[3] *= 0.16                  # squeeze
        csub = fig.add_axes(pos, frameon=False)
        cbar = plt.colorbar(mesh, cax=csub, orientation='horizontal')
        if dolog:
            # upper labels
            min_cbar = cvals[1]
            max_cbar = cvals[-1]
            lmin = np.log10(min_cbar)
            assert lmin/int(lmin) == 1., '\nmin_cbar must be power of 10.'
            nticks = np.floor(9*(np.log10(max_cbar) - lmin)) +1 # -1
            base = (np.arange(nticks) % 9) + 1 # 1,2,...,9,1,2,...9,1,2,...
            expo = np.arange(nticks) // 9      # 0,0,...,0,1,1,...1,2,2,...
            cticks = base * 10**expo * 10**lmin
            cticknames = [ str(i) for i in cticks ]
            cticknames = [ i if ('1' in i or '3' in i or '5' in i) else '' for i in cticknames ]
            cticknames = str2tex(cticknames, usetex=usetex)
            cbar.set_ticks(cticks)
            cbar.set_ticklabels(cticknames)
            cbar.ax.tick_params(direction='out', length=3, width=1, colors=fgcolor, labelsize=cbartsize, pad=-1.0,
                                bottom='off', top='on', labelbottom='off', labeltop='on')
        else:
            cticks = cvals[:-1] # remove last tick
            cticks = cvals      # dont remove last tick
            cticknames = [ str2tex('$C_{'+str(i)+'}$') if i>0 else str2tex('no data') for i in cticks ]
            cticknames = str2tex(cticknames, usetex=usetex)
            cbar.set_ticks(cticks[:-1]+np.diff(cticks)/2.) # draw ticks in the middle
            cbar.set_ticklabels(cticknames)
            cbar.ax.tick_params(direction='out', length=3, width=1, colors=fgcolor, labelsize=cbartsize, pad=-1.0,
                                bottom=False, top=True, labelbottom=False, labeltop=True)
            for iitick,itick in enumerate(cticks[:-1]):
                tick_location = cticks 
                csub.text(1.0/(len(cvals)-1)*(iitick+0.5), 0.5, str2tex(astr(percent_in_cat[iitick],prec=1)+'%',usetex=usetex),
                              fontsize='xx-small', va='center', ha='center',rotation=0, transform=csub.transAxes)


        # categorizer name
        if categorizer_name == 'categorize_meteo_decision_tree_EEE':
            clab = "EEE-based\n category"
        elif categorizer_name == 'categorize_meteo_decision_tree_expert':
            clab = "Expert-based\n category"
        else:
            clab = 'Not known\n category label'   
        clab = str2tex(clab, usetex=usetex)
        cbar.ax.yaxis.set_label_text(clab, fontsize=cbartsize, va='center', ha='right', rotation=0)
        cbar.ax.yaxis.set_label_coords(-0.03, 0.5, transform=csub.transAxes)

        dy = -0.5
        for irecipe,recipe_cat in enumerate(recipes[icategorizer]):

            # skip if this recipe is renames with None
            if ( categorizer_name == 'categorize_meteo_decision_tree_expert'):
                if rename_recipes_expert[recipe_cat] is None:
                    continue
            elif ( categorizer_name == 'categorize_meteo_decision_tree_EEE'):
                if rename_recipes_EEE[recipe_cat] is None:
                    continue
            else:
                raise ValueError('Categorizer not known!')

            idx = recipe_names.index(recipe_cat)

            vars_recipe = recipe_vars[idx]
            vars_recipe =  [ [ rename_vars[jj] for jj in ii ] for ii in vars_recipe ] # rename vars
            nvars = [len(vars_recipe[iitick]) for iitick,itick in enumerate(cticks[1:-1]) ]

            # Recipe name
            if ( categorizer_name == 'categorize_meteo_decision_tree_expert'):
                lbl    = rename_recipes_expert[recipe_cat].replace('_', ' ').title()+':'
            elif ( categorizer_name == 'categorize_meteo_decision_tree_EEE'):
                lbl    = rename_recipes_EEE[recipe_cat].replace('_', ' ').title()+':'
            else:
                raise ValueError('Categorizer not known!')
            
            iitick = 0
            csub.text(1.0/(len(cvals)-1)*(iitick+1.0), dy, str2tex(lbl,usetex=usetex), fontsize='xx-small', va='top', ha='right',rotation=0, transform=csub.transAxes)

            for iitick,itick in enumerate(cticks[1:-1]):
                
                if nvars[iitick] > 8:
                    cut = nvars[iitick]//3 + 1
                    lbl = (    ','.join([ ('$'+ii+'}$').replace('v','v_{') for ii in vars_recipe[iitick][0*cut:1*cut] ])+'\n'+
                               ','.join([ ('$'+ii+'}$').replace('v','v_{') for ii in vars_recipe[iitick][1*cut:2*cut] ])+'\n'+
                               ','.join([ ('$'+ii+'}$').replace('v','v_{') for ii in vars_recipe[iitick][2*cut:] ]) )
                elif nvars[iitick] > 4:
                    cut = nvars[iitick]//2 + 1
                    lbl = (    ','.join([ ('$'+ii+'}$').replace('v','v_{') for ii in vars_recipe[iitick][0*cut:1*cut] ])+'\n'+
                               ','.join([ ('$'+ii+'}$').replace('v','v_{') for ii in vars_recipe[iitick][1*cut:] ]) )
                else:
                    lbl = ','.join([ ('$'+ii+'}$').replace('v','v_{') for ii in vars_recipe[iitick] ])

                csub.text(1.0/(len(cvals)-1)*(iitick+1.5), dy, str2tex(lbl,usetex=usetex), fontsize='xx-small', va='top', ha='center',rotation=0, transform=csub.transAxes)

            # set dy for next recipe
            dy -= 0.6
            if any(np.array(nvars) > 4):
                dy -= 0.6

            
            

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


