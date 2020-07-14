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
# python figure_1.py -g figure_1_

# -------------------------------------------------------------------------
# General settings
#
dobw      = False # True: black & white
docomp    = True  # True: Print classification on top of modules
dosig     = False # True: add signature to plot
dolegend  = False # True: add legend to each subplot
doabc     = True  # True: add subpanel numbering
dotitle   = True  # True: add catchment titles to subpanels

# -------------------------------------------------------------------------
# Command line arguments
#

import argparse

pngbase   = ''
pdffile   = ''
usetex    = False

parser  = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                  description='''Plot basin shape.''')
parser.add_argument('-g', '--pngbase', action='store',
                    default=pngbase, dest='pngbase', metavar='pngbase',
                    help='Name basis for png output files (default: open screen window).')
parser.add_argument('-p', '--pdffile', action='store',
                    default=pdffile, dest='pdffile', metavar='pdffile',
                    help='Name of pdf output file (default: open screen window).')
parser.add_argument('-t', '--usetex', action='store_true', default=usetex, dest="usetex",
                    help="Use LaTeX to render text in pdf.")

args      = parser.parse_args()
pngbase   = args.pngbase
pdffile   = args.pdffile
usetex    = args.usetex

if (pdffile != '') & (pngbase != ''):
    print('\nError: PDF and PNG are mutually exclusive. Only either -p or -g possible.\n')
    parser.print_usage()
    import sys
    sys.exit()

del parser, args

# -----------------------
# add subolder scripts/lib to search path
# -----------------------
import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append('lib')

import color                            # in lib/
from position      import position      # in lib/
from abc2plot      import abc2plot      # in lib/
from brewer        import get_brewer    # in lib/
from autostring    import astr          # in lib/
from str2tex       import str2tex       # in lib/
from fread         import fread         # in lib/

# import fiona          # some shapefile coordinate stuff
import numpy as np
import xarray as xr
import pandas as pd
import copy                       # deep copy objects, arrays etc
import time
import os
import shapefile
from matplotlib.patches import Polygon, Ellipse
t1 = time.time()

# Projection information are in "GLERL_Map.prj":

# Data could not be read with Python's readshapefile:
#    	ValueError: shapefile must have lat/lon vertices  - it looks like this one has vertices
#	    in map projection coordinates. You can convert the shapefile to geographic
#	    coordinates using the shpproj utility from the shapelib tools
#	    (http://shapelib.maptools.org/shapelib-tools.html)
#
#     ogr2ogr -t_srs EPSG:4326 Contour22_ARCGIS_UTM_EPSG_4326.shp Contour22_ARCGIS_UTM.shp


catchfile_shp_subbasins = ','.join([os.path.join('../data/shapefiles/Contour22_ARCGIS_UTM_EPSG_4326')])
catchfile_shp_merged    = ','.join([os.path.join('../data/shapefiles/Contour22_ARCGIS_UTM_EPSG_4326')])


# -------------------------------------------------------------------------
# Customize plots
#

if (pdffile == ''):
    if (pngbase == ''):
        outtype = 'x'
    else:
        outtype = 'png'
else:
    outtype = 'pdf'
    
# Main plot
nrow        = 5           # # of rows of subplots per figure
ncol        = 2           # # of columns of subplots per figure
hspace      = 0.02         # x-space between subplots
vspace      = 0.05        # y-space between subplots
right       = 0.9         # right space on page
textsize    = 8           # standard text size
dxabc       = 1.0         # % of (max-min) shift to the right from left y-axis for a,b,c,... labels
# dyabc       = -13       # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dyabc       = 0.0         # % of (max-min) shift up from lower x-axis for a,b,c,... labels
dxsig       = 1.23        # % of (max-min) shift to the right from left y-axis for signature
dysig       = -0.05       # % of (max-min) shift up from lower x-axis for signature
dxtit       = 0           # % of (max-min) shift to the right from left y-axis for title
dytit       = 1.3         # % of (max-min) shift up from lower x-axis for title

lwidth      = 1.5         # linewidth
elwidth     = 1.0         # errorbar line width
alwidth     = 0.5         # axis line width
glwidth     = 0.5         # grid line width
msize       = 3.0         # marker size
mwidth      = 1.0         # marker edge width
mcol1       = color.colours('blue')      # primary marker colour
mcol2       = color.colours('red')       # secondary
mcol3       = color.colours('red')       # third
mcols       = ['0.0', '0.4', '0.4', '0.7', '0.7', '1.0']
lcol1       = color.colours('blue')   # primary line colour
lcol2       = '0.0'
lcol3       = '0.0'
lcols       = ['None', 'None', 'None', 'None', 'None', '0.0']
hatches     = [None, None, None, None, None, '//']

# Legend
llxbbox     = 0.0         # x-anchor legend bounding box
llybbox     = 1.0         # y-anchor legend bounding box
llrspace    = 0.          # spacing between rows in legend
llcspace    = 1.0         # spacing between columns in legend
llhtextpad  = 0.4         # the pad between the legend handle and text
llhlength   = 1.5         # the length of the legend handles
frameon     = False       # if True, draw a frame around the legend. If None, use rc

# PNG
dpi         = 600         # 150 for testing
transparent = False
bbox_inches = 'tight'
pad_inches  = 0.035

# Clock options
ymax = 0.6
ntextsize   = 'medium'       # normal textsize
# modules
bmod        = 0.5            # fraction of ymax from center to start module colours
alphamod    = 0.7            # alpha channel for modules
fwm         = 0.05           # module width to remove at sides
ylabel1     = 1.15           # position of module names
ylabel2     = 1.35           # position of class names
mtextsize   = 'large'        # 1.3*textsize # textsize of module labels
# bars
bpar        = 0.4            # fraction of ymax from center to start with parameter bars
fwb         = [0.7,0.4,0.3]  # width of bars
plwidth     = 0.5
# parameters in centre
bplabel     = 0.1            # fractional distance of ymax of param numbers in centre from 0-line
ptextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of param numbers in centre
# yaxis
space4yaxis = 2              # space for y-axis (integer)
ytextsize   = 'medium'       # 'small' # 0.8*textsize # textsize of y-axis
sig         = 'J Mai' # sign the plot

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
        mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        #mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
elif (outtype == 'png'):
    mpl.use('Agg') # set directly after import matplotlib
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(8.27,11.69)) # a4 portrait
    if usetex:
        mpl.rc('text', usetex=True)
    else:
        mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        #mpl.rc('font',**{'family':'serif','serif':['times']})
    mpl.rc('text.latex', unicode=True)
    mpl.rc('savefig', dpi=dpi, format='png')
else:
    import matplotlib.pyplot as plt
    mpl.rc('figure', figsize=(4./5.*8.27,4./5.*11.69)) # a4 portrait
mpl.rc('font', size=textsize)
mpl.rc('lines', linewidth=lwidth, color='black')
mpl.rc('axes', linewidth=alwidth, labelcolor='black')
mpl.rc('path', simplify=False) # do not remove

from matplotlib.patches import Rectangle, Circle, Polygon
from mpl_toolkits.basemap import Basemap


# Read catchments as shape files
#  --> during plotting

if (catchfile_shp_subbasins != ''):
    catchfile_shp_subbasins  = catchfile_shp_subbasins.split(',')
    catchfile_shp_subbasins  = [ i.replace(' ','') for i in catchfile_shp_subbasins ]
    catchfile_shp_merged     = catchfile_shp_merged.split(',')
    catchfile_shp_merged     = [ i.replace(' ','') for i in catchfile_shp_merged ]
    ncatchfile_shp = len(catchfile_shp_subbasins)
else:
    ncatchfile_shp = 0

# colors
if dobw:
    c = np.linspace(0.2, 0.85, nmod)
    c = np.ones(nmod)*0.7
    c = [ str(i) for i in c ]
    ocean_color = '0.1'
else:
    # c = [(165./255.,  0./255., 38./255.), # interception
    #      (215./255., 48./255., 39./255.), # snow
    #      (244./255.,109./255., 67./255.), # soil moisture
    #      (244./255.,109./255., 67./255.), # soil moisture
    #      (253./255.,174./255., 97./255.), # direct runoff
    #      (254./255.,224./255.,144./255.), # Evapotranspiration
    #      (171./255.,217./255.,233./255.), # interflow
    #      (116./255.,173./255.,209./255.), # percolation
    #      ( 69./255.,117./255.,180./255.), # routing
    #      ( 49./255., 54./255.,149./255.)] # geology
    c  = get_brewer('rdylbu11', rgb=True)
    tmp = c.pop(5)   # rm yellow
    np.random.shuffle(c)
    
    #c.insert(2,c[2]) # same colour for both soil moistures
    ocean_color = (151/256., 183/256., 224/256.)
    # ocean_color = color.get_brewer('accent5', rgb=True)[-1]

    cc = color.get_brewer('dark_rainbow_256', rgb=True)
    cc = cc[::-1] # reverse colors
    cmap = mpl.colors.ListedColormap(cc)


    # colors for each sub-basin from uwyellow4 (228,180,42) to gray
    graylevel = 0.2
    uwyellow = [251,213,79]
    cc = [ (    (uwyellow[0]+ii/(22.-1)*(256*graylevel-uwyellow[0]))/256.,
                (uwyellow[1]+ii/(22.-1)*(256*graylevel-uwyellow[1]))/256.,
                (uwyellow[2]+ii/(22.-1)*(256*graylevel-uwyellow[2]))/256.) for ii in range(22) ]
    cmap = mpl.colors.ListedColormap(cc)



    
    cc = color.get_brewer('Paired8', rgb=True)
    cmap = mpl.colors.ListedColormap(cc)




    cc = color.get_brewer('RdYlBu7', rgb=True)    # need to be at least 7 colors
    cmap = mpl.colors.ListedColormap(cc)

# selected catchments
head_names = [ "MBLANC",  "PERIB", "LM", "PD" ]
head_catch = [       16,       19,    9,   17 ]   # which subcatchment in shapefile
head_color = [        7,        6,    1,   0  ]   # make sure they receive the same color as in other plots   # Paired8
head_color = [        3,        2,    1,   0  ]   # make sure they receive the same color as in other plots   # RdYlBu7

    
# -------------------------------------------------------------------------
# Plot
#

if (outtype == 'pdf'):
    print('Plot PDF ', pdffile)
    pdf_pages = PdfPages(pdffile)
elif (outtype == 'png'):
    print('Plot PNG ', pngbase)
else:
    print('Plot X')
# figsize = mpl.rcParams['figure.figsize']

ifig = 0

# -------------------------------------------------------------------------
# Fig 1
#
ifig += 1
iplot = 0
print('Plot - Fig ', ifig)
fig = plt.figure(ifig)

# -----------------------------------------
# Quebec - overview
# -----------------------------------------
iplot += 1

#     [left, bottom, width, height]
pos = [0.1,0.8,0.45,0.15]
sub = fig.add_axes(pos)
#sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)) #, axisbg='none')

# Map: Europe
# m = Basemap(projection='lcc', llcrnrlon=-9, llcrnrlat=35.6, urcrnrlon=25.3, urcrnrlat=53,
#             lat_1=50, lat_2=70, lon_0=0, resolution='i') # Lambert conformal
# Map: USA
# m = Basemap(projection='lcc',
#             llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64, urcrnrlat=49,
#             lat_1=33, lat_2=45, lon_0=-95,
#             resolution='i') # Lambert conformal
# Map: Canada - Saguenay-Lac-St-Jean region
llcrnrlon =  -95.0
urcrnrlon =  -49.0
llcrnrlat =   37.0
urcrnrlat =   63.0
lat_1     =   50.0  # first  "equator"
lat_2     =   50.0  # second "equator"
lat_0     =   50.0  # center of the map
lon_0     =  -72.0  # center of the map
# m = Basemap(projection='lcc',
#             llcrnrlon=-80, llcrnrlat=43, urcrnrlon=-75, urcrnrlat=47,
#             lon_0=-77.5, lat_0=43, 
#             lat_1=44, lat_2=44, 
#             resolution='i') # Lambert conformal
map4 = Basemap(projection='lcc',
            llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
            lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
            area_thresh=1000., # only large lakes
            resolution='i') # Lambert conformal
          
# draw parallels and meridians.
# labels: [left, right, top, bottom]
map4.drawparallels(np.arange(-80.,81.,7.),  labels=[1,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
map4.drawmeridians(np.arange(-180.,181.,10.),labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

# draw cooastlines and countries
map4.drawcoastlines(linewidth=0.3)
map4.drawmapboundary(fill_color=ocean_color, linewidth=0.3)
map4.drawcountries(color='black', linewidth=0.3)
map4.fillcontinents(color='white', lake_color=ocean_color)

coord_catch   = []
subbasin_info = []
for i, catch in enumerate(catchfile_shp_subbasins[0:1]):
    print('Read ', catch)        
    s = map4.readshapefile(catch, 'cequeau_subbasins', linewidth=0.0, color='none', drawbounds=False)
    nshapes = len(map4.cequeau_subbasins)
    print('   --> has: ',nshapes,' sub-basins')

    for info, shape in zip(map4.cequeau_subbasins_info, map4.cequeau_subbasins):
        subbasin_info.append(info)
        coord_catch.append(shape)

    print('   --> has: ',len(coord_catch),' filtered sub-basins')

# Catchments
xmin =  99999999999.0
xmax = -99999999999.0
ymin =  99999999999.0
ymax = -99999999999.0
nshapes = np.shape(coord_catch)[0]
for ishape in range(nshapes):
    xy = np.array(coord_catch[ishape])

    xxmin = np.min(xy[:,0])
    xxmax = np.max(xy[:,0])
    yymin = np.min(xy[:,1])
    yymax = np.max(xy[:,1])
    if xmin > xxmin:
        xmin = xxmin
    if xmax < xxmax:
        xmax = xxmax
    if ymin > yymin:
        ymin = yymin
    if ymax < yymax:
        ymax = yymax
    
    # add catchment shape to plot
    if ishape in range(nshapes):
        if ishape in head_catch:
            icolor = cc[head_color[head_catch.index(ishape)]]
        else:
            icolor = cc[ishape%len(cc)]
        sub.add_patch(Polygon(xy[::1], facecolor=cc[0], edgecolor='black', linewidth=0.0,zorder = 300, alpha=0.6))  #linewidth=0.3

print("   --> lon range = [",map4(xmin,ymin,inverse=True)[0],",",map4(xmax,ymax,inverse=True)[0],"]")
print("   --> lat range = [",map4(xmin,ymin,inverse=True)[1],",",map4(xmax,ymax,inverse=True)[1],"]")

# Fake subplot for numbering
if doabc:
    
    #  [left, bottom, width, height]
    # lsub = fig.add_axes([pos[0], pos[1]-0.02, pos[2], 0.05])
    lsub = fig.add_axes([pos[0]-0.36, pos[1]+pos[3]-0.01, pos[2]-0.01, pos[3]])

    lsub.set_xlim([0,1])
    lsub.set_ylim([0,1])

    # subplot numbering
    abc2plot(lsub, dxabc, dyabc, iplot, lower=False,
                 bold=True, large=True,
                 mathrm=True, usetex=usetex,
                 horizontalalignment='left', verticalalignment='bottom')

    lsub.set_title('')
    lsub.set_xlabel('')
    lsub.set_ylabel('')
    lsub.set_xticks([])
    lsub.set_yticks([])
    lsub.set_axis_off()
    
# -----------------------------------------
# Subcatchments of Saguenay-Lac-St-Jean region - highlight only 4 head catchments with labels
# -----------------------------------------
iplot += 1

#     [left, bottom, width, height]
pos = [0.4,0.8,0.45,0.15]
sub = fig.add_axes(pos)
#sub = fig.add_axes(position(nrow,ncol,iplot,hspace=hspace,vspace=vspace)) #, axisbg='none')

# Map: Europe
# m = Basemap(projection='lcc', llcrnrlon=-9, llcrnrlat=35.6, urcrnrlon=25.3, urcrnrlat=53,
#             lat_1=50, lat_2=70, lon_0=0, resolution='i') # Lambert conformal
# Map: USA
# m = Basemap(projection='lcc',
#             llcrnrlon=-119, llcrnrlat=22, urcrnrlon=-64, urcrnrlat=49,
#             lat_1=33, lat_2=45, lon_0=-95,
#             resolution='i') # Lambert conformal
# Map: Canada - Saguenay-Lac-St-Jean region
llcrnrlon =  -77.0
urcrnrlon =  -67.0
llcrnrlat =   47.0
urcrnrlat =   53.0
lat_1     =   (llcrnrlat+urcrnrlat)/2.0  # first  "equator"
lat_2     =   (llcrnrlat+urcrnrlat)/2.0  # second "equator"
lat_0     =   (llcrnrlat+urcrnrlat)/2.0  # center of the map
lon_0     =   (llcrnrlon+urcrnrlon)/2.0  # center of the map
# m = Basemap(projection='lcc',
#             llcrnrlon=-80, llcrnrlat=43, urcrnrlon=-75, urcrnrlat=47,
#             lon_0=-77.5, lat_0=43, 
#             lat_1=44, lat_2=44, 
#             resolution='i') # Lambert conformal
map4 = Basemap(projection='lcc',
            llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat,
            lat_1=lat_1, lat_2=lat_2, lat_0=lat_0, lon_0=lon_0,
            area_thresh=100., # only large lakes
            resolution='i') # Lambert conformal
          
# draw parallels and meridians.
# labels: [left, right, top, bottom]
map4.drawparallels(np.arange(-80.,81.,2.),  labels=[1,0,0,0], dashes=[1,1], linewidth=0.25, color='0.5')
map4.drawmeridians(np.arange(-180.,181.,5.),labels=[0,0,0,1], dashes=[1,1], linewidth=0.25, color='0.5')

# draw cooastlines and countries
map4.drawcoastlines(linewidth=0.3)
map4.drawmapboundary(fill_color=ocean_color, linewidth=0.3)
map4.drawcountries(color='black', linewidth=0.3)
map4.fillcontinents(color='white', lake_color=ocean_color)

coord_catch   = []
subbasin_info = []
for i, catch in enumerate(catchfile_shp_subbasins[0:1]):
    print('Read ', catch)        
    s = map4.readshapefile(catch, 'cequeau_subbasins', linewidth=0.0, color='none', drawbounds=False)
    nshapes = len(map4.cequeau_subbasins)
    print('   --> has: ',nshapes,' sub-basins')

    for info, shape in zip(map4.cequeau_subbasins_info, map4.cequeau_subbasins):
        subbasin_info.append(info)
        coord_catch.append(shape)

    print('   --> has: ',len(coord_catch),' filtered sub-basins')

# Catchments
xmin =  99999999999.0
xmax = -99999999999.0
ymin =  99999999999.0
ymax = -99999999999.0
nshapes = np.shape(coord_catch)[0]
count = 0
for ishape in range(nshapes):
    xy = np.array(coord_catch[ishape])

    xxmin = np.min(xy[:,0])
    xxmax = np.max(xy[:,0])
    yymin = np.min(xy[:,1])
    yymax = np.max(xy[:,1])
    if xmin > xxmin:
        xmin = xxmin
    if xmax < xxmax:
        xmax = xxmax
    if ymin > yymin:
        ymin = yymin
    if ymax < yymax:
        ymax = yymax
    
    # add catchment shape to plot
    if ishape in range(nshapes):
        if ishape in head_catch:
            icolor = cc[head_color[head_catch.index(ishape)]]
        else:
            icolor = cc[ishape%len(cc)]
        if ishape in head_catch:
            sub.add_patch(Polygon(xy[::1], facecolor=icolor, edgecolor='black', linewidth=0.3,zorder = 300, alpha=0.6))

            xpt, ypt  = [ np.mean(xy[:,0]), np.mean(xy[:,1]) ]   # center of shape
            x2,  y2   = map4(-74.2, 52.8-0.33*head_catch.index(ishape)) #[ np.mean(xy[:,0])*0.9, np.mean(xy[:,1])*1.0 ]   # position of text
            sub.annotate(head_names[head_catch.index(ishape)],
                xy=(xpt, ypt),   xycoords='data',
                xytext=(x2, y2), textcoords='data',
                fontsize=6,
                verticalalignment='center',horizontalalignment='right',
                arrowprops=dict(arrowstyle="->",relpos=(1.0,0.5),linewidth=0.6),
                zorder=400
                )
            count += 1
        else:
            sub.add_patch(Polygon(xy[::1], facecolor='None', edgecolor='black', linewidth=0.3,zorder = 300, alpha=0.6))

print("   --> lon range = [",map4(xmin,ymin,inverse=True)[0],",",map4(xmax,ymax,inverse=True)[0],"]")
print("   --> lat range = [",map4(xmin,ymin,inverse=True)[1],",",map4(xmax,ymax,inverse=True)[1],"]")

# set title as time step
# sub.set_title(timestep,fontsize=textsize)

# Fake subplot for numbering
if doabc:
    
    #  [left, bottom, width, height]
    # lsub = fig.add_axes([pos[0], pos[1]-0.02, pos[2], 0.05])
    lsub = fig.add_axes([pos[0]-0.35, pos[1]+pos[3]-0.01, pos[2]-0.01, pos[3]])

    lsub.set_xlim([0,1])
    lsub.set_ylim([0,1])

    # subplot numbering
    abc2plot(lsub, dxabc, dyabc, iplot, lower=False,
                 bold=True, large=True,
                 mathrm=True, usetex=usetex,
                 horizontalalignment='left', verticalalignment='bottom')

    lsub.set_title('')
    lsub.set_xlabel('')
    lsub.set_ylabel('')
    lsub.set_xticks([])
    lsub.set_yticks([])
    lsub.set_axis_off()
    

if (outtype == 'pdf'):
    pdf_pages.savefig(fig)
    plt.close(fig)
elif (outtype == 'png'):
    pngfile = pngbase+"{0:04d}".format(ifig)+".png"
    fig.savefig(pngfile, transparent=transparent, bbox_inches=bbox_inches, pad_inches=pad_inches)
    plt.close(fig)


# -------------------------------------------------------------------------
# Finished
#

if (outtype == 'pdf'):
    pdf_pages.close()
elif (outtype == 'png'):
    pass
else:
    plt.show()

t2    = time.time()
strin = '[m]: '+astr((t2-t1)/60.,1) if (t2-t1)>60. else '[s]: '+astr(t2-t1,0)
print('Time ', strin)




    
