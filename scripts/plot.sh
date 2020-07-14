#!/bin/bash
#
# Produces plots for publication:
#
#      J. Mai,  R. Arsenault, B.A. Tolson, M. Latraverse, and K. Demeester (2020):
#      Application of Parameter Screening To Derive Optimal Initial State Adjustments for Streamflow Forecasting
#      Water Resources Research, ??, ???-???.
#
#
# Copyright 2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
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

set -e
#
prog=$0
pprog=$(basename ${prog})
dprog=$(dirname ${prog})
isdir=${PWD}
pid=$$

dotex=1      # LaTeX for font setting

# -----------------------------------------
# Set figure(s) you want to plot to "1" and run ./plot.sh
# -----------------------------------------

dofig1=0     # figure_1:               (A) Overview of the Saguenay–Lac-Saint-Jean watershed in Northern Quebec,
#                                      Canada. (B) Highlighted four subbasins used for this study. Details about the basins can be
#                                      found in Table 1.

dofig2=1     # figure_2:               Flowchart of the hydrologic model CEQUEAU. The adjustment variables vi in this
#                                      study are indicated by the labels in circles. The state variables considered in this study are the
#                                      level of the soil reservoir (adjusted through v1 in [mm]), the level of the groundwater reservoir
#                                      (adjusted through v2 in [mm]), the amount of snow in the snowpack (adjusted through v3 in
#                                      [mm]), and the state of snow in the snowpack (adjusted through v4 in [deg C]). The input vari-
#                                      ables adjusted in this study are precipitation (adjusted through v5 in [mm]) and temperature
#                                      (adjusted through v6 in [deg C]). The temperature is influencing evapotranspiration and snowmelt
#                                      processes. Hence adjustment v6 appears twice in the flow chart. More details about ranges, units,
#                                      and initial values of these adjustments can be found in Table 2.

dofig3=0    # figure_3:                Hydrographs with original and corrected forcings (red and blue lines, respectively)
#                                      for two head catchments (A) Peribonka and (B) Montagnes Blanches compared to the observed
#                                      discharge (gray circles) for two example years 2012 and 2013. Several performance metrics of
#                                      the hydrographs are given as labels above each panel (colors corresponding to the discharge
#                                      curves) demonstrating the performance increase achieved by correcting the meteorologic forcings
#                                      following the approach explained in Section 2.3.2.

dofig4=1     # figure_4:               Flowchart describing the experiments of this work. (A) We identify recipes, i.e. a
#                                      set of (sensitive) parameters to optimize based on meteorologic conditions. A proposed recipe
#                                      is identified using the Efficient Elementary Effects method (EEE); another recipe is determined
#                                      by a forecasting analyst based on expert knowledge. (B) These two recipes are then tested re-
#                                      garding their ability to yield optimal initial states (through automatic optimization) and the
#                                      performance of forecasts issued using these optimal initial states. The sections describing details
#                                      of the methods are indicated left of each box of the flowchart.
#                                      a) For these analyses and model runs adjusted forcings as explained in Sec. 2.3.2 are used.
#                                      b) For these analyses and model runs raw, non-adjusted forcings are used.
#                                      c) Instead of analyzing the full year only the summer period is considered here. This is matching
#                                      the time period where the adjusted forcings are used in the operational setup at Rio Tinto. The
#                                      winter period is not using adjusted forcings.
#                                      d) Instead of using random basins, we used the four basins described in Sec. 2.1 assuming that
#                                      the results are general enough to be applied in any other basin of the Saguenay–Lac-Saint-Jean
#                                      watershed.

dofig5=0     # figure_5:               Model variable screening using Efficient Elementary Effects (EEE) for the four
#                                      example subbasins (Tab. 1 and Fig. 1B) over the course of five years. The screening is based on
#                                      8-day simulations (T − 7 to T − 0) of streamflow starting at the indicated date (x-axis) using pre-
#                                      processed meteorological forcings and an open-loop run as initial condition. Dark red indicates
#                                      that the variable is informative (sensitive) while yellow indicates that the variable is identified to
#                                      be insensitive using the EEE analysis. Details about the variables can be found in Table 2.

dofig6_S1=0  # figure_6 and figure_S1: Sensitivity EEE for each model variable v1 to v6 (rows) and basin (column) de-
#                                      pending on day of the year (DOY) and temperature of simulation start day T − 7 (circles). The
#                                      size of the markers indicates the sensitivity while the color indicates if the variable was identified
#                                      to be insensitive (yellow) or informative (red) using the Efficient Elementary Effects (EEE). The
#                                      median sensitivity EEE for each parameter and basin is given above each panel. To identify the
#                                      temperature range a parameter should be considered for optimization, the first and 99th weighted
#                                      percentiles of temperature p^T_1 and p^T_99 (horizontal lines ) are used (added as labels below and
#                                      above the horizontal lines). To derive ranges for the DOY, a Density-Based Spatial Clustering of
#                                      Applications with Noise (DBSCAN) is applied (vertical lines). The limits of each DOY range are
#                                      the first and 99th percentile of the cluster k identified by DBSCAN (p^{C_k}_1 and p^{C_k}_99 ) and are added
#                                      as labels right and left of the vertical lines.

dofig7=0     # figure_7:               The classification of days based on (A) expert and (B) EEE-based classification
#                                      schemes (Tab. 3) for the data period with available meteorology (1953 - 2019) for one example
#                                      watershed (PD; Passes Dangereuses). The colors indicate different classes. The percentage of
#                                      days in each class Ci are added as labels to the color bar. Below the color bar the model vari-
#                                      ables adjusted for different recipes tested are added as labels. Details about the variables can be
#                                      found in Tab. 2.

dofig8=0     # figure_8:               Initial state adjustment performance for summer months (Jul to Nov). The (A)
#                                      nominal and (B) relative improvement of the initial state performance metric (WAE; Eq. 2)
#                                      compared to an open-loop simulation using six different optimization strategies in four basins
#                                      (Tab. 1, Fig. 1). The first three recipes (warm colors) are based on the Efficient Elementary Ef-
#                                      fects analysis determining which parameters are adjusted/ optimized under which meteorological
#                                      conditions at the simulation start day T − 7 while the three other recipes (cold colors) use expert
#                                      knowledge using both the same categorization of meteorological conditions (temperature and
#                                      precipitation) but different sets of parameters optimized in each meteorologic category. The me-
#                                      dian performance improvements are added as labels to each boxplot. The whiskers of the boxplot
#                                      indicate the 10th and 90th percentile, the box itself envelops the 25th and 75th percentile of the
#                                      data points. In total 1826 data points (simulation start days) are represented by each box. The
#                                      percentage of successful initial state adjustments and a summary of which inputs/states were
#                                      adjusted for the different recipes can be found in Table 4.

dofig9=0     # figure_9:               Forecast performance for summer months (Jul to Nov). The success rate is count-
#                                      ing if a forecast with initial state adjustments following different recipes (colored lines) is yielding
#                                      volume errors that are at least 5% better than the open-loop simulation for lead times up to 52
#                                      weeks for four subbasins (A-D). Success rates above 50% indicate strategies that reliably outper-
#                                      form the open-loop and are hence desirable.

dofig10=0    # figure_10:              Forecast performance for summer months (Jul to Nov). The volume error E_V(L)
#                                      (Eq. 7) between observed streamflow (hindcast period) and the open-loop simulations (gray
#                                      boxes) and forecasts issued using different recipes to adjust the initial states (colored boxes) for
#                                      four catchments (A-D) over lead times L up to one year. The optimal volume error of zero is
#                                      added to ease comparison (black line). The whiskers of the boxplot indicate the 10th and 90th
#                                      percentile, the box itself envelops the 25th and 75th percentile of the data points. In total up
#                                      to 1826 data points (simulation start days) are represented by each box. The dashed horizontal
#                                      lines are added for reference to indicate half of the mean annual daily discharge Q for each basin.
#                                      Details about the basins are given in Table 1. Details about the recipes and which variables v are
#                                      adjusted can be found in Table 4.

verbose=0 # 0: pipe stdout and stderr to /dev/null
          # 1: pipe stdout to /dev/null
          # 2: print everything

# ---------------------------------------------------------------------
# Treat switches
# ---------------------------------------------------------------------
if [[ ${dotex} -eq 1 ]] ; then
    texit='-t'
else
    texit=''
fi

pipeit=''
if [[ ${verbose} -eq 0 ]] ; then pipeit=' > /dev/null 2>&1' ; fi
if [[ ${verbose} -eq 1 ]] ; then pipeit=' > /dev/null' ; fi

# Figures
if [[ ${dofig1} -eq 1 ]] ; then
    echo ''
    echo 'Figure 1 in progress...'
    # ---------------------------------------------------------------------
    # Study domain and highlighted, annotated subwatersheds
    # ---------------------------------------------------------------------
    python figure_1.py -g figure_1_
    convert -trim figure_1_0001.png figure_1_0001-crop.png
    mv figure_1_0001-crop.png ../figures/figure_1.png
    rm figure_1_0001.png
fi

if [[ ${dofig2} -eq 1 ]] ; then
    echo ''
    echo 'Figure 2 in progress...'
    # --------------------------------------------------------------------
    # Flowchart CEQUEAU
    # --------------------------------------------------------------------
    pdflatex ../scripts/figure_2/figure_2.tex  ${pipeit}
    pdfcrop figure_2.pdf
    mv figure_2-crop.pdf ../figures/figure_2.pdf
    rm figure_2.pdf
    rm figure_2.log
    rm figure_2.aux
fi

basins="[Peribonka, Montagnes Blanches]"
if [[ ${dofig3} -eq 1 ]] ; then
    echo ''
    echo 'Figure 3 in progress...'
    # ---------------------------------------------------------------------
    # Simulated streamflow before and after meteo correction 
    # ---------------------------------------------------------------------
    python figure_3.py -b "${basins}" -m ../data/open_loop/resultats_1953-01-01_original.nc -n ../data/open_loop/resultats_1953-01-01_corrected.nc -t 2012-01-01_2013-12-31 -i ../data/observations/ObservedData.nc -o figure_3.pdf
    pdfcrop figure_3.pdf
    mv figure_3-crop.pdf ../figures/figure_3.pdf
    rm figure_3.pdf
fi

if [[ ${dofig4} -eq 1 ]] ; then
    echo ''
    echo 'Figure 4 in progress...'
    # --------------------------------------------------------------------
    # Flowchart: Experimental setup paper
    # --------------------------------------------------------------------
    pdflatex ../scripts/figure_4/figure_4.tex  ${pipeit}
    pdfcrop figure_4.pdf
    mv figure_4-crop.pdf ../figures/figure_4.pdf
    rm figure_4.pdf
    rm figure_4.log
    rm figure_4.aux
fi

if [[ ${dofig5} -eq 1 ]] ; then
    echo ''
    echo 'Figure 5 in progress...'
    # ---------------------------------------------------------------------
    # Results of EEE analysis of 7-day simulations regarding streamflow
    # ---------------------------------------------------------------------
    # python figure_5.py -i ../../CRD-DA/EEE_results_2011-2015/ -o ../data/eee/eee_results.nc -p figure_5.pdf -t pdf -u
    python figure_5.py -i ../data/eee/ -o ../data/eee/eee_results.nc -p figure_5.pdf -t pdf 
    pdfcrop figure_5.pdf
    mv figure_5-crop.pdf ../figures/figure_5.pdf
    rm figure_5.pdf
fi

if [[ ${dofig6_S1} -eq 1 ]] ; then
    echo ''
    echo 'Figure 6 and S1 in progress...'
    # ---------------------------------------------------------------------
    # Results of EEE analysis correlated with DOY, temp, and precip to identify recipes
    # ---------------------------------------------------------------------
    python figure_6_S1.py  -i ../data/eee/ -o ../data/eee/eee_results.nc -p figure_6_S1_ -t png 
    mv figure_6_S1_0001.png ../figures/figure_S1.png
    mv figure_6_S1_0002.png ../figures/figure_6.png
fi

categorizer='[categorize_meteo_decision_tree_expert, categorize_meteo_decision_tree_EEE]'
recipes="[[recipe_5,recipe_3,recipe_1],[recipe_1,recipe_4,recipe_3]]"
if [[ ${dofig7} -eq 1 ]] ; then
    echo ''
    echo 'Figure 7 in progress...'
    # ---------------------------------------------------------------------
    # Classification of basins (EEE and expert classes)
    # ---------------------------------------------------------------------

    python figure_7.py  -i ../data/observations/meteo.nc -a -p figure_7.pdf -t pdf -c "${categorizer}" -r "${recipes}"
    pdfcrop figure_7.pdf
    mv figure_7-crop.pdf ../figures/figure_7.pdf
    rm figure_7.pdf
fi

basins="[Passes Dangereuses, Peribonka, Lac Manouane, Montagnes Blanches]"
if [[ ${dofig8} -eq 1 ]] ; then
    echo ''
    echo 'Figure 8 in progress...'
    # ---------------------------------------------------------------------
    # Calibration results (EEE and expert classes)
    # ONLY calibrated starting days are considered
    # ---------------------------------------------------------------------
    #python figure_8.py -i "[../data/initial_state_adjustments/initial_state_adjustments_EEE.nc, ../data/initial_state_adjustments/initial_state_adjustments_expert.nc]" -b "${basins}" -p figure_8.pdf -t pdf  --only_summer
    python figure_8.py -i "[../data/forecast/forecast_performance_EEE.nc, ../data/forecast/forecast_performance_expert.nc]" -b "${basins}" -p figure_8.pdf -t pdf --only_summer
    pdfcrop figure_8.pdf
    mv figure_8-crop.pdf ../figures/figure_8.pdf
    rm figure_8.pdf
fi


basins="[Passes Dangereuses, Peribonka, Lac Manouane, Montagnes Blanches]"
if [[ ${dofig9} -eq 1 ]] ; then
    echo ''
    echo 'Figure 9 in progress...'

    # ---------------------------------------------------------------------
    # Forecast performance (EEE and expert classes)
    #    EEE:    recipe_1, recipe_2, recipe_3, recipe_4
    #    expert: recipe_1, recipe_2, recipe_3, recipe_5
    # ---------------------------------------------------------------------
    
    python figure_9.py -i "[../data/forecast/forecast_performance_EEE.nc, ../data/forecast/forecast_performance_expert.nc]" -b "${basins}" -p figure_9.pdf -t pdf
    pdfcrop figure_9.pdf
    mv figure_9-crop.pdf ../figures/figure_9.pdf  
    rm figure_9.pdf
    
fi

basins="[Passes Dangereuses, Peribonka, Lac Manouane, Montagnes Blanches]"
if [[ ${dofig10} -eq 1 ]] ; then
    echo ''
    echo 'Figure 10 in progress...'
    
    # ---------------------------------------------------------------------
    # Forecast performance volume error  (EEE and expert classes; different recipes)
    #
    # ALL starting days considered (not just calibrated or successfully calibrated)
    # ---------------------------------------------------------------------
    
    python figure_10.py -i "[../data/forecast/forecast_performance_EEE.nc, ../data/forecast/forecast_performance_expert.nc]" -b "${basins}" -p figure_10.pdf -t pdf 
    pdfcrop figure_10.pdf
    mv figure_10-crop.pdf ../figures/figure_10.pdf
    rm figure_10.pdf

fi
