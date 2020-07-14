# Application of Parameter Screening To Derive Optimal Initial State Adjustments for Streamflow Forecasting
*by J. Mai<sup> 1</sup>,  R. Arsenault<sup> 2</sup>, B.A. Tolson<sup> 1</sup>, M. Latraverse<sup> 3</sup>, and K. Demeester<sup> 3</sup>*<br><br>
*<sup> 1</sup> Dept. Civil and Environmental Engineering, University of Waterloo, Waterloo, ON, Canada.*<br>
*<sup> 2</sup> Dépt. de génie de la construction, École de Technologie Supérieure, Montreal, QC, Canada.*<br>
*<sup> 3</sup> Aluminium division, Rio Tinto, Jonquière, QC, Canada.*<br>

## Abstract
Streamflow forecasting is essential in many applications, including the management of hydropower reservoirs and of flood risks. To increase the performance of these hydrologic forecasts, a hydrological model is typically updated with optimal initial model states. These states are usually selected a priori and adjusted to obtain the optimal initial setup. The model states to adjust are generally selected by experts with deep knowledge of the model's behaviour. These relationships are rarely documented and formalized leading to difficulties in reproducibility and transferability to other models.
The method here provides a fully automatized system to find the model states based on the parameter-free sensitivity method of Efficient Elementary Effects (EEE), which identifies the most informative variables of the hydrologic model CEQUEAU. The analysis is applied to 1826 model setups for four Canadian basins. The identified states are formalized using clustering techniques to obtain adjustment "recipes" conditioned on given meteorological conditions. Several recipes are tested and compared to expert-based recipes. The recipes are tested \replaced{regarding}{according to} their ability to (a) obtain optimal initial states and (b) their forecast performance. The results show: 1) The EEE-based recipes perform similarly or better than the expert-based recipes in both objectives. 2) The recipes are similar across the four basins highlighting their potential of spatial transferability. 3) The adjustment of hydrologic model inputs significantly decreases forecast performance. The proposed approach to use EEE-based recipes to obtain optimal initial states for forecasting is reproducible, objective, and suitable to be applied to other models and basins.  The full paper can be found [here](https://doi.org/10.1002/???). 

## Creating Plots
This GitHub contains all scripts and data to reproduce the plots in the paper and Supplementary Material. Please see the main plotting script [here](https://github.com/julemai/EEE-DA/scripts/plot.sh) and select the figures you want to plot. The final figures will then be found in the folder `figures`. All the data used to produce those figures can be found in the folder `data`.

## Citation

### Journal Publication
J. Mai,  R. Arsenault, B.A. Tolson, M. Latraverse, and K. Demeester (2020) .<br>
Application of Parameter Screening To Derive Optimal Initial State Adjustments for Streamflow Forecasting.<br>
*Water Resources Research*, ??, ???-???.<br>
https://doi.org/10.1002/???.

### Code Publication
J. Mai,  R. Arsenault, B.A. Tolson, M. Latraverse, and K. Demeester (2020) .<br>
Application of Parameter Screening To Derive Optimal Initial State Adjustments for Streamflow Forecasting (v1.0).<br>
*Zenodo*<br>
https://doi.org/10.5281/zenodo.3945225
