# SAM_SHIP_TRACK_STATS
Source code to compute conditional stats in the ship track region as a function of the aerosol threshold determined immediately before the aerosol injection begins. Code altered from version 6.10.9. 

Affected files:

grid.f90
vars.f90
init.f90
hbuffer.f90
setparm.f90
printout.f90
hbuf_conditionals_init.f90
restart.f90
MICRO_M2005_PA/microphysics.f90
RAD_RRTM4PBL/rad_full.f90
statistics.f90

Overview: Conditionally sample statistics based on deviations from background aerosol concentration before an aerosol injection. A ship track column is present if the weighted average of the bottom N height levels (N controlled by n_avg_lev, with a default value of the bottom 30 grid levels) exceeds the standard deviation threshold (default value: std_t = 3.0). The standard deviation (std_aero) is the maximum deviation in either the horizontal or vertical directions. At each time step, each model grid column is assigned one of the 4 conditions based on the average aerosol concentration through n_avg_lev and the presence or absence of liquid water (threshold is anything greater than 0 g/kg). 

List of unique statistics calculated when doShipTrackConditionals = .true.

‘ST1’ - Fraction of domain covered by cloudy, ship track region.
‘ST2’ - Fraction of domain covered by clear, ship track region.
‘ST3’ - Fraction of domain covered by clear, no-ship track region.
‘ST4’ - Fraction of domain covered by cloudy, no-ship track region.

Mass, energy, and moisture budget variables/residuals

‘WEM’ - Entrainment rate calculated as a residual of the mass budget [m/s].
‘ZB’ - Inversion base height (boundary-layer depth) as calculated using an absolute temperature threshold of 0.03 K/m. Note that the threshold was determined using 5 m vertical grid spacing [m].
‘INVT’ - Inversion base height tendency [m/s].
‘WLO’ - Local mean vertical velocity in the sub-sampled region, calculated by averaging the vertical velocity at the inversion base height [m/s].
‘WLA’  - Large-scale vertical velocity across the entire domain (prescribed profile). Again, the average large-scale velocity in a sub-sampled region is found by averaging the inversion base height large-scale velocity value in each column [m/s].
‘RDIV’ - The net radiative flux is calculated similarly to WLO/WLA by evaluating the flux at the inversion base height in each column and then averaged across all columns. The net radiative flux divergence is then calculated by taking the average radiative flux at the inversion base minus the surface radiative flux [W/m2]
‘LDIV’ - The long-wave component of the radiative flux divergence [W/m2].
‘SDIV’ - The short-wave component of the radiative flux divergence [W/m2].
'TCB' - cloud-base temperature [K]
'QT_A' - Mean boundary-layer total water mixing ratio 
‘TL_T’ - Boundary-layer average (weighted by density and grid spacing) liquid water static energy tendency [K/s]. 
‘SH_F’ - The surface sensible heat flux, normally output in the 2D stats [W/m2].
‘QT_T’ - Boundary-layer average (weighted by density and grid spacing) total water mixing ratio tendency [1/s].
‘LH_F’ - Surface latent heat flux [W/m2].
‘PR_S’ - Surface rain rate [m/s].
‘PR_Z’ - Rain rate (precipitation flux) at the inversion base [m/s].
‘CL_T’ - Cloud-thickness tendency [m/s].
‘CB_T’ - Cloud-base tendency [m/s].
‘LWPT’ - Liquid water path tendency [kg/kg/m2/s].


Vertical profiles of conditionally sampled regions (variables not included in original version of statistics.f90)

‘SW_U’ - Vertical profile of the average upward short-wave radiative flux [W/m2].
‘SW_D’ - Vertical profile of the average downward short-wave flux [W/m2].
‘LW_U’ - Vertical profile of the average upward long-wave flux [W/m2].
‘LW_D’ - Vertical profile of the average downward long-wave flux [W/m2].
‘TS2’ - Liquid potential temperature variance [K2].
‘QS2’ - Total water variance.
‘WS2’ - Vertical velocity variance [m2/s2].
‘US2’ - Zonal velocity variance [m2/s2].
‘VS2’ - Meridional velocity variance [m2/s2].
‘WTL’ - Vertical liquid water static energy flux.
‘WQT’ - Vertical total water mixing ratio flux.
‘WTV’ - Vertical buoyancy flux.
‘PFS’ - Vertical flux of precipitation.

