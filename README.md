# SAM_SHIP_TRACK_STATS
Source code to compute conditional stats in the ship track region as a function of the aerosol threshold determined immediately before the aerosol injection begins. Code altered from version 6.10.9. 

Affected files:

grid.f90
vars.f90
hbuffer.f90
setparm.f90
printout.f90
hbuf_conditionals_init.f90
restart.f90
MICRO_M2005_PA/microphysics.f90 
RAD_RRTM4PBL/rad_full.f90
statistics.f90

Overview: Conditionally sample statistics based on deviations from background aerosol concentration before an aerosol injection. A ship track column is present if the weighted average of the bottom N height levels (N controlled by n_avg_lev, with a default value of the bottom 30 grid levels) exceeds the standard deviation threshold (default value: std_t = 3.0). The standard deviation (std_aero) is the maximum deviation in either the horizontal or vertical directions. At each time step, each model grid column is assigned one of the 4 conditions based on the average aerosol concentration through n_avg_lev and the presence/absence of liquid water (threshold is anything greater than 0 g/kg). Current version will only work with RAD_RRTM4PBL and MICRO_M2005_PA, but this can be easily ported to other radiation/microphysics packages.

List of unique statistics calculated when doShipTrackConditionals = .true.

‘ST1’ - Fraction of domain covered by cloudy, ship track region.
‘ST2’ - Fraction of domain covered by clear, ship track region.
‘ST3’ - Fraction of domain covered by clear, no-ship track region.
‘ST4’ - Fraction of domain covered by cloudy, no-ship track region.

Mass, energy, and moisture budget variables/residuals

‘ZB’ - Inversion base height (boundary-layer depth) as calculated using an absolute temperature threshold of 0.03 K/m. Note that the threshold was determined using 5 m vertical grid spacing [m].
‘ZBCT’ - Cloud top height (~boundary-layer depth) for cases when inversion base height may fail due to decoupling
‘WLO’ - Local mean vertical velocity in the sub-sampled region, calculated by averaging the vertical velocity at the inversion base height [m/s].
‘WLA’  - Large-scale vertical velocity across the entire domain (prescribed profile). Again, the average large-scale velocity in a sub-sampled region is found by averaging the inversion base height large-scale velocity value in each column [m/s].
'WLOT' - Local vertical velocity at cloud top height
'WLAT' - Large-scael vertical velocity at cloud top height
‘RDIV’ - The net radiative flux is calculated similarly to WLO/WLA by evaluating the flux at the inversion base height in each column and then averaged across all columns. The net radiative flux divergence is then calculated by taking the average radiative flux at the inversion base minus the surface radiative flux [W/m2]
‘LDIV’ - The long-wave component of the radiative flux divergence [W/m2].
‘SDIV’ - The short-wave component of the radiative flux divergence [W/m2].
'RDVT' - Net radiative flux divergence calcuated at cloud top height in each column
'LDVT' - Same as LDIV, but evaluated at cloud top height instead of inversion base 
'SDVT' - Same as SDIV, but evaluated at cloud top height instead of inversion base 
'QT_A' - Mean boundary-layer total water mixing ratio 
'TL_A' - Mean boundary-layer liquid potential temperature 
‘SH_F’ - The surface sensible heat flux, normally output in the 2D stats [W/m2].
‘LH_F’ - Surface latent heat flux [W/m2].
‘PR_S’ - Surface rain rate [m/s].
‘PR_Z’ - Rain rate (precipitation flux) at the inversion base [m/s].
'PRZT' - Rain rate at cloud top height [m/s].
'CBH' - Cloud-base height [m].
'TCB' - Cloud-base temperature [K].
'TLCL' - thetal cloud-layer average [K].
'TLSC' - thetal subcloud-layer average [K].
'QTCL' - Total water mixing ratio cloud-layer average [kg/kg].
'QTSC' - Total water mixing ratio subcloud-layer average [kg/kg].
'TCLT' - LWSE cloud-layer tendency [K/s].
'TSCT' - LWSE subcloud-layer tendency [K/s].
'QCLT' - Total water mixing ratio cloud-layer tendency [1/s].
'QSCT' - Total water mixing ratio subcloud-layer tendency [1/s].
‘LWP’ - Liquid water path [kg/m2].
'TADV' - Boundary-layer average thetal temperature advection [K/s].
'QADV' - Boundary-layer average moisture advection [kg/kg/s].
'TCLA' - Cloud-layer average thetal temperature advection [K/s].
'QCLA' - Cloud-layer average moisture advection [kg/kg/s].
'TSCA' - Subcloud-layer average thetal temperature advection [K/s].
'QSCA' - Subcloud-layer average moisture advection [kg/kg/s].

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

