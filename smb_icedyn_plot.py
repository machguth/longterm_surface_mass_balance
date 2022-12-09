"""
    *** Numerical calculation of surface mass balance profiles using a simple degree-day approach. ***

Purpose of the calculation is to investigate the influence of different annual temperature amplitudes on the
surface mass balance gradient in the abaltion area.

Surface mass balance is used to drive an ice flow model and thereby to study the influence
of changing surface mass balance gradients on advance and reatreat behaviour of glaciers during the last glacial period

This code plots the calculated ice volume and glacier length, as simulated by ELMER-ICE in 2D mode.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
import xarray as xr

dynfolder = r'M:/_temp_modelling/modeloutput/'
icedyn_s1 = dynfolder + 'S1/ts.nc'
icedyn_s2 = dynfolder + 'S2/ts.nc'
icedyn_s3 = dynfolder + 'S3/ts.nc'
icedyn_s4 = dynfolder + 'S4/ts.nc'

smbfolder = r'M:/_temp_modelling/modeloutput/SMB_47300_10100_v2/'  # is also the outfolder
smb = smbfolder + 'out_table_annual_parameters_yrs_47300-10100_v2.xlsx'

# annual climate data table - - only needed for subplots 1 and 2
climate = r'M:/_temp_modelling/modelanalysis/dbdz/T_lastGlacial_kindler_et_al_2014_annual.xlsx'

# (Mean annual air temperature (K) at elevation z_stat, present-day value
T_zpcl_pd = 281.9

# Mean annual air temperature (K) at elevation z_stat and coldest phase (LGM) - only needed for subplots 1 and 2
T_zpcl_lgm = [271.9, 265.65, 267.9, 267.9]

# Factor Psi in: p_scale = exp(psi * delta_T_lgm), Huybrechts (2002)
psi = [0.0704, 0.0704, 0.0704, 0.028]

# correct for polar amplification (i.e. polar temperature variability > mid-latitudal)
T_climate_pd = -29  # (°C) Temperature in 'climate' that corresponds to present-day T at T_zpcl(delta_T_lgm == 0)
T_climate_lgm = -49  # (°C) Temperature in 'climate' that corresponds to LGM T at T_zpcl(maximal delta_T_lgm)

plots = [1, 2, 3, 4, 5]  # which plots to be plotted

# define time range of meteorological conditions to be shown in the plot
min_yr = 10100
max_yr = 47300

# specify xlims for plot
xlims = [10000, 47500]

# *********************************************************************************************************************
# prepare figure
fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(nrows=5, figsize=(15, len(plots)*3), constrained_layout=True)

fig.suptitle('ELA, glacier area and glacier volume', fontsize=14)

# read icedynamics data
ds1 = xr.open_dataset(icedyn_s1)
ds2 = xr.open_dataset(icedyn_s2)
ds3 = xr.open_dataset(icedyn_s3)
ds4 = xr.open_dataset(icedyn_s4)

# *********************************************************************************************************************
# preparations for plots 1 and 2: calculate Delta T and precipitation scaling

if 1 or 2 in plots:

    T_zpcl_lgm = np.array(T_zpcl_lgm)
    psi = np.array(psi)

    delta_T_lgm = T_zpcl_lgm - T_zpcl_pd
    clim = pd.read_excel(climate)
    clim.set_index('year', inplace=True)
    clim = clim.loc[min_yr:max_yr]  # slicing using loc includes both first and last element

    # t_years = year_start - year_end
    T_scale = delta_T_lgm / (T_climate_lgm - T_climate_pd)

    clim['T_s1'] = T_zpcl_lgm[0] + (clim['t(C)'] - T_climate_lgm) * T_scale[0] - 273.16
    clim['T_s2'] = T_zpcl_lgm[1] + (clim['t(C)'] - T_climate_lgm) * T_scale[1] - 273.16
    clim['T_s3'] = T_zpcl_lgm[2] + (clim['t(C)'] - T_climate_lgm) * T_scale[2] - 273.16
    clim['T_s4'] = T_zpcl_lgm[3] + (clim['t(C)'] - T_climate_lgm) * T_scale[3] - 273.16

    clim['DT_s1'] = clim['T_s1'] - T_zpcl_pd + 273.16
    clim['DT_s2'] = clim['T_s2'] - T_zpcl_pd + 273.16
    clim['DT_s3'] = clim['T_s3'] - T_zpcl_pd + 273.16
    clim['DT_s4'] = clim['T_s4'] - T_zpcl_pd + 273.16

    clim['p_scale_1'] = np.exp(psi[0] * clim['DT_s1'])
    clim['p_scale_2'] = np.exp(psi[1] * clim['DT_s2'])
    clim['p_scale_3'] = np.exp(psi[2] * clim['DT_s3'])
    clim['p_scale_4'] = np.exp(psi[3] * clim['DT_s4'])

# *********************************************************************************************************************
print('preparations done')
# *********************************************************************************************************************
# Plot 1: variation of ELA over time

if 1 in plots:

    ymin, ymax = -22, 0

    for item in ([ax0.title, ax0.xaxis.label, ax0.yaxis.label] +
                 ax0.get_xticklabels() + ax0.get_yticklabels()):
        item.set_fontsize(14)

    ax0.yaxis.grid(True, linestyle='--')  # Here we add the grid
    ax0.minorticks_on()  # use a simple method to turn minor ticks on and let matplotlib decide how to space them

    ax0.set_ylim(ymin, ymax)
    ax0.set_xlim(xlims[0], xlims[1])
    xmin, xmax = ax0.get_xlim()

    ax0.plot(clim.index, clim['DT_s1'], label='$\Delta T_{1}$', linewidth=1, color='r')
    ax0.plot(clim.index, clim['DT_s2'], label='$\Delta T_{2}$', linewidth=1, color='b')
    ax0.plot(clim.index, clim['DT_s3'], label='$\Delta T_{3}$', linewidth=1, color='orange')
    ax0.plot(clim.index, clim['DT_s4'], label='$\Delta T_{4}$', linewidth=1, color='g', linestyle=':')

    ax0.set_xlabel('Years BP')
    ax0.set_ylabel('$\Delta T$ (K)')

    mis1 = Rectangle((xmin, ymin), 14500 - xmin, ymax - ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    mis3 = Rectangle((29000, ymin), 28000, ymax - ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    # mis5 = Rectangle((71000,ymin), xmax-71000, ymax-ymin,
    #             linewidth=0, edgecolor='none', facecolor='black',
    #             alpha=0.07, zorder = 2)

    ax0.add_patch(mis1)
    ax0.add_patch(mis3)
    # ax.add_patch(mis5)

    # ax2.text(7500, ymax-(ymax-ymin)*0.95, 'MIS1')
    ax0.text(13000, ymax - (ymax - ymin) * 0.95, 'MIS1')
    ax0.text(19500, ymax - (ymax - ymin) * 0.95, 'MIS2')
    ax0.text(41000, ymax - (ymax - ymin) * 0.95, 'MIS3')
    # ax2.text(62000, ymax-(ymax-ymin)*0.95, 'MIS4')
    # ax2.text(97500, ymax-(ymax-ymin)*0.95, 'MIS5')

    ax0.legend(loc='lower left', bbox_to_anchor=(0.905, 0.05))  # for plotting ELA

    # ax2.invert_xaxis()
    print('subplot 1 done')

# *********************************************************************************************************************
# Plot 2: variation of ELA over time

if 2 in plots:

    ymin, ymax = 0, 1

    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
                 ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(14)

    ax1.yaxis.grid(True, linestyle='--')  # Here we add the grid
    ax1.minorticks_on()  # use a simple method to turn minor ticks on and let matplotlib decide how to space them

    ax1.set_ylim(ymin, ymax)
    ax1.set_xlim(xlims[0], xlims[1])
    xmin, xmax = ax1.get_xlim()

    ax1.plot(clim.index, clim['p_scale_1'], label='$S_{pa}$ 1', linewidth=1, color='r')
    ax1.plot(clim.index, clim['p_scale_2'], label='$S_{pa}$ 2', linewidth=1, color='b')
    ax1.plot(clim.index, clim['p_scale_3'], label='$S_{pa}$ 3', linewidth=1, color='orange')
    ax1.plot(clim.index, clim['p_scale_4'], label='$S_{pa}$ 4', linewidth=1, color='g')

    ax1.set_xlabel('Years BP')
    ax1.set_ylabel('Scaling factor P (-)')

    mis1 = Rectangle((xmin, ymin), 14500 - xmin, ymax - ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    mis3 = Rectangle((29000, ymin), 28000, ymax - ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    # mis5 = Rectangle((71000,ymin), xmax-71000, ymax-ymin,
    #             linewidth=0, edgecolor='none', facecolor='black',
    #             alpha=0.07, zorder = 2)

    ax1.add_patch(mis1)
    ax1.add_patch(mis3)
    # ax.add_patch(mis5)

    # ax2.text(7500, ymax-(ymax-ymin)*0.95, 'MIS1')
    ax1.text(13000, ymax - (ymax - ymin) * 0.95, 'MIS1')
    ax1.text(19500, ymax - (ymax - ymin) * 0.95, 'MIS2')
    ax1.text(41000, ymax - (ymax - ymin) * 0.95, 'MIS3')
    # ax2.text(62000, ymax-(ymax-ymin)*0.95, 'MIS4')
    # ax2.text(97500, ymax-(ymax-ymin)*0.95, 'MIS5')

    ax1.legend(loc='lower left', bbox_to_anchor=(0.905, 0.05))  # for plotting ELA

    # ax2.invert_xaxis()
    print('subplot 2 done')

# *********************************************************************************************************************
# Plot 3: variation of ELA over time

if 3 in plots:

    df1 = pd.read_excel(smb, sheet_name='Scenario_1')
    df2 = pd.read_excel(smb, sheet_name='Scenario_2')
    df3 = pd.read_excel(smb, sheet_name='Scenario_3')
    df4 = pd.read_excel(smb, sheet_name='Scenario_4')

    ymin, ymax = 750, 2500

    for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
                 ax2.get_xticklabels() + ax2.get_yticklabels()):
        item.set_fontsize(14)

    ax2.yaxis.grid(True, linestyle='--')  # Here we add the grid
    ax2.minorticks_on()  # use a simple method to turn minor ticks on and let matplotlib decide how to space them

    ax2.set_ylim(ymin, ymax)
    ax2.set_xlim(xlims[0], xlims[1])
    xmin, xmax = ax2.get_xlim()

    ax2.plot(df1['year'], df1['ELA (m a.s.l.)'], label='ELA$_1$', linewidth=1, color='r')
    ax2.plot(df2['year'], df2['ELA (m a.s.l.)'], label='ELA$_{2}$', linewidth=1, color='b')
    ax2.plot(df3['year'], df3['ELA (m a.s.l.)'], label='ELA$_{3}$', linewidth=1, color='orange')
    ax2.plot(df4['year'], df4['ELA (m a.s.l.)'], label='ELA$_{4}$', linewidth=1, color='g')

    ax2.set_xlabel('Years BP')
    ax2.set_ylabel('ELA (m a.s.l.)')

    mis1 = Rectangle((xmin, ymin), 14500-xmin, ymax-ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    mis3 = Rectangle((29000, ymin), 28000, ymax-ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    # mis5 = Rectangle((71000,ymin), xmax-71000, ymax-ymin,
    #             linewidth=0, edgecolor='none', facecolor='black',
    #             alpha=0.07, zorder = 2)

    ax2.add_patch(mis1)
    ax2.add_patch(mis3)
    # ax.add_patch(mis5)

    # ax2.text(7500, ymax-(ymax-ymin)*0.95, 'MIS1')
    ax2.text(13000, ymax-(ymax-ymin)*0.95, 'MIS1')
    ax2.text(19500, ymax-(ymax-ymin)*0.95, 'MIS2')
    ax2.text(41000, ymax-(ymax-ymin)*0.95, 'MIS3')
    # ax2.text(62000, ymax-(ymax-ymin)*0.95, 'MIS4')
    # ax2.text(97500, ymax-(ymax-ymin)*0.95, 'MIS5')

    ax2.legend(loc='lower left', bbox_to_anchor=(0.905, 0.05))  # for plotting ELA

    # ax2.invert_xaxis()
    print('subplot 3 done')

# *********************************************************************************************************************
# Plot 4: variation of glacier length over time

if 4 in plots:

    for item in ([ax3.title, ax3.xaxis.label, ax3.yaxis.label] +
                 ax3.get_xticklabels() + ax3.get_yticklabels()):
        item.set_fontsize(14)

    ax3.yaxis.grid(True, linestyle='--')  # Here we add the grid
    ax3.minorticks_on()  # use a simple method to turn minor ticks on and let matplotlib decide how to space them

    ax3.set_xlim(xlims[0], xlims[1])
    xmin, xmax = ax3.get_xlim()

    ax3.plot(ds1.time * (-1), ds1.area, label='Area$_1$', linewidth=1, color='r')
    ax3.plot(ds2.time * (-1), ds2.area, label='Area$_{2}$', linewidth=1, color='b')
    ax3.plot(ds3.time * (-1), ds3.area, label='Area$_{3}$', linewidth=1, color='orange')
    ax3.plot(ds4.time * (-1), ds4.area, label='Area$_{4}$', linewidth=1, color='g')

    ax3.set_xlabel('Years BP')
    ax3.set_ylabel('Area (km$^2$)')

    mis1 = Rectangle((xmin, ymin), 14500 - xmin, ymax - ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    mis3 = Rectangle((29000, ymin), 28000, ymax - ymin,
                     linewidth=0, edgecolor='none', facecolor='black',
                     alpha=0.07, zorder=2)

    # mis5 = Rectangle((71000,ymin), xmax-71000, ymax-ymin,
    #             linewidth=0, edgecolor='none', facecolor='black',
    #             alpha=0.07, zorder = 2)

    ax3.add_patch(mis1)
    ax3.add_patch(mis3)
    # ax3.add_patch(mis5)

    # ax3.text(7500, ymax-(ymax-ymin)*0.95, 'MIS1')
    ax3.text(13000, ymax - (ymax - ymin) * 0.95, 'MIS1')
    ax3.text(19500, ymax - (ymax - ymin) * 0.95, 'MIS2')
    ax3.text(41000, ymax - (ymax - ymin) * 0.95, 'MIS3')
    # ax3.text(62000, ymax-(ymax-ymin)*0.95, 'MIS4')
    # ax3.text(97500, ymax-(ymax-ymin)*0.95, 'MIS5')

    ax3.legend(loc='lower left', bbox_to_anchor=(0.905, 0.5))

    # ax3.invert_xaxis()
    print('subplot 4 done')

# *********************************************************************************************************************
# Plot 5: variation of glacier volume over time

if 5 in plots:

    for item in ([ax4.title, ax4.xaxis.label, ax4.yaxis.label] +
                 ax4.get_xticklabels() + ax4.get_yticklabels()):
        item.set_fontsize(14)

    ax4.yaxis.grid(True, linestyle='--')  # Here we add the grid
    ax4.minorticks_on()  # use a simple method to turn minor ticks on and let matplotlib decide how to space them

    ax4.set_xlim(xlims[0], xlims[1])
    xmin, xmax = ax4.get_xlim()

    ax4.plot(ds1.time * (-1), ds1.vol, label='Volume$_1$', linewidth=1, color='r')
    ax4.plot(ds2.time * (-1), ds2.vol, label='Volume$_{2}$', linewidth=1, color='b')
    ax4.plot(ds3.time * (-1), ds3.vol, label='Volume$_{3}$', linewidth=1, color='orange')
    ax4.plot(ds4.time * (-1), ds4.vol, label='Volume$_{4}$', linewidth=1, color='g')

    ax4.set_xlabel('Years BP')
    ax4.set_ylabel('Volume (km$^3$)')

    # mis1 = Rectangle((xmin, ymin), 14500 - xmin, ymax - ymin,
    #                  linewidth=0, edgecolor='none', facecolor='black',
    #                  alpha=0.07, zorder=2)

    ax4.axvspan(xmin, 14500, color='tab:gray', alpha=0.1, zorder=2)

    # mis3 = Rectangle((29000, ymin), 28000, ymax - ymin,
    #                  linewidth=0, edgecolor='none', facecolor='black',
    #                  alpha=0.07, zorder=2)

    ax4.axvspan(29000, xmax, color='tab:gray', alpha=0.1, zorder=2)

    # mis5 = Rectangle((71000,ymin), xmax-71000, ymax-ymin,
    #             linewidth=0, edgecolor='none', facecolor='black',
    #             alpha=0.07, zorder = 2)

    # ax4.add_patch(mis1)
    # ax4.add_patch(mis3)
    # ax.add_patch(mis5)

    # ax4.text(7500, ymax-(ymax-ymin)*0.95, 'MIS1')
    ax4.text(13000, ymax - (ymax - ymin) * 0.95, 'MIS1')
    ax4.text(19500, ymax - (ymax - ymin) * 0.95, 'MIS2')
    ax4.text(41000, ymax - (ymax - ymin) * 0.95, 'MIS3')
    # ax4.text(62000, ymax-(ymax-ymin)*0.95, 'MIS4')
    # ax4.text(97500, ymax-(ymax-ymin)*0.95, 'MIS5')

    ax4.legend(loc='lower left', bbox_to_anchor=(0.905, 0.5))  # for plotting ELA

    # ax4.invert_xaxis()
    print('subplot 5 done')

# *********************************************************************************************************************
# Finalize plot

out_fig = smbfolder + 'GlacierDyn_time_series_comparison' + str(int(max_yr)) + '-' + str(int(min_yr)) + '.pdf'

plt.tight_layout()

plt.savefig(out_fig)
print('*** all done ***')
