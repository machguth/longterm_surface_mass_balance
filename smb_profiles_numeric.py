"""
    *** Numerical calculation of surface mass balance profiles using a simple degree-day approach. ***

Purpose of the calculation is to investigate the influence of different annual temperature amplitudes on the
surface mass balance gradient in the abaltion area.

Eventually, the surface mass balance should be used to drive an ice flow model and thereby to study the influence
of changing surface mass balance gradients on advance and reatreat behaviour of glaciers during the last glacial period

The profiles are calculated for two different annual amplitudes of surface temperature. In addition, both simulations
can be perturbed by an offset in air temperature and precipitation (totalling four scenarios if this option is used).

*** Execution of the code is partially multi threaded. Multi threading is used if overhead from multi threading is
smaller than time gain due to parallel computing (the threshold lies at around a total of 100 mass balance years). To
disable parallel computations, just increase 'parallel_ts' (in years) to a very large value. ***
"""

#import sys
#sys.path.insert(0, r'C:/Users/machguth/switchdrive/src/py/')

import os
import time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing as mp

import smb_profiles_numeric_funcs as smbf
import smb_profiles_numeric_plots as smbp
import smb_profiles_numeric_tables as smbt

if __name__=='__main__': # required to use parallel computation under Windows

    start = time.time()
    ##########################################  input, output, parameter definitions  #####################################
    # --------------------------
    T_zpcl_pd = 281.9  # (K) Mean annual air temperature at elevation z_stat, present-day value
    TAa_pd = 8.75  # (K) Annual amplitude of air temperature, present-day value
    #TAa = [8.75, 15]  # (K) Air temperature, annual amplitude, LGM conditions
    TAa = [8.75, 15, 15, 15]  # (K) Air temperature, annual amplitude, LGM conditions
    #TAa = [15, 15]  # (K) Air temperature, annual amplitude, LGM conditions
    Tsd = [3.5] * 4  # (K) Standard deviation daily air temperature ==> set to zero to suppress daily cycle of T
    Tg = [0.006] * 4  # (K m-1) Temperature lapse rate
    Ts = 1  # (°C) threshold temperature snowfall and rain

    T_off = [0] * 4  # (K) offset for simulating impact of climate change
    p_off = [0] * 4 # (m yr-1) offset for simulating impact of climate change

    p_a = [0.0002857] * 4 # () factor a in: p = a*z + b, where z is elevation and p is present-day annual precip
    p_b = [1.1412] * 4 # (m) factor b in: p = a*z + b, where z is elevation and p is present-day annual precip
    #psi = [0.028, 0.028]  # () psi in: p_scale = exp(psi * delta_T_lgm), Huybrechts (2002); set 0 for p_scale = 1.
    #psi = [0.0704, 0.028]  # () psi in: p_scale = exp(psi * delta_T_lgm), Huybrechts (2002); set 0 for p_scale = 1.
    psi = [0.0704, 0.0704, 0.0704, 0.028]  # () psi in: p_scale = exp(psi * delta_T_lgm), Huybrechts (2002); set 0 for p_scale = 1.

    # Calculating refreezing of meltwater. Two parameterizations can be used, either 'Reeh91' or 'Pfeffer91-Reeh05'
    # 'Reeh91' uses fixed fraction of snow w.e. to calculate amount of meltwater that refreezes before runoff starts
    # 'Pfeffer91-Reeh05' calculates fraction of snow w.e. based on firn T (using MAAT)
    refreeze = [True] * 4 # specify whether refreezing of meltwater is simulated
    #refreeze = [False, True]

    #refreeze_parameterization = ['Reeh91'] * 2 # parameterization to be used
    refreeze_parameterization = ['Pfeffer91-Reeh05'] * 4 # parameterization to be used

    pddf = [3.297, 6, 8.791]  # (mm K-1 d-1) degree day factors, first for snow, second for firn, third for ice

    zmm = [3250, 50]    # (m a.s.l.) max and min elevation z for calculation, identical to "hypsometry"
    step = 100          # (m) elevation step for calculation of B, identical to "hypsometry"

    z_paleoclim = 400   # (m a.s.l.) elevation for which quantitative paleoclimate dTAa exist

    t_start = 274       # (julian days) always start calculation October 1 (day of year 274)

    # Table of annual air temperature. Comment out to trigger linear run over a duration of "t_years"
    #climate = r'C:/Horst/modeling/modelanalysis/dbdz/T_lastGlacial_kindler_et_al_2014_annual.xlsx'
    climate = r'M:/_temp_modelling/modelanalysis/dbdz/T_lastGlacial_kindler_et_al_2014_annual.xlsx'

    # Note that delta_T_lgm (MAAT difference between present-day and LGM (coldest) conditions)
    # is not explicitly defined. Instead, it is specified implicitly by chosing values for T_zpcl_pd and
    # T_zpcl_lgm.

    # ---------------------------------------- Specify paleoclimate data --------------------------------------------
    t_years = 100  # (years) number of years to calculate - used if 'climate' table is not given.

    # Start and end year BP of model run. Model runs from higher to lower year numbers. Used if climate table.
    #year_start = 122900  # (years BP)
    year_start = 47300 # (years BP)
    #year_start = 23820  # (years BP)
    year_end = 10100  # (years BP)
    #year_end = 23720  # (years BP)

    # Mean annual air temperature at elevation z_stat and coldest phase (LGM), in Kelvin
    #T_zpcl_lgm = [270.9, 264.65]  # (K)
    #T_zpcl_lgm = [266.9, 266.9]  # (K)
    #T_zpcl_lgm = [270.9, 264.65, 266.9, 266.9]  # (K)
    T_zpcl_lgm = [271.9, 265.65, 267.9, 267.9]  # (K)

    # correct for polar amplification (i.e. polar temperature variability > mid-latitudinal T variability)
    T_climate_pd = -29  # (°C) Temperature in 'climate' that corresponds to present-day T at T_zpcl (delta_T_lgm == 0)
    T_climate_lgm = -49  # (°C) Temperature in 'climate' that corresponds to LGM T at T_zpcl(maximal delta_T_lgm)

    # ----------------------------------------- Specify hypsometry -------------------------------------------------
    #hypsometry = r'C:/Horst/modeling/modelanalysis/dbdz/s1_hypsometry_mod.xlsx'
    #hypsometry = r'C:/Users/Horst/switchdrive/_temp_modelling/modelanalysis/dbdz/s1_hypsometry_mod.xlsx'
    hypsometry = r'M:/_temp_modelling/modelanalysis/dbdz/s1_hypsometry_mod_v2.xlsx'
    #hypsometry = r'C:/Horst/modeling/modelanalysis/dbdz/s2_hypsometry.xlsx'

    parallel_ts = 100. # (mass balance years) threshold for use of parallel computing (recommended around 100).

    # ---------------------------------------------   Output    ----------------------------------------------------
    #outfolder = r'C:/Horst/modeling/modelanalysis/dbdz/'
    #outfolder = r'C:/Users/Horst/switchdrive/_temp_modelling/modeloutput/dbdz/'
    #outfolder = r'C:/Users/machguth/switchdrive/_temp_modelling/modeloutput/test/'
    outfolder = r'M:/_temp_modelling/modeloutput/SMB_47300_10100_v2/'
    # --------------------------

    # #############################################  preparations  ####################################################
    # check if output folder exists, if no create
    isdir = os.path.isdir(outfolder)
    if not isdir:
        os.mkdir(outfolder)

    # check if climate table is given. If not, create a dummy variable
    if 'climate' not in locals():
        climate = 'none'

    # Initiate the climate data arrays
    T_raw, T_zpcl, TAa, delta_T, t_years, year_end, year_start = smbf.make_clim_arrs(climate, T_zpcl_lgm,
                                                                                     TAa, TAa_pd, year_end, year_start,
                                                                                     T_zpcl_pd, T_climate_lgm,
                                                                                     T_climate_pd, t_years)

    # read the hypsometry
    df_hypso = pd.read_excel(hypsometry)  # , index_col = 'elevation')
    df_hypso['area'] = np.flip(df_hypso['area'].values, axis=0)  # flip the order of the entries to agree with order of dem and pdds arrays
    df_hypso['elevation'] = np.flip(df_hypso['elevation'].values, axis=0)  # flip the order of the entries to agree with order of dem and pdds arrays
    # --------------------------
    dem = np.linspace(zmm[0], zmm[1], int((zmm[0] - zmm[1])/step+1))

    pddf = np.array(pddf) # convert to numpy array

    # establish an array of temperature at sea level for each of the scenarios
    T0m = [0] * len(T_zpcl)
    for ind, i in enumerate(T_zpcl):
        T0m[ind] = i - 273.15 + z_paleoclim * Tg[ind]

    df_out = pd.DataFrame(columns=['z (m a.s.l.)', 'A (km2)', 's_Tpos (jul.d)', 'e_Tpos (jul.d)', 'dur_Tpos (d)',
                                   'Tmax (°C)', 'Tmin (°C)', 'pdd (°C d)', 'pdd_eff (°C d)',
                                   'b (m w.e)', 'Ba (m3 w.e.)', 'Retention (m w.e.)', 'r_max ()', 'rho_firn (kg m-3)'])

    df_out['z (m a.s.l.)'] = dem
    df_out['A (km2)'] = df_hypso['area']

    # table of ending dates and duration of months
    months = np.empty((2,12), dtype=np.int32)
    months[0, :] = [30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364]  # end julian day numbers of all 12 months
    months[1, :] = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    # set up a few output arrays
    t_arr = np.arange(t_start, t_years*365+t_start)
    Ba = np.zeros((len(T_zpcl), len(dem), t_years))
    param_arr = np.zeros((len(T_zpcl), 17, t_years))
    Bm = [] # table of monthly mass balance at all grid cells
    Tm = [] # table of monthly air temperature at all grid cells
    pdds = np.zeros((len(T_zpcl), 13, len(dem), t_years))
    clim_info = np.zeros((len(T_zpcl), 4))
    for ind, i in enumerate(T_zpcl):
        pdds[ind, 0, :, 0] = dem

    # ************* Overview of array content *************
    # pdds[case#, var#, dem, year_mb];
    #        variables = [start_abl, end_abl, dur_abl, maxT, minT, pdd_sum, pdd_sum_eff, B, BB];
    #        per DEM-gridcell and annual
    # Ba[case#, dem, year_mb]; variable = [annual mass balance] per DEM-gridcell
    # Bm[case#, months, dem]; variable = monthly mass balance per DEM-gridcell
    # Tm[case#, months, dem]; variable = monthly air temperature per DEM-gridcell
    # param_arr[case#, var#, year_mb]; variables = [B_glw, maxB, minB, dbdz, ELA, AAR]; glacier wide and annual

    ##############################################   computations   ###################################################

    # display some basic informations on the chosen climate
    for ind, i in enumerate(T_zpcl):
        clim_info[ind, :] = smbf.climate_info(T_zpcl[ind][0], TAa[ind][0], Tg[ind], delta_T[ind][0],
                         p_a[ind], p_b[ind], psi[ind], z_paleoclim, z_paleoclim, ind, T_off[ind], p_off[ind])

    # decide whether parallel processing saves time or not
    if len(T_zpcl) * t_years < parallel_ts:

        # numerically calculate surface mass balance for all DEM grid cells (= elevation classes)
        for ind, i in enumerate(T_zpcl):
            pdds[ind, 1:, :, :], Ba[ind, :, :], Bmt, param_arr[ind, :, :], mb_years, Tmt = \
                smbf.smb_numerical(dem, t_arr, T0m[ind], Tg[ind],
                                   TAa[ind], delta_T[ind], p_a[ind], p_b[ind], psi[ind], Ts,
                                   pddf, Tsd[ind], df_hypso, step, t_years, months, z_paleoclim,
                                   T_off[ind], p_off[ind], refreeze[ind], refreeze_parameterization[ind])

            Bm.append(Bmt)
            Tm.append(Tmt)
            print('finished numerical calculation #' + str(ind) + ' of a total of ' + str(
                len(T_zpcl)) + ' SMB calculations.')

    else:

        # numerically calculate surface mass balance for all DEM grid cells (= elevation classes)
        num_cores = mp.cpu_count()
        if num_cores > len(T_zpcl):
            num_cores = len(T_zpcl)
        print('\n *** Using {:}'.format(num_cores) + ' cores in parallel processing. ***')
        print('')

        results = Parallel(n_jobs=num_cores)(delayed(smbf.smb_numerical)(dem, t_arr, T0m[ind], Tg[ind],
                                             TAa[ind], delta_T[ind], p_a[ind], p_b[ind], psi[ind], Ts, pddf,
                                             Tsd[ind], df_hypso, step, t_years, months, z_paleoclim,
                                             T_off[ind], p_off[ind], refreeze[ind], refreeze_parameterization[ind])
                                             for ind, i in enumerate(T_zpcl))

        a1, a2, a3, a4, a5, a6 = zip(*results)

        for ind, i in enumerate(T_zpcl):
            pdds[ind, 1:, :, :] = a1[ind]
            Ba[ind, :, :] = a2[ind]
            Bm.append(a3[ind])
            param_arr[ind, :, :] = a4[ind]
            mb_years = a5[ind]
            Tm.append(a6[ind])

    # clean-up: check how many hydrological years there were and if less than t_years, remove superfluous years
    if mb_years < t_years:
        Ba = Ba[:, :, :mb_years]
        param_arr = param_arr[:, :, :mb_years]
        pdds = pdds[:, :, :, :mb_years]

    print('finished '+ str(len(T_zpcl)) +' numerical calculations.')

    # convert monthly smb list to numpy array
    Bm = np.array(Bm)
    Tm = np.array(Tm)

    # Modify Ba for the purpose of writing output, also add time axis
    t_axis = np.arange(year_start-1, year_end, -1)  # -1 as first year is not a complete SMB year (autumn start)
    Ba = Ba.transpose((0,2,1))
    Ba = np.insert(Ba, 0, t_axis, axis=2)

    print('---')
    print('diff. B glacier-wide (of first two scenarios, last year): %s' % (param_arr[0, 0, -1] - param_arr[1, 0, -1]))
    print('---')

    for ind, i in enumerate(T_zpcl):
        smbt.write_output_tables_csv_vertical(outfolder, Bm[ind, :, :], ind, dem, 'monthly', 'mass-balance')
        smbt.write_output_tables_csv_vertical(outfolder, Ba[ind, :, :], ind, dem, 'annual', 'mass-balance')

    for ind, i in enumerate(T_zpcl):
        smbt.write_output_tables_csv_vertical(outfolder, Tm[ind, :, :], ind, dem, 'monthly', 'air-temperature')

    # Plotting - three different plots are created
    smbp.plot_fig_1(outfolder, T_zpcl, dem, pdds[:, :, :, -1])
    smbp.plot_fig_2(outfolder, pdds[:, 3, :, -1], pdds[:, 4, :, -1], pdds[:, :, :, -1], TAa)
    smbp.plot_fig_3(clim_info, TAa, Tsd, T_zpcl, Tg, p_a, pddf, T_off, p_off, outfolder,
                   param_arr[:, 1, -1], param_arr[:, 2, -1], zmm, dem, Ba[:, -1, 1:], param_arr[:, 0, -1],
                    df_hypso, param_arr[:, 3, -1], param_arr[:, 4, -1], param_arr[:, 5, -1], z_paleoclim)

    # Write output tables
    smbt.write_output_table_xls1(outfolder, T_zpcl, df_out, pdds[:, :, :, -1])
    smbt.write_output_table_xls2(outfolder, param_arr, T_raw, T_zpcl, year_start, year_end, mb_years)

    print('Finished')

    end = time.time()
    print('')
    print('Process took {:4f}'.format(end-start) + ' seconds.')
