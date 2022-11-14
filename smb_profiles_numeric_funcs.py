"""
    *** Numerical calculation of surface mass balance profiles using a simple degree-day approach. ***

Purpose of the calculation is to investigate the influence of different annual temperature amplitudes on the
surface mass balance gradient in the abaltion area.

Eventually, the surface mass balance should be used to drive an ice flow model and thereby to study the influence
of changing surface mass balance gradients on advance and reatreat behaviour of glaciers during the last glacial period

The profiles are calculated for two different annual amplitudes of surface temperature. In addition, both simulations
can be perturbed by an offset in air temperature and precipitation (totalling four scenarios if this option is used).
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import math
import numpy as np
import csv
import pandas as pd
from scipy import stats
import scipy.special as sp


# ###################################################   functions   ###################################################
def make_clim_arrs(climate, T_zpcl_lgm, TAa, TAa_pd, year_end, year_start, T_zpcl_pd,
                   T_climate_lgm, T_climate_pd, t_years):

    T_zpcl_lgm = np.array(T_zpcl_lgm)
    delta_T_lgm = T_zpcl_lgm - T_zpcl_pd
    n = len(T_zpcl_lgm)  # number of scenarios

    if climate == 'none':  # Check whether an array of annual air temperatures is given
        # # below +1 year as, e.g., calcul. for one year (365 days) stretches over yr. 0 and yr. 1 (two years)
        # temp = np.zeros((2, t_years + 1))
        # delta_T = np.zeros((2, t_years + 1))
        # delta_T[0, :], delta_T[1, :] = delta_T_lgm[0], delta_T_lgm[1]
        # temp[0, :], temp[1, :] = T_zpcl_lgm[0], T_zpcl_lgm[1]
        # T_zpcl = temp
        # clim_rel = np.full(len(temp[0, :]), np.NAN)
        # year_end, year_start = 0, t_years
        #
        # # no need to modify seasonality if no climatology specified (in this case TAa is constant)
        # TAa = [np.full(t_years + 1, TAa[0]), np.full(t_years + 1, TAa[1])]

        # below +1 year as, e.g., calcul. for one year (365 days) stretches over yr. 0 and yr. 1 (two years)
        temp = np.zeros((n, t_years + 1))
        delta_T = np.zeros((n, t_years + 1))

        for i in np.arange(n):

            # THIS NEEDS MODIFICATION
            delta_T[0, :], delta_T[1, :] = delta_T_lgm[0], delta_T_lgm[1]
            temp[0, :], temp[1, :] = T_zpcl_lgm[0], T_zpcl_lgm[1]
            T_zpcl = temp
            clim_rel = np.full(len(temp[0, :]), np.NAN)
            year_end, year_start = 0, t_years

            # no need to modify seasonality if no climatology specified (in this case TAa is constant)
            TAa = [np.full(t_years + 1, TAa[0]), np.full(t_years + 1, TAa[1])]

    else:
        # clim = pd.read_excel(climate)
        # clim.set_index('year', inplace=True)
        # clim_rel = clim.loc[year_end:year_start, 't(C)']  # slicing using loc includes both first and last element
        # clim_rel = clim_rel[
        #        ::-1].values  # + 273.16 # reverse order of elements in the array, convert to numpy array
        #
        # t_years = year_start - year_end
        # T_scale = [delta_T_lgm[0] / (T_climate_lgm - T_climate_pd), delta_T_lgm[1] / (T_climate_lgm - T_climate_pd)]
        # T_zpcl = [T_zpcl_lgm[0] + (clim_rel - T_climate_lgm) * T_scale[0],
        #           T_zpcl_lgm[1] + (clim_rel - T_climate_lgm) * T_scale[1]]
        #
        # delta_T = [(clim_rel - T_climate_pd) * T_scale[0], (clim_rel - T_climate_pd) * T_scale[1]]
        #
        # TAa = [(delta_T_lgm[0] - (T_zpcl_lgm[0] - T_zpcl[0])) / delta_T_lgm[0] * (TAa[0] - TAa_pd) + TAa_pd,
        #        (delta_T_lgm[1] - (T_zpcl_lgm[1] - T_zpcl[1])) / delta_T_lgm[1] * (TAa[1] - TAa_pd) + TAa_pd]

        clim = pd.read_excel(climate)
        clim.set_index('year', inplace=True)
        clim_rel = clim.loc[year_end:year_start, 't(C)']  # slicing using loc includes both first and last element
        clim_rel = clim_rel[
               ::-1].values  # + 273.16 # reverse order of elements in the array, convert to numpy array

        t_years = year_start - year_end

        for i in np.arange(n):

            # THIS NEEDS MODIFICATION
            T_scale = [delta_T_lgm[0] / (T_climate_lgm - T_climate_pd), delta_T_lgm[1] / (T_climate_lgm - T_climate_pd)]
            T_zpcl = [T_zpcl_lgm[0] + (clim_rel - T_climate_lgm) * T_scale[0],
                      T_zpcl_lgm[1] + (clim_rel - T_climate_lgm) * T_scale[1]]

            delta_T = [(clim_rel - T_climate_pd) * T_scale[0], (clim_rel - T_climate_pd) * T_scale[1]]

            TAa = [(delta_T_lgm[0] - (T_zpcl_lgm[0] - T_zpcl[0])) / delta_T_lgm[0] * (TAa[0] - TAa_pd) + TAa_pd,
                   (delta_T_lgm[1] - (T_zpcl_lgm[1] - T_zpcl[1])) / delta_T_lgm[1] * (TAa[1] - TAa_pd) + TAa_pd]

    return clim_rel, T_zpcl, TAa, delta_T, t_years, year_end, year_start


# *********************************************************************************************************************
# calculate annual means/sums of T and P at elevation z
def climate_info(T_z1, TAa, Tg, delta_T, p_a, p_b, psi, z1, z2, ind, T_off, p_off):
    MAATz2 = T_z1 - (z2 - z1) * Tg + T_off
    Tsz2 = MAATz2 + TAa
    Twz2 = MAATz2 - TAa
    Pz2 = (p_a * z2 + p_b) * np.exp(psi * delta_T) + p_off

    if ind != -1:
        scenario = ind + 1
        print('\n Scenario {:}'.format(scenario) + ' at {:}'.format(int(z2)) + ' m a.s.l., at start of time axis:')
        print('MAAT: {:.2f}'.format(MAATz2 - 273.15) + ' °C')
        print('Tmax: {:.2f}'.format(Tsz2 - 273.15) + ' °C')
        print('Tmin: {:.2f}'.format(Twz2 - 273.15) + ' °C')
        print('P: {:.2f}'.format(Pz2) + ' m/yr')
        print('Delta T: {:.2f}'.format(delta_T) + ' °C')

    return MAATz2, Tsz2, Twz2, Pz2


# *********************************************************************************************************************
# the actual calculation of the surface mass balance
def smb_numerical(dem, t_arr, T0m, Tg, TAa, delta_T, p_a, p_b, psi, Ts, pddf,
                  Tsd, df_hypso, step, t_years, months, z_paleoclim, T_off, p_off,
                  refreeze, refreeze_parameterization):

    pddf = pddf/1000.  # convert PDDFs from mm d-1 °C-1 to m d-1 °C-1

    pdds = np.zeros((9, len(dem), t_years))

    year = 0
    year_mb = 0
    snow = np.zeros_like(dem)
    firn = np.zeros_like(dem)
    ice = np.zeros_like(dem)
    snow_1 = np.zeros_like(dem)
    snow_2 = np.zeros_like(dem)
    smin = np.zeros_like(dem)

    la = 500  # maximum number of days to be kept in arrays b_arr and bt_arr
    b_arr = np.zeros((la, len(dem), 7))  # daily array - includes snow, firn, ice, bt, Ttz, PDD_simple, PDD_eff
    bt_arr = np.zeros((la, len(dem)))  # array for daily cumulative mass balances
    Ba_arr = np.zeros((len(dem), t_years))  # array for annual mass balance profiles
    param_arr = np.zeros((16, t_years))  # array for annual and glacier-wide smb parameters
    Bm_arr = []  # array for monthly mass balance per grid-cell
    Tm_arr = []  # array for monthly air temperature per grid-cell

    t_arr_jul = np.zeros(la)

    for num, day_abs in enumerate(t_arr):  # FOR num = jstart_c, jend_c DO BEGIN

        # calculate day of the year (1-365 or 366)
        day = day_abs - year * 365

        # write to array of julian day numbers (from 1 to 365)
        t_arr_jul = np.roll(t_arr_jul, 1)
        t_arr_jul[0] = day

        if num == 0:  # initialize refreezing
            r_max = refreezing(refreeze, refreeze_parameterization, num, T0m, Tg, dem)

        mn = np.where(months[0,:] == day)
        if len(mn[0]) == 1 and num > 31:  # test if a month has just ended (and if enough total days have passed)
            Bm_arr.append(daily_to_monthly_v2(bt_arr[:, :], year, t_arr_jul, months, mn[0], 'sum'))
            Tm_arr.append(daily_to_monthly_v2(b_arr[:, :, 4], year, t_arr_jul, months, mn[0], 'mean'))

        if day_abs/365 == np.floor(day_abs/365):  # To be carried out every end of a year
            year += 1
            r_max = refreezing(refreeze, refreeze_parameterization, num, T0m, Tg, dem)  # update refreezing

        if day_abs/365 == np.floor(day_abs/365) and num > 365:  # only to be executed after a full year has passed
            Ba_arr[:, year_mb] = b_arr[0, :, 3] - b_arr[364, :, 3]  # diff. between 1. and last day of the last 365 days
            pdds[:, :, year_mb] = smb_param_dem(t_arr_jul, b_arr, dem, df_hypso)
            param_arr[0, year_mb], param_arr[1, year_mb], param_arr[2, year_mb] = \
                glacier_wide_smb(Ba_arr[:, year_mb], df_hypso)
            param_arr[3, year_mb], param_arr[4, year_mb], param_arr[5, year_mb] = \
                smb_param_glacier_wide(Ba_arr[:, year_mb], dem, df_hypso, step)
            param_arr[6, year_mb] = mass_input(dem, param_arr[4, year_mb], pdds[8, :, year_mb])
            param_arr[7, year_mb], param_arr[8, year_mb], param_arr[9, year_mb], param_arr[10, year_mb] = \
                climate_info(T0m[year], TAa[year], Tg, delta_T[year],
                             p_a, p_b, psi, 0, z_paleoclim, -1, T_off, p_off)
            param_arr[11, year_mb], param_arr[12, year_mb], param_arr[13, year_mb], param_arr[14, year_mb] = \
                climate_info(T0m[year], TAa[year], Tg, delta_T[year],
                             p_a, p_b, psi, 0, param_arr[4, year_mb], -1, T_off, p_off)
            param_arr[15, year_mb] = TAa[year]

            year_mb += 1

        # ********************************** temperature *************************************************
        # Calculate temperature distribution over the entire array
        Ttz = T0m[year] - dem * Tg + TAa[year] * (- np.cos((2. * math.pi * (day - 30.)) / 365.)) + T_off

        # Use Tsd to calculate effective temperature
        # compute positive part of temperature everywhere
        positivepart = np.greater(Ttz, 0)*Ttz

        # compute Calov and Greve (2005) integrand, ignoring division by zero (copied from code by J. Seguinot)
        with np.errstate(divide='ignore', invalid='ignore'):
            normTtz = Ttz / (np.sqrt(2)*Tsd)
        calovgreve = (Tsd/np.sqrt(2*np.pi)*np.exp(-normTtz**2) +
                      Ttz/2*sp.erfc(-normTtz))

        # use positive part if sigma == zero, else use Calov and Greve (copied from code by J. Seguinot)
        Ttz_eff = np.where(Tsd == 0., positivepart, calovgreve)

        # ********************************** precipitation ************************************************
        # Calculate precipitation at elevation z and for 5-days period
        # accumulation if air temperature lower Ts (in °C)

        # annual precip over dem, including scaling as f(Delta_T)
        p_dem = (p_a * dem + p_b) * np.exp(psi * delta_T[year]) + p_off

        if np.floor(num / 5.0) == num / 5.0: # precipitation every fifth day
            ptz = (p_dem / 365.) * 5
            acc = (Ttz <= Ts) * ptz
        else:
            acc = 0.

        # -------- accumulation  ----------
        snow = snow + acc

        ls = Ttz_eff * pddf[0]

        pot_smelt = (ls > 0.) * ls

        # --------------- central element ----------------
        # calculates mass balance for day num (Bt) using smelt, fmelt and imelt

        # effective snow melt can't exceed amount of snow
        smelt = np.minimum(snow, pot_smelt)

        # firn melt is proportional to excess snow melt
        pot_fmelt = (pot_smelt - smelt) * pddf[1] / pddf[0]

        # effective firn melt can't exceed amount of firn
        fmelt = np.minimum(firn, pot_fmelt)

        # ice melt is proportional to excess firn melt
        imelt = (pot_fmelt - fmelt) * pddf[2] / pddf[1]

        snow -= smelt
        firn -= fmelt
        ice -= imelt

        # daily mass balance of day num
        bt_arr = np.roll(bt_arr, 1, axis=0)
        bt_arr[0, :] = acc - smelt - fmelt - imelt

        # cumulative mass balance of day num
        btc = snow + firn + ice

        # ------------------- minima --------------------
        # snow is tested if minima has ocurred

        asmin = smin
        smin = (snow_1 < smin) * (snow > snow_1) * (snow_1 <= snow_2) * snow_1

        # in case that none of the above conditions is true,keep old smin
        # but if snow == 0, then smin must be 0
        smin = smin + (smin == 0.) * (snow > 0) * asmin

        # end of the first year (tstart[1])
        if day_abs / 365 == np.floor(day_abs / 365) and (year == 1): # year == 1 bcs was already increased by 1
            smin = snow  # snow minima is set to the actual snow surface

        # end of every following year
        if day_abs / 365 == np.floor(day_abs / 365) and (year > 1):
            ice = ice + firn           # firn turns into ice
            firn = smin                # last winters snow which was not melted during ablation season turns into firn
            snow = snow - smin         # snow that has now turned into firn is substracted from the snow layer
            smin = snow                # new smin is set to the actual snow surface

        # --------------------- course of the snow cover  ---------------------
        # arrays of the snow surface of the three last calculation steps are required to determine minima
        snow_2 = snow_1
        snow_1 = snow

        b_arr = np.roll(b_arr, 1, axis=0)
        b_arr[0, :, 0], b_arr[0, :, 1], b_arr[0, :, 2], b_arr[0, :, 3], \
        b_arr[0, :, 4], b_arr[0, :, 5], b_arr[0, :, 6] = snow, \
                   firn, ice, btc, Ttz, (Ttz > 0.) * Ttz, Ttz_eff

        if (day_abs / 3650000 == np.floor(day_abs / 3650000)):
            print('Calculated ' + str(int(year)) + ' years.')

    return pdds, Ba_arr, Bm_arr, param_arr, year_mb, Tm_arr


# *********************************************************************************************************************
# calculate monthly SMB
def daily_to_monthly_v2(bt_arr, year, t_arr_jul, months, mn, task):

    a = (months[1, mn])[0]
    if task == 'sum':
        monthly = np.sum(bt_arr[0:a, :], axis=0)
    if task == 'mean':
        monthly = np.mean(bt_arr[0:a, :], axis=0)

    monthly = np.insert(monthly, 0, (mn/12+year))

    return monthly


# *********************************************************************************************************************
# calculate glacier-wide SMB
def glacier_wide_smb(B, df_hypso):
    B_glw = 0
    for ind2, i2 in enumerate(B):
        B_glw += i2 * df_hypso['area'][ind2] * 10 ** 6
    B_glw = B_glw / (df_hypso['area'].sum() * 10 ** 6)
    maxB = np.amax(B, axis=0)
    minB = np.amin(B, axis=0)

    return B_glw, maxB, minB


# *********************************************************************************************************************
# calculate mass input above the ELA (assuming a static glacier!, hence only a rough indicator)
def mass_input(dem, ELA, BB):
    pos = np.where(BB > 0)
    BBacc = np.sum(BB[pos]) / 10**9  # and also convert to km3 (= GT as density is always 1 in the model)
    return BBacc


# *********************************************************************************************************************
# calculate annual mass balance parameters per DEM pixel
def smb_param_dem(t_arr_jul, b_arr, dem, df_hypso):

    # extract the last 365 days from arrays and conduct mass balance analysis
    t_arr_jul_yr = t_arr_jul[0:365]
    b_arr_yr = b_arr[0:365, :, :]

    # write to output arrays of size identical to dem
    pdd_sum = np.sum(b_arr_yr[:, :, 5], axis=0)  # sum up positive daily mean T over the last 365 days
    pdd_sum_eff = np.sum(b_arr_yr[:, :, 6], axis=0)  # sum up positive T (incl. daily cycle) over the last 365 days
    B = b_arr_yr[0, :, 3] - b_arr_yr[-1, :, 3]  # use difference between first and last day of the last 365 days
    maxT = np.amax(b_arr_yr[:, :, 4], axis=0)
    minT = np.amin(b_arr_yr[:, :, 4], axis=0)
    start_abl = np.zeros_like(dem)
    end_abl = np.zeros_like(dem)
    dur_abl = np.zeros_like(dem)
    for ind, i in enumerate(dem):
        if np.sum(b_arr_yr[:, ind, 5]) > 0:  # check that there is melt
            melt_period = np.where(b_arr_yr[:, ind, 4] > 0)
            start_abl[ind] = np.amin(t_arr_jul_yr[melt_period[0]])
            end_abl[ind] = np.amax(t_arr_jul_yr[melt_period[0]])
            dur_abl[ind] = end_abl[ind] - start_abl[ind]
            # else: not needed because 0 (default value in the three arrays) stands for not melt

    BB = B * df_hypso['area'] * 10 ** 6

    pdds = [start_abl, end_abl, dur_abl, maxT, minT, pdd_sum, pdd_sum_eff, B, BB]

    return pdds


# *********************************************************************************************************************
# use linear regression to calculate the ELA and surface mass balance gradient of the ablation area
def smb_param_glacier_wide(B, dem, df_hypso, step):

    abl_area = (np.where(B[:] < 0))[0]
    acc_area = (np.where(B[:] >= 0))[0]
    if abl_area.size > 0:  # abaltion area needs to exist
        slope, intercept, r_value, p_value, std_err = stats.linregress(dem[abl_area], B[abl_area])
        dbdz = slope * 100
    else: # if no ablation area exists
        dbdz = 0

    if abl_area.size > 0 and acc_area.size > 0:  # both ablation and accumulation area need to exists
        ELA = df_hypso['elevation'][acc_area[-1]] - step * (B[acc_area[-1]]) / (
            B[abl_area[0]] * (-1) + B[acc_area[-1]])
        AAR = (df_hypso['area'][acc_area[:-1]].sum() + (
            df_hypso['area'][acc_area[-1]] + df_hypso['area'][abl_area[0]]) *
            (B[acc_area[-1]]) / (B[abl_area[0]] * (-1) + B[acc_area[-1]])) / df_hypso['area'].sum()
    elif abl_area.size == 0:
        AAR = 1
        ELA = 0
    elif acc_area.size == 0:
        AAR = 0
        ELA = 0

    return dbdz, ELA, AAR


# ******************************************************************************************************************
# calculate refreezing
def refreezing(refreeze, refreeze_parameterization, num, T0m, Tg, dem):

    """ This module calculates the retention factor Cmax (see Reeh, 1991). Cmax is either (i) assumed 0.6 according
    to Reeh (1991) or (ii) is calculated as a function of firn temperature according to Pfeffer et al. (1991) and with
    help of a firn density fucntion by Reeh et al. (2005).
    Reeh (1991) result in a fixed retention fraction over all years and elevations.
    Pfeffer et al (1991) / Reeh et al. (2005) result in a retention fraction that varies as a function of
    annual air temperatures, thus varying over time and elevation.
    In the combined Pfeffer et al. (1991) / Reeh et al. (2005) approach, firn temperature is approximated based on
    annual mean air temperature. Firn temperature here is not the average annual firn temperature but
    represents more the spring firn temperature at relatively shallow depth. It can thus be approximated by
    annual mean air temperature.
    In case of using approach (ii), an error in the equation by Pfeffer et al. (1991) was corrected (STILL
    NEEDS TO BE IMPLEMENTED). """

    if refreeze:

        # required parameters
        retention_reeh91 = 0.6  # Fixed value
        Lm = 334000  # [J/kg] latent heat of water
        Ci = 1950  # [J K-1 kg-1] volumetric heat capacity ice at approx. -5°C
        rho_pc = 0.83  # [kg dm-3]  Porel close-off density

        # *************************************** start of model run ***************************************************
        # Here the initial retention fraction is defined. Is simply set to 0 as already
        # beginning next year a retention fraction is defined based on either a fixed value (Reeh_1991) or
        # based on air temperature (Pfeffer_1991+Reeh_2005)
        if num == 0:
            r_max = dem * 0  # Array of zeros - set retention fraction to 0

        # ***************************************** start of each year ************************************************
        # update annual mean temperature, ice temperature (cold content) and refreezing potential at the beginning
        # of the year:  originally update was placed at end of winter balance, this, however, could lead to issues
        # in the calculation of refreezing potential if the new r_max is lower and the refreezing potential is
        # currently exhausted. in this case excess refreezing will be converted to runoff.

        # calculate r_max
        if num > 0:  # do in all other cases. Note that r_max is calculated only every end of a year

            # fixed fraction of the water equivalent of the snow cover can be refrozen
            if refreeze_parameterization == 'Reeh_1991':
                r_max = retention_reeh91 + dem * 0

            # refreezing potential depends on firn T and firn rho. firn density calculation from Reeh et al. (2005)
            if refreeze_parameterization == 'Pfeffer_1991+Reeh_2005':
                # determine the relevant values
                ann_mean_T = T0m[year] - dem * Tg  # Ta, in "year", along DEM axis
                # Following Reeh et al. (2005): calculate firn density and convert from kg m-3 to kg dm-1
                rho_firn = (625 + 18.7 * ann_mean_T + 0.293 * ann_mean_T ** 2.) / 1000
                # Calculate the ratio of melt/accumulation (= Cmax) according to Pfeffer et al. (1991)
                r_max = (Ci/Lm * (-1) * ann_mean_T + ((rho_pc - rho_firn) / rho_firn)) * \
                                (1 + ((rho_pc - rho_firn) / rho_firn))**(-1.)
                r_max = np.where(r_max < 0, 0, r_max)

            print('======================================================')
            print('mean ann_mean_T:           ', np.mean(ann_mean_T))
            print('mean rho_firn:             ', np.mean(rho_firn))
            print('mean r_max:        ', np.mean(r_max))
            print('======================================================')

    else:
        r_max = dem * 0  # array of zeros

    return r_max
