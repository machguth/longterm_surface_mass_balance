"""calculating surface mass balance profiles using a most simple degree-day approach
for GG.0444 "Alpine Cryosphere and Geomorphology" The profiles are calculated for current
climate and for +1°C climate and plotted together in one plot"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import math
import numpy as np
import csv
import pandas as pd
from scipy import stats

#################################################################################################
#                           input, output, parameter definitions
#################################################################################################

# --------------------------
T_ampl = [8.6, 15]
#T_ampl = [7, 9]    # (K) amplitude is half the annual temperature range
zmm = [3600, 300]  # (m a.s.l.) max and min elevation z for calculation
step = 100          # (m) elevation step for calculation of B

Tg = [0.006, 0.006]

#T_zstat = [265.77, 259.37] # (K)
T_zstat = [266.77, 260.37] # (K)
z_stat = 1000 # (m a.s.l.) elevation for which climate input is specified

z_paleoclimate = 400 # (m a.s.l.) elevation for which quantitative paleoclimate data exist

pddf = [3.297, 8.791] # degree day factors (mm K-1 d-1), first for snow, second for ice

#p = [0.167, 0.126] # precipitation at z_stat (m yr-1)
#p = [0.126, 0.126] # precipitation at z_stat (m yr-1)
p = [0.23, 0.23] # precipitation at z_stat (m yr-1)
pg = np.array([12, 12]) # precipitation gradient (percent per 100m-1 yr-1)

T_offset = 0 # (K) offset for simulating impact of climate change
p_offset = 0.0 # (m yr-1) offset for simulating impact of climate change

hypsometry = r'C:/Horst/modeling/modelanalysis/dbdz/s1_hypsometry.xlsx'
#hypsometry = r'C:/Horst/modeling/modelanalysis/dbdz/s2_hypsometry.xlsx'

# --------------------------
outfolder = r'C:/Horst/modeling/modelanalysis/dbdz/'
# --------------------------

# font definitions for comments in output figures
font_red = {'family': 'sans-serif',
        'color':  'r',
        'weight': 'normal',
        'size': 10,
        }
font_blue = {'family': 'sans-serif',
        'color':  'b',
        'weight': 'normal',
        'size': 10,
        }
font_black = {'family': 'sans-serif',
        'color':  'k',
        'weight': 'normal',
        'size': 10,
        }

#################################################################################################
#                                      preparations
#################################################################################################

specifier = 'Ta{:}'.format(T_ampl[0]) +' {:}'.format(T_ampl[1]) +'__T{:.2f}'.format(T_zstat[0]) + ' {:.2f}'.format(T_zstat[1]) \
            + '__P{:.3f}'.format(p[0]) +' {:.3f}'.format(p[1]) +'__Tg{:.2f}'.format(Tg[0]*1000) + ' {:2f}'.format(Tg[1]) +'__Pg{:}'.format(pg[0]) \
            + ' {:}'.format(pg[1]) +'__PDD{:.1f}'.format(pddf[0]) +' {:.1f}'.format(pddf[1])+'__T_offset{:.1f}'.format(T_offset)

# read the hypsometry
df_hypso = pd.read_excel(hypsometry)#, index_col = 'elevation')
df_hypso['area'] = np.flip(df_hypso['area'].values) # flip the order of the entries to agree with order of dem and pdds arrays
df_hypso['elevation'] = np.flip(df_hypso['elevation'].values) # flip the order of the entries to agree with order of dem and pdds arrays
# --------------------------
dem = np.linspace(zmm[0], zmm[1], (zmm[0] - zmm[1])/step+1)

# specify the time array, usually 365 elements for each day of the year
# this aray is only used to ANALYTICALLY calculate the annual sum of PDDs, no numerical calculation is performed
t_array = np.linspace(1, 365, 365)

pg = pg/10000  # convert precip gradient to scaling factor per m instead of % per 100m

# if climate change simulations are done, simply extend some arrays to include four instead of two scenarios
# the additional scenarios are identical to the first two, except for the imposed change in climate
if T_offset != 0 or p_offset != 0:
    T_zstat = [T_zstat[0], T_zstat[1], T_zstat[0] + T_offset, T_zstat[1] + T_offset]
    p = [p[0], p[1], p[0]+p_offset, p[1]+p_offset]
    pg = [pg[0], pg[1], pg[0], pg[1]]
    Tg = [Tg[0], Tg[1], Tg[0], Tg[1]]
    T_ampl = [T_ampl[0], T_ampl[1], T_ampl[0], T_ampl[1]]

tTza = np.zeros((len(T_zstat), 365))
maxT = []
minT = []
maxB = []
minB = []
# establish an array of temperature at sea level for each of the scenarios
T0 = [0] * len(T_zstat)
for ind, i in enumerate(T_zstat):
    T0[ind] = i - 273.15 + z_stat * Tg[ind]

# output array, contains: elevation, start_abl, end_abl, dur_abl, maxt, mint, pdd, B (m w.e.), B (m3 w.e.)
pdds = np.zeros((len(T_zstat), 9, len(dem)))

df_out = pd.DataFrame(columns=['elev', 'area', 's_abl', 'e_abl', 'maxt', 'mint', 'pdd', 'b (m.w.e)', 'B (m3 w.e.)'])

df_out['elev'] = dem
df_out['area'] = df_hypso['area']



#################################################################################################
#                                        computations
#################################################################################################

# --------------------------------------------------------
# some basic informations on the chosen climate
for ind, i in enumerate(T_zstat):
    MAAT_paleoclimate = T_zstat[ind] + (z_stat - z_paleoclimate) * Tg[ind]
    Ts_paleoclimate = MAAT_paleoclimate + T_ampl[ind]
    Tw_paleoclimate = MAAT_paleoclimate - T_ampl[ind]
    p_paleoclimate = (p[ind] * 1000) * (1 - (z_stat - z_paleoclimate) * pg[ind])
    scenario = ind + 1
    print('Scenario {:}'.format(scenario) +' at {:}'.format(int(z_paleoclimate)) + ' m a.s.l.')
    print('MAAT: {:.1f}'.format(MAAT_paleoclimate - 273.15) + ' °C')
    print('Tmax: {:.1f}'.format(Ts_paleoclimate - 273.15) + ' °C')
    print('Tmin: {:.1f}'.format(Tw_paleoclimate - 273.15) + ' °C')
    print('P: {:.2f}'.format(p_paleoclimate/1000) + ' m/yr')
    print('')
print('--------------------------')

# --------------------------------------------------------
# calculate PDDs for all elevation classes
for ind0, i0 in enumerate(dem):

    for ind1, i1 in enumerate(T_zstat):

        zeros_t = []  # little array to store the daynumbers of transition from neg to pos T

        # orig. form.: tTza[ind] = T0[ind] - dem * Tg + T_ampl[ind] * ((-1)*cos((2. * math.pi * (x - 30.)) / 365.))
        # note that T0 is the temperature at sea level
        k1 = T0[ind1] - i0 * Tg[ind1]
        k2 = T_ampl[ind1]
        k3 = 2. * math.pi
        k4 = k3 * 30. / 365.
        k5 = k3 / 365
        tTza[ind1] = k1 - k2 * np.cos(k5 * t_array - k4)

        # find zeros (numerically, temp. solution)
        for ind2, i2 in enumerate(tTza[ind1][:-1]):
            if tTza[ind1][ind2] < 0 and tTza[ind1][ind2 + 1] > 0:
                zeros_t = np.append(zeros_t, ind2 + 1)
            if tTza[ind1][ind2] > 0 and tTza[ind1][ind2 + 1] < 0:
                zeros_t = np.append(zeros_t, ind2 + 1)

        # integrate to calculate the total number of positive degree days
        if len(zeros_t) == 2:
            # integral k1 - k2 * cos(k5 * x - k4) dx = k1*x + k2 * sin(k4 - k5 * x)/k5
            pdd = k1 * zeros_t[1] + k2 * np.sin(k4 - k5 * zeros_t[1]) / k5 - (
            k1 * zeros_t[0] + k2 * np.sin(k4 - k5 * zeros_t[0]) / k5)
        else:
            pdd = 0
            zeros_t = [0, 0]

        pdds[ind1, 0:7, ind0] = [i0, zeros_t[0], zeros_t[1], zeros_t[1] - zeros_t[0], max(tTza[ind1]), min(tTza[ind1]),
                               pdd]

        maxT = np.append(maxT, max(tTza[ind1])) # July maximum temperatures for each elevation band and scenario are: pdds[ind1, 4, :]
        minT = np.append(minT, min(tTza[ind1]))

# --------------------------------------------------------
# calculate PDDs at z_stat
pdds_zstat = np.zeros((len(T_zstat), 7))

for ind1, i1 in enumerate(T_zstat):

    zeros_t = []  # little array to store the daynumbers of transition from neg to pos T

    # orig. form.: tTza[ind] = T0[ind] - dem * Tg + T_ampl[ind] * ((-1)*cos((2. * math.pi * (x - 30.)) / 365.))
    # note that T0 is the temperature at sea level
    k1 = T0[ind1] - z_stat * Tg[ind1]
    k2 = T_ampl[ind1]
    k3 = 2. * math.pi
    k4 = k3 * 30. / 365.
    k5 = k3 / 365
    tTza[ind1] = k1 - k2 * np.cos(k5 * t_array - k4)

    # find zeros (numerically, temp. solution)
    for ind2, i2 in enumerate(tTza[ind1][:-1]):
        if tTza[ind1][ind2] < 0 and tTza[ind1][ind2 + 1] > 0:
            zeros_t = np.append(zeros_t, ind2 + 1)
        if tTza[ind1][ind2] > 0 and tTza[ind1][ind2 + 1] < 0:
            zeros_t = np.append(zeros_t, ind2 + 1)

    # integrate to calculate the total number of positive degree days
    if len(zeros_t) == 2:
        # integral k1 - k2 * cos(k5 * x - k4) dx = k1*x + k2 * sin(k4 - k5 * x)/k5
        pdd = k1 * zeros_t[1] + k2 * np.sin(k4 - k5 * zeros_t[1]) / k5 - (
                k1 * zeros_t[0] + k2 * np.sin(k4 - k5 * zeros_t[0]) / k5)
    else:
        pdd = 0
        zeros_t = [0, 0]

    pdds_zstat[ind1, :] = [z_stat, zeros_t[0], zeros_t[1], zeros_t[1] - zeros_t[0], max(tTza[ind1]), min(tTza[ind1]),
                       pdd]

print('--------------------------')
for ind1, i1 in enumerate(T_zstat):
    swe = pdds_zstat[ind1, 6] * pddf[0]
    scenario = ind1 + 1
    print('Scenario %s' %scenario)
    print('PDDs at elevation z_stat: {:.3f}'.format(pdds_zstat[ind1, 6]))
    print('Snow w.e. needed to place ELA at z_stat: {:.3f}'.format(swe))
print('--------------------------')

# ----------------------------------------------------------
# calculate SMB

# output array, contains: smb for both test cases
B = np.zeros((len(T_zstat), len(dem)))
for ind, i in enumerate(T_zstat):
    for ind1, i1 in enumerate(pdds[ind][6]):
        if (pdds[ind][6][ind1] * pddf[0]) < ((p[ind] * 1000) * (1 + (dem[ind1] - z_stat) * pg[ind])) :
            snowmelt = pdds[ind][6][ind1] * pddf[0]
            ice_pdds = 0
            icemelt = 0
        else:
            snowmelt = (p[ind] * 1000) * (1 + (dem[ind1] - z_stat) * pg[ind])
            ice_pdds = pdds[ind][6][ind1] - snowmelt / pddf[0]
            icemelt = ice_pdds * pddf[1]

        B[ind][ind1] = ((p[ind] * 1000) * (1 + (dem[ind1] - z_stat) * pg[ind]) - snowmelt - icemelt) / 1000
        pdds[ind][7][ind1] = B[ind][ind1]
        pdds[ind][8][ind1] = B[ind][ind1] * df_hypso['area'][ind1] * 10**6

    maxB = np.append(maxB, max(B[ind]))
    minB = np.append(minB, min(B[ind]))


# ----------------------------------------------------------
# calculate glacier-wide SMB
B_glw = np.zeros(len(T_zstat))
for ind, i in enumerate(T_zstat):
    for ind2, i2 in enumerate(B[ind]):
        B_glw[ind] += i2 * df_hypso['area'][ind2] * 10**6
    B_glw[ind] = B_glw[ind] / (df_hypso['area'].sum() * 10**6)

for i in B_glw:
 print('B glacier-wide: %s' %i)
print('---')
print('difference in B glacier-wide (of first two scenarios): %s' %(B_glw[0] - B_glw[1]))

# ----------------------------------------------------------
# use linear regression to calculate the ELA and surface mass balance gradients (only ablation area)

dbdz = np.zeros(len(T_zstat))
ELA = np.zeros(len(T_zstat))
AAR = np.zeros(len(T_zstat))

for ind, i in enumerate(T_zstat):
    abl_area = (np.where(B[ind][:] < 0))[0]
    acc_area = (np.where(B[ind][:] >= 0))[0]
    slope, intercept, r_value, p_value, std_err = stats.linregress(dem[abl_area], B[ind][abl_area])
    dbdz[ind] = slope * 100

    ELA[ind] = df_hypso['elevation'][acc_area[-1]] - step * (B[ind][acc_area[-1]])/(B[ind][abl_area[0]]*(-1) + B[ind][acc_area[-1]])
    AAR[ind] = (df_hypso['area'][acc_area[:-1]].sum() + (df_hypso['area'][acc_area[-1]] + df_hypso['area'][abl_area[0]])* \
               (B[ind][acc_area[-1]])/(B[ind][abl_area[0]]*(-1) + B[ind][acc_area[-1]]))/ df_hypso['area'].sum()
    print(AAR[ind])
    #ELA[ind] = df_hypso['elevation'][abl_area[0]]
    #ELA[ind] = intercept/slope * (-1)

#################################################################################################
#                                        plot figures
#################################################################################################

# -----------------------------------------------------------------------------------------------
# Plot 1: PDDs vs. elevation
out_fig = outfolder + 'pdds_as_func_of_z.pdf'

for ind, i in enumerate(T_zstat):
    plt.plot(dem, pdds[ind][6][:])

plt.xlabel('elevation (m a.s.l.)')
plt.ylabel('positive degree days')
plt.title("Positive degree days as function of elevation and climate")
plt.savefig(out_fig)

# -----------------------------------------------------------------------------------------------
# Plot 2: PDDs vs. July temperatures
out_fig = outfolder + 'pdds_as_func_of_T.pdf'

maxT = math.ceil(max(maxT)/2)*2 # get axis range
minT = math.floor(min(minT)/2)*2 # get axis range

# plot only two lines because PPD per Tjuly do not depend on T_offset
plt.plot(pdds[0][4][:], pdds[0][6][:], color='r', linewidth=1.0)
plt.plot(pdds[1][4][:], pdds[1][6][:], color='b', linewidth=1.0)

# place the text annotations
plt.text(math.ceil(max(pdds[0][4][:]))-((maxT+1)/5), math.ceil(max(pdds[0][6][:])), r'$T_{A}$'+ ': {:} K'.format(T_ampl[0]), fontdict=font_red)
plt.text(math.ceil(max(pdds[1][4][:])), math.ceil(max(pdds[1][6][:])), r'$T_{A}$'+ ': {:} K'.format(T_ampl[1]), fontdict=font_blue)

plt.axis([-1, maxT, 0, 2000])
plt.xlabel(r'T$_{July}$ (°C)')
plt.ylabel('positive degree days (count)')
plt.title('Positive degree days as function of July temperature')
# plt.show()
plt.savefig(out_fig)

# -----------------------------------------------------------------------------------------------
# Plot 3: mass balance profiles
out_fig = outfolder + 'smb_prof_'+specifier+'.pdf'

# Create figure and axes
fig = plt.figure()

gs=GridSpec(1,4) # 2 rows, 3 columns

ax1=fig.add_subplot(gs[0,:-1]) # First row, first column
ax2=fig.add_subplot(gs[0,3]) # First row, second column


#fig,ax = plt.subplots(1,2, sharey=True, sharex=False)

# get range of X-axis
maxB = math.ceil(max(maxB)/2)*2 # get axis range
minB = math.floor(min(minB)/2)*2 # get axis range

# define the axes
if maxB - minB < 6:
    smb_axis = np.linspace(maxB, minB, (maxB - minB)/0.25 + 1)
else:
    smb_axis = np.linspace(maxB, minB, (maxB - minB) / 0.5 + 1)
#print('smb axis %s' % smb_axis)

# # draw a custom grid
# if zmm[0] - zmm[1] < 1200:
#     dem_axis = np.linspace(zmm[0], zmm[1], (zmm[0] - zmm[1])/step + 1)
# else:
#     dem_axis = np.linspace(zmm[0], zmm[1], (zmm[0] - zmm[1]) / (2*step) + 1)

for ind, i in enumerate(dem):
    ax1.plot([minB, maxB], [i, i], color='k', linewidth=0.2, linestyle='dotted', zorder = 0)
for ind, i in enumerate(smb_axis):
    ax1.plot([i, i], [zmm[1], zmm[0]], color='k', linewidth=0.2, linestyle='dotted', zorder = 0)


# make vertical zero line
ax1.plot([0, 0], [zmm[1], zmm[0]], color='k', linewidth=0.8)

# plot the data
colours = ['r', 'b', 'r', 'b']
for ind, i in enumerate(T_zstat):
    if ind > 1:
        ax1.plot(B[ind], dem, color=colours[ind], linewidth=1.0, linestyle='dashed')
    else:
        ax1.plot(B[ind], dem, color=colours[ind], linewidth=1.6)

if T_offset != 0 or p_offset != 0:
    deltaB = np.zeros(len(dem)) + maxB
    ax1.plot((deltaB + B[2] - B[0]), dem, color='r', linewidth=1.0, linestyle='dotted')
    ax1.plot((deltaB + B[3] - B[1]), dem, color='b', linewidth=1.0, linestyle='dotted')

ax2.barh(dem, df_hypso['area']*0.001, height=50, align='edge')#align='center',
ax2.set_ylim([zmm[1],zmm[0]])
#ax2.set_xlim([0,2])

# place legend
yTop = zmm[0] - (zmm[0]-zmm[1])/20
xLeft = minB + (maxB - minB)/28
xShift = (maxB - minB)/24
yShift = (zmm[0]-zmm[1])/20

for ind, i in enumerate(T_zstat):
    yi = yTop - ind*yShift
    if ind > 1:
        ax1.plot([xLeft, xLeft+xShift], [yi,yi], color=colours[ind], linewidth=1.0, linestyle='dashed')
        ax1.text(xLeft+xShift+xShift/3, yi-yShift/4,
                 r'$db/dz$: {:.2f} m w.e./100m'.format(dbdz[ind])+r' | ELA: {:} m a.s.l.'.format(int(ELA[ind])),
                 fontdict=font_black, size=6) #  (100 m)$^{-1}$
    else:
        ax1.plot([xLeft, xLeft+xShift], [yi,yi], color=colours[ind], linewidth=1.6)
        ax1.text(xLeft+xShift+xShift/3, yi-yShift/4,
                 r'$db/dz$: {:.2f} m w.e./100m'.format(dbdz[ind])+r' | ELA: {:} m a.s.l.'.format(int(ELA[ind])),
                 fontdict=font_black, size=6) #  (100 m)$^{-1}$

if T_offset != 0 or p_offset != 0:
    TOffStr = str(T_offset)
    pOffStr = str(p_offset)
    yi = yTop - (ind+1)*yShift
    ax1.plot([xLeft, xLeft+xShift], [yi,yi], color='r', linewidth=1.0, linestyle='dotted')
    yi = yTop - (ind+2)*yShift
    ax1.plot([xLeft, xLeft+xShift], [yi,yi], color='b', linewidth=1.0, linestyle='dotted')
    ax1.text(xLeft + xShift + xShift/5, yi + yShift/5,
             r'}', fontdict=font_black, size=16)
    ax1.text(xLeft + 2*xShift, yi + yShift/4,
             r'$\Delta b$ (m w.e.) for $T_{offset}$='+TOffStr+r' K, $p_{offset}$='+pOffStr+r' m yr$^{-1}$',
             fontdict=font_black, size=6)


# now plot all the information on the applied variable combinations
# for better readability first only the variable names and units
ax1.text(xLeft, yTop - (ind+3)*yShift, r'$T_A$ (K): ', fontdict=font_black, size=6)
ax1.text(xLeft, yTop - (ind+3.7) * yShift, r'$MAAT$ (K) @{:}'.format(z_stat) + ' m a.s.l.:', fontdict=font_black, size=6)
ax1.text(xLeft, yTop - (ind+4.4)*yShift, r'$P$ (m yr$^{-1}$): ', fontdict=font_black, size=6)
ax1.text(xLeft, yTop - (ind+5.1)*yShift, r'$\gamma_{T}$ (K km$^{-1}$): ', fontdict=font_black, size=6)
ax1.text(xLeft, yTop - (ind+5.8)*yShift, r'$\gamma_{P}$ (% (100 m)$^{-1}$ yr$^{-1}$): ', fontdict=font_black, size=6)
ax1.text(xLeft, yTop - (ind+6.5)*yShift, r'PDDF$_{snow}$ (mm K$^{-1}$ d$^{-1}$): ', fontdict=font_black, size=6)
ax1.text(xLeft, yTop - (ind+7.2)*yShift, r'PDDF$_{ice}$ (mm K$^{-1}$ d$^{-1}$): ', fontdict=font_black, size=6)

# ax1.text(xLeft, yTop - (ind+9.2)*yShift, r'B$_1$ (m w.e.): ', fontdict=font_black, color='r', size=6)
# ax1.text(xLeft, yTop - (ind+9.9)*yShift, r'B$_2$ (m w.e.): ', fontdict=font_black, color='b', size=6)

dv = 8.3
fs = 5

# then plot the values, thereby plotting a value only once if identical for both cases
if T_ampl[0] == T_ampl[1]:
    ax1.text(xLeft + dv * xShift, yTop - (ind + 3) * yShift, '{:}'.format(T_ampl[0]), fontdict=font_black, size=fs)
else:
    ax1.text(xLeft+dv*xShift, yTop - (ind+3)*yShift, '{:}'.format(T_ampl[0])+'; '+'{:}'.format(T_ampl[1]), fontdict=font_black, size=fs)
if T_zstat[0] == T_zstat[1]:
    ax1.text(xLeft + dv * xShift, yTop - (ind + 3.7) * yShift, '{:.2f}'.format(T_zstat[0]), fontdict=font_black, size=fs)
else:
    ax1.text(xLeft + dv * xShift, yTop - (ind+3.7) * yShift, '{:.2f}'.format(T_zstat[0]) + '; ' + '{:.2f}'.format(T_zstat[1]), fontdict=font_black, size=fs)
if p[0] == p[1]:
    ax1.text(xLeft + dv * xShift, yTop - (ind + 4.4) * yShift, '{:.3f}'.format(p[0]), fontdict=font_black, size=fs)
else:
    ax1.text(xLeft+dv*xShift, yTop - (ind+4.4)*yShift, '{:.3f}'.format(p[0])+'; '+'{:.3f}'.format(p[1]), fontdict=font_black, size=fs)
if Tg[0] == Tg[1]:
    ax1.text(xLeft+dv*xShift, yTop - (ind+5.1)*yShift, '{:.2f}'.format(Tg[0]*1000), fontdict=font_black, size=fs)
else:
    ax1.text(xLeft + dv * xShift, yTop - (ind + 5.1) * yShift, '{:.2f}'.format(Tg[0] * 1000)+'; '+'{:.2f}'.format(Tg[1]), fontdict=font_black, size=fs)
if pg[0] == pg[1]:
    ax1.text(xLeft + dv * xShift, yTop - (ind + 5.8) * yShift, '{:}'.format(pg[0]), fontdict=font_black, size=fs)
else:
    ax1.text(xLeft+dv*xShift, yTop - (ind+5.8)*yShift, '{:}'.format(pg[0])+'; '+'{:}'.format(pg[1]), fontdict=font_black, size=fs)
ax1.text(xLeft+dv*xShift, yTop - (ind+6.5)*yShift, '{:.1f}'.format(pddf[0]), fontdict=font_black, size=fs)
ax1.text(xLeft+dv*xShift, yTop - (ind+7.2)*yShift, '{:.1f}'.format(pddf[1]), fontdict=font_black, size=fs)

for ind2, i2 in enumerate(T_zstat): # here add information on specific net balance and AAR
    if ((ind2 + 1) % 2) == 0:
        colo = 'b'
    else:
        colo = 'r'
    ax1.text(xLeft, yTop - (ind + 8+0.7*ind2) * yShift, r'B'+str(ind2+1)+' (m w.e.): ', fontdict=font_black, color=colo, size=fs)
    ax1.text(xLeft+dv*xShift, yTop - (ind + 8+0.7*ind2)*yShift, '{:.3f}'.format(B_glw[ind2]), fontdict=font_black, color=colo, size=fs)
    ax1.text(xLeft, yTop - (ind + 8+0.7*(ind2+len(T_zstat))) * yShift, r'AAR'+str(ind2+1)+' (-): ', fontdict=font_black, color=colo, size=fs)
    ax1.text(xLeft+dv*xShift, yTop - (ind + 8+0.7*(ind2+len(T_zstat)))*yShift, '{:.3f}'.format(AAR[ind2]), fontdict=font_black, color=colo, size=fs)

rect = patches.Rectangle((xLeft, yTop - (ind+2.5)*yShift), (maxB - minB)*2.8/8, -4.9*yShift,
                         linewidth=0, edgecolor='none', facecolor='white',
                         alpha=0.8, zorder = 2)

# Add the patch to the Axes
ax1.add_patch(rect)


ax1.axis([minB, maxB, zmm[1], zmm[0]])
ax1.set_ylabel('elevation $z$ (m a.s.l.)')
ax1.set_xlabel('surface mass balance (m w.e.)')
ax2.set_xlabel('A (in 1000 km$^2$)')
ax1.set_title('surface mass balance', fontdict=font_black)
ax2.set_title('hypsometry', fontdict=font_black)
plt.tight_layout()
plt.savefig(out_fig)


#################################################################################################
#                                    write output table
#################################################################################################


# --------------------------
writer = pd.ExcelWriter(outfolder + 'out_table.xlsx')
for ind, i in enumerate(T_zstat):
    for i2 in range(2, 9):
       df_out.iloc[:,i2] = pdds[ind,i2]
    df_out.to_excel(writer,'Scenario_'+str(ind + 1))
writer.save()