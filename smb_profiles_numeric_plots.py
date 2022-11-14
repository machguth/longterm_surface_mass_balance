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
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import math
import numpy as np


# *********************************************************************************************************************
def def_fonts():
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
    return(font_red, font_blue, font_black)


# *********************************************************************************************************************
# Plot 1: PDDs vs. elevation
def plot_fig_1(outfolder, T_zpcl, dem, pdds):
    fonts = def_fonts()
    out_fig = outfolder + 'pdds_as_func_of_z_num.pdf'

    for ind, i in enumerate(T_zpcl):
        plt.plot(dem, pdds[ind][6][:])

    plt.xlabel('elevation (m a.s.l.)')
    plt.ylabel('positive degree days')
    plt.title("Positive degree days as function of elevation and climate")
    plt.savefig(out_fig)


# *********************************************************************************************************************
# Plot 2: PDDs vs. July temperatures
def plot_fig_2(outfolder, maxT, minT, pdds, TAa):
    fonts = def_fonts()
    out_fig = outfolder + 'pdds_as_func_of_T_num.pdf'

    maxT = math.ceil(np.amax(maxT) / 2) * 2  # get axis range
    minT = math.floor(np.amin(minT) / 2) * 2  # get axis range

    # plot only two lines because PPD per Tjuly do not depend on T_offset
    plt.plot(pdds[0][4][:], pdds[0][6][:], color='r', linewidth=1.0)
    plt.plot(pdds[1][4][:], pdds[1][6][:], color='b', linewidth=1.0)

    # place the text annotations
    plt.text(math.ceil(max(pdds[0][4][:]))-((maxT + 1) / 5), math.ceil(max(pdds[0][6][:])), r'$T_{A}$' +
             ': {:} K'.format(TAa[0][0]), fontdict=fonts[0])
    plt.text(math.ceil(max(pdds[1][4][:])), math.ceil(max(pdds[1][6][:])), r'$T_{A}$' +
             ': {:} K'.format(TAa[1][0]), fontdict=fonts[0])

    plt.axis([-1, maxT, 0, 2000])
    plt.xlabel(r'T$_{July}$ (°C)')
    plt.ylabel('positive degree days (count)')
    plt.title('Positive degree days as function of July temperature')
    plt.savefig(out_fig)


# *********************************************************************************************************************
# Plot 3: mass balance profiles
def plot_fig_3(clim_info, TAa, Tsd, T_zpcl, Tg, p_a, pddf, T_off, p_off,
               outfolder, maxB, minB, zmm, dem, B, B_glw, df_hypso,
               dbdz, ELA, AAR, z_paleoclim):

    fonts = def_fonts()

    matx, mntx = 1, 0.5  # major and minor ticks (in m w.e.) for x-axis of main (ax1) plot
    maty, mnty = 500, 100  # major and minor ticks (in m) for y-axis of all three plots

    p_a = np.array(p_a)*1000

    specifier = 'Taa{:.2f}'.format(TAa[0][0]) + ' {:.2f}'.format(TAa[1][0]) + '__T{:.2f}'.format(
                clim_info[0, 0]) + ' {:.2f}'.format(clim_info[1, 0]) + '__Tsd{:.2f}'.format(
                Tsd[0]) + ' {:.2f}'.format(Tsd[1]) \
                + '__P{:.3f}'.format(clim_info[0, 3]) + ' {:.3f}'.format(clim_info[1, 3]) + '__Tg{:.2f}'.format(
                Tg[0] * 1000) + ' {:2f}'.format(Tg[1]) + '__Pg{:}'.format(p_a[0]) \
                + ' {:}'.format(p_a[1]) + '__PDD{:.1f}'.format(pddf[0]) + ' {:.1f}'.format(
                pddf[1]) # + '__T_offset{:.1f}'.format(T_offset) + '__p_offset{:.2f}'.format(p_offset)

    out_fig = outfolder + 'smb_prof_num_'+specifier+'.pdf'

    # Create figure and axes
    fig = plt.figure(figsize=[9, 5])
    gs = GridSpec(1, 5)  # number of rows and columns
    ax1 = fig.add_subplot(gs[0, :-2])  # First row, first column
    ax2 = fig.add_subplot(gs[0, 4])  # First row, second column
    ax3 = fig.add_subplot(gs[0, 3])  # First row, second column

    # get range of X-axis
    maxB = math.ceil(np.amax(maxB) / 2) * 2  # get axis range
    minB = math.floor(np.amin(minB) / 2) * 2  # get axis range

    # # define the axes
    # if maxB - minB < 6:
    #     smb_axis = np.linspace(maxB, minB, int((maxB - minB)/0.25 + 1))
    # else:
    #     smb_axis = np.linspace(maxB, minB, int((maxB - minB) / 0.5 + 1))

    # Define major and minor ticks
    ymajor_ticks = np.arange(np.ceil(zmm[1] / maty) * maty, zmm[0] + 1, maty)
    yminor_ticks = np.arange(np.ceil(zmm[1] / mnty) * mnty, zmm[0] + 1, mnty)

    xmajor_ticks = np.arange(minB, maxB + 0.1, matx)
    xminor_ticks = np.arange(minB, maxB + 0.1, mntx)

    ax1.set_xticks(xmajor_ticks)
    ax1.set_xticks(xminor_ticks, minor=True)
    ax1.set_yticks(ymajor_ticks)
    ax1.set_yticks(yminor_ticks, minor=True)

    ax1.grid(which='minor', color='k', linewidth=0.2, linestyle='dotted', zorder=0)
    ax1.grid(which='major', color='k', linewidth=0.2, linestyle='dotted', zorder=0)

    # make vertical zero line
    ax1.plot([0, 0], [zmm[1], zmm[0]], color='k', linewidth=0.8)

    # *** First plot: plot the dTAa
    colours = ['r', 'b', 'r', 'b']
    for ind, i in enumerate(T_zpcl):
        if ind > 1:
            ax1.plot(B[ind], dem, color=colours[ind], linewidth=1.0, linestyle='dashed')
        else:
            ax1.plot(B[ind], dem, color=colours[ind], linewidth=1.6)

    # *** third plot: plot the hypsometry
    ax2.barh(dem, df_hypso['area']*0.001, height=50, align='edge')  # align='center',
    ax2.set_ylim([zmm[1], zmm[0]])
    ax2.set_yticks(ymajor_ticks)
    ax2.set_yticks(yminor_ticks, minor=True)
    minorLocator2 = AutoMinorLocator()
    ax2.xaxis.set_minor_locator(minorLocator2)

    # *** Secondd plot: plot the delta b
    ax3.plot(B[1] - B[0], dem, color='black', linewidth=1.0, linestyle='dashed')
    # if T_offset != 0 or p_offset != 0:
    #     ax3.plot(B[2] - B[0], dem, color='r', linewidth=1.0, linestyle='dotted')
    #     ax3.plot(B[3] - B[1], dem, color='b', linewidth=1.0, linestyle='dotted')
    # else:
    #     ax3.plot(B[1] - B[0], dem, color='black', linewidth=1.0, linestyle='dashed')

    ax3.set_ylim([zmm[1], zmm[0]])
    ax3.set_yticks(ymajor_ticks)
    ax3.set_yticks(yminor_ticks, minor=True)
    minorLocator3 = AutoMinorLocator()
    ax3.xaxis.set_minor_locator(minorLocator3)
    ax3.xaxis.grid(which='minor', color='k', linewidth=0.2, linestyle='dotted', zorder=0)
    ax3.xaxis.grid(which='major', color='k', linewidth=0.2, linestyle='dotted', zorder=0)

    # *** First plot: create and place legend
    yTop = zmm[0] - (zmm[0]-zmm[1])/20
    xLeft = minB + (maxB - minB)/28
    xShift = (maxB - minB)/24
    yShift = (zmm[0]-zmm[1])/23

    for ind, i in enumerate(T_zpcl):
        yi = yTop - ind*yShift/1.3
        if ind > 1:
            ax1.plot([xLeft, xLeft + xShift], [yi, yi], color=colours[ind], linewidth=1.0, linestyle='dashed')
            ax1.text(xLeft + xShift + xShift / 3, yi - yShift / 4,
                     r'B' + str(ind + 1)+': $db/dz$: {:.2f} m w.e./100m'.format(dbdz[ind]) +
                     r' | ELA: {:} m a.s.l.'.format(int(ELA[ind])),
                     fontdict=fonts[2], size=6)  # (100 m)$^{-1}$
        else:
            ax1.plot([xLeft, xLeft + xShift], [yi, yi], color=colours[ind], linewidth=1.6)
            ax1.text(xLeft + xShift + xShift / 3, yi - yShift / 4,
                     r'B'+str(ind+1)+': $db/dz$: {:.2f} m w.e./100m'.format(dbdz[ind]) +
                     r' | ELA: {:} m a.s.l.'.format(int(ELA[ind])),
                     fontdict=fonts[2], size=6)  # (100 m)$^{-1}$

    # if T_offset != 0 or p_offset != 0:
    #     TOffStr = str(T_offset)
    #     pOffStr = str(p_offset)
    #     yi = yTop - (ind + 1) * yShift / 1.3
    #     ax1.plot([xLeft, xLeft + xShift], [yi, yi], color='r', linewidth=1.0, linestyle='dotted')
    #     yi = yTop - (ind+2)*yShift/1.3
    #     ax1.plot([xLeft, xLeft + xShift], [yi, yi], color='b', linewidth=1.0, linestyle='dotted')
    #     ax1.text(xLeft + xShift + xShift/5, yi + yShift / 5,
    #              r'}', fontdict=fonts[2], size=12)
    #     ax1.text(xLeft + 2 * xShift, yi + yShift / 4,
    #              r'$\Delta b$ (m w.e.) for $T_{offset}$=' + TOffStr + r' K, $p_{offset}$=' + pOffStr + r' m yr$^{-1}$',
    #              fontdict=fonts[2], size=6)

    # now plot all the information on the applied variable combinations
    # for better readability first only the variable names and units
    ax1.text(xLeft, yTop - (ind + 2.0) * yShift, r'$T_A$ (K): ',
             fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 2.7) * yShift, r'$\sigma_{daily}$ (K): ',
             fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 3.4) * yShift, r'$MAAT$ (K) @{:}'.format(z_paleoclim) +
             ' m a.s.l.:', fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 4.1) * yShift, r'$P$ (m yr$^{-1}$) ' + r'@{:}'.format(z_paleoclim) +
             ' m a.s.l.:', fontdict=fonts[2], size=6)
    # ax1.text(xLeft, yTop - (ind+4.1)*yShift, r'$P$ (m yr$^{-1}$): ', fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 4.8) * yShift, r'$\gamma_{T}$ (K km$^{-1}$): ',
             fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 5.5) * yShift, r'$\gamma_{P}$ ((1000 m)$^{-1}$ yr$^{-1}$): ',
             fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 6.2) * yShift, r'PDDF$_{snow}$ (mm K$^{-1}$ d$^{-1}$): ',
             fontdict=fonts[2], size=6)
    ax1.text(xLeft, yTop - (ind + 6.9) * yShift, r'PDDF$_{ice}$ (mm K$^{-1}$ d$^{-1}$): ',
             fontdict=fonts[2], size=6)

    dv = 8.3
    fs = 5

    # then plot the values, thereby plotting a value only once if identical for both cases
    if TAa[0][0] == TAa[1][0]:
        ax1.text(xLeft + dv * xShift, yTop - (ind + 2.0) * yShift, '{:.2f}'.format(TAa[0][0]),
                 fontdict=fonts[2], size=fs)
    else:
        ax1.text(xLeft+dv*xShift, yTop - (ind+2.0)*yShift, '{:.2f}'.format(TAa[0][0]) + '; ' +
                 '{:.2f}'.format(TAa[1][0]), fontdict=fonts[2], size=fs)
    if Tsd[0] == Tsd[1]:
        ax1.text(xLeft + dv * xShift, yTop - (ind + 2.7) * yShift, '{:}'.format(Tsd[0]),
                 fontdict=fonts[2], size=fs)
    else:
        ax1.text(xLeft+dv*xShift, yTop - (ind+2.7)*yShift, '{:}'.format(Tsd[0])+'; '+'{:}'.format(Tsd[1]),
                 fontdict=fonts[2], size=fs)
    if clim_info[0, 0] == clim_info[1, 0]:
        ax1.text(xLeft + dv * xShift, yTop - (ind + 3.4) * yShift, '{:.2f}'.format(clim_info[0, 0]),
                 fontdict=fonts[2], size=fs)
    else:
        ax1.text(xLeft + dv * xShift, yTop - (ind+3.4) * yShift, '{:.2f}'.format(clim_info[0, 0]) +
                 '; ' + '{:.2f}'.format(clim_info[1, 0]), fontdict=fonts[2], size=fs)
    if clim_info[0, 3] == clim_info[1, 3]:
        ax1.text(xLeft + dv * xShift, yTop - (ind + 4.1) * yShift, '{:.3f}'.format(clim_info[0, 3]),
                 fontdict=fonts[2], size=fs)
    else:
        ax1.text(xLeft+dv*xShift, yTop - (ind+4.1)*yShift, '{:.3f}'.format(clim_info[0, 3]) +
                 '; ' + '{:.3f}'.format(clim_info[1, 3]), fontdict=fonts[2], size=fs)
    if Tg[0] == Tg[1]:
        ax1.text(xLeft+dv*xShift, yTop - (ind+4.8)*yShift, '{:.2f}'.format(Tg[0]*1000),
                 fontdict=fonts[2], size=fs)
    else:
        ax1.text(xLeft + dv * xShift, yTop - (ind + 4.8) * yShift, '{:.2f}'.format(Tg[0] * 1000) +
                 '; ' + '{:.2f}'.format(Tg[1]), fontdict=fonts[2], size=fs)
    if p_a[0] == p_a[1]:
        ax1.text(xLeft + dv * xShift, yTop - (ind + 5.5) * yShift, '{:}'.format(p_a[0]),
                 fontdict=fonts[2], size=fs)
    else:
        ax1.text(xLeft+dv*xShift, yTop - (ind+5.5)*yShift, '{:}'.format(p_a[0]) +
                 '; ' + '{:}'.format(p_a[1]), fontdict=fonts[2], size=fs)
    ax1.text(xLeft+dv*xShift, yTop - (ind+6.2)*yShift, '{:.1f}'.format(pddf[0]),
             fontdict=fonts[2], size=fs)
    ax1.text(xLeft+dv*xShift, yTop - (ind+6.9)*yShift, '{:.1f}'.format(pddf[2]),
             fontdict=fonts[2], size=fs)

    for ind2, i2 in enumerate(T_zpcl):  # here add information on specific net balance and AAR
        if ((ind2 + 1) % 2) == 0:
            colo = 'b'
        else:
            colo = 'r'
        ax1.text(xLeft, yTop - (ind + 8 + 0.7 * ind2) * yShift, r'B'+str(ind2 + 1) +
                 ' (m w.e.): ', fontdict=fonts[2], color=colo, size=fs)
        ax1.text(xLeft+dv*xShift, yTop - (ind + 8 + 0.7 * ind2) * yShift, '{:.3f}'.format(B_glw[ind2]),
                 fontdict=fonts[2], color=colo, size=fs)
        ax1.text(xLeft, yTop - (ind + 8 + 0.7 * (ind2+len(T_zpcl))) * yShift, r'AAR'+str(ind2+1)+' (-): ',
                 fontdict=fonts[2], color=colo, size=fs)
        ax1.text(xLeft+dv*xShift, yTop - (ind + 8+0.7*(ind2+len(T_zpcl))) * yShift, '{:.3f}'.format(AAR[ind2]),
                 fontdict=fonts[2], color=colo, size=fs)

    # rect = patches.Rectangle((xLeft, yTop - (ind+2.5)*yShift), (maxB - minB)*2.8/8, -4.9*yShift,
    #                      linewidth=0, edgecolor='none', facecolor='white',
    #                      alpha=0.8, zorder = 2)

    # Add the patch to the Axes
    # ax1.add_patch(rect)

    # All three plots: Labelling
    ax1.axis([minB, maxB, zmm[1], zmm[0]])
    ax1.set_ylabel('elevation $z$ (m a.s.l.)')
    ax1.set_xlabel('$b$ (m w.e.)')
    ax2.set_xlabel('A (in 1000 km$^2$)')
    ax1.set_title('surface mass balance', fontdict=fonts[2])
    ax2.set_title('hypsometry', fontdict=fonts[2])
    # if T_offset != 0 or p_offset != 0:
    #     ax3.set_title('$\Delta b$ (B4-B2; B3-B1)', fontdict=fonts[2])
    # else:
    #     ax3.set_title('$\Delta b$ (B2-B1)', fontdict=fonts[2])
    ax3.set_title('$\Delta b$ (B2-B1)', fontdict=fonts[2])
    ax3.set_xlabel('$\Delta b$ (m w.e.)')
    plt.tight_layout()
    plt.savefig(out_fig)
