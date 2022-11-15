"""
    *** Numerical calculation of surface mass balance profiles using a simple degree-day approach. ***

Purpose of the calculation is to investigate the influence of different annual temperature amplitudes on the
surface mass balance gradient in the abaltion area.

Eventually, the surface mass balance should be used to drive an ice flow model and thereby to study the influence
of changing surface mass balance gradients on advance and reatreat behaviour of glaciers during the last glacial period

The profiles are calculated for two different annual amplitudes of surface temperature. In addition, both simulations
can be perturbed by an offset in air temperature and precipitation (totalling four scenarios if this option is used).
"""

import numpy as np
import csv
import pandas as pd


# *********************************************************************************************************************
# write output table nr. 1 (as XLSX)
def write_output_table_xls1(outfolder, T_zpcl, df_out, pdds):
    with pd.ExcelWriter(outfolder + 'out_table_num.xlsx') as writer:
        for ind, i in enumerate(T_zpcl):
            for i2 in range(2, 13):
                df_out.iloc[:, i2] = pdds[ind, i2 - 1]
                df_out.iloc[:, i2] = pdds[ind, i2 - 1]
            df_out.to_excel(writer, sheet_name='Scenario_' + str(ind + 1))


# *********************************************************************************************************************
# write output table nr. 2 (as XLSX)
def write_output_table_xls2(outfolder, param_arr, T_raw, T_zpcl, year_start, year_end, mb_years):
    with pd.ExcelWriter(outfolder + 'out_table_annual_parameters_yrs_' +
                        str(int(year_start))+'-'+str(int(year_end))+'.xlsx') as writer:
        df_out = pd.DataFrame(columns=['year', 'B (m w.e.)', 'B_max (m w.e.)', 'B_min (m w.e.)', 'db/dz (m (100 m)-1 yr-1)',
                                       'ELA (m a.s.l.)', 'AAR (-)', 'B_acc (GT)', 'T_GISP_raw (C)', 'MAAT@400 m asl.(°C)',
                                       'Tjuly@400 m asl.(°C)', 'Tjan@400 m asl.(°C)', 'P@400 m asl.(m yr-1)',
                                       'MAAT@ELA (°C)', 'Tjuly@ELA (°C)',
                                       'Tjan@ELA (°C)', 'P@ELA (m yr-1)', 'T_amplitude (K)', 'mean_of_r_max ()'])

        years = np.arange(year_end, year_start, 1)
        years = years[::-1]  # reverse array
        # check length, usually there is one mb year less
        # than total years (bcs. of mb years corresponding to hydrol. years)
        if len(years) > mb_years:
            years = years[:mb_years]

        for ind, i in enumerate(T_zpcl):
            df_out.iloc[:, 0] = years
            for i2 in range(0, 7):
                df_out.iloc[:, i2+1] = param_arr[ind, i2, :]
            df_out.iloc[:, 8] = T_raw[:mb_years]
            for i2 in range(7, len(param_arr[ind, :, :])):
                df_out.iloc[:, i2+2] = param_arr[ind, i2, :]
            df_out.to_excel(writer, sheet_name='Scenario_'+str(ind + 1))


# *********************************************************************************************************************
# write output tables nr. 2 to len(T_zpcl) + 1 (as individual CSV files)
def write_output_tables_csv_vertical(outfolder, b_arr, ind, dem, t_resolution, content):

    # create x-axis label
    header = dem.tolist()
    header.insert(0, 'time')

    b_arr = np.round(b_arr, 5)

    outfile_fin = outfolder + t_resolution + '_' + content + '_case' + str(ind + 1) + '.txt'

    with open(outfile_fin, "w") as outcsv:
        # configure writer to write standard csv file
        writer = csv.writer(outcsv, delimiter='\t', quotechar='|',
                            quoting=csv.QUOTE_NONE, lineterminator='\n')  # ,escapechar='\\'
        writer.writerow(header)
        writer.writerows(b_arr)


# # *******************************************************************************************************************
# # write output tables nr. 2 to len(T_zpcl) + 1 (as individual CSV files)
# def write_output_tables_csv_horizontal(outfolder, bm_arr, dates, ind, content):
#
#     # create x-axis label
#     header = []
#     for num, i in enumerate(dates[0]):
#         header.append(str(i) + '/' + str(dates[1][num]))
#
#     header.insert(0, 'elevation')
#
#     outfile_fin = outfolder + 'monthly_'+content+'_case' + str(ind + 1) + '.txt'
#     bm_arr = np.transpose(bm_arr)
#
#     with open(outfile_fin, "w") as outcsv:
#         # configure writer to write standard csv file
#         writer = csv.writer(outcsv, delimiter='\t', quotechar='|', quoting=csv.QUOTE_NONE, lineterminator='\n')
#         writer.writerow(header)
#         writer.writerows(bm_arr)
#

# # *******************************************************************************************************************
# # convert table of daily b to monthly b
# # IS THIS PIECE OF CODE STILL USED? APPEARS TO BE CALLED NOWHERE
# def daily_to_monthly(bt_arr, t_years, t_start):
#     ml = 365/12  # length of a month in days
#     months_s = np.arange(0, 365, ml)
#     months = []
#     for m in months_s:
#         months.append(np.arange(np.ceil(m), m + ml - 0.001, 1, dtype=int))  # -0.001, otherwise one element too many
#
#     years = np.arange(365, t_years*365, 365, dtype=int)
#     years_s = np.insert(years, 0, t_start)  # start julday of every year
#     years_e = np.arange(365, (t_years-1)*365+0.1, 365, dtype=int)  # end julday of every year
#     years_e = np.append(years_e, (t_years-1)*365+t_start)
#     years = np.array((years_s, years_e))
#     years_c = years - t_start  # same as 'years' but shifted by t_start to be able to read directly from bt_arr
#
#     bm_arr = []
#     years_monthly = []
#     months_monthly = []
#
#     for ind1, i1 in enumerate(years[0]):
#
#         if ind1 == 0: # search starting month, only first year, otherwise always 1. month
#             for ind2, i2 in enumerate(months):
#                 c = np.where(years[0,0] >= i2[0] and years[0,0] <= i2[-1], True, False)
#                 if c:
#                     s_month = ind2
#         else:
#             s_month = 0
#
#         year = ind1 + 1
#         bt_rel = bt_arr[:, years_c[0, ind1]: years_c[1, ind1]]
#
#         for ind3, i3 in enumerate(months):
#             try:
#                 bt_rel_m = bt_rel[:, i3[:]]
#                 bm_arr.append(np.sum(bt_rel_m, axis=1))
#                 years_monthly.append(year)
#                 months_monthly.append(ind3 + 1 + s_month)
#             except:
#                 try:
#                     # some months have 30, some 31 days, hence possible that last month does not fit by one day
#                     # thus: try with one day less
#                     bt_rel_m = bt_rel[:, i3[:-1]]
#                     bm_arr.append(np.sum(bt_rel_m, axis=1))
#                     years_monthly.append(year)
#                     months_monthly.append(ind3 + 1 + s_month)
#                 except:
#                     pass
#
#     bm_arr = np.array(bm_arr)
#
#     dates = (years_monthly, months_monthly)
#
#     return bm_arr, dates
