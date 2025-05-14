import numpy as np
import pandas as pd
#==================================SRS======================================887=============================================


def SRS(atten_a_1, atten_a_2, atten_a_3, atten_b_1, atten_b_2, atten_b_3):
    
    print('atten_a_1', atten_a_1)
    print('atten_a_2', atten_a_2)
    print('atten_a_3', atten_a_3)

    print('atten_b_1', atten_b_1)
    print('atten_b_2', atten_b_2)
    print('atten_b_3', atten_b_3)
    wavelengths = np.array([782, 801, 808, 828, 848, 887 ])

    srs_extinction_coefficients = pd.read_csv(
    "src/concentrations_ucln_srs/defaults.csv")

    srs_ex_co = []
    for wavelength in wavelengths:
        for i, j in enumerate(srs_extinction_coefficients.iloc[:, 0]):
            if wavelength == j:
                srs_ex_co.append(srs_extinction_coefficients.iloc[i, :])

    srs_ex_co = np.array(srs_ex_co)
    srs_ext_coeffs = srs_ex_co[:, 1:3]  # just take oxy and deoxy for now
    print('ext_coeffs', srs_ext_coeffs)
    #ext_coeffs_t = ext_coeffs.T
    srs_ext_coeffs_inv = np.linalg.pinv(srs_ext_coeffs)
    print('ext_coeffs', srs_ext_coeffs_inv)
    atten_a_1.columns = ['LED_A_782_DET1', 'LED_A_801_DET1', 'LED_A_808_DET1', 'LED_A_828_DET1', 'LED_A_848_DET1',
                     'LED_A_887_DET1']
    atten_a_2.columns = ['LED_A_782_DET2', 'LED_A_801_DET2', 'LED_A_808_DET2', 'LED_A_828_DET2', 'LED_A_848_DET2',
                     'LED_A_887_DET2']
    atten_a_3.columns = ['LED_A_782_DET3', 'LED_A_801_DET3', 'LED_A_808_DET3', 'LED_A_828_DET3', 'LED_A_848_DET3',
                     'LED_A_887_DET3']

    atten_b_1.columns = ['LED_B_782_DET1', 'LED_B_801_DET1', 'LED_B_808_DET1', 'LED_B_828_DET1', 'LED_B_848_DET1',
                     'LED_B_887_DET1']
    atten_b_2.columns = ['LED_B_782_DET2', 'LED_B_801_DET2', 'LED_B_808_DET2', 'LED_B_828_DET2', 'LED_B_848_DET2',
                     'LED_B_887_DET2']
    atten_b_3.columns = ['LED_B_782_DET3', 'LED_B_801_DET3', 'LED_B_808_DET3', 'LED_B_828_DET3', 'LED_B_848_DET3',
                     'LED_B_887_DET3']

    Group_A_782 = [atten_a_1['LED_A_782_DET1'], atten_a_2['LED_A_782_DET2'], atten_a_3['LED_A_782_DET3']]
    Group_A_801 = [atten_a_1['LED_A_801_DET1'], atten_a_2['LED_A_801_DET2'], atten_a_3['LED_A_801_DET3']]
    Group_A_808 = [atten_a_1['LED_A_808_DET1'], atten_a_2['LED_A_808_DET2'], atten_a_3['LED_A_808_DET3']]
    Group_A_828 = [atten_a_1['LED_A_828_DET1'], atten_a_2['LED_A_828_DET2'], atten_a_3['LED_A_828_DET3']]
    Group_A_848 = [atten_a_1['LED_A_848_DET1'], atten_a_2['LED_A_848_DET2'], atten_a_3['LED_A_848_DET3']]
    Group_A_887 = [atten_a_1['LED_A_887_DET1'], atten_a_2['LED_A_887_DET2'], atten_a_3['LED_A_887_DET3']]

    Group_B_782 = [atten_b_1['LED_B_782_DET1'], atten_b_2['LED_B_782_DET2'], atten_b_3['LED_B_782_DET3']]
    Group_B_801 = [atten_b_1['LED_B_801_DET1'], atten_b_2['LED_B_801_DET2'], atten_b_3['LED_B_801_DET3']]
    Group_B_808 = [atten_b_1['LED_B_808_DET1'], atten_b_2['LED_B_808_DET2'], atten_b_3['LED_B_808_DET3']]
    Group_B_828 = [atten_b_1['LED_B_828_DET1'], atten_b_2['LED_B_828_DET2'], atten_b_3['LED_B_828_DET3']]
    Group_B_848 = [atten_b_1['LED_B_848_DET1'], atten_b_2['LED_B_848_DET2'], atten_b_3['LED_B_848_DET3']]
    Group_B_887 = [atten_b_1['LED_B_887_DET1'], atten_b_2['LED_B_887_DET2'], atten_b_3['LED_B_887_DET3']]

    detector_A_distance = np.array([[3], [4], [5]])
    detector_B_distance = np.array([[5], [4], [3]])

    ##Do group A first##
    A_groupA = np.hstack((detector_A_distance,np.ones_like(detector_A_distance) ))

    # 782nm
    atten_wldep_A_782 = np.array(Group_A_782, dtype=float)  # change atten_wldp with group from above
    m_782_A = np.zeros((1, len(atten_wldep_A_782[0])))
    c_782_A = np.zeros((1, len(atten_wldep_A_782[0])))
    for x in range(0, len(atten_wldep_A_782[0])):
        m_782_A[0, x], c_782_A[0, x] = np.linalg.lstsq(A_groupA, atten_wldep_A_782[:, x], rcond=None)[0]

    # 801nm
    atten_wldep_A_801 = np.array(Group_A_801, dtype=float)  # change atten_wldp with group from above
    m_801_A = np.zeros((1, len(atten_wldep_A_801[0])))
    c_801_A = np.zeros((1, len(atten_wldep_A_801[0])))
    for x in range(0, len(atten_wldep_A_801[0])):
        m_801_A[0, x], c_801_A[0, x] = np.linalg.lstsq(A_groupA, atten_wldep_A_801[:, x], rcond=None)[0]

    # 808nm
    atten_wldep_A_808 = np.array(Group_A_808, dtype=float)  # change atten_wldp with group from above
    m_808_A = np.zeros((1, len(atten_wldep_A_808[0])))
    c_808_A = np.zeros((1, len(atten_wldep_A_808[0])))
    for x in range(0, len(atten_wldep_A_808[0])):
        m_808_A[0, x], c_808_A[0, x] = np.linalg.lstsq(A_groupA, atten_wldep_A_808[:, x], rcond=None)[0]

    # 828nm
    atten_wldep_A_828 = np.array(Group_A_828, dtype=float)  # change atten_wldp with group from above
    m_828_A = np.zeros((1, len(atten_wldep_A_828[0])))
    c_828_A = np.zeros((1, len(atten_wldep_A_828[0])))
    for x in range(0, len(atten_wldep_A_808[0])):
        m_828_A[0, x], c_828_A[0, x] = np.linalg.lstsq(A_groupA, atten_wldep_A_828[:, x], rcond=None)[0]

    # 848nm
    atten_wldep_A_848 = np.array(Group_A_848, dtype=float)  # change atten_wldp with group from above
    m_848_A = np.zeros((1, len(atten_wldep_A_848[0])))
    c_848_A = np.zeros((1, len(atten_wldep_A_848[0])))
    for x in range(0, len(atten_wldep_A_848[0])):
        m_848_A[0, x], c_848_A[0, x] = np.linalg.lstsq(A_groupA, atten_wldep_A_848[:, x], rcond=None)[0]


    # 887nm
    atten_wldep_A_887 = np.array(Group_A_887, dtype=float)  # change atten_wldp with group from above
    m_887_A = np.zeros((1, len(atten_wldep_A_887[0])))
    c_887_A = np.zeros((1, len(atten_wldep_A_887[0])))
    for x in range(0, len(atten_wldep_A_887[0])):
        m_887_A[0, x], c_887_A[0, x] = np.linalg.lstsq(A_groupA, atten_wldep_A_887[:, x], rcond=None)[0]

    ##Do group B##
    A_groupB = np.hstack((detector_B_distance,np.ones_like(detector_B_distance) ))

    # 782nm
    atten_wldep_B_782 = np.array(Group_B_782, dtype=float)  # change atten_wldp with group from above
    m_782_B = np.zeros((1, len(atten_wldep_B_782[0])))
    c_782_B = np.zeros((1, len(atten_wldep_B_782[0])))
    for x in range(0, len(atten_wldep_B_782[0])):
        m_782_B[0, x], c_782_B[0, x] = np.linalg.lstsq(A_groupB, atten_wldep_B_782[:, x], rcond=None)[0]

    # 801nm
    atten_wldep_B_801 = np.array(Group_B_801, dtype=float)  # change atten_wldp with group from above
    m_801_B = np.zeros((1, len(atten_wldep_B_801[0])))
    c_801_B = np.zeros((1, len(atten_wldep_B_801[0])))
    for x in range(0, len(atten_wldep_B_801[0])):
        m_801_B[0, x], c_801_B[0, x] = np.linalg.lstsq(A_groupB, atten_wldep_B_801[:, x], rcond=None)[0]

    # 808nm
    atten_wldep_B_808 = np.array(Group_B_808, dtype=float)  # change atten_wldp with group from above
    m_808_B = np.zeros((1, len(atten_wldep_B_808[0])))
    c_808_B = np.zeros((1, len(atten_wldep_B_808[0])))
    for x in range(0, len(atten_wldep_B_808[0])):
        m_808_B[0, x], c_808_B[0, x] = np.linalg.lstsq(A_groupB, atten_wldep_B_808[:, x], rcond=None)[0]

    # 828nm
    atten_wldep_B_828 = np.array(Group_B_828, dtype=float)  # change atten_wldp with group from above
    m_828_B = np.zeros((1, len(atten_wldep_B_828[0])))
    c_828_B = np.zeros((1, len(atten_wldep_B_828[0])))
    for x in range(0, len(atten_wldep_B_828[0])):
        m_828_B[0, x], c_828_B[0, x] = np.linalg.lstsq(A_groupB, atten_wldep_B_828[:, x], rcond=None)[0]


    # 848nm
    atten_wldep_B_848 = np.array(Group_B_848, dtype=float)  # change atten_wldp with group from above
    m_848_B = np.zeros((1, len(atten_wldep_B_848[0])))
    c_848_B = np.zeros((1, len(atten_wldep_B_848[0])))
    for x in range(0, len(atten_wldep_B_848[0])):
        m_848_B[0, x], c_848_B[0, x] = np.linalg.lstsq(A_groupB, atten_wldep_B_848[:, x], rcond=None)[0]

    # 887nm
    atten_wldep_B_887 = np.array(Group_B_887, dtype=float)  # change atten_wldp with group from above
    m_887_B = np.zeros((1, len(atten_wldep_B_887[0])))
    c_887_B = np.zeros((1, len(atten_wldep_B_887[0])))
    for x in range(0, len(atten_wldep_B_887[0])):
        m_887_B[0, x], c_887_B[0, x] = np.linalg.lstsq(A_groupB, atten_wldep_B_887[:, x], rcond=None)[0]

    h = 6.3e-4
    # group A
    k_mua_782_A = (1 / (3 * (1 - (h * wavelengths[0])))) * (np.log(10) * m_782_A - (2 / np.mean(detector_A_distance)))**2
    k_mua_801_A = (1 / (3 * (1 - (h * wavelengths[1])))) * (np.log(10) * m_801_A - (2 / np.mean(detector_A_distance)))**2
    k_mua_808_A = (1 / (3 * (1 - (h * wavelengths[2])))) * (np.log(10) * m_808_A - (2 / np.mean(detector_A_distance)))**2
    k_mua_828_A = (1 / (3 * (1 - (h * wavelengths[3])))) * (np.log(10) * m_828_A - (2 / np.mean(detector_A_distance)))**2
    k_mua_848_A = (1 / (3 * (1 - (h * wavelengths[4])))) * (np.log(10) * m_848_A - (2 / np.mean(detector_A_distance)))**2
    k_mua_887_A = (1 / (3 * (1 - (h * wavelengths[5])))) * (np.log(10) * m_887_A - (2 / np.mean(detector_A_distance)))**2

    #groupB
    k_mua_782_B = (1 / (3 * (1 - (h * wavelengths[0])))) * (np.log(10) * m_782_B - (2 / np.mean(detector_B_distance)))**2
    k_mua_801_B = (1 / (3 * (1 - (h * wavelengths[1])))) * (np.log(10) * m_801_B - (2 / np.mean(detector_B_distance)))**2
    k_mua_808_B = (1 / (3 * (1 - (h * wavelengths[2])))) * (np.log(10) * m_808_B - (2 / np.mean(detector_B_distance)))**2
    k_mua_828_B = (1 / (3 * (1 - (h * wavelengths[3])))) * (np.log(10) * m_828_B - (2 / np.mean(detector_B_distance)))**2
    k_mua_848_B = (1 / (3 * (1 - (h * wavelengths[4])))) * (np.log(10) * m_848_B - (2 / np.mean(detector_B_distance)))**2
    k_mua_887_B = (1 / (3 * (1 - (h * wavelengths[5])))) * (np.log(10) * m_887_B - (2 / np.mean(detector_B_distance)))**2


    k_mua_A = np.vstack([k_mua_782_A, k_mua_801_A,k_mua_808_A,k_mua_828_A,k_mua_848_A, k_mua_887_A])
    k_mua_B = np.vstack([k_mua_782_B, k_mua_801_B,k_mua_808_B,k_mua_828_B,k_mua_848_B, k_mua_887_B])

    C_A = np.matmul(srs_ext_coeffs_inv, k_mua_A)
    C_B = np.matmul(srs_ext_coeffs_inv, k_mua_B)

    oxy_a = C_A[0]
    deoxy_a = C_A[1]

    oxy_b = C_B[0]
    deoxy_b = C_B[1]

    sto2_a = (oxy_a / (oxy_a + deoxy_a)) * 100
    sto2_b = (oxy_b / (oxy_b + deoxy_b)) * 100

    print( "k_mua_A", k_mua_A)
    print( "k_mua_B", k_mua_B)
    print( "C_A", C_A)
    print( "C_B", C_B)
    print( "oxy_a", oxy_a)
    print( "deoxy_a", deoxy_a)
    print( "oxy_a", oxy_b)
    print( "deoxy_a", deoxy_b)

    return sto2_a, sto2_b, srs_ext_coeffs_inv


#==================================Dual Slope===================================================================================


def dual_slope_wavelength(data, wavelengths, LED_A_det_seps, LED_B_det_seps, srs_ext_coeffs_inv,  dsf=8):
    ds_mua_results = {}

    
    for wavelength in wavelengths:
        col_A_det1 = f'LED_A_{wavelength}_DET1'
        col_A_det3 = f'LED_A_{wavelength}_DET3'
        col_B_det1 = f'LED_B_{wavelength}_DET1'
        col_B_det3 = f'LED_B_{wavelength}_DET3'
        
        LED_A_data = data[[col_A_det1, col_A_det3]].values
        LED_B_data = data[[col_B_det3, col_B_det1]].values
        
        ss_LED_A = []
        ss_LED_B = []
        
        for row in LED_A_data:
            num_detectors = len(LED_A_det_seps)
            avg_det = np.mean(LED_A_det_seps)
            var_det = np.var(LED_A_det_seps, ddof=1)
            slope_sum = 0
            for i in [num_detectors // 2]:
                r_N = LED_A_det_seps[num_detectors - (i + 1)]
                r_1 = LED_A_det_seps[i]
                I_N = row[num_detectors - (i + 1)]
                I_1 = row[i]  # Added missing I_1
                log_ratio_r = 2 * np.log(r_N / r_1)
                log_ratio_I = np.log(I_N / I_1)
                slope_sum += (r_N - avg_det) * (log_ratio_r + log_ratio_I)
            ss_LED_A.append(slope_sum / (num_detectors * var_det))

        for row in LED_B_data:
            num_detectors = len(LED_B_det_seps)
            avg_det = np.mean(LED_B_det_seps)
            var_det = np.var(LED_B_det_seps, ddof=1)
            slope_sum = 0
            for i in [num_detectors // 2]:
                r_N = LED_B_det_seps[num_detectors - (i + 1)]
                r_1 = LED_B_det_seps[i]
                I_N = row[num_detectors - (i + 1)]
                I_1 = row[i]  # Added missing I_1
                log_ratio_r = 2 * np.log(r_N / r_1)
                log_ratio_I = np.log(I_N / I_1)
                slope_sum += (r_N - avg_det) * (log_ratio_r + log_ratio_I)
            ss_LED_B.append(slope_sum / (num_detectors * var_det))
        
        ss_LED_A = np.array(ss_LED_A) / dsf
        ss_LED_B = np.array(ss_LED_B) / dsf


        # Ensure the results are stored properly for each wavelength
        ds_mua_results[wavelength] = -(ss_LED_A + ss_LED_B) / (2 * dsf)

    # Now this line will work correctly since ds_mua_results contains all wavelengths
    mua_dual_slope = np.vstack([ds_mua_results[wl] for wl in [782, 801, 808, 828, 848, 887]])


    conc_dual_slope = np.matmul(srs_ext_coeffs_inv, mua_dual_slope)

    # Separate the oxygenated and deoxygenated concentrations.
    oxy_dual_slope = conc_dual_slope[0, :]
    deoxy_dual_slope = conc_dual_slope[1, :]

    # Compute StO2 for each time point.
    sto2_dual_slope = (oxy_dual_slope / (oxy_dual_slope + deoxy_dual_slope)) * 100

    return sto2_dual_slope
