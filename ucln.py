import numpy as np
import pandas as pd

def process_time_series(data):
    data_df = pd.DataFrame(data)
    print('data_df head 5 rows', data_df.head(5))
    data = data_df
    print('data head 5 rows', data.head(5))
    data.columns = data.columns = ['Time', 'System Time (s)', 'Sample Time (s)',
                    'LED_A_782_DET1', 'LED_A_782_DET2', 'LED_A_782_DET3', 
                    'LED_A_801_DET1', 'LED_A_801_DET2', 'LED_A_801_DET3', 
                    'LED_A_808_DET1', 'LED_A_808_DET2', 'LED_A_808_DET3', 
                    'LED_A_828_DET1', 'LED_A_828_DET2', 'LED_A_828_DET3',
                    'LED_A_848_DET1', 'LED_A_848_DET2', 'LED_A_848_DET3', 
                    'LED_A_887_DET1', 'LED_A_887_DET2', 'LED_A_887_DET3', 
                    'LED_A_DARK_DET1', 'LED_A_DARK_DET2', 'LED_A_DARK_DET3',
                    'LED_B_782_DET1', 'LED_B_782_DET2', 'LED_B_782_DET3', 
                    'LED_B_801_DET1', 'LED_B_801_DET2', 'LED_B_801_DET3', 
                    'LED_B_808_DET1', 'LED_B_808_DET2', 'LED_B_808_DET3', 
                    'LED_B_828_DET1', 'LED_B_828_DET2', 'LED_B_828_DET3',
                    'LED_B_848_DET1', 'LED_B_848_DET2', 'LED_B_848_DET3', 
                    'LED_B_887_DET1', 'LED_B_887_DET2', 'LED_B_887_DET3', 
                    'LED_B_DARK_DET1', 'LED_B_DARK_DET2', 'LED_B_DARK_DET3']

    Time = data['Time']
    wavelengths = np.array([782, 801, 808, 828, 848, 887 ])
    dpf = 4.99
    extinction_coefficients_path = "src/concentrations_ucln_srs/defaults.csv"
    extinction_coefficients = pd.read_csv(extinction_coefficients_path)
    
    ex_co = []
    for wavelength in wavelengths:
        for i, j in enumerate(extinction_coefficients.iloc[:, 0]):
            if wavelength == j:
                ex_co.append(extinction_coefficients.iloc[i, :])
    
    ex_co = np.array(ex_co)
    ext_coeffs = ex_co[:, 1:4]  
    ext_coeffs_inv = np.linalg.pinv(ext_coeffs)
    
    wavelength_dependency = ex_co[:, 4]
    
    GroupA_Detector1 = [f'LED_A_{w}_DET1' for w in wavelengths]
    GroupA_Detector2 = [f'LED_A_{w}_DET2' for w in wavelengths]
    GroupA_Detector3 = [f'LED_A_{w}_DET3' for w in wavelengths]
    GroupB_Detector1 = [f'LED_B_{w}_DET1' for w in wavelengths]
    GroupB_Detector2 = [f'LED_B_{w}_DET2' for w in wavelengths]
    GroupB_Detector3 = [f'LED_B_{w}_DET3' for w in wavelengths]
    
    # Extracting the data for each group
    group_a_1 = data[GroupA_Detector1].apply(pd.to_numeric, errors='coerce').fillna(0).values
    group_a_2 = data[GroupA_Detector2].apply(pd.to_numeric, errors='coerce').fillna(0).values
    group_a_3 = data[GroupA_Detector3].apply(pd.to_numeric, errors='coerce').fillna(0).values
    group_b_1 = data[GroupB_Detector1].apply(pd.to_numeric, errors='coerce').fillna(0).values
    group_b_2 = data[GroupB_Detector2].apply(pd.to_numeric, errors='coerce').fillna(0).values
    group_b_3 = data[GroupB_Detector3].apply(pd.to_numeric, errors='coerce').fillna(0).values
    
    x, y = group_a_1.shape
    atten_a_1 = np.log10(group_a_1[0, :] / group_a_1)
    atten_a_2 = np.log10(group_a_2[0, :] / group_a_2)
    atten_a_3 = np.log10(group_a_3[0, :] / group_a_3)
    atten_b_1 = np.log10(group_b_1[0, :] / group_b_1)
    atten_b_2 = np.log10(group_b_2[0, :] / group_b_2)
    atten_b_3 = np.log10(group_b_3[0, :] / group_b_3)
    
    atten_a_1 = pd.DataFrame(atten_a_1, columns=GroupA_Detector1)
    atten_a_2 = pd.DataFrame(atten_a_2, columns=GroupA_Detector2)
    atten_a_3 = pd.DataFrame(atten_a_3, columns=GroupA_Detector3)
    atten_b_1 = pd.DataFrame(atten_b_1, columns=GroupB_Detector1)
    atten_b_2 = pd.DataFrame(atten_b_2, columns=GroupB_Detector2)
    atten_b_3 = pd.DataFrame(atten_b_3, columns=GroupB_Detector3)
    
    groups = {
        'Group_A_782': [atten_a_1[GroupA_Detector1[0]], atten_a_2[GroupA_Detector2[0]], atten_a_3[GroupA_Detector3[0]]],
        'Group_A_801': [atten_a_1[GroupA_Detector1[1]], atten_a_2[GroupA_Detector2[1]], atten_a_3[GroupA_Detector3[1]]],
        'Group_A_808': [atten_a_1[GroupA_Detector1[2]], atten_a_2[GroupA_Detector2[2]], atten_a_3[GroupA_Detector3[2]]],
        'Group_A_828': [atten_a_1[GroupA_Detector1[3]], atten_a_2[GroupA_Detector2[3]], atten_a_3[GroupA_Detector3[3]]],
        'Group_A_848': [atten_a_1[GroupA_Detector1[4]], atten_a_2[GroupA_Detector2[4]], atten_a_3[GroupA_Detector3[4]]],
        'Group_A_887': [atten_a_1[GroupA_Detector1[5]], atten_a_2[GroupA_Detector2[5]], atten_a_3[GroupA_Detector3[5]]],
    }

    print('group_a_1', group_a_1)
    print('group_a_2', group_a_2)
    print('group_a_3', group_a_3)
    print('group_b_1', group_b_1)
    print('group_b_2', group_b_2)
    print('group_b_3', group_b_3)
    
    return atten_a_1, atten_a_2, atten_a_3, atten_b_1, atten_b_2, atten_b_3, wavelengths, dpf, wavelength_dependency, ext_coeffs_inv


def UCLN(data):
        # Using the formula function to extract the necessary variables
        atten_a_1, atten_a_2, atten_a_3, atten_b_1, atten_b_2, atten_b_3, wavelengths, dpf, wavelength_dependency, ext_coeffs_inv = process_time_series(data)
    
        # Apply wavelength dependency
        atten_a_1_wldep = atten_a_1 / wavelength_dependency
        atten_a_2_wldep = atten_a_2 / wavelength_dependency
        atten_a_3_wldep = atten_a_3 / wavelength_dependency
        atten_b_1_wldep = atten_b_1 / wavelength_dependency
        atten_b_2_wldep = atten_b_2 / wavelength_dependency
        atten_b_3_wldep = atten_b_3 / wavelength_dependency

        print('atten_a_1_wldep', atten_a_1_wldep)
        print('atten_a_2_wldep', atten_a_2_wldep)
        print('atten_a_3_wldep', atten_a_3_wldep)
        print('atten_b_1_wldep', atten_b_1_wldep)
        print('atten_b_2_wldep', atten_b_2_wldep)
        print('atten_b_3_wldep', atten_b_3_wldep)
        
        # Calculate concentrations using the inverse of extinction coefficients
        conc_a_1 = np.transpose(np.matmul(ext_coeffs_inv, atten_a_1_wldep.T) * (1 / (3 * dpf)))
        conc_a_2 = np.transpose(np.matmul(ext_coeffs_inv, atten_a_2_wldep.T) * (1 / (4 * dpf)))
        conc_a_3 = np.transpose(np.matmul(ext_coeffs_inv, atten_a_3_wldep.T) * (1 / (5 * dpf)))
        conc_b_1 = np.transpose(np.matmul(ext_coeffs_inv, atten_b_1_wldep.T) * (1 / (3 * dpf)))
        conc_b_2 = np.transpose(np.matmul(ext_coeffs_inv, atten_b_2_wldep.T) * (1 / (4 * dpf)))
        conc_b_3 = np.transpose(np.matmul(ext_coeffs_inv, atten_b_3_wldep.T) * (1 / (5 * dpf)))

        print('conc_a_1', conc_a_1)
        print('conc_a_2', conc_a_2)
        print('conc_a_3', conc_a_3)
        print('conc_b_1', conc_b_1)
        print('conc_b_2', conc_b_2)
        print('conc_b_3', conc_b_3)

        conc_a_1_df = pd.DataFrame(conc_a_1)
        conc_a_2_df = pd.DataFrame(conc_a_2)
        conc_a_3_df = pd.DataFrame(conc_a_3)
        conc_b_1_df = pd.DataFrame(conc_b_1)
        conc_b_2_df = pd.DataFrame(conc_b_2)
        conc_b_3_df = pd.DataFrame(conc_b_3)


        conc_a_1_df.columns = ['HbO', 'HHb', 'oxCCO']
        conc_a_2_df.columns = ['HbO', 'HHb', 'oxCCO']
        conc_a_3_df.columns = ['HbO', 'HHb', 'oxCCO']
        conc_b_1_df.columns = ['HbO', 'HHb', 'oxCCO']
        conc_b_2_df.columns = ['HbO', 'HHb', 'oxCCO']
        conc_b_3_df.columns = ['HbO', 'HHb', 'oxCCO']

        # Print the first 5 rows
        print('conc_a_1_df:\n', conc_a_1_df.head())
        print('conc_a_2_df:\n', conc_a_2_df.head())
        print('conc_a_3_df:\n', conc_a_3_df.head())
        print('conc_b_1_df:\n', conc_b_1_df.head())
        print('conc_b_2_df:\n', conc_b_2_df.head())
        print('conc_b_3_df:\n', conc_b_3_df.head())

        return conc_a_1_df, conc_a_2_df, conc_a_3_df, conc_b_1_df, conc_b_2_df, conc_b_3_df, atten_a_1, atten_a_2, atten_a_3, atten_b_1, atten_b_2, atten_b_3, wavelengths

