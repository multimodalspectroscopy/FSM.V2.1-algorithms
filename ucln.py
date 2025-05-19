import numpy as np
import pandas as pd

def UCLN(data: pd.DataFrame):
    # Define wavelengths used in FSM-V2.1
    wavelengths = [782, 801, 808, 828, 848, 887]

    # Load extinction coefficients
    extinction_path = "src/concentrations_ucln_srs/defaults.csv"
    epsilon = pd.read_csv(extinction_path)

    ext_coeffs = []
    ext_coeffs_cor = []

    for wl in wavelengths:
        row = epsilon[epsilon['wavelength'] == wl]
        if not row.empty:
            ext = row[['HbO2', 'HHb', 'CCO']].values.flatten()
            corr = row['wl_dep'].values[0]
            ext_coeffs.append(ext)
            ext_coeffs_cor.append(ext * corr)
        else:
            raise ValueError(f"Wavelength {wl} not found in extinction coefficient data.")

    ext_coeffs_cor = np.array(ext_coeffs_cor)
    ext_coeffs_inv = np.linalg.pinv(ext_coeffs_cor)

    # Optode distances in mm
    optode_distance_A = np.array([30, 40, 50])
    optode_distance_B = np.array([50, 40, 30])
    DPF = 4.99

    # Parse intensity values from data
    def get_group_data(prefix, det):
        return np.stack([
            pd.to_numeric(data[f'{prefix}_{wl}_DET{det}'], errors='coerce').fillna(1e-10).values
            for wl in wavelengths
        ], axis=0)  # shape: (6, n_samples)

    group_a_1 = get_group_data('LED_A', 1)
    group_a_2 = get_group_data('LED_A', 2)
    group_a_3 = get_group_data('LED_A', 3)
    group_b_1 = get_group_data('LED_B', 1)
    group_b_2 = get_group_data('LED_B', 2)
    group_b_3 = get_group_data('LED_B', 3)

    # Compute attenuation
    def calc_attenuation(group):
        return np.log10(100 / np.clip(group, 1e-10, None))  # shape: (6, n_samples)

    att_a_1 = calc_attenuation(group_a_1)
    att_a_2 = calc_attenuation(group_a_2)
    att_a_3 = calc_attenuation(group_a_3)
    att_b_1 = calc_attenuation(group_b_1)
    att_b_2 = calc_attenuation(group_b_2)
    att_b_3 = calc_attenuation(group_b_3)

    # Delta attenuation
    def delta_atten(att):
        return att[:, 0:1] - att  # broadcast over all samples

    att_a_1 = delta_atten(att_a_1)
    att_a_2 = delta_atten(att_a_2)
    att_a_3 = delta_atten(att_a_3)
    att_b_1 = delta_atten(att_b_1)
    att_b_2 = delta_atten(att_b_2)
    att_b_3 = delta_atten(att_b_3)

    # Concentration calculation
    def get_conc(delta_att, distance):
        conc = 10000 * ext_coeffs_inv @ (delta_att / (distance * DPF))
        return pd.DataFrame(conc.T, columns=['HHb', 'HbO', 'oxCCO'])

    conc_a_1_df = get_conc(att_a_1, optode_distance_A[0])
    conc_a_2_df = get_conc(att_a_2, optode_distance_A[1])
    conc_a_3_df = get_conc(att_a_3, optode_distance_A[2])
    conc_b_1_df = get_conc(att_b_1, optode_distance_B[0])
    conc_b_2_df = get_conc(att_b_2, optode_distance_B[1])
    conc_b_3_df = get_conc(att_b_3, optode_distance_B[2])

    # Convert attenuation arrays to DataFrames
    def att_to_df(att):
        return pd.DataFrame(att.T, columns=[f'{wl}' for wl in wavelengths])

    atten_a_1 = att_to_df(att_a_1)
    atten_a_2 = att_to_df(att_a_2)
    atten_a_3 = att_to_df(att_a_3)
    atten_b_1 = att_to_df(att_b_1)
    atten_b_2 = att_to_df(att_b_2)
    atten_b_3 = att_to_df(att_b_3)

    return (
        conc_a_1_df, conc_a_2_df, conc_a_3_df,
        conc_b_1_df, conc_b_2_df, conc_b_3_df,
        atten_a_1, atten_a_2, atten_a_3,
        atten_b_1, atten_b_2, atten_b_3, wavelengths
    )

