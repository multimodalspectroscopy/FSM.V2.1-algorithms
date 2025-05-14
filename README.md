# FSM.V2.1-algorithms

This repository provides implementations of different physiological signal processing algorithms designed to work with multi-wavelength LED data for estimating tissue HHB, HBO2, oxCCO concentrations and oxygen saturation (StO2). The repo is structured into three branches, each containing one of the following algorithms:

- **UCLN Algorithm**
- **SRS Algorithm**
- **Dual Slope Algorithm**

---

## 📁 Branch Overview

### 🔷 `UCLN_algorithm` Branch

#### Files:
- `ucln.py` — Main UCLN algorithm implementation.
- `defaults.csv` — Extinction coefficients for various wavelengths.
- `__init__.py` — Exports the main `UCLN` function.

#### Purpose:
The UCLN algorithm processes raw LED signal data to chromophore concentrations (HbO, HHb, and oxCCO) using the Modified Beer-Lambert Law and wavelength-dependent attenuation correction.

#### How to Use:

```python
from ucln import UCLN
import pandas as pd

# Load your raw data into a DataFrame
data = pd.read_csv("your_led_data.csv")

# Run the UCLN algorithm
results = UCLN(data)

# Unpack the results
(conc_a1, conc_a2, conc_a3,
 conc_b1, conc_b2, conc_b3,
 atten_a1, atten_a2, atten_a3,
 atten_b1, atten_b2, atten_b3,
 wavelengths) = results
````

---

### 🔶 `SRS_algorithm` Branch

#### Files:

* `srs_dualSlope.py` — Contains the `SRS` function.
* `defaults.csv` — Shared extinction coefficient file.
* `__init__.py` — Exports `SRS` and `dual_slope_wavelength`.

#### Purpose:

The SRS algorithm estimates StO2 using the slope of attenuation versus distance for each wavelength.

#### How to Use:

```python
from srs_dualSlope import SRS

# Provide attenuation DataFrames for each detector group (precomputed using UCLN's `process_time_series`)
sto2_a, sto2_b, inv_ext_coeffs = SRS(atten_a1, atten_a2, atten_a3, atten_b1, atten_b2, atten_b3)

print("StO2 for Group A:", sto2_a)
print("StO2 for Group B:", sto2_b)
```

---

### 🔷 `Dual_Slope_algorithm` Branch

#### Files:

* `srs_dualSlope.py` — Includes both `SRS` and `dual_slope_wavelength`.
* `defaults.csv` — Extinction coefficients shared across algorithms.
* `__init__.py` — Exports both main functions.

#### Purpose:

The Dual Slope algorithm uses intensity readings at different detector separations to compute the absorption coefficient, which is then used to estimate StO2. It's a more robust variant compared to traditional SRS.

#### How to Use:

```python
from srs_dualSlope import dual_slope_wavelength
import pandas as pd
import numpy as np

# Load preprocessed data
data = pd.read_excel("cleaned_data.xlsx")

# Define known inputs
wavelengths = np.array([782, 801, 808, 828, 848, 887])
LED_A_det_seps = [3, 5]  # e.g., detector distances in cm
LED_B_det_seps = [5, 3]
inv_ext_coeffs = np.linalg.pinv(pd.read_csv("defaults.csv").iloc[:, 1:3].to_numpy())  # Example inverse

# Compute StO2
sto2 = dual_slope_wavelength(data, wavelengths, LED_A_det_seps, LED_B_det_seps, inv_ext_coeffs)

print("Dual Slope StO2:", sto2)
```

---

## 🧪 Data Requirements

All algorithms expect input data in a consistent format:

* Raw LED signal CSV or Excel file with columns such as:

  * `Time`, `System Time (s)`, `Sample Time (s)`
  * LED signal columns like `LED_A_808_DET1`, `LED_B_887_DET3`, etc.
  * Dark signal columns like `LED_A_DARK_DET1`, etc.
* `defaults.csv` file must be present in `src/concentrations_ucln_srs/`, containing extinction coefficients for each wavelength.

---

## 🔁 Interoperability

You can use the output from `UCLN` as input to both `SRS` and `dual_slope_wavelength'.

---

## 📄 License

MIT License
