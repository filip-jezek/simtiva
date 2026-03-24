# Eleveld 2018 PK/PD Dataset Dictionary

This dictionary documents the structural headers natively mapped across the raw Propofol non-linear mixed-effects (NONMEM) modeling datasets. 

## Dataset 1: Pharmacokinetic (PK) Model Correlates
*(Source: [supplementary_digital_content_1.txt](file:///c:/home/git/simtiva/data/supplementary_digital_content_1.txt) mapped against [supplementary_digital_content_2.txt](file:///c:/home/git/simtiva/data/supplementary_digital_content_2.txt))*

### Identifiers & Events
* **`ID`**: Unique Subject/Patient explicitly mapped across the entire database (1,033 unique individuals).
* **`DID`**: Likely a secondary individual identification tag or underlying sub-dataset tracking ID. 
* **`TIME`**: Time of the documented event (dosing or blood sampling).
* **`EVID`**: Event ID standard to NONMEM. A value of `0` denotes an **observation** (blood sample taking), and `1` denotes a **dosing event** (bolus or infusion administration).
* **`STDY`**: Study Identifier flag (`1` through `32`). Indicates the original source trial the patient data was gathered from.
* **`GRP`**: A grouping factor (`1`, `2`, `3`), typically separating sub-cohorts like neonates, children, and adults.

### Dose & Pharmacokinetic Measurements
* **`AMT`**: Amount of propofol administered during a dosing event (`EVID = 1`).
* **`RATE`**: Infusion rate for continuous administration (`EVID = 1`). A rate of `0` typically denotes an instantaneous bolus.
* **`DV`**: Dependent Variable (the actual measured Propofol concentration during sampling). Negative values denote measurements Below the Limit of Quantification (BLQ) or LLOQ flags natively bounded by logarithmic offsets.
* **`DVTY`**: Dependent Variable Type separator (`1` vs `2`). Distinguishes types of assay measurements.
* **`P1V2`**: Sample constitution: Plasma (`1`) vs Whole Blood (`2`) assays. 

### Core Demographic Covariates
* **`AGE`**: Chronological age in fractional years, spanning from 0.0027 (~1 day old) to 88 years old.
* **`WGT`**: Total Body Weight in kilograms (0.68 kg to 160 kg).
* **`HGT`**: Height in centimeters (33 cm to 200 cm).
* **`BMI`**: Body Mass Index derived from WGT and HGT.
* **`PMA`**: Post-Menstrual Age in fraction of years. Multiplying `PMA * 52.` transforms it cleanly into weeks, which the structural model uses to compute the sigmoidal clearance maturation curve (`PMW=PMA*52.`).

### Physiologically Mapped Flags
* **`M1F2` (Gender/Sex)**: `1` = Male, `2` = Female.
  * *NONMEM Proof*: Mapped algebraically: `FFMM*(2-M1F2) + FFMF*(M1F2-1)`. Setting female (`2`) zeros out male Fat-Free Mass (`FFMM*0`) and enables the female calculation (`FFMF*1`).
* **`A1V2` (Sampling Site)**: `1` = Arterial, `2` = Venous. 
  * *NONMEM Proof*: Mapped algebraically as [(A1V2-1)](file:///c:/home/git/simtiva/pharmacology.js#5438-5449). Applying `2` introduces venous delay penalties (`THETA(17)` mathematically adds an extra 1.419L max volume to `V1`, and `THETA(18)` reduces compartmental diffusion `Q2`).
* **`TECH` (Opiate Co-administration)**: `1` = Propofol Alone, `2` = With Opiates. 
  * *NONMEM Proof*: Bound to [(TECH-1)](file:///c:/home/git/simtiva/pharmacology.js#5438-5449) in `KV3` and `KCL`. When opiates are administered (`TECH=2`), it triggers exponential decay metrics natively annotated as *"CL declines with age with opiates"* and *"V3 declines with age with opiates"*.
* **`FFMZ` (Normalized Fat-Free Mass)**: Matches the `NFFM=FFM/FFMR` calculation inside the model, bounding the actual mass index normalized proportionally to standard healthy targets spanning `0.012` to `1.54`.

## Dataset 2: Pharmacodynamic (PD) Model Correlates
*(Source: [supplementary_digital_content_3.txt](file:///c:/home/git/simtiva/data/supplementary_digital_content_3.txt) mapped against [supplementary_digital_content_4.txt](file:///c:/home/git/simtiva/data/supplementary_digital_content_4.txt))*

The PD dataset preserves the identical core demographics (AGE, WGT, HGT, M1F2, etc.) but fundamentally shifts the dependent variables and structural tracking to map Pharmacodynamic effect (specifically, Bispectral Index or **BIS** tracking).

### Core PD Modifiers
* **`DV` (Bispectral Index / BIS)**: Unlike Dataset 1 where DV represents Propofol plasma concentration, the `DV` in the PD dataset represents the measured **BIS** score (spanning ~`0.3` to `98.0`).
* **`DVTY` (Assay Type)**: Exclusively flags `4` internally representing BIS observation rows uniquely.
* **`EVID`**: Same as PK (`0` for BIS observation, `1` for Propofol dosing).

### Individual Bayesian PK Pre-Calculations
Instead of forcing the NONMEM PD model to blindly re-evaluate population-level Pharmacokinetics purely from covariates repeatedly, the PD dataset natively **injects** patient-specific, individually tuned (Bayesian *post-hoc*) pharmacokinetic constants directly into the data matrix.

* **`EV1`**: Patient's exact central volume ($V_1$).
* **`EV2`**: Patient's exact rapid peripheral volume ($V_2$).
* **`EV3`**: Patient's exact slow peripheral volume ($V_3$).
* **`ECL`**: Patient's exact elimination clearance ($CL$).
* **`EQ2`**: Patient's exact inter-compartmental clearance to $V_2$ ($Q_2$).
* **`EQ3`**: Patient's exact inter-compartmental clearance to $V_3$ ($Q_3$).

*(NONMEM Proof: The script explicitly states `V1=EV1`, `V2=EV2`, `CL=ECL` bypassing THETA covariate algebra and linking directly to the dataset structural block).*

### Explicit PD Output
The accompanying [supplementary_digital_content_4.txt](file:///c:/home/git/simtiva/data/supplementary_digital_content_4.txt) fits the PD mathematical response structure bounding:
* **`E50`**: Concentration at 50% maximal effect (baseline ~`3.08` mg/L).
* **`KE0`**: Effect-site equilibration rate mapping plasma delay.
* **`EMAX`**: The true baseline awake BIS index explicitly derived by the model per patient (population mean ~`93.0`).
