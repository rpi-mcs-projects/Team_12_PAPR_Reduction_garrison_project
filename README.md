# PAPR Reduction for OFDM Under PA Nonlinearity (MATLAB)

This project simulates **PAPR reduction methods for OFDM** in the presence of a **nonlinear power amplifier (PA)** and evaluates performance using:

- **CCDF(PAPR)** (pre-PA)
- **EVM and BER** vs **SNR** across **PA output back-off (BO)**
- **ACLR** vs **BO**
- **Target-based view**: EVM improvement vs BASE at **BO = 0 dB** for **SNR = 20/25 dB**
- **Sensitivity**: SLM **U=2 vs 4** and PTS **B=2 vs 4** at **BO=0 dB, SNR=20 dB**

All figures and results are saved to `./results`.

---

## Repository / File Overview

Typical structure:
├── mid_term_check.m # main entry point (single-file project)
└── results/ # auto-created output folder (PNG/CSV/MAT)


The main script/function is:
- `garrison_final_proj.m` (contains all helper functions at the bottom)

---

## Requirements

### MATLAB
- MATLAB **R2021a or newer** recommended (earlier versions may work, but not guaranteed)

### Required Toolboxes
This code uses the following MATLAB functions that require toolboxes:
- **Communications Toolbox**
  - `qammod`, `qamdemod`, `awgn`, `de2bi`
- **Signal Processing Toolbox**
  - `pwelch`, `hann` (for ACLR/PSD)

If you get “undefined function” errors for any of the above, install the corresponding toolbox.

### Hardware
- **No external hardware required.**
- Simulation runs entirely in MATLAB.

### Dataset
- **No dataset required.**
- OFDM bits are generated internally using `randi(...)` (reproducible via `rng(...)`).

---

## How to Run

1. Open MATLAB and set the **Current Folder** to the project directory containing `mid_term_check.m`.

2. Run the main function from the Command Window:
```matlab
garrison_final_proj
