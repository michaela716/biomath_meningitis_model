# biomath_meningitis_model
# Bacterial meningitis model

This project models bacteria (B) and the immune response (IL-6, I) in
the cerebrospinal fluid using two differential equations.

Authors: [Michaela Willis, Ashlee Pretorius]

## Files

- mcalculations.m : sets the parameters, finds the equilibria, and checks their stability
- msimulations.m : [one sentence on what this file does, e.g. simulates B and I over time and plots them]

## How to run

1. Download all files from this page (green Code button, then Download ZIP) and unzip them.
2. Open MATLAB and go to the unzipped folder.
3. Run mcalculations.m, then run msimulations.m.

## Requirements

MATLAB R2019b or newer. No extra toolboxes needed.

## Expected output

- Calibrated killing rate k is about 0.0042
- R0 is about 17.21
- The disease-free equilibrium is unstable
- The persistent equilibrium is stable (a spiral)
