# Fugacity Analysis of a Water–Ammonia Binary System

> **Note:** This project accompanies a detailed technical report with experimental data, derivations, and theoretical background.

## Project Summary

This project investigates the non-ideal vapor-phase behavior of a binary water–ammonia system using the Peng–Robinson equation of state (PR EOS). Fugacity coefficients and fugacities are computed across varying temperature, pressure, and composition conditions, with results compared against generalized correlations and the Lewis fugacity rule.

The study emphasizes both **rigorous thermodynamic modeling** and **model validation**, highlighting when simplifying assumptions break down under non-ideal conditions.

---

## Thermodynamic Framework

- Equation of State: Peng–Robinson (PR EOS)
- Phase considered: Vapor phase
- Mixing rule: Quadratic mixing with binary interaction parameter
- Binary interaction parameter:
  - Temperature-dependent
  - Derived from published experimental data (2012 study)

Pure-component parameters are calculated using critical properties and acentric factors for water and ammonia.

---

## Modeling Approach

### 1. Binary Interaction Parameter Analysis

- Computation of the cross-parameter `a₁₂(T)` with and without the binary interaction correction
- Visualization of temperature dependence
- Comparison illustrating the impact of non-ideal interactions

### 2. Fugacity Coefficient Calculations

Fugacity coefficients are calculated using the PR EOS as functions of:
- Pressure (constant temperature and composition)
- Temperature (constant pressure and composition)
- Vapor-phase composition (constant temperature and pressure)

The cubic EOS is solved numerically at each condition, with the vapor-phase compressibility factor selected for fugacity calculations.

### 3. Generalized Correlation Comparison

- Fugacity coefficients are independently estimated using generalized correlations
- Double interpolation is performed in reduced temperature and pressure space
- Results are used to assess the validity of simplified correlation-based methods

### 4. Lewis Fugacity Rule Evaluation

- Fugacities computed using:
  - Rigorous PR EOS mixture calculations
  - Lewis-rule approximations with generalized correlations
- Percent error quantified for:
  - Near-ideal conditions
  - Strongly non-ideal conditions

---

## Results & Insights

- Binary interaction effects significantly influence predicted fugacity behavior
- Lewis fugacity rule performs well under near-ideal conditions
- Errors increase substantially at elevated pressures, highlighting limits of ideal-mixture assumptions
- PR EOS provides consistent and physically reasonable predictions across the full parameter space

Plots and numerical outputs are generated automatically by the MATLAB script.

---

## Numerical & Computational Methods

- Cubic equation of state root solving
- Vapor-phase root selection
- Parametric sweeps over temperature, pressure, and composition
- Double interpolation for tabulated correlation data
- Error quantification and comparison metrics

---

## Tools & Technologies

- MATLAB
- Numerical root finding
- Thermodynamic property estimation
- Data visualization and parametric analysis

---

## Documentation

A full technical report accompanies this codebase and includes:
- Theoretical background
- Equation derivations
- Assumptions and limitations
- Literature references
- Detailed discussion of results

---

## Portfolio Relevance

This project demonstrates:
- Applied chemical thermodynamics
- EOS-based mixture modeling
- Validation of simplifying assumptions
- Research-oriented coding and documentation

It reflects skills directly applicable to process simulation, separations, energy systems, and graduate-level chemical engineering work.
