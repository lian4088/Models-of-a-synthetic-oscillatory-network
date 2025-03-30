
# Synthetic Oscillatory Network: Repressilator Simulation

This project simulates the dynamics of a synthetic gene regulatory network called the **Repressilator**, originally proposed by Elowitz and Leibler (2000), using a deterministic system of ordinary differential equations (ODEs).

The Repressilator is a synthetic cyclic network composed of three transcriptional repressors — LacI, TetR, and λ cI — each repressing the next in a loop. This negative feedback system generates periodic oscillations in both mRNA and protein concentrations, forming a simple synthetic biological clock.

## 📁 Project Structure

```
repressilator-oscillation/
│
├── project code.py               # Python code for ODE simulation and plotting
├── oscillation analysis project.pdf  # Project report describing model construction and interpretation
├── 35002125.pdf                  # Original Nature article (Elowitz & Leibler, 2000)
├── README.md                     # This file
```

## 🔬 Model Overview

- The network consists of three coupled pairs of mRNA-protein ODEs (6 total).
- The model uses:
  - Hill repression function for transcription
  - First-order degradation for both mRNA and proteins
  - Empirical parameters drawn from biological estimates
- Simulation conducted using `scipy.integrate.odeint` in Python.

## 🧪 Simulation Highlights

- Initial conditions set for symmetrical or asymmetrical analysis
- Oscillations observed in both mRNA and protein concentrations
- Time series visualizations of LacI, TetR, and cI regulators

## 📊 Key Parameters

| Parameter      | Value                        |
|----------------|------------------------------|
| Transcription rate (α) | 0.5 transcripts/second |
| Baseline transcription (α₀) | 0.0005 transcripts/second |
| Hill coefficient (n)   | 2 |
| Translation efficiency (β) | 20 proteins/mRNA |
| Km (half-max repressor) | 40 monomers/cell |
| mRNA half-life | 2 min |
| Protein half-life | 10 min |

## 📎 Citation

If you use this project or adapt it for your work, please cite:
```
Lian, Lian (2024). Modeling a Synthetic Oscillatory Network Based on the Repressilator. GitHub Repository.
Original study: Elowitz, M. B., & Leibler, S. (2000). Nature, 403(6767), 335–338.
```

## 🧠 Insights

- Demonstrates how even a minimal genetic network can generate oscillations.
- Showcases how symmetry in initial conditions can suppress oscillatory behavior.
- Provides a reproducible Python implementation of a classic synthetic biology model.

## 📬 Contact

For questions, feedback, or collaboration, reach out to: xlian8@wisc.edu
