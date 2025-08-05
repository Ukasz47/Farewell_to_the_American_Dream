# Farewell to the American Dream – How Endogenous Feedback Loops Shape Wealth Accumulation
**Master Thesis Code Guide**

## Code Overview

The replication package consists of two main parts:

1. **Part 1** solves the stylized version of the model, used to generate Figures 1 and 2 from the thesis.  
   These provide basic intuition about the model’s core mechanisms.  
   - Run `Figure1.jl` and `Figure2.jl` to reproduce the corresponding plots.  
   - The file `functions_stylized.jl` contains the main functions for solving the model presented in Section 3.1 using the Endogenous Grid Method.

2. **Part 2** solves the full model, used to generate Figure 3 and the results presented in Tables 1–7.  
   - Run `Model_full.jl` for the baseline model described in Section 3.2.  
   - Files `Model_homothetic.jl` and `Model_fixed_phi.jl` generate results for alternative model specifications (Section 5.3).  
   - `functions.jl` contains the main EGM-based functions for solving the models.  
   - `analysis.jl` includes auxiliary routines for post-processing and analysis.

Both parts use shared files:  
- `toolbox.jl` – general-purpose functions  
- `plotting.jl` – generates the visual outputs

---

## Requirements and Packages

The code was written in **Julia** and relies on publicly available packages for numerical computation and visualization.  
Core packages include:

- `Interpolations.jl`
- `Plots.jl`
- `FastGaussQuadrature.jl`
- `Parameters.jl`
- `LinearAlgebra`
- `Distributions.jl`
- `Roots.jl`
- `SparseArrays`
- `Measures.jl`
- `StatsPlots.jl`
- `ProgressMeter.jl`

All packages are open source and can be installed via the Julia package manager:

```julia
import Pkg
Pkg.add([
    "Interpolations",
    "Plots",
    "FastGaussQuadrature",
    "Parameters",
    "LinearAlgebra",
    "Distributions",
    "Roots",
    "SparseArrays",
    "Measures",
    "StatsPlots",
    "ProgressMeter"
])
