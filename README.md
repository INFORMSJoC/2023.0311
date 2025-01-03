# Activated Benders Decomposition for Paratransit Itinerary Planning

This repository supplements the above mentioned paper, published in INFORMS Journal of Computing. 

## Authors

This codebase was developed by Kayla Kummerlowe (née Cummings) (kaylakummerlowe \<at\> gmail \<dot\> com).

Coauthors A Jacquillat and V Vaze contributed to models, methodology, and experiment design underlying the code.

## Utilization

- data: Simulate trips and build a SIPPAR case study - see dedicated README in that directory.
- src: Package containing models and algorithms from the paper. 
- notebooks: Jupyter notebook with examples of how to execute the main algorithms and acceleration strategies.
- results: Synthetic results from Appendix C.1 Table 9. To replicate this Table, run notebooks/solve.ipynb. The output is also in ./results/summary.csv. Further instructions are below.

## Installation

This repo includes code written in [Python](https://www.python.org/downloads/), [R](https://cran.r-project.org/mirrors.html), and [Julia](https://julialang.org/downloads/), using [Gurobi](https://www.gurobi.com/downloads/gurobi-software/#GO1003) and [Jupyter notebook](https://jupyter.org/) software.

### Julia and Gurobi Installation Notes 

The model and algorithm implementations are based on the JuMP package in the Julia programming language, as well as the Gurobi Optimizer. At the time of writing (January 2, 2025), the JuMP package does not yet support the v11 or v12 solver, so v9 or [v10](https://www.gurobi.com/downloads/gurobi-software/#GO1003) must be installed.


1. Install [Julia](https://julialang.org/downloads/).
2. Obtain a [Gurobi license](https://www.gurobi.com/free-trial/).
3. Download [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-software/#GO1003) software.
4. Install the [Julia wrapper for Gurobi](https://github.com/jump-dev/Gurobi.jl) - their README is a valuable resource for proper installation, to ensure the license is visible to the Optimizer.
5. Install other Julia packages:
  - JuMP
  - ArgParse
  - CSV
  - DataFrames
  - Dates
  - Glob
  - JSON
  - Printf
  - StatsBase