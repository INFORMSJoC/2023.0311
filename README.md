[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Activated Benders Decomposition for Paratransit Itinerary Planning

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Activated Benders Decomposition for Day-ahead Paratransit Planning](https://doi.org/10.1287/ijoc.2023.0311) by Kayla Kummerlowe née Cummings, Alexandre Jacquillat, and Vikrant Vaze. This codebase was developed by Kummerlowe. Jacquillat and Vaze contributed to models, methodology, and experiment design underlying the code.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2023.0311

https://doi.org/10.1287/ijoc.2023.0311.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{ActivatedBendersDecomposition,
  author =        {K. Cummings and A. Jacquillat and V. Vaze},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Activated Benders Decomposition for Day-ahead Paratransit Itinerary Planning}},
  year =          {2025},
  doi =           {10.1287/ijoc.2023.0311.cd},
  url =           {https://github.com/INFORMSJoC/2023.0311},
  note =          {Available for download at \url{https://github.com/INFORMSJoC/2023.0311}}
}  
```

## Description

The purpose of this repository is to build itineraries for employee paratransit drivers based on trip requests finalized on the day prior to service. This project also provides Benders-based decomposition algorithms with acceleration methods to solve the underlying models.

## Building

This repo includes code written in [Python](https://www.python.org/downloads/), [R](https://cran.r-project.org/mirrors.html), and [Julia](https://julialang.org/downloads/), using [Gurobi](https://www.gurobi.com/downloads/gurobi-software) and [Jupyter notebook](https://jupyter.org/) software.

### Julia and Gurobi Installation Notes 

The model and algorithm implementations are based on the JuMP package in the Julia programming language, as well as the Gurobi Optimizer. 

1. Install [Julia](https://julialang.org/downloads/).
2. Obtain a [Gurobi license](https://www.gurobi.com/free-trial/).
3. Download [Gurobi Optimizer](https://www.gurobi.com/downloads/gurobi-software) software.
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

### Execution

To build the case study, follow the directions in data/README.md.

To build the model and solve with Activated Benders Decomposition, and to obtain Appendix C.1 Table 9, execute the cells in scripts/solve.ipynb.

## Results

Synthetic results from Appendix C.1 Table 9 are available in results/summary.csv. 

## Replicating

To generate the results in results/summary.csv, do the following:

* Set seed "4747" in data/simulate_rides.R.
* Set the following parameters in data/configure_parameters.py:
  - num_trips - 30
  - hourly_wage - 20
  - penalty - 100
  - num_scenarios - 100
  - num_itin - 10000

Then follow the steps under Execution header.

Based on the specs of your machine, the exact solve times will vary.

## Support

For support in using this software, submit an [issue](https://github.com/INFORMSJoC/2023.0311/issues/new).