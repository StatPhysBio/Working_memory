# Working_memory


This repository contains the code and analysis associated with the manuscript:

Schnaack, Peliti, and Nourmohammad, [Learning and organization of memory for evolving patterns, arXiv:2106.02186](https://arxiv.org/abs/2106.02186).

It allows reproduction and a simplified visualization of the numerical results reported in the manuscript.

## Dependencies

The code is written in [Julia](https://julialang.org) and depends on several packages. Below we give a list of version numbers of the packages for which the code is known to run.
- Julia 1.5.3
- Distributions 0.24.6
- SpecialFunctions 1.1.0
- QuadGK 2.4.1
- FileIO 1.4.5
- IJulia 1.23.1
- Plots 1.9.1

## Usage

The repository includes the source code to generate and analyze data as well as jupiter notebooks that can be used to produce all figures used in the manuscript. To generate your own figures, you can, for example, install [IJulia](https://github.com/JuliaLang/IJulia.jl) and run
```bash
using IJulia
notebook()
```
in your Julia terminal to launch the IJulia notebook in your browser. Navigate to the directory of the notebook folder ```bash  Working_memory/Notebooks_Hopfield_vs_Compartments ```. Then select one of the notebooks to reproduce the figures of the manuscript.
In order to reduce the computation time, the default simulation parameters are set to fewer repetitions and a smaller resolution compared to the manuscript, but the results capture the reported behavior. Further, data for the parameters set in the notebooks is already stored in ```bash  Working_memory/Data ```. The system will load this data and perform the analysis and produce the corresponding figures. To produce your own data eigther remove the corresponding data files or change the simulation parameters.
To get comparable results set the parameters to the values given in the manuscript. Note: As most of the simulations are stochastic you generally do not expect precisely equivalent plots.


## Notebooks

There are three notebooks in the folder each will produce figures that show the analysis presented in the manuscript but using smaller systems.

### Figure_2_and_3.ipynb

Calculates the optimal learning rate and optimal perfromance for a distributed memory (C=1) and a fully compartemtelized 1-to-1 strategy (C=N). These results are similar to figures 2 and 3 in the manuscript.

###

## Contact

Any issues or questions should be addressed to [me](mailto:oskar.schnaack@ds.mpg.de).
