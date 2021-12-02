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

The repository includes the source code to generate and analyze data as well as Jupiter notebooks that can be used to produce all figures used in the manuscript. To generate your own figures, you can, for example, install [IJulia](https://github.com/JuliaLang/IJulia.jl) and run
```bash
using IJulia
notebook()
```
in your Julia terminal to launch the IJulia notebook in your browser. Navigate to the directory of the notebook folder ```bash  Working_memory/Notebooks_Hopfield_vs_Compartments ```. Then select one of the notebooks to reproduce the figures of the manuscript.
In order to reduce the computation time, the default simulation parameters are set to fewer repetitions and a smaller resolution compared to the manuscript, but the results capture the reported behavior. Further, data for the parameters set in the notebooks is already stored in ```bash  Working_memory/Data ```. The system will load this data and perform the analysis and produce the corresponding figures. To produce your own data either remove the corresponding data files or change the simulation parameters.
To get comparable results set the parameters to the values given in the manuscript. Note: As most of the simulations are stochastic you generally do not expect precisely equivalent plots.


## Notebooks

There are three notebooks in the folder.  Each notebook will produce figures that show the analysis presented in the manuscript but using smaller systems.


please be aware, that we use the learning rate λ as well as the memory retention rate (γ = 1 - λ) throughout these notebooks and the source code. Further, we might use x as being equivalent to the mutation rate μ.

### Figure_2_and_3.ipynb

Calculates the optimal learning rate and optimal performance for a distributed memory (C=1) and a fully compartmentalized 1-to-1 strategy (C=N). These results are similar to figures 2 and 3 in the manuscript.

### Temperature_transitions.ipynb

Shows the dependence of the performance on the inverse temperature <img src="https://render.githubusercontent.com/render/math?math=\beta_{\rm H}">. And the correct allocation of the patterns to the compartments as a function of the inverse temperature <img src="https://render.githubusercontent.com/render/math?math=\beta_{\rm S}">.
We assume a constant learning rate during these simulations, and to calculate the corresponding heat maps in figure 4 of the manuscript an optimization for each parameter is needed. While this is not done in this reduced analysis, the results do show the transitions presented in the manuscript.

###  Plots_SI_reconstructed_vs_not_reconstructed_patterns.ipynb
Presents the analysis to produce Figures SXX,XXX,SXX.

The mode of data production relies on the existence of optimizations for the chosen parameters. Alternatively, it is possible to set the optimal learning rate λ by hand.


## Source code

the Julia code for the simulation and analysis presented in the manuscript is split into 4 files which are located in the directory ```bash Working_memory/src```

### Hopdield-model.jl

The functions defines in this file constitute the basic dynamics of the extended Hopfield model we introduce in the manuscript. Among others, these include the initialization of patterns, the interaction matrix, the mutations of patterns, and continuous updates of the interactions.

### Analysis-tools.jl

Gives functions that are used for the data analysis such as the protocol to find the optimal learning rate from performance on a learning rate grid.


### Self_recognition-Compartments.jl


Gives the necessary functions to run simulations. Such as the protocols for data products that include the initialization of patterns and networks and the reconstruction of evolved patterns.


### Self_recognition-Simulations.jl


The functions prested here can be used to produce data that is stored in the directory ```bash  Working_memory/Data```. These functions call the simulations from the file ```bash Self_recognition-Compartments.jl``` and then store the data. Before new data is produced the system will always check whether this data has already been produced. If you wish to produce new data with the same parameters, please delete or move the old data from its original storage.


## Contact

Any issues or questions should be addressed to [me](mailto:oskar.schnaack@ds.mpg.de).
