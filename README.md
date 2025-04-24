# Homework Assignment 19

![Assignment 19](https://github.com/PGE383-HPC/assignment19/actions/workflows/main.yml/badge.svg)

During a single-phase one-dimensional core flooding experiment, pressure readings were observed at 99 evenly spaced points in the interior of a core. The pressure as a function of $x$ observations are stored in the data file [data.csv](./data/data.csv) and shown in the figure below

<img src="images/data.png" width=500>

The left essential boundary condition is $p(0) = 15$ and the right essential boundary condition is $p(1) = 5$. 

In [https://github.com/PGE383-HPC/assignment18](Assignment18) we used an assumed polynomial form of the mobility function 

$$
\lambda(x) = x^{\theta_1} + \theta_2
$$

where $\theta_1$ and $\theta_2$ where unknown parameters that we were trying to learn by solving the *inverse* PDE problem with [Gripap.jl](https://gridap.github.io/Gridap.jl/stable/) and [Optimisers.jl](https://fluxml.ai/Optimisers.jl/dev/).

For this assignment, we will not assume any particular form of the mobility function, but instead use a *universal function approximator* in the form of an unknown neural network, i.e.  

$$
\lambda(x) = \mathcal{NN}_{\vec{\theta}}(x)
$$

where ${\vec{\theta}}$ are the weights and biases of a neural network.  The template code uses the Julia package [https://github.com/FluxML/Flux.jl](Flux.jl) for assigning the architecture of the neural network (do not change).  Your job is the write the training loop to learn the weights and biases.

Implement your solution in the [assignment19.jl](src/assignment19.jl) file; specifically, your job is to complete the function `train(filename::String, nepochs::Integer=5000, tolerance::Real=10, rng=MersenneTwister(1234))` where `filename` is the path to [data.csv](data/data.csv), `nepochs` is the maximum number of training iterations, `tolerance` is the loss function tolerance to stop the training iterations, and `rng` is the random number generator used to seed the inital weights and biases of the network.  The loss function has already been provided in the function `loss()`.


## Testing

To see if your answer is correct, run the following command at the Terminal
command line from the repository's root directory

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

the tests will run and report if passing or failing.
