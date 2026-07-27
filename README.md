# Connectivity

## Run the Julia movement-gradient experiment

Open the project root (this folder) in Positron, then select the project Julia
environment when prompted. In the Julia console, install the dependencies once:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

To work interactively, open `src/movement_gradient_experiments.jl` and run the
file. This loads all functions without automatically starting the lengthy full
analysis. Start it explicitly in the Julia console:

```julia
main()
```

Alternatively, run the complete analysis from a terminal at the project root:

```sh
julia --project=. src/movement_gradient_experiments.jl
```

Progress is printed for each of the six scans. Results are written to
`src/movement_gradient_output/`.
