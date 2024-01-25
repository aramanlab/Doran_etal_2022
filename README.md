# Doran_etal_2022

*current state: prerelease, This will be curated into a more interpretable and user friendly version soon.*

**See [https://github.com/aramanlab/Doran_etal_2023](https://github.com/aramanlab/Doran_etal_2023) for the curated version**

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project for the paper

> Personalized signatures of human gut bacteria discovered from proteome coevolution
> Doran et al.
>, ((paper year)) 

To (locally) reproduce this project, do the following:

1. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
   ```
   git clone git@github.com:aramanlab/Doran_etal_2022.git
   cd Doran_etal_2022
   git submodule init
   git submodule update
   ```
2. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
