# HW0

## Getting Started
All the homeworks are distributed as Jupyter notebooks. Follow these instructions to get everything set up.

1. Install Julia. Download v1.6.5 from the [Julia website](https://julialang.org/downloads/#long_term_support_release). Extract the binaries onto a stable location on your computer, and add the `/bin` folder to your system path. Alternatively you can follow these [platform specific instructions](https://julialang.org/downloads/platform/). 
2. Clone this repo, and open a terminal in the root directory
2. Start a Julia REPL in your terminal using `julia`. This should work as long as the binaries are on your system path.
3. Install the [IJulia](https://github.com/JuliaLang/IJulia.jl) using the Julia package manager. In the REPL, enter the package manager using `]`, then enter `add IJulia` to add it to your system path.
4. In the REPL (hit backspace to exit the package manager), enter `using IJulia`
5. Launch the notebook using `notebook()` or `jupyterlab()`

## Updating `studentinfo()`
Update your name and Andrew ID in the file `src/HW0.jl`. 

## Doing the homework 
Fill in `src/Q1.ipynb` and `src/Q2.ipynb`. There are tests in these notebooks that ensure you have completed the exercises correctly. If all of the tests pass, you can consider the assignment completed. 

## Submitting your homework
Make sure your repo lives under the Class Organization. This will be done automatically when you use the GitHub Classrooms link we send provide. To submit, simply create a release. Follow [these instructions](https://github.com/Optimal-Control-16-745/JuliaIntro/blob/main/docs/Submission%20Instructions.md) for more info on creating the release.

## Adding the Upstream Repo
We may release changes to the homework periodically if errors or bugs are found. Follow these instructions for linking your repo to the original template and pulling changes. It's always a good idea to branch your code before pulling from the upstream repo in case something goes wrong or the merge is particularly nasty. Do the right away after creating your repo.
```
git remote add upstream https://github.com/Optimal-Control-16-745/HW0
git pull upstream main --allow-unrelated-histories
```
