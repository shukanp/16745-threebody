
# HW3: Intro to Nonlinear Trajectory Optimization 
## Due date: Friday, March 26

In this homework you will implement some foundational methods for solving simple nonlinear trajectory optimization problems. Here's an overview of the problems:
1. Implement iLQR to get a planar quadrotor to do a simple flip. Track with TVLQR to make it robust to wind.
2. Write a sequential-quadratic-programming (SQP) solver to find a solution to the canonical cartpole problem. Generate a TVLQR controller to make it robust to model mismatch.

For Q2, all your solutions can be implemented in regular Julia files. This will allow you to use any Julia IDE you want (if you prefer to not use Jupyter). Make sure your Jupyter notebook still runs and passes the tests before you submit, though, since that is what we'll be checking with the autograder. 

## Submitting your homework
Make sure your repo lives under the Class Organization. This will be done automatically when you use the GitHub Classrooms link we send provide.

1. Create a release. Follow [these instructions](https://github.com/Optimal-Control-16-745/JuliaIntro/blob/main/docs/Submission%20Instructions.md) for more info on creating the release.

## Adding the Upstream Repo
We may release changes to the homework periodically if errors or bugs are found. Follow these instructions for linking your repo to the original template and pulling changes. It's always a good idea to branch your code before pulling from the upstream repo in case something goes wrong or the merge is particularly nasty. Do the right away after creating your repo. 
```
git remote add upstream https://github.com/Optimal-Control-16-745/hw3
git pull upstream main --allow-unrelated-histories
```
