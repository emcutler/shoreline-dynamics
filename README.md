# shoreline-dynamics
Approximates output from event-driven lifecycle shoreline management models

Reference for accompanying paper:
Cutler, E.M.; Albert, M.R. and White, K.D. A Low-cost shoreline dynamic model with sea level rise and storm impacts. Journal of Coastal Research. In review

The model can be run from main.m
shorelineDynamics.m performs the calculations for changes in shoreline position
parameters.m contains parameter values tuned to three locations in Florida: Hutchinson Island, Vilano Beach, and Gasparilla Island under three future sea level rise scenarios
sensitivityAnalysis.m displays output from global sensitivity analalysis for these three locations

The model outputs mean nourishment interval necessary to maintain a specified level of risk reduction and first order sensitivity scores for uncertain parameters. To store these outputs, it is necessary to have the following directories in the home directory: "~/sensiOutputs/", "~/outputs/"
