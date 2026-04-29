# Project Implementation

## Context

The MC simulation I am trying to implement model the behaviour of a dense melt of ring polymers with sparse linear polymers in solution. The linear polymers have reactive end and can merge together forming longer reactive linear polymers. In this process the ring polymers can be linked together forming a cluster of linked rings. I am interested in simulating the network evolution of this system.

## DSMC Implementation

Use the previous implementation of the DSMC algorithm (in the /DSMC folder) to simulate the gelation process of the linear component of the melt by transforming this part of the code to python. This would mean that you will have an array of linear polymer lengths that will evolve with those rules.

## Network Growth Implementation

The network growth is a poissonian process, meaning that at every cyclisation event as modelled by the DSMC algorithm, there is a probability that the newly formed ring is linked to one or more rings in the system. This number follows a poisson distribution with given average. This average is modelled from MD simulations and given in this code as a function of both the linear and ring polymer lengths. So the linear will link a number of rings drawn from that poisson distribution. 

## Stages implementation

Follow the structure of the previous implemenation in simulate_network_growth_progressive.py but make the following changes.

- After one cyclisation the linear becomes part of the ring list hence can get linked by the linear polymers
- mlin new linear polymers are added at each stage and the DSMC controlled ligation process follows until all the linear polymers have cyclisised.
- in this polydisperse system each population have a probability of being selected to be linked to. This, similarly to what discussed before must be a function of the concentration, ring length and linear length that I still have to work out. So for now assume that the returned number is indeed drawn from poissonian with certain average and I'll come back to you in the future with the actual function.