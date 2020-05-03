To learn more details about the features, check out "Solving Travelling Salesman Problem with Genetic Algorithms.pdf"
To run the genetic algorithm yourself, run <<< python3 tsp.py -f cities1000.txt

The best solution is generated with 20 iterations using the following default configuration for GA (Genetic Algorithm instance) in 6 minutes.
size=100,                   # size of the population
growth_factor=1.5,          # number of offsprings per organism
point_mutation_p=0.1,       # probability of point_mutation
frame_mutation_p=0.1,       # probability of frame_mutation
elitism=0.05,               # probability of suriving without reproduction
pmx_p=0.5,                  # probability of partially_mapped_crossover
live_learn_p=0.5,           # probability of live and learn reproduction
mutation_factor=1,          # number of basepairs mutated in point_mutation
mutation_scheduler=True,    # adapt mutation level based on fitness level.
greedy_init=False,          # initalize population with greedy solutions of tsp
epoch_limit=20

More infomation about the above configuartions can be found in the accompanying paper.
Thank you.
