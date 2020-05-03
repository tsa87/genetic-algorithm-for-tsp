##################################################
## Author: Tony Shen
##################################################

from utils import *
from random import randrange, uniform, sample, choice
from statistics import mean
from copy import deepcopy
from argparse import ArgumentParser
import pickle
import numpy as np

# traveling salesperson problem
class TSP:
    def __init__(self, file_path):
        self.cord_mat = []
        self.dist_mat = [[]]
        self.length = 0
        self.load_file(file_path)

    def load_file(self, file_path):
        with open(file_path) as f:
            for line in f:
                line = line.split()
                self.cord_mat.append((int(line[1]), int(line[2])))

        self.length = len(self.cord_mat)

        # compute distance matrix
        self.dist_mat = [[0 for j in range(self.length)] for i in range(self.length)]
        for i in range(self.length):
            for j in range(i, self.length):
                dist = euclidean_distance(
                    self.cord_mat[i][0],
                    self.cord_mat[i][1],
                    self.cord_mat[j][0],
                    self.cord_mat[j][1]
                )
                self.dist_mat[i][j] = dist
                self.dist_mat[j][i] = dist


class GA:
    def __init__(
        self,
        tsp,                        # traveling salesperson problem instance
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
        epoch_limit=20,
        output=False
    ):

        self.tsp = tsp
        self.organisms = []

        #utils
        self.output=output
        self.min_history = []
        self.max_history = []
        self.avg_history = []
        self.epoch = 0
        self.epoch_limit = epoch_limit


        #config
        self.size = size
        self.growth_factor = growth_factor
        self.mutation_factor = mutation_factor
        self.mutation_scheduler = mutation_scheduler

        #probability
        self.pmx_p=pmx_p
        self.elitism=elitism
        self.live_learn_p=live_learn_p
        self.point_mutation_p=point_mutation_p
        self.frame_mutation_p=frame_mutation_p

        if greedy_init:
            self.greedy_init()
        else:
            self.normal_init()

        scores = [self.fitness_level(node) for node in self.organisms]
        self.display_stats(scores)

    def __str__(self):
        return f"s={self.size}gf={self.growth_factor}epoch={self.epoch}"

    def evolve(self):
        while self.epoch < self.epoch_limit:
            self.run_epoch()

        if self.output:
            self.write_result()

    def write_result(self):
        write_progress(self.min_history, 'min')
        write_progress(self.max_history, 'max')
        write_progress(self.avg_history, 'avg')
        with open("solution.txt", 'w+') as f:
            levels = [self.fitness_level(organism) for organism in self.organisms]
            i = levels.index(min(levels))
            for city in self.organisms[i]:
                f.write(str(city) + " ")

    def run_epoch(self):
        self.off_springs = []

        for _ in range(round(self.growth_factor*self.size)):
        # reproduce
            fitness = [1/self.fitness_level(organism) for organism in self.organisms]
            total = sum(fitness)
            select_p = [score/total for score in fitness]
            organism = self.organisms[np.random.choice(self.size, p=select_p)]

            p = uniform(0, 1)
            if p < self.elitism:
                off_spring = organism
            else:
                p = uniform(0, 1)
                if p < self.pmx_p:
                    partner = choice(self.organisms)
                    off_spring = self.partially_mapped_crossover(organism, partner)
                else:
                    off_spring = self.live_learn(organism)

            self.off_springs.append(off_spring)

        self.off_springs.sort(key = lambda off_spring : self.fitness_level(off_spring))
        self.organisms = deepcopy(self.off_springs[:self.size])

        # mutate
        p = uniform(0, 1)
        for organism in self.organisms:
            if p < self.point_mutation_p:
                self.point_mutation(organism)
            if p < self.frame_mutation_p:
                self.frame_mutation(organism)

        scores = [self.fitness_level(node) for node in self.organisms]

        min_score = round(min(scores))
        high_score = round(max(scores))
        avg_score = round(mean(scores))
        self.min_history.append(min_score)
        self.max_history.append(high_score)
        self.avg_history.append(avg_score)

        self.display_stats(scores)

        if self.mutation_scheduler:
            self.mutation_factor = max(round((min_score-600000)/100000), 0) + 1

        self.epoch += 1


    def point_mutation(self, organism):
        for _ in range(self.mutation_factor):
            # pick a gene that is most costly

            dists = [0 for _ in range(self.tsp.length)]

            for i in range(self.tsp.length-1):
                cost = self.tsp.dist_mat[organism[i]][organism[i+1]]
                dists[i] += cost
                dists[i+1] += cost

            cost = self.tsp.dist_mat[organism[-1]][organism[0]]
            dists[-1] += cost
            dists[0] += cost

            total = sum(dists)
            select_p = [dist/total for dist in dists]

            a = np.random.choice(self.tsp.length, p=select_p)
            b = randrange(0, tsp.length)
            swap(organism, a, b)

    def frame_mutation(self, organism):
        start = randrange(0, tsp.length)
        end = randrange(0, tsp.length)

        if end >= start:
            segment = organism[start:end]

            # select a point that is not in [start, end)
            point = randrange(0, tsp.length-(start-end))
            if point < start: #insertion point before segment
                organism = organism[:point] + segment + organism[point:start] + organism[end:]
            else: #insertion point behind segment
                point = point - start + end
                organism = organism[:start] + organism[end:point] + segment + organism[point:]

        #wrap around segment
        else:
            segment = organism[start:] + organism[:end]
            point = randrange(end, start)
            organism = organism[end:point] + segment + organism[point:start]


    # http://user.ceng.metu.edu.tr/~ucoluk/research/publications/tspnew.pdf
    def partially_mapped_crossover(self, organism1, organism2):
        point = randrange(0, self.tsp.length)   # off_spring perserves [0, point) of organism1
        off_spring = deepcopy(organism2)
        for i in range(point):
            edit = organism1[i]
            j = off_spring.index(edit)
            swap(off_spring, i, j)
        return off_spring


    def live_learn(self, organism):
        point = randrange(1, self.tsp.length)   # off_spring perserves [0, point) of organism

        off_spring = deepcopy(organism[:point])
        visited = set(organism[:point])

        curr_node = off_spring[-1]
        for _ in range(self.tsp.length - point):
            curr_node = self.greedy_next(curr_node, visited)
            off_spring.append(curr_node)

        return off_spring

    def greedy_init(self):
        for _ in range(self.size):
            organism = [randrange(0, self.tsp.length)]
            visited = set(organism)

            curr_node = organism[-1]
            for _ in range(tsp.length-1):
                curr_node = self.greedy_next(curr_node, visited)
                organism.append(curr_node)

            self.organisms.append(organism)


    def greedy_next(self, curr_node, visited):
        neighbor_dist = self.tsp.dist_mat[curr_node]
        sorted_indices = np.argsort(neighbor_dist)

        for node in sorted_indices:
            if node not in visited:
                visited.add(node)
                return node

    def normal_init(self):
        genes = [i for i in range(self.tsp.length)]
        for _ in range(self.size):
            organism = sample(genes, self.tsp.length)
            self.organisms.append(organism)

    def fitness_level(self, organism, debug=False):
        if debug:
            assert(len(organism) == self.tsp.length)
            assert(no_rep(organism))

        total_dist = 0
        for i in range(self.tsp.length-1):
            total_dist += self.tsp.dist_mat[organism[i]][organism[i+1]]
        total_dist += self.tsp.dist_mat[organism[-1]][organism[0]]

        return total_dist


    def display_stats(self, scores, display_all=False):
        print(f"\nEpoch: #{self.epoch}")

        if display_all:
            for score in scores:
                print(round(score))

        print(f"Min distance: {round(min(scores))}")
        print(f"Max distance: {round(max(scores))}")
        print(f"Avg distance: {round(mean(scores))}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--file', '-f', required=True)
    args = parser.parse_args()

    tsp = TSP(args.file)
    ga = GA(tsp)
    ga.evolve()
