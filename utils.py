from math import sqrt
import matplotlib.pyplot as plt

#https://stackoverflow.com/questions/13214809/pretty-print-2d-python-list
def pretty_print(mat):
    s = [[str(e) for e in row] for row in mat]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))

def euclidean_distance(x1, y1, x2, y2):
    return sqrt((x1-x2)**2 + (y1-y2)**2)

# check if elements of a list are all unique
def no_rep(list):
    visited = set()
    for node in list:
        if node not in visited:
            visited.add(node)
        else:
            return False
    return True

def swap(list, a, b):
    tmp = list[a]
    list[a] = list[b]
    list[b] = tmp

def gen_plot(list, fname):
    fig = plt.plot(list)
    plt.ylabel("distance")
    plt.ylim(0, 1200000)
    plt.xlim(0, len(list))
    plt.xlabel("iterations")
    plt.savefig(fname + '.png')
    plt.show()

def write_progress(list, fname):
    with open(fname + '.txt', 'w+') as f:
        for item in list:
            f.write(str(item) + ',')
