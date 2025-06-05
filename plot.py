#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys


def do_plot(nodes, elems):
    xs = []
    ys = []
    for pt in nodes.values():
        xs.append(pt[0])
        ys.append(pt[1])
        pass
    plt.plot(xs, ys, 'bo')
    for ns in elems.values():
        plt.plot([nodes[ns[0]][0], nodes[ns[1]][0], nodes[ns[2]][0], nodes[ns[0]][0]],
                 [nodes[ns[0]][1], nodes[ns[1]][1], nodes[ns[2]][1], nodes[ns[0]][1]], 'k-')
    plt.axis([min(xs)-.5, max(xs)+.5, min(ys)-.5, max(ys)+.5])
    plt.show()
    pass


def get_nodes(f):
    nodes = {}
    line = f.readline()
    if line != '$nodes\n':
        print('incorrect format')
        pass
    line = f.readline()
    while line != '$elements\n' and line:
        nos = line.split()
        idno = int(nos[0])
        pt = [float(nos[1]), float(nos[2])]
        nodes[idno] = pt
        line = f.readline()
        pass
    return nodes


def get_elems(f):
    elems = {}
    line = f.readline()
    while line:
        nos = line.split()
        elems[int(nos[0])] = [int(nos[1]), int(nos[2]), int(nos[3])]
        line = f.readline()
        pass
    return elems


def main():
    if len(sys.argv) != 2:
        print('Usage: filename')
        sys.exit(1)
        pass
    fname = sys.argv[1]
    with open(fname) as f:
        nodes = get_nodes(f)
        elems = get_elems(f)
        pass
    do_plot(nodes, elems)  
    f.close()
    pass


if __name__ == '__main__':
    main()
    pass
