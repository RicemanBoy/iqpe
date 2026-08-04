import time
import gates as gate
from utils import create_tree, save_tree

#basis = [gate.su2(gate.H), gate.su2(gate.T)]
basis = [gate.Toffoli, gate.H1(), gate.H2(), gate.H3()]
depth = 10

t = time.time()
tree = create_tree(basis, max_depth=depth)
print('Elapsed time:', time.time() - t)
print('Number of gates:', len(tree['names']))

save_tree(tree, 'SolovayKitaev/trees/HToff_{}.pkl'.format(depth))
