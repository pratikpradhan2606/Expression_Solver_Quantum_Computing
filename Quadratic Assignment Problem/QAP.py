
"""### **Quadartic Assignment Problem (QAP)**"""

pip install pyqubo

from pyqubo import Array, Placeholder, Constraint
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

n_city = 4
x = Array.create('c', (n_city, n_city), 'BINARY')

# Constraint not to visit more than two cities at the same time.
time_const = 0.0
for i in range(n_city):
    # If you wrap the hamiltonian by Const(...), this part is recognized as constraint
    time_const += Constraint((sum(x[i, j] for j in range(n_city)) - 1)**2, label="time{}".format(i))

# Constraint not to visit the same city more than twice.
city_const = 0.0
for j in range(n_city):
    city_const += Constraint((sum(x[i, j] for i in range(n_city)) - 1)**2, label="city{}".format(j))

# distance of route
distance = 0.0
#example 1
F = np.array([
    [0, 30, 5, 10],
    [30, 0, 5, 20],
    [5, 5, 0, 50],
    [10, 20, 50, 0]
])

D = np.array([
  [0, 10, 40, 30],
  [10, 0, 20, 50],
  [40, 20, 0, 20],
  [30, 50, 20, 0]
])
n_city = 4
for i in range(n_city):
  for j in range(n_city):
          for k in range(n_city):
              for l in range(n_city):
                      distance += D[k][l] * F[i][j] * x[i][k] *x[j][l]

# Construct hamiltonian
A = Placeholder("A")
max_bound = np.max(D) *  np.max(F) + 1
print(max_bound)
H = distance + A * (time_const + city_const)

# Compile model
model = H.compile()

# Generate QUBO
feed_dict = {'A': max_bound}
bqm = model.to_bqm(feed_dict=feed_dict)
#print(bqm)

"""**1. Using Simmulated Annealing**"""

import neal
sa = neal.SimulatedAnnealingSampler()
sampleset = sa.sample(bqm,num_reads=10000)

# Decode solution
decoded_samples = model.decode_sampleset(sampleset, feed_dict=feed_dict)
best_sample = min(decoded_samples, key=lambda x: x.energy)
num_broken = len(best_sample.constraints(only_broken=True))
if num_broken == 0:
    print(best_sample)
else:
    print(best_sample.constraints())

# Optimal solution as a 2D matrix
opt = np.zeros((n_city, n_city))
for i in range(n_city):
    for j in range(n_city):
        if best_sample.array('c', (i, j)) == 1:
            opt[i][j] = 1

print(opt)
energy = sampleset.first.energy
print(energy)

"""**2. Using Hybrid Solver**"""

pip install dwave.system

from dwave.system import  LeapHybridSampler
sa = LeapHybridSampler(token = "DEV-79b52af764c9cc6cefaf61d21dd2820ef46f5f72")
sampleset = sa.sample(bqm,time_limit=3)

# Decode solution
decoded_samples = model.decode_sampleset(sampleset, feed_dict=feed_dict)
best_sample = min(decoded_samples, key=lambda x: x.energy)
num_broken = len(best_sample.constraints(only_broken=True))

if num_broken == 0:
    print(best_sample)
else:
    print(best_sample.constraints())
  
# Optimal solution as a 2D matrix
opt = np.zeros((n_city, n_city))
for i in range(n_city):
    for j in range(n_city):
        if best_sample.array('c', (i, j)) == 1:
            opt[i][j] = 1

print(opt)
energy = sampleset.first.energy
print(energy)

""" **3. Using Quantum Solver**"""

from dwave.system import DWaveSampler, EmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler(token = "DEV-79b52af764c9cc6cefaf61d21dd2820ef46f5f72"))
sampleset = sampler.sample(bqm, num_reads=1000)

# Decode solution
decoded_samples = model.decode_sampleset(sampleset, feed_dict=feed_dict)
best_sample = min(decoded_samples, key=lambda x: x.energy)
num_broken = len(best_sample.constraints(only_broken=True))

if num_broken == 0:
    print(best_sample)
else:
    print(best_sample.constraints())
  
# Optimal solution as a 2D matrix
opt = np.zeros((n_city, n_city))
for i in range(n_city):
    for j in range(n_city):
        if best_sample.array('c', (i, j)) == 1:
            opt[i][j] = 1

print(opt)
energy = sampleset.first.energy
print(energy)

"""**4. Using Classical Solver**"""

import dimod
from dimod import ExactSolver

sampler = dimod.ExactSolver()
sampleset = sampler.sample(bqm)

# Decode solution
decoded_samples = model.decode_sampleset(sampleset, feed_dict=feed_dict)
best_sample = min(decoded_samples, key=lambda x: x.energy)
num_broken = len(best_sample.constraints(only_broken=True))

if num_broken == 0:
    print(best_sample)
else:
    print(best_sample.constraints())
  
# Optimal solution as a 2D matrix
opt = np.zeros((n_city, n_city))
for i in range(n_city):
    for j in range(n_city):
        if best_sample.array('c', (i, j)) == 1:
            opt[i][j] = 1

print(opt)
energy = sampleset.first.energy
print(energy)