# metaUKPstu: solving the 0/1 unidimensional knapsack problem
Implementation in Julia (compliant Julia v1.2) of well-known optimization algorithms applied on the 0/1 unidimensional knapsack problem (01UKP) for pedagogical purposes.

------

Elements de soutien en Julia (v1.2 et suivants) des exercices des cours "optimisation discrète et combinatoire" et "métaheuristiques" en master 1 informatique parcours "Optimisation en Recherche Opérationnelle (ORO)", année 2019-2020.

------

Algorithms available: 
  
-  Generator of instances for the 01UKP
-  Upper bound given by linear relaxation of the 01UKP
-  Lower bound given by a random feasible solution for the 01UKP
-  Lower bound given by a solution computed by descend constructive method for the 01UKP
-  Simple exploration heuristic (heuExplore)
-  Simulated Annealing metaheuristic (metaSA)
-  swap move for the 01UKP (swap)
-  add_or_drop move for the 01UKP (addOrDrop)
-  Plotting solutions and the cooling schedule
