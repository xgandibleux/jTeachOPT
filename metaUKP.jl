# ============================================================
# solveUKP 
#
#  Bounds, heuristics and metaheuristics applied on the 01UKP implemented in Julia (compliant v1.0.0)
#
#  Xavier Gandibleux (Xavier.Gandibleux@univ-nantes.fr)
#  Universite de Nantes, Faculty of sciences and technologies
#
#  Ver 0.2.0 - Released : 27 March 2016 - Modified : 24 August 2018
#  Objectif de ce package (English version will follow soon):
#
#  Ensemble de primitives julia a vocation pedagogique en support a l'enseignement
#  des metaheuristisques et du sac-a-dos mono-objectif unidimensionnel en variable 01 (01UKP).
#
#  Developpe pour aider a (1) l'illustration des concepts en optimisation et -plus modestement-
#  a la mise en oeuvre de julia, (2) base au travail de developpement de composants additionnels
#  par les etudiants.
#
#  Utilise depuis 2016-2017 en M1 informatique parcours optimisation en recherche operationnelle :
#  course Metaheuristics
#
# Avertissement :
#
#  Etant aussi aficionado du C et de ADA, on retrouvera un style c-like dans l'implementation.

# Algorithms/procedures implemented :
#
#  Linear relaxation of the 01UKP
#  Random feasible solution for the 01UKP
#  Descend constructive method for the 01UKP
#  Swap move for the 01UKP (swap)
#  Add_or_drop move for the 01UKP (addOrDrop)
#  Simple exploration heuristic (heuExplore)
#  Simulated Annealing metaheuristic (metaSA)
#  Computing and plotting of a cooling schedule (coolingIllustration)

# ============================================================

# Compulsory to be compliant with Julia v1.0.0:

using Random
using Printf
using LinearAlgebra

# ------------------------------------------------------------
# Global constants for reporting or not the activity of algorithms

const verbose = false

# ------------------------------------------------------------
# Global vectors storing the solutions for graphical purposes

global zBest    = [] # each best values of f(x) collected 
global zBestAll = [] # best values of f(x) collected for each iteration
global zAll     = [] # all values of f(x) collected for each iteration
global tAll     = [] # values of the cooling schedule for each iteration

# ------------------------------------------------------------
# Datastructure of a 01UKP instance (single objective unidimensional binary knapsack problem)

mutable struct instance

    n ::Int64      # instance size (integer)
    c              # cost of items (integer)
    w              # weight of items (integer)
    W ::Int64      # rhs (integer)

end

# ------------------------------------------------------------
# Datastructure of a 01UKP solution 

mutable struct solution

    # for all problems defined by one vector of binary variables  1..n

    x              # variables (binary)
    v0             # index of els of the variable equal to 0 (integer)
    v1             # index of els of the variable equal to 1 (integer)
    z    ::Int64   # performance (integer)

    # for the UKP

    r    ::Int64   # residual capacity of the constraint (integer)
    somX ::Int64   # sum of x_i = 1 (integer)

end

# ------------------------------------------------------------
# split of a solution into 2 index vectors V0 and V1 according respectively if x_i=0 or x_i=1

function splitX(s::solution)
    s.v0 = [] ; s.v1 = []
    for i=1:length(ukp.c)
        if s.x[i] == 0
            push!(s.v0,i)
        else
            push!(s.v1,i)
        end
    end
end

# ------------------------------------------------------------
# evaluate the objective function for a given solution

function fctUKP(c,x)
    return dot(c,x)   # sum_{i=1}^{n} c_i x_i
end

# ------------------------------------------------------------
# Build randomly a feasible neighboor with a swap (random exchange)
# Discussed in Chapter 2, course Metaheuristics

function swap(s::solution)
    i0::Int64=-1;  i1::Int64=-1   #astuces
    i01::Int64=-1; i10::Int64=-1  #astuces

    sVoisin = deepcopy(s)

    # repeter jusqu'a obtenir un swap produisant une solution realisable
    # (ne traite pas les cas pathologiques ie swap pas possible et donc cyclage)
    while true
        i0  = rand(1:length(sVoisin.v0)) ; i1  = rand(1:length(sVoisin.v1))
        i01 = sVoisin.v0[i0] ;             i10 = sVoisin.v1[i1]

        if (verbose == true)
            @printf "swap) i01: %d - i10: %d \n" i01 i10
        end

        ((dot(ukp.w,sVoisin.x) - ukp.w[i10] + ukp.w[i01]) <= ukp.W)  && break # repeat...until

    end

    # l'echange induit une solution realisable
    sVoisin.x[i01]=1;  sVoisin.x[i10]=0
    sVoisin.v0[i0]=i10; sVoisin.v1[i1]=i01;
    sVoisin.z = fctUKP(ukp.c, sVoisin.x)
    sVoisin.r = ukp.W - dot(ukp.w,sVoisin.x)
    if (verbose == true)
        @printf "swap) z=%d x=%s r=%d \n" dot(ukp.c,sVoisin.x) sVoisin.x ukp.W-dot(ukp.w,sVoisin.x)
    end
    return sVoisin
end

# ------------------------------------------------------------
# Build randomly a feasible neighboor with an add_ou_drop (flip 0to1 or 1to0 randomly)
# Discussed in Chapter 2, course Metaheuristics

function addOrDrop(s::solution)
    i0::Int64=-1;  i1::Int64=-1   #astuces
    i01::Int64=-1; i10::Int64=-1

    sVoisin = deepcopy(s)

    i0  = rand(1:length(sVoisin.v0)) ; i1  = rand(1:length(sVoisin.v1))
    i01 = sVoisin.v0[i0] ;             i10 = sVoisin.v1[i1]

    if (verbose == true)
        @printf "swap) i01: %d - i10: %d \n" i01 i10
    end

    if (dot(ukp.w,sVoisin.x) + ukp.w[i01]) <= ukp.W
        # add possible => applique
        sVoisin.x[i01]=1;
    else
        if length(sVoisin.v1) > 1 # il existe au moins une variable x_i=1
            # drop possible => applique
            sVoisin.x[i10]=0
        end
    end

    # faineantise, peut mieux faire sans reconstruire le split
    splitX(sVoisin)

    sVoisin.z = fctUKP(ukp.c, sVoisin.x)
    sVoisin.r = ukp.W - dot(ukp.w,sVoisin.x)
    if (verbose == true)
        @printf "addDrop) z=%d x=%s r=%d \n" sVoisin.z sVoisin.x sVoisin.r
    end
    return sVoisin
end

# ------------------------------------------------------------
# Generate randomly an instance for the UKP with c and w uniformly distributed

function generateRandomlyInstanceUKP(n = 100, max_ci = 100, max_wi = 30)

    verboseUtility = false # rapporte (ou pas) les items par ordre decroissant

    # --- creation de l'instance
    rnd_c = rand(1:max_ci,n); # c_i \in [1,max_ci]
    rnd_w = rand(1:max_wi,n) # w_i \in [1,max_wi]
    
    # rank the items according the decreasing values u_i = c_i/w_i
    utilite = rnd_c ./ rnd_w
    reord = sortperm(utilite, rev=true)
    ukp   = instance(n, zeros(n), zeros(n), 0)
    for i = 1:n
        ukp.c[i] = rnd_c[reord[i]]
        ukp.w[i] = rnd_w[reord[i]]
        if (verboseUtility == true)
            @printf "(%d %d %.2f) \n " ukp.c[i] ukp.w[i] utilite[reord[i]]
        end
    end
            
    ukp.W = round(Int64, sum(ukp.w)/2)
                
    return ukp
end

# ------------------------------------------------------------
# Descend method for computing a greedy solution for the UKP
# Discussed in Chapter 2, course Metaheuristics

function computeGreedySolutionUKP(ukp)
    
    # ---
    # Calcule la solution gloutonne avec pour utilite : u(i) = c(i)/p(i)
    sGreedy = solution(zeros(Int64, ukp.n), [], [], 0, 0, 0)
    sommeW = 0
    for i = 1:ukp.n
        #@printf "%f\n" ukp.c[i]/ukp.w[i]
        if (sommeW + ukp.w[i] <= ukp.W)
            sGreedy.x[i] = 1
            sommeW = sommeW + ukp.w[i]
        end
    end
    sGreedy.z = fctUKP(ukp.c, sGreedy.x)
    sGreedy.r = ukp.W - dot(ukp.w, sGreedy.x)
    sGreedy.somX=sum(sGreedy.x) # somme des x_i = 1

    # ---
    # scinde le vecteur de x binaire en 2 vecteurs d'indices
    splitX(sGreedy)
                        
    return sGreedy
end

# ------------------------------------------------------------
# Compute the linear relaxation for the 01UKP
# Discussed in Chapter x, course Integer Programming

function computeLinearRelaxationUKP(ukp)
    
    # identify the last item integrally selected
    sommeW = 0 ; s = 1
    while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
        sommeW = sommeW + ukp.w[s]
        s = s + 1
    end
    s = s - 1
                
    # compute the upper bound
    if ((ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) > 0) && (s<ukp.n)
        # contrainte non saturee => ajout de la partie fractionnaire de l'item bloquant
        zRelax = sum(ukp.c[1:s]) + (ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) * (ukp.c[s+1] ./ ukp.w[s+1])
    else
        # contrainte saturee => rien a faire
        zRelax = sum(ukp.c[1:s])
    end
    return zRelax
end

# ------------------------------------------------------------
# Compute randomly a feasible solution for the 01UKP
# Discussed in Chapter 2, course Metaheuristics

function computeRandomSolutionUKP(ukp)
    
    # ---
    # Construit aleatoirement une solution initiale realisable
    s0 = solution(zeros(Int64, ukp.n), [], [], 0, 0, 0)
    semence = randperm(ukp.n) # une permutation aleatoire de (1..n)
    sommeW = 0
    for i = 1:ukp.n
        if (sommeW + ukp.w[semence[i]] <= ukp.W)
            s0.x[semence[i]] = 1
            sommeW = sommeW + ukp.w[semence[i]]
        end
    end
    s0.z = fctUKP(ukp.c,s0.x)
    s0.r = ukp.W - dot(ukp.w,s0.x)

    # ---
    # scinde le vecteur de x binaire en 2 vecteurs d'indices
    splitX(s0)

    return s0
end

# ------------------------------------------------------------
# Elementary exploration of a neighborhood of a solution with a random 
# Discussed in Chapter 2, course Metaheuristics
# Call :
#          heuExplore(s0, zBest)
#          @printf "\n\nzBest=%s \n" zBest

function heuExplore(s::solution, zBest)
    push!(zBest,s.z)

    essais = 10 # nombre de voisins construits
    for essai = 1:essais
        if (verbose == true)
            @printf "\n \nv0=%s v1=%s\n" s.v0 s.v1
        end
        sVoisin = deepcopy(swap(s))
        if (verbose == true)
            @printf "\nv0=%s v1=%s\n" sVoisin.v0 sVoisin.v1
            @printf "z=%d  zVoisin=%d " s.z sVoisin.z
        end
        if sVoisin.z > s.z
            push!(zBest,sVoisin.z)
            s = deepcopy(sVoisin)
        end
    end
end

# ------------------------------------------------------------
# Simulated Annealing metaheuristic for a function to maximize
# Discussed in Chapter 4, course Metaheuristics

function metaSA(s0::solution, t0::Float64, lPalier::Int64, α::Float64, move, tLow::Float64, nogoodMax::Int64, verbose)

    verboseDecision = false # rapporte (ou pas) les 4 decisions de l'algo (A++,A+,A-,R)

    if (verbose == true)
        @printf "t0= %d, lPalier= %d, α= %f \n" t0 lPalier α
    end

    sCur  = deepcopy(s0) # solution courante
    sBest = deepcopy(s0) # meilleure solution

    push!(zAll,s0.z) #faineantise : des structures globales pour les graphiques
    push!(zBest,s0.z)
    push!(zBestAll,s0.z) ; zMax = s0.z
    push!(tAll,t0)

    t = t0
    nogood = 0

    iter=1
    while true

        for palier = 1:lPalier
            if (verbose == true)
                @printf "\n \nv0=%s v1=%s\n" sCur.v0 sCur.v1
            end
            sVois = deepcopy(move(sCur))
            push!(zAll,sCur.z)
            if (verbose == true)
                @printf "\nv0=%s v1=%s\n" sVois.v0 sVois.v1
                @printf "zCur=%d  zVois=%d || " sCur.z sVois.z
            end

            Δz = sVois.z - sCur.z
            prob = rand()
            if (verboseDecision == true)
                @printf "%4d) Δz=%5d  t=%6.2f Pa=%5.3f p=%5.3f || " iter Δz t exp(-(sCur.z - sVois.z)/t) prob
            end

            if ((Δz >= 0) || (exp(-(sCur.z - sVois.z)/t)>prob))

                sCur = deepcopy(sVois)
                nogood = 0

                if (sVois.z > sBest.z)
                    sBest = deepcopy(sVois)
                    @printf "A++ => zBest=%d\n" sBest.z
                    push!(zBest,sBest.z)
                    zMax = sBest.z
                else
                    if (Δz >= 0)
                        if (verboseDecision == true)
                            @printf "A+ \n"
                        end
                    else
                        if (verboseDecision == true)
                            @printf "A- \n"
                        end
                    end
                end
            else
                nogood +=1
                if (verboseDecision == true)
                    @printf "R  \n"
                end
            end
            push!(zBestAll,zMax)
            push!(tAll,t)
            iter = iter +1
        end

        t = t * α

        ((t<tLow) || (nogood >= nogoodMax)) && break # repeat...until
    end
    if (t<tLow) @printf("stopping SA because : t<tLow\n") end
    if (nogood >= nogoodMax) @printf("stopping SA because : nogood >= nogoodMax\n") end
    return sBest
end

# ------------------------------------------------------------
# Plot a graph with the results obtained by the SA (z for all solutions, best, greedy)

function plotSolutionsSA()
    title("01UKP | n=" * string(ukp.n)) # * " | α=" * string(α) * " | lPalier=" * lPalier )
    xlabel("Algorithm's Iterations")
    ylabel("Function Value")
    grid()
    #xlim(0, length(zBest))
    #ylim(0, length(zAll))
    #ylim(0, int64(floor(maximum(zBest)/10)+1)*10)
    #plot(zBest, marker = "o")

    plot(zAll,label="all solutions")
    plot(zBestAll, linewidth=2.0, color="green", label="best solutions")

    vGreedy=fill(sGreedy.z,length(zAll))
    vRelax=fill(zRelax,length(zAll))
    plot(vGreedy, linewidth=1.0, color="green",linestyle="--", label="greedy solution")
    plot(vRelax, linewidth=1.0, color="red",linestyle="-", label="UB (linear relaxation)")
    legend(loc=4, fontsize ="small")
end

# ------------------------------------------------------------
# Plot a graph of the cooling scheduled computed

function plotTemperatureSA()
    title("Cooling Schedule")
    xlabel("iterations")
    ylabel("temperature")
    plot(tAll)
end

# ------------------------------------------------------------
# Compute and plot a cooling schedule for illustration purposes
# call :
#         coolingIllustration(t0=100 , α=0.7,  lgPalier=6 , tLow = 1)

function coolingIllustration(t0=100.0 , α=0.7,  lgPalier=6 , tLow = 1.0)
    # t0=100; alpha=0.9; lgPalier = 5 ; tLow = 1

    @printf "t0 = %d, lgPalier = %d, α=%f \n" t0 lgPalier α

    t=t0; i=1;
    ay= [] #; push!(ay,t)

    while true
        for j in 1:lgPalier
            push!(ay,t)
            i=i+1
        end
        t = t * α
        t < tLow && break # repeat...until
    end

    ax = [ i for i=1:length(ay) ] # [1:1:length(ay)]

    title("Cooling Schedule")
    xlabel("iterations")
    ylabel("temperature")
    b = bar(ax,ay,color="#0f87bf",align="center", width=1.0, edgecolor ="lightgrey")

    # sortie ecran du tableau de valeurs
    @printf "   i   k      t \n"
    for i = 1:length(ay)
        @printf "%4d %3d %6.2f \n" i ceil(Int64,i/lgPalier) ay[i]
    end
end


# ============================================================
# MAIN ENTRY POINT
# ============================================================

# ------------------------------------------------------------
# Generate an instance for the 01UKP (n, max_ci, max_wi)

ukp = generateRandomlyInstanceUKP(100, 100, 30)

# ------------------------------------------------------------
# Compute a lower and upper bound on the optimal solution

sGreedy = computeGreedySolutionUKP(ukp)

zRelax = computeLinearRelaxationUKP(ukp)

# ------------------------------------------------------------
# Compute a feasible solution

s0 = computeRandomSolutionUKP(ukp)

# ------------------------------------------------------------

if (verbose == true)
    @printf "init) c=%s\n" ukp.c
    @printf "init) w=%s W=%d\n" ukp.w ukp.W
    @printf "init) z=%d x=%s r=%d\n" s0.z s0.x s0.r
end

# ------------------------------------------------------------
# Numerical experiment => perform nbrRuns

nbrRuns = 1
sBest = solution(zeros(Int64, ukp.n), [], [], 0, 0, 0)
allzBest = []

for run = 1:nbrRuns

    # ----------------------------------------------------------
    # Search a good solution with the simulated annealing

    # move 1: improve a feasible solution by modifying the cardinality of variables equal to 1
    t0 = 300.0 ; lPalier = ceil(Int64, 1.5 * ukp.n) ; α = 0.95 ; tLow =0.5 ; nogoodMax = 5 * ukp.n
    @printf "t0 = %d, lPalier = %d, α=%f \n" t0 lPalier α
    sBest=metaSA(s0, t0, lPalier, α, addOrDrop, tLow, nogoodMax, verbose)

    # move 2: improve a feasible solution in maintaining the cardinality of variables equal to 1
    t0 = 50.0 ; lPalier = ceil(Int64, 2.5 * ukp.n) ; α = 0.8 ; tLow =0.05 ; nogoodMax = 10 * ukp.n
    @printf "t0 = %d, lPalier = %d, α=%f \n" t0 lPalier α
    sBest=metaSA(sBest, t0, lPalier, α, swap, tLow, nogoodMax, verbose)

    push!(allzBest,sBest.z)

end #nbrRuns

# ------------------------------------------------------------
# Plotting the results

#Pkg.add("PyPlot") # Mandatory before the first use of this package
using PyPlot
plotSolutionsSA()
#plotTemperatureSA()


# Quantitative summary of results collected
@printf "\n zS0=%d (%3d)  zMaxSA=%d (%3d)   zGreedy=%d (%3d)  zRelax=%7.2f" s0.z sum(s0.x)  sBest.z sum(sBest.x)  sGreedy.z sum(sGreedy.x)  zRelax

# ============================================================

