# ============================================================
# jTeachOPT
#
#  Heuristics, metaheuristics and exact algorithms applied on the 01UKP implemented in Julia
#
#  Xavier Gandibleux (Xavier.Gandibleux@univ-nantes.fr)
#  Universite de Nantes, Faculty of sciences and technologies - IRCCyN UMR CNRS 6597
#
#  Ver 0.1.3 - Released : 23 Mars 2016

# Objectif de ce package (English version will follow soon):
#
#  Ensemble de primitives julia a vocation pedagogique en support a l'enseignement
#  des metaheuristisques et du sac-a-dos mono-objectif unidimensionnel en variable 01 (01UKP).
#  Developpe pour aider a (1) l'illustration des concepts en optimisation et -plus modestement-
#  a la mise en oeuvre de julia, (2) base au travail de developpement de composants additionnels
#  par les etudiants.

# Utilise (des 2016-2017) en M1 informatique parcours optimisation en recherche operationnelle :
#
#  cours Metaheuristics

# Avertissement :
#
#  S'agissant d'une premiere experience personnelle avec Julia, et etant aficionado du C,
#  on verra un style c-like dans l'implementation et plusieurs astuces sont appliquees
#  pour repondre a des traitements que je n'ai pas vu comment faire davantage dans l'esprit de julia.

# Algorithms/procedures implemented :
#
#  Descend constructive method for the 01UKP
#  Random feasible solution for the 01UKP
#  Linear relaxation of the 01UKP
#  Simple exploration heuristic (heuExplore)
#  Simulated Annealing metaheuristic (metaSA)
#  swap move for the 01UKP (swap)
#  add_or_drop move for the 01UKP (addOrDrop)
#  Computing and plotting of a cooling schedule (examenRefroidissement)

# ============================================================

# ------------------------------------------------------------
# constantes globales (besoin au bavardage de l'algo)

verbose = false

# ------------------------------------------------------------
# structure d'une instance de 01UKP (sac a dos unidimmensionel en variables 01)

type instance
  n  # taille de l'instance
  c  # couts des items
  w  # poids des items
  W  # rhs
end

# ------------------------------------------------------------
# structure d'une solution mono-objectif

type solution
  x  # variable
  v0 # indices des els de la variable a 0
  v1 # indices des els de la variable a 1
  z  # performance
  r  # capacite residuelle # propre au 01UKP
  somX # somme des x_i = 1
end

# ------------------------------------------------------------
# evalue la fonction objectif en une solution realisable

function fctUKP(c,x)
  return dot(c,x) # sum_{i=1}^{n} c_i x_i
end

# ------------------------------------------------------------
# split de la solution en 2 vecteurs d'indices de variables a 0 et a 1

function splitX(s)
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
# construit un voisin aleatoire realisable par swap (echange aleatoire)

function swap(s)
  i0::Int64=-1;  i1::Int64=-1   #astuces
  i01::Int64=-1; i10::Int64=-1

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
# construit un voisin aleatoire realisable par add_ou_drop (flip 01 ou 10 aleatoire)

function addOrDrop(s)
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
# procede a l'exploration elementaire du voisinage d'une solution a partir d'un swap aleatoire
# Appel :
#          heuExplore(s0, zBest)
#          @printf "\n\nzBest=%s \n" zBest

function heuExplore(s, zBest)
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
# metaheuristique "Recuit simule" pour une fonction a maximiser

function metaSA(s0, t0, lPalier, α, move, tLow, verbose)

  verboseDecision = false # rapporte (ou pas) les 4 decisions de l'algo (A++,A+,A-,R)

  if (verbose == true)
    @printf "t0= %d, lPalier= %d, α= %f \n" t0 lPalier α
  end

  sCur  = deepcopy(s0) # solution courante
  sBest = deepcopy(s0) # meilleure solution

  push!(zAll,s0.z) #faineantise : des structures globales pour les graphiques
  push!(zBest,s0.z)
  push!(zMaxAll,s0.z) ; zMax = s0.z
  push!(vtemp,t0)

  t = t0

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
        if (verboseDecision == true)
          @printf "R  \n"
        end
      end
      push!(zMaxAll,zMax)
      push!(vtemp,t)
      iter = iter +1
    end

    t = t * α

    (t<tLow) && break # repeat...until
  end
  return sBest
end

# ------------------------------------------------------------
# Trace un graphique des resultats du SA (z pour toutes solutions, meilleures, glouton)

function traceSolutionsSA()
title("01UKP | n=" * dec(ukp.n)) # * " | α=" * dec(α) * " | lPalier=" * lPalier )
  xlabel("itérations de l'algorithme")
  ylabel("valeurs de f(x)")
  grid()
  #xlim(0, length(zBest))
  #ylim(0, length(zAll))
  #ylim(0, int64(floor(maximum(zBest)/10)+1)*10)
  #plot(zBest, marker = "o")

  plot(zAll,label="all solutions")
  plot(zMaxAll, linewidth=2.0, color="green", label="best solutions")

  vGreedy=fill(sGreedy.z,length(zAll))
  vRelax=fill(zRelax,length(zAll))
  plot(vGreedy, linewidth=1.0, color="green",linestyle="--", label="greedy solution")
  plot(vRelax, linewidth=1.0, color="red",linestyle="-", label="UB (linear relaxation)")
  legend(loc=4, fontsize ="small")
end

# ------------------------------------------------------------
# Trace la courbe de refroidissement du SA

function traceTemperatureSA()
  title("Schéma de refroidissement")
  xlabel("itérations")
  ylabel("température")
  plot(vtemp)
end

# ------------------------------------------------------------
# calcule et trace un schema de refroidissement a finalite d'etude
# appel :
#         examenRefroidissement(t0=100 , α=0.7,  lgPalier=6 , tLow = 1)

function examenRefroidissement(t0=100 , α=0.7,  lgPalier=6 , tLow = 1)
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

  title("Schéma de refroidissement")
  xlabel("itérations")
  ylabel("température")
  b = bar(ax,ay,color="#0f87bf",align="center", width=1.0, edgecolor ="lightgrey")

  # sortie ecran du tableau de valeurs
  @printf "   i   k      t \n"
  for i = 1:length(ay)
    @printf "%4d %3d %6.2f \n" i ceil(Int8,i/lgPalier) ay[i]
  end
end

# ------------------------------------------------------------
# Generate randomly an instance for the UKP with c and w unformly distributed
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
                        
    return sGreedy
end

# ------------------------------------------------------------
# Compute the linear relaxation for the 01UKP
function computeLinearRelaxationUKP(ukp)
    
    # identify the last item integrally selected
    sommeW = 0 ; s = 1
    while (s <= ukp.n) && (sommeW + ukp.w[s] <= ukp.W)
      sommeW = sommeW + ukp.w[s]
      s = s + 1
    end
    s = s - 1
                
    # compute the upper bound
    if (ukp.W - dot(ukp.w[1:s], sGreedy.x[1:s])) > 0
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


# ============================================================
# POINT d'ENTREE PRINCIPAL
# ============================================================

# ------------------------------------------------------------
# Global vectors storing the solutions for graphical purposes

zBest   = [] # meilleures obtenues
zAll    = [] # toutes
zMaxAll = [] # max de toutes
vtemp   = [] # courbe refroidissement

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
# Search a good solution with the simulated annealing

# move 1: pour flirter avec le glouton a sommeX variable
t0 = 300 ; lPalier = ceil(Int64, 1.5 * ukp.n) ; α = 0.95 ; tLow =0.5
@printf "t0 = %d, lPalier = %d, α=%f \n" t0 lPalier α
sBest=metaSA(s0, t0, lPalier, α, addOrDrop, tLow, verbose)

# move 2: pour creuser autour du glouton a sommeX fixe
t0 = 50 ; lPalier = ceil(Int64, 2.5 * ukp.n) ; α = 0.8 ; tLow =0.05
@printf "t0 = %d, lPalier = %d, α=%f \n" t0 lPalier α
sBest=metaSA(sBest, t0, lPalier, α, swap, tLow, verbose)


# ------------------------------------------------------------
# Plotting the results

#Pkg.add("PyPlot") # Mandatory before the first use of this package
using PyPlot
traceSolutionsSA()
#traceTemperatureSA()


# Quantitative summary of results collected
@printf "\n zS0=%d (%3d)  zMaxSA=%d (%3d)   zGreedy=%d (%3d)  zRelax=%7.2f" s0.z sum(s0.x)  pop!(zMaxAll) sum(sBest.x)  sGreedy.z sum(sGreedy.x)  zRelax

# ============================================================

