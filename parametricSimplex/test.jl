include("simplex_XG.jl")
using SimplexMethod

include("simplexparam_XG.jl")
#using ParametricSimplexMethod

# Example 1 

c1 = [ 3;  1; 0; 0]
c2 = [-1; -2; 0; 0]
T  = [ 0  1  1  0 ;
       3 -1  0  1 ]
d = [3; 6]

c1 = Array{Float64}(c1)
c2 = Array{Float64}(c2)
T  = Array{Float64}(T)
d  = Array{Float64}(d)

println("Objective 1: ")
@time opt_x1, obj1 = simplex_method(c1, T, d)
println(" ")
println("Objective 2: ")
@time opt_x2, obj2 = simplex_method(c2, T, d)
println(" ")
println("Objective 1&2: ")
@time  parametric_simplex_method(c1, c2, T, d)
println(" ")
