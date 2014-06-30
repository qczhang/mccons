# EXAMPLE 1: the IREs

#in this particular case, there is enough
#similarity in the structures to get good
#sets of similar structures easily


#------------------------------------------------------------------------------
#load the fasta reader and get the secondary structures
require("../utility/fastaReader.jl")
data = fastaRead("data.txt")


#convert the strings to RNA 2D structure representation
require("../../src/RNA_2D.jl")
alleles = Vector{RNA_2D.structure}[]

#since most of the IREs don't have the same size, and we want to
#use distance functions that depend on it, we'll use only the IREs of size 30...

for (name, structures) in data
    if length(structures[1]) == 76
        RNA2Dstructures = RNA_2D.structure[]
        for struct in structures
            push!(RNA2Dstructures, RNA_2D.structure(name, struct))
        end
        push!(alleles, RNA2Dstructures)
    end
end

#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#choose the distance functions, we will use naive ones already
#implemented in the RNA_2D module.
#We modify them to compare all against all structures in the vector of
#structures that will constitute a solution and square them 
#to penalize in a nonlinear fashion.

#we square the distance functions to get a nonlinear error

#the hausdorff distance
squaredHausdorff(s1, s2) = (RNA_2D.compareHausdorff(s1, s2))^2
#the levenshtein distance
squaredLevenshtein(s1, s2) = (RNA_2D.levenshteinDistance(s1, s2))^2
#the mountain distsance
squaredMountain(s1, s2) = (RNA_2D.compareMountainDistance(s1, s2))^2


#we compare all structure against all structure
#(but assume that the distance is symmetric, so we only calculate
#the lower triangle of the matrix)
function compareLowerTriangle(v::Vector, evalFunction::Function)
    distance = 0
    for i = 1:length(v)-1
        for j = i+1:length(v)
            distance += evalFunction(v[i], v[j])
        end
    end
    distance
end

ltSquaredHausdorff(v::Vector) = compareLowerTriangle(v::Vector, squaredHausdorff)
ltSquaredLevenshtein(v::Vector) = compareLowerTriangle(v::Vector, squaredLevenshtein)
ltSquaredMountain(v::Vector) = compareLowerTriangle(v::Vector, squaredMountain)


evalDistance(v::Vector) = [-ltSquaredHausdorff(v),
                           -ltSquaredLevenshtein(v),
                           -ltSquaredMountain(v)]
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
# a small function to write the result to a file
function writeResult{T<:String}(fileName::T, P, alleles::Vector)
  #use it to write to specified file name
  #give alleles to know what rank the structure was in the suboptimals
  listDir = readdir()
  name = string(fileName)
  if name in listDir
    error("$name is already present in the local directory, delete it if you want to write")
  end
  #extract the dotbrackets fron the structures in alleles
  alleles2 = Vector{String}[]
  for i in alleles
    toAdd = String[]
    for j in i
      push!(toAdd, j.dotBracket)
    end
    push!(alleles2, toAdd)
  end
  f = open(name, "w")

  ind = P.individuals
  fitnessLen = length(ind[1].fitness)

  for i = 1:length(ind)
    #write the id
    println(f, ">$i")
    fit = ind[i].fitness
    s = fit[1]
    for j = 2:fitnessLen
      s = string(s,",", fit[j])
    end
    #write the fitnesses on the same line
    println(f, s)

    #write the dot brackets associated
    for j=1:length(ind[i].genes)
      struct = ind[i].genes[j].dotBracket
      index = findfirst(alleles2[j], struct)
      println(f, string(struct, " ", index))
    end
  end
  close(f)
end



#we load the NSGA_II algorithm and let it do its magic (we could choose 
#more refined parameters 
require("../../src/NSGA_II.jl")
function main(popSize::Int, numGenerations::Int)
    result = NSGA_II.main(alleles, #the structure alternatives
                          evalDistance, # the evaluation function
                          popSize, #the size of the population, choosen on whim in this case
                          numGenerations, #the number of generations for the genetic algorithm
                          0.2, #the probabiltiy of crossover for the genetic algorithm
                          0.2) #the probabilty of mutation for the genetic algorithm
end



#to launch, write these for example
result = main(500, 50)
writeResult("result.txt", result[1], alleles)