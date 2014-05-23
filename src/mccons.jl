



#------------------------------------------------------------------------------
#BEGIN readme



#original MC-CONS
#Input : n RNAs sequences with each m associated subotimal secondary structures
#Output: alignment of n secondary structures, giving a "consensus" structure
#(one secondary structure per RNA sequence)

#this version
#Input : same as the original + one or many RNA secondary structure distance functions
#Output: a number of nondominated sets of n secondary structures, according
#to the chosen distance functions. The output can then be filtered according to
#desired usage
#example: if the desired output is an alignment, one can run the algorithm,
#get n sets of similar structures and then run n alignments and learn
#which distance function makes more sense for the kind of alignement used.

#END   readme
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN imports


require("NSGA_II.jl")
require("RNA_2D.jl")



#END   imports
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN data




#mock testing
mock = Dict{String, String}()
mock["a"] = "AAAAAAAAAATTTTTTTTTT"
mock["b"] = "TTTTTTTTTTAAAAAAAAAA"
mock["c"] = "ATATATATATATATATATAT"
mock["d"] = "TATATATATATATATATATA"

#same length IREs
IREs30 = Dict{String,String}()
for i in keys(IREs)
  if length(IREs[i])==30
    IREs30[i] = IREs[i]
  end
end

#same length tRNA
yeastTRNAs76 = Dict{String, String}()
for i in keys(yeastTRNAs)
  if length(yeastTRNAs[i])==76
    yeastTRNAs76[i] = yeastTRNAs[i]
  end
end

#END
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN folding functions

#BEGIN callFlashFold
function callFlashFold(sequence::String, ft::Int)
  data = split(readall(`./ff -seq $sequence -ft $ft`), "\n")
  if(data[end] == "")
    data = data[1:end-1]
  end
  data = map(x->split(x), data)
  data2 = (String, FloatingPoint)[]
  for x in data
    push!(data2, ((convert(String, x[1])), float(x[2])))
  end
  return unique(data2)
end
#END callFlashFold



#BEGIN foldAll
function foldAll(dict::Dict, n::Int)
  #uses 
  #get the keys
  k = collect(keys(dict))
  vals = Vector{RNA_2D.structure}[]
  for i = 1:length(k)
    folded = callFlashFold(dict[k[i]], n)
    folds = map(x -> RNA_2D.structure(i, x[1], x[2]), folded)
    push!(vals, folds)
  end
  return vals
end
#END

#END
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN optimizing functions (precalculate or memoize)


#BEGIN precalculate
function precalculate{T}(distanceFunction::Function, alleles::Vector{Vector{T}})
  #allows memoization of the results of all possible distance function 
  #comparisons for the given alleles
  
  #O(n^2) on number of alleles in space and time (so be careful)
  
  #we require that two properties are respected
  # -distance function is symmetric
  # -fitness function is the sum of the distance function across all alleles
  
  allAlleles = reduce(vcat, alleles)
  
  comparisons = Dict{(T,T), Number}()

  len = length(allAlleles)
  
  for elem1 = 1:len
    for elem2 = elem1:len
      #add the result of the comparison to the dict
      v1 = allAlleles[elem1]
      v2 = allAlleles[elem2]
      elems = (v1, v2)
      if !(elems in keys(comparisons))
        #since it is symmetrical...
        if (v2, v1) in keys(comparisons)
          comparisons[elems] = comparisons[(v2, v1)]
        else
          distance = distanceFunction(allAlleles[elem1], allAlleles[elem2])
          comparisons[elems] = distance
        end
      end
    end
  end
  
  #create the memoized function
  function memoized{T}(genes::Vector{T})
    fitness = 0
    for i1 = 1:length(genes)
      for i2 = i1:length(genes)
        elems = (genes[i1], genes[i2])
        fitness += comparisons[elems]
      end
    end
    return fitness
  end
  
  return memoized
end
#END


#BEGIN memoize distance
function memoizeDist(distanceFunction::Function, toEvaluate::Vector{RNA_2D.structure}, memory::Dict{(String,String), Number})
  summation = 0
  len = length(toEvaluate)
  for i=1:len
    for j=i:len
      struct1 = toEvaluate[i]
      d1 = struct1.dotBracket
      struct2 = toEvaluate[j]
      d2 = struct2.dotBracket
      #dist(x,x) = 0
      if d1 != d2
        #
        if !((d1,d2) in keys(memory))
          distance = distanceFunction(struct1, struct2)
          memory[(d1,d2)] = distance
          memory[(d2,d1)] = distance
        end
        summation += memory[(d1,d2)]
      end
    end
  end
  return summation
end
#END


#END
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN main function
function main(moleculeDict = IREs, popSize = 250, numIterations = 100, alleleSize = 30)
  mutationOperator = geneticAlgorithmOperators.uniformMutate
  crossoverOperator = geneticAlgorithmOperators.uniformCrossover
  
  alleles = foldAll(moleculeDict, alleleSize)
  
  #memoize the results
  #const bpset = Dict{(String,String), Number}()
  const hausdorff = Dict{(String,String), Number}()
  const levenshtein = Dict{(String,String),Number}()
  const mountain = Dict{(String, String), Number}()
  
  #do squared error function
  squaredHausdorff(s1, s2) = (RNA_2D.compareHausdorff(s1, s2))^2
  squaredLevenshtein(s1, s2) = (RNA_2D.levenshteinDistance(s1, s2))^2
  squaredMountain(s1, s2) = (RNA_2D.compareMountainDistance(s1, s2))^2
  
  #memoizedBPSetDist(v::Vector{RNA_2D.structure}) = memoizeDist(RNA_2D.compareBPSet, v, bpset)
#   memoizedHausdorff(v::Vector{RNA_2D.structure}) = memoizeDist(RNA_2D.compareHausdorff, v, hausdorff)
#   memoizedLevenshtein(v::Vector{RNA_2D.structure}) = memoizeDist(RNA_2D.levenshteinDistance, v, levenshtein)
#   memoizedMountain(v::Vector{RNA_2D.structure}) = memoizeDist(RNA_2D.compareMountainDistance, v, mountain)
  
  memoizedHausdorff(v::Vector{RNA_2D.structure}) = memoizeDist(squaredHausdorff, v, hausdorff)
  memoizedLevenshtein(v::Vector{RNA_2D.structure}) = memoizeDist(squaredLevenshtein, v, levenshtein)
  memoizedMountain(v::Vector{RNA_2D.structure}) = memoizeDist(squaredMountain, v, mountain)
  
  #create the evaluation vector
  evalDistance(v::Vector) = [-memoizedLevenshtein(v), -memoizedHausdorff(v), -memoizedMountain(v)]
  println("memoized functions done")
  
  r = NSGA_II.main(alleles,
                   evalDistance,
                   popSize,
                   numIterations,
                   0.2,
                   0.2,
                   crossoverOperator,
                   mutationOperator)
  println("comparisons calculated : $(length(hausdorff)/2)")
  return (r, alleles)
end
#END
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN io

function writeResult{T<:String}(fileName::T, P::NSGA_II.population, alleles::Vector)
  #will not rewrite files in theory
  #assert that the file is not already in the local directory
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
  
  #format is:
  #>id
  #fitness
  #dotBrackets
  #end
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



function extractFitness{T<:String}(fileName::T)
  fileName = string(fileName)
  if !(fileName in readdir())
    error("$fileName is not in the local directory")
  end
  
  f = open(fileName, "r")
  lines = readlines(f)
  #println(lines)
  close(f)
  
  f2 = open(string("fitness_",fileName), "w")
  #println(lines[1])
  i=1
  while i < length(lines)
    #println(i)
    #println(lines[i])
    if beginswith(lines[i], ">")
      index = chomp(lines[i][2:end])
      fitnesses = chomp(lines[i+1])
      println(f2, string(index,",", fitnesses))
      i+=2
    else
      i+=1
    end
  end
  close(f2)

end



function writeAndExtract{T}(fileName::T, P::NSGA_II.population, alleles::Vector)
  writeResult(fileName, P, alleles)
  extractFitness(fileName)

end

#END
#------------------------------------------------------------------------------



