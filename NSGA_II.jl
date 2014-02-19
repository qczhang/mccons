#module NSGA_II
#-------------------------DEFINITION-------------------------
#Simple implementation of the NSGA-II multiobjective
#genetic algorithm.

#-------------------------imports-------------------------



#-------------------------exported methods-------------------------
#export



#-------------------------type definitions-------------------------
immutable solution
  units::Vector
  fitness::Vector
end



type population
  individuals::Vector{solution}
  distances::Dict{Vector, FloatingPoint} #this stores the distances result
end



#-------------------------genetic algorithm methods-------------------------
function nonDominatedSort(pop::population) 
  #multi-objective optimization using evolutionnary algorithms p.43
  #input: 2N population
  #output: at least N individuals, the output likely has one front too much
  
  #get cutoff
  len = length(pop.individuals)
  cutoff = div(len, 2) + len%2
  
  
  #1- evaluate the whole population
  values = (Int, Int, Array{Int,1})[]
  for i = 1:len
    push!(values, evaluateAgainstOthers(pop, i, nonDominatedCompare))
  end
  
  #2- nondominated sort
  fronts = Array{Int,1}[]
  while length(values) > cutoff
    #find the "dominators"
    front = filter(x->x[2] == 0, values)
    frontIndices = map(x->x[1], front)
    push!(fronts, frontIndices)
    
    #exclude the dominators from the values (could reuse)
    values = filter(x->x[2] != 0, values)
    
    #decrement the count based on the latest front
    for i in 1:length(values)
      #delete the last front from the indices
      substracted = fastDelete(values[i][3], frontIndices)
      #substract the difference of cardinality
      cardinality = length(substracted)
      values[i] = (values[i][1], cardinality, substracted)
    end
  end
  return fronts
end



function crowdingDistance(pop::Vector{solution})
  #calculate the crowding distance 
  #step 1
  #assing index to keep track after sorting
  values = (Int, Vector)[]
  for i = 1:length(pop)
    push!(values, (i, pop.individuals[i].fitness))
  end
  
  #step 2 
  #get the unique fitness vectors
  fitnessMapping = uniqueFitness(values)
  fitnesses = keys(fitnessMapping)
  
  #step 3 
  #do the crowding distance
  popSize = length(fitnesses) #we need the number of uniques
  vectorSize = length(fitnesses[1]) #all same size
  
  #we use only the fitnesses
  sorts = {} #think about the type...
  for i = 1:vectorSize
    push!(sorts, sort(fitnesses, by = x->x[i]), rev =true)
  end
  
  #get the max and min of each objective
  minFitnesses = map(x->x[end], sorts)
  maxFitnesses = map(x->x[1],   sorts)
  rangeSize = map(x->x[1]-x[2], zip(maxFitnesses, minFitnesses))
  
  #test
  map(x-> @assert x>=0, rangeSize )
  
  #assign infinite distance to extremities
  distances = Dict{Vector{Int}, Float}()
  #initialize at zero
  for i in fitnesses
    distances[i] = 0
  end
  
  map(x->distances[x[end]] = Inf, sorts)
  map(x->distances[x[1]]   = Inf, sorts)
  
  #calculate the other ones
  for i = 1:vectorSize
    for j = 2:(popSize-1) #first and last are already calculated
      distances[sorts[i][j]] += (sorts[i][j-1] - sorts[i][j+1])/rangeSize[i]
    end
  end
  
  return distances
end



function lastFrontSelection!(pop::population, fronts::Array{Array{Int,1},1}, k::Int)
   #input: 2N population, non dominated fronts, individuals to select from last set
   #calculate crowding distance for all but last front
   if length(fronts)>1
    for i = 1:length(fronts)-1
      p = pop.individuals[fronts[i]]
      dist = crowdingDistance(p)
      merge!(pop.distances, dist)
    end
  end
  
  #last crowding distance is kept local (it will be recalculated after selecting k solutions)
  #keep in mind that the order of fronts[end] is fixed, keep correspondance between lastfront and it
  const lastFront = map(x->x.fitness, pop.individuals[fronts[end]])
  const dist = crowdingDistance(lastFront)
  
  #create mapping fitness => index
  fitnessToIndex = Dict{Vector{Int}, Vector{Int}}()
  for i = 1:(length(lastFront))
    fitnessToIndex[lastFront[i]] = push!(get(fitnessToIndex, lastFront[i], Int[]), i)
  end
  
  #F is a list of fitness sorted by decreasing crowding distance
  F = sort(lastFront, rev = true)
  chosenOnes = Int[]
  j = 1
  while length(chosenOnes) != k
    len = length(fitnessToIndex(F[j]))
    #more than one solution with same fitness
    if len > 1
      sample = rand([1:len])
      index = fitnessToIndex[F[j]][sample]
      push!(consenOnes, index)
      #individuals can be picked only once
      deleteat!(fitnessToIndex[F[j]], sample)
    #only one solution with this fitness
    else
      index = fitnessToIndex[F[j]][1]
      push!(chosenOnes, index)
      deleteat!(F, j)
    end
    j= (j+1)%length(F)
  end
  
  #assign the crowding distance to the now confirmed last front
  p = pop.individuals[chosenOnes]
  dist = crowdingDistance(p)
  merge!(pop.distances, dist)
#   return 
end




function UFTournSelection(L::population)

end


 

#-------------------------helper methods-------------------------
function nonDominatedCompare (a::Vector, b::Vector, comparator = >)
  # a > b --> 0: a==b, 1: a>b, -1: b>a, pairwise vector comparison
  #nonDominatedCompare([0, 0, 2], [0, 0, 1]) = 1
  #nonDominatedCompare([0, 0, 1], [0, 1, 0]) = 0 
  #nonDominatedCompare([1, 0, 1], [1, 1, 1]) =-1
  @assert length(a) == length(b)
  AdomB=false
  BdomA=false
  for i in zip(a,b)
    if i[1] != i[2]
      if(comparator(i[1], i[2]))
	AdomB = true
      else
	BdomA = true
      end
    end
    #immediate return if nondominated
    if AdomB && BdomA
      return 0
    end
  end
  #return result
  if AdomB
    return 1
  end
  if BdomA
    return -1
  end
  if !AdomB && !BdomA
    return 0
  end
end



function evaluateAgainstOthers(pop::population, index::Int, compare_method = nonDominatedCompare)
  #compare the object at index with the rest of the vector
  count = 0
  dominatedby = Int[]
  indFit = pop.individuals[index].fitness
  #before the index
  if index!= 1
    for i = 1: (index-1)
      if compare_method(pop.individuals[i].fitness, indFit) == 1
        count += 1
        push!(dominatedby, i)
      end
    end
  end
  #after the index
  if index != length(pop.individuals)
    for i = (index+1):length(pop.individuals) #exclude the index
      if compare_method(pop.individuals[i].fitness, indFit) == 1
        count += 1
        push!(dominatedby, i)
      end
    end
  end
  #output is sorted
  return (index, count, dominatedby)
end



function fastDelete(values::Vector, deletion::Vector)
  #helper, inside non dominated sorting
  #both vectors are sorted, start in natural, start > 0, nonempty
  @assert values[1] > 0
  @assert deletion[1] > 0
  result = Int[]
  indexDel = 1
  for i in values
    #iterate to the next "good" index, value > or = to i
    while (deletion[indexDel] < i) && (indexDel < length(deletion))
      indexDel += 1
    end
    if i!=deletion[indexDel]
      push!(result, i)
    end
  end
  return result
end



function uniqueFitness(a::Array{(Int,Array{Int,1}),1})
  #helper used in crowdingDistance()
  #map fitness to index (get unique fitness)
  equal = Dict{Vector{Int}, Vector{Int}}()
  for i in a
    equal[i[2]] = push!(get(equal, i[2], Int[]), i[1])
  end
  return equal
end

  
  
#--module end
#end