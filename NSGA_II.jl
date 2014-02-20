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


#every fitness has an associated crowding distance too
#the key is not the individual but its fitness
#	-same fitness implies same front
#	-same fitness implies same crowding (but is treated correctly in the last front selection)

#\todo: add an assertion that solutions must have the same size and type of fitness vector 
type population
  solutions::Vector{solution}
  distances::Dict{Vector, (Int, FloatingPoint)} #distances[fitness] = (Front, CrowdingDistance)
end



typealias hallOfFame population



#-------------------------genetic algorithm methods-------------------------
function nonDominatedSort(pop::population) 
  #multi-objective optimization using evolutionnary algorithms p.43
  #input: 2N population
  #output: at least N solutions, the output likely has one front too much
  
  #get cutoff
  len = length(pop.solutions)
  cutoff = div(len, 2) + len%2
  
  #step 1
  #evaluate the whole population
  values = (Int, Int, Array{Int,1})[]
  for i = 1:len
    push!(values, evaluateAgainstOthers(pop, i, nonDominatedCompare))
  end
  
  #step 2 
  #nondominated sort
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



function addToHallOfFame(wholePopulation::population, firstFrontIndices::Vector{Int}, bestPop::hallOfFame)
  #get the first front from the whole population
  firstFront = wholePopulation.solutions[firstFrontIndices]
  
  #add the first front to the hall of fame
  for i in firstFront
    push!(bestPop.solutions, pop[i])
  end
  
  #acquire the domination values
  values = (Int, Int, Array{Int,1})[]
  for i=1:length(bestPop.solutions)
    push!(values, evaluateAgainstOthers(bestPop, i, nonDominatedCompare))
  end
  
  #get the first front
  firstFront = filter(x->x[2]==0, values)
  
  #get first front indices
  firstFrontIndices = map(x->x[1], firstFront)
  
  #reassign the hall of fame
  bestPop.individuals = bestPop[firstFrontIndices]
end




function crowdingDistance(wholePopulation::population, frontIndices::Vector{Int}, frontID::Int, update::Bool)
  #calculate and modify the crowding and front value in the whole population for a certain front
  #don't use if on the last front, it doesn't make any sense
  #step 1
  #fetch the individual fitness from the front
  front = map(x->x.fitness, wholePopulation.solutions[frontIndices])
  
  #step 2
  #create a dict {fitness => crowdingDistance}
  fitnessToCrowding = Dict{Vector, (Int, FloatingPoint)}()
  for i in front
    fitnessToCrowding[i] = (frontID, 0.0)
  end

  fitKeys = collect(keys(fitnessToCrowding))
  
  #step 3 
  #do the crowding distance
  numUniqueFitness = length(fitKeys)
  fitnessVectorSize = length(fitKeys[1])
  
  #step 4
  #reverse sort all the unique fitness vectors per each objective value
  sorts = {} #todo: static type this one
  for i = 1:vectorSize
    push!(sorts, sort(fitKeys, by = x->x[i]), rev = true)
  end
  
  #step 5
  #get the max and min of each objective
  maxFitnesses = map(x->x[1],   sorts)
  minFitnesses = map(x->x[end], sorts)
  rangeSize = map(x->x[1]-x[2], zip(maxFitnesses, minFitnesses))
  
  #step 6
  #assign infinite crowding distance to maximum and minimum fitness vectors for each objective 
  map(x->fitnessToCrowding[x[end]] = (frontId, Inf), sorts)
  map(x->fitnessToCrowding[x[1]]   = (frontId, Inf), sorts)
  
  #step 7
  #assign crowding distances to the other fitness vectors for each objectives
  for i = 1:fitnessVectorSize
    for j = 2:(numUniqueFitness-1)
      fitnessToCrowding[sorts[i][j]] = (frontID, (fitnessToCrowding[sorts[i][j]][2] + (((sorts[i][j-1] - sorts[i][j+1]))/rangeSize[i])))
    end
  end
  
  #step 8
  #merge the front crowding distance dictionary with the one of the population
  #this assigns both the crowding distance and the front value
  if update == true
    merge!(wholePopulation.distances, fitnessToCrowding)
  end
  
  #step 9
  #return
  return fitnessToCrowding
end




function lastFrontSelection(wholePopulation::population, lastFrontIndices::Vector{Int}, lastFrontId::Int, k::Int)
  #select k solutions from the last front 
  
  #create mapping fitness => crowding
  fitnessToCrowding = crowdingDistance(wholePopulation, lastFrontIndices, -1, false)
  
  #create mapping fitness => lastFrontIndices index
  fitnessToIndex = Dict{Vector{Int}, Vector{Int}}()
  for i = 1:length(lastFrontIndices)
    fitnessAtIndex = wholePopulation.solutions[lastFrontIndices[i]].fitness
    fitnessToIndex[fitnessAtIndex] = push!(get(fitnessToIndex, fitnessAtIndex, Int[]), lastFrontIndices[i])
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
      #solutions can be picked only once
      deleteat!(fitnessToIndex[F[j]], sample)
    #only one solution with this fitness
    else
      index = fitnessToIndex[F[j]][1]
      push!(chosenOnes, index)
      deleteat!(F, j)
    end
    j= (j+1)%length(F)
  end
  
  #get the new crowding distance values for the last front and push it to the whole population
  crowdingDistance(wholePopulation, chosenOnes, lastFrontID, true)
  
   return chosenOnes
end



function sample(L::Vector, k::Int)
  #take k elements from L without replacing
  L2 = copy(L)
  result = {}
  for i = 1:k
    randIndex = rand(1:length(L2))
    push!(result, L2[randIndex])
    deleteat!(L2, randIndex)
  end
  return result
end  



function UFTournSelection(L::population)
  #unique fitness based tournament selection
  N = length(L.solutions)
  
  #map fitnesses to indices
  values = (Int, Vector)[]
  for i=1:length(population.solutions)
    push!(values, (i, population.solutions[i].fitness))
  end
  
  fitnessToIndex = uniqueFitness(values)
  
  #extreme case : all fitnesses are equal
  if length(fitnessToIndex) == 1
    return L
  end
  
  #
  S = {}
  #while the size of S is not equal to N...
  while length(S) != N
    k = min(2*(N - length(S)), length(fitnessToIndex))
  end
  
  
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
  indFit = pop.solutions[index].fitness
  #before the index
  if index!= 1
    for i = 1: (index-1)
      if compare_method(pop.solutions[i].fitness, indFit) == 1
        count += 1
        push!(dominatedby, i)
      end
    end
  end
  #after the index
  if index != length(pop.solutions)
    for i = (index+1):length(pop.solutions) #exclude the index
      if compare_method(pop.solutions[i].fitness, indFit) == 1
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




  
  
#--module end
#end