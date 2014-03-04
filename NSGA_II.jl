

#BEGIN readme
#Simple implementation of the NSGA-II multiobjective
#genetic algorithm.

#Based on "Revisiting the NSGA-II crowding-distance computation"
#END



#BEGIN imports

#END



#BEGIN exports

#END



#BEGIN type definitions
immutable solution
  units::Vector
  fitness::Vector
  
  function solution(units::Vector, fitnessValues::Vector)
    @assert length(units) != 0
    @assert length(fitnessValues) != 0
    self = new(units, fitnessValues)
  end
  
  function solution(units::Vector, fitnessFunction)
    @assert length(units) != 0
    self = new(units, fitnessFunction(units))
  end
end



type population
  solutions::Vector{solution}
  # {fitness => (Front, CrowdingDistance)}
  distances::Dict{Vector, (Int, FloatingPoint)}
  
  function population()
    #initialize empty
    self = new(solution[], Dict{Vector, (Int, FloatingPoint)}())
  end
  
  function population(solutions::Vector{solution})
    #initialize with solutions but no distances
    @assert length(solutions) != 0
    d = Dict{Vector, (Int, FloatingPoint)}()
    self = new(solutions, d)
  end
  
  function population(solutions::Vector{solution}, distances::Dict{Vector, (Int, FloatingPoint)})
    #initialize with solutions and distances
    @assert length(solutions) != 0
    @assert length(distances) != 0
    self = new(solutions, distances)
  end
end



#Hall of Fame is a special population
typealias hallOfFame population



#END



#BEGIN genetic algorithm main methods


#BEGIN nonDominatedSort
function nonDominatedSort(pop::population)
  #sort a population into m nondominating fronts (1 = best, m = worst)
  #until at least half the original number of solutions are added
  
  #get number of solutions to keep
  populationSize = length(pop.solutions)
  cutoff = ceil(populationSize/2)

  #values = (index, domination count, identity of dominators)
  values = (Int, Int, Vector{Int})[]
  for i = 1:populationSize
    push!(values, evaluateAgainstOthers(pop, i, nonDominatedCompare))
  end
  
  #hold indices of the n solutions into m 
  result = Vector{Int}[]
  
  while length(values) > cutoff
    #find the nondominated solutions
    front = filter(x->x[2] == 0, values)
    #get indices
    frontIndices = map(x->x[1], front)
    #add them to the current front
    push!(result, frontIndices)
    
    
    #find the dominated solutions
    values = filter(x->x[2] != 0, values)
    
    #remove the latest front indices from the values
    for i = 1:length(values)
      #delete the last front from the indices
      substracted = fastDelete(values[i][3], frontIndices)
      #substract the difference of cardinality
      values[i] = (values[i][1], length(substracted), substracted)
    end
    
  end
  
  return result
end
#END 



#BEGIN addToHallOfFame
function addToHallOfFame(wholePopulation::population, firstFrontIndices::Vector{Int}, bestPop::hallOfFame)
  #solutions from first front
  firstFront = wholePopulation.solutions[firstFrontIndices]
  
  #add them to the hall of fame
  for i in firstFront
    push!(bestPop.solutions, pop[i])
  end
  
  #acquire the domination values
  values = (Int, Int, Array{Int,1})[]
  for i=1:length(bestPop.solutions)
    push!(values, evaluateAgainstOthers(bestPop, i, nonDominatedCompare))
  end
  
  #get first front solutions
  firstFront = filter(x->x[2]==0, values)
  
  #get indices
  firstFrontIndices = map(x->x[1], firstFront)
  
  #substitute in hall of fame
  bestPop.individuals = bestPop[firstFrontIndices]
end
#END



#BEGIN crowdingDistance
function crowdingDistance(wholePopulation::population, frontIndices::Vector{Int}, frontID::Int, update::Bool)
  #calculate and modify the crowding and front value in the whole population for a certain front
  #fetch fitness from the front
  front = map(x->x.fitness, wholePopulation.solutions[frontIndices])
  
  #map {fitness => crowdingDistance}
  fitnessToCrowding = Dict{Vector, (Int, FloatingPoint)}()
  for v in front
    fitnessToCrowding[v] = (frontID, 0.0)
  end

  #get information about unique fitness
  fitKeys = collect(keys(fitnessToCrowding))
  numUniqueFitness = length(fitKeys)
  fitnessVectorSize = length(fitKeys[1])
  
  #reverse sort fitness vectors for each fitness value
  sorts = Vector{Vector}[]
  maxFitnesses = Number[]
  minFitnesses = Number[]
  
  #get max, min and range for each objective
  for i = 1:fitnessVectorSize
    tmp = sort(fitKeys, by = x->x[i], rev = true)
    push!(maxFitnesses, tmp[1][i])
    push!(minFitnesses, tmp[end][i])
    push!(sorts, tmp)
  end
  rangeSize = map(x->(x[1]-x[2]), zip(maxFitnesses, minFitnesses))
  
  #assign infinite crowding distance to maximum and minimum fitness vectors for each objective 
  map(x->fitnessToCrowding[x[end]] = (frontID, Inf), sorts)
  map(x->fitnessToCrowding[x[1]]   = (frontID, Inf), sorts)

  #assign crowding distances to the other fitness vectors for each objectives
  for i = 1:fitnessVectorSize
    for j = 2:(numUniqueFitness-1)
      fitnessToCrowding[sorts[i][j]] = (frontID, (fitnessToCrowding[sorts[i][j]][2] + (((sorts[i][j-1][i] - sorts[i][j+1][i]))/rangeSize[i])))
    end
  end
  
  #add computed front id and crowding distances to population dict
  if update == true
    merge!(wholePopulation.distances, fitnessToCrowding)
  end
  
  return fitnessToCrowding
end
#END



#BEGIN lastFrontSelection
function lastFrontSelection(wholePopulation::population, lastFrontIndices::Vector{Int}, lastFrontId::Int, k::Int)
  @assert 0 < k <= length(lastFrontIndices)

  #map {fitness => crowding distance}
  fitnessToCrowding = crowdingDistance(wholePopulation, lastFrontIndices, -1, false)

  #map {fitness => indices}
  fitnessToIndex = Dict{Vector, Vector{Int}}()
  for i = 1:length(lastFrontIndices)
    fitnessAtIndex = wholePopulation.solutions[lastFrontIndices[i]].fitness
    fitnessToIndex[fitnessAtIndex] = push!(get(fitnessToIndex, fitnessAtIndex, Int[]), lastFrontIndices[i])
  end
  
  #sort fitness by decreasing crowding distance
  fitnessByCrowding = sort(collect(keys(fitnessToCrowding)), by = k->fitnessToCrowding[k], rev = true)
  
  #choose solutions by iterating through unique fitness list in decreasing order of crowding distance
  chosenOnes = Int[]
  j = 1
  while length(chosenOnes) != k
    len = length(fitnessToIndex[fitnessByCrowding[j]])
    
    if len > 1 #multiple solutions with same fitness
      sample = rand(1:len)
      index = fitnessToIndex[fitnessByCrowding[j]][sample]
      push!(chosenOnes, index)
      #solutions can be picked only once
      deleteat!(fitnessToIndex[fitnessByCrowding[j]], sample)
    
    else #single solution with this fitness
      index = fitnessToIndex[fitnessByCrowding[j]][1]
      push!(chosenOnes, index)
      deleteat!(fitnessByCrowding, j)
    end
    
    j += 1
    #wrap around
    if j>length(fitnessByCrowding)
      j = 1
    end
  end
  
  #get the new crowding distance values for the last front and push it to the whole population
  crowdingDistance(wholePopulation, chosenOnes, lastFrontID, true)
end
#END



#BEGIN UFTournSelection
function UFTournSelection(pop::population)
  #unique fitness based tournament selection
  cardinality = length(pop.solutions)
  result = solution[]
  
  #map {fitness => indices}
  fitnessToIndex = Dict{Vector, Vector{Int}}()
  for i = 1:cardinality
    fitnessToIndex[pop.solutions[i]] = push!(get(fitnessToIndex, pop.solutions[i], Int[]), i)
  end
  
  #keys of map {fitness => indices}
  fitnessKeys = collect(keys(fitnessToIndex))
  
  #edge case : only one fitness
  if length(fitnessToIndex) == 1
    return pop
  end
  
  #solutions selected to parent new population
  chosen = Vector[]
  
  while length(chosen) != cardinality
    k = min((2*(cardinality - length(chosen))), length(fitnessToIndex))
    
    #get k fitnesses and associate to their (front, crowding)
    candidateFitnesses = SAMPLE(fitnessKeys, k)
    vals = map(x->pop.distances[x], candidateFitnesses)
    
    #push the chosen fitnesses in the chosen array
    for i in Range(1, 2, k)
      push!(chosen,  candidateFitnesses[i + crowdedCompare(vals[i], vals[i+1])])
    end
    
  end
  
  #select a random 
  S = map(x -> fitnessToIndex[x][rand(1 : length(fitnessToIndex[x]))], chosen)

  return S
end
#END



#BEGIN generateOffsprings
function generateOffsprings(parents::Vector{solution}, 
                            probabilityOfCrossover::FloatingPoint,
                            probabilityOfMutation::FloatingPoint,
                            evaluationFunction,
                            alleles,
                            mutationOperator = uniformMutate,
                            mutationStrength = 0.05,
                            crossoverOperator = halfUniformCrossover)
  
  #parents -> children
  children = population()
  popSize = length(parents)
  
  #choice vectors
  mut   = map(x->x <= probabilityOfMutation,  rand(length(parents)))
  cross = map(x->x <= probabilityOfCrossover, rand(length(parents)))
  
  #zip them together
  choices = zip(cross, mut)
  
  for i = 1:popSize  
    #no change implies no re-evaluation
    if choices[i] == [false, false]
      push!(children.solutions, parents[i])
      
    else
      #else proceed through operators and evaluate new solution
      newUnits = parents[i].units
      
      #crossover operation
      if choices[i][1] == true
        #find second parent to mix with
        secondParentIndex = rand(1:(popSize-1))
        
        if secondParentIndex >= i
          secondParentIndex += 1
        end
        
        newUnits = halfUniformCrossover(newUnits, parents[secondParentIndex].units)
      end
    
      #mutate operation
      if choices[i][2] == true
        newUnits = mutationOperator(newUnits, mutationStrength, alleles)
      end
      
      push!(children.solutions, solution(newUnits, evaluationFunction))
    end
  end
  
  return children
end
#END

#END



#BEGIN helper methods



#BEGIN SAMPLE
function SAMPLE(L::Vector, k::Int)
  #take k elements from L without replacing
  result = {}
  len = length(L)
  Lcopied = deepcopy(L)
  
  if k == len
    return Lcopied
  end
  
  for i = 1:k
    randIndex = rand(1:len)
    push!(result, Lcopied[randIndex])
    deleteat!(Lcopied, randIndex)
    len -= 1
  end
  return result
end
#END



#BEGIN crowdedCompare
function crowdedCompare(valueA::(Int, FloatingPoint), valueB::(Int, FloatingPoint))
  @assert valueA[1]>0
  @assert valueB[1]>0
  @assert valueA[2]>=0
  @assert valueA[2]>=0
  # A rank < B rank
  if valueA[1] < valueB[1]
    return 0
  # B rank < A rank
  elseif valueA[1] > valueB[1]
    return 1
  # A dist > B dist
  elseif valueA[2] > valueB[2]
    return 0
  # B dist > A dist
  elseif valueA[2] < valueB[2]
    return 1
  # A == B, choose either
  else
    return rand(0:1)
  end
end  
#END



#BEGIN nonDominatedCompare
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
#END



#BEGIN evaluateAgainstOthers
function evaluateAgainstOthers(pop::population, 
                               index::Int, 
                               compare_method = nonDominatedCompare)
  #compare object at index with rest of the vector
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
#END



#BEGIN fastDelete
function fastDelete(values::Vector, deletion::Vector)
  #helper, inside non dominated sorting
  #both vectors are sorted, start in natural, start > 0, nonempty
  @assert values[1] > 0
  @assert deletion[1] > 0
  @assert issorted(values)
  @assert issorted(deletion)
  result = Int[]
  indexDel = 1
  for i in values
    #iterate to the next valid index, value >= to i
    while (deletion[indexDel] < i) && (indexDel < length(deletion))
      indexDel += 1
    end
    if i!=deletion[indexDel]
      push!(result, i)
    end
  end
  return result
end
#END



#BEGIN  uniformCrossover
function uniformCrossover(units1::Vector, units2::Vector)
  @assert length(units1) == length(units2) != 0
  
  newUnits = deepcopy(units1)
  units2 = deepcopy(units2)
  
  for i = 1:length(newUnits)
    if rand() < 0.5
      newUnits[i] = units2[i]
    end
  end
  
  return newUnits
end
#END



#BEGIN halfUniformCrossover
function halfUniformCrossover(units1::Vector, units2::Vector)
  #from wikipedia:
  #In the half uniform crossover scheme (HUX), exactly half of the nonmatching
  #bits are swapped. Thus first the Hamming distance (the number of differing bits) 
  #is calculated. This number is divided by two. The resulting number is how many 
  #of the bits that do not match between the two parents will be swapped.
  @assert length(units1) == length(units2) != 0
  result = deepcopy(units1)
  units2 = deepcopy(units2)
  
  #find how many are matching
  matching = map(x->x[1]==x[2] ? 1 : 0, zip(result, units2))
  sumDifferent = length(matching) - reduce(+, matching)
  
  #need to swap at least half of the number of non matching
  toSwap = ceil(sumDifferent / 2)
  
  #if not matching, 0.5 probability of exchange
  for i = 1:length(matching)
    if matching[i] == 0
      if rand() < 0.5
        result[i] = units2[i]
      end
    end
  end
  
  return result
end
#END



#BEGIN onePointCrossover
function onePointCrossover(units1::Vector, units2::Vector)
  #swap at one point
  @assert length(units1) == length(units2) != 0
  result = deepcopy(units1)
  units2 = deepcopy(units2)
  
  #find point
  p = rand(1:length(result))
  
  #swap the values after the point
  for i = p:length(result)
    result[i] = units2[i]
  end
  
  return result
end
#END



#BEGIN uniformMutate
function uniformMutate(individualUnits::Vector, probability::FloatingPoint, alleles::Vector{Vector})
  #copy the individual units
  newUnits = deepcopy(individualUnits)
  
  #for each unit, mutate if random is inferior to given probability
  for i = 1:length(newUnits)
    if(rand() < probability)
      newUnits[i] = alleles[i][rand(1:length(alleles[i]))]
    end
  end
  
  return newUnits
end
#END



#END


