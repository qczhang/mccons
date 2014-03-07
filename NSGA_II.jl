module NSGA_II




#------------------------------------------------------------------------------
#BEGIN readme



#Implementation of the NSGA-II multiobjective
#genetic algorithm as described in:

# Revisiting the NSGA-II crowding-distance computation
# Felix-Antoine Fortin
# Marc Parizeau
# Universite Laval, Quebec, PQ, Canada
# GECCO '13 Proceeding of the fifteenth annual conference on Genetic and evolutionary computation conference
# Pages 623-630



#END   readme
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN imports



require("geneticAlgorithmOperators")



#END   imports
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN type definitions



immutable solution
  #unit, individual, basic block of the solution
  units::Vector
  fitness::Vector
  
  function solution(units::Vector, fitnessValues::Vector)
    #fitness value is precomputed
    @assert length(units) != 0
    @assert length(fitnessValues) != 0
    new(units, fitnessValues)
  end
  
  function solution(units::Vector, fitnessFunction::Function)
    #fitness value is to be computed
    @assert length(units) != 0
    new(units, fitnessFunction(units))
  end
end



type population
  #the compound of all solutions
  #includes a mapping of fitness values to crowding distance
  
  solutions::Vector{solution}
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
  
  function population(solutions::Vector{solution}, 
                      distances::Dict{Vector, (Int, FloatingPoint)})
    #initialize with solutions and distances
    @assert length(solutions) != 0
    @assert length(distances) != 0
    self = new(solutions, distances)
  end
end



#Hall of Fame is a special population
typealias hallOfFame population



#END   type definitions
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN NSGA-II main methods



#BEGIN nonDominatedSort
function nonDominatedSort(pop::population, 
                          comparisonMethod = nonDominatedCompare)
  #sort a population into m nondominating fronts (1 = best, m = worst)
  #until at least half the original number of solutions are added
  
  #get number of solutions to keep
  populationSize = length(pop.solutions)
  cutoff = ceil(populationSize/2)

  #values = (index, domination count, identity of dominators)
  values = (Int, Int, Vector{Int})[]
  for i = 1:populationSize
    push!(values, evaluateAgainstOthers(pop, i, comparisonMethod))
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
function addToHallOfFame(wholePopulation::population, 
                         firstFrontIndices::Vector{Int}, 
                         HallOfFame::hallOfFame)
  #solutions from first front
  firstFront = wholePopulation.solutions[firstFrontIndices]
  
  #add them to the hall of fame
  for i in firstFront
    push!(HallOfFame.solutions, pop[i])
  end
  
  #acquire the domination values
  values = (Int, Int, Array{Int,1})[]
  for i=1:length(HallOfFame.solutions)
    push!(values, evaluateAgainstOthers(HallOfFame, i, nonDominatedCompare))
  end
  
  #get first front solutions
  firstFront = filter(x->x[2]==0, values)
  
  #get indices
  firstFrontIndices = map(x->x[1], firstFront)
  
  #substitute in hall of fame
  HallOfFame.individuals = HallOfFame[firstFrontIndices]
end
#END



#BEGIN crowdingDistance
function crowdingDistance(P::population, 
                          frontIndices::Vector{Int}, 
                          frontID::Int, 
                          update::Bool)
  #calculate and modify the crowding and front value in 
  #the whole population for a certain front
  
  #fetch fitness from the front
  front = map(x->x.fitness, P.solutions[frontIndices])
  
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
  
  #assign infinite crowding distance to maximum and 
  #minimum fitness vectors for each objective 
  map(x->fitnessToCrowding[x[end]] = (frontID, Inf), sorts)
  map(x->fitnessToCrowding[x[1]]   = (frontID, Inf), sorts)

  #assign crowding distances to the other 
  #fitness vectors for each objectives
  for i = 1:fitnessVectorSize
    for j = 2:(numUniqueFitness-1)
      previousValue = fitnessToCrowding[sorts[i][j]][2]
      toAdd = ((sorts[i][j-1][i] - sorts[i][j+1][i])/rangeSize[i])
      fitnessToCrowding[sorts[i][j]] = (frontID, (previousValue + toAdd))
    end
  end
  
  #add computed front id and crowding distances to population dict
  if update == true
    merge!(P.distances, fitnessToCrowding)
  end
  
  return fitnessToCrowding
end
#END



#BEGIN lastFrontSelection
function lastFrontSelection(P::population, 
                            lastFrontIndices::Vector{Int}, 
                            lastFrontID::Int, 
                            k::Int)
  @assert 0 < k <= length(lastFrontIndices)

  #map {fitness => crowding distance}
  fitnessToCrowding = crowdingDistance(P, lastFrontIndices, -1, false)

  #map {fitness => indices}
  fitnessToIndex = Dict{Vector, Vector{Int}}()
  for i = 1:length(lastFrontIndices)
    fitnessAtIndex = wholePopulation.solutions[lastFrontIndices[i]].fitness
    value = get(fitnessToIndex, fitnessAtIndex, Int[])
    fitnessToIndex[fitnessAtIndex] = push!(value, lastFrontIndices[i])
  end
  
  #sort fitness by decreasing crowding distance
  fitnessByCrowding = sort(collect(keys(fitnessToCrowding)), 
                           by = k->fitnessToCrowding[k], 
                           rev = true)
  
  #choose solutions by iterating through unique fitness list 
  #in decreasing order of crowding distance
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
  
  #get the new crowding distance values for the 
  #last front and push it to the whole population
  crowdingDistance(P, chosenOnes, lastFrontID, true)
  
  #return the indices of the chosen solutions on the last front
  return chosenOnes
end
#END



#BEGIN UFTournSelection
function UFTournSelection(pop::population)
  #Unique Fitness based Tournament Selection
  
  #select across entire range of fitnesses to avoid 
  #bias by reoccuring fitnesses
  
  #size of the population
  cardinality = length(pop.solutions)
  result = solution[]
  
  #map {fitness => indices}
  fitnessToIndex = Dict{Vector, Vector{Int}}()
  for i = 1:cardinality
    value = get(fitnessToIndex, pop.solutions[i], Int[])
    fitnessToIndex[pop.solutions[i]] = push!(value, i)
  end
  
  #keys of map {fitness => indices}
  fitnessKeys = collect(keys(fitnessToIndex))
  
  #edge case : only one fitness
  if length(fitnessToIndex) == 1
    return pop
  end
  
  #solutions selected to parent new population
  \fiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiixxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx vals
  vals = 
  while length(chosen) != cardinality
    k = min((2*(cardinality - length(chosen))), length(fitnessToIndex))
    
    #get k fitnesses and associate to their (front, crowding)
    candidateFitnesses = SAMPLE(fitnessKeys, k)
    vals = map(x->pop.distances[x], candidateFitnesses)
    
    #Choose n fittest out of 2n
    #by comparing pairs of neighbors
    
  end
  
  
  chosen = Vector[]
  for i in Range(1, 2, k)
    push!(chosen, candidateFitnesses[i + crowdedCompare(vals[i], vals[i+1])])
  end
    
  #randomly choose a solution from the solutions sharing the fitness
  return map(x -> fitnessToIndex[x][rand(1 : length(fitnessToIndex[x]))], chosen)
end
#END



#BEGIN initializePopulation
function initializePopulation{T}(alleles::Vector{Vector{T}}, 
                                 fitnessFunction::Function, 
                                 n::Int)
  #
  @assert n>0
  result = population()
  index = 1
  while index <= n
    units = {}
    for i in alleles
      push!(units, i[rand(1:length(i))])
    end
    sol = solution(units, fitnessFunction)
    push!(result.solutions, sol)
    index+=1
  end
  return result
end
#END



#BEGIN generateOffsprings
function generateOffsprings(parents::Vector{solution}, 
                            probabilityOfCrossover::FloatingPoint,
                            probabilityOfMutation::FloatingPoint,
                            evaluationFunction::Function,
                            alleles,
                            mutationOperator = geneticAlgorithmOperators.uniformMutate,
                            mutationStrength = 0.05,
                            crossoverOperator = geneticAlgorithmOperators.halfUniformCrossover)
  
  #parents -> children
  children = population()
  popSize = length(parents)
  
  #deciding who is mutating and having crossovers
  isMutated   = map(x->x <= probabilityOfMutation,  rand(length(parents)))
  hasCrossover= map(x->x <= probabilityOfCrossover, rand(length(parents)))
  
  evolutionaryEvents = zip(isMutated, hasCrossover)
  
  for i = 1:popSize
  \fiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiixxxxxxxxxxxxxxxxxxxxxxxxxxxxxx do it with one push
    #no change implies no re-evaluation
    if evolutionaryEvents[i] == [false, false]
      push!(children.solutions, parents[i])
      
    else
      #else proceed through operators and evaluate new solution
      
      #fitness function will evaluate the units
      newUnits = deepcopy(parents[i].units)
      
      #recombination (crossover)
      if evolutionaryEvents[i][1] == true
        #randomly choose second parent
        secondParentIndex = rand(1:(popSize-1))
        
        #leave a gap to not select same parent
        if secondParentIndex >= i
          secondParentIndex += 1
        end
        
        #combine two parents units (on which the fitness is based)
        newUnits = crossoverOperator(newUnits, parents[secondParentIndex].units)
      end
    
      #mutation
      if evolutionaryEvents[i][2] == true
        newUnits = mutationOperator(newUnits, mutationStrength, alleles)
      end
      
      
      push!(children.solutions, solution(newUnits, evaluationFunction))
    end
    
    newSolution = solution(newUnits, fitness)
  end
  
  return children
end
#END


#BEGIN NSGA_II_main
function NSGA_II_main{T}(alleles::Vector{Vector{T}},
                         fitnessFunction::Function,
                         populationSize::Int,
                         iterations::Int,
                         probabilityOfCrossover::FloatingPoint,
                         probabilityOfMutation::FloatingPoint)
  @assert populationSize > 0
  @assert iterations > 0
  #main loop of the NSGA-II algorithm
  #executes selection -> breeding until number of iterations is reached
  
  HallOfFame = hallOfFame()
  
  #initialize
  P1 = initializePopulation(alleles, fitnessFunction, n)
  Q1 = initializePopulation(alleles, fitnessFunction, n)
  
  #merge two initial parents
  Pt = population(vcat(P1.solutions, Q1.solutions))
  
  #iterate selection -> breeding
  for i = 1:iterations
    
    #sort into non dominated fronts
    F = nonDominatedSort(Pt)
    
    addToHallOfFame(Pt, F[1], HallOfFame)
    
    #get the last front indices
    lastFront = F[end]
    
    #separate last front from rest
    F = F[1:end-1]
    Pnext = population()
    
    #calculate the crowding distances for all but the last front
    for j = 1:length(F)
      front = F[j]
      crowdingDistance(Pt, front, j, true)
    end
    
    #------------------------
    #how many solutions left to choose from the last front 
    toChoose = length(Pt.solutions) - length(reduce(vcat, F))
    
    #perform last front selection (side effect: add crowding distance for last front fitnesses)
    indicesLastFront = lastFrontSelection(Pt, lastFront, length(F) + 1, toChoose)
    
    #-------------------------
    
    #concatenate the list of all individuals to be parents
    chosenIndices = vcat(reduce(vcat, F), indicesLastFront)
    
    
    
    #create the Pt+1 population 
    Pt2 = population(Pt.solutions[chosenIndices], Pt.distances)
    
    
    #select based on crowding comparison
    Pt2prime = UFTournSelection(Pt2)
    
    #generate offsprings
    Pt = generateOffsprings(Pt2prime, probabilityOfCrossover, probabilityOfMutation,  )
    
    #merge earlier and actual pr
  end
  
  return HallOfFame
end
#END



#END   NSGA-II main methods
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
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
function crowdedCompare(valueA::(Int, FloatingPoint), 
                        valueB::(Int, FloatingPoint))
  #crowded comparison operator
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
  # non domination operator
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



#END   helper methods
#------------------------------------------------------------------------------




end

