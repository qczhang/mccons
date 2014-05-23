module NSGA_II
using ProgressMeter

#------------------------------------------------------------------------------
#BEGIN readme



# Implementation of genetic algorithm mutation and crossover operators
# consult:

# Wikipedia...

# Genetic Algorithms in Search, Optimization, and Machine Learning
# David E. Goldberg



#END   readme
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#BEGIN crossover operators



#BEGIN  uniformCrossover
function uniformCrossover(genes1::Vector, genes2::Vector)
  @assert length(genes1) == length(genes2) != 0

  newGenes = deepcopy(genes1)
  genes2 = deepcopy(genes2)

  for i = 1:length(newGenes)
    if rand() < 0.5
      newGenes[i] = genes2[i]
    end
  end
  return newGenes
end
#END




#BEGIN onePointCrossover
function onePointCrossover(genes1::Vector, genes2::Vector)
  #exchange genes after a certain point in the chromosome
  @assert length(genes1) == length(genes2) != 0
  result = deepcopy(genes1)
  genes2 = deepcopy(genes2)
  
  #find point
  point = rand(1:length(result))
  
  #swap the values after the point
  for i = point:length(result)
    result[i] = genes2[i]
  end
  
  return result
end
#END



#END   crossover operators
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#BEGIN mutation operators



#BEGIN uniformMutate
function uniformMutate(originalGenes::Vector, alleles::Vector, probability = 0.05)
  #copy the individual genes
  newGenes = deepcopy(originalGenes)
  
  #for each unit, mutate if random is inferior to given probability
  for i = 1:length(newGenes)
    if(rand() < probability)
      newGenes[i] = alleles[i][rand(1:length(alleles[i]))]
    end
  end
  
  return newGenes
end
#END



#END   mutation operators
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN readme



#Implementation of the NSGA-II multiobjective
#genetic algorithm as described in:

# Revisiting the NSGA-II crowding-distance computation
# Felix-Antoine Fortin
# Marc Parizeau
# Universite Laval, Quebec, PQ, Canada
# GECCO '13 Proceeding of the fifteenth annual conference on Genetic 
# and evolutionary computation conference
# Pages 623-630



#END   readme
#------------------------------------------------------------------------------







#------------------------------------------------------------------------------
#BEGIN type definitions


immutable individual
  #unit, individual, basic block of the solution
  genes::Vector
  fitness::Vector
  
  function individual(genes::Vector, fitnessValues::Vector)
    #fitness value is precomputed
    @assert length(genes) != 0
    @assert length(fitnessValues) != 0
    new(genes, fitnessValues)
  end
  
  function individual(genes::Vector, fitnessFunction::Function)
    #fitness value is to be computed
    @assert length(genes) != 0
    new(genes, fitnessFunction(genes))
  end
end



type population
  #the compound of all individuals
  #includes a mapping of fitness values to crowding distance
  
  individuals::Vector{individual}
  distances::Dict{Vector, (Int, FloatingPoint)}
  
  function population()
    #initialize empty
    self = new(individual[], Dict{Vector, (Int, FloatingPoint)}())
  end
  
  function population(individuals::Vector{individual})
    #initialize with individuals but no distances
    @assert length(individuals) != 0
    d = Dict{Vector, (Int, FloatingPoint)}()
    self = new(individuals, d)
  end
  
  function population(individuals::Vector{individual},
                      distances::Dict{Vector, (Int, FloatingPoint)})
    #initialize with individuals and distances
    @assert length(individuals) != 0
    @assert length(distances) != 0
    self = new(individuals, distances)
  end
end


#hall of fame is a special population to keep
#the best individuals of all generations
typealias hallOfFame population



#END   type definitions
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN NSGA-II main methods



#BEGIN initializePopulation
function initializePopulation{T}(alleles::Vector{Vector{T}}, 
                                 fitnessFunction::Function, 
                                 n::Int)
  #from alleles and a fitness function, initialize a population
  #of n individuals
  @assert n>0
  P = population()
  index = 1
  while index <= n
    genes = T[]
    for i in alleles
      push!(genes, i[rand(1:length(i))])
    end
    sol = individual(genes, fitnessFunction)
    push!(P.individuals, sol)
    index+=1
  end
  return P
end
#END



#BEGIN nonDominatedSort
function nonDominatedSort(P::population, 
                          comparisonMethod = nonDominatedCompare)
  #sort a population into m nondominating fronts (1 = best, m = worst)
  #until at least half the original number of individuals are added
  
  #get number of individuals to keep
  populationSize = length(P.individuals)
  cutoff = populationSize/2

  #evaluate all against all individuals to get domination count
  #and identity of dominators (for dominated individuals)
  #values = (index, domination count, identity of dominators)
  values = (Int, Int, Vector{Int})[]
  for i = 1:populationSize
    push!(values, evaluateAgainstOthers(P, i, comparisonMethod))
  end
  
  #result is a vectors of vectors of indices (direct mapping to individuals)
  result = Vector{Int}[]
  
  #add fronts until there are at least n individuals in the fronts
  while length(values) > cutoff
    #first, we must find which individuals are not dominated and 
    #separate them from the rest
    currentFront = Int[]
    dominatedValues = (Int, Int, Vector{Int})[]
    
    for i in values
      if i[2] == 0
        #the individual is dominating, we add its index to frontIndices
        push!(currentFront, i[1])
      else
        #the individual is dominated
        push!(dominatedValues, i)
      end
    end
    
    #push the current front to the result
    push!(result, currentFront)
    
    #update the values, removing the dominating
    values = dominatedValues
    
    #remove the latest front indices from the values
    for i = 1:length(values)
      #remove indices from the current front
      substracted = fastDelete(values[i][3], currentFront)
      #substract the difference of cardinality
      values[i] = (values[i][1], length(substracted), substracted)
    end
  end
  
  return result
end
#END 



#BEGIN crowdingDistance

#the crowding distance measures the proximity of a
#solution to its immediate neighbors of the same front. it is used
#to preserve diversity, later in the algorithm.

#this particular version uses fitness directly to avoid bias
#problems introduced by many solutions sharing a same fitness
#(see article referred in the readme for more details).
function crowdingDistance(P::population,
                          frontIndices::Vector{Int},
                          frontID::Int,
                          update::Bool)
  #update flag signals to add or not the distances calculated to P.distances
  
  #get all the fitnesses from the individuals in the front
  front = map(x->x.fitness, P.individuals[frontIndices])
  
  #add these fitnesses to a dict (non unique fitnesses will disappear)
  #and assign an initial crowding distance of 0
  #{fitness => crowdingDistance}
  fitnessToCrowding = Dict{Vector, (Int, FloatingPoint)}()
  for v in front
    fitnessToCrowding[v] = (frontID, 0.0)
  end

  #get how many fitnesses we have and how many objectives they each have
  fitKeys = collect(keys(fitnessToCrowding))
  numFitness = length(fitKeys)
  numOfObjectives = length(fitKeys[1])
  
  #sort in decreasing order the fitness vectors for each objective
  sortedByObjective = Vector{Vector}[]
  maxOfObjective = Number[]
  minOfObjective = Number[]
  
  for i = 1:numOfObjectives
    tmp = sort(fitKeys, by = x->x[i], rev = true)
    push!(maxOfObjective, tmp[1][i])
    push!(minOfObjective, tmp[end][i])
    push!(sortedByObjective, tmp)
  end
  
  #\todoooo fix error Inf - Inf = NaN or 0 / 0
  #get the range of the fitnesses for each objective
  rangeOfObjective = map(x->(x[1]-x[2]), zip(maxOfObjective, minOfObjective))

  #assign infinite crowding distance to maximum and
  #minimum fitness of each objective
  map(x->fitnessToCrowding[x[end]] = (frontID, Inf), sortedByObjective)
  map(x->fitnessToCrowding[x[1]]   = (frontID, Inf), sortedByObjective)

  #assign crowding distances to the other
  #fitness vectors for each objectives
  for i = 1:numOfObjectives
    #edge case here! if range == 0, 0 / 0 will give NaN
    if rangeOfObjective[i] != 0
      for j = 2:(numFitness-1)
        previousValue = fitnessToCrowding[sortedByObjective[i][j]][2]
        toAdd = ((sortedByObjective[i][j-1][i] - sortedByObjective[i][j+1][i])/rangeOfObjective[i])
        @assert toAdd != NaN
        fitnessToCrowding[sortedByObjective[i][j]] = (frontID, (previousValue + toAdd))
      end
    end
  end
  
  #test that all are valid numbers
  for i in keys(fitnessToCrowding)
    val = fitnessToCrowding[i]
    if val[1]==NaN
      println(fitnessToCrowding)
    elseif val[2] == NaN
      println("###################################################BUG BUG#######")
      println(fitnessToCrowding)
      println("#################################################################")
    end
  end
  
  #if we wish to update, the dict of computed distances is merged to the main one
  #in P.distances
  if update == true
    merge!(P.distances, fitnessToCrowding)
  end
  
  return fitnessToCrowding
end
#END



#BEGIN lastFrontSelection

#the last selected front usually contains more individuals
#than needed for the parent population. supposing
#that the previous front contain n-k individuals, we
#need to select k individuals to reach a population
#of n individuals.

#since they do not dominate each other, we must
#select them by crowding distance
function lastFrontSelection(P::population,
                            lastFrontIndices::Vector{Int},
                            k::Int)
  @assert 0 < k <= length(lastFrontIndices)

  #map {fitness => crowding distance}
  fitnessToCrowding = crowdingDistance(P, lastFrontIndices, -1, false)

  #map {fitness => indices}
  fitnessToIndex = Dict{Vector, Vector{Int}}()
  for i = 1:length(lastFrontIndices)
    fitnessAtIndex = P.individuals[lastFrontIndices[i]].fitness
    value = get(fitnessToIndex, fitnessAtIndex, Int[])
    fitnessToIndex[fitnessAtIndex] = push!(value, lastFrontIndices[i])
  end
  
  #sort fitness by decreasing crowding distance
  fitnessByCrowding = sort(collect(keys(fitnessToCrowding)), 
                           by = k->fitnessToCrowding[k], 
                           rev = true)
  
  #choose individuals by iterating through unique fitness list
  #in decreasing order of crowding distance
  chosenOnes = Int[]
  j = 1
  while length(chosenOnes) != k
    len = length(fitnessToIndex[fitnessByCrowding[j]])
    
    if len > 1 #multiple individuals with same fitness
      sample = rand(1:len)
      index = fitnessToIndex[fitnessByCrowding[j]][sample]
      push!(chosenOnes, index)
      #individuals can be picked only once
      deleteat!(fitnessToIndex[fitnessByCrowding[j]], sample)
    
    else #single individual with this fitness
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
  
  #return the indices of the chosen individuals on the last front
  return chosenOnes
end
#END



#BEGIN UFTournSelection
function UFTournSelection(P::population)
  #Unique Fitness based Tournament Selection
  
  #select across entire range of fitnesses to avoid 
  #bias by reoccuring fitnesses
    
  sizeOfPopulation = length(P.individuals)
  
  #associate fitness to indices of individuals
  #map {fitness => indices}
  fitnessToIndex = Dict{Vector, Vector{Int}}()
  for i = 1:sizeOfPopulation
    value = get(fitnessToIndex, P.individuals[i].fitness, Int[])
    fitnessToIndex[P.individuals[i].fitness] = push!(value, i)
  end
  fitnesses = collect(keys(fitnessToIndex))
  
  
  #edge case : only one fitness, return the population as it was
  if length(fitnessToIndex) == 1
    return P.individuals
  end
  
  #else we must select parents
  newParents = individual[]
  
  while length(newParents) != sizeOfPopulation
    #we either pick all the fitnesses and select a random individual from them
    #or select a subset of them. depends on how many new parents we still need to add
    k = min((2*(sizeOfPopulation - length(newParents))), length(fitnessToIndex))

    #sample k fitnesses and get their (front, crowing) from P.distances
    candidateFitnesses = SAMPLE(fitnesses, k)
    frontAndCrowding = map(x->P.distances[x], candidateFitnesses)
    
    #Choose n fittest out of 2n
    #by comparing pairs of neighbors
    chosenFitnesses = Vector[]
    i = 1
    while i < k
      #crowdedCompare returns an offset (0 if first solution is better, 1 otherwise)
      selectedIndex = i + crowdedCompare(frontAndCrowding[i], frontAndCrowding[i+1])
      push!(chosenFitnesses, candidateFitnesses[selectedIndex])
      i += 2
    end
    
    #we now randomly choose an individual from the indices associated with the chosen fitnesses
    for i in chosenFitnesses
      chosenIndex = fitnessToIndex[i][rand(1:length(fitnessToIndex[i]))]
      push!(newParents, P.individuals[chosenIndex])
    end
    
  end
  
  
  return newParents
end
#END



#BEGIN generateOffsprings

#final step of the generation, creates the next population
#from the children templates
function generateOffsprings(childrenTemplates::Vector{individual}, 
                            probabilityOfCrossover::FloatingPoint,
                            probabilityOfMutation::FloatingPoint,
                            evaluationFunction::Function,
                            alleles,
                            mutationOperator,
                            crossoverOperator)
  
  #initialize
  childrenPopulation = population()
  popSize = length(childrenTemplates)
  
  #deciding who is mutating and having crossovers
  willMutate   = map(x->x <= probabilityOfMutation,  rand(length(childrenTemplates)))
  willRecombine= map(x->x <= probabilityOfCrossover, rand(length(childrenTemplates)))
  
  evolutionaryEvents = collect(zip(willMutate, willRecombine))
  
  for i = 1:popSize
    #initialize new genes and a new fitness from childrenTemplates genes and fitness
    newGenes = deepcopy(childrenTemplates[i].genes)
    newFitness = deepcopy(childrenTemplates[i].fitness)
    modified = false

    if evolutionaryEvents[i][1] == true
      modified = true
      #recombination (crossover)
      
      #randomly choose second parent
      secondParentIndex = rand(1:(popSize-1))

      #leave a gap to not select same parent
      if secondParentIndex >= i
        secondParentIndex += 1
      end

      #combine two childrenTemplates genes (on which the fitness is based)
      newGenes = crossoverOperator(newGenes, childrenTemplates[secondParentIndex].genes)
    end
    
    if evolutionaryEvents[i][2] == true
      modified = true
      #mutation
      
      newGenes = mutationOperator(newGenes, alleles)
    end
      
    #if modified, re-evaluate
    if modified
      newFitness = evaluationFunction(newGenes)
    end
    
    #add newly created individuals to the children population
    push!(childrenPopulation.individuals, individual(newGenes, newFitness))
  end

  return childrenPopulation
end
#END



#BEGIN addToHallOfFame
function addToHallOfFame(P::population, 
                         firstFrontIndices::Vector{Int},
                         HallOfFame::hallOfFame,
                         maxSize=400)
  #add the best individuals to the Hall of Fame population to save them for
  #further examination. we merge the first front of the actual population
  #with the rest of the hall of fame to then select the first front of it.
  
  #we know from previous calculation the indices of the best individuals
  firstFront = P.individuals[firstFrontIndices]
  #println("num in hall of fame $(length(HallOfFame.individuals))")
  #println("num external first front $(length(firstFront))")

  
  
  #we add add the best individuals to the Hall of Fame
  for i in firstFront
    push!(HallOfFame.individuals, i)
  end

  
  #elmiminate duplicates (since it is elitist, same individuals may reappear)
  HallOfFame.individuals = unique(HallOfFame.individuals)
  
  #find the first non dominated front, to select the best individuals of the new Hall of Fame
  
  #acquire the domination values
  values = (Int, Int, Array{Int,1})[]
  for i=1:length(HallOfFame.individuals)
    push!(values, evaluateAgainstOthers(HallOfFame, i, nonDominatedCompare))
  end
  
  #get first front individuals
  firstFront2 = filter(x->x[2]==0, values)
  
  #get indices
  firstFrontIndices2 = map(x->x[1], firstFront2)
  firstFrontIndividuals = HallOfFame.individuals[firstFrontIndices2]
  
  fitnesses = unique(map(x->x.fitness, firstFrontIndividuals))
  
  
  #unique genes
  selected = individual[]
  allGenes = Set{Vector}()
  for i in firstFrontIndividuals
    if !(i.genes in allGenes)
      push!(allGenes, i.genes)
      push!(selected, i)
    end
  end

  HallOfFame.individuals = selected
end
#END



#BEGIN main

function main(alleles::Vector,
              fitnessFunction::Function,
              populationSize::Int,
              iterations::Int,
              probabilityOfCrossover = 0.1,
              probabilityOfMutation = 0.05,
              crossoverOperator = uniformCrossover,
              mutationOperator = uniformMutate)
  @assert populationSize > 0
  @assert iterations > 0
  
  #progress bar stuff
  p = Progress(iterations, 1, "Generating solutions", 50)
  
  
  #main loop of the NSGA-II algorithm
  
  #create hall of fame to save the best individuals
  HallOfFame = hallOfFame()
  
  #initialize with two randomly initialized populations
  kickstartingPopulation = initializePopulation(alleles, fitnessFunction, populationSize)
  previousPopulation = initializePopulation(alleles, fitnessFunction, populationSize)
  
  #merge two initial parents
  mergedPopulation = population(vcat(kickstartingPopulation.individuals, previousPopulation.individuals))
  
  
  #|selection -> offsprings| -> |selection -> offsprings| -> ...
  for i = 1:iterations
    #sort the merged population into non dominated fronts
    fronts = nonDominatedSort(mergedPopulation)
    
    #add the best individuals to the hall of fame
    addToHallOfFame(mergedPopulation, fronts[1], HallOfFame)
    
    
    if length(fronts) == 1 || length(fronts[1]) >= populationSize
        crowdingDistance(mergedPopulation, fronts[1], 1, true)
        selectedFromLastFront = lastFrontSelection(mergedPopulation,
                                                   fronts[1],
                                                   populationSize)
        #put the indices of the individuals in all 
        #fronts that were selected as parents
        parentsIndices = selectedFromLastFront
    else
        #separate last front from rest, it is treated differently with
        #lastFrontSelection function
        indexOfLastFront = length(fronts)
        lastFront = fronts[indexOfLastFront]
        fronts = fronts[1: (indexOfLastFront - 1)]
        
        #calculate the crowding distances for all but the last front and 
        #update the mergedPopulation.distance
        for j = 1:length(fronts)
            front_j = fronts[j]
            crowdingDistance(mergedPopulation, front_j, j, true)
        end
        
        #calculate how many individuals are left to 
        #select (there's n-k in the previous fronts)
        k = populationSize - length(reduce(vcat, fronts))

        
        #find the indices of the k individuals we need from the last front
        selectedFromLastFront = lastFrontSelection(mergedPopulation, 
                                                lastFront,
                                                k)
        
        #update the crowding distance on the last front
        crowdingDistance(mergedPopulation, selectedFromLastFront, indexOfLastFront, true)
        
        #put the indices of the individuals in all 
        #fronts that were selected as parents
        parentsIndices = vcat(reduce(vcat, fronts), selectedFromLastFront)
    end

    #--------------------------------------------------------------------
    
    #at this point, we have all we need to create the next population:
    #   -indices of n individuals that were selected from the merged population
    #   -crowding distance and front information (in actualPop.distance)
    
    parentPopulation = population(mergedPopulation.individuals[parentsIndices],
                                  mergedPopulation.distances)
    
    
    #print the dict of distances
#     for dist in keys(parentPopulation.distances)
#       print("[")
#       print(dist)
#       print("]")
#       print(" : ")
#       print(parentPopulation.distances[dist])
#       println("")
#     end
    #we make a tournament selection to select children
    #the templates are actual parents
    childrenTemplates = UFTournSelection(parentPopulation)
    
    #we apply genetic operators (recomination and mutation) to obtain next pop
    nextPopulation = generateOffsprings(childrenTemplates,
                                        probabilityOfCrossover,
                                        probabilityOfMutation,
                                        fitnessFunction,
                                        alleles,
                                        mutationOperator,
                                        crossoverOperator)
    #we now have a new population, we must now
    # -create a new merged population
    # -assign this newly produce population as being the previous of the next loop
    mergedPopulation = population(vcat(nextPopulation.individuals, previousPopulation.individuals))
    previousPopulation = nextPopulation
    
    next!(p)
  end
  
  return [HallOfFame, previousPopulation]
end
#END



#END   NSGA-II main methods
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN helper methods

#BEGIN SAMPLE
function SAMPLE(L::Vector, k::Int)
  #helper local function
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



#BEGIN printResult
function printResult(P::population)
  mi1 = minimum(map(x->x.fitness[1], P.individuals))
  ma1 = maximum(map(x->x.fitness[1], P.individuals))
  
  mi2 = minimum(map(x->x.fitness[2], P.individuals))
  ma2 = maximum(map(x->x.fitness[2], P.individuals))
  
end
#END



#END   helper methods
#------------------------------------------------------------------------------



end


