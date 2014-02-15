#module NSGA_II
#-------------------------DEFINITION-------------------------
#Simple implementation of the NSGA-II multiobjective
#genetic algorithm.

#-------------------------imports-------------------------
using Base.Test


#-------------------------exported methods-------------------------
#export 


#-------------------------type definitions-------------------------
immutable solution
  units::Vector
  fitness::Vector
end

type population
  individuals::Vector{solution}
end

#-------------------------genetic algorithm methods-------------------------
function nonDominatedCompare (a::Vector, b::Vector, comparator = >)
  # a > b --> 0: a==b, 1: a>b, -1: b>a, pairwise vector comparison
  #A solution is called nondominated, Pareto optimal, 
  #Pareto efficient or noninferior, if none of the objective 
  #functions can be improved in value without degrading some 
  #of the other objective values.
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
    
    if AdomB && BdomA
      return 0
    end
  end
 
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



function evaluate(pop::population, index::Int, compare_method = nonDominatedCompare)
  #compare the object at index with the rest of the vector
  #output must be sorted
  count = 0
  dominatedby = (Int)[]
  indFit = pop[index].fitness
  #before the index
  if index!= 1
    for i = 1: (index-1)
      if compare_method(pop[i].fitness, indFit) == 1
        count += 1
        push!(dominatedby, i)
      end
    end
  end
  #after the index
  if index != length(pop)
    for i = (index+1):length(pop) #exclude the index
      if compare_method(pop[i].fitness, indFit) == 1
        count += 1
        push!(dominatedby, i)
      end
    end
  end
  return (index, count, dominatedby)
end



function fastDelete(values::Vector, deletion::Vector)
  #both are sorted, we take that into account
  #deletion vector cannot be empty by definition
  #start at naturals, not whole integer range
  @assert values[1] > 0
  @assert deletion[1] > 0
  result = Int[]
  indexDel = 1
  for i in values
    #iterate to the next "good" index
    #the value of deletion is either > or = to i
    while (deletion[indexDel] < i) && (indexDel < length(deletion))
      indexDel += 1
    end

    if i!=deletion[indexDel]
      push!(result, i)
    end
  end
  return result
end
  

function nonDominatedSort(pop::population) 
  #TO REPLACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #multi-objective optimization using evolutionnary algorithms p.43
  
  
  #step 1 of non dominated sorting
  #sort the pop in non dominated fronts
  
  #could also require population to be pair...
  cutoff = div(length(pop), 2) + length(pop)%2
  
  values = {} #to improve...
  len = length(pop)
  
  #1- evaluate the whole population
  for i = 1:len
    push!(values, evaluate(pop, i, nonDominatedCompare))
  end
  
  fronts = {}
  #2- get the fronts
  while length(values) > cutoff #we stop at half
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

function uniqueFitness(a::Array{(Int,Array{Int,1}),1})
  #map fitness to index (get unique fitness)
  equal = Dict{Vector{Int}, Vector{Int}}()
  for i in a
    equal[i[2]] = push!(get(equal, i[2], Int[]), i[1])
  end
  return equal
end
      
   


function crowdingDistance(pop::population)
  #calculate the crowding distance for the 2N population
  #this is based on "Revisiting the NSGA-II crowding distance computation"
  
  #step 1
  #assing index to keep track after sorting
  values = (Int, Vector)[]
  for i = 1:length(pop)
    push!(values,(i,pop[i].fitness))
  end
  
  #step 2 -----------------------------------------------------------
  #get the unique fitness vectors
  fitnessMapping = uniqueFitness(values)
  fitnesses = keys(fitnessMapping)
  
  #step 3 -----------------------------------------------------------
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


  
  


#--module end
#end