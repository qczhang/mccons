#module NSGA_II
#-------------------------DEFINITION-------------------------
#Simple implementation of the NSGA-II multiobjective
#genetic algorithm.

#-------------------------imports-------------------------
using Base.Test


#-------------------------exported methods-------------------------
#export 


#-------------------------type definitions-------------------------
type arrangement
  #represents combination of structures
  parts::Vector
  fitness::Vector
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



function evaluate(population::Vector{arrangement}, index::Int, compare_method = nonDominatedCompare)
  #compare the object at index with the rest of the vector
  #output must be sorted
  count = 0
  dominatedby = (Int)[]
  indFit = population[index].fitness
  #before the index
  if index!= 1
    for i = 1: (index-1)
      if compare_method(population[i].fitness, indFit) == 1
        count += 1
        push!(dominatedby, i)
      end
    end
  end
  #after the index
  if index != length(population)
    for i = (index+1):length(population) #exclude the index
      if compare_method(population[i].fitness, indFit) == 1
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
  



function nonDominatedSort(population::Vector{arrangement})
  #multi-objective optimization using evolutionnary algorithms p.43
  #step 1 of non dominated sorting
  #sort the population in non dominated fronts
  
  #could also require population to be pair...
  cutoff = div(length(population), 2) + length(population)%2
  
  values = {} #to improve...
  len = length(population)
  
  #1- evaluate the whole population
  for i = 1:len
    push!(values, evaluate(population, i, nonDominatedCompare))
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

function findEqualVectors(a::(Int, Vector))
  #triangular numbers, O(n^2)
  #a :: index, vector)
  #use shift!
  equal = (Int, Int)[]
  toVisit = trues(length(a))
  for i = 1:(length(a)-1)
    if toVisit[i] == true
      elem = a[i][2]
      for j = (i+1):(length(a)-1)
	if toVisit[j] == 1.0
	  if a[j][2] == elem
	    push!(equal, (a[i][1], a[j][1]))
	    toVisit[j] = false
	  end
	end
      end
    end
  end  
  return result
end
  
  

function findIdentical(a::Array{(Int,Array{Int,1}),1})
  #a :: (index, int vector)
  #the vector must be sorted by first index
  equivalents = (Int, Int)[]
  
  
  i = 1
  while i < length(a)
    equivClass = Vector{Int, T}[]
    val = a[i][2][1]
    #get all the elements in the equivalence class
    while a[i][2][1] == val
      push!(equivClass, a[i])
      i+=1
    end
    #if the class is not empty
    if !(length(equivClass) == 1)
      
      
    


function crowdingDistance(population::Vector{arrangement})
  #calculate the crowding distance for the 2N population
  #this is based on "Revisiting the NSGA-II crowding distance computation"
  #step 1
  values = (Int, Vector)[]
  for i = 1:length(population)
    push!(values,(i,population[i].fitness))
  end
  
  #step 2
  #sort by first fitness value
  sort!(values, by=x->x[2][1], 2)
  
  #step 3
  #find uniques and concatenate
  #go in ascending order, verify if equal, add (first, equalx) until equalx is different
  nonUniques = (Int, Int)[]
  
  #step 4
  
  



function lastFrontSelection(input::Vector{arrangement}, n::Int, maxObjs::Vector, minObjs::Vector)
  #multi-objective optimization using evolutionnary algorithms p.248
  #n is the number of individuals to select from this population
  
  #step 1 
  l = length(input)
  m = length(input[1].fitness)
  lastfront = (Int, arrangement)[]
  for i=1:l
    push!(lastfront,(i, input[i]))
  end
  
  values = zeros(l)

    
  #step 2
  sorts = {}
  for i = 1:m
    push!(sorts, map(x->x[1],sort(lastfront, by = x->x[2].fitness[i], rev = true)))
  end
  
  #step 3
  #for each objective 
  for i = 1:m
    ar = sorts[i]
    #first and last have infinite value
    divider =  maxObjs[i] - minObjs[i]
    
    values[ar[1]] = Inf
    values[ar[end]] = Inf
    #rest have value according to multi-objective optimization using evolutionnary algorithms p.250
    for j = 2:l-1
      values[ar[j]] += (input[j+1].fitness[i]  - input[j-1].fitness[i])/divider
    end 
  end
  
  #step 4
  #reassign indices, sort and output the n first indices...
  
  result = (Int, Float64)[]
  for i=1:l
    push!(result, (i, values[i]))
  end
  
  result = sort(result, by = x->x[2], rev=true)[1:n]
  return map(x->x[1], result)
end
  
  


#--module end
#end