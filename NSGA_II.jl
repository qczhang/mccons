module NSGA_II

#-------------------------imports-------------------------
using Base.Test


#-------------------------exported methods-------------------------

export 


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



function fastDelete!(values::Vector, deletion::Vector)
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
  
  values = (Int, Int, Int[])  
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
    push!(result, frontIndices)
    
    #exclude the dominators from the values (could reuse)
    values = filter(x->x[2] != 0, values)
    
    #decrement the count based on the latest front
    for i in values
      #delete the last front from the indices
      substracted = fastDelete(i[2], frontIndices)
      #substract the difference of cardinality
      diff = length(i[2]) - length(substracted)
      i[3] = substracted
      i[2] = diff
    end
  end
  return result
end

function crowdingSort(lastfront::Vector{arrangement})
  #multi-objective optimization using evolutionnary algorithms p.248
  #input is the last front that is over the population limit
  
  #typeof(Inf) == Float64...
  
  #step 1 
  l = length(lastfront)
  m = length(lastfront[1].fitness)
  
  #step 2
  sorts = {}
  for i = 1:m
    push!(sort(lastfront, by = x->x.fitness[i], rev = true))
  end
  
  #step 3
  fkakakakakakkakkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk
  

end







#-------------------------unit test methods-------------------------
function test_nonDominatedCompare()
  #unit test
  a = [1,2,3,4]
  b = [1,2,3,4]
  c = [1,3,3,4]
  d = [0,2,3,4]
  e = [1,2,3,5]
  f = [1,2,3,3]
  #equal
  @test nonDominatedCompare(a,b) == 0
  @test nonDominatedCompare(b,a) == 0
  
  #smaller
  @test nonDominatedCompare(a,c) == -1
  @test nonDominatedCompare(d,a) == -1
  @test nonDominatedCompare(f,a) == -1
  @test nonDominatedCompare(a,e) == -1
  
  #bigger
  @test nonDominatedCompare(c,a) == 1
  @test nonDominatedCompare(a,d) == 1
  @test nonDominatedCompare(a,f) == 1
  @test nonDominatedCompare(e,a) == 1
  
  return true
end

function test_evaluate()
  #unit test
  a = arrangement([], [1,2,3,4,5])
  b = arrangement([], [0,0,4,2,0])
  c = arrangement([], [0,0,1,1,1])
  d = arrangement([], [0,2,2,1,5])
  pop = [a,b,c,d]
  @test evaluate(pop, 1, nonDominatedCompare) == (1,0,[])
  @test evaluate(pop, 2, nonDominatedCompare) == (2,0,[])
  @test evaluate(pop, 3, nonDominatedCompare) == (3,2,[1,4])
  @test evaluate(pop, 4, nonDominatedCompare) == (4,1,[1])
  return true
end


function slowDelete(values::Vector, deletion::Vector)
  return filter(x->!(x in deletion), values)
end 

function generatePosRandInt(n::Int)
  return filter(x->x>0, unique(sort(rand(Int16, n))))
end

function test_fastDelete(repet::Int, size::Int)
  #unit test
  #as far as can tell, equivalent to slowDelete (warning, the rand(Int, n) is too sparse to test )
  #use rand(Int16,n) or something like that
  for i= 1:repet
    values = generatePosRandInt(size)
    deletion = generatePosRandInt(size)
    @test slowDelete(values, deletion) == fastDelete(values, deletion)
  end
  return true
end


function test_all()
  #unit test all
  test_evaluate()
  test_fastDelete(1000,1000)
  test_nonDominatedCompare()
  return true
end



#--module end
end