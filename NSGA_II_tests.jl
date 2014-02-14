#Unit tests for the NSGA_II module

require("NSGA_II")
using Base.Test


function test_nonDominatedCompare(n::Int, fitnessSize::Int)
  #exhaustive unit test
  tests = {}
  for i =1:n
    push!(tests, randomFitnessArray(fitnessSize))
  end
  function all_compare(x,y, op)
    #helper
    for i in zip(x,y)
      if !(op(i[1],i[2]))
	return false
      end
    end
    return true
  end
  for i in tests
    for j in tests
      v = nonDominatedCompare(i,j)
      if v == 1
	@test all_compare(i,j, >=) == true
      elseif v== -1
	@test all_compare(i,j, <=) == true
      end
    end
  end
  return true
end

function randomFitnessArray(fitnessLen::Int)
  #helper
  return map(abs, rand(Int, fitnessLen))
end

function test_nonDominatedSort(cardinality::Int, fitnessLen::Int)
  #unit test
  #exhaustive
  individuals = arrangement[]
  for i =  1: cardinality
    push!(individuals, arrangement([],randomFitnessArray(fitnessLen)))
  end
  sorts = nonDominatedSort(individuals)
  #no domination within the same front
  for i = 1:length(sorts)
    ar = sorts[i]
    for j in ar
      for k in ar
	@test nonDominatedCompare(individuals[j].fitness, individuals[k].fitness) == 0
      end
    end
  end
  #domination or equivalence of all for greater front
  #all in 1 dominate all in 2
  if(length(sorts)>1)
    for i = 1:length(sorts)-1
      a = sorts[i]
      b = sorts[i+1]
      for j in a
	for k in b
	  @test nonDominatedCompare(individuals[j].fitness, individuals[k].fitness) in (0,1)
	end
      end
    end
  end
    
  return true
end
  
  
function test_evaluate(cardinality::Int, fitnessLen::Int, compare_method = nonDominatedCompare)
  #exhaustive unit test
  #generate the population
  population = arrangement[]
  for i =  1: cardinality
    push!(population, arrangement([],randomFitnessArray(fitnessLen)))
  end
  
  #evaluate all individuals
  result = {}
  for i = 1:cardinality
    push!(result, evaluate(population, i, compare_method))
  end
  
  #verify validity
  for i = 1:cardinality
    if !(isempty(result[i][3]))
      for j in result[i][3]
	@test compare_method(population[result[i][1]].fitness, population[j].fitness) == -1
      end
    end
  end
  
  return true
end

function slowDelete(values::Vector, deletion::Vector)
  #helper
  return filter(x->!(x in deletion), values)
end 

function generatePosRandInt(n::Int)
  #helper
  return filter(x->x>0, unique(sort(rand(Int16, n))))
end

function test_fastDelete(repet::Int, size::Int)
  #unit test, exhaustive
  for i= 1:repet
    values = generatePosRandInt(size)
    deletion = generatePosRandInt(size)
    @test slowDelete(values, deletion) == fastDelete(values, deletion)
  end
  return true
end

function test_lastFrontSelection()
  a = arrangement([], [1,2,3,4,5])
  b = arrangement([], [0,0,4,2,0])
  c = arrangement([], [2,0,1,1,1])
  d = arrangement([], [0,2,2,1,5])
  e = arrangement([], [2,1,1,1,1])
  f = arrangement([], [2,2,2,2,2])
  g = arrangement([], [2,2,2,2,2])
  h = arrangement([], [2,3,2,2,2])
  pop = [a,b,c,d,e,f,g,h]
  maxObjs = [5,5,5,5,5]
  minObjs = [0,0,0,0,0]
  
  

end


function test_findEqualVectors(n::Int)
  vals = Array{(Int,Array{Int,1}),1}[]
  for i = 1:n
    push!(vals, (i, randomFitnessArray(5)))
  end
  
  
    
end

function test_all()
  #exhaustive
  test_nonDominatedCompare(1000,3)
  test_nonDominatedSort(2000, 5)
  test_evaluate(1000,5)
  test_fastDelete(2000,2000)
  
  
  
  
  #non exhaustive
  test_lastFrontSelection()
  test_findEqualVectors()
  return true
end


