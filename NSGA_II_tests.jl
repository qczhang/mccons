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
  return map(abs, rand(Int16, fitnessLen))
end

function test_nonDominatedSort(cardinality::Int, fitnessLen::Int)
  #unit test
  #exhaustive
  pop = population(Vector{solution}[],Dict{Vector, FloatingPoint}())
  for i =  1: cardinality
    push!(pop.individuals, solution([],randomFitnessArray(fitnessLen)))
  end
  

  sorts = nonDominatedSort(pop)
  #no domination within the same front
  for i = 1:length(sorts)
    ar = sorts[i]
    for j in ar
      for k in ar
	@test nonDominatedCompare(pop.individuals[j].fitness, pop.individuals[k].fitness) == 0
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
	  @test nonDominatedCompare(pop.individuals[j].fitness, pop.individuals[k].fitness) in (0,1)
	end
      end
    end
  end
    
  return true
end
  
  
function test_evaluateAgainstOthers(cardinality::Int, fitnessLen::Int, compare_method = nonDominatedCompare)
  #exhaustive unit test
  #generate the population
  pop = population(Vector{solution}[],Dict{Vector, FloatingPoint}())
  for i =  1: cardinality
    push!(pop.individuals, solution([],randomFitnessArray(fitnessLen)))
  end
  #evaluate all individuals
  result = {}
  for i = 1:cardinality
    push!(result, evaluateAgainstOthers(pop, i, compare_method))
  end
  #verify validity
  for i = 1:cardinality
    if !(isempty(result[i][3]))
      for j in result[i][3]
	@test compare_method(pop.individuals[result[i][1]].fitness, pop.individuals[j].fitness) == -1
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




function test_uniqueFitness(n::Int)
  #exhaustive unit test
  #generate data
  vals = [(1, randomFitnessArray(5))]
  for i = 2:n
    push!(vals, (i, randomFitnessArray(5)))
  end
  #generate doubles
  i = 0
  while i < n/2
    push!(vals, (n+i+1, vals[rand(1:n)][2]))
    i+=1
  end
  #go through the dict and test for equality
  d = uniqueFitness(vals)
  for i in keys(d)
    if(length(d[i])>1)
      v = vals[d[i][1]][2]
      for j in d[i]
	@test vals[j][2] == v
      end
    end
  end
  return true  
end




  
  
function test_all()
  #exhaustive
  test_nonDominatedCompare(1000,3)
  test_evaluateAgainstOthers(1000,5)
  test_fastDelete(2000,2000)
  test_nonDominatedSort(2000, 5)
  test_uniqueFitness(1000)
  
  
  
  #non exhaustive
  return true
end


