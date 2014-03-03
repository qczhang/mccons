

#BEGIN readme
#Unit tests for the NSGA_II module
#END



#BEGIN imports
require("NSGA_II")
using Base.Test
#END



#BEGIN unit tests

#BEGIN test_solution
function test_solution()
  #tests constuctors for solution
  s1 = solution([2], [4])
  s2 = solution([2], x->map(y->2*y, x))
  @test s1 == s2
end
#END

#BEGIN test_population
function test_population()
  #tests constructors for population
  
end

#END


#BEGIN test_nonDominatedCompare
function test_nonDominatedCompare(n::Int, fitnessSize::Int)
  #test by property 
  tests = {}
  
  #create vector of random fitness arrays
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
      #dominating
      if v == 1
	@test all_compare(i,j, >=) == true
      #dominated
      elseif v == -1
	@test all_compare(i,j, <=) == true
      #non dominated and non dominating
      elseif v == 0
	@test all_compare(i,j, >) == false
	@test all_compare(i,j, <) == false
      end
    end
  end
  return true
end
#END test_nonDominatedCompare



#BEGIN randomFitnessArray
function randomFitnessArray(fitnessLen::Int)
  #helper
  return rand(1:10000, fitnessLen)
end
#END randomFitnessArray



#BEGIN test_nonDominatedSort
function test_nonDominatedSort(cardinality::Int, fitnessLen::Int)
  #unit test
  #exhaustive
  pop = population()
  for i =  1: cardinality
    push!(pop.solutions, solution([42],randomFitnessArray(fitnessLen)))
  end

  sorts = nonDominatedSort(pop)
  #no domination within the same front
  for i = 1:length(sorts)
    ar = sorts[i]
    for j in ar
      for k in ar
	@test nonDominatedCompare(pop.solutions[j].fitness, pop.solutions[k].fitness) == 0
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
	  @test nonDominatedCompare(pop.solutions[j].fitness, pop.solutions[k].fitness) in (0,1)
	end
      end
    end
  end
    
  return true
end
#END test_nonDominatedSort



#BEGIN test_evaluateAgainstOthers
function test_evaluateAgainstOthers(cardinality::Int, fitnessLen::Int, compare_method = nonDominatedCompare)
  #exhaustive unit test
  #generate the population
  pop = population()
  for i =  1: cardinality
    push!(pop.solutions, solution([42],randomFitnessArray(fitnessLen)))
  end
  #evaluate all solutions
  result = {}
  for i = 1:cardinality
    push!(result, evaluateAgainstOthers(pop, i, compare_method))
  end
  #verify validity
  for i = 1:cardinality
    if !(isempty(result[i][3]))
      for j in result[i][3]
	@test compare_method(pop.solutions[result[i][1]].fitness, pop.solutions[j].fitness) == -1
      end
    end
  end
  
  return true
end
#END test_evaluateAgainstOthers



#BEGIN slowDelete
function slowDelete(values::Vector, deletion::Vector)
  #helper
  return filter(x->!(x in deletion), values)
end 
#END slowDelete



#BEGIN generatePosRandInt
function generatePosRandInt(n::Int, minInt = 1, maxInt = 10000)
  #helper
  @assert minInt < maxInt
  @assert n > 0
  return sort(rand(minInt:maxInt, n))
end
#END generatePosRandInt



#BEGIN test_fastDelete
function test_fastDelete(repet::Int, size::Int)
  #unit test, exhaustive
  for i= 1:repet
    values = generatePosRandInt(size)
    deletion = generatePosRandInt(size)
    @test slowDelete(values, deletion) == fastDelete(values, deletion)
  end
  return true
end
#END test_fastDelete



#BEGIN generatePop
function generatePop(size::Int, fitnessLen::Int)
  p = solution[]
  for i = 1:size
    push!(p, solution({42}, randomFitnessArray(fitnessLen)))
  end
  
  return population(p)
end
#END



#BEGIN test_crowdingDistance
function test_crowdingDistance()
  #create population
  p = solution[]
  push!(p, solution({42}, [0,5]))
  push!(p, solution({42}, [0,5]))
  push!(p, solution({42}, [2,2]))
  push!(p, solution({42}, [3,1]))
  push!(p, solution({42}, [3,1]))
  push!(p, solution({42}, [3,1]))
  push!(p, solution({42}, [5,0]))
  pop = population(p)
  
  #sort into fronts
  fronts = nonDominatedSort(pop)
  
  #calculate crowding distances
  crowdingDistance(pop, fronts[1], 1, true)
  
  #test against manually calculated values
  @test pop.distances[[0,5]] == (1,Inf)
  @test pop.distances[[2,2]] == (1,1.4)
  @test pop.distances[[3,1]] == (1,1.0)
  @test pop.distances[[5,0]] == (1,Inf)
  
end
#END



#BEGIN test_lastFrontSelection


#END



#BEGIN test_all
function test_all()
  #exhaustive
  test_nonDominatedCompare(1000,3)
  test_evaluateAgainstOthers(1000,5)
  test_fastDelete(2000,2000)
  test_nonDominatedSort(2000, 5)
  test_crowdingDistance()
  println("All unit tests succeeded")
  
  return true
end
#END test_all

test_all()


#END


