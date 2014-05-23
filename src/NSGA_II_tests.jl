


#------------------------------------------------------------------------------
#BEGIN readme



#Unit tests for the NSGA_II module



#END   readme
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#BEGIN imports



require("NSGA_II.jl")
using Base.Test



#END   imports
#------------------------------------------------------------------------------




#------------------------------------------------------------------------------
#BEGIN unit tests



#BEGIN test_nonDominatedSort
function test_nonDominatedSort(cardinality::Int, fitnessLen::Int)
  #verify property of non dominated sorting results
  pop = NSGA_II.population()
  for i =  1: cardinality
    push!(pop.individuals, NSGA_II.individual([42],randomFitnessArray(fitnessLen)))
  end

  sorts = NSGA_II.nonDominatedSort(pop)
  @test length(sorts) != 0
  
  #no domination within the same front
  for i = 1:length(sorts)
    ar = sorts[i]
    for j in ar
      for k in ar
        fit1 = pop.individuals[j].fitness
        fit2 = pop.individuals[k].fitness
        @test NSGA_II.nonDominatedCompare(fit1, fit2) == 0
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
    fit1 = pop.individuals[j].fitness
    fit2 = pop.individuals[k].fitness
    @test NSGA_II.nonDominatedCompare(fit1, fit2) in (0,1)
  end
      end
    end
  end
    
  return true
end
#END



#BEGIN test_individual
function test_individual()
  #tests constuctors for individual
  s1 = NSGA_II.individual([2], [4])
  s2 = NSGA_II.individual([2], x->map(y->2*y, x))
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
      v = NSGA_II.nonDominatedCompare(i,j)
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
#END



#BEGIN randomFitnessArray
function randomFitnessArray(fitnessLen::Int)
  #helper
  return rand(1:10000, fitnessLen)
end
#END



#BEGIN test_evaluateAgainstOthers
function test_evaluateAgainstOthers(cardinality::Int, 
                                    fitnessLen::Int, 
                                    compare_method = NSGA_II.nonDominatedCompare)
  #test evaluation of individual against everyone else in the population (helper)
  
  #generate the population
  pop = NSGA_II.population()
  for i =  1: cardinality
    push!(pop.individuals, NSGA_II.individual([42],randomFitnessArray(fitnessLen)))
  end
  #evaluate all individuals
  result = {}
  for i = 1:cardinality
    push!(result, NSGA_II.evaluateAgainstOthers(pop, i, compare_method))
  end
  #verify validity
  for i = 1:cardinality
    if !(isempty(result[i][3]))
      for j in result[i][3]
        a = pop.individuals[result[i][1]].fitness
        b = pop.individuals[j].fitness
        @test compare_method(a, b) == -1
      end
    end
  end
  
  return true
end
#END



#BEGIN generatePosRandInt
function generatePosRandInt(n::Int, minInt = 1, maxInt = 10000)
  #helper
  @assert minInt < maxInt
  @assert n > 0
  return sort(rand(minInt:maxInt, n))
end
#END



#BEGIN test_fastDelete
function test_fastDelete(repet::Int, size::Int)

  function slowDelete(values::Vector, deletion::Vector)
    #helper, used to compare with fastDelete
    return filter(x->!(x in deletion), values)
  end
  
  #unit test, exhaustive
  for i= 1:repet
    values = generatePosRandInt(size)
    deletion = generatePosRandInt(size)
    @test slowDelete(values, deletion) == NSGA_II.fastDelete(values, deletion)
  end
  return true
end
#END



#BEGIN generatePop
function generatePop(size::Int, fitnessLen::Int)
  p = NSGA_II.individual[]
  for i = 1:size
    push!(p, NSGA_II.individual({42}, randomFitnessArray(fitnessLen)))
  end
  
  return NSGA_II.population(p)
end
#END



#BEGIN test_crowdingDistance
function test_crowdingDistance()
  #create population
  p = NSGA_II.individual[]
  push!(p, NSGA_II.individual({42}, [0,5]))
  push!(p, NSGA_II.individual({42}, [0,5]))
  push!(p, NSGA_II.individual({42}, [2,2]))
  push!(p, NSGA_II.individual({42}, [3,1]))
  push!(p, NSGA_II.individual({42}, [3,1]))
  push!(p, NSGA_II.individual({42}, [3,1]))
  push!(p, NSGA_II.individual({42}, [5,0]))
  pop = NSGA_II.population(p)

  #sort into fronts
  fronts = NSGA_II.nonDominatedSort(pop)

  #calculate crowding distances
  NSGA_II.crowdingDistance(pop, fronts[1], 1, true)
  
  #test against manually calculated values
  @test pop.distances[[0,5]] == (1,Inf)
  @test pop.distances[[2,2]] == (1,1.4)
  @test pop.distances[[3,1]] == (1,1.0)
  @test pop.distances[[5,0]] == (1,Inf)

end
#END



#BEGIN test_lastFrontSelection


#END

#BEGIN test_main
function test_main(n::Int)
  #uses the 0-1 sum to check it does indeed optimize
  allele = [0,1]
  ALLELES = Vector{Int}[]
  for i=1:50
    push!(ALLELES, allele)
  end
  
  function f(x)
    #we maximize the sum of the genes
    v = 0
    for i in x
      v+=i[1]
    end
    return v
  end
  
  function g(x)
    #we minimize the sum of the genes
    v= 0 
    for i in x
      v-=i[1]
    end
    return v
  end
  
  evalF(x) = [f(x), g(x)]
  
  mutationOperator = NSGA_II.uniformMutate
  crossoverOperator = NSGA_II.uniformCrossover
  x =  NSGA_II.main(ALLELES,
                    evalF,
                    100,
                    n,
                    0.1,
                    0.05,
                    crossoverOperator,
                    mutationOperator)
  


end



#BEGIN test_all
function test_all()
  #exhaustive
  test_nonDominatedCompare(500,3)
  test_evaluateAgainstOthers(500,5)
  test_fastDelete(1000,1000)
  test_nonDominatedSort(500, 3)
  test_crowdingDistance()
  println("All unit tests succeeded")
  
  return true
end
#END test_all

test_all()


#END   unit tests
#------------------------------------------------------------------------------

