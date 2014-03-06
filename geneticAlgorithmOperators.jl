#genetic operators
module geneticAlgorithmOperators



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



#END   crossover operators
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#BEGIN mutation operators



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



#END   mutation operators
#------------------------------------------------------------------------------



end