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



end