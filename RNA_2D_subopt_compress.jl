#reduce the number of suboptimal based on the notion of probabilistic shape
#but not on a whole space of folding, simply the input space

#we assume that the input will be a list of suboptimal structures with
#associated free energy (kCal / mol)

push!(LOAD_PATH, chomp(readall(`pwd`)))
import RNA_2D



function boltzmann(energy::FloatingPoint, T = 295)
  #temperature is in kelvin
  R = BigFloat(0.00198717)
  eu^(-energy / (R * T))
end



function sumFreeEnergy(suboptList::Vector{(String, FloatingPoint)})
  reduce(+, map(boltzmann, map(x->x[2], suboptList)))
end



function probabilisticPartition(suboptList::Vector{(String, FloatingPoint)})
  #sum the free energy
  Qs = sumFreeEnergy(suboptList)
  #create the dicts {shape => free energy}
  l5 = Dict{String, FloatingPoint}()
  l3 = Dict{String, FloatingPoint}()
  l1 = Dict{String, FloatingPoint}()
  structureToShape = Dict{String, Vector{String}}()
  #iterate over suboptimals
  for i in suboptList
    shapes = RNA_2D.RNAshapes(i[1])
    #add to the mapping {2D structure => shapes}
    structureToShape[i[1]] = shapes
    #add to the energy sum
    l5[shapes[1]] = get(l5, shapes[1], 0.0) + i[2]
    l3[shapes[2]] = get(l3, shapes[2], 0.0) + i[2]
    l1[shapes[3]] = get(l1, shapes[3], 0.0) + i[2]
  end
  return {l5, l3, l1, structureToShape}
end



function RMSD()

end



function testCompressionPerformance(initialsuboptList::Vector{(String, FloatingPoint)}, finalsuboptList::Vector{(String, FloatingPoint)})
  #estimator is the RMSD
  p1 = probabilisticPartition(initialsuboptList)
  p2 = probabilisticPartition(finalsuboptList)
  l1_keys = keys(p1[1])
  l3_keys = keys(p1[2])
  l5_keys = keys(p1[3])
  l1_bias = 0.0
  l3_bias = 0.0
  l5_bias = 0.0
  for i in l5_keys
    l5_bias += ((l5[i] - get(p2[1], i, 0.0))^2) 
  end
  l5_bias /= length(l5_keys)
  for i in l3_keys
    l3_bias += ((l3[i] - get(p2[2], i, 0.0))^2)
  end
  l3_bias /= length(l3_keys)
  for i in l1_keys
    l1_bias += ((l1[i] - get(p2[3], i, 0.0))^2)
  end
  l1_bias /= length(l1_keys)
  return {l5_bias, l3_bias, l1_bias}
end



function stochasticUniversalSampling(suboptList::Vector{(String, FloatingPoint)}, n::Int)
  #assign boltzmann to structure
  values = collect(zip(map(x->x[1], suboptList), map(x->boltzmann(x[2]), suboptList)))
  #sum the energy
  Qs = reduce(+, map(x->x[2], values))
  #sort by reverse order of energy
  sort!(values, by = x->x[2], rev = true)
  #acquire the step size
  step = Qs / n
  #choose a small value within the first elem (is it really how?)
  target = rand()* values[1][2]
  index = 1
  accumulatedValue = 0.0
  #add the chosen ones to the result
  chosenOnes = (String, FloatingPoint, Int)[]
  for i=1:n #revise the boundary conditions to be sure
    target += step
    while (accumulatedValue + value[index][2])  < target
      accumulatedValue+= value[index][2]
      index = ((index+1) %n) # we wrap around
    end
    push!(chosenOnes, values[index])
  end
  
  return suboptList[map(x->x[3], chosenOnes)]
end
  
  

function probabilisticChoose(n::Int, suboptList::Vector{(String, FloatingPoint)})
  @assert n > 0
  #from "Reducing bias and inefficiency in the selection algorithm"
  #http://en.wikipedia.org/wiki/Stochastic_universal_sampling
  
  #choose n elements from the partition
  #there is a big risk of losing much of the fidelity if there is too much
  #granularity (n too small)
  partition = probabilisticPartition(suboptList)
  shapeToStructure
  #use some kind of localized universal sampling
  probab = probabilisticPartition(suboptList)
  
  #do a 3 level sampling
  #lvl5
  for i = 1:n
    \todo: select per class
    
  end    
end


