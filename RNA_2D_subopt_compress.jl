#reduce the number of suboptimal based on the notion of probabilistic shape
#but not on a whole space of folding, simply the input space

#we assume that the input will be a list of suboptimal structures with
#associated free energy (kCal / mol)

#BEGIN setup
push!(LOAD_PATH, chomp(readall(`pwd`)))
import RNA_2D
#END setup



#BEGIN boltzmann
function boltzmann(energy::FloatingPoint, T = 295)
  #temperature is in kelvin
  R = BigFloat(0.00198717)
  eu^(-energy / (R * T))
end
#END boltzmann



#BEGIN callFlashFold
function callFlashFold(sequence::String, ft::Int)
  data = split(readall(`./f32 -seq $sequence -ft $ft`), "\n")
  if(data[end] == "")
    data = data[1:end-1]
  end
  data = map(x->split(x), data)
  data = map(x->(x[1], float(x[2])), data)
  return data
end
#END callFlashFold



#BEGIN test_comparePartitions
function test_comparePartitions(proportion::FloatingPoint)
  @assert 0<= proportion <= 1
  a = getObject("test_compression/test_compress_tRNA-GLY_10k")
  b = suboptCompress(a, proportion)
  
  #get all the info from the partitions
  partition1 = abstractShape_energy_Partition(a)
  partition2 = abstractShape_energy_Partition(b)
  
  #look at level5 first (arbitrary)
  p1 = partition1[1]
  p2 = partition2[1]
  
  compared = comparePartitions(p1,p2)
  l5_1 = partition1[4]
  l5_2 = partition2[4]
  
  cardinalities = Dict{String, (Int, Int)}()
  for i in keys(compared)
    cardinalities[i] = (length(get(l5_1, i, [])), length(get(l5_2, i, [])))
  end
  
  result = {}
  for i in keys(compared)
    push!(result, (compared[i], cardinalities[i][1], cardinalities[i][2]))
  end
  sort!(result, by=x->x[1], rev=true)
  
  return result
end
#END test_comparePartitions



#BEGIN comparePartitions
function comparePartitions(partition1, partition2)
  #look at level5 first (arbitrary)
  #the input is directly from abstractShape_energy_Partition
  @assert length(p1)==length(p2)==6
  
  #get the energies per abstract shape
  p1 = partition1[1]
  p2 = partition2[1]
  
  #get the shape=>suboptimals to evaluate the number of structures per abstract shape
  l5_1 = partition1[4]
  l5_2 = partition2[4]
  
  #
  cardinalities = Dict{String, (Int, Int)}()
  for i in keys(compared)
    cardinalities[i] = (length(get(l5_1, i, [])), length(get(l5_2, i, [])))
  end
  
  result = {}
  for i in keys(compared)
    push!(result, (compared[i], cardinalities[i][1], cardinalities[i][2]))
  end
  sort!(result, by=x->x[1], rev=true)
  #
  Qs1 = BigFloat(0)
  Qs2 = BigFloat(0)
  
  #get the Qs for both of them
  for i in keys(p1)
    Qs1 += p1[i] 
  end
  
  for i in keys(p2)
    Qs2 += p2[i]
  end
  
  #divide the values
  for i in keys(p1)
    p1[i] /= Qs1
  end
  for i in keys(p2)
    p2[i] /= Qs2
  end
  
  #assign the absolute difference of probability
  #do not assume that they have the all the same keys
  result = Dict{String, BigFloat}()
  for i in keys(p1)
    result[i] = p1[i]
  end
  
  for i in keys(p2)
    result[i] = abs(get(result, i, BigFloat(0)) - p2[i])
  end
  return result
end
#END comparePartitions



#BEGIN fixedPointFlashFold
function fixedPointFlashFold(structure::String,initialFT::Int, maxFT::Int, incrementSize::Int, epsilon::FloatingPoint)
  #run flashfold until there is a fixed point reached in the output
  #using a small increment will automatically stop it (will create so little change)
  #will only run if f32 in the path...
  dataBefore = abstractShape_energy_Partition(callFlashFold(structure, initialFT))
  while initialFT < maxFT
    #fetch the next suboptimals and calculate their partition
    initialFT += incrementSize
    dataAfter = abstractShape_energy_Partition(callFlashFold(structure, initialFT))
    #compare the previous and the actual
   
  end
end
#END fixedPointFlashFold



#BEGIN abstractShape_energy_Partition
function abstractShape_energy_Partition(suboptList::Vector)
  #{shape => suboptimals}
  l5structures = Dict{String, Vector}()
  l3structures = Dict{String, Vector}()
  l1structures = Dict{String, Vector}()
  
  #{shape => boltzmann probability}
  l5energy = Dict{String, BigFloat}()
  l3energy = Dict{String, BigFloat}()
  l1energy = Dict{String, BigFloat}()
  
  for i in suboptList
    #acquire the abstract shapes
    shapes = RNA_2D.RNAshapes(i[1])

    #add to the energy sum
    l5energy[shapes[1]] = get(l5energy, shapes[1], BigFloat(0.0)) + boltzmann(i[2])
    l5structures[shapes[1]] = push!(get(l5structures, shapes[1], {}), i)
    
    l3energy[shapes[2]] = get(l3energy, shapes[2], BigFloat(0.0)) + boltzmann(i[2])
    l3structures[shapes[2]] = push!(get(l3structures, shapes[2], {}), i)
    
    l1energy[shapes[3]] = get(l1energy, shapes[3], BigFloat(0.0)) + boltzmann(i[2])
    l1structures[shapes[3]] = push!(get(l1structures, shapes[3], {}), i) 
  end

  return {l5energy, l3energy, l1energy, l5structures, l3structures, l1structures}
end
#END abstractShape_energy_Partition



#BEGIN rouletteWheelSel
function rouletteWheelSel(l::Vector, proportion::FloatingPoint)
  @assert 0 <= proportion <= 1
  #its granular, not much to do about it
  n = ceil(length(l)*proportion)
  
  #assign boltzmann probability
  l2 = map(x->(x[1], boltzmann(x[2])), l)
  
  #get probability sum
  Qs = reduce(+, map(x->x[2], l2))
  
  #initialize
  result = (String, FloatingPoint)[]
  i = 0 
  while i < n
    #target is between 0 and Qs
    r = rand()
    target = r * Qs
    @assert 0<=target<=Qs
#     println("rand = $r")
#     println("target = $target")
    index = 1
    cumul = BigFloat(0)
    
    while (cumul + (l2[index][2])) < target
      cumul += l2[index][2]
      index+=1
    end
    
    push!(result, l[index])
    i+=1
  end
  
  return result
end
#END rouletteWheelSel



#BEGIN stochasticUniversalSampling
function stochasticUniversalSampling(suboptListOriginal::Vector, proportion::FloatingPoint)
  #assert the proportion makes sense
  @assert 0<=proportion<=1
  
  #calculate the number of sequences to extract
  n = ceil(proportion*length(suboptListOriginal))
  
  #convert floating point to BigFloat to avoid overflow and add the index
  suboptList = (String, BigFloat, Int)[]
  for i = 1:length(suboptListOriginal)
    push!(suboptList, (suboptListOriginal[i][1], boltzmann(suboptListOriginal[i][2]), i))
  end
  
  #sum the energy
  Qs = reduce(+, map(x->x[2], suboptList))
  
  #sort by reverse order of Boltzmann energy (state energy)
  #sort!(suboptList, by = x->x[2], rev = true)
  
  #acquire the step size
  step = BigFloat(Qs / n)
  
  #choose a small value within the first elem (is it really how?)
  target = rand()* suboptList[1][2]
  index = 1
  accumulatedValue = BigFloat(0.0)
  
  #add the chosen ones to the result
  chosenOnes = (String, FloatingPoint)[]
  
  #choose like a comb
  #from "Reducing bias and inefficiency in the selection algorithm"
  #http://en.wikipedia.org/wiki/Stochastic_universal_sampling
#   println("len suboptlist = $(length(suboptList))")
#   println("stepsize = $step")
  for i=1:n
#     println("i = $i")
#     println("index = $index")
#     println("accumulated val = $accumulatedValue")
#     println("target = $target")
#     
    target += step
    while (accumulatedValue + suboptList[index][2])  < target
      accumulatedValue += suboptList[index][2]
      index += 1
      if index == (length(suboptList)+1)
	index = 1
      end
#       println("index = $index")
    end
    push!(chosenOnes, suboptListOriginal[suboptList[index][3]])
  end
  
  return chosenOnes
end
#END stochasticUniversalSampling  
  


#BEGIN suboptCompress
function suboptCompress(suboptList::Vector, proportion::FloatingPoint)
  #there is a big risk of losing much of the fidelity if there is too much
  #granularity (the proportion is too small in relation to the original)
  @assert 0 <= proportion <= 1
 
  #calculate the abstract shape clustering and cumulative boltzmann
  partition = abstractShape_energy_Partition(suboptList)
  
  #shape => boltzmann probability
  shapeToBoltzmann = partition[1:3]
  
  #shape => suboptimals
  shapeToSubopts = partition[4:6]
  
  #selection vectors for different levels of shape
  selection_lvl_5 = (String, FloatingPoint)[]
  selection_lvl_3 = Vector[]
  selection_lvl_1 = Vector[]
  
  #level 5 sampling
  l5Keys = collect(keys(shapeToSubopts[1]))
  for k in l5Keys
    selection_lvl_5 = vcat(selection_lvl_5, stochasticUniversalSampling(shapeToSubopts[1][k], proportion))
  end
  
  #only level5 is implemented yet
  return selection_lvl_5
end
#END suboptCompress



#BEGIN tests



#BEGIN getObject
function getObject(name::String)
  #open
  f = open(name, "r")
  #fetch
  data = deserialize(f)
  #close
  close(f)
  #correct the last element
  data = data[1:end-1]
  #split the strings
  data = map(x->split(x, " "), data)
  #fetch only the dotbracket and convert the energy value to float
  data = map(x->(x[1], float(x[2])), data)
  return data
end
#END getObject



#BEGIN test1
function test1()
  #first test on tRNA-ASN
  data = getObject("/tests_compression/test_compress_tRNA-ASN")
end
#END test1



#BEGIN test2
function test2()
  data = getObject("test_compression/test_compress_tRNA-GLY_10k")
  
end
#END test2
 
 
  
#BEGIN testCompressionPerformance
function testCompressionPerformance(initialsuboptList::Vector, proportion::FloatingPoint)

  #compress the suboptimals
  finalsuboptList = suboptCompress(initialsuboptList, proportion)
  
  #calculate their respective shape / energy partitions
  p1 = abstractShape_energy_Partition(initialsuboptList)
  p2 = abstractShape_energy_Partition(finalsuboptList)
  
#   println("keys of initial = $(collect(keys(p1[1])))")
#   println("keys of final = $(collect(keys(p2[1])))")
  #initialize the result dict
  result = Dict{String, (BigFloat, BigFloat)}()
  
  #calculate the sums of boltzmann energy 
  s1 = BigFloat(0)
  for i in keys(p1[1])
    s1 += p1[1][i]
  end
  
  s2 = BigFloat(0)
  for i in keys(p2[1])
    s2 += p2[1][i]
  end
  
  #calculate the relative proportion
  for i in keys(p1[1])
    result[i] = ((p1[1][i])/s1, (p2[1][i])/s2)
  end
  #\todo debug la somme semble pas etre a 100............
  
  #get the dicts to arrays
  ar = (String, BigFloat, BigFloat)[]
  for i in keys(result)
    push!(ar, (i, result[i][1], result[i][2]))
  end
  sort!(ar, by = x->x[2], rev = true)
  return ar
end
  
#   
#   for i in l5_keys
#     l5_bias += ((l5[i] - get(p2[1], i, 0.0))^2) 
#   end
#   l5_bias /= length(l5_keys)
#   for i in l3_keys
#     l3_bias += ((l3[i] - get(p2[2], i, 0.0))^2)
#   end
#   l3_bias /= length(l3_keys)
#   for i in l1_keys
#     l1_bias += ((l1[i] - get(p2[3], i, 0.0))^2)
#   end
#   l1_bias /= length(l1_keys)
#   return {l5_bias, l3_bias, l1_bias}
# end
#END testCompressionPerformance



#END tests
