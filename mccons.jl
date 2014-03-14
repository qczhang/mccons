#MCCONS rework



#get NSGA-II and RNA2D 
require("NSGA_II.jl")
require("RNA_2D.jl")
require("geneticAlgorithmOperators")


#BEGIN callFlashFold
function callFlashFold(sequence::String, ft::Int)
  data = split(readall(`./f32 -seq $sequence -ft $ft`), "\n")
  if(data[end] == "")
    data = data[1:end-1]
  end
  data = map(x->split(x), data)
  data2 = (String, FloatingPoint)[]
  for x in data
    push!(data2, ((convert(String, x[1])), float(x[2])))
  end
  return data2
end
#END callFlashFold



function pairwiseAll(f::Function, v::Vector)
  #helper, compares all to all, including self to self,
  #since distances have 0, doesn't make a difference
  result = {}
  for i=1:length(v)
    for j = i:length(v)
      push!(result, f(v[i], v[j]))
    end
  end  
  return result
end



#BEGIN evaluation functions
function evalBPSetDistance(v::Vector)
  v2 = map(x->x.base_pair_set, v)
  reduce(+, pairwiseAll(RNA_2D.compareBPSet, v2))
end



function evalMountainDistance(v::Vector)
  v2 = map(x->x.mountain, v)
  reduce(+, pairwiseAll(RNA_2D.compareMountainDistance, v2))
end



function evalHausdorff(v::Vector)
  v2 = map(x->x.base_pair_set, v)
  reduce(+, pairwiseAll(RNA_2D.compareHausdorff, v2))
end



function evalDistance(v::Vector)
  return [-evalBPSetDistance(v), -evalHausdorff(v)]
end

function evalDistanceEqualLength(v::Vector)
  #uses the mountain distance in addition to base pair set and hausdorff 
  return [evalBPSetDistance(v), evalHausdorff(v), evalMountainDistance(v)]
end

#END



#BEGIN data
IREs = {
#length = 
"AY112742_1_12_41" => "GUcCUGCUUCAACAGUGCUUGAACGGaAC",
"BC019840_1_11_40" => "GUcUUGCUUCAACAGUGUUUGAACGGaAC",
"AF086786_1_2_28" => "UcGUUCGUCCUCAGUGCAGGGCAACaG",
"S57280_1_391_417" => "UcGUUCGUCCUCAGUGCAGGGCAAUaG",
"AF171078_1_1416_1442" => "UgGUUCGUCCUCAGUGCAGGGCAACaG",
"AJ426432_1_1593_1619" => "AUUAUCGGGAGCAGUGUCUUCCAUAAU",
"X13753_1_1371_1397" => "AUUAUCGGGGGCAGUGUCUUCCAUAAU",
"X01060_1_3482_3508" => "AUUAUCGGAAGCAGUGCCUUCCAUAAU",
"BC001188_1_3791_3817" => "AUUAUCGGGAACAGUGUUUCCCAUAAU",
"X13753_1_1481_1507" => "AUUAUCGGGGACAGUGUUUCCCAUAAU",
"AJ426432_1_1658_1684" => "UAUAUCGGAGACAGUGAUCUCCAUAUG" ,
"M58040_1_3309_3335" => "UAUAUCGGAGaCAGUGAcCUCCAUAUG" ,
"X13753_1_1434_1460" => "UAUAUCGGAGGCAGUGACCUCCAUAUG" ,
"X01060_1_3950_3976" => "UGUAUCGGAGACAGUGAUCUCCAUAUG" ,
"X01060_1_3432_3458" => "UUUAUCAGUGACAGAGUUCACUAUAAA" ,
"X13753_1_830_856" => "UUUAUCAGUGACAGCGUUCACUAUAAA" ,
"AY120878_1_50_76" => "GgUCGCGUCAACAGUGUUUGAUCGAaC" ,
#length = 
"AB062402_1_11_40" => "uUUCCUGCUUCAACAGUGCUUGGACGGAAc" ,
"AB073371_1_5_34" => "uCUCCUGCUUCAACAGUGCUUGGACGGAGc" ,
"AF285177_1_3_32" => "GUUCCUGCUUCAACAGUGCUUGGACGGAAC" ,
"M16343_1_1306_1335" => "GUUCCUGCgUCAACAGUGCUUGGaCGGAAC",
"AF338763_1_11_40" => "UUaCCUGCUUCAACAGUGCUUGAACGGcAA",
"AF117958_1_132_161" => "ucUCUUGUUUCAACAGUGUUUGGACGGAac",
"AF266195_1_14_43" => "AUUCUUGCUUCAACAGUGUUUGAACGGAAU" ,
"D86625_1_6_35" => "GUUCUUGUUUCAACAGUGAUUGAACGGAAC" ,
"S77386_1_28_57" => "GUUCUUGCUUCAACAGUGAUUGAACGGAAC" ,
"M12120_1_24_53" => "GUUCUUGCUUCAACAGUGUUUGAACGGAAC" ,
"L39879_1_1190_1219" => "GUaCUUGCUUCAACAGUGUUUGAACGGaAC",
"J02741_1_400_429" => "aUCUUGCUUCAACAGUGUUUGGACGGAa" }
#END



function generateAlleles(RNAs::Vector{String}, n::Int)
  @assert n > 0
  #we use flash fold to generate n suboptimals per RNA
  alleles = Vector{(String, FloatingPoint)}[]
  for i in RNAs
    push!(alleles, callFlashFold(i, n))
  end
  return alleles
end



function foldAllToStructure(dict::Dict, n::Int)
  #get the keys
  k = collect(keys(dict))
  vals = Vector{RNA_2D.structure}[]
  
  for i = 1:length(k)
    folded = callFlashFold(dict[k[i]], n)
    folds = map(x -> RNA_2D.structure(i, x[1], x[2]), folded)
    push!(vals, folds)
  end
  return vals
end


function dotBracketsToShapes(dotB::Vector, lvl::Int)
  #
  @assert lvl in [1,2,3]
  shapeToDotB = Dict()
  for i in dotB
    shape = RNA_2D.RNAshapes(i)[lvl]
    println(shape)
    actualValue = get(shapeToDotB, shape, Any[])
    push!(actualValue, i)
    shapeToDotB[shape] = actualValue
  end
  return shapeToDotB
end

function separateAllelesByAbstractShapes(alleles::Vector, lvl::Int)
  #returns a set of similar alleles (according to abstract shape) 
  #to be tested for similarity
  @assert lvl in [1,2,3]
  #lvl: 1->lvl5, 2->lvl3, 3->lvl1
  differentShapes = Set{String}()
  shapeToStructures = Dict{String, Vector{(Int, String)}}()
  
  #find unique shapes and map structures 
  #to their abstract shapes (lvl5, lvl3, lvl1)
  for moleculeType = 1:length(alleles)
    for struct = 1:length(alleles[moleculeType])
      println(length(moleculeType))
      dotBracket = alleles[moleculeType][struct].dotBracket
      println(dotBracket)
      shape = RNA_2D.RNAshapes(dotBracket)[lvl]
      actualMapping = get(shapeToStructures, shape, (Int, String)[])
      println(actualMapping)
      println(typeof(actualMapping))
      push!(actualMapping, (moleculeType, dotBracket))
      println(actualMapping)
      shapeToStructures[shape] = actualMapping
    end
  end
  
  println(length(differentShapes))
  
  #assemble sets of shapes
  return shapeToStructures
end

function printResultDotB(P::NSGA_II.population, i::Int)
  @assert i in [0:length(P.individuals)]
  solution = P.individuals[i]
  println("fitness \n$(solution.fitness)")
  for i in solution.genes
    println(i.dotBracket)
  end
end


function main(popSize = 50,numIterations=50, alleleSize = 50)
  mutationOperator = geneticAlgorithmOperators.uniformMutate
  crossoverOperator = geneticAlgorithmOperators.uniformCrossover
  
  ALLELES = foldAllToStructure(IREs, alleleSize)
  
  r = NSGA_II.main(ALLELES,
                   evalDistance,
                   popSize,
                   numIterations,
                   0.2,
                   0.15,
                   crossoverOperator,
                   mutationOperator)

end



function writeResult{T<:String}(fileName::T, extension::T, P::NSGA_II.population, alleles::Vector)
  #will not rewrite files in theory
  
  #assert that the file is not already in the local directory
  listDir = readdir()
  name = string(fileName, extension)
  if name in listDir
    error("$name is already present in the local directory, delete it if you want to write")
  end
  
  f = open(name, "w")
  
  #format is:
  #>id
  #fitness
  #dotBrackets
  #end
  ind = P.individuals
  fitnessLen = length(ind[1].fitness)
  
  for i=1:length(ind)
    #write the id
    println(f, ">$i")
    
    fit = ind[i].fitness
    s = fit[1]
    for j = 2:fitnessLen
      s = string(s," ", fit[j])
    end
    #write the fitnesses on the same line
    println(f, s)
    
    #write the dot brackets associated
    for j=1:length(ind[i].genes)
    \todddddddddddddddddddddddooooooooooooooooooooooooooooooo append subopt number
#       dotBracket = 
#       indexAt = findfirst(alleles[j], )
      println(f, j.dotBracket)
    end
    #write end symbol
    println(f, "end")
  end
  
  close(f)
end


function extractFitness{T<:String}(fileName::T)
  fileName = string(fileName)
  if !(fileName in readdir())
    error("$fileName is not in the local directory")
  end
  
  f = open(fileName, "r")
  lines = readlines(f)
  println(lines)
  close(f)
  
  f2 = open(string("fitness_",fileName), "w")
  #println(lines[1])
  i=1
  while i < length(lines)
    #println(i)
    #println(lines[i])
    if beginswith(lines[i], ">")
      #println(lines[i])
      index = chomp(lines[i][2:end])
      fitnesses = chomp(lines[i+1])
      println(f2, string(index," ", fitnesses))
      i+=2
    else
      i+=1
    end
  end
  close(f2)

end


