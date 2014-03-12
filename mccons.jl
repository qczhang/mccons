#MCCONS rework



#get NSGA-II and RNA2D 
require("NSGA_II.jl")
require("RNA_2D.jl")



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
  return [evalBPSetDistance(v), evalHausdorff(v)]
end

function evalDistanceEqualLength(v::Vector)
  #uses the mountain distance in addition to base pair set and hausdorff 
  return [evalBPSetDistance(v), evalHausdorff(v), evalMountainDistance(v)]
end

#END



#BEGIN data
yeastTRNAs = {
"tRNA-ASN" =>
"GACUCCAUGGCCAAGUUGGUUAAGGCGUGCGACUGUUAAUCGCAAGAUCGUGAGUUCAACCCUCACUGGGGUCGCCA",

"tRNA-GLY" =>
"GCGCAAGUGGUUUAGUGGUAAAAUCCAACGUUGCCAUCGUUGGGCCCCGGUUCGAUUCCGGGCUUGCGCACCA",

"tRNA-ILE" =>
"GGUCUCUUGGCCCAGUUGGUUAAGGCACCGUGCUAAUAACGCGGGGAUCAGCGGUUCGAUCCCGCUAGAGACCACCA",

"tRNA-LYS" =>
"UCCUUGUUAGCUCAGUUGGUAGAGCGUUCGGCUUUUAACCGAAAUGUCAGGGGUUCGAGCCCCCUAUGAGGAGCCA",

"tRNA-MET" =>
"GCUUCAGUAGCUCAGUAGGAAGAGCGUCAGUCUCAUAAUCUGAAGGUCGAGAGUUCGAACCUCUCCUGGAGCACCA",

"tRNA-THR" =>
"GCUUCUAUGGCCAAGUUGGUAAGGCGCCACACUAGUAAUGUGGAGAUCAUCGGUUCAAAUCCGAUUGGAAGCACCA",

"tRNA-TRP" =>
"GAAGCGGUGGCUCAAUGGUAGAGCUUUCGACUCCAAAUCGAAGGGUUGCAGGUUCAAUUCCUGUCCGUUUCACCA",

"tRNA-ALA" =>
"GGGCGUGUGGCGUAGUCGGUAGCGCGCUCCCUUAGCAUGGGAGAGGUCUCCGGUUCGAUUCCGGACUCGUCCACCA",

"tRNA-ARG" =>
"UUCCUCGUGGCCCAAUGGUCACGGCGUCUGGCUACGAACCAGAAGAUUCCAGGUUCAAGUCCUGGCGGGGAAGCCA",

"tRNA-ASP" =>
"UCCGUGAUAGUUUAAUGGUCAGAAUGGGCGCUUGUCGCGUGCCAGAUCGGGGUUCAAUUCCCCGUCGCGGAGCCA",

"tRNA-GLU" =>
"UCCGAUAUAGUGUAACGGCUAUCACAUCACGCUUUCACCGUGGAGACCGGGGUUCGACUCCCCGUAUCGGAGCCA",

"tRNA-HIS" =>
"GGCCAUCUUAGUAUAGUGGUUAGUACACAACAUUGUGGCUGUUGAAACCCUGGUUCGAUUCUAGGAGGUGGCACCA",

"tRNA-PHE" =>
"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",

"tRNA-VAL" =>
"GGUUUCGUGGUCUAGUCGGUUAUGGCAUCUGCUUAACACGCAGAACGUCCCCAGUUCGAUCCUGGGCGAAAUCACCA"}

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



function foldYeasttRNAs(n::Int)
  vals = Vector{RNA_2D.structure}[]
  k = collect(keys(yeastTRNAs))
  
  for i = 1:length(k)
    folded = callFlashFold(yeastTRNAs[k[i]], n)
    folds = map(x -> RNA_2D.structure(i, x[1], x[2]), folded)
    push!(vals, folds)
  end
  return vals
end


ALLELES = foldYeasttRNAs(50)
populationSize = 50

kickstartingPopulation = NSGA_II.initializePopulation(ALLELES, evalDistance, populationSize)
previousPopulation = NSGA_II.initializePopulation(ALLELES, evalDistance, populationSize)

mergedPopulation = NSGA_II.population(vcat(kickstartingPopulation.individuals, previousPopulation.individuals))
fronts = NSGA_II.nonDominatedSort(mergedPopulation)


indexOfLastFront = length(fronts)
lastFront = fronts[indexOfLastFront]
fronts = fronts[1: (indexOfLastFront - 1)]

for j = 1:length(fronts)
  front_j = fronts[j]
  NSGA_II.crowdingDistance(mergedPopulation, front_j, j, true)
end

k = populationSize - length(reduce(vcat, fronts))



selectedFromLastFront = NSGA_II.lastFrontSelection(mergedPopulation, 
                                               lastFront,
                                               k)
                                               
NSGA_II.crowdingDistance(mergedPopulation, selectedFromLastFront, indexOfLastFront, true)

parentsIndices = vcat(reduce(vcat, fronts), selectedFromLastFront)

parentPopulation = NSGA_II.population(mergedPopulation.individuals[parentsIndices],
                                  mergedPopulation.distances)

childrenTemplates = NSGA_II.UFTournSelection(parentPopulation)
