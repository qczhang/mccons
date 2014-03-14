#dumb cons is dumb



function callFlashFold(sequence::String, ft::Int)
  data = split(readall(`./f32 -seq $sequence -ft $ft`), "\n")
  if(data[end] == "")
    data = data[1:end-1]
  end
  data = map(x->split(x), data)
  data2 = String[]
  for x in data
    push!(data2, convert(String, x[1]))
  end
  return data2
end



function foldAll(sequences::Dict, ft::Int)
  result = Vector{String}[]
  for i in keys(sequences)
    seq = sequences[i]
    push!(result, callFlashFold(seq, ft))
  end
  return result
end


#BEGIN data
# EMBL Number Sequence / Structure Rank
# Consensus ((((((.((((((...).)))))))))))
# AY112742_1_12_41 GUcCUGCUUCAACAGUGCUUGAACGGaAC (1st)
# BC019840_1_11_40 GUcUUGCUUCAACAGUGUUUGAACGGaAC (1st)
# AF086786_1_2_28 -UcGUUCGUCCUCAGUGCAGGGCAACaG- (1st)
# S57280_1_391_417 -UcGUUCGUCCUCAGUGCAGGGCAAUaG- (5th)
# AF171078_1_1416_1442 -UgGUUCGUCCUCAGUGCAGGGCAACaG- (1st)
# AJ426432_1_1593_1619 -AUUAUCGGGAGCAGUGUCUUCCAUAAU- (1st)
# X13753_1_1371_1397 -AUUAUCGGGGGCAGUGUCUUCCAUAAU- 
# X01060_1_3482_3508 -AUUAUCGGAAGCAGUGCCUUCCAUAAU- (1st)
# BC001188_1_3791_3817 -AUUAUCGGGAACAGUGUUUCCCAUAAU- 
# X13753_1_1481_1507 -AUUAUCGGGGACAGUGUUUCCCAUAAU- 
# AJ426432_1_1658_1684 -UAUAUCGGAGACAGUGAUCUCCAUAUG- (1st)
# M58040_1_3309_3335 -UAUAUCGGAGaCAGUGAcCUCCAUAUG- (1st)
# X13753_1_1434_1460 -UAUAUCGGAGGCAGUGACCUCCAUAUG- (1st)
# X01060_1_3950_3976 -UGUAUCGGAGACAGUGAUCUCCAUAUG- (1st)
# X01060_1_3432_3458 -UUUAUCAGUGACAGAGUUCACUAUAAA- (1st)
# X13753_1_830_856 -UUUAUCAGUGACAGCGUUCACUAUAAA- (1st)
# AY120878_1_50_76 -GgUCGCGUCAACAGUGUUUGAUCGAaC- (1st)
# Consensus (((((.(.((((((...).)))))))))))
# AB062402_1_11_40 uUUCCUGCUUCAACAGUGCUUGGACGGAAc 
# AB073371_1_5_34 uCUCCUGCUUCAACAGUGCUUGGACGGAGc 
# AF285177_1_3_32 GUUCCUGCUUCAACAGUGCUUGGACGGAAC 
# M16343_1_1306_1335 GUUCCUGCgUCAACAGUGCUUGGaCGGAAC 
# AF338763_1_11_40 UUaCCUGCUUCAACAGUGCUUGAACGGcAA 
# AF117958_1_132_161 ucUCUUGUUUCAACAGUGUUUGGACGGAac (21th)
# BC016354_1_30_59 ucUCUUGCUUCAACAGUGUUUGGACGGAac 
# AF266195_1_14_43 AUUCUUGCUUCAACAGUGUUUGAACGGAAU 
# D86625_1_6_35 GUUCUUGUUUCAACAGUGAUUGAACGGAAC 
# S77386_1_28_57 GUUCUUGCUUCAACAGUGAUUGAACGGAAC 
# M12120_1_24_53 GUUCUUGCUUCAACAGUGUUUGAACGGAAC 
# L39879_1_1190_1219 GUaCUUGCUUCAACAGUGUUUGAACGGaAC 
# J02741_1_400_429 -aUCUUGCUUCAACAGUGUUUGGACGGAa- 



# EMBL Number Sequence / Structure Rank
# Consensus ((((((.((((((...).)))))))))))
IREs = {
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


function clusterByLength(alleles::Vector)
  #the suboptimals within each nested vectors have the same length
  lengthToCategory = Dict{Int, Vector{Int}}()
  for i=1:length(alleles)
    l = length(alleles[i][1])
    actualValue = get(lengthToCategory,l, Int[])
    push!(actualValue, i)
    lengthToCategory[l] = actualValue
  end
  return lengthToCategory
end



function setify{T}(vector::Vector{T})
  #put all the elements of a vector into a set
  result = Set{T}()
  for i in vector
    push!(result, i)
  end
  return result
end



function intersectAtIndices(alleles::Vector, indices::Vector{Int})
  i1 = indices[1]
  result = setify(alleles[i1])
  for i = 2:(length(indices))
    index = indices[i]
    result = intersect(result, setify(alleles[index]))
  end
  return result
end



function dumbCons{T<:String}(alleles::Vector{Vector{T}})
  #cluster them by length
  byLength = clusterByLength(alleles)
  
  #no gap free alignment if they all have different length
  if length(byLength) == length(alleles)
    return []
  end
  
  #convert to sets to do intersections
  allelesSets = Set[]
  for i=1:length(alleles)
    set = Set{T}()
    for j in alleles[i]
      push!(set, j)
    end
    push!(allelesSets, set)
  end
  
  intersectDict = Dict{Int, Set{T}}()
  #do the intersections
  for i in keys(byLength)
    #we get the vector of indices
    indices = byLength[i]
    #
    intersectDict[i] = intersectAtIndices(alleles, indices)
  end
  
  return intersectDict
end


function main(n::Int)
  @assert n>0
  dotBrackets = foldAll(IREs, n)
  
  dumbCons(dotBrackets)
end

