#filter folding output by shape

require("RNA_2D.jl")

tRNAs = {
"ASN" =>
"GACUCCAUGGCCAAGUUGGUUAAGGCGUGCGACUGUUAAUCGCAAGAUCGUGAGUUCAACCCUCACUGGGGUCGCCA",
"GLY" =>
"GCGCAAGUGGUUUAGUGGUAAAAUCCAACGUUGCCAUCGUUGGGCCCCGGUUCGAUUCCGGGCUUGCGCACCA",
"ILE" =>
"GGUCUCUUGGCCCAGUUGGUUAAGGCACCGUGCUAAUAACGCGGGGAUCAGCGGUUCGAUCCCGCUAGAGACCACCA",
"LYS" =>
"UCCUUGUUAGCUCAGUUGGUAGAGCGUUCGGCUUUUAACCGAAAUGUCAGGGGUUCGAGCCCCCUAUGAGGAGCCA",
"MET" =>
"GCUUCAGUAGCUCAGUAGGAAGAGCGUCAGUCUCAUAAUCUGAAGGUCGAGAGUUCGAACCUCUCCUGGAGCACCA",
"THR" =>
"GCUUCUAUGGCCAAGUUGGUAAGGCGCCACACUAGUAAUGUGGAGAUCAUCGGUUCAAAUCCGAUUGGAAGCACCA",
"TRP" =>
"GAAGCGGUGGCUCAAUGGUAGAGCUUUCGACUCCAAAUCGAAGGGUUGCAGGUUCAAUUCCUGUCCGUUUCACCA",
"ALA" =>
"GGGCGUGUGGCGUAGUCGGUAGCGCGCUCCCUUAGCAUGGGAGAGGUCUCCGGUUCGAUUCCGGACUCGUCCACCA",
"ARG" =>
"UUCCUCGUGGCCCAAUGGUCACGGCGUCUGGCUACGAACCAGAAGAUUCCAGGUUCAAGUCCUGGCGGGGAAGCCA",
"ASP" =>
"UCCGUGAUAGUUUAAUGGUCAGAAUGGGCGCUUGUCGCGUGCCAGAUCGGGGUUCAAUUCCCCGUCGCGGAGCCA",
"GLU" =>
"UCCGAUAUAGUGUAACGGCUAUCACAUCACGCUUUCACCGUGGAGACCGGGGUUCGACUCCCCGUAUCGGAGCCA",
"HIS" =>
"GGCCAUCUUAGUAUAGUGGUUAGUACACAACAUUGUGGCUGUUGAAACCCUGGUUCGAUUCUAGGAGGUGGCACCA",
"PHE" =>
"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA",
"VAL" =>
"GGUUUCGUGGUCUAGUCGGUUAUGGCAUCUGCUUAACACGCAGAACGUCCCCAGUUCGAUCCUGGGCGAAAUCACCA"}



tRNAconsensus = {
"ASN" =>
"(((((((..((((.........))))((((((((...))))))))....((((((...)..))))))))))))....",
"GLY" =>
"(((((((..(((((.....)))))((((((((...))))))))..((((((...)..))))))))))))....",
"ILE" =>
"(((((((..(((((.......)))))((((((((...))))))))....((((((...)..))))))))))))....",
"LYS" =>
"(((((((..(((((......)))))((((((((...))))))))....((((((...)..))))))))))))....",
"MET" =>
"(((((((..(((((......)))))((((((((...))))))))....((((((...)..))))))))))))....",
"THR" =>
"(((((((..((((........))))((((((((...))))))))....((((((...)..))))))))))))....",
"TRP" =>
"(((((((..(((((.....)))))((((((((...))))))))....((((((...)..))))))))))))....",
"ALA" =>
"(((((((..(((((......)))))(((((((.....)))))))....((((((...)..))))))))))))....",
"ARG" =>
"(((((((..(((((......)))))((((((((...))))))))....((((((...)..))))))))))))....",
"ASP1" =>
"(((((((..(((((......)))))((((((.......))))))...((((((...)..))))))))))))....",
"ASP2" =>
"(((((((..(((((......)))))(((((((.....)))))))...((((((...)..))))))))))))....",
"GLU" =>
"(((((((..(((((......)))))((((((((...))))))))...((((((...)..))))))))))))....",
"HIS" =>
"((((((((..(((((......)))))(((((((.....)))))))...((((((...)..)))))))))))).)..",
"PHE1" =>
"(((((((..((((........))))(((((((.....)))))))....((((((...)..))))))))))))....",
"PHE2" =>
"(((((((..(((((......)))))((((((.(....)))))))....((((((...)..))))))))))))....",
"VAL" =>
"(((((((..(((((.......)))))((((((((...))))))))....((((((...)..))))))))))))...."}


rRNA5S = {
"Ecoli1" =>
"UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU",
"Ecoli2" =>
"UGCCUGGCGGCAGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCUCAUGCGAGAGUAGGGAACUGCCAGGCAU",
"Ecoli3" =>
"UGCCUGGCGGCAGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU",
"Ecoli4" =>
"UGUCUGGCGGCAGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGACAU",
"Ecoli5" =>
"UGCCUGGCGGCCUUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU",
"Ecoli6" =>
"UGUCUGGCGGCAGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGACUCCCCAUGCGAGAGUAGGGAACUGCCAGACAU",
"Ecoli10" =>
"UGUCUGGCGGCAGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCUCAUGCGAGAGUAGGGAACUGCCAUGCAU",
"Ecoli11" =>
"UGCCUGGCGGCAGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAUCA",
"Ecoli14" =>
"UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACGCCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGACAU",
}

rRNA5Sconsensus = {
"Ecoli1" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli2" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli3" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli4" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli5" =>
"(((((((((((.(.((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli6" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli10" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).)))))))))))).",
"Ecoli11" =>
"((((((((((.(.(((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...))))))))))))))))).))))))))))))...",
"Ecoli14" =>
"((((((((((.(..((((((((((((((((((((((((.....))))))))))).)))))))))).)))(((((((((((((((((...)))))))))))))))))..))))))))))).",
}


function callFlashFold(sequence::String, ft::Int)
  #calls flash fold
  data = split(readall(`./f32 -seq $sequence -ft $ft`), "\n")
  if(data[end] == "")
    data = data[1:end-1]
  end
  #
  data = map(x->split(x), data)
  data2 = (String, FloatingPoint)[]
  for x in data
    push!(data2, ((convert(String, x[1])), float(x[2])))
  end
  #
  unique(data2)
end


function callFlashFoldMask(sequence::String, mask::String, ft::Int)
  #calls flash fold
  data = split(readall(`./f32 -seq $sequence -ft $ft -m $mask`), "\n")
  if(data[end] == "")
    data = data[1:end-1]
  end
  #
  data = map(x->split(x), data)
  data2 = (String, FloatingPoint)[]
  for x in data
    push!(data2, ((convert(String, x[1])), float(x[2])))
  end
  #
  unique(data2)
end


function computeShapes(V::Vector)
  result = {}
  for shape in V
    push!(result, RNA_2D.RNAshapes(shape))
  end
  result
end

function convertToFFMask(mask::String)
  #
  mask = replace(mask, ".","u")
  mask = replace(mask, "x", ".")
  mask = replace(mask, "u", "x")
end


function calculateCumulative(distribution::Vector,lvl::Int)
  #calculates the cumulative distribution of shapes for the
  #desired level
  #             1 = lvl 1 (less abstract)
  #             2 = lvl 3
  #             3 = lvl 5 (most abstract)
  @assert lvl in [1,2,3]
  shapes = map(x->x[lvl], distribution)
  shapeDict = Dict{String, Int}()
  for shape in shapes
    shapeDict[shape] = get(shapeDict, shape, 0) + 1
  end
  shapeDict
end


function calculateL1(V::Vector)
  a = computeShapes(map(x->x[1],V))
  b = calculateCumulative(a, 1)
end


function calculateAllCumulativeL1(yeastTRNAs, ft::Int)
  @assert ft >= 0
  k = keys(yeastTRNAs)
  result = {}
  cumul = Dict{String, Int}()
  for n in k
    D = calculateL1(callFlashFold(yeastTRNAs[n], ft))
    L = (String, Int)[]
    for (k,v) in D
      cumul[k] = get(cumul, k, 0) + v
      push!(L, (k,v))
    end
    sort!(L, by=x->x[2], rev=true)
    println(n)
    println(L)
    push!(result, (n, L))
  end
  (result, cumul)
end


function filterByShape(sequence::String, foldSize::Int, numberWantedBrackets::Int, shape::String, shapeLvl::Int, cutoff = 100000)
  @assert shapeLvl in [1,2,3]

  function numFitting(dotBrackets, shape, shapeLvl)
    reduce(+, map(x->x==shape, map(x->RNA_2D.RNAshapes(x[1])[shapeLvl], dotBrackets)))
  end
  
  function fitting(dotBrackets, shape, shapeLvl)
    result = (String, FloatingPoint)[]
    shapes = map(x->RNA_2D.RNAshapes(x[1])[shapeLvl], dotBrackets)
    filter(x->RNA_2D.RNAshapes(x[1])[shapeLvl] == shape, dotBrackets)
  end
  
  dotBrackets = callFlashFold(sequence, foldSize)
  numFittingBrackets = numFitting(dotBrackets, shape, shapeLvl)
  while numFittingBrackets < numberWantedBrackets && foldSize < cutoff
    foldSize = foldSize * 2
    dotBrackets = callFlashFold(sequence, foldSize)
    numFittingBrackets = numFitting(dotBrackets, shape, shapeLvl)
    if numFittingBrackets > numberWantedBrackets
      break
    end
    println(numFittingBrackets, " at ", foldSize)
  end
  
  fitting(dotBrackets, shape, shapeLvl)
end