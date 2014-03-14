#dumb cons is dumb
require("mccons")


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
  ALLELES = foldYeasttRNAs(n)
  alleles = Vector{String}[]
  #fetch the dot brackets from the structure type
  for i in ALLELES
    toAdd = String[]
    for j in i
      push!(toAdd, j.dotBracket)
    end
    push!(alleles, toAdd)
  end
  #
  
  return dumbCons(alleles)
end

