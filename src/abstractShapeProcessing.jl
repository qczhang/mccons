#preprocessing to relate abstract shapes within predicted secondary structures

#try 1:
#use the longest common subsequence to find similar abstract shapes
# function longestCS{T}(s1::T, s2::T)
#   #only allows deletions and insertions
#   #http://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Computing_the_length_of_the_LCS
#   m = length(s1)
#   n = length(s2)
#   C = zeros(m+1, n+1)
#   #the array is already filled with zeros
#   for i = 2:(m+1)
#     for j=2:(n+1)
#       if s1[i-1] == s2[j-1]
#         C[i,j] = C[i-1, j-1] + 1
#       else
#         C[i,j] = max(C[i, j-1], C[i-1, j])
#       end
#     end
#   end
#   #println(C)
#   return C[m+1, n+1]
# end

# oneBracketDistance(s1, s2) = (max(length(s1), length(s2)) - longestCS(s1, s2)) <= 2
#this approach is erroneous because of the "[[]]" and "[][]" which work according to the function
#but are not at a distance of 1 insertion :(

#try 2: use deletion
function bracketToPairs{T<:String}(s::T)
  #does the same as finding the base pair set
  #but on a level 5 abstract shape
  stack = Int[]
  result = (Int,Int)[]
  for i = 1:length(s)
    if s[i] == '['
      push!(stack, i)
    elseif s[i] == ']'
      push!(result, (pop!(stack), i))
    else
      error("wrong format for the abstract shape")
    end
  end  
  return sort(result)
end



function pairsToBrackets(pairs::Vector{(Int,Int)})
  #base pair -> brackets
  indexBrackets = (Int, Char)[]
  for i = 1:length(pairs)
    push!(indexBrackets, (pairs[i][1], '['))
    push!(indexBrackets, (pairs[i][2], ']'))
  end
  return CharString(map(x->x[2], sort(indexBrackets, by = x->x[1])))
end



function getAtOneDeletion{T<:String}(s::T)
  #returns all abstract shapes at a distance of one bracket deletion
  pairs = bracketToPairs(s)
  #
  newSets = Vector{(Int,Int)}[]
  #
  for i = 1:length(pairs)
    toAdd= (Int,Int)[]
    for j = 1:length(pairs)
      if j!= i
        push!(toAdd, pairs[j])
      end
    end
    push!(newSets, toAdd)
  end
  #
  result = String[]
  for i in newSets
    bracket = pairsToBrackets(i)
    if !(bracket in result) && bracket != ""
      push!(result, bracket)
    end
  end
  #
  return result
end


