#preprocessing to relate abstract shapes within predicted secondary structures
require("RNA_2D.jl")

#use the longest common subsequence to find sequences 
function longestCS{T}(s1::T, s2::T)
  #only allows deletions and insertions
  #http://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Computing_the_length_of_the_LCS
  m = length(s1)
  n = length(s2)
  C = zeros(m+1, n+1)
  #the array is already filled with zeros
  for i = 2:(m+1)
    for j=2:(n+1)
      if s1[i-1] == s2[j-1]
        C[i,j] = C[i-1, j-1] + 1
      else
        C[i,j] = max(C[i, j-1], C[i-1, j])
      end
    end
  end
  #println(C)
  return C[m+1, n+1]
end


oneBracketDistance(s1, s2) = (max(length(s1), length(s2)) - longestCS(s1, s2)) <= 2