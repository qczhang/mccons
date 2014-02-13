module RNA_2D
#-------------------------DEFINITION-------------------------
#RNA 2D module is used to compare secondary structures
#It contains type representation, methods to convert representation,
#compare different structures (some restricted to structures of same length).

#-------------------------imports-------------------------
using Base.Test


#-------------------------exported methods-------------------------
export 



#-------------------------type definition-------------------------
immutable structure
  family::Int #needed in mccons
  dotBracket::ASCIIString
  mountain::Vector{Int}
  base_pair_set::Vector{(Int,Int)}
  
  function structure(family::Int, dotBracketInput::ASCIIString)
    @test testDotBracket(dotBracketInput) 
    mountain = dotBracketToMountain(dotBracketInput)
    base_pair_set = dotBracketToBPSet(dotBracketInput)
    self = new(family, dotBracketInput, mountain, base_pair_set)
  end
end


#-------------------------verification methods-------------------------
function testDotBracket(dotBracket::ASCIIString)
  #verifies Vienna dot-bracket for "()" and unbalanced structure
  #could be done with regex... to investigate
  counter = 0
  lastchar = '('
  for i in dotBracket
    #add to structure
    if(i=='(')
      counter+=1
    elseif(i==')')
      counter-=1
    elseif(i!='.')
      return false
    end
    #catch illegal structure
    if(counter < 0)
      return false
    end
    #catch illegal pair
    if(lastchar == '(' && i == ')')
      return false
    end
    lastchar=i
  end
  #catch unbalanced structure
  if(counter!= 0)
    return false
  end
  return true
end


#-------------------------transformation methods-------------------------
function dotBracketToMountain(dotBracket::ASCIIString)
  #transforms Vienna dotbracket to mountain representation
  #e.g. "(.)" -> [0,1,1,0]
  counter = 0
  lastchar = '('
  val::Vector{Int} = [0]
  for i in dotBracket
    #add to structure
    if(i=='(')
      counter+=1
    elseif(i==')')
      counter-=1
    end
    append!(val, [counter])
    lastchar=i
  end
  return val
end

function dotBracketToBPSet(dotBracket::ASCIIString)
  #transforms Vienna dotbracket to base pair set (sorted list by first base of the pair)
  # "((..))" -> [(1,6), (2,5)]
  bpset= (Int,Int)[]
  accumulator = Int[]
  count = 0
  for i in dotBracket
    if(i == '(')
      push!(accumulator, count)
    elseif(i ==')')
      push!(bpset, (pop!(accumulator), count))
    end
    count += 1
  end
  return sort(bpset)
end




#-------------------------comparison methods-------------------------
function compareMountainDistance(m1::Vector{Int}, m2::Vector{Int})
  #lp1 mountain distance on two mountains representation of same length
  #e.g. [1,2,2,2,1], [1,2,3,2,1] = 1
  @assert length(m1) == length(m2)
  absdiff(x::(Int,Int))= abs(x[1]-x[2])
  return mapreduce(absdiff, +, zip(m1, m2))
end


function fastCompareBPSet(bp1::Vector{(Int,Int)}, bp2::Vector{(Int,Int)})
  #THIS ONE STILL NEEDS DEBUGGING AND FURTHER TESTING
  #naive base pair distance (cardinality of symmetric difference, |(A\B)U(B\A)|)
  #
  i = 0
  id1 = 1
  id2 = 1
  result = 0
  while i < min(bp1[end][1], bp2[end][1])
    res = Int[]
    #add the 2nd of tuple if first is equal to index
    #println(bp1[id1])
    if(bp1[id1][1] == i)
      push!(res, bp1[id1][2])
      id1 += 1
    end
    #println(bp2[id2])
    if(bp2[id2][1] == i)
      push!(res, bp2[id2][2])
      id2 += 1
    end
    #println(res)
    #if the tuple has 1, add 1, if it has 2 and 2 are diff, add 2
    if(length(res) == 1)
      result += 1
    elseif(length(res) == 2)
      if(res[1] != res[2])
      result += 2
      end
    end
    #println("cumul = $result")
    i+=1
  end
  #end processing
  if(!(bp1[id1] == bp2[id2]))
    result += 2
  end
  if(id1 == length(bp1))
    result += length(bp2) - id2
    return result
  end
  if(id2 == length(bp2))
    result += length(bp1) - id1
    return result
  end
end




#-------------------------unit test methods-------------------------
function test_testDotBracket()
  #unit test
  @test testDotBracket("((((.)))") == false #missing brackets on the right
  @test testDotBracket("(((.))))") == false #missing brackets on the left
  @test testDotBracket("(((())))") == false
end


function slowCompareBPSet(bp1::Vector{(Int,Int)}, bp2::Vector{(Int,Int)})
  #O(n^2) used for debugging
  bp12 =(Int, Int)[]
  bp21 =(Int, Int)[]
  for i in bp1
    if !(i in bp2)
    push!(bp12, i)
    end
  end
  for i in bp2
    if !(i in bp1)
    push!(bp21, i)
    end
  end
  return length(bp12) + length(bp21)
end


function test_compareBPSet()
  #unit test
  a = dotBracketToBPSet("(((((....)))))")
  b = dotBracketToBPSet("(((....)(((....))(..))))")
  @test slowCompareBPSet(a,b) == 12 == fastCompareBPSet(a,b)

  c = a
  @test slowCompareBPSet(a,c) == 0 == fastCompareBPSet(a,c)
  
  d = slowCompareBPSet("(((....)(((....))(.)(..))))")
  return true
end



#-------------------------export all the methods!-------------------------

export compareMountainDistance, fastCompareBPSet



#--module end
end