

#BEGIN readme
#unit tests for the RNA_2D module
#END readme



#BEGIN imports
require("RNA_2D.jl")
using Base.Test
#END



#BEGIN unit tests
function test_testDotBracket()
  #unit test
  @test RNA_2D.testDotBracket("((((.)))") == false #missing brackets on the right
  @test RNA_2D.testDotBracket("(((.))))") == false #missing brackets on the left
end


function test_init_structure()
  a = RNA_2D.structure(1, "(.)(((((....)))))")
  b = RNA_2D.structure(2, "(.)(.)(...........)", -42.0)
  return true
end
  
  
  
function slowCompareBPSet(bp1::Vector{(Int,Int)}, bp2::Vector{(Int,Int)})
  #helper method
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



function test_compareBPSet(n::Int)
  #exhaustive unit test
  @assert n > 0
  for i = 1:n
    a = RNA_2D.dotBracketToBPSet(randomDotBracketPlus())
    b = RNA_2D.dotBracketToBPSet(randomDotBracketPlus())
    @test slowCompareBPSet(a,b) == RNA_2D.compareBPSet(a,b)
    @test RNA_2D.compareBPSet(a,b) == RNA_2D.compareBPSet(b,a)
  end
  return true
end



function randomDotBracket()
  #helper method
  #generate random valid dot bracket
  #choice of three moves
  # 1- (
  # 2- )
  # 3- .
  opening = ['(', '.', '.']
  closing = [')', '(', '.','.']
  function addSymbol(choices, values, stack)
    sym = choices[rand(1:length(choices))]
    if sym == '('
      push!(values, sym)
      return stack + 1
      
    elseif sym == ')' && stack != 0
      push!(values, sym)
      return stack -1
      
    elseif sym == '.'
      push!(values, sym)
    
    else return -1
    end
    
    return stack
  end
  
  partialAddSymbol = s->addSymbol(s, val, stack)
  
  stack = 0
  val = Char[]
  stack = partialAddSymbol(opening)
  
  while true && length(val) < 40
    #println(stack)
    #println(CharString(val))
    if val[end] == '('
      stack = partialAddSymbol(opening)
      
    else #symbol is either ')' or '.'
      stack = partialAddSymbol(closing)
      
      if stack == -1
	return CharString(val)
      end
    end
  end
    
  #just to avoid the problem of ()
  push!(val, '.')
  while stack > 0
    push!(val, ')')
    stack -= 1
  end
  return CharString(val)
end



function randomDotBracketPlus()
  #generate non empty dotbrackets (ie. not "." or ".."...)
  x = randomDotBracket()
  while true
    for i in x
      if i !='.'
	return x
      end
    end
    x = randomDotBracket()
  end
end



function test_randomDotBracket(n::Int)
  #unit test
  @assert n > 0
  for i = 1:n
    @test RNA_2D.testDotBracket(randomDotBracket())==true
  end
  return true
end



function test_RNAshapes(n::Int)
  #use only on computers where RNAshapes is installed
  @assert n > 0
  for i = 1:n
    dotB = randomDotBracket()
    t5 = chomp(readall(`RNAshapes -D $dotB -t5`))
    t3 = chomp(readall(`RNAshapes -D $dotB -t3`))
    t1 = chomp(readall(`RNAshapes -D $dotB -t1`))
    result = RNA_2D.RNAshapes(dotB)
    @test t5 == result[1]
    @test t3 == result[2]
    @test t1 == result[3]
  end
  return true
end



function test_compareHausdorff(n::Int)
  @assert n>0
  #simple debug example
  S1 = "........((((...))))."
  S2 = ".......((((...)))).."
  
  B1 = RNA_2D.dotBracketToBPSet(S1)
  B2 = RNA_2D.dotBracketToBPSet(S2)
  
  @test RNA_2D.compareHausdorff(B1,B2) == 1
  
  for i = 1:n
    a = RNA_2D.dotBracketToBPSet(randomDotBracketPlus())
    b = RNA_2D.dotBracketToBPSet(randomDotBracketPlus())
    RNA_2D.compareHausdorff(a,b)
  end
end


function test_all()
  test_testDotBracket()
  test_randomDotBracket(10000)
  test_compareBPSet(1000)
  test_compareHausdorff(100)
  test_init_structure()
  println("All unit tests succeeded")
  return true
end
#END



test_all()




