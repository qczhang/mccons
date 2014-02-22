module RNA_2D
#-------------------------DEFINITION-------------------------
#RNA 2D module is used to compare secondary structures
#It contains type representation, methods to convert representation,
#compare different structures (some restricted to structures of same length).

#-------------------------imports-------------------------



#-------------------------exported methods-------------------------
export structure, testDotBracket, dotBracketToMountain, dotBracketToBPSet,
compareMountainDistance, fastCompareBPSet, compareHausdorff, RNAshapes



#-------------------------type definition-------------------------
immutable structure
  family::Int #needed in mccons
  dotBracket::String
  mountain::Vector{Int}
  base_pair_set::Vector{(Int,Int)}
  energy::FloatingPoint
  
  function structure(family::Int, dotBracketInput::String)
    @assert testDotBracket(dotBracketInput) 
    mountain = dotBracketToMountain(dotBracketInput)
    base_pair_set = dotBracketToBPSet(dotBracketInput)
    self = new(family, dotBracketInput, mountain, base_pair_set, -Inf)
  end
  
  function structure(family::Int, dotBracketInput::String, energy::FloatingPoint)
    @assert testDotBracket(dotBracketInput) 
    mountain = dotBracketToMountain(dotBracketInput)
    base_pair_set = dotBracketToBPSet(dotBracketInput)
    self = new(family, dotBracketInput, mountain, base_pair_set, energy)
  end    
end



#-------------------------verification methods-------------------------
function testDotBracket(dotBracket::String)
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
function dotBracketToMountain(dotBracket::String)
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



function dotBracketToBPSet(dotBracket::String)
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



function compareBPSet(bp1::Vector{(Int,Int)}, bp2::Vector{(Int,Int)})
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



function compareHausdorff(bp1::Vector{(Int,Int)}, bp2::Vector{(Int,Int)})
  function distanceBP(a::(Int, Int), b::(Int, Int))
    return maximum((abs(a[1]-b[1])), abs(a[2]-b[2]))
  end
  
  function distanceBPtoSet(a::(Int, Int), b::Vector{(Int,Int)})
    return minimum(map(x->compareHausdorff(a,x), b))
  end
  
  hausdorffLefttoRight = maximum(map(x->distanceBPtoSet(x,bp2), bp1))
  hausdorffRightToLeft = minimum(map(x->distanceBPtoSet(x,bp1), bp2))
  return maximum(hausdorffLefttoRight, hausdorffRightToLeft)
end



function RNAShape_fetchStems(structure::String)
  #assert given structure is adequate
  @assert testDotBracket(structure) == true
  list_opener = {}
  list_stems = {}
  
  i = 1
  list_stem_end = {}
    
  # separate in stems
  while i <= length(structure)
  
    #add to opening list  
    if structure[i] == '('
      push!(list_opener, i)
    
    #find the closing in the opening list
    elseif structure[i] == ')'
      current_stem = Dict()
      current_stem["opener"] = Int[]
      current_stem["closer"] = Int[]
      current_stem["open_dict"] = Dict{Int,Int}()
      
      while i <= length(structure)
	if structure[i] == ')'
	  opener = pop!(list_opener)
	  push!(current_stem["opener"], opener)
	  push!(current_stem["closer"], i)
	  current_stem["open_dict"][opener] = i
	  
	  if (!isempty(list_opener)) && (list_opener[end] in list_stem_end)
	    push!(list_stem_end, list_opener[end])
	    break
	  end
	  
	elseif structure[i] =='('
	  if (!isempty(list_opener))
	    push!(list_stem_end, list_opener[end])
	  end  
	  i -= 1
	  break
	  
	end 
        i += 1 
      end #inner while end
      push!(list_stems, current_stem)
    end
    i += 1
  end
  return list_stems
end



function RNAshapes(structure::String)
  #returns leve 1, 3 and 5 or RNAshapes (from Bielefeld)
  #fetch the stems
  list_stems = RNAShape_fetchStems(structure)
  
  #for convenience
  add(a::String, b::String) = string(a,b)

  # build the level1 for each stems
  range_occupied = {}
  dict_lvl1 = Dict()
  for stem in list_stems
    #println("for $stem in $list_stems")
    range_open  = collect([ (minimum(stem["opener"])) : (maximum(stem["opener"])) ])
    range_close = collect([ (minimum(stem["closer"])) : (maximum(stem["closer"])) ])
    range_occupied = vcat(range_occupied, range_open, range_close)
    
    #println("range_open = $range_open")
    #println("range_close = $range_close")
    #println("range_occupied = $range_occupied")
    temp_lvl1_open = ""
    temp_lvl1_close = ""
    last_opener = None
    last_closer = None
    
    for opener in sort(stem["opener"])
      if last_opener == None
	temp_lvl1_open = add(temp_lvl1_open, "[")
	temp_lvl1_close = add("]", temp_lvl1_close)
      else
	if abs(opener - last_opener) != 1
	  temp_lvl1_open = add(temp_lvl1_open, "_")
	end
	if abs(stem["open_dict"][opener] - last_closer) != 1
	  temp_lvl1_close = add("_", temp_lvl1_close)
	end
	if (endswith(temp_lvl1_open , "_")) || (beginswith(temp_lvl1_close, "_"))
	  temp_lvl1_open = add(temp_lvl1_open, "[")
	  temp_lvl1_close = add("]", temp_lvl1_close)
	end
      end
      last_opener = opener
      last_closer = stem["open_dict"][opener]
    end
    
    dict_lvl1[ minimum(stem["opener"]) ] = {"elem" => temp_lvl1_open,  "lvl5" => "["}
    dict_lvl1[ minimum(stem["closer"]) ] = {"elem" => temp_lvl1_close, "lvl5" => "]"}
  end
  
  # assemble level1
  level1 = ""
  level5 = ""
  #println("dict lv1 = $dict_lvl1")
  for i= 1:length(structure)
  
    if i in keys(dict_lvl1)
      level1 = add(level1, dict_lvl1[i]["elem"])
      level5 = add(level5, dict_lvl1[i]["lvl5"])
    end
    
    if structure[i] == '.' && (! endswith(level1, "_")) && (!(i in range_occupied))
      level1 = add(level1, "_")
    end
    
  end
  #println("lv 1 = $level1")
  level1 = replace(level1, "[_]", "[]")
  level1 = replace(level1, " ", "")
  level3 = replace(level1, "_", "")
  
  #particular edge case
  if(level5 == "")
    level5 = "_"
    level3 = "_"
  end
  return {level5, level3, level1}
end



#--module end
end