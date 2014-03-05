#random string generator from context free grammar
#based on Generating Strings at Random from a Context Free Grammar
#the grammar still needs to be non ambiguous

using Base.Test


function separateIntoSymbols(s::String, nonterminals::Set{Char}, terminals::Set{Char})
  result_terminals = Set{Char}()
  result_nonterminals = Set{Char}()
  
  for i in s
  
    if i in nonterminals 
      push!(result_nonterminals, i)
    elseif i in terminals
      push!(result_terminals, i)
    else
      error("symbol $i is not in the grammar")
    end
  
  end
  
  return (result_nonterminals, result_terminals)
end






type productionRule
  input::Char
  output::Vector
  
  function productionRule{S<:String}(input::Char, output::Vector{S})
    unique_output = unique(output)
    if length(unique_output) != length(output)
      error("invalid duplicate production in the rule")
    end
    self = new(input, output)
  end
  
end



function prod(p::productionRule, s::Char)
  if s == p.input
    return p.output
  else
    error("this rule does not take $s as argument")
  end
end



function test_productionRule()
  p = productionRule('S', ["A", "a"])
  @test prod(p, 'S') == ["A", "a"]
end



#the grammar must not contain epsilon productions
type grammar
  nonterminals::Set{Char}
  terminals::Set{Char}
  start::Char
  productionRules ::Vector{productionRule}
  
  function grammar(nonterminals::Set{Char}, terminals::Set{Char}, start::Char, productionRules::Vector{productionRule})
    #the intersection of terminal and non terminal symbols must be empty
    @assert intersect(nonterminals, terminals) == Set()
    
    for i in nonterminals
      @assert int(i) in 65:90
    end

    #assert terminal and non terminal symbols used in the production rules make sense
    prod_nonterms = Set{Char}()
    prod_terms = Set{Char}()
    for rule in productionRules
      #push the nonterminal acting as input
      push!(prod_nonterms, rule.input)
      
      #separate the outputs using separateIntoSymbols, will throw an error if a symbol is not either in terminals or nonterminals
      for str in rule.output
	separateIntoSymbols(str, nonterminals, terminals)
      end
    end
    
    #assert at least one rule uses the start symbol
    @assert start in prod_nonterms
    
    #initialize the object
    self = new(nonterminals, terminals, start, productionRules)
  end
end


function grammmar_generate(G::grammar, l::Int)
  #used to generate words up to a certain length
  @assert l>0
  
  nonterminated = Dict{Int, Set{String}}()
  
  terminated = Dict{Int, Set{String}}()
  
  #we use the guarantee of non epsilon characters
  productionRules = G.productionRules
  
  
end



function isTerminated(s::String, nonterminals::Set{Char}, terminals::Set{Char})
  r = separateIntoSymbols(s, nonterminals, terminals)
  if length(r[1]) == 0
    return true
  end
  return false
end

#well formed nested parentheses
#     S -> SS
#     S -> ()
#     S -> (S)

function test_generate()
  nont = Set('S')
  term = Set('(', ')')
  st = 'S'
  prods = [productionRule('S', ["SS"])]
  push!(prods, productionRule('S', ["()"]))
  push!(prods, productionRule('S', ["(S)"]))
  g = grammar(nont, term, st, prods)
  
  nonterminated = Dict{Int, Set{String}}()
  terminated = Dict{Int, Set{String}}()
  
  #start 
  l1 = String[]
  for i in g.productionRules
  
    if i.input == g.start
      append!(l1, i.output)
    end
  end
  
  #println(l1)
  nt = filter(x->!(isTerminated(x, g.nonterminals, g.terminals)), l1)
  t =  filter(x->  isTerminated(x, g.nonterminals, g.terminals),  l1)
  
  #
  function addNonterminated(x, nonterminated::Dict{Int, Set{String}})
    nonterminated[length(x)] = push!(get(nonterminated, length(x), Set{String}()), x)
  end
  
  function addTerminated(x, terminated::Dict{Int, Set{String}})
    terminated[length(x)] =  push!(get(terminated, length(x), Set{String}()), x)
  end
  
  #
  map(x->addNonterminated(x, nonterminated), nt)
  map(x->addTerminated(x, terminated), t)
  
  #
  return(nonterminated, terminated)
end

function findNonterminals(G::grammar, s::String)
  

end


function generate(G::grammar, s::String)
  (s, G.nonterminals, G.terminals)
  
end


function generateAll(G::grammar,  level::Int)
  #
  function addNonterminated(x, nonterminated::Dict{Int, Set{String}})
    nonterminated[length(x)] = push!(get(nonterminated, length(x), Set{String}()), x)
  end
  
  function addTerminated(x, terminated::Dict{Int, Set{String}})
    terminated[length(x)] =  push!(get(terminated, length(x), Set{String}()), x)
  end
  
  
  @assert level > 0
  l = 0
  nonterminated = Dict{Int, Set{String}}()
  terminated = Dict{Int, Set{String}}()
  while l <= level
  if l == 0
    #we use the starting symbol and generate
    toExpand = Vector{String}
    
  
end


function test_grammar()
  nont = Set('S','A', 'B')
  term = Set('a')
  st = 'S'
  prods = [productionRule('S', ["A","a"])]
  push!(prods, productionRule('S', ["A","B"]))
  g = grammar(nont, term, st, prods)
  
end

#valid RNA folding grammar
#S -> A
#A -> AA | 





    
