#random string generator from context free grammar
#based on Generating Strings at Random from a Context Free Grammar




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
end


#the grammar must not contain epsilon productions

type grammar
  nonterminals::Set{Char}
  terminals::Set{Char}
  start::Char
  productionRules ::Set{productionRule}
  
  function grammar(nonterminals::Set{Char}, terminals::Set{Char}, start::Char, productionRules::Set{productionRule})
    #the intersection of terminal and non terminal symbols must be empty
    @assert intersect(nonterm, term) = Set()
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
    
    #asser that at least one rule uses the start symbol
    @assert start in prod_nonterms
    
    self = new(nonterminals, terminals, start, productionRules)
  end
end
    
    
#valid RNA folding grammar
#S -> A
#A -> AA | 
#
  

    
