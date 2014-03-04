#random string generator from context free grammar
#based on Generating Strings at Random from a Context Free Grammar



type nonterminal
  #uppercase char for terminals
  v::Char
  function nonterminal(v::Char)
    @assert int(v) in 65:90
    self = new(v)
  end
end



type terminal
  #lowercase char for terminals
  v::Char
  function terminal(v::Char)
    @assert int(v) in 97:122
    self = new(v)
  end
end



typealias nonterminals Set{nonterminal}



typealias terminals Set{terminal}



function separateIntoSymbols(s::String)
  term = terminals()
  nonterm = nonterminals()
  for i in s
    unicodeVal = int(i)
    if unicodeVal in 65:90
      push!(nonterm, nonterminal(i))
    elseif unicodeVal in 97:122
      push!(term, terminal(i))
    else
      throw(DomainError())
    end
  end
  
  return (term, nonterm)
end




type productionRule
  input::nonterminal
  output::Vector
  
  function productionRule{S<:String}(input::nonterminal, output::Vector{S})
    self = new(input, output)
  end
  
end



function prod(p::productionRule, n::nonterminal)
  println(p)
#   println(n)
  if n == p.input
    return p.output
  else
    throw(DomainError())
  end
end



function test_productionRule()
  p = productionRule(nonterminal('S'), ["A", "a"])
end


#the grammar must not contain epsilon productions

type grammar
  nonterm :: nonterminals
  term :: terminals
  start :: nonterminal
  productions ::Set{productionRule}
  
  function grammar(syms::nonterminals, term::terminals, start::nonterminal, productions::Vector{productionRule})
    #the intersection of terminal and non terminal symbols must be empty
    @assert intersect(nonterm, term) = Set()
    
    #assert terminal and non terminal symbols used in the production rules make sense
    nonterms = Set()
    terms = Set()
    for i in productions
      inp = Set()
      outp = Set()
      for j in i.input
	push!(inp, collect(j))
      end
      
      for j in i.output
	push!(outp, collect(j))
      end
      
      nonterms = union(nonterms, inp)
      terms = union(terms, outp)
    end
    @assert union(intersect(nonterms, nonterm), intersect(terms, term))
    
    self = new(nonterm, term, start, productions)
  end
end
    
    
#valid RNA folding grammar
#S -> A
#A -> AA | 
#
  

    
