#random string generator from context free grammar
#based on Generating Strings at Random from a Context Free Grammar


type nonterminal
  v::Char
end

type terminal
  v::Char
end

typealias nonterminals Set{nonterminal}

typealias terminals Set{terminal}


type productionRule
  input::nonterminal
  output::Set
  
  function produce(n::sym)
    if n == input
      return output
    else
      return None
    end
  end  
end


#the grammar must not contain epsilon productions

type grammar
  nonterm :: nonterminals
  term :: terminals
  start :: sym
  productions ::Set{productionRule}
  
  function grammar(syms::nonterminals, term::terminals, start::sym, productions::Set{productionRule})
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
    @assert union(intersect(nonterms, nonterm), intersect(terms, term)
    
    self = new(nonterm, term, start, productions)
  end
end
    
  

    
