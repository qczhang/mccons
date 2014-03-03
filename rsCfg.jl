#random string generator from context free grammar
#based on Generating Strings at Random from a Context Free Grammar


typealias sym = Char

typealias terminal = Char

typealias nonterminals = Set{sym}

typealias terminals =  Set{terminal}


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
    for i in productions
      @assert i.input in nonterm
      @assert union(intersect(i.output, term), intersect(i.output,nonterm)) == output

    
  
type productionRule
  input::sym
  output::Set{sym}
  
  function produce(n::nonterminal)
    @assert n == input
    return output
  end
end
    
