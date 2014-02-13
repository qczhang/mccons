# SIDENOTES
# general observations or remarks




#On modules and packaging... http://docs.julialang.org/en/latest/manual/modules/

# The statement using Lib means that a module called Lib will be available 
# for resolving names as needed. When a global variable is encountered that
# has no definition in the current module, the system will search for it in 
# Lib and import it if it is found there. This means that all uses of that 
# global within the current module will resolve to the definition of that 
# variable in Lib.
# 
# The statement import BigLib: bar, baz means that the names bar and baz from
# the BigLib module will be available as needed (but no other names).


#Conventions
#f!(x) implies that f modifies x
#don't put () around conditional statement evaluation



