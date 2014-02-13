#MCCONS rework

#Since we work with modules in the local path
#add the local path to the load_path
push!(LOAD_PATH, chomp(readall(`pwd`)))