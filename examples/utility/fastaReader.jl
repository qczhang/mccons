function fastaRead(fileName)
    #reads a fasta file and returns Vector{(String, Vector{String})}
    f = open(fileName, "r")
    lines = readlines(f)
    close(f)
    
    result = (String, Vector{String})[]
    name = ""
    data = String[]
    for line in lines
        if beginswith(line, ";")
            continue
        elseif beginswith(line, ">")
            push!(result, (name, data))
            name = chomp(string(split(line, ">")[2]))
            data =  String[]
        elseif line != "\n"
            push!(data, chomp(line))
        end
    end
    push!(result, (name, data))
    result[2:end]
end