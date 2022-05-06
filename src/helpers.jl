function lastline(io)
    local line
    for l in eachline(io)
        line = l
    end
    line
end