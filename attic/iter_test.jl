struct Squares
    count::Int
end

function Base.iterate(S::Squares, state=1)
    println("state", state)
    if state > S.count
        return nothing
    else
        return (state*state, state+1)
    end
end

for i in Squares(7)
    println(typeof(i), i)
end
