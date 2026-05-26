module StructOperator

function map(a::T,op) where T 
    values = [op(getfield(a,name)) for name in fieldnames(T)]
    return T(values...)
end

function map(a::T, op, state) where T
    values = Array{Any}(undef, fieldcount(T))

    for (i, name) in enumerate(fieldnames(T))
        values[i], state = op(getfield(a, name), state)
    end

    return T(values...)
end

function reduce(a::T, op, state) where T
    for name in fieldnames(T)
        state = op(getfield(a, name), state)
    end

    return state
end

function combine(a::T,b::T,op) where T
    values = [op(getfield(a,name),getfield(b,name)) for name in fieldnames(T)]
    return T(values...)
end

function safe(op)
    return (args...) -> Base.reduce((a,b) -> a||b ,args .=== nothing,init=false) ? nothing : op(args...)
end

end