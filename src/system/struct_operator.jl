module StructOperator

function map(a::T,op) where T 
    values = [op(getfield(a,name)) for name in fieldnames(T)]
    return T(values...)
end

function combine(a::T,b::T,op) where T
    values = [op(getfield(a,name),getfield(b,name)) for name in fieldnames(T)]
    return T(values...)
end

function safe(op)
    return (args...) -> reduce((a,b) -> a||b ,args .=== nothing,init=false) ? nothing : op(args...)
end

end