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

function set_values(default,values)
    return combine(values,default,something_with_default(nothing))
end

function something_with_default(default)
    return (a,b) -> a !== nothing ? a : b !== nothing ? b : default
end

function reset(default,flags)
    return combine(default,flags, (d,f) -> f !== nothing ? nothing : d)
end

end