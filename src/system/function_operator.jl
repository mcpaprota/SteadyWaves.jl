module FunctionOperator

struct F
    func

    # Preserve existing F
    F(f::F) = f

    # Wrap callable functions
    F(f::Function) = new(f)

    # Promote numbers to constant functions
    F(value::Number) = new((args...) -> value)

    F(value) = value
end

# Call overload
(f::F)(args...) = f.func(args...)

# -------------------
# Arithmetic operators
# -------------------

# Addition
Base.:+(f::F, g::F) = F((args...) -> f(args...) + g(args...))
Base.:+(f::F, c::Number) = F((args...) -> f(args...) + c)
Base.:+(c::Number, f::F) = F((args...) -> c + f(args...))

# Subtraction
Base.:-(f::F, g::F) = F((args...) -> f(args...) - g(args...))
Base.:-(f::F, c::Number) = F((args...) -> f(args...) - c)
Base.:-(c::Number, f::F) = F((args...) -> c - f(args...))

# Unary minus
Base.:-(f::F) = F((args...) -> -f(args...))

# Multiplication
Base.:*(f::F, g::F) = F((args...) -> f(args...) * g(args...))
Base.:*(f::F, c::Number) = F((args...) -> f(args...) * c)
Base.:*(c::Number, f::F) = F((args...) -> c * f(args...))

# Division
Base.:/(f::F, g::F) = F((args...) -> f(args...) / g(args...))
Base.:/(f::F, c::Number) = F((args...) -> f(args...) / c)
Base.:/(c::Number, f::F) = F((args...) -> c / f(args...))

# Power
Base.:^(f::F, g::F) = F((args...) -> f(args...) ^ g(args...))
Base.:^(f::F, c::Number) = F((args...) -> f(args...) ^ c)
Base.:^(c::Number, f::F) = F((args...) -> c ^ f(args...))

end