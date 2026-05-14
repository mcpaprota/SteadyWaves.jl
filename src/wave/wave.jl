module Wave

using ..Index:IndexStruct
using ..Surface:SurfaceStruct, Surface
using ..Velocity:VelocityStruct, velocity_struct_factory
using ..Params
using ..StructOperator: combine, map, safe

"""
Structure with wave properties:
- `eta`: structure with properties related to free surface
- `v`: structure with properties related to velocity
- `D`: depth
- `C`: celerity
- `R`: wave-related parameter
- `H`: wave height
- `U`: average velocity
- `Q`: flow-related parameter
- `N`: number of control points along a wave
- `L`: wavelength
- `T`: wave period
- `F`: wave flux
- `P`: wave pressure function `P(x, z)`
- `sigma`: surface tension coefficient
- `raw`: internal representesion of equation system variables
"""
struct WaveStruct
    eta
    v
    D
    C
    R
    H
    U
    Q
    N
    L
    T
    F
    P
    sigma
    raw
    WaveStruct(eta,v,D,C,R,H,U,Q,N,L,T,F,P,sigma,raw) = new(eta,v,D,C,R,H,U,Q,N,L,T,F,P,sigma,raw)

    # create struct that replaces values given as key words in default
    WaveStruct(;
        eta = nothing,
        v = nothing,
        D = nothing,
        C = nothing,
        R = nothing,
        H = nothing,
        U = nothing,
        Q = nothing,
        N = nothing,
        L = nothing,
        T = nothing,
        F = nothing,
        P = nothing,
        sigma = nothing,
        raw = nothing,
    )= new(eta,v,D,C,R,H,U,Q,N,L,T,F,P,sigma,raw)

    # creates compiler X = (w_c, u) -> u[idx.X] from IndexStruct
    WaveStruct(idx::IndexStruct,config::Params.ConfigStruct) = new(
        SurfaceStruct(idx,config.eta_type),
        velocity_struct_factory(idx,config),
        (w_c, u) -> u[idx.D],
        (w_c, u) -> idx.C == 0 ? nothing : u[idx.C],
        (w_c, u) -> u[idx.R],
        (w_c, u) -> u[idx.H],
        (w_c, u) -> u[idx.U],
        (w_c, u) -> u[idx.Q],
        (w_c, u) -> idx.N,
        (w_c, u) -> nothing,
        (w_c, u) -> nothing,
        (w_c, u) -> nothing,
        (w_c, u) -> nothing,
        (w_c, u) -> 0,
        (w_c, u) -> u,
    )

    # create structure from array u and compiler X = compiler.X(compiler,u)
    WaveStruct(u,compiler::WaveStruct,inner_compiler::WaveStruct) = begin
        return map(compiler,v->v(inner_compiler,u))
    end

    WaveStruct(u,compiler::WaveStruct) = WaveStruct(u,compiler,compiler)
end

function set_compilator_values(default::WaveStruct,dless::WaveStruct)
    wrapped = map(dless,safe(dl -> typeof(dl) <: Real ? (w_c,u) -> dl : dl))

    return combine(wrapped,default,something)
end

function set_compilator_values(default::WaveStruct,dim::WaveStruct,df::WaveStruct)
    dless = combine(dim,df,safe((dim,df)-> (w_c, u) -> dim * df(w_c,u)))

    return combine(dless,default,something)
end

function set_values(default,values)
    return combine(values,default,something_with_default(nothing))
end

function something_with_default(default)
    return (a,b) -> a !== nothing ? a : b !== nothing ? b : default
end

end