module Current

using ..Params
using ..Index

function celerity_from_eulerian_current(w_c,u)
    return w_c.U(w_c,u) + w_c.c_e(w_c,u)
end
function celerity_from_eulerian_current(w)
    return w.U + w.c_e
end

function eulerian_current_net_zero_transport(w_c,u)
    return w_c.Q(w_c,u)/w_c.D(w_c,u)
end

function eulerian_current_factory(c_e,config)
    if config.cc == Params.CC_EULER

        return (w_c,u) -> 0

    elseif config.cc == Params.CC_STOKES

        return eulerian_current_net_zero_transport

    else
        throw(error("Unknown current criterion $(config.cc)"))
    end
end

end