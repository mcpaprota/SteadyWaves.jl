module NonlinearSystem

using ..Wave:WaveStruct

using NonlinearSolve

struct ConditionStruct
    condition
    range
    is_singular

    ConditionStruct(condition) = new(condition,nothing,true)

    ConditionStruct(condition,range) = new(condition,range,false)

end

function nonlinear_system_base!(du, u, conditions,compiler)

    @assert sum(u) !== NaN "$u"

    w = WaveStruct(u,compiler)

    index = 1

    for con_struct in conditions
        if con_struct.is_singular
            du[index] = con_struct.condition(w)
            index += 1
        else
            for m in con_struct.range
                du[index] = con_struct.condition(w, m)
                index += 1
            end
        end

    end

    return nothing
end

function fourier_approx_base(u,compiler,conditions)

    _nonlinear_system!(du,u,p) = nonlinear_system_base!(du,u,conditions,compiler)

    problem = NonlinearProblem(_nonlinear_system!, u)
    solution = solve(problem, RobustMultiNewton())
    u[:] = solution.u

    return WaveStruct(u,compiler)
end

end