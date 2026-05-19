# SPDX-License-Identifier: MIT

# Functions for calculating the output
module Output
using ..Physics

using ..Index
using ..Wave: WaveStruct
using ..Indirect

"""
    wave_period(w)

Calculate wave period `T` from solution `w`.

Dimensionality of result depend on dimensionality of struct `w`

"""
function wave_period(w)
    return w.T
end

"""
    wavelength(w)

Calculate wavelength `L` from solution `w`.

Dimensionality of result depend on dimensionality of struct `w`

"""
function wavelength(w)
    return w.L
end

"""
    wave_number(w)

Calculate wave number `K` from solution `w`.

Dimensionality of result depend on dimensionality of struct `w`

"""
function wave_number(w)
    return 2pi/w.L
end

"""
    wave_height(w)

Calculate wave height `H` from solution `w`.

Dimensionality of result depend on dimensionality of struct `w`

"""
function wave_height(w)
    return w.H
end

"""
    vertical_velocity(w, kx, kz)

Calculate vertical velocity from solution `w`
at coordinates (`kx`, `kz`) or (`x`, `z`).

If `w` is dimensional provide dimensional (`x`, `z`).
If not provide (`k*x`, `k*z`).

"""
function vertical_velocity(w,kx,kz)
    return w.v.z(kx,kz)
end

"""
    horizontal_velocity(w, kx, kz)

Calculate horizontal velocity from solution `w`
at coordinates (`kx`, `kz`) or (`x`, `z`).

If `w` is dimensional provide dimensional (`x`, `z`).
If not provide (`k*x`, `k*z`).

"""
function horizontal_velocity(w,kx,kz)
    return w.v.x(kx,kz)
end

"""
    pressure(w, kx, kz)

Calculate pressure from solution `w`
at coordinates (`kx`, `kz`) or (`x`, `z`).

If `w` is dimensional provide dimensional (`x`, `z`).
If not provide (`k*x`, `k*z`).

"""
function pressure(w,kx,kz)
    return w.P(kx,kz)
end

"""
    wave_period(w::WaveStruct)

Calculate wave period `T` from wave struct `w`.

"""
function wave_period(w::WaveStruct)
    return w.T
end

"""
    wavelength(w::WaveStruct)

Calculate wavelength `L` from wave struct `w`.

"""
function wavelength(w::WaveStruct)
    return w.L
end

"""
    wave_power(w::WaveStruct)

Calculate wave power `F` from wave struct `w`.

"""
function wave_power(w::WaveStruct)
    return w.F
end

"""
    elevation(w::WaveStruct, kx)

    elevation(w::WaveStruct, kx, t)
    elevation(w::WaveStruct, kx, t, c)

Calculate free surface elevation from solution `w`
at coordinate `kx` or `x`.

If only `x` is provided returns value from wave perspective.

If `t` is provided returns value from static perspective.

If `t, s` is provided returns value from an arbitrary perspective.

If `w` is dimensional provide dimensional `x`.
If not provide `k*x`.

"""
function elevation(w::WaveStruct,kx)
    return w.eta.z(kx)
end

function elevation(w::WaveStruct,kx,t)
    return elevation(w,kx+w.C*t)
end

function elevation(w::WaveStruct,kx,t,c)
    return elevation(w,kx+c*t)
end

"""
    surface_tension(w, kx)

    surface_tension(w, kx, t)
    surface_tension(w, kx, t, c)

Calculate surface tension from solution `w`
at coordinate `kx` or `x`.

If only `x` is provided returns value from wave perspective.

If `t` is provided returns value from static perspective.

If `t, s` is provided returns value from an arbitrary perspective.

If `w` is dimensional provide dimensional `x`.
If not provide `k*x`.

"""
function surface_tension(w,kx)
    return Indirect.indirect_surface_tension(
        w.sigma,
        w.eta.z.dz_dx_1(kx),
        w.eta.z.dz_dx_2(kx)
    )
end

function surface_tension(w,kx,t)
    return surface_tension(w,kx+w.C*t)
end

function surface_tension(w,kx,t,c)
    return surface_tension(w,kx+c*t)
end

end