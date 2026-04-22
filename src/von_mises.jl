# algorithm from
#     DJ Best & NI Fisher (1979). Efficient Simulation of the von Mises
#     Distribution. Journal of the Royal Statistical Society. Series C
#     (Applied Statistics), 28(2), 152-157.
#     from https://github.com/JuliaStats/Distributions.jl/blob/f852803aa1a90c330e66a70e8d16847744cb8338/src/samplers/vonmises.jl
function rand_von_mises(μ::Float64, κ::Float64)
    τ::Float64 = 1.0 + sqrt(1.0 + 4 * abs2(κ))
    ρ::Float64 = (τ - sqrt(2.0 * τ)) / (2.0 * κ)
    r::Float64 = (1.0 + abs2(ρ)) / (2.0 * ρ)

    f::Float64 = 0.0
    local x::Float64
    if κ > 700.0
        x = μ + randn(Float64) / sqrt(κ)
    else
        while true
            t::Float64, u::Float64 = 0.0, 0.0
            while true
                d::Float64 = abs2(rand(Float64) - 0.5)
                e::Float64 = abs2(rand(Float64) - 0.5)
                if d + e <= 0.25
                    t = d / e
                    u = 4 * (d + e)
                    break
                end
            end
            z::Float64 = (1.0 - t) / (1.0 + t)
            f = (1.0 + r * z) / (r + z)
            c::Float64 = κ * (r - f)
            if c * (2.0 - c) > u || log(c / u) + 1 >= c
                break
            end
        end
        acf::Float64 = acos(f)
        x = μ + (rand(Bool) ? acf : -acf)
    end
    return x
end