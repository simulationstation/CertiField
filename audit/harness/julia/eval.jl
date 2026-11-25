# Evaluate a SpecialFunctions.jl / Distributions.jl function on rows of hexadecimal floats.
#   julia --project=<this folder> eval.jl <function> <in.csv> <out.txt>
# Input and output are hexadecimal floats ("%a"), so every Float64 passes exactly.
# An exception is written as NaN (the audit samples only documented domains).
using SpecialFunctions, Distributions, Printf

const F = Dict{String,Function}(
    "gamma_inc_p" => (a, x) -> gamma_inc(a, x)[1],
    "gamma_inc_q" => (a, x) -> gamma_inc(a, x)[2],
    "gamma_inc_inv" => (a, p) -> gamma_inc_inv(a, p, 1 - p),
    "beta_inc_p" => (a, b, x) -> beta_inc(a, b, x)[1],
    "beta_inc_q" => (a, b, x) -> beta_inc(a, b, x)[2],
    "beta_inc_inv" => (a, b, p) -> beta_inc_inv(a, b, p)[1],
    "logabsgamma" => x -> logabsgamma(x)[1],
    "digamma" => x -> digamma(x),
    "trigamma" => x -> trigamma(x),
    "erfcinv" => y -> erfcinv(y),
    "erfcx" => x -> erfcx(x),
    "zeta" => (s, q) -> zeta(s, q),
    "besselix" => (v, x) -> besselix(v, x),
    "besselkx" => (v, x) -> besselkx(v, x),
    "besselj" => (v, x) -> besselj(v, x),
    "normal_logccdf" => x -> logccdf(Normal(), x),
    "t_logccdf" => (x, df) -> logccdf(TDist(df), x),
    "chisq_logccdf" => (x, k) -> logccdf(Chisq(k), x),
    "f_logccdf" => (x, d1, d2) -> logccdf(FDist(d1, d2), x),
    "beta_logccdf" => (x, a, b) -> logccdf(Beta(a, b), x),
    "gamma_logccdf" => (x, a) -> logccdf(Gamma(a), x),
    "binom_logccdf" => (k, n, p) -> logccdf(Binomial(Int(n), p), k),
    "t_invlogccdf" => (lp, df) -> invlogccdf(TDist(df), lp),
    "beta_invlogccdf" => (lp, a, b) -> invlogccdf(Beta(a, b), lp),
)

f = F[ARGS[1]]
open(ARGS[3], "w") do io
    for line in eachline(ARGS[2])
        s = strip(line)
        isempty(s) && continue
        r = [parse(Float64, t) for t in split(s, ",")]
        y = try
            Float64(f(r...))
        catch
            NaN
        end
        println(io, @sprintf("%a", y))
    end
end
