using Plots
using ForwardDiff: derivative

f = x -> exp(x)cos(8x^2)

x = range(0,1,100)
plot(x, f.(x); linewidth=2, legend=false)

r = [sqrt(π/2)/(2sqrt(2)),sqrt(-π/2+2π)/(2sqrt(2)),sqrt(π/2 + 2π)/(2sqrt(2))]
scatter!(r, zero(r))
savefig("slides/figures/slides4_funct.pdf")

