# Consider the following SDE:
# 𝑑𝑋(𝑡) = 𝜇𝑋(𝑡)𝑑𝑡 + 𝜎𝑋(𝑡)𝑑𝑊(𝑡) , 𝑋(0) = 3, 𝜇 = 2 , 𝜎 = 0.10
# a) Simulate (over the interval [0,20]) this stochastic process using an implicit method of the form
# 𝑋𝑛+1 = 𝑋𝑛 + (1 − 𝜃)𝛥𝑡𝑓(𝑋𝑛) + 𝜃𝛥𝑡𝑓(𝑋𝑛+1) + √𝛥𝑡𝛼𝑛𝑔(𝑋𝑛)

# using Plots
using Colors, Plots;gr()
using Random
using Statistics
using LaTeXStrings
using DifferentialEquations

function f(x, μ)
    return μ*x
end

function g(x, σ)
    return σ*x
end

function implicit_euler(μ, σ, θ, Δt, N)
    x = zeros(N)
    x[1] = 3
    for i in 2:N
        α = randn()
        x[i] = x[i-1] + (1-θ)*Δt*f(x[i-1], μ) + θ*Δt*f(x[i], μ) + sqrt(Δt)*g(x[i-1], σ)*α
    end
    return x
end

μ = 2
σ = 0.10
θ = 0.5
Δt = 0.01
N = 2000
x = implicit_euler(μ, σ, θ, Δt, N)
t = 0:Δt:20-Δt
plot(t, x, label = "Implicit Euler", xlabel = L"t", ylabel = L"X(t)", title = "Implicit Euler Method", legend = :topleft, dpi = 1000)
savefig("imgs/4implicit_euler.png")

# b) Repeat the simulation using the explicit Euler method
function explicit_euler(μ, σ, Δt, N)
    x = zeros(N)
    x[1] = 3
    for i in 2:N
        α = randn()
        x[i] = x[i-1] + Δt*f(x[i-1], μ) + sqrt(Δt)*g(x[i-1], σ)*α
    end
    return x
end

x = explicit_euler(μ, σ, Δt, N)
plot(t, x, label = "Explicit Euler", xlabel = L"t", ylabel = L"X(t)", title = "Explicit Euler Method", legend = :topleft, dpi = 1000)
plot!(t, implicit_euler(μ, σ, θ, Δt, N), label = "Implicit Euler", dpi = 1000)
savefig("imgs/4explicit_euler.png")

# c) For what values of 𝜇 𝑎𝑛𝑑 𝜎 is the SDE mean-square stable.
# The SDE is mean-square stable if the following condition is satisfied:
# It is called mean-square stable if, for every ε > 0, there exists δ > 0 such that
# E[|X(t)|^2] ≤ ε for all t ≥ 0 whenever E[|X(0)|^2] ≤ δ.

# Getting the values of μ and σ that make the SDE mean-square stable
# ranges of μ and σ
μ = 0:0.1:75
σ = 0:0.1:20
stable = []
for μ in μ
    for σ in σ
        y = implicit_euler(μ, σ, θ, Δt, N)
        if mean(y.^2) < 1
            push!(stable, (μ, σ))
        end
    end
end

# plotting the values of μ and σ that make the SDE mean-square stable
μ = [x[1] for x in stable]
σ = [x[2] for x in stable]

edges_x = 0:0.1:20
edges_y = 0:0.1:75
histogram2d(μ, σ, bins = (edges_y, edges_x), xlabel = L"\mu", ylabel = L"\sigma", title = "Mean Square Stable", dpi = 1000, c = :blues)
savefig("imgs/4mean_square_stable.png")

# d) For what values of 𝜃 is the implicit method mean-square stable
theta = 0:0.1:1
stable = []
# set μ and σ to 2 and 0.10 respectively
μ = 2
σ = 0.10
for t in theta
    z = implicit_euler(μ, σ, t, Δt, N)
    if mean(z.^2) < 1
        push!(stable, t)
    end
end

edges_x = 0:0.1:1
histogram(stable, bins = edges_x, xlabel = L"\theta", ylabel = "Frequency", title = "Mean Square Stable", dpi = 1000, c = :blues)
savefig("imgs/4mean_square_stability_theta.png")

# For what values of 𝜇 𝑎𝑛𝑑 𝜎 is the SDE asymptotically stable.

# set μ and σ to 2 and 0.10 respectively
# stable if limn→∞ |Yn| = 0, with probability one, for any X0.

# Getting the values of μ and σ that make the SDE asymptotically stable
μ = 0:0.1:75
σ = 0:0.1:20
stable = []
probabilities = []
for μ in μ
    for σ in σ
        y = implicit_euler(μ, σ, θ, Δt, N)
        prob = 0
        for i in 1:N
            if abs(y[i]) < 1
                prob += 1
            end
        end
        if prob == N
            push!(stable, (μ, σ))
        end
    end
end

# plotting the values of μ and σ that make the SDE asymptotically stable
μ = [x[1] for x in stable]
σ = [x[2] for x in stable]

edges_x = 0:0.1:20
edges_y = 0:0.1:75

histogram2d(μ, σ, bins = (edges_y, edges_x), xlabel = L"\mu", ylabel = L"\sigma", title = "Asymptotically Stable", dpi = 1000, c = :blues)
savefig("imgs/4asymptotically_stable.png")


# For what values of 𝜃 is the implicit method asymptotically stable
theta = 0:0.1:1
stable = []
# set μ and σ to 2 and 0.10 respectively
μ = 2
σ = 0.10
for t in theta
    z = implicit_euler(μ, σ, t, Δt, N)
    if abs(mean(z)) < 1
        push!(stable, t)
    end
end

edges_x = 0:0.1:1
histogram(stable, bins = edges_x, xlabel = L"\theta", ylabel = "Frequency", title = "Asymptotically Stable", dpi = 1000, c = :blues)
savefig("imgs/4asymptotically_stability_theta.png")
