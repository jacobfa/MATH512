using Random
using Statistics
using Plots
using LaTeXStrings
using DifferentialEquations


# Consider the following SDE:
# 𝑑𝑋(𝑡) = 𝜇𝑋(𝑡)𝑑𝑡 + 𝜎𝑋(𝑡)𝑑𝑊(𝑡) , 𝑋(0) = 2, 𝜇 = 0.1, 𝜎 = 0.15

# Let 𝑎 = 0.5 𝑎𝑛𝑑 𝑏 = 3.
# Compute the mean exit time function 𝑣(𝑥) for 𝑥 ∈ [0.5, 3]

# We will use the Euler Maruyama method to solve the SDE

x0 = 2
μ = 0.1
σ = 0.15
a = 0.5
b = 3
t = 1
dt = 0.01
n = Int(t / dt)
number_of_samples = 10000

function euler_maruyama(x, dt, dw)
    return x + μ * x * dt + σ * x * dw
end

function exact_solution(x0, μ, σ, t, dw)
    return x0 * exp((μ - 0.5 * σ^2) * t + σ * dw)
end

# Choose a step size Δt
# Choose a number of paths, M
# for s = 1to M
    #  Set tn = 0 and Xn = X0 While Xn > aand Xn < b  Compute a N(0,1) sample ξn 
    #  Replace Xn by Xn+∆tf(Xn)+√∆tξng(Xn) 
    #  Replace tn by tn+∆t 
    # end 
    # set Ts exit=tn−1/2 ∆t 
# end set aM= 1M M s=1Ts exit 
# setb2 M= 1 M−1 M s=1(Ts exit−aM)2

function mean_exit_time(x0, μ, σ, a, b, t, dt, number_of_samples)
    aM = 0
    b2M = 0
    vals = []
    for s = 1:number_of_samples
        Xn = x0
        tn = 0
        while Xn > a && Xn < b
            ξn = randn()
            Xn = euler_maruyama(Xn, dt, ξn)
            tn += dt
        end
        Ts_exit = tn - 0.5 * dt
        push!(vals, Ts_exit)
        aM += Ts_exit
    end
    aM /= number_of_samples
    for s = 1:number_of_samples
        b2M += (vals[s] - aM)^2
    end
    b2M /= number_of_samples - 1
    return vals, aM, b2M
end

vals, aM, b2M = mean_exit_time(x0, μ, σ, a, b, t, dt, number_of_samples)
println("Mean exit time: ", aM)
println("Variance of exit time: ", b2M)

# Plot the mean exit time function 𝑣(𝑥) for 𝑥 ∈ [0.5, 3]
x = 0.5:0.01:3
v = []
for i in x
    local vals = []
    local aM = 0
    local b2M = 0
    vals, aM, b2M = mean_exit_time(i, μ, σ, a, b, t, dt, number_of_samples)
    push!(v, aM)
end

plot(x, v, label = "Mean exit time", xlabel = L"x", ylabel = L"v(x)", title = "Mean exit time function", legend = :topleft, dpi = 1000)
savefig("imgs/5mean_exit_time.png")

# plot the histogram of the exit times
histogram(vals, label = "Exit times", xlabel = L"t", ylabel = "Frequency", title = "Histogram of exit times", legend = :topleft, dpi = 1000)
savefig("imgs/5exit_times_histogram.png")
