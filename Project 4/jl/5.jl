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

function compute_exit_time(x0, μ, σ, t, dt, n, number_of_samples, a, b)
    exit_times = zeros(number_of_samples)
    for i in 1:number_of_samples
        dw = sqrt(dt) * randn()
        x = x0
        for j in 1:n
            x = euler_maruyama(x, dt, dw)
            dw = sqrt(dt) * randn()
            if x < a || x > b
                exit_times[i] = j * dt
                break
            end
        end
    end
    return mean(exit_times)
end

exit_times = [compute_exit_time(x0, μ, σ, t, dt, n, number_of_samples, a, b) for i in 1:1000]

histogram(exit_times, bins=50, label=L"\tau", xlabel="Value", ylabel="Frequency", title=L"Histogram of $\tau$")
savefig("./imgs/5.png")
println("Mean exit time: ", mean(exit_times))
