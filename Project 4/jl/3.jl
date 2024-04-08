
using Random
using Statistics
using Plots
using LaTeXStrings
using DifferentialEquations

# Consider the following SDE
# 𝑑𝑋(𝑡) = 𝜇𝑋(𝑡)𝑑𝑡 + 𝜎𝑋(𝑡)𝑑𝑊(𝑡) , 𝑋(0) = 3, 𝜇 = 2 , 𝜎 = 0.10
# Where 𝑡 ∈ [0,1]

# a) Show that the Euler Maruyama method has weak order of convergence equal to one. That is
# |𝐸[𝑋1] − 𝐸[𝑋(1)]| = 𝐶Δ𝑡. Here 𝑋(1) is the exact solution at time 1 and 𝑋1 is the computed solution at time
# 1

function euler_maruyama(f, g, x0, t0, tf, Δt)
    t = t0:Δt:tf
    x = zeros(length(t))
    x[1] = x0
    for i in 1:length(t)-1
        x[i+1] = x[i] + f(x[i], t[i])*Δt + g(x[i], t[i])*randn()*sqrt(Δt)
    end
    return t, x
end

function f(x, t)
    return 2*x
end

function g(x, t)
    return 0.10*x
end

function exact_solution(x0, t0, tf, Δt)
    t = t0:Δt:tf
    x = zeros(length(t))
    x[1] = x0
    for i in 1:length(t)-1
        x[i+1] = x[i]*exp(2*Δt)
    end
    return t, x
end

t0 = 0
tf = 1
Δt = 0.01
x0 = 3

t, x = euler_maruyama(f, g, x0, t0, tf, Δt)
t_exact, x_exact = exact_solution(x0, t0, tf, Δt)

plot(t, x, label="Euler Maruyama", xlabel=L"t", ylabel=L"X(t)", title="Euler Maruyama vs Exact Solution")
plot!(t_exact, x_exact, label="Exact Solution")
savefig("./imgs/3a_comparison.png")

# show that the Euler Maruyama method has weak order of convergence equal to one
# |𝐸[𝑋1] − 𝐸[𝑋(1)]| = 𝐶Δ𝑡
# where 𝑋(1) is the exact solution at time 1 and 𝑋1 is the computed solution at time 1
# 𝐸[𝑋1] = mean(x) and 𝐸[𝑋(1)] = x_exact[end]
C = abs(mean(x) - x_exact[end])/Δt
println("C = ", C)

# b) Show that the Euler Maruyama method has strong order of convergence equal to one half. That is
# 𝐸|𝑋1− 𝑋(1)| = 𝐶Δ𝑡0.5. Here 𝑋(1) is the exact solution at time 1 and 𝑋1 is the computed solution at time 1.

# show that the Euler Maruyama method has strong order of convergence equal to one half
# 𝐸|𝑋1− 𝑋(1)| = 𝐶Δ𝑡0.5
C = abs(mean(x .- x_exact))/Δt^0.5
println("C = ", C)

# plot the strong, weak and exact solutions
strong_error = [abs(x[i] - x_exact[i]) for i in 1:length(x)]
weak_error = [abs(mean(x[1:i]) - x_exact[i]) for i in 1:length(x)]
strong_solution = [x_exact[i] for i in 1:length(x)]
weak_solution = [mean(x[1:i]) for i in 1:length(x)]

plot(t, strong_error, label="Strong Error", xlabel=L"t", ylabel="Error", title="Strong Error vs Weak Error")
plot!(t, weak_error, label="Weak Error")
savefig("./imgs/3b_error.png")

# Plot the strong, weak and exact solutions
plot(t, x, label="Euler Maruyama", xlabel=L"t", ylabel=L"X(t)", title="Strong Solution vs Weak Solution vs Exact Solution")
plot!(t, strong_solution, label="Strong Solution")
plot!(t, weak_solution, label="Weak Solution")
savefig("./imgs/3b_comparison.png")
