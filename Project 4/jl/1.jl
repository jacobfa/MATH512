# Let W_t be a standard Wiener Process, that is the drift parameter is zero and the Variance parameter 𝜎2 = 1.
# Suppose that we divide the interval [0,2] into L subintervals [t_i, t_{i+1}] , with 𝑡𝑖 = 𝑖𝛿𝑡 and 𝛿𝑡 = 2/𝐿
# Let 𝑊𝑖 = 𝑊(𝑡𝑖 ) 𝑎𝑛𝑑 𝛿𝑊𝑖 = 𝑊𝑖+1 − 𝑊𝑖. Verify numerically that
# a) sum_{i=0}^{L-1} |𝛿𝑊𝑖| is unbounded as 𝛿t goes to zero.
# b) sum_{i=0}^{L-1} 𝛿𝑊𝑖^2 converges to 2 in probability as 𝛿t goes to zero

using Random
using Statistics
using Plots
using LaTeXStrings

function simulate(L)
    δt = 2/L
    W = 0.0
    sum1 = 0.0
    sum2 = 0.0
    for i in 1:L
        δW = sqrt(δt) * randn()
        sum1 += abs(δW)
        sum2 += δW^2
    end
    return sum1, sum2
end

L = 1000
N = 1000
sum1 = zeros(N)
sum2 = zeros(N)
for i in 1:N
    sum1[i], sum2[i] = simulate(L)
end

println("Mean of sum1: ", mean(sum1))
println("Mean of sum2: ", mean(sum2))

L = 100:100:10000
sum_dW = zeros(length(L))
sum_dW2 = zeros(length(L))

for i in 1:length(L)
    δt = 2 / L[i]
    W = 0.0
    dW = zeros(L[i])
    for j in 1:L[i]
        dW[j] = sqrt(δt) * randn()
    end
    sum_dW[i] = sum(abs.(dW))
    sum_dW2[i] = sum(dW.^2)
end

plot(L, sum_dW, label=L"\sum_{i=0}^{L-1} |\delta W_i|", xlabel="L", ylabel="Value", title=L"Sums of $|\delta W_i|$ and $\delta W_i^2$")
plot!(L, sum_dW2, label=L"\sum_{i=0}^{L-1} \delta W_i^2")
savefig("./imgs/convergence1.png")