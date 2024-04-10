using Random
using Statistics
using Plots
using LaTeXStrings


# Evaluate numerically the stochastic integrals
# a) It𝑜̂ ∫ 𝑊(𝑡)𝑑𝑊(𝑡) on the interval [0,2]

function Ito_integral(N)
    t = 0.0
    W = 0.0
    sum = 0.0
    for i in 1:N
        δt = 2/N
        dW = sqrt(δt) * randn()
        sum += W * dW
        W += dW
    end
    return sum
end


N = 1000
M = 1000
sum = zeros(M)
for i in 1:M
    sum[i] = Ito_integral(N)
end

println("Mean of sum: ", mean(sum))

# histogram of results

itos = [Ito_integral(N) for i in 1:M]
histogram(itos, bins=50, label=L"\int W(t) dW(t)", xlabel="Value", ylabel="Frequency", title=L"Histogram of $\int W(t) dW(t)$")
savefig("./imgs/2a.png")


# b) Stratonovich ∫ 𝑊(𝑡) ∘ 𝑑𝑊(𝑡) on the interval [0,2]

function Stratonovich_integral(N)
    t = 0.0
    W = 0.0
    sum = 0.0
    for i in 1:N
        δt = 2/N
        dW = sqrt(δt) * randn()
        sum += (W + 0.5 * dW) * dW
        W += dW
    end
    return sum
end

N = 1000
M = 1000

sum = zeros(M)
for i in 1:M
    sum[i] = Stratonovich_integral(N)
end

println("Mean of sum: ", mean(sum))

# histogram of results

strats = [Stratonovich_integral(N) for i in 1:M]
histogram(strats, bins=50, label=L"\int W(t) \circ dW(t)", xlabel="Value", ylabel="Frequency", title=L"Histogram of $\int W(t) \circ dW(t)$")
savefig("./imgs/2b.png")

# E[∫ 𝑊(𝑡)𝑑𝑊(𝑡)] on the interval [0,2]

function E_integral(N)
    t = 0.0
    W = 0.0
    sum = 0.0
    for i in 1:N
        δt = 2/N
        dW = sqrt(δt) * randn()
        sum += W * dW
        W += dW
    end
    return sum
end

N = 1000
M = 1000

sum = zeros(M)
for i in 1:M
    sum[i] = E_integral(N)
end

println("Mean of sum: ", mean(sum))

# histogram of results

es = [E_integral(N) for i in 1:M]

histogram(es, bins=50, label=L"\mathbb{E} \left[ \int W(t) dW(t) \right]", xlabel="Value", ylabel="Frequency", title=L"Histogram of $\mathbb{E} \left[ \int W(t) dW(t) \right]$")
savefig("./imgs/2c.png")

# E[∫ 𝑊(𝑡)𝑑𝑊(𝑡)]^2 on the interval [0,2]
 
function E2_integral(N)
    t = 0.0
    W = 0.0
    sum = 0.0
    for i in 1:N
        δt = 2/N
        dW = sqrt(δt) * randn()
        sum += W * dW
        W += dW
    end
    return sum^2
end

N = 1000
M = 1000

sum = zeros(M)
for i in 1:M
    sum[i] = E2_integral(N)
end

println("Mean of sum: ", mean(sum))

# histogram of results

es2 = [E2_integral(N) for i in 1:M]

histogram(es2, bins=50, label=L"\mathbb{E} \left[ \int W(t) dW(t) \right]^2", xlabel="Value", ylabel="Frequency", title=L"Histogram of $\mathbb{E} \left[ \int W(t) dW(t) \right]^2$")

savefig("./imgs/2d.png")

# e) 𝐸 [∫ 𝑊2(𝑡)𝑑𝑡 on the interval [0, 2]

function E3_integral(N)
    t = 0.0
    W = 0.0
    sum = 0.0
    for i in 1:N
        δt = 2/N
        dW = sqrt(δt) * randn()
        sum += W^2 * δt
        W += dW
    end
    return sum
end

N = 1000
M = 1000

sum = zeros(M)
for i in 1:M
    sum[i] = E3_integral(N)
end

println("Mean of sum: ", mean(sum))

# histogram of results

es3 = [E3_integral(N) for i in 1:M]

histogram(es3, bins=50, label=L"\mathbb{E} \left[ \int W^2(t) dt \right]", xlabel="Value", ylabel="Frequency", title=L"Histogram of $\mathbb{E} \left[ \int W^2(t) dt \right]$")
savefig("./imgs/2e.png")

# For 𝑡 ∈ [0,2] evaluate ∫ 𝑊(𝑡)𝑑𝑊(𝑡) in the range [0,t], ∫ 𝑊(𝑡) ∘ 𝑑𝑊(𝑡) in the range [0,t] and 1/2 ∫ 𝑑𝑡 in the range [0,t]. What do you observe
L = 500000

function integrals(L)
    delta_t = 2 / L
    t = range(0, 2, length=L+1)
    W = zeros(L + 1)
    dW = zeros(L)
    for i in 1:L
        dW[i] = sqrt(delta_t) * randn()
        W[i + 1] = W[i] + dW[i]
    end
    return t, W, dW
end

t, W, dW = integrals(L)

plot(t[1:end-1], cumsum(W[1:end-1] .* dW), label=L"\int W(t) dW(t)", xlabel="t", ylabel="Integral", title="Stochastic integrals")
plot!(t[1:end-1], cumsum((W[1:end-1] + W[2:end]) / 2 .* dW), label=L"\int W(t) \circ dW(t)")
plot!(t[1:end-1], t[1:end-1] / 2, label=L"\frac{1}{2} \int dt")
savefig("./imgs/2f.png")