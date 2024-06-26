using Random, Distributions, DataFrames, StringDistances, LinearAlgebra

function f(n::Int, beta::Float64, gamma::Float64, i::Int, t_max::Float64)
    t::Float64 = 0
    df = DataFrame(S=Int[], I=Int[], t=Float64[])
    push!(df, (n-i, i, t)) 
    
    while t < t_max
        current_total_rate::Float64 = i*(n-i)*beta + i*gamma
        t_draw = rand(Exponential(1/current_total_rate), 1)
        t = t + t_draw[1]
        r_draw = rand(Bernoulli((i*(n-i)*beta)/(i*(n-i)*beta + i*gamma)), 1)
        if r_draw[1] == 1
            i = i + 1
        else
            i = i - 1
        end
        push!(df, (n-i, i, t))
        if i == 0
            break
        end
        print(t)
        print("\n")
    end   
    return df
end

function generate_seqs(nucleotides::Int)
    if nucleotides < 1
        throw(DomainError(nucleotides, "argument must be a positive"))
    end
    alphabet = ["A","C","G","T"]
    return join.(collect(Iterators.product(ntuple(_ -> alphabet, nucleotides)...))[:])
end

function mutated_string(seq::String, seqs::Vector{String}, dists::Matrix{Int64}, nucleotides::Int)
    return rand(seqs[vec(dists[:,seqs.==seq].==(nucleotides-1))])
end

function hamming_similarity_mat(nucleotides::Int, seqs::Vector{String})
    return nucleotides .- pairwise(Hamming(), seqs)
end

function landscape_sim_MVN(U::Matrix{Float64}, D::Vector{Float64}) ##using SVD to MVN trick
    return U * sqrt(diagm(D)) * rand(Normal(), size(U, 1))
end

function copula(draws::Vector{Float64}, target_distribution::String, p = 0.2)
    if target_distribution == "Uniform"
        return cdf(Normal(), draws)
    elseif target_distribution == "Gamma"
        return quantile(Gamma(), cdf(Normal(), draws))
    elseif target_distribution == "ZIGamma"
        zigamma_q = cdf(Normal(), draws)
        zigamma_q[vec(zigamma_q .< p)] = fill(0.0, length(zigamma_q[vec(zigamma_q .< p)]))
        zigamma_q[vec(zigamma_q .>= p)] = cdf(Gamma(), (zigamma_q[vec(zigamma_q .>= p)] .-p ) ./ (1-p))
        return zigamma_q
    else
        throw(DomainError(target_distribution, "argument must be an implemented distribution"))
    end
end