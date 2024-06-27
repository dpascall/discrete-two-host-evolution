using Random, Distributions, DataFrames, StringDistances, LinearAlgebra

function f(n::Int, beta::Float64, gamma::Float64, mu::Float64, i::Int, t_max::Float64)
    t::Float64 = 0 ##initialise time
    m::Int = 0 ##initialise mutations
    rates = [0.0, 0.0, 0.0] ##initalise rates
    df_recording = DataFrame(S = Int[], I = Int[], Evolution = Int[], t = Float64[]) ##data frame for storing events
    df_tracking = DataFrame(ID = 1:n, Status = fill(0, n), Genotype = fill("", n)) ##data frame for storing individuals
    push!(df_recording, (n-i, i, 0, t))  ##store initial conditions
    
    while t < t_max ##until pre-defined end point
        rates[1] = i*(n-i)*beta ##update infection rate
        rates[2] = i*gamma ##update recovery rate
        rates[3] = i*mu ##update evolutionary rate
        t_draw = rand(Exponential(1/sum(rates)), 1) ##sample time to next event
        t = t + t_draw[1] ##step forward time
        r_draw = rand(Categorical(rates./sum(rates)), 1) ##determine event type
        if r_draw[1] == 1 ##if infection...
            i = i + 1 ##...increase number of infected individuals
        elseif r_draw[1] == 2 ##if clearance...
            i = i - 1 ##...decrease number of infected individuals
        else ##otherwise mutate
            m = m + 1 ##increase count of evolution
        end
        push!(df_recording, (n-i, i, m, t)) ##update recording
        if i == 0 ##end if extinction
            break
        end
    end   
    return df_recording ##return results
end

function generate_seqs(nucleotides::Int) ##generate initial sequence set
    if nucleotides < 1
        throw(DomainError(nucleotides, "argument must be a positive"))
    end
    alphabet = ["A","C","G","T"] ##alphabet of nucleotides
    return join.(collect(Iterators.product(ntuple(_ -> alphabet, nucleotides)...))[:]) ##create all possible orderings of n nucleotides
end

function mutated_string(seq::String, seqs::Vector{String}, dists::Matrix{Int64}, nucleotides::Int) ##get candidate mutation
    return rand(seqs[vec(dists[:,seqs.==seq].==(nucleotides-1))]) ##sample from all sequences of hamming distance 1 (precomputed)
end

function hamming_similarity_mat(nucleotides::Int, seqs::Vector{String})
    return nucleotides .- pairwise(Hamming(), seqs) ##generate pairwise hamming distances between all inputed sequences
end

function landscape_sim_MVN(U::Matrix{Float64}, D::Vector{Float64}) ##using SVD to MVN trick
    return U * sqrt(diagm(D)) * rand(Normal(), size(U, 1)) ##x ~ Normal() then UD^(1/2)x is MVN with given correlation structure
end

function copula(draws::Vector{Float64}, target_distribution::String, p = 0.2) ##using copulas to convert from MVN to target structure
    if target_distribution == "Uniform"
        return cdf(Normal(), draws) ##marginal uniformity under cdf transform
    elseif target_distribution == "Gamma"
        return quantile(Gamma(), cdf(Normal(), draws)) ##marginal uniformity under cdf transform -> Gamma quantile
    elseif target_distribution == "ZIGamma" ##bespoke implementation of Zero-inflated Gamma quantile function
        zigamma_q = cdf(Normal(), draws) ##marginal uniformity under cdf transform
        zigamma_q[vec(zigamma_q .< p)] = fill(0.0, length(zigamma_q[vec(zigamma_q .< p)])) ##those with random draws less than the probability of 0 equal 0
        zigamma_q[vec(zigamma_q .>= p)] = cdf(Gamma(), (zigamma_q[vec(zigamma_q .>= p)] .-p ) ./ (1-p)) ##others used rescaled Gamma quantile function
        return zigamma_q
    else
        throw(DomainError(target_distribution, "argument must be an implemented distribution"))
    end
end


f(200, 0.005, 0.5, 0.5, 1, 200.0)