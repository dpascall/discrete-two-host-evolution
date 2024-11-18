#!/usr/bin/env julia
using Random, Distributions, DataFrames, StringDistances, LinearAlgebra, CSV, Tables, StatsBase
Random.seed!(1234)

function f(h1::Int, h2::Int, beta::Float64, gamma::Float64, mu::Float64, i0::Int, t_max::Float64, nucleotides::Int, fitness::DataFrame, Hamming::Matrix, Ne::Int)
    t::Float64 = 0 ##initialise time
    m::Int = 0 ##initialise mutations
    rates = [0.0, 0.0, 0.0] ##initalise rates
    seqs = generate_seqs(nucleotides) ##get sequences
    i = i0
    n = h1 + h2
    df_recording = DataFrame(S = Int[], I = Int[], t = Float64[]) ##data frame for storing events
    df_tracking = DataFrame(ID = 1:n, Status = fill(0, n), Genotype = fill("", n), Host = [fill(1, h1); fill(2, h2)]) ##data frame for storing individuals
    df_tracking.Genotype[sample(1:n, i0, replace = false)] .= sample(seqs, 1) ##random single genotype to start -- can currently start from 0 fitness genotypes
    df_tracking.Status[df_tracking.Genotype .!= ""] .= 1 ##randomly set initial infections
    df_genotypes = DataFrame(Genotypes = Vector{Vector{String}}(), t = Vector{Float64}()) ## store genotypes over time
    push!(df_genotypes, (df_tracking.Genotype[df_tracking.Genotype .!= ""], t)) ##store starting genotypes

    push!(df_recording, (n-i0, i0, t))  ##store initial conditions
    
    while t < t_max ##until pre-defined end point
        rates[1] = i*(n-i)*beta ##update infection rate
        rates[2] = i*gamma ##update recovery rate
        rates[3] = i*mu ##update evolutionary rate
        t_draw = rand(Exponential(1/sum(rates)), 1) ##sample time to next event
        t = t + t_draw[1] ##step forward time
        if t > t_max ##end if event pushes past max time
            break
        end
        r_draw = rand(Categorical(rates./sum(rates)), 1) ##determine event type
        if r_draw[1] == 1 ##if infection...
            i = i + 1 ##...increase number of infected individuals
            i_ID = sample(findall(==(1), df_tracking.Status), 1, replace = false) ##randomly select infector candidate
            s_ID = sample(findall(==(0), df_tracking.Status), 1, replace = false) ##randomly select infection candidate
            df_tracking.Status[s_ID] .= 1 ##update infection status
            df_tracking.Genotype[s_ID] .= df_tracking.Genotype[i_ID] ##update genotype
        elseif r_draw[1] == 2 ##if clearance...
            i = i - 1 ##...decrease number of infected individuals
            i_ID = sample(findall(==(1), df_tracking.Status), 1, replace = false)
            df_tracking.Status[i_ID] .= 0 ##randomly select clearance candidate
            df_tracking.Genotype[i_ID] .= "" ##update genotype
        else ##otherwise mutate
            m_ID = sample(findall(==(1), df_tracking.Status), 1, replace = false)
            df_tracking.Genotype[m_ID] .= consensus_test(df_tracking.Genotype[m_ID][1], mutated_string(df_tracking.Genotype[m_ID][1], fitness.Sequence, Hamming), fitness, df_tracking.Host[m_ID][1], Ne)
            push!(df_genotypes, (df_tracking.Genotype[df_tracking.Genotype .!= ""], t))
        end
        push!(df_recording, (n-i, i, t)) ##update recording
        if i == 0 ##end if extinction
            break
        end
    end   
    return df_recording, df_tracking, df_genotypes ##return results
end

function consensus_test(resident_seq::String, replacement_seq::String, fitness::DataFrame, host_type::Int, Ne::Int) ##test if new mutant takes over
    if host_type == 1
        within_fitness_1 = fitness.Host1Within[fitness.Sequence .== resident_seq][1]
        within_fitness_2 = fitness.Host1Within[fitness.Sequence .== replacement_seq][1]
        return selection_algorithm(resident_seq, replacement_seq, within_fitness_1, within_fitness_2, Ne)
    else
        within_fitness_1 = fitness.Host2Within[fitness.Sequence .== resident_seq][1]
        within_fitness_2 = fitness.Host2Within[fitness.Sequence .== replacement_seq][1]
        return selection_algorithm(resident_seq, replacement_seq, within_fitness_1, within_fitness_2, Ne)
    end
end

function selection_algorithm(resident_seq::String, replacement_seq::String, within_fitness_1::Float64, within_fitness_2::Float64, Ne::Int) ##calculation for replacement
    if within_fitness_2 ≈ 0 ## a float is never equal to 0
        return resident_seq ##if defective variant, return resident
    elseif (within_fitness_2 - within_fitness_1)/within_fitness_1 > (1/Ne)
        return replacement_seq ##if old seq fitter and greater than threshold return new seq
    elseif (within_fitness_2 - within_fitness_1)/within_fitness_1 < -(1/Ne)
        return resident_seq ##if new seq fitter and greater than threshold return old seq
    else
        s = (within_fitness_2 - within_fitness_1)/within_fitness_1
        return sample([replacement_seq, resident_seq], Weights([((1-exp(-s))/(1-exp(-Ne*s))), 1-((1-exp(-s))/(1-exp(-Ne*s)))])) 
        ##if (1/Ne) > s > -(1/Ne) sample probability from under nearly neutral regime
    end
end

function mutated_string(seq::String, seqs::Vector{String}, dists::Matrix{Int64}) ##get candidate mutation
    return sample(seqs[vec(dists[:,seqs.==seq].==1)]) ##sample from all sequences of hamming distance 1 (precomputed)
end

function landscape_sim_MVN(factor::Matrix{Float64}) ##using SVD to MVN trick
    return factor * rand(Normal(), size(factor, 1)) ##x ~ Normal() then UD^(1/2)x is MVN with given correlation structure
end

function generate_seqs(nucleotides::Int) ##generate initial sequence set
    if nucleotides < 1
        throw(DomainError(nucleotides, "argument must be a positive"))
    end
    alphabet = ["A","C","G","T"] ##alphabet of nucleotides
    return join.(collect(Iterators.product(ntuple(_ -> alphabet, nucleotides)...))[:]) ##create all possible orderings of n nucleotides
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

function fitness_calc(factor::Matrix{Float64}, target_distribution::String)
    fitness = hcat(generate_seqs(Int.(log(size(factor)[1]/4)/log(4))), reshape(collect(copula(landscape_sim_MVN(factor), target_distribution)), Int.(size(factor)[1]/4), 4))
    fitness = DataFrame(fitness, :auto)
    rename!(fitness, ["Sequence", "Host1Within", "Host2Within", "Host1Between", "Host2Between"])
    fitness[!,:"Sequence"] = convert.(String,fitness[!,:"Sequence"])
    fitness[!,:"Host1Within"] = convert.(Float64,fitness[!,:"Host1Within"])
    fitness[!,:"Host2Within"] = convert.(Float64,fitness[!,:"Host2Within"])
    fitness[!,:"Host1Between"] = convert.(Float64,fitness[!,:"Host1Between"])
    fitness[!,:"Host2Between"] = convert.(Float64,fitness[!,:"Host2Between"])
    return fitness
end

Hamming = CSV.File("6HammingDistance.csv", header = 0) |> Tables.matrix
factor = CSV.File("factor_6_1.0_0.5_1.0.csv", header = 0) |> Tables.matrix

fitness = fitness_calc(factor, "ZIGamma")

results = f(100, 100, 0.005, 0.5, 2.0, 10, 200.0, Int.(log(size(Hamming)[1])/log(4)), fitness, Hamming, 40)

fitness_trajectory = DataFrame(repeat = Vector{Int}(), t = Vector{Float64}(), MeanHost1Fitness = Vector{Float64}(), MeanHost2Fitness = Vector{Float64}())

for i in 1:nrow(results[[3]][1])
    push!(fitness_trajectory, (3, results[[3]][1][i,2], mean(fitness.Host1Within[indexin(results[[3]][1][i,1], fitness.Sequence)]), mean(fitness.Host2Within[indexin(results[[3]][1][i,1], fitness.Sequence)])))
end

CSV.write("trajectories.csv",fitness_trajectory)

###generate data for poster figure

factor_paths = ["factor_6_1.0_-1.0_1.0.csv", "factor_6_1.0_-0.75_1.0.csv", "factor_6_1.0_-0.5_1.0.csv", "factor_6_1.0_-0.25_1.0.csv", "factor_6_1.0_0.0_1.0.csv", "factor_6_1.0_0.25_1.0.csv", "factor_6_1.0_0.5_1.0.csv", "factor_6_1.0_0.75_1.0.csv", "factor_6_1.0_1.0_1.0.csv"]
correlations = [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]

timetopercentile = DataFrame(correlation = Vector{Float64}(), Time1Host = Vector{Union{Missing, Float64}}(), Time2Host = Vector{Union{Missing, Float64}}())

t_twohosts = missing
t_onehost = missing

for i in 1:9
    factor = CSV.File(factor_paths[i], header = 0) |> Tables.matrix
    for j in 1:50
        fitness = fitness_calc(factor, "ZIGamma")
        results_twohosts = f(190, 10, 0.005, 0.5, 2.0, 10, 200.0, Int.(log(size(Hamming)[1])/log(4)), fitness, Hamming, 40)
        results_onehost = f(200, 0, 0.005, 0.5, 2.0, 10, 200.0, Int.(log(size(Hamming)[1])/log(4)), fitness, Hamming, 40)
        targets = fitness.Sequence[fitness.Host1Within .> cdf(Gamma(), (0.99-0.2)/(1-0.2))]
        twohoststep1 = map(x ->  indexin(x, targets), results_twohosts[[3]][1].Genotypes)
        twohoststep2 = map(x ->  filter(!isnothing, x), twohoststep1)
        twohoststep3 = map(x -> length(x), twohoststep2)
        if sum(twohoststep3) > 0
            t_twohosts = results_twohosts[[3]][1][findfirst(x -> x>0, twohoststep3),2]
        else
            t_twohosts = missing
        end
        onehoststep1 = map(x ->  indexin(x, targets), results_onehost[[3]][1].Genotypes)
        onehoststep2 = map(x ->  filter(!isnothing, x), onehoststep1)
        onehoststep3 = map(x -> length(x), onehoststep2)
        if sum(onehoststep3) > 0
            t_onehost = results_onehost[[3]][1][findfirst(x -> x>0, onehoststep3),2]
        else
            t_onehost = missing
        end
        push!(timetopercentile, (correlations[i], t_onehost, t_twohosts))
        t_twohosts = missing
        t_onehost = missing   
    end
end