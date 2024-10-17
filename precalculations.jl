#!/usr/bin/env julia
##takes length of nucleotide string as an argument
using Random, Distributions, StringDistances, DelimitedFiles, LinearAlgebra

function generate_seqs(nucleotides::Int) ##generate initial sequence set
    if nucleotides < 1
        throw(DomainError(nucleotides, "argument must be a positive"))
    end
    alphabet = ["A","C","G","T"] ##alphabet of nucleotides
    return join.(collect(Iterators.product(ntuple(_ -> alphabet, nucleotides)...))[:]) ##create all possible orderings of n nucleotides
end

function hamming_similarity_mat(nucleotides::Int, seqs::Vector{String})
    return nucleotides .- pairwise(Hamming(), seqs) ##generate pairwise hamming distances between all inputed sequences
end

function calculate_covariance_factor(decomp_mu::Eigen, decomp_rho::Eigen, decomp_hamming::Eigen)
    decomp_vectors = kron(decomp_mu.vectors, kron(decomp_rho.vectors, decomp_hamming.vectors)) ##generate final covariance matrix eigenvectors
    decomp_values = kron(diagm(decomp_mu.values), kron(diagm(decomp_rho.values), diagm(decomp_hamming.values))) ##generate final covariance matrix eigenvalues
    decomp_values[decomp_values .< 0] .= 0 ##fix negative eigenvalues to 0 - assumes true matrix PSD
    return decomp_vectors * sqrt.(decomp_values) ##return factor
end

nucleotides = parse(Int, ARGS[1])

seqs = generate_seqs(nucleotides)
sim_mat = hamming_similarity_mat(nucleotides, seqs)

name = string(nucleotides) * "HammingDistance.csv"

open(joinpath(dirname(@__FILE__), name), "w") do io
    writedlm(io, sim_mat, ',') ##write out hamming mat
end

for index in [1.0, 5.0, Inf] ##perform SVDs for each of the conditions
    decomp_hamming = eigen((sim_mat./nucleotides).^index) ##use decomposability to calculate each separately
    for rho in [-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
        decomp_rho = eigen([[1,rho] [rho,1]]) ##use decomposability to calculate each separately
        for mu in [-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
            decomp_mu = eigen([[1,mu] [mu,1]]) ##use decomposability to calculate each separately
            factor = calculate_covariance_factor(decomp_mu, decomp_rho, decomp_hamming)
            targetname = "factor_" * string(nucleotides) * "_" * string(mu) * "_" * string(rho) * "_" * string(index) * ".csv"
            open(joinpath(dirname(@__FILE__), targetname), "w") do io
                writedlm(io, factor, ',') ##write out 
            end
        end
    end
#end