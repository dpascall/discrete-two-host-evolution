#!/usr/bin/env julia
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

function calculate_covarience_matrix(nucleotides::Int, hamming_mat::Matrix{Int64}, rho::Float64, index::Float64)
    return kron((hamming_mat./nucleotides).^index, ([[1, rho] [rho, 1]])) ##kronecker of seq correlation matrix and cross landscape correlation matrix
end

seqs = generate_seqs(8)
sim_mat = hamming_similarity_mat(8, seqs)

open(joinpath(dirname(@__FILE__), "8HammingDistance.csv"), "w") do io
    writedlm(io, sim_mat, ',') ##write out hamming mat
end

Threads.@threads for rho in [-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1] ##perform SVDs for each of the conditions
    for i in [1.0, 5.0, 10.0]
        decomp = svd(calculate_covarience_matrix(8, sim_mat, rho, i), full = true, alg = LinearAlgebra.QRIteration())
        targetname = "factor_" * string(rho) * "_" * string(index) * ".csv"
        open(joinpath(dirname(@__FILE__), targetname), "w") do io
            writedlm(io, decomp, ',') ##write out 
        end
    end
end