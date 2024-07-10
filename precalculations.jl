#!/usr/bin/env julia
##takes length of nucleotide string as an argument
using Random, Distributions, StringDistances, DelimitedFiles, LinearAlgebra, MPI, MPIClusterManagers, Distributed

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
    return kron(([[1, rho] [rho, 1]]), (hamming_mat./nucleotides).^index) ##kronecker of seq correlation matrix and cross landscape correlation matrix
end

man = MPIManager(np = 7)
addprocs(man)

nucleotides = parse(Int, ARGS[1])

seqs = generate_seqs(nucleotides)
sim_mat = hamming_similarity_mat(nucleotides, seqs)

name = string(nucleotides) * "HammingDistance.csv"

open(joinpath(dirname(@__FILE__), name), "w") do io
    writedlm(io, sim_mat, ',') ##write out hamming mat
end

covariance = Elemental.Matrix(Float64)

for rho in [-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1] ##perform SVDs for each of the conditions
    for i in [1.0, 5.0, 10.0]
        covariance = calculate_covarience_matrix(nucleotides, sim_mat, rho, i)
        @mpi_do man decomp = svd(covariance, full = true, alg = LinearAlgebra.QRIteration())
        factor = decomp.U * diagm(sqrt.(decomp.S))
        targetname = "factor_" * string(rho) * "_" * string(index) * ".csv"
        open(joinpath(dirname(@__FILE__), targetname), "w") do io
            writedlm(io, factor, ',') ##write out 
        end
    end
end