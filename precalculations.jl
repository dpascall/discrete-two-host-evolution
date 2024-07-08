using Random, Distributions, StringDistances, LinearAlgebra, CSV, Tables

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

seqs = generate_seqs(4)
sim_mat = hamming_similarity_mat(4, seqs)

CSV.write(file = "8HammingDistance.csv", Tables.table(sim_mat))

for rho in [-1.0, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1] ##perform SVDs for each of the conditions
    for i in [1.0, 5.0, 10.0]
        decomp = svd(calculate_covarience_matrix(8, sim_mat, rho, i), full = true, alg = LinearAlgebra.QRIteration())
        CSV.write(file = "factor_" * string(rho) * "_" * string(index) * ".csv", Tables.table(decomp.U  * sqrt(diagm(decomp.S))))
    end
end
