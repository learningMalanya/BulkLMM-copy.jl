##################################################################
# Fast linear mixed models
##################################################################
#
# We implement linear mixed models for data which has covariance of
# the form tau2*K + sigma2*I, where sigma2 and tau2 are positive
# scalars, K is a symmetric positive definite "kinship" matrix and I
# is the identity matrix.
#
# Models of this type are used for GWAS and QTL mapping in structured
# populations.
#
# ################################################################

##################################################################
# rotateData: rotate by orthogonal transformation
##################################################################
"""
rotateData: Rotates data with respect to the kinship matrix

y = phenotype matrix
X = predictor matrix
K = kinship matrix, expected to be symmetric and positive definite
"""
function rotateData(y::AbstractArray{Float64,2},X::AbstractArray{Float64,2},
                    K::Array{Float64,2})

    # check dimensions
    n = size(y,1)
    if( ( size(X,1) != n ) | ( size(K,1) != n ))
        error("Dimension mismatch.")
    end

#    # check symmetry and positive definiteness of K
#    if( !(issym(K)) )
#        error("K is not symmetric.")
#    end

#    if( !(isposdef(K)) )
#        error("K is not positive definite.")
#    end

    # spectral decomposition of a symmetric matrix
    EF = eigen(K)

    # return rotated phenotype, covariates, and eigenvalues
    return EF.vectors'y, EF.vectors'X, EF.values

end

"""
rotateData: Rotates data with respect to the kinship matrix

y = phenotype matrix
X = predictor matrix
K = kinship matrix, expected to be symmetric and positive definite
weights = vector of sample sizes (or weights inversely proportional to
    error variance)
"""

function rotateData(y::AbstractArray{Float64,2},X::AbstractArray{Float64,2},
                    K::Array{Float64,2}, weights::Array{Float64,1})

    # check dimensions
    n = size(y,1)
    if( ( size(X,1) != n ) | ( size(K,1) != n ))
        error("Dimension mismatch.")
    end

    # make vector of square root of the sample sizes
    w = sqrt.(n)
    # transform kinship
    scale!(K, w)
    scale!(w, K)
    # transform phenotype
    scale!(w, y)
    # transform predictor
    scale!(w, X)

    # pass to old function
    return rotateData(y,X,K)

end

function makeweights(h2::Float64,lambda::Array{Float64,1})
    return h2*lambda .+ (1.0-h2)
end


##################################################################
# function to fit linear mixed model by optimizing heritability
##################################################################
mutable struct LMM_estimates
    b::Array{Float64,2}
    sigma2::Float64
    h2::Float64
    ell::Float64
end


"""
fitlmm: fit linear mixed model using grid of values

y: 2-d array of (rotated) phenotypes
X: 2-d array of (rotated) covariates
lambda: 1-d array of eigenvalues
ngrid: number of grid values to consider
"""
function fitlmm(y::Array{Float64,2}, X::Array{Float64,2}, lambda::Array{Float64,1}, ngrid::Int64; 
                reml::Bool = false, loglik::Bool = true)
    h2vec = convert(Vector, (0:ngrid)/ngrid)

    function logLik0(h2::Float64)
        out = wls(y, X, 1.0./makeweights(h2, lambda); reml = reml, loglik = loglik)
        return -out.ell
    end

    ell = map(logLik0,h2vec)
    h2 = h2vec[argmin(ell)]
    est = wls(y, X, 1.0./makeweights(h2, lambda); reml = reml, loglik = loglik)

    return LMM_estimates(est.b, est.sigma2, h2, est.ell)
end



"""
fitlmm: fit linear mixed model

y: 2-d array of (rotated) phenotypes
X: 2-d array of (rotated) covariates
lambda: 1-d array of eigenvalues
reml: boolean indicating ML or REML estimation

"""
function fitlmm(y::Array{Float64, 2}, X::Array{Float64, 2}, lambda::Array{Float64, 1};
                reml::Bool = false, loglik::Bool = true, h20::Float64 = 0.5, d::Float64 = 1.0)
    function logLik0(h2::Float64)
        out = wls(y, X, 1.0./makeweights(h2, lambda); reml = reml, loglik = loglik)
        return -out.ell
    end

    opt = optimize(logLik0, max(h20-d, 0.0), min(h20+d, 1.0))
    h2 = opt.minimizer
    est = wls(y, X, 1.0./makeweights(h2, lambda); reml = reml, loglik = loglik)
    return LMM_estimates(est.b, est.sigma2, h2, est.ell)
end
