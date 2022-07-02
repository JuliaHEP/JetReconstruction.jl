# try using LorentzVectorHEP directly here

@inline energy(p) = p[1]

@inline px(p) = p[2]

@inline py(p) = p[3]

@inline pz(p) = p[4]

@inline pt(p) = @fastmath sqrt(p[2]^2 + p[3]^2)
kt = pt

@inline phi(p) = @fastmath atan(p[3], p[2])
ϕ = phi

@inline mass(p) = @fastmath sqrt(p[0]^2 - p[2]^2 - p[3]^2 - p[4]^2)

@inline eta(p) = atanh(sqrt(p[2]^2 + p[3]^2 + p[4]^2)/p[1]) # WARNING: possibly incorrect
η = eta
