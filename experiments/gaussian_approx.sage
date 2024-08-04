from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
from sage.probability.probability_distribution import RealDistribution

n = 50
sigma = 3.2
bound = 5 * sigma
repeat = 100000
  
C = RealDistribution('gaussian', sigma)
D = DiscreteGaussianDistributionIntegerSampler(sigma = sigma)

np_succ = 0.0

# Random orthogonal matrix, for coordinate change
B = MatrixSpace(RDF,n).random_element().QR()[0] 
Binv = B.inverse()

for _ in range(repeat):  
  v_std = vector(RDF, n)
  for i in range(n):
    v_std[i] = D()

  v_d = Binv * v_std
  
  v_d_abs = vector([abs(vv) for vv in v_d])
  if max(v_d_abs) < bound:
    np_succ += 1

print("p_np for discrete & rotate = ", np_succ / repeat)
p_np_theory = RR(erf(bound / (sqrt(2)*sigma)))**n
print("v.s. p_np for conti = ", p_np_theory)

U = RealDistribution('uniform', [-bound, bound])
admissibility_succ = 0.0

CG = RealDistribution('gaussian', sigma)

for _ in range(repeat):  
  v_std = vector(RDF, n)
  x = vector(RDF, n)
  for i in range(n):
    v_std[i] = CG.get_random_element()
    x[i] = U.get_random_element()  
  
  v_d = Binv * v_std

  check = v_d + x
  check_abs = vector([abs(vv) for vv in check])

  if max(check_abs) < bound:
    admissibility_succ += 1
  
X = RR(2*bound/(sqrt(2)*sigma))
p_adm_theory = RR(erf(X) + (exp(-X**2) - 1)/(sqrt(pi)*X))**n

print("p_adm for discrete & rotate = ", admissibility_succ / repeat)
print("v.s. p_adm for conti = ", p_adm_theory)

