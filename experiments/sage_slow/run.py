from basic import *

n = 20; m = 20; stddev = 3.2; q = 1024; h = 8

def run_test(n, stddev, q, h, m, w, l, iteration = 1):
  GSnorm = [q] * m

  # Set Parameters
  R = ambiguity(n, h, w)
  r, p_adm = constraint_area(R, stddev, l, GSnorm, unif=False)

  print("------- Params -------")
  print("- # of reps (Ambiguity):", R)
  print("- Projection dim:", r)
  print("- log(p_admissible):", log(p_adm, 2))
  print("- Expected reps (R * p_adm):", R * p_adm)

  print("------- Actual Tests -------")
  # Actual List Construction & Collision Finding
  for _ in range(iteration):
    M, s = gen_instance(n-1, stddev, q, h-1, m=m)
    L = build_list(M, w, r, l, q)
    sol = near_collision(L, stddev, q, h)
    e = vector(ZZ, M*s)
    e = [modq(e[i], q) for i in range(len(e))]
    print("- |L_0| =", len(sol), "/ |L_1| =", len(L), "/ e[-r:] =", e[-r:])