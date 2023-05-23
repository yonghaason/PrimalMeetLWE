from itertools import combinations


hw = 10  
dim = 16
r1 = 2*(1.01)^(dim-1)
delta = 1.01
l = 0.5


## R matrix construct
R = matrix(QQ,dim)
R[0,0] = r1
for i in range(dim-1):
	R[1+i,1+i] = R[i,i]/delta
for i in range(dim-1):
	for j in range(i+1,dim):
		R[i,j] = (random()-1/2) * R[i,i]


## Volume ratio
ratio = RR((2*l)^dim/ prod([R[i,i] for i in range(dim)]))


## Check
Card = binomial(dim,hw)
Card2 = binomial(hw,hw/2)
RR(Card*Card2/prod([R[i,i] for i in range(dim)]))


RR(log(Card*Card2,2))


## Mod operation
def Modmatrix(v, R):
	for i in range(dim):
		v -= round(v[dim-1-i]/R[dim-1-i,dim-1-i])*R.transpose()[dim-1-i]
	return v

## matrix M setup
M = matrix(QQ,dim)
for i in range(dim):
	for j in range(dim):
		M[i,j] = (random()-1/2) * R[i,i]

## low hamming weight vector setup
S = []
for i in range(dim):
	S += [i+1]

subsets = list(combinations(S,hw))
Card = binomial(dim,hw)
Card2 = binomial(hw,hw/2)

## Count the number of Ms, which is lying on the cube.

Test = 0
for i in range(Card):
	listpm1 = subsets[i]
	subsubsets = list(combinations(listpm1,hw/2))
	for j in range(Card2):
		temps = vector(ZZ,dim)
		listp1 = subsubsets[j]
		listm1 = list(set(listpm1) - set(listp1))
		for k1 in range(hw/2):
			temps[listp1[k1]-1]=1
			temps[listm1[k1]-1]=-1
		tempevaluated = Modmatrix(M*temps, R)
		if vector(tempevaluated).norm(infinity) < l:
			Test += 1


print(RR(Test/Card*Card2))
print(RR(ratio))