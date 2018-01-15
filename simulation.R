
#---Load library
library(transport)
#change to local location
source('/functions.R')

#---Set parameters for one simulation
# p - mixture proportion
# mu - mean in first coordinate of second population
# dim - dimension of samples
# N - sample size of each population
# sample 1 is N(0,1); sample 2 is (1-p) N(0,1) + p N(mu,1), for mu = (mu,0,0,...,0)
p = 0
mu = 0
dim = 20
N = 100000

#---Generate samples 
x = generate_samples(p,mu,dim,N);
x_1 = x[[1]];
x_2 = x[[2]];

#---Random binning of combined sample, compute center of each bin and frequency in each bin, distance matrix of centers
result = random_bin(x_1,x_2);
locations = result$center;
freq = result$freq;

distance = distance_comp(locations)

#---compute distances
#EMD
A = wpp(locations,freq[,1])
B = wpp(locations,freq[,2])
emd = wasserstein(A, B, p=1, prob=FALSE)

#QF
m = max(distance)
qd <- sqrt(t(freq[,2]-freq[,1])%*%((1-distance/m)%*%(freq[,2]-freq[,1])))

#chi-square
chi2d <- sum((freq[,2]-freq[,1])^2 / (freq[,2]+freq[,1])) 

#---Print result
print(c(emd,qd,chi2d))

#---notification
library(beepr)
beep()
