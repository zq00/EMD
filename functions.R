
#-----------------------------
#generate_samples(p,mu,dim,N) : generate sample 1 from N(0,1) and sample 2 from (1-p) N(0,1) + p N(mu,1), for mu = (mu,0,0,...,0), each of size N 
#INPUT: p - mixture proportion
#       mu - mean in first coordinate of second population
#       dim - dimension of samples
#       N - sample size of each population
#OUTPUT: list of length 2, each member is a matrix of dimenstion N*dim corresponding to sample 1 and sample 2 

generate_samples = function(p,mu,dim,N){
  x_1 = matrix(rnorm(dim*N,0,1),N,dim);
  x_2 = matrix(rnorm(dim*N,0,1),N,dim);
  
  n = rbinom(1,N,p);
  x_2[(1:n),1] = x_2[(1:n),1]+mu;
  
  return(list(x_1, x_2));
}

#-----------------------------
#random_bin(x_1, x_2) : performs random binning on combined sample
#INPUT: x_1, x_2 - (matrix) of dimension N*dim, corresponds to sample 1 and sample 2
#OUTPUT: center - (matrix) of dimension nbin*dim, center of mass of each bin, computed on the combined sample
#        freq - (matrix) of dimension nbin*2, frequency in each bin for each sample

random_bin = function(x_1, x_2){
  level = 0;
  nbin = 2^level;
  
  N_1 = dim(x_1)[1];
  N_2 = dim(x_2)[1];
  N = N_1 + N_2;
  
  dim = dim(x_1)[2];
  
  x = rbind(x_1 , x_2);
  label = matrix(1,N,1);
  
  while (nbin <= (N/(2*log(N)))){
    for(i in 1:nbin){
      s = which(label==i);
      
      splitdim = which.max(apply(x[s,],2,var)); 
      
      newsplit = quantile(x[s,splitdim],0.5);
      
      label[s[which(x[s,splitdim]>newsplit)]] = nbin + i;
    }
    
    level = level + 1;
    nbin = 2^level;
  }
  
  center = matrix(0,nbin,dim);
  freq_1 = matrix(0,nbin,1);
  freq_2 = matrix(0,nbin,1);
  
  for(i in 1:nbin){
    s = label==i;
    
    center[i,]=colMeans(x[s,]);
    
    freq_1[i] = sum(label[1:N_1] == i) / N_1;
    freq_2[i] = sum(label[(N_1+1):N] == i)/N_2;
  }
  
  result = list(center,cbind(freq_1,freq_2));
  names(result) = c("center","freq");

  return(result);
}

#-----------------------------
#distance_comp(center) : computes distance matrix
#INPUT: center - (matrix) of dimension nbin*dim, center of mass of each bin
#OUTPUT: distance - (matrix) of dimension nbin*nbin, Euclidean distance between each pair of centers of each bin

distance_comp = function(center){
  nbin = dim(center)[1];
  
  distance = matrix(0,nbin,nbin);
  for(i in 1:nbin){
    for(j in 1:nbin){
      #2-norm
      distance[i,j] = sqrt(sum((center[i,]-center[j,])^2))
    }
  }
  
  return(distance)
}











