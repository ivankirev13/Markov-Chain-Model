library(seqinr)
dat<-read.fasta("sequences.fasta")
genome = dat[[1]]
#calculate the number of each basis in the genome, then divide by the length of the whole sequence to get the proportion
table_1 = count(genome,1)/getLength(genome)
print(table_1)

#select the fist 100000 bases
n = genome[1:100000]

#for each window of 1000 bases calculate the propotion of G + C
proportion = c()
for (i in c(1:99000)){
  proportion = c(proportion,GC(n[i:(i+1000)]))
}

#plot these proportions agains position
plot(c(1:99000), proportion,xlab="Position", ylab = "Proportion of G+ C")

#put the count of bases for a, c, g, t in a vector q_0
q_0 = as.vector(table_1)

#calculate each $q_{ij}$
qij = c(q_0[1]*q_0[1],q_0[1]*q_0[2],q_0[1]*q_0[3],q_0[1]*q_0[4],q_0[2]*q_0[1],q_0[2]*q_0[2],q_0[2]*q_0[3],q_0[2]*q_0[4],q_0[3]*q_0[1],q_0[3]*q_0[2],q_0[3]*q_0[3],q_0[3]*q_0[4],q_0[4]*q_0[1],q_0[4]*q_0[2],q_0[4]*q_0[3],q_0[4]*q_0[4])
print(qij)

#calculate the number of each pair in the genome
table_2 = count(genome,2)
table_2

#calculate $n_{ij}$
nij = as.vector(table_2)

#calculate $p_{ij}$
pij = nij/(getLength(genome)-1)

#set initial value for the sum
log_likelihood = 0

#calculate the sum from above
for (i in 1:length(pij)){
  log_likelihood = log_likelihood + nij[i]*(log(pij[i]) - log(qij[i]))
}
print(log_likelihood)

#function that returns the log-likelihood ratio
likelihood_function = function(genome){
  
  #here q_0 and qij will be the same since they don't change on permtations of the initial genome
  q_0 = as.vector(count(genome,1)/getLength(genome))
  qij = c(q_0[1]*q_0[1], q_0[1]*q_0[2], q_0[1]*q_0[3], q_0[1]*q_0[4], q_0[2]*q_0[1], q_0[2]*q_0[2], q_0[2]*q_0[3], q_0[2]*q_0[4], q_0[3]*q_0[1], q_0[3]*q_0[2], q_0[3]*q_0[3], q_0[3]*q_0[4], q_0[4]*q_0[1], q_0[4]*q_0[2], q_0[4]*q_0[3], q_0[4]*q_0[4])
  
  #calculate nij and pij using the same method as before
  nij = as.vector(count(genome,2))
  pij = nij/(getLength(genome)-1)
  
  #calculate the log-likelihood ratio as before
  log_likelihood = 0
  for (i in 1:length(pij)){
    log_likelihood = log_likelihood + nij[i]*(log(pij[i]/qij[i]))}
  
#return log-likelihood ratio
return(log_likelihood)
  
#initialize vector in which we collect the ratios
ratios = c()
  
#calculate log-likelihood ratio for 100 permutations
for (i in 1:100){
  new_genome = permutation(genome)
  ratios = c(ratios,likelihood_function((new_genome)))
}
print(ratios)
  
#qqplot
qqplot(2*ratios, rchisq(100,df = 9))
abline(a = 0, b = 1, col = 'red', lty = 2)
  
#first calculate ni, we have repeated each sum four times since our vector needs to have 16 entries and each i takes up four of them.
ni = c(rep(sum(nij[1:4]),4), rep(sum(nij[5:8]),4), rep(sum(nij[9:12]),4), rep(sum(nij[13:16]),4))
  
#now we calculate and print the maximum likelihood estimates
phat = nij/ni
print(phat)
  
#create the transition matrix from the maximum likelihood data
P = matrix(phat,nrow = 4, byrow = TRUE)
print(P)
  
#create a function that returns a sample of size m
simulation = function(m,initial_value = 1){
  x = vector(length = m)
  x[1] = initial_value
  for (i in 2:m){
    x[i] = sample(x = 4, size = 1,prob = P[x[i-1],])
  }
    
  #since our result is a vector with numbers 1,2,3,4, we each number into letters from the genome
  x = plyr::mapvalues(x, from = c(1,2,3,4), to = c("a", "c", "g", "t"))
  return(x)
    
#create n different samples of size m and store them in a list
m = 100000
n = 100
realizations = list(rep(0),n)
for (i in 1:n){
  realizations[i] = list(simulation(m, initial_value = 1))
}


#calculate the maximum likelihood estimates for each sample
mle = vector(length = n)
for (i in 1:length(realizations)){
  newNij = as.vector(count(realizations[[i]],2))
  newNi = c(rep(sum(newNij[1:4]),4), rep(sum(newNij[5:8]),4), rep(sum(newNij[9:12]),4), rep(sum(newNij[13:16]),4))
  mle[i] = list(newNij/newNi)
}
  
#calculate the difference between the mle for each sample and the mle based on the DNA data
error = list(rep(c(0),n))
for (i in 1:n){
  error[i] = list(mle[[i]] - phat)
}
#create a vector containing all the errors from the 100 samples
errors = c()
for (i in 1:n){
  errors = c(errors, error[[i]])
}
#create a histogram of all the errors
hist(errors)  

#create a qqplot to show that the errors are approximately small and follow a normal distribution with standart deviation sd = 0.003

qqplot(rnorm(1600, sd = sd(errors)), errors)
abline(a = 0, b = 1, col = 'red', lty = 2)

correlation = c()
for (i in 1:15){
  for (j in (i + 1):16){
    vec1 = c()
    vec2 = c()
    for (k in 1:50){
      vec1 = c(vec1, mle[[k]][i])
      vec2 = c(vec2, mle[[k]][j])
    }
    correlation = c(correlation, cor(vec1,vec2))
  }
}
hist(correlation)

  
