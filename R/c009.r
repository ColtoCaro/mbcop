# c009.r : 
# simple bootstrap

#*********************************************************************************;

mycluster1 = function(x, ignored1,  ignored2) { 
  # random clustering into 2 clusters
  # x : n by 1 vector
  n1 = length(x) / 2;
  n2 = length(x) - n1;
  cl = c(rep(1, n1), rep(2, n2));
  return(list(cluster=cl));
}

#*********************************************************************************;

mycluster2 = function(x, ignored1,  ignored2) {
  # below/above median
  # x : n by 1 vector
  cl = 1 + (x > median(x));
  return(list(cluster=cl));
}

#*********************************************************************************;

ch2 = function(n) {n * (n - 1) / 2}

#*********************************************************************************;

smry2x2 = function(a, g) {
  # a : true cluster id
  # g : assigned cluster id
  # return the 2x2 table of pairs (as a vector)
  n   = table(a, g)
  mdd = ch2(sum(n))  # total #pairs
  m22 = sum(ch2(n))  # #links by both
  m2d = sum(ch2(apply(n, 1, sum)))  #true links
  md2 = sum(ch2(apply(n, 2, sum)))  #assigned links
  m12 = md2 - m22
  m21 = m2d - m22
  m11 = mdd - m2d - md2 + m22
  return(c(m11, m12, m21, m22))
}

#*********************************************************************************;

get_indices = function(freq)
{
  sens = freq[,4] / (freq[,3] + freq[,4])
  spec = freq[,1] / (freq[,1] + freq[,2])
  return(list(sens=sens, spec=spec))
}

#*********************************************************************************;

get_kappa = function(sens, spec, prev)
{
  p11 = spec * (1 - prev)
  p12 = 1 - prev - p11
  p22 = sens * prev
  p21 = prev - p22
  e11 = (1 - prev) * (p11 + p21)
  e22 = prev * (p12 + p22)
  kappa = 2 * (p11 - e11) / (1 - (e11 + e22)) # (O-E)/(1-E)
  return(kappa) 
}
#*********************************************************************************;

const_n = 100;
const_G = 2;
const_delta = 2;
const_B = 1e3;
const_m = const_n / const_G;
const_mu = c(0, const_delta);

true_cluster = c(rep(1, const_m), rep(2, const_m));
x = rnorm(const_n, const_mu[true_cluster]);
kmx = kmeans(x, 2, nstart = 10);
freq = matrix(smry2x2(true_cluster, kmx$cluster), ncol=4);
get_indices (freq);
# for delta=(1,2,3,4) and the optimal classifier: sens=spec = (0.573, 0.733, 0.875, 0.956)


# Bootstrap sampling;
a = matrix(NA, const_n, const_B); # each column is a clustering (a partition);
for (j in 1:const_B) {
  s = sample(const_n, rep = T);         # bootstrap sample;
  kms    = kmeans(x[s], 2, nstart = 10);  
  a[,j][s] = kms$cluster;   # credit assignement to the original data points;
}

# look at pairs of data points, rows of a[,];
k = 0;
linked = rep(NA, const_n * (const_n - 1) / 2);
row1id = linked;
row2id = linked;
for (row1 in 1:(const_n - 1)) {
  for (row2 in (row1 + 1):const_n) {
     k = k + 1;
     linked[k] = mean(a[row1,] == a[row2,], na.rm=T);  # prop(linked);
     row1id[k] = row1;
     row2id[k] = row2;
  }
}

summary(linked);
hist(linked);   # look for peaks near 0 and 1;

# look for problematic data points;
q = (abs(linked - 0.5) < 0.3);
table(q);
sort(table(c(row1id[q],row2id[q])), decreasing = T);
 
      
# compare to the true partition ;
row1class = true_cluster [row1id];
row2class = true_cluster [row2id];
linked_same = linked[row1class == row2class];
linked_different = linked[row1class != row2class];
summary(linked_same);  # sens?
summary(1 - linked_different);  # spec?

    
# compare to original estimated partition;
row1class = kmx$cluster [row1id];
row2class = kmx$cluster [row2id];
linked_same     = linked[row1class == row2class] ;
linked_different = linked[row1class != row2class] ;
summary(linked_same) ;  # sens?
summary(1 - linked_different);   # spec?


