#File containing internal helpful functions


ldr = function(x)
{
  # x is p * n, n points in p dimensions, p > n-1
  # find coordinates of the n points in a hyperplane in n-1 dimensions
  n = ncol(x)
  q = qr(x - x[,n])   # QR.  Origin = the last point, arbitrary
  r = qr.R(q)[,q$pivot]
  return(r[-n,])   # n points (columns) in n-1 dimensions (rows)
}


getDist <- function(x1, x2, mergeTable, distVec){
  #Function to find the distance at which two points (x1, x2) are merged
  i <- 1
  loc1 <- -1*x1
  loc2 <- -1*x2
  while(loc1 != loc2){
    if(mergeTable[i, 1] == loc1 | mergeTable[i, 2] == loc1){loc1 <- i}
    if(mergeTable[i, 1] == loc2 | mergeTable[i, 2] == loc2){loc2 <- i}
    if(loc1 == loc2){return(distVec[i])}
    i <- i + 1
  }
}

ch2 = function(n) {n * (n - 1) / 2}

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
