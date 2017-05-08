# ldr.r : lossless dimension reduction using QR
# timing comparisons

#------------------------------------------------------------------



ldr = function(x)
{
  # x is p * n, n points in p dimensions, p > n-1
  # find coordinates of the n points in a hyperplane in n-1 dimensions
  n = ncol(x)
  q = qr(x - x[,n])   # QR.  Origin = the last point, arbitrary
  r = qr.R(q)[,q$pivot]
  return(r[-n,])   # n points (columns) in n-1 dimensions (rows)
}

