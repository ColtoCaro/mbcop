#File containing internal helpful functions

#-----------------------Functions used in creating HROC plots
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

#--------------------------Functions used for creating pairwise indices

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


#-----------------------Functions used for creating the triplet indices
# triplib.r :  Compare partitions of triplets and pairs
#              Ordered and unordered versions
#              Fast version of triplet counting, O(m^3),
#                 m := #sets(p1) * #sets(p2)
#              This version should be faster than the O(n^3) method
#                whenever the average cluster size > #(clusters),
#                i.e. #(clusters) < sqrt(n)
#              The slower version runs in O(n^3)
#              Function library only
#--------------------------------------------------------

#========================================================

u2o_triplets <- function(x) {
  # Convert counts from unordered to ordered
  # x is a 5*5 matrix tabulating unordered triplets
  # Returns the 5*5 matrix tabulating the ordered triplets
  # Sum(matrix returned) = 6 * sum(x)

  if (!is.matrix(x)) return(NA)
  d = dim(x)
  if (d[1] != 5 | d[2] != 5) return(NA)
  y = matrix(NA, 5, 5)
  diag(y) = 6 * diag(x)
  y[1,1] = 6 * x[1,1]
  y[1,5] = 6 * x[1,5]
  y[5,1] = 6 * x[5,1]
  y[5,5] = 6 * x[5,5]
  y[2,2] = y[3,3] = y[4,4] =
    2 * (x[2,2] + x[3,3] + x[4,4])
  y[2,3] = y[2,4] = y[3,2] = y[3,4] = y[4,2] = y[4,3] =
    x[2,3] + x[2,4] + x[3,2] + x[3,4] + x[4,2] + x[4,3]
  y[1,2:4] = 2 * sum(x[1, 2:4])
  y[5,2:4] = 2 * sum(x[5, 2:4])
  y[2:4,1] = 2 * sum(x[2:4,1])
  y[2:4,5] = 2 * sum(x[2:4,5])
  dimnames(y) = dimnames(x)
  return(y)
}

#--------------------------------------------------------
xtab_otriplets <- function(p1, p2) {
  # Returns: 5*5 Matrix of xtab frequencies of how ORDERED triplets are
  # partitioned in p1 vs in p2
  # Counts ordered triplets
  x = xtab_utriplets (p1, p2)
  y = u2o_triplets(x)
  return(y)
}
#--------------------------------------------------------

xtab_utriplets <- function(p1, p2) {
  # Fastest version.
  # Returns 5*5 Matrix of xtab frequencies of how UNORDERED triplets are
  # partitioned in p1 vs in p2
  # Works by processing a crosstab of p1*p2
  # O(K1 K2)

  n  = length(p1)
  if (n < 3) return(NA)
  if (n != length(p2)) return(NA)

  f   = as.matrix(table(p1, p2))
  nrows = nrow(f)
  ncols = ncol(f)
  rowid = rep(1:nrows, times=ncols) # row number
  colid = rep(1:ncols, each =nrows) # column number
  rt    = apply(f, 1, sum)  # row totals
  ct    = apply(f, 2, sum)  # col totals
  rtpc  = rt[rowid]         # row total per cell
  ctpc  = ct[colid]         # col total per cell
  f     = as.vector(f)

  A1 = choose(rt, 2) %*% (n - rt)
  A3 = sum(choose(rt, 3))
  A0   = choose(n, 3) - A1 - A3
  B1 = choose(ct, 2) %*% (n - ct)
  B3 = sum(choose(ct, 3))
  fc2   = choose(f, 2)
  A3B1 =  fc2 %*% (rtpc - f)
  A1B3 =  fc2 %*% (ctpc - f)
  A1B1d=  fc2 %*% (n - rtpc - ctpc + f)    # diag
  A1B1o=  sum(f * (rtpc - f) * (ctpc - f)) # off-diag
  A1B1 =  A1B1d + A1B1o
  A3B3 =  sum(choose(f, 3))
  A0B1 = B1 - A1B1 - A3B1
  A0B3 = B3 - A1B3 - A3B3
  A1B0 = A1 - A1B1 - A1B3
  A3B0 = A3 - A3B1 - A3B3
  A0B0 = A0 - A0B1 - A0B3

  ret = c(A0B0, A1B0,  0, 0, A3B0,
          A0B1, A1B1d, 0, 0, A3B1,
          0   , A1B1o, 0, 0, 0,
          rep(0, 5),
          A0B3, A1B3,  0, 0, A3B3)
  ret = matrix(ret, nrow=5)

  return(ret)
}

#--------------------------------------------------------

xtab_utriplets2 <- function(p1, p2) {
  # Fast version.
  # Returns 5*5 Matrix of xtab frequencies of how UNORDERED triplets are
  # partitioned in p1 vs in p2
  # Works by processing a crosstab of p1*p2
  # A triplet can belong to:
  #    * one cell
  #    * two different cells
  #    * three different cells
  # The three cases are processed separately
  # O(K1^3 K2^3)

  coder <- function(y1, y2, y3) {
    z = (y1 == y2) + 2 * (y1 == y3) + 4 * (y2 == y3)
    factor(z, levels = c(0,1,2,4,7))
  }

  if (length(p1) != length(p2)) return(NA)
  if (length(p1) < 3) return(NA)

  ret = matrix(0, 5, 5)
  f   = as.matrix(table(p1, p2))
  nrows = nrow(f)
  ncols = ncol(f)
  rowid = rep(1:nrows, times=ncols) # row number
  colid = rep(1:ncols, each =nrows) # column number
  f     = as.vector(f)

  # Case 1:
  # All members of the triplet come from one cell
  code1 = coder(1, 1, 1)
  ret[code1, code1] = sum(choose(f, 3))

  # Case 2:
  # Members of the triplet come from two different cells
  if (length(f) >= 2) {
    cell_pairs = combn(length(f), 2)
    row1 = rowid[cell_pairs[1,]]
    col1 = colid[cell_pairs[1,]]
    n1   =     f[cell_pairs[1,]]
    row2 = rowid[cell_pairs[2,]]
    col2 = colid[cell_pairs[2,]]
    n2   =     f[cell_pairs[2,]]
    # The two cells are (row1,col1) and (row2,col2)
    # Cell counts are n1 and n2
    # Two obs from cell 1
    code1 = coder(row1, row1, row2)
    code2 = coder(col1, col1, col2)
    m = choose(n1, 2) * n2
    a = tapply(m, list(code1,code2), sum)
    a = ifelse(is.na(a), 0, a)
    ret = ret + a

    # Two obs from cell 2
    code1 = coder(row1, row2, row2)
    code2 = coder(col1, col2, col2)
    m = n1 * choose(n2, 2)
    a = tapply(m, list(code1,code2), sum)
    a = ifelse(is.na(a), 0, a)
    ret = ret + a
  }

  # Case 3:
  # Members of the triplet come from three different cells
  if (length(f) >= 3) {
    cell_trips = combn(length(f), 3)
    row1 = rowid[cell_trips[1,]]
    col1 = colid[cell_trips[1,]]
    n1   =     f[cell_trips[1,]]
    row2 = rowid[cell_trips[2,]]
    col2 = colid[cell_trips[2,]]
    n2   =     f[cell_trips[2,]]
    row3 = rowid[cell_trips[3,]]
    col3 = colid[cell_trips[3,]]
    n3   =     f[cell_trips[3,]]
    # The three cells are (row1,col1), (row2,col2) and (row3,col3)
    # Cell counts are n1, n2 and n3
    code1 = coder(row1, row2, row3)
    code2 = coder(col1, col2, col3)
    m = n1 * n2 * n3
    a = tapply(m, list(code1, code2), sum)
    a = ifelse(is.na(a), 0, a)
    ret = ret + a
  }

  return(ret)
}

#--------------------------------------------------------
xtab_otriplets2 <- function(p1, p2) {
  x = xtab_utriplets2 (p1, p2)
  y = u2o_triplets(x)
  return(y)
}
#--------------------------------------------------------

xtab_otriplets0 <- function(p1, p2) {
  # Matrix of xtab frequencies of how ORDERED triplets are
  # partitioned in p1 vs in p2
  # SUPER SLOW, use only for checking the other functions
  # SUPER SLOW, do NOT use in production code
  # O(n^3)

  coder <- function(y1, y2, y3) {
    z = (y1 == y2) + 2*(y1 == y3) + 4*(y2 == y3)
  }

  if (length(p1) != length(p2)) return(NA)
  if (length(p1) < 3) return(NA)
  n = length(p1)
  a1 = integer(n*(n-1)*(n-2))
  a2 = integer(n*(n-1)*(n-2))
  ijk = 0
  for (i in 1:n) {
    for (j in 1:n) {
      for (k in 1:n) {
        if (i != j & i != k & j != k) {
          ijk = ijk + 1
          a1[ijk] = coder(p1[i], p1[j], p1[k])
          a2[ijk] = coder(p2[i], p2[j], p2[k])
        }}}}
  a1 = factor(a1, levels = c(0,1,2,4,7))
  a2 = factor(a2, levels = c(0,1,2,4,7))
  return(table(a1, a2))
}

#--------------------------------------------------------

xtab_utriplets1 <- function(p1, p2) {
  # 5*5 Matrix of xtab frequencies of how UNORDERED triplets are
  # partitioned in p1 vs in p2
  # This is a slow version, O(n^3)
  # O(n^3)

  coder <- function(y) {
    z = (y[1,]==y[2,]) + 2*(y[1,]==y[3,]) + 4*(y[2,]==y[3,])
    factor(z, levels = c(0,1,2,4,7))
  }
  if (length(p1) != length(p2)) return(NA)
  if (length(p1) < 3) return(NA)
  a1 = coder(combn(p1, 3))
  a2 = coder(combn(p2, 3))
  return(table(a1, a2))
}

#--------------------------------------------------------
xtab_otriplets1 <- function(p1, p2) {
  x = xtab_utriplets1 (p1, p2)
  y = u2o_triplets(x)
  return(y)
}
#========================================================

xtab_upairs <- function(p1, p2) {
  # Returns 2*2 Matrix of xtab frequencies of how UNORDERED pairs are
  # partitioned in p1 vs in p2
  # Works by processing a crosstab of p1*p2
  # O(K1 K2)

  if (length(p1) != length(p2)) return(NA)
  if (length(p1) < 2) return(NA)

  ret = matrix(0, 2, 2)
  colnames(ret) = rownames(ret) = c(0, 1)
  f   = as.matrix(table(p1, p2))
  rt  = apply(f, 1, sum)
  ct  = apply(f, 2, sum)
  ret[2, 2] = sum(choose(f, 2))
  ret[2, 1] = sum(choose(rt, 2)) - ret[2, 2]
  ret[1, 2] = sum(choose(ct, 2)) - ret[2, 2]
  ret[1, 1] = choose(length(p1), 2) - sum(ret)

  return(ret)
}

#--------------------------------------------------------
xtab_opairs <- function(p1,p2) {
  # Matrix of xtab frequencies of how ORDERED pairs are
  # partitioned in p1 vs in p2
  # Counts ordered pairs

  return(2 * xtab_upairs(p1,p2))
}

#========================================================

kappa_stat <- function(x) {
  # x is a square matrix
  if (!is.matrix(x)) return(NA)
  d = dim(x)
  if(d[1] != d[2]) return(NA)
  n = sum(x)
  if (n <= 0) return(NA)
  x = x / n
  rt = apply(x, 1, sum)
  ct = apply(x, 2, sum)
  po = sum(diag(x))
  pe = sum(diag(rt %*% t(ct)))
  num = po - pe
  den = 1 - pe
  if (num == 0 & den == 0) return(1)
  return(num / den)
}

#--------------------------------------------------------

four_kappas <- function(x) {
  # x is a 3*3 matrix, rows and cols correspond to 0,1,3

  if (!is.matrix(x)) return(NA)
  d = dim(x)
  if(d[1] != 3 | d[2] != 3) return(NA)

  k013 = kappa_stat(x)
  k01 = kappa_stat(x[1:2,1:2])
  k13 = kappa_stat(x[2:3,2:3])
  k03 = kappa_stat(x[c(1,3), c(1,3)])
  return(c(k01,k03,k13,k013))
  # return(list(k01=k01,k03=k03,k13=k13,k013=k013))
}
#--------------------------------------------------------

five_kappas <- function(x) {
  # x is a 5*5 matrix, rows and cols correspond to 0,1,1,1,3
  # This relies on the special structure of x

  if (!is.matrix(x)) return(NA)
  d = dim(x)
  if(d[1] != 5 | d[2] != 5) return(NA)

  i = x[2,2]
  j = x[2,3]
  k1 = (i - j) / (i + j + j) # kappa in the inner 3*3
  y = cbind(x[,1], x[,2]+x[,3]+x[,4], x[,5])
  z = as.matrix(rbind(y[1,], y[2,]+y[3,]+y[4,], y[5,]))

  return(c(k1, four_kappas(z)))
}

#--------------------------------------------------------

comp_part <- function(p1, p2) {
  # Compare two partitions
  # For triplets, display ordered and unordered

  print("Ordered pairs:")
  xt2o = xtab_opairs(p1,p2)
  print(xt2o)

  print("Unordered Triplets:")
  xt3u = xtab_utriplets(p1,p2)
  print(xt3u)

  if (!is.matrix(xt3u)) return()

  print("Ordered Triplets:")
  xt3o = u2o_triplets(xt3u) #same as: xt3o = xtab_otriplets(p1,p2)
  print(xt3o)

  print("Kappa: pairs, ordered triplets:")
  print(c(kappa_stat(xt2o), kappa_stat(xt3o) ))

  print("Margins of ordered triplets:")
  m =  sum(xt3o)
  print("Row totals:")
  print(apply(xt3o, 1, sum) / m)
  print("Column totals:")
  print(apply(xt3o, 2, sum) / m)

  # Inner 3x3
  print("Inner 3x3, diagonal, kappa:")
  i33 = sum(xt3o[2:4,2:4])  # inner 3x3
  d33 = sum(diag(xt3o)[2:4]) / i33  # diagonal fraction
  k33 = 1.5 * (d33 - 1/3)
  print(c(i33 / m, d33, k33))

  print("Five_kappas (triplets): k1, k01, k03, k13, k013:")
  print(five_kappas(xt3o))

}

#--------------------------------------------------------

