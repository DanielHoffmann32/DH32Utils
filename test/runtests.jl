using DH32Utils
using Base.Test
using GZip

# write your own tests here
x = vec(readcsv(gzopen("testdist.csv.gz")))
w = vec(readcsv(gzopen("testweights.csv.gz")))
q1,q2 = block_average_spread(x, dimin=1000, ddi=1000, weighted=true, w=w)
@test abs(q1[1]-1.7432) < 0.001
@test abs(q2[1]-2.0877) < 0.001
