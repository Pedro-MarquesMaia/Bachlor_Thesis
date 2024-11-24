import numpy as np
import matplotlib.pyplot as pl

N = 3

x = np.zeros(shape=(N,3),dtype=np.int64)
v = 0*x

x[0][0],x[0][1] = -970004360000, 243087530000
x[1] = 0*x[0]
x[2] = -x[0]

v[0][0],v[0][1] = 466203685, 432365730
v[1] = -2*v[0]
v[2] = 1*v[0]

hdt = 50
for s in range(20):
    x += v*hdt
    acc = 0*v
    for i in range(N):
        for j in range(i):
            dr = x[i] - x[j]
            dr = 1.0*dr  # otherwise may overflow
            denom = np.sum(dr**2)**(3/2)
            eqop = dr/denom
            eqop *= 1e30
            eqop = (np.rint(eqop)).astype(int)
            acc[i] -= eqop
            acc[j] += eqop
    v += acc*2*hdt
    x += v*hdt
    pl.plot(x[0][0],x[0][1],'.',color='red')
    pl.plot(x[1][0],x[1][1],'.',color='green')
    pl.plot(x[2][0],x[2][1],'.',color='blue')
    vs = 0*v[0]
    lz = 0
    for i in range(N):
        vs += v[i]
        lz += x[i][0]*float(v[i][1]) - x[i][1]*float(v[i][0])
    print(vs,'%e' % lz)

pl.gca().set_aspect(1)
pl.show()
