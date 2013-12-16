from numpy import *

__author__ = 'pheodor'


class CmCmGenerator():
    def __init__(self, N, x0, y0, nTrans=1000):
        self.N = N
        self.x0 = x0
        self.y0 = y0
        self.coefficients = []
        self.nTrans = nTrans

    def generate(self, kappa):
        n = self.N + self.nTrans
        xx = empty(n)
        yy = empty(n)
        xx[0] = self.x0
        yy[0] = self.y0
        twoPi = 2 * pi
        for i in xrange(1, n, 1):
            xx[i] = (xx[i - 1] + self.coefficients[0] * twoPi + self.coefficients[1] * sin(xx[i - 1])) % twoPi + kappa * yy[i - 1]
            yy[i] = (yy[i - 1] + self.coefficients[2] * twoPi + self.coefficients[3] * sin(yy[i - 1])) % twoPi
        radX = xx[self.nTrans:]
        radY = yy[self.nTrans:]
        return radX, radY
