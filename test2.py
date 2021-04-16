# coding:UTF-8
'''
Created on 2015年5月19日
@author: zhaozhiyong
'''

from numpy import *


# fun
def fun(x):
    return 10*x[0,0]**2+x[1,0]**2


# gfun
def gfun(x):
    result = zeros((2, 1))
    result[0, 0] = 20*x[0,0]
    result[1, 0] = 2*x[1,0]
    return result


from numpy import *


def dfp(fun, gfun, x0):
    result = []
    maxk = 500
    rho = 0.55
    sigma = 0.4
    m = shape(x0)[0]
    Hk = eye(m)
    k = 0
    while (k < maxk):
        gk = mat(gfun(x0))  # 计算梯度
        dk = -mat(Hk) * gk
        m = 0
        mk = 0
        while (m < 20):
            newf = fun(x0 + rho ** m * dk)
            oldf = fun(x0)
            if (newf < oldf + sigma * (rho ** m) * (gk.T * dk)[0, 0]):
                mk = m
                break
            m = m + 1

        # DFP校正
        x = x0 + rho ** mk * dk
        sk = x - x0
        yk = gfun(x) - gk
        if (sk.T * yk > 0):
            Hk = Hk - (Hk * yk * yk.T * Hk) / (yk.T * Hk * yk) + (sk * sk.T) / (sk.T * yk)

        k = k + 1
        x0 = x
        result.append(x0)

    return result

if __name__ == '__main__':
    x0 = mat([[0.1], [1.]])
    result = dfp(fun, gfun, x0)
    print(result)