

import numpy as np
def golden_section(f, section, delta,flag=1):
    '''
    黄金分割法
    :param f: 单谷函数
    :param section: 初始区间
    :param delta: 区间精度
    :return: 目标函数的最小值
    '''
    t = 0.618  # 黄金分割比
    # 当前区间
    alpha = section[0]
    beta = section[1]

    lam = alpha + (1 - t) * (beta - alpha)
    mu = alpha + t * (beta - alpha)

    alpha_list = [alpha]
    beta_list = [beta]
    lam_list = [lam]
    mu_list = [mu]

    while (beta - alpha >= delta):
        if f(lam) - f(mu) > 0:
            alpha = lam
            lam = mu
            mu = alpha + t * (beta - alpha)
        else:
            beta = mu
            mu = lam
            lam = alpha + (1 - t) * (beta - alpha)
        alpha_list.append(alpha)
        beta_list.append(beta)
        lam_list.append(lam)
        mu_list.append(mu)
    if flag==1:
        print('alpha_list:' + str(alpha_list) + '\n')
        print('beta_list:' + str(beta_list) + '\n')
        print('lambda_list:' + str(lam_list) + '\n')
        print('mu_list:' + str(mu_list) + '\n')

    return (alpha + beta) / 2


def fibonacci_search(f, section, delta):
    '''
    斐波那契法
    '''
    arr = [1, 1]  # Fibonacci数组
    stop_condition = (section[1] - section[0]) / delta  # 求n
    i = 2
    alpha = section[0]
    beta = section[1]
    alpha_list = [alpha]
    beta_list = [beta]
    lam_list = []
    mu_list = []

    while True:
        tmp = arr[i - 1] + arr[i - 2]
        arr.append(tmp)
        if tmp >= stop_condition:
            break
        i += 1

    k = len(arr) - 1
    lam = alpha + (arr[k - 2] / arr[k]) * (beta - alpha)
    mu = alpha + (arr[k - 1] / arr[k]) * (beta - alpha)
    lam_list.append(lam)
    mu_list.append(mu)

    while k >= 2:

        if f(lam) - f(mu) > 0:
            alpha = lam
            lam = mu
            mu = alpha + (arr[k - 1] / arr[k]) * (beta - alpha)
        else:
            beta = mu
            mu = lam
            lam = alpha + (arr[k - 2] / arr[k]) * (beta - alpha)

        alpha_list.append(alpha)
        beta_list.append(beta)
        lam_list.append(lam)
        mu_list.append(mu)
        k -= 1
    print('Fibonacci:' + str(arr) + '\n')
    print('alpha_list:' + str(alpha_list) + '\n')
    print('beta_list:' + str(beta_list) + '\n')
    print('lambda_list:' + str(lam_list) + '\n')
    print('mu_list:' + str(mu_list) + '\n')

    return (alpha + beta) / 2


def bin_split(f_dao, section, delta):
    '''
    二分法
    :param f_dao: 函数导数
    :param section:
    :param delta:
    :return:
    '''
    alpha = section[0]
    beta = section[1]
    alpha_list = [alpha]
    beta_list = [beta]
    lam_list = []
    dao = []
    while (beta - alpha >= delta):
        lam = (alpha + beta) / 2
        lam_list.append(lam)
        dao.append(f_dao(lam))
        if f_dao(lam) == 0:
            break
        elif f_dao(lam) > 0:
            beta = lam
        else:
            alpha = lam
        alpha_list.append(alpha)
        beta_list.append(beta)

    print('alpha_list:' + str(alpha_list) + '\n')
    print('beta_list:' + str(beta_list) + '\n')
    print('lambda_list:' + str(lam_list) + '\n')
    print('dao_list:' + str(dao) + '\n')
    return lam_list[-1]


def dichotomous(f, section, delta, eps):
    '''
    dichotomous法
    :param section:
    :param delta:
    :param eps:
    :return:
    '''
    alpha = section[0]
    beta = section[1]
    alpha_list = [alpha]
    beta_list = [beta]
    lam = (alpha + beta) / 2 - eps
    mu = (alpha + beta) / 2 + eps
    lam_list = [lam]
    mu_list = [mu]

    while (beta - alpha >= delta):

        if f(lam) > f(mu):
            alpha = lam
        else:
            beta = mu
        lam = (alpha + beta) / 2 - eps
        mu = (alpha + beta) / 2 + eps

        lam_list.append(lam)
        mu_list.append(mu)
        alpha_list.append(alpha)
        beta_list.append(beta)

    print('alpha_list:' + str(alpha_list) + '\n')
    print('beta_list:' + str(beta_list) + '\n')
    print('lambda_list:' + str(lam_list) + '\n')
    print('mu_list:' + str(mu_list) + '\n')

    return (alpha + beta) / 2



def goldstein(f,gradf,x,d,lamb,alpha,beta,p):
    '''
    :param f: 目标函数
    :param gradf: 函数梯度
    :param x: 初始点
    :param d: 初始方向
    :param lamb: 初始步长
    :param alpha: 步长增大系数
    :param beta: 步长缩短系数
    :return: 最佳lambda
    '''

    lamb_list=[lamb]
    while True:
        x_next=x+lamb*d
        # print(p)
        # print(gradf(x))
        # print(lamb)
        # print(d)
        # print(gradf(x).T.dot(lamb))
        if f(x_next)-f(x)>p*lamb*float(gradf(x).T.dot(d)):
            lamb*=beta
            lamb_list.append(lamb)
            continue
        if f(x_next)-f(x)<(1-p)*lamb*float(gradf(x).T.dot(d)):
            lamb*=alpha
            lamb_list.append(lamb)
            continue
        print('lambda的变化情况：'+str(lamb_list))
        return lamb
def goldstein_price(f,gradf,x,d,lamb,alpha,beta,p,sigma):
    '''
    :param f: 目标函数
    :param gradf: 函数梯度
    :param x: 初始点
    :param d: 初始方向
    :param lamb: 初始步长
    :param alpha: 步长增大系数
    :param beta: 步长缩短系数
    :param sigma: 条件2的参数
    :return: 最佳lambda
    '''

    lamb_list=[lamb]
    while True:
        x_next=x+lamb*d
        # print(p)
        # print(gradf(x))
        # print(lamb)
        # print(d)
        # print(gradf(x).T.dot(lamb))
        if f(x_next)-f(x)>p*lamb*float(gradf(x).T.dot(d)):
            lamb*=beta
            lamb_list.append(lamb)
            continue
        if f(x_next)-f(x)<sigma*lamb*float(gradf(x).T.dot(d)):
            lamb*=alpha
            lamb_list.append(lamb)
            continue
        print('lambda的变化情况：'+str(lamb_list))
        return lamb



def wolf_powell(f,gradf,x,d,lamb,alpha,beta,p):
    '''
    :param f: 目标函数
    :param gradf: 函数梯度
    :param x: 初始点
    :param d: 初始方向
    :param lamb: 初始步长
    :param alpha: 步长增大系数
    :param beta: 步长缩短系数
    :return: 最佳lambda
    '''

    lamb_list=[lamb]
    while True:
        x_next=x+lamb*d
        # print(p)
        # print(gradf(x))
        # print(lamb)
        # print(d)
        # print(gradf(x).T.dot(lamb))
        if f(x_next)-f(x)>p*lamb*float(gradf(x).T.dot(d)):
            lamb*=beta
            lamb_list.append(lamb)
            continue
        if gradf(x_next).T.dot(lamb*d)<(1-p)*lamb*float(gradf(x).T.dot(d)):
            lamb*=alpha
            lamb_list.append(lamb)
            continue
        print('lambda的变化情况：'+str(lamb_list))
        return lamb


def Hfunc(H, p, q):
    '''
    DFP算法的修正公式
    :param H:
    :param p:
    :param q:
    :return:
    '''
    return H+(1/(p.T.dot(q)))*p.dot(p.T)-(1/(q.T.dot(H).dot(q)))*(H.dot(q).dot(q.T).dot(H))

def Hfunc_bfgs(H,p,q):
    '''
    BFGS的修正公式
    :param H:
    :param p:
    :param q:
    :return:
    '''
    return H+(1+1/p.T.dot(q)*q.T.dot(H).dot(q))*(1/p.T.dot(q))*p.dot(p.T)-(1/p.T.dot(q))*(p.dot(q.T).dot(H)+H.dot(q).dot(p.T))

def quasi_Newton(x,eps,gradf,f,hfunc,**lists):
    '''
    拟牛顿法算法
    :param x: 初始点
    :param eps: 允许误差
    :param gradf: 梯度函数
    :param hfunc: 修正公式
    :return:
    '''



    # 设定初始值
    n=x.shape[0]
    g=gradf(x)
    H=np.identity(n)
    k=1
    d = -g.dot(H)

    def f_lambda(lamb_x):
        return f(x + lamb_x * d)
    lam = golden_section(f_lambda, [-100, 100], 0.03, flag=0)

    lists['x_list'].append(x)
    # lists['grad_list'].append(gradf(x))
    # lists['H_list'].append(H)
    # lists['d_list'].append(d)
    # lists['lambda_list'].append(lam)
    # lists['k_list'].append(k)

    while True:
        d = -g.dot(H)

        def f_lambda(lamb_x):
            return f(x + lamb_x * d)
        # 使用黄金分割法确定步长
        lam = golden_section(f_lambda, [-100, 100], 0.03,flag=0)
        x_next=x+lam*d
        if abs(np.linalg.norm(gradf(x_next)))<=eps:
            x=x_next
            break
        if k==n:
            lists['x_list'].append(x_next)
            lists['grad_list'].append(gradf(x_next))
            lists['H_list'].append(H)
            lists['d_list'].append(d)
            lists['lambda_list'].append(lam)
            lists['k_list'].append(k)
            return quasi_Newton(x_next,eps,gradf,f,hfunc,
                       x_list=lists['x_list'],
                       grad_list=lists['grad_list'],
                       H_list=lists['H_list'],
                       d_list=lists['d_list'],
                       lambda_list=lists['lambda_list'],
                       k_list=lists['k_list'])
        else:
            g_next=gradf(x_next)
            p=x_next-x
            q=g_next-g
            H=hfunc(H,p.reshape((2,1)),q.reshape((2,1)))
            x=x_next
            g=g_next
            k+=1

            lists['x_list'].append(x_next)
            # lists['grad_list'].append(gradf(x_next))
            # lists['H_list'].append(H)
            # lists['d_list'].append(d)
            # lists['lambda_list'].append(lam)
            # lists['k_list'].append(k)
    # print('k    x        grad_list             H_list                d_list      lambda_list')
    # print('----------------------------------------------------------------------------')
    lists['x_list'].append(x_next)
    print('x的取值变化：')
    for i in range(len(lists['x_list'])):
        print(
              # str(lists['k_list'][i])+'  '+
              str(lists['x_list'][i])+'  '
              # str(lists['grad_list'][i])+'  '+
              # str(lists['H_list'][i])+'  '+
              # str(lists['d_list'][i])+'  '+
              # str(lists['lambda_list'][i])
              )
    return x

if __name__ == '__main__':
    x=np.array([0.1,1.])
    eps=1e-6
    def gradf(x):
        return np.array([20*x[0],2*x[1]])

    def f(x):
        return 10*x[0]**2+x[1]**2

    res_x=quasi_Newton(x,eps,gradf,f,Hfunc,k_list=[],x_list=[],grad_list=[],H_list=[],d_list=[],lambda_list=[])
    print('最终x='+str(res_x))
    print('此时f='+str(f(res_x)))

    # x=np.array([2.,1.])
    # eps=1e-6
    # def gradf(x):
    #     return np.array([4*(x[0]-1),2*x[1]])
    #
    # def f(x):
    #     return 2*x[0]**2+x[1]**2-4*x[0]+2
    #
    #
    # print(dfp(x, eps, gradf, f, x_list=[], grad_list=[], H_list=[], d_list=[], lambda_list=[]))

