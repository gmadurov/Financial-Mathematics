

from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import inv


def inital_Wt(cov):
    ones = np.ones((len(cov), 1))
    Wt = np.matmul(np.transpose(ones), inv(cov)) / \
        np.matmul((np.matmul(np.transpose(ones), inv(cov))),    ones)
    return Wt


def getMin_μσ2(μt, cov, show=True):
    Wt = inital_Wt(cov)
    σ2_min = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))
    μ_min = np.matmul(Wt, np.transpose(μt))
    if show:return f'μ_min  = {μ_min[0][0]}  \nσ2_min = { σ2_min[0][0]}'
    return μ_min[0][0],σ2_min[0][0]

def getMin_μσ2W(μt, cov, W, show=True):
    Wt = np.transpose(W)
    σ2_min = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))
    μ_min = np.matmul(Wt, np.transpose(μt))
    if show:return f'μ_min  = {μ_min[0][0]}  \nσ2_min = { σ2_min[0][0]}'
    return μ_min[0][0],σ2_min[0][0]

def getμσ2(μ, cov):
    varlist = []
    mewplotlist = []
    for i in range(100000):
            W = np.random.random((len(cov),1))
            W /= np.sum(W)
            Wt= np.transpose(W)
            mewplot = np.matmul(Wt, μ)[0][0]
            sig2 = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))[0][0]
            mewplotlist.append(mewplot)
            varlist.append(sig2)
    plt.title('Feasible Set')
    plt.xlabel('σ^2 (volatility)')
    plt.ylabel('μ (return)')
    plt.scatter(varlist, mewplotlist, c='none',s=2, edgecolors='b')
    plt.scatter(getMin_μσ2(μt, cov,False)[1],getMin_μσ2(μt, cov,False)[0], c='red', edgecolors='r')
    plt.show()


cov = np.array([
    [0.2, 0.1, 0.1],
    [0.1, 0.3, 0.1],
    [0.1, 0.1, 0.4]
])
μt = np.array([[0.1, 0.2, 0.3]])
μ = np.transpose(μt)
# W = np.transpose(inital_Wt(cov))


# cov = np.array([
#     [0.2,-0.1],
#     [-0.1,0.5],
# ])
# μt = np.array([[0.1, 0.4]])
# μ = np.transpose(μt)

# W = np.array([[0.5],[0.5]])
print(
    getMin_μσ2(μt, cov),
    # getMin_μσ2W(μt,cov,W),
    getμσ2(μ, cov),
        # μ 
)
