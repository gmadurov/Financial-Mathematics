
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.linalg import inv


DF = pd.read_csv('Assignment3-data.csv', delimiter=';')
AR = np.array(DF)


def inital_W(C):
    '''finds the initial weight matrix W, given a covariance matrix C'''
    ones = np.ones((len(C), 1))
    W = np.transpose(np.matmul(np.transpose(ones), inv(C)) /
                     np.matmul((np.matmul(np.transpose(ones), inv(C))),    ones))
    return W


def getMin_mewsigma2(mu, C, show=True):
    '''gets the minimum value of mu and sigma**2 for a given expected return matrix mu\n
    and a covariance matrix C, if show is true function returns a string'''
    W = inital_W(C)
    sigma2_min = np.matmul(np.matmul(np.transpose(W), C), W)
    mew_min = np.matmul(np.transpose(W), mu)
    if show:
        return f'mew_min  = {mew_min[0][0]}  \nsigma2_min = { sigma2_min[0][0]}'
    return mew_min[0][0], sigma2_min[0][0]


def get_mewsigma2(mu, C, W, show=True):
    '''gets the value of mu and sigma**2, for a given expected return matrix mu\n and a covariance matrix C, and a weight matrix W.\nif show is true function returns a string'''
    sigma2 = np.matmul(np.matmul(np.transpose(W), C), W)
    mu = np.matmul(np.transpose(W), mu)
    if show:
        return f'mu  = {mu[0][0]}  \n sigma2 = { sigma2[0][0]}'
    return mu[0][0], sigma2[0][0]


def Markowitz(muval, mu, C):
    '''Solve the Markowitz problem given a value specific value of return- muval
    the excpected return vector (mu) and the covariance matrix (C)'''
    if muval == 0:
        muval = getMin_mewsigma2(mu, C, False)[0]
    # set up vector for all the weights on the Right Hand Side (RHS)
    RHS = np.array([[0] for i in range(len(mu))])
    # adds the value of mu and 1 to the RHS vector (forming the last 2 constraints)
    RHS = np.vstack((RHS, np.array([[muval], [1]])))
    # the following gets to the covariance matrix and adds the labmda variable coefficients
    # to create the Left Hand Side of the equation
    LHS = np.append(C, -mu, axis=1)
    LHS = np.append(LHS, -np.ones((len(C), 1)), axis=1)
    LHS = np.vstack((LHS, np.append(mu, [0, 0])))
    LHS = np.vstack((LHS, np.append(np.ones(len(C)), [0, 0])))
    # this line solves it and returns (ONLY) the values of the weights
    # return pd.DataFrame(LHS), pd.DataFrame(RHS)
    return (np.matmul(inv(LHS), RHS)[:len(mu)])


def plotEfficientFrontier(mu, C, lb=-1, ub=2, show=True):
    varlist = []
    mewplotlist = []
    colour = []
    for i in range(lb*1000, ub*1000):  # change this range to what you want to show
        i /= 1000
        W = Markowitz(i, mu, C)
        if i < 0:
            colour.append('red')
        elif i > 0 and i <= 1:
            colour.append('green')
        elif i > 1:
            colour.append('blue')
        else:
            colour.append('none')
        mewplot = np.matmul(np.transpose(W), mu)[0][0]
        sig2 = sqrt(np.matmul(np.matmul(np.transpose(W), C), (W))[0][0])
        mewplotlist.append(mewplot)
        varlist.append(sig2)
    if show:
        plt.title('Efficient Frontier')
        plt.xlabel('sigma (volatility)')
        plt.ylabel('mu (return)')
        plt.scatter(varlist, mewplotlist, c=colour, s=2, )
        plt.show()
        return
    return list(zip(varlist, mewplotlist))


DF_R = pd.DataFrame()
for col in DF.columns[1:]:
    # DF_R.loc[0, 'r_'+col] = 0 # taken out because you only use sum of 1 to T
    for i in range(1, len(DF[col])):
        try:
            DF_R.loc[i, 'r_'+col] = ((DF.loc[i, col] -
                                     DF.loc[i-1, col]) / DF.loc[i-1, col])
        except IndexError:
            print(f'error in {i}, in {col}')
R_array = np.array(DF_R)
# R_array
mu = np.array([[i] for i in np.sum(R_array, axis=0)])
C = np.cov(R_array, rowvar=False)*365
# C, mu


def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv += [r'\end{bmatrix}']
    return '\n'.join(rv)


# print(bmatrix(np.transpose(mu)))


def two_fund_theorem_from_matrix(mu, C, r, show=True, fig=None):
    if not fig:
        fig = plt.figure(figsize=(6, 5))
    ones = np.ones((len(C), 1))
    W_min = inital_W(C)
    W_market = np.transpose(np.matmul(np.transpose(mu-r*ones), inv(C)) /
                            np.matmul((np.matmul(np.transpose(mu-r*ones), inv(C))), ones))
    # W_market = ((mu-r*ones).T@ inv(C))/ (mu-r*ones).T@inv(C)@ones
    mew1, sigma21 = get_mewsigma2(mu, C, W_min, False)
    mew2, sigma22 = get_mewsigma2(mu, C, W_market, False)
    # if mew1 < mew2 and sigma21 < sigma22:
    #     mew1, sigma21, mew2, sigma22 = mew2, sigma22, mew1, sigma21
    mewpoints, sigmapoints, colours = [], [], []
    covMinMarket = np.matmul(np.matmul(np.transpose(W_min), C), W_market)

    for w in range(-10, 130):
        w /= 100
        if w < 0:
            colours.append('red')
        elif w > 0 and w < 1:
            colours.append('green')
        elif w > 1:
            colours.append('blue')
        else:
            colours.append('black')
        sigmav = sqrt(((w**2)*(sigma21)) + (2*w*(1-w) *
                  covMinMarket) + (((1-w)**2)*(sigma22)))
        mewv = w*mew1+(1-w)*mew2
        mewpoints.append(mewv)
        sigmapoints.append(sigmav)
    axarr = fig.add_subplot(1, 1, 1)
    plt.plot(sigmapoints, mewpoints)
    if show:
        plt.title('Efficient Frontier(two fund theorem)')
        plt.xlabel('sigma (volatility)')
        plt.ylabel('mu (return)')
        plt.scatter((sigma1), mew1,  c='b')
        plt.scatter((sigma2), mew2,  c='r')
        plt.show()
        return
    return sigmapoints, mewpoints


def OneFundTheorem(mu, C, rf, end=3.5):
    ones = np.ones((len(C), 1))
    W_market = np.transpose(np.matmul(np.transpose(mu-rf*ones), inv(C)) /
                            np.matmul((np.matmul(np.transpose(mu-rf*ones), inv(C))),    ones))
    W_min = inital_W(C)
    market_mew, sigma2p1 = get_mewsigma2(mu, C, W_market, False)
    mew_min, sigma_min = get_mewsigma2(mu, C, W_min, False)
    market_sigma = sqrt(sigma2p1)
    sigma_min = sqrt(sigma_min)
    # print(market_mew, market_sigma)
    print()
    sigma = np.linspace(0, end)
    mew = rf+sigma*(market_mew-rf)/market_sigma
    fig = plt.figure(figsize=(10, 8))
    # two_fund_theorem_from_points(
    #     market_mew, market_sigma, mew_min, sigma_min, covMinMarket, show=False, fig=fig)
    two_fund_theorem_from_matrix(mu, C, rf, show=False, fig=fig)
    plt.title('Efficient Frontier(one fund theorem)')
    plt.xlabel('sigma (volatility)')
    plt.ylabel('mu (return)')
    plt.scatter((0), rf,  c='green', s=30)
    plt.plot(sigma, mew, c='black', linestyle='dashed')
    plt.scatter((sigma_min), mew_min,  c='red')
    plt.scatter((market_sigma), market_mew,  c='green', s=30)
    plt.show()
    # min = 'min'
    print(
        f'$\mu_{{min}} = {round(mew_min, 6)}$ and  $\sigma_{{min}}^2 = {round(sigma_min, 6)}^2$')
    print(
        f'$\mu_{{M}} = {round(market_mew, 6)}$ and  $\sigma_{{M}}^2 = {round(market_sigma, 6)}^2$')
    print('min weights ', bmatrix(np.transpose(W_min)))
    print('market weights ', bmatrix(np.transpose(W_market)))


OneFundTheorem(mu, C, 1/100)

R_array2 = R_array[:, :4]
mew2 = np.array([[i] for i in np.sum(R_array2, axis=0)])
C2 = np.cov(R_array2, rowvar=False)*365

OneFundTheorem(mew2, C2, 1/100, 6.5)
