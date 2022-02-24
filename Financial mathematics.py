import random
import re
from math import ceil, exp, sqrt
from tkinter import NE

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
from numpy.linalg import inv


def linear_mortgage(interest, borrow, installment):
    arr = np.zeros((4, 2))
    arr[0, 0] = borrow
    redemption = 250000/360
    arr[2, 0] = 1072 - redemption
    i = 0
    continu = True
    while continu:
        try:
            arr[0, i+1]
        except IndexError:
            z = np.array([[0] for i in range(len(arr))])
            arr = np.append(arr, z, axis=1)
        arr[1, i] = arr[0, i]*interest
        arr[2, i] = min(arr[0, i], redemption)
        arr[3, i] = arr[1, i]+arr[2, i]
        try:
            arr[0, i+1] = arr[0, i] - arr[2, i]
            if arr[0, i] <= redemption:
                continu = False
        except:
            pass
        i += 1
    array = np.append(np.array(arr[:, :5]), np.array(arr[:, -5:]), axis=1)
    array = arr
    df = pd.DataFrame(
        array, columns=['1', '2', '3', '4', '5', '357', '358', '359', '360', '361'])
    df.index = ['Remainder', 'Interest', 'Redemption', 'Payment']
    return (df)
    # with open('linear_mortgage.txt', 'w') as f:
    #     f.write(df.to_latex())


def annuity_mortgage(interest, borrow, annuity):
    arr = np.zeros((4, 2))
    arr[0, 0] = borrow
    i = 0
    continu = True
    while continu:
        try:
            arr[0, i+1]
        except IndexError:
            z = np.array([[0] for i in range(len(arr))])
            arr = np.append(arr, z, axis=1)
        arr[1, i] = round(arr[0, i]*interest, 9)
        arr[3, i] = min(arr[0, i]+arr[1, i], annuity)
        arr[2, i] = arr[3, i]-arr[1, i]
        try:
            arr[0, i+1] = arr[0, i]-arr[2, i]
            if arr[0, i] <= annuity:
                continu = False
        except:
            pass
        i += 1
    df = pd.DataFrame(arr, columns=[f'{i}' for i in range(1, len(arr[0])+1)])
    return df


def simple_interest(pv, r, periods):
    return (1+r*periods)*pv


def coumpoundInterestAtEndOfYear(pv, r, periods):
    return ((1+r*periods)**periods)*pv


def compoundInterestAtVariousIntervals(pv, r, m,  periods):
    return ((1+r/m)**(m*periods))*pv


def continousCompoundInterest(pv, r, periods):
    return exp(r*periods)*pv


def PV_annuity(A, r, periods):
    return sum([A/(1+r)**i for i in range(1, periods+1)]), (A/r)*(1-(1/(1+r)**periods))


def xrange(start=0, end=0, interval=1):
    return [i*interval for i in range(start, end)]


def NetPV(periods, r, m=1, interval=1):
    c = 1/(1+r/m)
    return round(sum([periods[ceil(i)]*c**(i*m) for i in xrange(0, len(periods), interval)]), 3)


def IRR(xi):
    n = len(xi)
    for i in range(0, 100):
        for j in range(0, 100):
            for k in range(0, 100):
                k /= 10
                r = (i+j+k)/100
                c = 1/(1+r)
                net_pv = round(sum({xi[i]*c**i for i in range(n)}), 1)

                if net_pv <= 0.001 and net_pv >= -0.001:

                    return r
    else:
        print('notfound')


def annuity_mortgage(r):
    PV = (904/(r/12))*(1-(1/(1+(r/12))**(30*12)))
    return PV - 250000  # net PV


def linear_mortgage(r):
    PV = sum([1072/(1+(r/12))**(12*(i))for i in range(1, 31)])
    return PV - 250000  # net PV


def PresentValueBONDS(FV, couponRate, maturity):
    # if not isinstance(FV,list):
    #     FV = [FV for i in range(maturity)]
    pv = sum([couponRate*FV/(1+(couponRate))**i for i in range(1,
             maturity+1)]) + FV/(1+couponRate**maturity)
    return pv


def DurationBONDS(FV, couponRate, maturity):
    if isinstance(couponRate, int):
        couponRate = [couponRate for i in range(maturity)]
    dur = (1/(PresentValueBONDS)) * \
        (couponRate[t-1]*FV*t/(1+couponRate)**t for t in range(1, maturity+1))
    return dur


def DurationContiuationsBONDS(FV, couponRate, maturity):
    if isinstance(couponRate, int):
        couponRate = [couponRate for i in range(maturity)]


print(
    # NetPV(x, 0.0184),
    #   PV_annuity(904, 0.0184, 30*12)
    # ((1+0.08/2)**2 -1),
    # NetPV([300,0,230,0,300],0.08,2, 0.5),
    # NetPV([0,100,120,100,500],((1+0.08/2)**2 -1),2)
    # linear_mortgage(interest1, borrow1, redemption1),
    # annuity_mortgage(interest, borrow, redemption)
)
mew = np.array([0.1, 0.2])
μ = np.array([[0.1], [0.2]])
cov = np.array([[0.2, 0.1],
                [0.1, 0.3]]
               )
std_range = [i for i in range(20)]
mew_range = [i for i in range(20)]

varlist = []
varlist1 = []
mewplotlist = []
mewplotlist1 = []
# for w in range(-100000, 100001):
for w in range(0, 11):
    break
    w /= 10
    W = np.array([[w, 1-w]])
    mewplot = w*mew[0]+(1-w)*mew[1]
    # print(W, μ)
    mewplot1 = np.matmul((W), μ)
    # print(mewplot, mewplot1)
    var = ((w**2)*cov[0, 0]+2*w*(1-w)*(cov[0, 1]) + ((1-w)**2)*cov[1, 1])
    sig2 = np.matmul(np.matmul(np.transpose(W), cov), W)
    mewplotlist.append(mewplot)
    varlist.append(var)
    mewplotlist1.append(mewplot1)
    varlist1.append(sig2)

plt.plot(varlist, mewplotlist)
plt.plot()
plt.plot(varlist1, mewplotlist)

# plt.show()
# print(varlist[:4], varlist1[:4])


def plot_efficient_frontier(cov, μt):
    μ = np.transpose(μt)
    varlist = []
    mewplotlist = []
    ε = 0
    for ε in range(-100, 101):
        # for w in range(-100000, 100001):
        ε /= 100000
        Wt = inital_Wt(cov)
        W = np.transpose(Wt)
        # print(Wt.shape, μ.shape)
        mewplot = np.matmul(Wt, μ)[0][0]
        sig2 = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))[0][0]
        # print(mewplot, sig2)
        mewplotlist.append(mewplot)
        varlist.append(sig2)

    plt.plot(varlist, mewplotlist)
    plt.plot()

    plt.show()


def inital_Wt(cov):
    ones = np.ones((len(cov), 1))
    Wt = np.matmul(np.transpose(ones), inv(cov)) / \
        np.matmul((np.matmul(np.transpose(ones), inv(cov))),    ones)
    return Wt


def getMin_μσ2(μt, cov, show=True):
    Wt = inital_Wt(cov)
    σ2_min = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))
    μ_min = np.matmul(Wt, np.transpose(μt))
    if show:
        return f'μ_min  = {μ_min[0][0]}  \nσ2_min = { σ2_min[0][0]}'
    return μ_min[0][0], σ2_min[0][0]


def plotFeasibleSet(μ, cov):
    varlist1, mewplotlist1= [],[]
    for i in range(100000):
        W = np.random.random((len(cov), 1))
        W /= np.sum(W)
        Wt = np.transpose(W)
        mewplot1 = np.matmul(Wt, μ)[0][0]
        sig21 = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))[0][0]
        mewplotlist1.append(mewplot1), varlist1.append(sig21)
    plt.scatter(varlist1, mewplotlist1, c='none', s=5, edgecolors='r')
    plt.title('Feasible Set')
    plt.xlabel('σ^2 (volatility)')
    plt.ylabel('μ (return)')
    plt.scatter(getMin_μσ2(μt, cov, False)[1], getMin_μσ2(
        μt, cov, False)[0], c='red', edgecolors='r')
    plt.show()


def plotEfficientFrontier(μ, cov):
    varlist = []
    mewplotlist = []
    muval = 0
    for i in range(300):
        i /= 1000
        W = solveLP(i, μ, cov)
        mewplot = np.matmul(np.transpose(W), μ)[0][0]
        sig2 = np.matmul(np.matmul(np.transpose(W), cov), (W))[0][0]
        mewplotlist.append(mewplot)
        varlist.append(sig2)
    # varlist1 = []
    # mewplotlist1 = []
    # for i in range(100000):
    #     W = np.random.random((len(cov), 1))
    #     W /= np.sum(W)
    #     Wt = np.transpose(W)
    #     mewplot1 = np.matmul(Wt, μ)[0][0]
    #     sig21 = np.matmul(np.matmul(Wt, cov), np.transpose(Wt))[0][0]
    #     mewplotlist1.append(mewplot1)
    #     varlist1.append(sig21)
    # plt.scatter(varlist1, mewplotlist1, c='none', s=2, edgecolors='r')
    plt.title('Feasible Set')
    plt.xlabel('σ^2 (volatility)')
    plt.ylabel('μ (return)')
    plt.scatter(varlist, mewplotlist, c='none', s=2, edgecolors='b')
    plt.scatter(getMin_μσ2(μt, cov, False)[1], getMin_μσ2(
        μt, cov, False)[0], c='red', edgecolors='r')
    plt.show()


def solveLP(muval, mu, cov):
    if muval == 0:
        muval = getMin_μσ2(μt, cov, False)[0]
    RHS = np.array([[0], [0], [0], [muval], [1]])
    LHS = np.append(cov, -mu, axis=1)
    LHS = np.append(LHS, -np.ones((len(cov), 1)), axis=1)
    LHS = np.vstack((LHS, np.append(mu, [0, 0])))
    LHS = np.vstack((LHS, np.append(np.ones(len(cov)), [0, 0])))

    return (np.matmul(inv(LHS), RHS)[:3])


if __name__ == '__main__':
    cov = np.array([
        [0.2, 0.1, 0.1],
        [0.1, 0.3, 0.1],
        [0.1, 0.1, 0.4]
    ])
    μt = np.array([[0.1, 0.2, 0.3]])
    μ = np.transpose(μt)
    W = np.transpose(inital_Wt(cov))

    print(
        solveLP(0, μ, cov),
        '\n',
        plotFeasibleSet(μ, cov)
    )
