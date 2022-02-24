

import enum
from math import exp

import numpy as np


def main():
    x = [0 for i in range(24)]+[250000]
    r = [i/1000 for i in range(1,26)]
    if len(r) != len(x):
        raise ValueError('vector r and x are not equal in length')

    def P(x, r): return round(
        sum([(x[t]*exp(-r[t]*(t+1))) for t, e in enumerate(x)]), 4)

    def D(x, r): return round(
        1/P(x, r)*sum([((t+1)*x[t]*exp(-r[t]*(t+1))) for t, e in enumerate(x)]), 4)

    def D1(x, r): return round(
        1/P(x, r)*sum([((t+1)*(r[t])*x[t]*exp(-r[t]*(t+1))) for t, e in enumerate(x)]), 4)

    def D2(x, r):
        return round(1/P(x, r)*sum([((t+1)*(r[t]**2)*x[t]*exp(-r[t]*(t+1))) for t, e in enumerate(x)]), 4)

    immune = [100, 0, 100, 0, 100]

    def immunize(immunne, x, r):
        # np.linalg.solve
        pvs = np.array([(s*exp(-r[t]*(t+1))) for t, s in enumerate(immune) if s!= 0])
        RHS = np.array([P(x, r), P(x, r)*D(x, r), P(x, r)*D1(x, r)])
        LHS = np.array([pvs, pvs*[1,  3,  5],
                        pvs*[0.01,  0.09,  0.25]])
        return (np.linalg.solve(LHS, RHS))
    print(
        # immunize(immune, x, r)
        # P(x, r), D(x, r), D1(x, r), D2(x, r)
    )
   
if __name__ == '__main__':
    main()
