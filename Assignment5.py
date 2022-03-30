import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt, exp, factorial
import pprint
from scipy.stats import norm


def range_(start=0, end=1, factor=1, step=1):
    return [i * factor for i in range(start, end, step)]


def bicoef(top, bottom):
    return factorial(top) / (factorial(bottom) * factorial(top - bottom))


class Option:
    def __init__(self, v, sigma, number_of_periods, T, K, S0):
        self.S = {0: S0}
        self.v = v
        self.sigma = sigma
        self.number_of_periods = number_of_periods
        self.deltaT = 1 / self.number_of_periods
        self.T = T
        self.K = K
        self.excercise_matrix = {}
        self.C = {}
        self.prob = 0

        self.get_values()

    def get_values(self):
        self.r = self.v + 1 / 2 * self.sigma**2
        self.R = exp(self.r * (self.deltaT))
        self.u = exp(self.sigma * sqrt(self.deltaT))
        self.d = 1 / self.u
        self.p = 1 / 2 + 1 / 2 * self.v / self.sigma * sqrt(self.deltaT)
        self.q = (self.R - self.d) / (self.u - self.d)

    def binomialMap(self):

        for j in range_(-self.number_of_periods, 1, -1):
            for i in range(0, j + 1):
                self.S["U" * i + "D" * (j - i)] = max(
                    0, (self.S[0] * (self.u**i) * (self.d ** (j - i)))
                )
        """initial RHS of the C tree"""
        for i in range_(0, self.number_of_periods + 1):
            self.C["U" * i + "D" * (self.number_of_periods - i)] = max(
                0, (self.S["U" * i + "D" * (self.number_of_periods - i)] - self.K)
            )
        """Calculate the rest of the C tree"""
        for j in range_(-self.number_of_periods + 1, 1, -1):
            for i in range_(0, j + 1):
                self.C["U" * i + "D" * (j - i)] = (
                    1
                    / self.R
                    * (
                        self.q * self.C["U" * (i) + "U" + "D" * (j - i)]
                        + (1 - self.q) * self.C["U" * (i) + "D" * (j - i) + "D"]
                    )
                )
        """making the excersice matrix for american option"""
        for j in range_(-self.number_of_periods, 1, -1):
            for i in range_(0, j + 1):
                if (
                    self.C["U" * i + "D" * (j - i)]
                    <= self.S["U" * i + "D" * (j - i)] - self.K
                ):
                    self.excercise_matrix["U" * i + "D" * (j - i)] = "excercise"
                else:
                    self.excercise_matrix["U" * i + "D" * (j - i)] = "not excercise"
        for i in range_(-self.number_of_periods, 0, -1):
            if (
                self.excercise_matrix["U" * i + "D" * (self.number_of_periods - i)]
                == "excercise"
            ):
                self.prob += (
                    self.p**i
                    * (1 - self.p) ** (self.number_of_periods - i)
                    * bicoef(self.number_of_periods, i)
                )
            else:
                pass
                # print('U'*i+'D'*(self.number_of_periods-i))
        dic = {}
        for k in self.C.keys():
            dic[k] = [self.C[k], self.S[k], self.excercise_matrix[k]]
        print(
            pd.DataFrame((dic.values()), index=dic.keys(), columns=["C", "S", "Option"])
        )
        return self.prob

    def Black_Scholes(self, time=1):
        self.D1 = {}
        self.D2 = {}

        for t in range(time):
            self.D1[t] = (
                1
                / (self.sigma * sqrt(self.T - t))
                * (
                    np.log(self.S[0] / self.K)
                    + ((self.r + self.sigma**2 / 2) * (self.T - t))
                )
            )
            self.D2[t] = (
                1
                / (self.sigma * sqrt(self.T - t))
                * (
                    np.log(self.S[0] / self.K)
                    + ((self.r - self.sigma**2 / 2) * (self.T - t))
                )
            )
        call_price = norm.cdf(self.D1[0]) * self.v - norm.cdf(self.D2[0]) * K * exp(
            -self.r * (t)
        )
        put_price = call_price - self.S[0] + self.K * exp(-self.r * self.T)

        return "call prize is", call_price, "Put prize is ", put_price


def make_binomial_tree(number_of_periods):
    columns = {}
    columns2 = {}
    if number_of_periods < 4:
        st = """\\begin{tikzpicture}[>=stealth,sloped]
            \\matrix (tree) [%
            matrix of nodes,
            minimum size=1cm,
            column sep=3.5cm,
            row sep=1cm,
            ]
        {"""
    else:
        st = """\\begin{tikzpicture}[>=stealth,sloped]
            \\matrix (tree) [%
            matrix of nodes,minimum size=0.5cm,
        column sep=1cm,
        row sep=0.5cm,]{"""
    st2 = ""
    for j in range_(-number_of_periods, 1, -1):
        for i in range_(0, j + 1):
            try:
                columns[j].append("")
                columns2[j].append("")
            except:
                columns[j] = ["" for i in range(number_of_periods - j)]
                columns2[j] = ["" for i in range(number_of_periods - j)]
            columns2[j].append(("U" * i + "D" * (j - i)))
            columns[j].append(
                (
                    f"""{'$S_0U' if i==1 else f'$S_0U^{i}'if i>0 else '$S_0'} {"D$" if j-i==1 else f'D^{(j-i)}$'if j-i>0 else '$'}"""
                    if not (i == j and j == 0)
                    else "$S_0$"
                )
            )
        for i in range(number_of_periods - j):
            columns[j].append("")
            columns2[j].append("")
    for i in range_(-2 * number_of_periods, 1, -1):
        for j in range(len(columns)):
            if j != len(columns) - 1:
                # print(columns2[j][i])
                try:
                    st += S[columns2[j][i]] + " & "
                except:
                    st += columns[j][i] + " & "
                if columns[j][i] != "":
                    U1, D1, U2, D2 = (
                        columns2[j + 1][i - 1].count("D"),
                        columns2[j + 1][i - 1].count("U"),
                        columns2[j + 1][i + 1].count("D"),
                        columns2[j + 1][i + 1].count("U"),
                    )
                    # adding the ligns and their values between nodes on the tree
                    st2 += (
                        (
                            f"\draw[->] (tree-{i+1}-{j+1}) -- (tree-{i}-{j+2}) node [midway,above] "
                        )
                        + f'{"{"}$ {"p" if U1 ==1 else f"p^{U1}" if U1 > 0 else ""} {"(1-p)" if D1 == 1 else f"(1-p)^{D1}" if D1 > 0 else ""} ${"}"};\n'
                        + (
                            f"\draw[->] (tree-{i+1}-{j+1}) -- (tree-{i+2}-{j+2}) node [midway,above] "
                        )
                        + f'{"{"}$ {"p" if U2 ==1 else f"p^{U2}" if U2 > 0 else ""} {"(1-p)" if D2 == 1 else f"(1-p)^{D2}" if D2 > 0 else ""} ${"}"};\n'
                    )
            else:
                st += columns[j][i]

        st += " \\\\\n"
    st += "};\n" + st2 + "\end{tikzpicture}"
    return st


S = {}
v = -0.133751
sigma = sqrt(0.554548)
number_of_periods = 12
deltaT = 1 / number_of_periods
T = 1
K = 42000
S = 42600.132813
Option(v, sigma, number_of_periods, T, K, S).binomialMap()
