import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt, exp, factorial
from scipy.stats import norm


def range_(start=0, end=1, factor=1, step=1):
    """allows for range to be a float and in reverse"""
    return [i * factor for i in range(start, end, step)]


def bicoef(top, bottom):
    """simple binomial coeficient function"""
    return factorial(top) / (factorial(bottom) * factorial(top - bottom))


class Option:
    """option pricing model h in the diviend(doesnt really work yet)"""

    def __init__(
        self,
        S0,
        sigma,
        number_of_periods,
        Time_horizon,
        strike_prices,
        log_variation=None,
        interest_rate=None,
        deltaT=None,
        h=0,
    ):
        """start the values and call get values function which calculates R,u,d,p,q"""
        self.S = {"": "", 0: S0}
        self.v = log_variation or interest_rate - 1 / 2 * sigma**2
        self.sigma = sigma
        self.number_of_periods = number_of_periods
        self.deltaT = deltaT or (1 / self.number_of_periods)
        self.T = Time_horizon
        self.K = strike_prices
        self.excercise_matrixCall = (
            {}
        )  # whether you should excersice a call option should it be american
        self.excercise_matrixPut = (
            {}
        )  # whether you should excersice a put option should it be american
        self.C = {}  # for calls tree
        self.P = {}  # for put tree
        self.probCall = 0  # probability of bieng in the money put option
        self.probPut = 0  # prob of being in the money for a put option
        self.h = h
        self.r = interest_rate or self.v + 1 / 2 * self.sigma**2
        self.get_values()

    def get_values(self):
        """calculates R,u,d,p,q"""
        self.R = exp(self.r * (self.deltaT))
        self.u = exp(self.sigma * sqrt(self.deltaT))
        self.d = 1 / self.u
        self.p = 1 / 2 + 1 / 2 * self.v / self.sigma * sqrt(self.deltaT)
        self.q = (self.R - self.d) / (self.u - self.d)

    def binomialMap(self, show=True):
        """makes a binomial map calculates all the values needed thus C and S and whether the option should be maximized
        show = True outputs a table with all the values and probability of being in the money
        false  outputs a tuple of lists of values"""

        """ the folowing lines: 
        for j in range_(-self.number_of_periods, 1, -1):
            for i in range_(0, j + 1):
                

        allow us to work thorugh a dictionary in reverse order so from the right hand side of the tree to the left hand side of tree
        this was done because the tree needs to be calculated from left to right and using the "range"function did suffice
                """

        for j in range_(-self.number_of_periods, 1, -1):
            for i in range_(0, j + 1):
                self.S["U" * i + "D" * (j - i) + "S0"] = max(
                    0, (self.S[0] * (self.u**i) * (self.d ** (j - i)))
                )
        """initial RHS of the C and P tree"""
        for i in range_(0, self.number_of_periods + 1):
            self.C["U" * i + "D" * (self.number_of_periods - i) + "S0"] = max(
                0,
                (self.S["U" * i + "D" * (self.number_of_periods - i) + "S0"] - self.K),
            )
            self.P["U" * i + "D" * (self.number_of_periods - i) + "S0"] = max(
                0,
                (self.K - self.S["U" * i + "D" * (self.number_of_periods - i) + "S0"]),
            )
        """Calculate the rest of the C and P tree"""
        for j in range_(-self.number_of_periods + 1, 1, -1):
            for i in range_(0, j + 1):
                self.C["U" * i + "D" * (j - i) + "S0"] = (
                    1
                    / self.R
                    * (
                        self.q * self.C["U" * (i) + "U" + "D" * (j - i) + "S0"]
                        + (1 - self.q) * self.C["U" * (i) + "D" * (j - i) + "D" + "S0"]
                    )
                )
                self.P["U" * i + "D" * (j - i) + "S0"] = (
                    1
                    / self.R
                    * (
                        self.q * self.P["U" * (i) + "U" + "D" * (j - i) + "S0"]
                        + (1 - self.q) * self.P["U" * (i) + "D" * (j - i) + "D" + "S0"]
                    )
                )
        """making the excersice matrix for american option"""
        for j in range_(-self.number_of_periods, 1, -1):
            for i in range_(0, j + 1):
                if (
                    self.C["U" * i + "D" * (j - i) + "S0"]
                    <= self.S["U" * i + "D" * (j - i) + "S0"] - self.K
                ):
                    self.excercise_matrixCall[
                        "U" * i + "D" * (j - i) + "S0"
                    ] = "excercise"
                else:
                    self.excercise_matrixCall[
                        "U" * i + "D" * (j - i) + "S0"
                    ] = "not excercise"
                if (
                    self.P["U" * i + "D" * (j - i) + "S0"]
                    >= self.S["U" * i + "D" * (j - i) + "S0"] - self.K
                ):
                    self.excercise_matrixPut[
                        "U" * i + "D" * (j - i) + "S0"
                    ] = "excercise"
                else:
                    self.excercise_matrixPut[
                        "U" * i + "D" * (j - i) + "S0"
                    ] = "not excercise"
        for i in range_(-self.number_of_periods, 0, -1):
            if (
                self.excercise_matrixCall[
                    "U" * i + "D" * (self.number_of_periods - i) + "S0"
                ]
                == "excercise"
            ):
                self.probCall += (
                    self.p**i
                    * (1 - self.p) ** (self.number_of_periods - i)
                    * bicoef(self.number_of_periods, i)
                )
            if (
                self.excercise_matrixPut[
                    "U" * i + "D" * (self.number_of_periods - i) + "S0"
                ]
                == "excercise"
            ):
                self.probPut += (
                    self.p**i
                    * (1 - self.p) ** (self.number_of_periods - i)
                    * bicoef(self.number_of_periods, i)
                )
        if show:
            dic = {}
            for k in self.C.keys():
                dic[k] = [self.C[k], self.P[k], self.S[k], self.excercise_matrixCall[k]]
            print(
                pd.DataFrame(
                    (dic.values()), index=dic.keys(), columns=["C", "P", "S", "Option"]
                ).tail()
            )
            return "call probability", self.probCall, "put probability", self.probPut
        else:
            return (
                self.S,
                self.C,
                self.P,
                self.excercise_matrixCall,
                self.excercise_matrixPut,
            )

    def Black_Scholes(self):
        """calculate the black scholes model"""
        self.D1 = (
            1
            / (self.sigma * sqrt(self.T))
            * (
                np.log(self.S[0] / self.K)
                + ((self.r - self.h + self.sigma**2 / 2) * (self.T))
            )
        )
        self.D2 = (
            1
            / (self.sigma * sqrt(self.T))
            * (
                np.log(self.S[0] / self.K)
                + ((self.r - self.h - self.sigma**2 / 2) * (self.T))
            )
        )
        # you might need to add smal t to the equation here
        self.call_price = norm.cdf(self.D1) * self.S[0] * exp(
            -self.h * (self.T)
        ) - norm.cdf(self.D2) * self.K * exp(-self.r * (self.T))
        self.put_price = self.call_price - self.S[0] + self.K * exp(-self.r * self.T)

        return (
            "Call prize is "
            + str(self.call_price)
            + " w.p. "
            + str(norm.cdf(self.D2))
            + "\nPut prize is "
            + str(self.put_price)
            + " w.p. "
            + str(1 - norm.cdf(self.D2))
        )

    def greeks_call(self, show=False):
        """returns --> deltaCall, thetaCall, gammaCall, vegaCall, rho
        this was only calculated for calls because of the assignment purpose"""
        self.Black_Scholes()
        self.deltaCall = norm.cdf(self.D1)
        self.thetaCall = (self.S[0] * norm.cdf(self.D1) * self.sigma) / (
            2 * sqrt(T)
        ) - self.r * self.K * exp(-self.r * self.T) * norm.cdf(self.D2)
        self.gammaCall = norm.cdf(self.D1) / (self.S[0] * self.sigma * sqrt(T))
        self.vegaCall = self.S[0] * sqrt(T) * norm.cdf(self.D1)
        self.rho = self.K * self.T * exp(-self.r * self.T) * norm.cdf(self.D2)
        if show:
            print(
                f"Delta = {self.deltaCall},\nTheta = {self.thetaCall},\nGamma = {self.gammaCall},\nVega = {self.vegaCall},\nRho = {self.rho}"
            )
        return self.deltaCall, self.thetaCall, self.gammaCall, self.vegaCall, self.rho

    def make_binomial_tree(self):
        """takes data from the binomial map fucntion and puts in into a binomial lattice tree in latex just copy paste the PRINTED output and it should plot a good graph in latex
        unfortuently doesnt work for dividends or split (yet)"""
        # THIS FUNCTION HONESTLY DOENST HAVE ANY CALCULATION PURPOSES
        # ITS SOLE PURPUSE IS PUTTING THE VALUES CALCULATED IN THE BINOMIAL MAP FUNCTION
        # INTO A TREE ON LATEX
        # AT THE TOP OF THE TABLES IS THE STOCK PRICE, MIDDLE IS THE CALL PRICE, BOTTOM IS THE PUT FUNCTION
        self.binomialMap(False)
        self.columns = {}
        self.columns2 = {}
        if self.number_of_periods < 4:
            self.st = """\\documentclass{article}\n\\usepackage{tikz}\n\\usetikzlibrary{matrix}\n\\usepackage{lscape}\n\\begin{document}\n\\begin{landscape}\n\\begin{tikzpicture}[>=stealth,sloped]\n\\matrix (tree) [%\nmatrix of nodes,\nminimum size=1cm,\ncolumn sep=3.5cm,\nrow sep=1cm,\n]\n{\n"""
        else:
            self.st = """\\documentclass{article}\n\\usepackage{tikz}\n\\usetikzlibrary{matrix}\n\\usepackage{lscape}\n\\begin{document}\n\\begin{landscape}\n\\begin{tikzpicture}[>=stealth,sloped]\n\\matrix (tree) [%\nmatrix of nodes,minimum size=0.5cm,\ncolumn sep=1cm,\nrow sep=0.5cm,]{\n"""
        self.st2 = ""
        for j in range_(-self.number_of_periods, 1, -1):
            for i in range_(0, j + 1):
                pass
        l = self.number_of_periods + 1
        for k in self.C.keys():
            i = k.count("U")
            j = k.count("D")
            try:
                self.columns[i + j].append("")
                self.columns2[i + j].append("")
            except:
                self.columns[i + j] = [
                    "" for l in range(self.number_of_periods - j - i)
                ]
                self.columns2[i + j] = [
                    "" for l in range(self.number_of_periods - j - i)
                ]
            self.columns2[i + j].append(self.C[k])
            self.columns[i + j].append(
                "U" * i + "D" * j + ("S0" if (j + i >= 0) else "")
            )
        self.C[""] = ""
        self.P[""] = ""
        for k in self.C.keys():
            i = k.count("U")
            j = k.count("D")
            for l in range(self.number_of_periods - j - i):
                self.columns[i + j].append("")
                self.columns2[i + j].append("")
        for i in range_(-2 * self.number_of_periods, 1, -1):
            for j in range(len(self.columns)):
                if j != len(self.columns) - 1:
                    try:
                        SSS = str(round(self.S[self.columns[j][i]], 3)) + "\\\\"
                        CCC = str(round(self.C[self.columns[j][i]], 3)) + "\\\\"
                        PPP = str(round(self.P[self.columns[j][i]], 3)) + "\\\\"
                        self.st += (
                            "\\begin{tabular}{ |c| } \n \hline\n"
                            + SSS
                            + CCC
                            + PPP
                            + "\\hline\\end{tabular} & "
                        )
                    except:
                        # print(self.columns[j][i])
                        self.st += self.columns[j][i] + " & "
                    if self.columns[j][i] != "":
                        U1, D1, U2, D2 = (
                            self.columns[j + 1][i - 1].count("D"),
                            self.columns[j + 1][i - 1].count("U"),
                            self.columns[j + 1][i + 1].count("D"),
                            self.columns[j + 1][i + 1].count("U"),
                        )
                        # adding the ligns and their values between nodes on the tree
                        self.st2 += (
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
                    try:
                        SSS = str(round(self.S[self.columns[j][i]], 3)) + "\\\\"
                        CCC = str(round(self.C[self.columns[j][i]], 3)) + "\\\\"
                        PPP = str(round(self.P[self.columns[j][i]], 3)) + "\\\\"
                        self.st += (
                            "\\begin{tabular}{ |c| } \n \hline\n"
                            + SSS
                            + CCC
                            + PPP
                            + "\\hline\\end{tabular} & "
                        )
                    except:
                        self.st += self.columns[j][i] + " & "

            self.st += " \\\\\n"
        self.st += (
            "};\n" + self.st2 + "\end{tikzpicture}\n\end{landscape}\n\end{document}"
        )
        return self.st

    @property
    def BS(self):
        return self.Black_Scholes()


sigma = sqrt(1 / 66)
number_of_periods = 4
deltaT = 0
T = 1
K = 38
S = 39
r = 1 / 100
v = r - 1 / 2 * sigma**2
op = Option(S, sigma, number_of_periods, T, K, v, r, h=0)
print(op.make_binomial_tree())


def P(st, k):
    return max(0, k - st)


def C(st, k):
    return max(st - k, 0)


def f(line):
    out = []
    for st in line:
        k1 = 40000
        k2 = 42000
        k3 = 44000
        out.append(C(st, k1) - 3 * C(st, k2) + 2 * C(st, k3))
    return out


line = np.linspace(35000, 50000, 10000)
plt.plot(line, f(line))
plt.show()
# f(line)

v = -0.133751
sigma = sqrt(0.550338)
number_of_periods = 12
T = 1
strike_prices = [40000, 42000, 44000]
S = 42600.132813
r = v + 1 / 2 * sigma**2
# v = r - 1 / 2 * sigma**2
ls = []
pb = []
w = [1, -3, 2]
greeks = []
for j in strike_prices:
    # print(j)
    op = Option(S, sigma, number_of_periods, T, j, log_variation=None, interest_rate=r)
    # make_binomial_tree(3)
    # print(op.BS)
    op.BS
    ls.append(op.call_price)
    pb.append(norm.cdf(op.D2))
    greeks.append(op.greeks_call())

delta, theta, gamma, vega, rho = 0, 0, 0, 0, 0
# for i,v in enumerate(greeks[0]):
for j, v in enumerate(greeks):
    delta += greeks[j][0] * w[j]
    theta += greeks[j][1] * w[j]
    gamma += greeks[j][2] * w[j]
    vega += greeks[j][3] * w[j]
    rho += greeks[j][4] * w[j]

delta, gamma, theta, vega, rho
# Part 2
num = 1
t = {
    0: 42600.132813,
    2: 37068.769531,
    4: 43194.503906,
    6: 47458.140625,
    8: 38059.902344,
    10: 39664.25,
    12: 43074.105469,
}
t


def delta_hedging(t, r, num, sigma, TF, K, div_factor=1):
    """THIS FUNCTION PERFORMS DELTA HEDGING AND RETURN IT IN A TABLE FORMAT"""
    data = []
    for T, S in t.items():
        try:
            # print((S, K, T / 12, r, sigma, TF))
            d = norm.cdf(D_10(S, K, T / 12, r, sigma, TF))
        except ZeroDivisionError:
            if S > K:
                d = 1
            else:
                d = 0
        try:
            change = d - data[-1][1]
            cost = round(change * num * S / div_factor, 3)
            tot_cost = round((data[-1][4] + data[-1][5] + cost), 2)
        except:
            change = d
            cost = round(d * num * S / div_factor, 3)
            tot_cost = cost
        ls = [
            S,
            round(d, 3),
            round(change * num, 3),
            cost,
            tot_cost / div_factor,
            (tot_cost * (r)) / div_factor,
        ]
        data.append(ls)

    df = pd.DataFrame(
        data,
        columns=[
            "Stock price",
            "Delta",
            "Shares purchased",
            "Costs",
            "Cumulative Costs",
            "Interest",
        ],
    )
    return df


def D_10(S0, K, T, r, s, TF):
    return 1 / (s * sqrt(TF - T)) * (np.log(S0 / K) + ((r + s**2 / 2) * (TF - T)))


r = 0.143523
TF = 1
num = 1

delta_hedging(t, r, num, sigma, TF, 40000)

delta_hedging(t, r, num, sigma, TF, 42000)

delta_hedging(t, r, num, sigma, TF, 44000)

df40000 = delta_hedging(t, r, num, sigma, TF, 40000)
CC40 = df40000["Cumulative Costs"][len(df40000) - 1]

df42000 = delta_hedging(t, r, num, sigma, TF, 42000)
CC42 = df42000["Cumulative Costs"][len(df42000) - 1]

df44000 = delta_hedging(t, r, num, sigma, TF, 44000)
CC44 = df44000["Cumulative Costs"][len(df44000) - 1]

CC40 - 3 * CC42 + 2 * CC44
