#%%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt, exp, factorial
import pprint
from scipy.stats import norm
from string import Template


def range_(start=0, end=1, factor=1, step=1):
    return [i * factor for i in range(start, end, step)]


def bicoef(top, bottom):
    return factorial(top) / (factorial(bottom) * factorial(top - bottom))


#%%
class Option:
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
        self.S = {"": "", 0: S0}
        self.v = log_variation or interest_rate - 1 / 2 * sigma**2
        self.sigma = sigma
        self.number_of_periods = number_of_periods
        self.deltaT = deltaT or (1 / self.number_of_periods)
        self.T = Time_horizon
        self.K = strike_prices
        self.excercise_matrixCall = {}
        self.excercise_matrixPut = {}
        self.C = {}
        self.P = {}
        self.probCall = 0
        self.probPut = 0
        self.h = h
        self.r = interest_rate or self.v + 1 / 2 * self.sigma**2
        self.get_values()

    def get_values(self):
        self.R = exp(self.r * (self.deltaT))
        self.u = exp(self.sigma * sqrt(self.deltaT))
        self.d = 1 / self.u
        self.p = 1 / 2 + 1 / 2 * self.v / self.sigma * sqrt(self.deltaT)
        self.q = (self.R - self.d) / (self.u - self.d)

    def binomialMap(self, show=True):
        for j in range_(-self.number_of_periods, 1, -1):
            for i in range(0, j + 1):
                self.S["U" * i + "D" * (j - i) + "S0"] = max(
                    0, (self.S[0] * (self.u**i) * (self.d ** (j - i)))
                )
        """initial RHS of the C tree"""
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
        """returns --> deltaCall, thetaCall, gammaCall, vegaCall"""
        self.Black_Scholes()
        self.deltaCall = norm.cdf(self.D1)
        self.thetaCall = (self.S[0] * norm.cdf(self.D1) * self.sigma) / (
            2 * sqrt(T)
        ) - r * self.K * exp(-self.r * self.T * norm.cdf(self.D2))
        self.gammaCall = norm.cdf(self.D1) / (self.S[0] * self.sigma * sqrt(T))
        self.vegaCall = self.S[0] * sqrt(T) * self.sigma * norm.cdf(self.D1)
        self.rho = self.K * self.T * exp(-self.r * self.T) * norm.cdf(self.D2)
        if show:
            print(
                f"Delta = {self.deltaCall},\nTheta = {self.thetaCall},\nGamma = {self.gammaCall},\nVega = {self.vegaCall},\nRho = {self.rho}"
            )
        return self.deltaCall, self.thetaCall, self.gammaCall, self.vegaCall, self.rho

    def make_binomial_tree(self):
        self.binomialMap(False)
        self.columns = {}
        self.columns2 = {}
        if not self.h > 0:
            if self.number_of_periods < 4 and self.h == 0:
                self.st = "\\documentclass{article}\n \\usepackage{tikz}\n \\usetikzlibrary{matrix}\n \\usepackage{lscape}\n \\begin{document}\n \\begin{landscape}\n \\begin{tikzpicture}[>=stealth,sloped]\n \\matrix (tree) [%\n matrix of nodes,\n minimum size=1cm,\n column sep=3.5cm,\n row sep=1cm,\n]\n{"
            else:
                self.st = "\\documentclass{article}\n \\usepackage{tikz}\n \\usetikzlibrary{matrix}\n \\usepackage{lscape}\n \\begin{document}\n \\begin{landscape}\n \\begin{tikzpicture}[>=stealth,sloped]\n \\matrix (tree) [%\n matrix of nodes,\n minimum size=0.5cm,\n column sep=1cm,\n row sep=0.5cm,\n]\n{\n"
            self.st2 = ""
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
        else:
            self.st = "\\documentclass{article}\n\\usepackage{tikz}\n\\begin{document}\n% End of automatically added code\n% Set the overall layout of the tree\n\\tikzstyle{level 1}=[level distance=4cm, sibling distance=3.5cm,->]\n\\tikzstyle{level 2}=[level distance=4cm, sibling distance=2cm,->]\n% Define styles for bags and leafs\n\\tikzstyle{bag} = [text width=2em, text centered]\n\\tikzstyle{end} = []\n% The sloped option gives rotated edge labels. Personally\n% I find sloped labels a bit difficult to read. Remove the sloped options\n% to get horizontal labels. \n\\begin{tikzpicture}[grow=right, sloped]\n"
            self.st2 = ""
            l = self.number_of_periods + 1
            for k in self.C.keys():
                i = k.count("U")
                j = k.count("D")
                try:
                    self.columns[i + j].append(k)
                except:
                    self.columns[i + j] = [k]
                    self.columns2[i + j] = [k]
            children = {}
            # for k in self.columns.keys():
            #     print(self.columns[k])
            #     for i, num in enumerate(self.columns[k]):
            #         if not k+1 == 3:
            #             for ind in range(i+2):
            #                 print(k+1,ind, self.columns[k + 1][ind])
            for k in self.columns.keys():
                for i, num in enumerate(self.columns[k]):
                    if k == self.number_of_periods:
                        # print(1, k)
                        U1, D1 = self.columns[k][i].count("U"), self.columns[k][
                            i
                        ].count("D")
                        children[self.columns[k][i]] = Template(
                            f"\n"
                            + "\t" * k
                            + f'node[$node $child $end {"{"}{"}"}\n'
                            + "\t" * k
                            + f"edge from parent\n"
                            + "\t" * k
                            + f'node[above] {"{"}{"}"}\n'
                            + "\t" * k
                            + f'node[below] {"{"}$$ {"p" if U1 ==1 else f"p^{U1}" if U1 > 0 else ""} {"(1-p)" if D1 == 1 else f"(1-p)^{D1}" if D1 > 0 else ""} $${"}"}'
                        )
                    else:
                        for ind in range(i, i + 1):
                            # children[num] = (ind, num, self.columns[k+1][ind],self.columns[k+1][ind+1])
                            # print(2,k, ind)
                            U1, D1 = self.columns[k + 1][ind].count("U"), self.columns[
                                k + 1
                            ][ind].count("D")
                            children[self.columns[k][ind]] = Template(
                                f'$node \nchild {"{"}'
                                + "\t" * k
                                + f"$child {'}'}"
                                + "$arrow"
                                # + "\t" * k
                                # + f"\nedge from parent\n"
                                # + "\t" * k
                                # + f'node[above] {"{"}{"}"}\n'
                                # + "\t" * k
                                # + f'node[below] {"{"}$$ {"p" if U1 ==1 else f"p^{U1}" if U1 > 0 else ""} {"(1-p)" if D1 == 1 else f"(1-p)^{D1}" if D1 > 0 else ""} $${"}"};\n\t{"}"}'
                            )
            childs = {}
            for k in self.columns.keys():
                for i, num in enumerate(self.columns[k]):
                    if not k + 1 == 3:
                        for ind in range(i, i + 2):
                            try:
                                # childs[self.columns[k][i]] += (children[self.columns[k+1][ind]].substitute(child=self.columns[k+1][ind]))
                                u += 1
                            except:
                                pass
                                # childs[self.columns[k][i]] = (children[self.columns[k+1][ind]].substitute(child=self.columns[k+1][ind]))
                            # print(u, k, i,num, ind)
            # print(self.columns)
            # for i in children.values():
            #     print(i.template)

            try:
                for k in self.columns.keys():
                    for i, num in enumerate(self.columns[k]):
                        if k == self.number_of_periods:
                            phrase = children[self.columns[k][i]].substitute(
                                node="end, label=right:{",
                                child=("$" + self.columns[k][i] + "$"),
                                end="}]",
                            )
                            children[self.columns[k][i]] = phrase
                        elif k == 0:
                            U1, D1 = self.columns[0][i].count("U"), self.columns[k][
                                i
                            ].count("D")
                            phrase = children[self.columns[0][0]].substitute(
                                node=("\\node[bag] {$" + self.columns[0][0] + "$}"),
                                child=(children[self.columns[1][0]] + " "),
                                arrow="",
                            )
                            phrase += children[self.columns[0][0]].substitute(
                                node="",
                                child=(children[self.columns[1][1]] + " "),
                                arrow=";",
                            )
                            children[self.columns[0][0]] = phrase
                        else:
                            for ind in range(i, i + 2):
                                # ind = 1
                                # print(k, i, ind)
                                U1, D1 = self.columns[k][i].count("U"), self.columns[k][
                                    i
                                ].count("D")
                                try:
                                    # print(

                                    #     self.columns[k][i])
                                    # print(
                                    #     children[self.columns[k][i]])
                                    # print(
                                    #     self.columns[k][ind])
                                    # print(self.columns[k + 1][ind])
                                    # print(children[self.columns[k + 1][ind]])
                                    # print(k,i, children[self.columns[k][i]])
                                    phrase = children[self.columns[k][ind]].substitute(
                                        node=(
                                            "node[bag] {$" + self.columns[k][ind] + "$}"
                                        ),
                                        child=(
                                            children[self.columns[k + 1][ind]] + " "
                                        ),
                                        arrow="",
                                    )
                                    phrase += children[
                                        self.columns[k][i + 1]
                                    ].substitute(
                                        node="",
                                        child=(
                                            children[self.columns[k + 1][ind + 1]] + " "
                                        ),
                                        arrow=f"\nedge from parent\n"
                                        # + "\t" * k
                                        + f'node[above] {"{"}{"}"}\n'
                                        # + "\t" * k
                                        + f'node[below] {"{"}$ {"p" if U1 ==1 else f"p^{U1}" if U1 > 0 else ""} {"(1-p)" if D1 == 1 else f"(1-p)^{D1}" if D1 > 0 else ""} ${"}"};\n\t',
                                    )
                                    children[self.columns[k][ind]] = phrase
                                    # print(phrase)
                                except:
                                    pass

                                    # phrase = children[
                                    #     self.columns[k][i - 1]
                                    # ].substitute(
                                    #     node=(
                                    #         "node[bag] {$" + self.columns[k][i] + "$}"
                                    #     ),
                                    #     child=(children[self.columns[k][i]] + " "),
                                    #     arrow="",
                                    # ) + children[
                                    #     self.columns[k][i]
                                    # ].substitute(
                                    #     node="",
                                    #     child=(children[self.columns[k][i]] + " "),
                                    #     arrow=f"\nedge from parent\n"
                                    #     # + "\t" * k
                                    #     + f'node[above] {"{"}{"}"}\n'
                                    #     # + "\t" * k
                                    #     + f'node[below] {"{"}$ {"p" if U1 ==1 else f"p^{U1}" if U1 > 0 else ""} {"(1-p)" if D1 == 1 else f"(1-p)^{D1}" if D1 > 0 else ""} ${"}"};\n\t{"}"}',
                                    # )
                                finally:
                                    pass
                                    # print(
                                    #     "error",
                                        # self.columns[k][i],
                                        # children[self.columns[k][i]].template,
                                    # )
                                # phrase = children[self.columns[k][ind]].substitute(
                                #     node=("node[bag] {$" + self.columns[k][ind] + "$}"),
                                #     child=children[self.columns[k + 1][ind]],
                                # )
                                # phrase += children[self.columns[k][ind + 1]].substitute(
                                #     node="",
                                #     child=children[self.columns[k + 1][ind + 1]],
                                # )
                                # print(phrase)
                                # children[self.columns[k][ind]] = phrase
                            # print(u, k, i,num, ind)

            finally:
                pass
                # for k, v in children.items():
                #     print(k)
                #     if type(v) == Template:
                #         print(v.template)
                #     else:
                #         print(v)

                # print("###############################")
            # print(children["S0"])
            # print(self.columns)
            self.st += children["S0"]
            self.st += (
                "\n\\end{tikzpicture}\n% Automatically added code\n\\end{document}"
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
# op.make_binomial_tree()
# print(op.BS)
# op = Option(S, sigma, v, number_of_periods, T, K, r, h=0)
print(op.make_binomial_tree())
# print(op.BS)

# DS0
# child {
# 	node[end, label=right:{ $DDS0$ }] {}
# 	edge from parent
# 	node[above] {}
# 	node[below] {$  (1-p)^2 $};
# 	}
# 	 }] {}
# 	edge from parent
# 	node[above] {}
# 	node[below] {$  (1-p) $};
# 	}
# child {
# 	node[bag, label=right:{ child {
# 	node[end, label=right:{ $UDS0$ }] {}
# 	edge from parent
# 	node[above] {}
# 	node[below] {$ p (1-p) $};
# 	}
# 	 }] {}
# 	edge from parent
# 	node[above] {}
# 	node[below] {$ p  $};
# 	}

# %%
