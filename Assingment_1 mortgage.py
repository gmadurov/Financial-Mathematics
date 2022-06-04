import numpy as np
import pandas as pd
from scipy.optimize import fsolve

interest = 0.0184 / 12
borrow = 250000
redemption = 904

interest1 = 0.0184 / 12
borrow1 = 250000
redemption1 = 1072


def linear_mortgage(interest, borrow, installment):
    arr = np.zeros((4, 2))
    arr[0, 0] = borrow
    redemption = 250000 / 360
    arr[2, 0] = 1072 - redemption
    i = 0
    continu = True
    while continu:

        try:
            arr[0, i + 1]
        except IndexError:
            z = np.array([[0] for i in range(len(arr))])
            arr = np.append(arr, z, axis=1)
        arr[1, i] = arr[0, i] * interest
        arr[2, i] = min(arr[0, i], redemption)
        arr[3, i] = arr[1, i] + arr[2, i]
        try:
            arr[0, i + 1] = arr[0, i] - arr[2, i]
            if arr[0, i] <= redemption:
                continu = False
        except:
            pass
        i += 1
    array = np.append(np.array(arr[:, :5]), np.array(arr[:, -5:]), axis=1)
    print(arr[:, :5], arr[:, -5:])
    df = pd.DataFrame(
        array, columns=["1", "2", "3", "4", "5", "357", "358", "359", "360", "361"]
    )
    df.index = ["Remainder", "Interest", "Redemption", "Payment"]
    print(df)
    # with open('linear_mortgage.txt', 'w') as f:
    #     f.write(df.to_latex())


def annuity_mortgage(interest, borrow, annuity):
    arr = np.zeros((4, 2))
    arr[0, 0] = borrow
    i = 0
    continu = True
    for r in range(360):
        try:
            arr[0, i + 1]
        except IndexError:
            z = np.array([[0] for i in range(len(arr))])
            arr = np.append(arr, z, axis=1)
        arr[1, i] = round(arr[0, i] * interest, 9)
        arr[3, i] = min(arr[0, i] + arr[1, i], annuity)
        arr[2, i] = arr[3, i] - arr[1, i]
        try:
            arr[0, i + 1] = arr[0, i] - arr[2, i]
            if arr[0, i] <= annuity:
                continu = False
        except:
            pass
        i += 1
    df = pd.DataFrame(arr, columns=[f"{i}" for i in range(1, len(arr[0]) + 1)])
    print(df)


def part_4(interest, borrow, annuity):
    # 500(1 +r/12)^12−i
    amount = 500 * sum([(1 + interest / 12) ** (12 - i) for i in range(1, 13)])
    print(amount)
    arr = np.zeros((4, 2))
    arr[0, 0] = borrow
    i = 0
    continu = True
    while continu:
        try:
            arr[0, i + 1]
        except IndexError:
            z = np.array([[0] for i in range(len(arr))])
            arr = np.append(arr, z, axis=1)
        arr[1, i] = round(arr[0, i] * interest, 2)
        arr[3, i] = min(arr[0, i] + arr[1, i], annuity)
        arr[2, i] = arr[3, i] - arr[1, i]
        if i % 12 == 0 and i > 0:
            arr[2, i] += amount
        try:
            arr[0, i + 1] = arr[0, i] - arr[2, i]
            if arr[0, i] <= annuity:
                continu = False
        except:
            pass
        i += 1
    df = pd.DataFrame(arr, columns=[f"{i}" for i in range(1, len(arr[0]) + 1)])
    print(sum(arr[3, :]) / borrow)
    print(df)
    # with open('ouput.txt', 'w') as file:
    #     file.write(df.to_latex())


print(
    linear_mortgage(interest1, borrow1, redemption1),
    annuity_mortgage(interest, borrow, redemption),
)


def annuity_mortgage(r):
    PV = (904 / (r / 12)) * (1 - (1 / (1 + (r / 12)) ** (30 * 12)))
    return PV - 250000  # net PV


print(fsolve(annuity_mortgage, 0.0184)[0])


def linear_mortgage(r):
    PV = sum([1072 / (1 + (r / 12)) ** (12 * (i)) for i in range(1, 31)])
    return PV - 250000  # net PV


print(fsolve(linear_mortgage, 0.0184)[0])


def part_4(interest, borrow, annuity, m):
    """dfd is the final tablau of calculations"""
    # 500(1 +r/12)^12−i
    amount = 500 * sum([(1 + interest / 12) ** (12 - i) for i in range(1, 13)])
    # print(amount)
    arr = np.zeros((4, 2))
    arr[0, 0] = borrow
    i = 0
    continu = True
    while continu:
        # for r in range(360):
        try:
            arr[0, i + 1]
        except IndexError:
            z = np.array([[0] for i in range(len(arr))])
            arr = np.append(arr, z, axis=1)
        arr[1, i] = round(arr[0, i] * interest / m, 2)
        arr[3, i] = min(arr[0, i] + arr[1, i], annuity)
        arr[2, i] = arr[3, i] - arr[1, i]
        if i % 12 == 0 and i > 0:
            arr[2, i] += amount
        try:
            arr[0, i + 1] = arr[0, i] - arr[2, i]
            if arr[0, i] <= annuity:
                continu = False
        except:
            pass
        i += 1
    array = np.append(np.array(arr[:, :5]), np.array(arr[:, -5:]), axis=1)
    dfd = pd.DataFrame(arr, columns=[f"{i}" for i in range(1, len(arr[0]) + 1)])
    # df = pd.DataFrame(array, columns=[f'{i}' for i in range(1, len(array[0])+1)])
    # df.index = ['Remainder','Interest','Redemption','Payment']

    dfd.index = ["Remainder", "Interest", "Redemption", "Payment"]
    # print(df)
    print(dfd)
    # print(sum(arr[3, :]))
    # print(amount)
    # with open('linear_mortgage.txt','w') as f:
    #     f.write(df.to_latex())


part_4(0.0184, 250000, 904, 12)
