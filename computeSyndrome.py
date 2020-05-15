import numpy as np
n = 4
GF = {
    0: 0b0001,
    1: 0b0010,
    2: 0b0100,
    3: 0b1000,
    4: 0b0011,
    5: 0b0110,
    6: 0b1100,
    7: 0b1011,
    8: 0b0101,
    9: 0b1010,
    10: 0b0111,
    11: 0b1110,
    12: 0b1111,
    13: 0b1101,
    14: 0b1001
}

def fieldMul(a, b):
    p = 0b0
    for _ in range(n):
        if b % 2 == 1:
            p ^= a
        b >>= 1
        a <<= 1
    return reducedForm(p)

def reducedForm(num):
    numGF = num & 0b1111
    binArr = list(bin(num)[2: ])
    binArr = binArr[0: len(binArr) - 4]
    binArr.reverse()
    for i in range(n, n + len(binArr)):
        if(binArr[n - i] == '1'):
            numGF ^= GF[i]
    return numGF

a = 0b0101
b = 0b1100
print(fieldMul(a, b))