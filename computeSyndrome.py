import numpy as np
n = 4   #GF(2^4)
t = 2   #Error Correcting Capability
GF = {}

r = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0] #Received Vector

def generateField(n):
    GF[0] = 0b1
    for i in range(1, n):
        GF[i] = GF[i - 1] << 1
    for i in range(n, 2**n - 1):
        GF[i] = GF[i - n] ^ GF[i - n + 1]   #Primitive Polynomial = x^4 = x + 1

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

def computeSyndrome(r, x):
    h = []
    for i in range(0, 2**n - 1):
        h.append((x * i) % (2**n - 1))
    
    syndrome = 0
    for i in range(0, len(r)):
        if r[i] == 1:
            syndrome = syndrome ^ GF[h[i]]
    return syndrome

def main():
    generateField(n)
    for i in range(2*t):
        print(bin(computeSyndrome(r, i + 1)))

if __name__ == "__main__":
    main()
