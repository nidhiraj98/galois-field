import numpy as np
n = 4   #GF(2^4)
t = 3   #Error Correcting Capability
GF = {}

r = [0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0] #Received Vector

def generateField(n):
    GF[0] = 0b1
    for i in range(1, n):
        GF[i] = GF[i - 1] << 1
    for i in range(n, 2**n - 1):
        GF[i] = GF[i - n] ^ GF[i - n + 1]   #Primitive Polynomial = x^4 = x + 1

def inverse(num):
    num = [str(i) for i in num]
    num = int("".join(num), 2)
    for key, val in GF.items():
        if val == num:
            inv_key = (15 - key) % 15
            break
    inv = GF[inv_key]
    return inv
    

def fieldMul(a, b):
    p = 0b0
    for _ in range(n):
        if b % 2 == 1:
            p ^= a
        b >>= 1
        a <<= 1
    # print(p)
    return reducedForm(p)

def reducedForm(num):
    if len(bin(num)[2: ]) <= 4:
        return num
    numGF = num & 0b1111
    binArr = list(bin(num)[2: ])
    binArr = binArr[0: len(binArr) - 4]
    binArr.reverse()
    for i in range(n, n + len(binArr)):
        if(binArr[n - i] == '1'):
            numGF ^= GF[i]
    # print("Field Mul = ", numGF)
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
    # print(fieldMul(6, 12))
    syndromes = []
    for i in range(2*t):
        syndromes.append(computeSyndrome(r, i + 1))
    print(syndromes)

    errLocTable = {
        'mu': [-0.5] + [i for i in range(0, t + 1)],
        'sigma': [[1], [1]],
        'd': [[1], [syndromes[0]]],
        'l': [0, 0],
        'diff': [-1, 0]
    }

    for i in range(2, t + 1):
        mu = errLocTable['mu'][i - 1]
        # print("mu = ", mu)
        if errLocTable['d'][i - 1] == 0:
            sigma = errLocTable['sigma'][i - 1]
        else:
            p = errLocTable['mu'][errLocTable['diff'].index(max(errLocTable['diff'][:i - 1]))]
            # print("p = ", p)
            d_u = errLocTable['d'][i - 1]
            d_u = [str(j) for j in d_u]
            d_u = int("".join(d_u), 2)
            # print(d_u)
            d_p_inv = inverse(errLocTable['d'][errLocTable['mu'].index(p)])
            # print(d_p_inv)
            sigma_p = errLocTable['sigma'][errLocTable['mu'].index(p)]
            # print(sigma_p)
            big_term = list(np.polymul([fieldMul(d_u, d_p_inv)], sigma_p)) + [0 for _ in range(int(2 * (mu - p)))]
            # print(big_term)
            sigma = list(np.polyadd(big_term, errLocTable['sigma'][i - 1])) #add x^2*(mu - p) #dry run here
        l = len(sigma) - 1
        diff = 2 * errLocTable['mu'][i] - l
        x = len(sigma)
        d = list(bin(syndromes[2 * mu + 3 - 1])[2:])
        d = [int(j) for j in d]
        for j in range(1, x):
            syn = list(bin(syndromes[2 * mu + 3 - 1 - l])[2: ])
            syn = [int(k) for k in syn]
            d = np.polyadd(d, np.polymul(syn, [sigma[x - j]]))
            d = [int(k % 2) for k in d]
        errLocTable['sigma'].append(sigma)
        errLocTable['d'].append(d)
        errLocTable['l'].append(l)
        errLocTable['diff'].append(diff)
    errLocPoly = errLocTable['sigma'][t]
    print(errLocPoly)

    # print(errLocTable)
if __name__ == "__main__":
    main()
