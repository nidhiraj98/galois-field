import numpy as np
n = 4   #GF(2^4)
t = 3   #Error Correcting Capability
GF = {}

r = [1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0] #Received Vector

def generateField(n):
    GF[0] = 0b1
    for i in range(1, n):
        GF[i] = GF[i - 1] << 1
    for i in range(n, 2**n - 1):
        GF[i] = GF[i - n] ^ GF[i - n + 1]   #Primitive Polynomial = x^4 = x + 1

def inverse(num):
    # num = [str(i) for i in num]
    # num = int("".join(num), 2)
    for key, val in GF.items():
        if val == num:
            inv_key = (15 - key) % 15
            break
    inv = GF[inv_key]
    return inv
    

def fieldMul(a, b):
    if a == 0 or b == 0:
        return 0
    # print(a, b)
    for key, val in GF.items():
        if val == a:
            alpha_a = key
        if val == b:
            alpha_b = key
    p = GF[(alpha_a + alpha_b) % (2**n - 1)]
    # print("Field Mul = ", p)
    return p

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
    try:
        generateField(n)
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

        for i in range(2, t + 2):
            mu = errLocTable['mu'][i - 1]

            if errLocTable['d'][i - 1] == 0: #Compute sigma
                sigma = errLocTable['sigma'][i - 1]
            else:
                p = errLocTable['mu'][errLocTable['diff'].index(max(errLocTable['diff'][:i - 1]))]
                # print("p = ", p)
                d_u = errLocTable['d'][i - 1]
                d_u = [str(j) for j in d_u]
                d_u = int("".join(d_u), 2)
                # print("du = ", d_u)
                d_p_inv = errLocTable['d'][errLocTable['mu'].index(p)]
                d_p_inv = [str(j) for j in d_p_inv]
                d_p_inv = inverse(int("".join(d_p_inv), 2))
                # print("dp = ", d_p_inv)
                sigma_p = errLocTable['sigma'][errLocTable['mu'].index(p)]
                # print("sigma_p = ", sigma_p)
                big_term = list(np.polymul([fieldMul(d_u, d_p_inv)], sigma_p)) + [0 for _ in range(int(2 * (mu - p)))]
                # print("bigterm = ", big_term)
                sigma_u = errLocTable['sigma'][i - 1]
                x = len(big_term) - len(sigma_u)
                sigma_u = [0 for _ in range(x)] + sigma_u
                sigma = list(np.bitwise_xor(big_term, sigma_u))

            l = len(sigma) - 1
            diff = 2 * errLocTable['mu'][i] - l

            errLocTable['sigma'].append(sigma)
            errLocTable['l'].append(l)
            errLocTable['diff'].append(diff)

            if(i == t + 1):
                break
            x = len(sigma)
            d = list(bin(syndromes[2 * mu + 3 - 1])[2:]) #initialize d
            d = [int(j) for j in d]
            for j in range(1, x):   #compute d
                syn = list(bin(syndromes[2 * mu + 3 - 1 - l])[2: ])
                syn = [int(k) for k in syn]
                d = np.polyadd(d, np.polymul(syn, [sigma[x - j]]))
                d = [int(k % 2) for k in d]
            errLocTable['d'].append(d)

        # print(errLocTable)

        errLocPoly = errLocTable['sigma'][t + 1]
        # print(errLocPoly)
        errLocPoly.reverse()
        errLoc = []
        for i in range(0, 2**n - 1):
            sum = 0
            for j in range(0, len(errLocPoly)):
                # print("i, j = ", i, j)
                sum ^= fieldMul(errLocPoly[j], GF[(i * j) % (2**n - 1)])
            if sum == 0:
                errLoc.append((2**n - 1) - i)

        print(errLoc)
        for i in errLoc:
            r[i] ^= 1
        print("Message Sent: ", r)
    except:
        print("Message corrupted. Request Retransmission")
if __name__ == "__main__":
    main()
