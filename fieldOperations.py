import generateField
n = 4 #GF(2^8)

r = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]

def main():
    GF = generateField.field(n)
    print(GF)
    print(computeSyndrome(GF, r, 2))

def inverse(GF, num):       #Compute inverse of a field element
    for key, val in GF.items():
        if val == num:
            inv_key = ((2**n - 1) - key) % (2**n - 1)
            break
    inv = GF[inv_key]
    return inv
    

def fieldMul(GF, a, b):     #Multiply two elements in the field
    if a == 0 or b == 0:
        return 0
    for key, val in GF.items():
        if val == a:
            alpha_a = key
        if val == b:
            alpha_b = key
    p = GF[(alpha_a + alpha_b) % (2**n - 1)]
    return p

def computeSyndrome(GF, r, x):      #Compute the syndrome of a given vector
    h = []
    for i in range(0, 2**n - 1):
        h.append((x * i) % (2**n - 1))
    syndrome = 0
    for i in range(0, len(r)):
        if r[i] == 1:
            syndrome = syndrome ^ GF[h[i]]
    return syndrome


if __name__ == "__main__":
    main()
