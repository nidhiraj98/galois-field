n = 4 #GF(2^4)
GF = {}

def field(n):
    GF[0] = 0b1
    for i in range(1, n):
        GF[i] = GF[i - 1] << 1
    for i in range(n, 2**n - 1):
        # GF[i] = GF[i - n] ^ GF[i - n + 2] ^ GF[i - n + 3] ^ GF[i - n + 4]   #GF(2^8)
        # GF[i] = GF[i - n] ^ GF[i - n + 3] #GF(2^7)
        # GF[i] = GF[i - n] ^ GF[i - n + 1] #GF(2^6)
        # GF[i] = GF[i - n] ^ GF[i - n + 2] #GF(2^5)
        GF[i] = GF[i - n] ^ GF[i - n + 1] #GF(2^4)
        # GF[i] = GF[i - n] ^ GF[i - n + 1] #GF(2^3)
    return GF

def main():
    print(field(n))

if __name__ == "__main__":
    main()
