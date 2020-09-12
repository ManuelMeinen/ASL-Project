import random
sk = bin(random.getrandbits(256))[2:]
pad = (256 - len(sk)) * '0'
sk = pad + sk
assert(len(sk)==256) 
if len(sk)==256:
    print(sk)

