import sys 

LIMB_SIZE = 17
NO_LIMBS = int(255/LIMB_SIZE)
def get_limbs17(x):
  l = []
  for i in range(0,NO_LIMBS-1):
    #print(LIMB_SIZE*i,LIMB_SIZE*(i+1))
    v = int(x[LIMB_SIZE*i:LIMB_SIZE*(i+1)],2)
    l.append(v)
  #print(len(x[LIMB_SIZE*(NO_LIMBS-1):256]))
  v = int(x[LIMB_SIZE*(NO_LIMBS-1):256],2)
  l.append(v)
  l.reverse()
  return l

def check(l, a):
  bstr = ""
  for i in range(NO_LIMBS-1,-1,-1):  
    bstr = bstr + (str(bin(l[i]))[2:]).zfill(LIMB_SIZE)
  print(int(bstr,2))
  print(int(a,2))


args = sys.argv
if len(args) < 2:
    print("No arguments passed.\nPass a list of space separated decimal numbers to convert to radix 17 limbs.")
    exit()

for i in range(1,len(args)):
    dec_num = args[i]
    print(get_limbs17((str(bin(int(dec_num,10)))[2:]).zfill(255)))