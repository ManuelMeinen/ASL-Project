import sys 

LIMB_SIZE_1 = 26
LIMB_SIZE_2 = 25
NO_LIMBS = int(255/LIMB_SIZE_1)
def get_limbs25(x):
  x = str(x)
  x = str(bin(int(x,10)))[2:].zfill(255)
  l = []
  i = 0
  v0 = int(x[0:25],2) #25
  v1 = int(x[25:51],2) #26
  v2 = int(x[51:76],2) #25
  v3 = int(x[76:102],2) #26
  v4 = int(x[102:127],2) #25
  v5 = int(x[127:153],2) #26
  v6 = int(x[153:178],2) #25
  v7 = int(x[178:204],2) #26
  v8 = int(x[204:229],2) #25
  v9 = int(x[229:],2) #26
  l.append(v0)
  l.append(v1)
  l.append(v2)
  l.append(v3)
  l.append(v4)
  l.append(v5)
  l.append(v6)
  l.append(v7)
  l.append(v8)
  l.append(v9)

  l.reverse()
  return l


def check(l, a):
  bstr = ""
  i = 0
  p = l[0] + (2**26) * l[1] + (2**51) * l[2] + (2**77)*l[3] + (2**102)*l[4] + (2**128)*l[5] + (2**153)*l[6] + (2**179)*l[7] + (2**204)*l[8] + (2**230)*l[9]
  print(p)
  print(a)

def get_bin_from_r17(l):
  bstr = ""
  for i in range(NO_LIMBS-1,-1,-1):  
    bstr = bstr + (str(bin(l[i]))[2:]).zfill(17)
  return bstr

if __name__ == '__main__':
  args = sys.argv
  if len(args) < 2:
    print("No arguments passed.\nPass a list of space separated decimal numbers to convert to radix 25.5 limbs.")
    exit()

  for i in range(1,len(args)):
    dec_num = args[i]
    print(get_limbs25(dec_num))

# a_bin = get_bin_from_r17([93848, 98738, 33478, 93614, 28872, 94227, 67689, 56183, 5696, 49829, 67138, 101447, 15544, 71231, 30135])
# l_a = get_limbs25(a_bin.zfill(255))

# b_bin = get_bin_from_r17([93378, 87592, 99224, 129887, 91690, 114010, 62066, 58051, 58184, 38652, 92295, 55055, 40636, 33242, 25333])
# l_b = get_limbs25(b_bin.zfill(255))

# resa_bin = get_bin_from_r17([56154, 55259, 1631, 92430, 120563, 77165, 129756, 114234, 63880, 88481, 28361, 25431, 56181, 104473, 55468])
# l_resa = get_limbs25(resa_bin.zfill(255))

# ress_bin = get_bin_from_r17([470, 11146, 65326, 94798, 68253, 111288, 5622, 129204, 78583, 11176, 105915, 46391, 105980, 37988, 4802])
# l_ress = get_limbs25(ress_bin.zfill(255))
# print("l_a", l_a)
# print("l_b", l_b)
# print("l_resa", l_resa)
# print("l_ress", l_ress)

# print(get_limbs25(str(bin(57896044618658097711785492504339381575729535664163146648190498831576564236739))[2:]))

# print(get_limbs25(str(bin(48674846124269517316669609744293450289806781891760030838356554460966925028204))[2:]))

# check(l_a, int(a_bin,2))
# check(l_b, int(b_bin,2))
# check(l_resa, int(resa_bin,2))
