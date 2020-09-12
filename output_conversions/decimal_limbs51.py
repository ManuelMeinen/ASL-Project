import sys 

def get_limbs51(x):
  l = []
  for i in range(0,4):
    v = int(x[51*i:51*(i+1)],2)
    l.append(v)
  v = int(x[51*4:256],2)
  l.append(v)
  l.reverse()
  return l

args = sys.argv
if len(args) < 2:
    print("No arguments passed.\nPass a list of space separated decimal numbers to convert to radix 51 limbs.\n\nSample input:\npython get_limbs51.py 23423465588619757183253169998465881899659595894318503666586588616682836581719 13311224122668955738146552185873637781933051739291637575097427091064598392472 11189996199370894582589062329635885928504840703623090769010532534865250839746 \n \nSample output:\n[664433579777367, 798145316342188, 1448017419764148, 1540225854329050, 911028651365221\n[575160602422936, 1618811318005166, 97863898957929, 1742854989071013, 517724694265016]\
\n[1704666820865218, 1958688903789407, 999601117524594, 945849795254012, 435221983174332]")
    exit()

for i in range(1,len(args)):
    dec_num = args[i]
    print(get_limbs51((str(bin(int(dec_num,10)))[2:]).zfill(255)))