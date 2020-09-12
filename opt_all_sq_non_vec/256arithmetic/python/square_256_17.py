#!/usr/bin/env python
# coding: utf-8

# In[49]:


NO_LIMBS = 15


# In[50]:


i
prod = []
for k in range(NO_LIMBS*2 - 1):
    prod.append([])

#First half
for k in range(NO_LIMBS-1,-1,-1):
    l_a = []
    l_2a = []
    i = 0
    j = k
    if k%2 == 0:
        l_a.append((int(k/2),int(k/2)))
    while (i < math.ceil(k/2)):
        l_2a.append((i,j))
        i += 1
        j -= 1
    prod[k] = (l_2a,l_a)

#Second half
for k in range(1, NO_LIMBS):
    i = k
    j = NO_LIMBS - 1
    l_a = []
    l_2a = []
    
    if (k+NO_LIMBS-1)%2 == 0:
        l_a.append((int((k+NO_LIMBS-1)/2), int((k+NO_LIMBS-1)/2)))
    
    while (i<j):
        l_2a.append((i,j))
        i += 1
        j -= 1
    prod[k+NO_LIMBS-1] = (l_2a,l_a)

print(prod)
    


# In[51]:


def gen_code(prod):
    code = ""
    for i in range(0, NO_LIMBS - 1):
        code += "_2a[%d] = a[%d] << 1;\n" % (i,i)
    code += "\n"
    
    i = 0
    for p in prod:
        code += gen_line(p[0],p[1],i)
        i += 1
    return code


# In[52]:


def gen_line(l_2a,l_a, i):
    l = "r[%d] = " %(i)
    for x in l_2a:
        l += "_2a[%d]*a[%d] + " % (x[0],x[1])
    
    for x in l_a:
        l += "a[%d]*a[%d] + " % (x[0],x[1])
    
    return l[:-3] + ";\n"


# In[53]:


print(gen_code(prod))

