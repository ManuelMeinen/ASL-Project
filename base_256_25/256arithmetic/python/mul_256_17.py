#!/usr/bin/env python
# coding: utf-8

# In[59]:


NO_LIMBS = 15


# In[61]:


prod = []
for k in range(NO_LIMBS*2 - 1):
    prod.append([])

#First half
for k in range(NO_LIMBS-1,-1,-1):
    i = 0
    j = k
    l = []
    while(j >= 0):
        l.append((i,j))
        i += 1
        j -= 1
    prod[k] = l

#second half
for k in range(1, NO_LIMBS):
    i = k
    j = NO_LIMBS - 1
    code_str = ""
    l = []
    while(j>=k):
        l.append((i,j))
        j -= 1
        i += 1
    prod[k+NO_LIMBS-1] = l
print(prod)
    


# In[62]:


def gen_code(prod):
    i = 0
    c = ""
    for tuples in prod:
        l = gen_line(tuples, i)
        c += l
        i += 1
    return c


# In[63]:


def gen_line(tuples, i):
    s = "r[%d] = " % (i)
    for t in tuples:
        s += " a[%d] * b[%d] +" % (t[0],t[1])
    return s[:-2] + ";\n"
        


# In[64]:


print(gen_code(prod))

