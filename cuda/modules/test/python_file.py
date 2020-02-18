import os
import sys
sys.path.append(os.path.join("..", "..", "..", "build"))
from pycliff.metric import EuclideanMetric
from pycliff.multivector import *
#
DIMS = 10
Multivector.set_N(DIMS)
generate_T(EuclideanMetric())
#
print((e(1)))

# print(CONJUGATE(1+e(1)+3*(e(1)^e(2))))
# print

# T1 = 0
# T2 = 0
# for i in range(DIMS):
    # T1 = T1 + e(i)
# for i in range(DIMS):
    # T2 = T2 + e(i)
# print(GP(T1,T2))
# T = get_T()
# T_OP = extract_OP_tensor(Multivector.get_N(), T)

# print(GP(e(0), e(1)^e(2)))

# l = [e(0), e(1), e(2), e(3), e(1)^e(2), e(1)^e(3), e(2)^e(3), e(1)^e(2)^e(3)]

# a = GP((e(1)+ e(3)),(e(1)+e(2)))
# print(a)

# for i in range(len(l)):
	# for j in range(len(l)):
		# print(l[i], " ", l[j], " -> ", GP(l[i], l[j]))


# A = 0
# B = 0
# # for i in range(DIMS+1):
# i = 3
# A = A + (e(i)*i)
# B = B + (e(i))
#
# #A = e(2)
# #B = e(1)
# print(A)
# print(-e(1))
# print(+e(1))
# print(LCONT(e(300)+10,e(300)))
# print(RCONT(e(300)+10,e(300)))
# print(GP(A,B))
# print(GP(e(1), e(2), T_OP));
#print(Multivector.get_N())
#print(GP(1+e(1), e(2)))
#for i in range(100):
#    print("LCONT: ", LCONT(3+e(1), 4+e(2)))

# for i in range(100):
# print("GP: ", GP(3+e(1), 4+e(2)))
# print("LCONT: ", LCONT(3+e(1), 4+e(2)))
# print("RCONT: ", RCONT(3+e(1), 4+e(2)))
# print("^: ", (3+e(1))^(4+e(2)))
# print("SCP: ", SCP((3+e(1)),(4+e(2))))
# print("\n\n\n\n\n\n")
