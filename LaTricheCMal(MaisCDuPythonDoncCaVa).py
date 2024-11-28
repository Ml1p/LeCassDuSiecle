from sympy import *
import numpy as np

#-----------------------------------------
# EXO 1
#-----------------------------------------
x,y,z,t=symbols("x y z t")

A=[1,1,3]
B=[7,5,1]
V=[-2*y,4*z,-2*x]

C=[A[0]+t*(B[0]-A[0]),A[1]+t*(B[1]-A[1]),A[2]+t*(B[2]-A[2])]
dl=[diff(C[0],t),diff(C[1],t),diff(C[2],t)]

for i in range(3):
    V[i]=V[i].subs(x,C[0])
    V[i]=V[i].subs(y,C[1])
    V[i]=V[i].subs(z,C[2])

Vdl=V[0]*dl[0]+V[1]*dl[1]+V[2]*dl[2]

print("\nQuestion 1 : \n\nI = "+str(integrate(Vdl,(t,0,1))))



#-----------------------------------------
# EXO 2
#-----------------------------------------
rho,phi=symbols("rho phi")


A=[4,8*pi/4]
B=[4,5*pi/4]
V=[2*rho*(1+4*phi/pi),-2*rho*(1-4*phi/pi)]
R=4


V[1]=V[1]*rho
V[1]=V[1].subs(rho,R)
print("\nQuestion 2-1 :\n\nI = "+str(integrate(V[1],(phi,A[1],B[1]))/pi))

A=[7,8*pi/4]
B=[2,8*pi/4]
phi0=8*pi/4


#------------------------- Partie 2


V[0]=V[0].subs(phi,phi0)

print("\nQuestion 2-2 :\n\nI = "+str(integrate(V[0],(rho,A[0],B[0]))))

#-----------------------------------------
# EXO 3
#-----------------------------------------
u=symbols("u")


V=[u/(-3*x+8*y)**2,-16/(-3*x+8*y)**2]


g=V[1]-diff(integrate(V[0],x),y)
u0=solve(g,u)[0]
print("\nQuestion 3-1: \nu = "+str(solve(g,u)[0]))


A=[1,0]
B=[0,1]

f=[integrate(V[0],x),integrate(V[0],x)]
f[0]=f[0].subs(u,u0)
f[1]=f[1].subs(u,u0)
f[0]=f[0].subs(x,B[0])
f[0]=f[0].subs(y,B[1])
f[1]=f[1].subs(x,A[0])
f[1]=f[1].subs(y,A[1])

I=f[0]-f[1]

print("\nQuestion 3-2 : \n\nI = "+str(I))


#-----------------------------------------
# EXO 4
#-----------------------------------------


f=2*x**2+16*x*y-y*z
V=[2*x-y,-4*y+3*z,-z+3*x]


grad=[diff(f,x),diff(f,y),diff(f,z)]
eq1=Eq(grad[0],0)
eq2=Eq(grad[1],0)
eq3=Eq(grad[2],1)
print("\nQuestion 4-1:\n"+str(solve([eq1,eq2,eq3],(x,y,z))))

print("\nQuestion 4-2:\n"+str(diff(V[0],x)+diff(V[1],y)+diff(V[2],z)))

c=diff(V[1],x)-diff(V[0],y)
a=diff(V[2],y)-diff(V[1],z)
b=diff(V[0],z)-diff(V[2],x)

print("\nQuestion 4-3:\na = "+str(a)+"\tb = "+str(b)+"\tc = "+str(c))

#-----------------------------------------
# EXO 5
#-----------------------------------------

f=abs(x-y)
bx=[0,7]
by=[0,4]

print("\nQuestion 5:\n"+str(integrate(f,(x,bx[0],bx[1]),(y,by[0],by[1]))))

#-----------------------------------------
# EXO 6
#-----------------------------------------
a,b,c=symbols("a b c")


A=[2,5,3]
B=[6,5,3]
C=[2,10,3]
D=[2,5,7]
c0=20
f=y-5

eq=Eq(a*(x-D[0])+b*(y-D[1])+c*(z-D[2]),0)

eq=eq.subs(c,c0)
eq1=eq
eq2=eq

eq1=eq1.subs(x,C[0])
eq1=eq1.subs(y,C[1])
eq1=eq1.subs(z,C[2])
eq2=eq2.subs(x,B[0])
eq2=eq2.subs(y,B[1])
eq2=eq2.subs(z,B[2])

a0=solve([eq1,eq2],(a,b)).get(a)
b0=solve([eq1,eq2],(a,b)).get(b)

print("\nQuestion 6-1:\n"+str(solve([eq1,eq2],(a,b))))

t0=min(A[2],B[2],C[2],D[2])
u0=max(A[2],B[2],C[2],D[2])
v=min(A[1],B[1],C[1],D[1])
w=min(A[0],B[0],C[0],D[0])
print("\nQuestion 6-2:\nt = "+str(t0)+"\tu = "+str(u0)+"\tv = "+str(v)+"\tw = "+str(w))


fTrident=solve(a0*(x-D[0])+b0*(y-D[1])+c0*(z-D[2]),x)

print("Trident: "+str(fTrident))

fphi=solve(b0*(y-D[1])+c0*(z-D[2]),y)

a,b,c=symbols("a b c")

eq1=Eq(a*(C[2]-D[2])+b,C[1])
eq2=Eq(a*(D[2]-D[2])+b,D[1])
s=solve([eq1,eq2],(a,b))
aphi=s.get(a)
bphi=s.get(b)

fphi=aphi*(z-D[2])+bphi

print("Phi = "+str(fphi))

print("\nQuestion 6-3:\nI = "+str(integrate(f,(x,w,fTrident),(y,v,fphi),(z,t0,u0))))

#-----------------------------------------
# EXO 7
#-----------------------------------------
r,Teta=symbols("r Teta")


f=(x-2)**2+y**2-4
V=[2*r*cos(Teta),-2*r*sin(Teta)]
coefx=2
r0=9

b0=Rational(coefx*2,r0)

print("\nQuestion 7-1:\na = "+str(1/2)+"\tb = "+str(b0))

V[0]=V[0]*r**2*sin(Teta)
V[0]=V[0].subs(r,r0)

I=V[0].subs(sin(Teta),u)
I=I/cos(Teta)
print("Uniquement le commencement "+str(I))
I=integrate(I,(u,pi/2,sin(cos(phi))),(phi,-pi/2,pi/2)) #Franchement pas sur des bornes sur u (sin(cos(phi)) ça résoud pas le pb)
print("La fin ? "+str(I))


I=integrate(V[0],(Teta,0,b0*cos(phi)))
I=I.subs(cos(phi),u)
I=I/-sin(phi)
I=integrate(I,(u,-1,1),(phi,-pi/2,pi/2))


print("\nQuestion 7-2:\nI = "+str(I))