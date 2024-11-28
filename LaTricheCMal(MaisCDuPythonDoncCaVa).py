from sympy import *
import numpy as np


# Modifiez que les trucs avant la logique (Les Champs de vecteurs, les points, les fonctions...)



#-----------------------------------------
# EXO 1
#-----------------------------------------
x,y,z,t=symbols("x y z t")


A=[5,2,3]   # Point A (séparez les coordonées par des virgules comme dans l'exemple)
B=[2,1,5]   # Point B
V=[4*y,-2*z,-x]     # Champs de Vecteur V (Coordonnées séparées par des virgules)
# [Vx , Vy ]

# Pas touche après ça
#-----------------------------------------


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


V=[-5*rho*(1-3*phi/pi),rho*(1+2*phi/pi)]

R=4            # Valeur de R Donée dans la première partie de l'exo
A=[R,2*pi/4]
B=[R,pi/4]

# Partie 2

phi0=2*pi/4 # Le phi qu'ils vous donnent dans la seconde partie de l'exo
A2=[6,phi0]
B2=[1,phi0]


# Pshhit
#-----------------------------------------

V[1]=V[1]*rho
V[1]=V[1].subs(rho,R)
print("\nQuestion 2-1 :\n\nI = "+str(integrate(V[1],(phi,A[1],B[1]))/pi))

V[0]=V[0].subs(phi,phi0)

print("\nQuestion 2-2 :\n\nI = "+str(integrate(V[0],(rho,A2[0],B2[0]))))

#-----------------------------------------
# EXO 3
#-----------------------------------------
u=symbols("u")


V=[u/(-7*x-8*y)**2,16/(-7*x-8*y)**2]  # Champs de vecteur (notez u, u simplement)
                                    # Et ² c'est **2


# Retirez vos mains après cette partie
#-----------------------------------------
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

f=-2*x**2-16*x*y+2*y*z         # Fonction qu'ils vous donnent
V=[-2*x+6*y,2*y-5*z,2*z+4*x]


# Toujours pareil, après ça c'est le cerveau
#-----------------------------------------

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

bx=[0,7] # Bornes de l'intégrale sur x, [borne négative, borne positive]
by=[0,4] # Bornes de l'intégrale sur y


# Approchez à vos risques et périls
#-----------------------------------------

f=abs(x-y)
print("\nQuestion 5:\n"+str(integrate(f,(x,bx[0],bx[1]),(y,by[0],by[1]))))


#-----------------------------------------
# EXO 6
#-----------------------------------------
a,b,c=symbols("a b c")

A=[5,3,6]
B=[6,3,6]
C=[5,7,6]
D=[5,3,9]
c0=4        # La valeur de c qu'ils vous donnent dans l'énoncé
f=y-3       # La fonction qu'il faut intégrer


# Bas les pattes
#-----------------------------------------

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

print("\nTrident (ou psi comme vous préférez): "+str(fTrident)+" (C'est pas factorisé comme ils le demandent mais vous pouvez le copier coller et ça passe crême)")

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

V=[-5*r*cos(Teta),4*r*sin(Teta)]

coefx=5     # C'est le 2 dans (x-2)² ou le 5 dans (x-5)² (et il est positif c'est pas -5)
r0=11       # Le r constant qu'ils vous donnent


# Fuyez pauvres fou
#-----------------------------------------

b0=Rational(coefx*2,r0)

print("\nQuestion 7-1:\na = "+str(1/2)+"\tb = "+str(b0))

V[0]=V[0]*r**2*sin(Teta)
V[0]=V[0].subs(r,r0)

I=V[0].subs(sin(Teta),u)
I=I/cos(Teta)
I=integrate(I,(u,0,b0*cos(phi)),(phi,-pi/2,pi/2))/pi

print("\nQuestion 7-2:\nI = "+str(I))