library(MASS)
library(pracma)

p=3
t=3
n=18

T_d1=matrix(c(1,0,0,0,0,1,0,1,0),nrow=p,ncol=t,byrow=T)
T_d2=matrix(c(1,0,0,0,0,1,0,1,0),nrow=p,ncol=t,byrow=T)
T_d3=matrix(c(1,0,0,0,0,1,0,1,0),nrow=p,ncol=t,byrow=T)
T_d4=matrix(c(1,0,0,0,0,1,0,1,0),nrow=p,ncol=t,byrow=T)
T_d5=matrix(c(1,0,0,0,0,1,0,1,0),nrow=p,ncol=t,byrow=T)
T_d6=matrix(c(1,0,0,0,0,1,0,1,0),nrow=p,ncol=t,byrow=T)
T_d7=matrix(c(0,1,0,1,0,0,0,0,1),nrow=p,ncol=t,byrow=T)
T_d8=matrix(c(0,1,0,1,0,0,0,0,1),nrow=p,ncol=t,byrow=T)
T_d9=matrix(c(0,1,0,1,0,0,0,0,1),nrow=p,ncol=t,byrow=T)
T_d10=matrix(c(0,1,0,1,0,0,0,0,1),nrow=p,ncol=t,byrow=T)
T_d11=matrix(c(0,1,0,1,0,0,0,0,1),nrow=p,ncol=t,byrow=T)
T_d12=matrix(c(0,1,0,1,0,0,0,0,1),nrow=p,ncol=t,byrow=T)
T_d13=matrix(c(0,0,1,0,1,0,1,0,0),nrow=p,ncol=t,byrow=T)
T_d14=matrix(c(0,0,1,0,1,0,1,0,0),nrow=p,ncol=t,byrow=T)
T_d15=matrix(c(0,0,1,0,1,0,1,0,0),nrow=p,ncol=t,byrow=T)
T_d16=matrix(c(0,0,1,0,1,0,1,0,0),nrow=p,ncol=t,byrow=T)
T_d17=matrix(c(0,0,1,0,1,0,1,0,0),nrow=p,ncol=t,byrow=T)
T_d18=matrix(c(0,0,1,0,1,0,1,0,0),nrow=p,ncol=t,byrow=T)

T_d=rbind(T_d1, T_d2, T_d3, T_d4, T_d5, T_d6, T_d7, T_d8, T_d9, T_d10, T_d11, T_d12, T_d13, T_d14, T_d15, T_d16, T_d17, T_d18)

T_d_ortho_1=matrix(c(1,0,0,0,1,0,0,0,1),nrow=p,ncol=t,byrow=T)

psi=matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0),nrow=p,ncol=p,byrow=T)

I_p=diag(p)
I_n=diag(n)
I_t=diag(t)

J_p=matrix(rep(1,p^2),nrow=p,ncol=p,byrow=T)
J_n=matrix(rep(1,n^2),nrow=n,ncol=n,byrow=T)
J_t=matrix(rep(1,t^2),nrow=t,ncol=t,byrow=T)

H_n=I_n-J_n/n
H_t=I_t-J_t/t

F_d=kronecker(I_n,psi)%*%T_d

r1=seq(-0.49,0.99,0.01)

Nr1=c()
Dr1=c()
A_efficiency=c()

Nr2=c()
Dr2=c()
D_efficiency=c()

Nr3=c()
Dr3=c()
E_efficiency=c()


for(i in 1:length(r1))
{
  V_1=matrix(c(1,r1[i],r1[i]^2,r1[i],1,r1[i],r1[i]^2,r1[i],1),nrow=p,ncol=p,byrow=T)
  V_2=matrix(c(1,r1[i],r1[i],r1[i],1,r1[i],r1[i],r1[i],1),nrow=p,ncol=p,byrow=T)
  
  x=2
  
  delta_1=1/sum(rowSums(solve(V_1)))
  V_1_star=solve(V_1)-delta_1*solve(V_1)%*%J_p%*%solve(V_1)
  A_1_star=kronecker(H_n,V_1_star)
  
  
  delta_2=1/sum(rowSums(solve(V_2)))
  V_2_star=solve(V_2)-delta_2*solve(V_2)%*%J_p%*%solve(V_2)
  A_2_star=kronecker(H_n,V_2_star)
  
  C_d11_1=t(T_d)%*%A_1_star%*%T_d
  C_d12_1=t(T_d)%*%A_1_star%*%F_d
  C_d21_1=t(C_d12_1)
  C_d22_1=t(F_d)%*%A_1_star%*%F_d
  
  C_d_1 = C_d11_1 - C_d12_1%*%pinv(C_d22_1)%*%C_d21_1
  
  C_d11_2=t(T_d)%*%A_2_star%*%T_d
  C_d12_2=t(T_d)%*%A_2_star%*%F_d
  C_d21_2=t(C_d12_2)
  C_d22_2=t(F_d)%*%A_2_star%*%F_d
  
  C_d_2 = C_d11_2 - C_d12_2%*%pinv(C_d22_2)%*%C_d21_2
  
  s1=sum(1/Re(eigen(C_d_1)$values[1:t-1]))
  s2=sum(1/Re(eigen(C_d_2)$values[1:t-1]))
  
  Dr1[i]=x*s1+s2
  
  p1=prod(1/Re(eigen(C_d_1)$values[1:t-1]))
  p2=prod(1/Re(eigen(C_d_2)$values[1:t-1]))
  
  Dr2[i]=x^(t-1)*p1*p2
  
  u1=max(1/Re(eigen(C_d_1)$values[1:t-1]))
  u2=max(1/Re(eigen(C_d_2)$values[1:t-1]))
  
  Dr3[i]=max(x*u1,u2)
    
  q_1_11 = sum(diag(t(T_d_ortho_1)%*%V_1_star%*%T_d_ortho_1))
  q_1_12 = sum(diag(t(T_d_ortho_1)%*%V_1_star%*%psi%*%T_d_ortho_1))
  q_1_22 = sum(diag(t(T_d_ortho_1)%*%t(psi)%*%V_1_star%*%psi%*%T_d_ortho_1)) - V_1_star[1,1]/t
  
  Q_1 = n/(t-1)*matrix(c(q_1_11,q_1_12,q_1_12,q_1_22),nrow=2,ncol=2,byrow=T)
  
  C_d_ortho_1 = (det(Q_1)/q_1_22)*H_t
  
  q_2_11 = sum(diag(t(T_d_ortho_1)%*%V_2_star%*%T_d_ortho_1))
  q_2_12 = sum(diag(t(T_d_ortho_1)%*%V_2_star%*%psi%*%T_d_ortho_1))
  q_2_22 = sum(diag(t(T_d_ortho_1)%*%t(psi)%*%V_2_star%*%psi%*%T_d_ortho_1)) - V_2_star[1,1]/t
  
  Q_2 = n/(t-1)*matrix(c(q_2_11,q_2_12,q_2_12,q_2_22),nrow=2,ncol=2,byrow=T)
  
  C_d_ortho_2 = (det(Q_2)/q_2_22)*H_t
  
  s1_ortho=sum(1/Re(eigen(C_d_ortho_1)$values[1:t-1])) 
  s2_ortho=sum(1/Re(eigen(C_d_ortho_2)$values[1:t-1]))
  
  Nr1[i]=x*s1_ortho+s2_ortho
  
  p1_ortho=prod(1/Re(eigen(C_d_ortho_1)$values[1:t-1]))
  p2_ortho=prod(1/Re(eigen(C_d_ortho_2)$values[1:t-1]))
  
  Nr2[i]=x^(t-1)*p1_ortho*p2_ortho
  
  u1_ortho=max(1/Re(eigen(C_d_ortho_1)$values[1:t-1]))
  u2_ortho=max(1/Re(eigen(C_d_ortho_2)$values[1:t-1]))
  
  Nr3[i]=max(x*u1_ortho,u2_ortho)
  
  A_efficiency[i]=Nr1[i]/Dr1[i]
  
  D_efficiency[i]=Nr2[i]/Dr2[i]
  
  E_efficiency[i]=Nr3[i]/Dr3[i]
  
}

data=data.frame(A_efficiency,D_efficiency,E_efficiency)
names(data)=c("A Efficiency", "D Efficiency", "E Efficiency")

#Calling library writexl
library(writexl)
#Storing dataframe in given path as excel file
write_xlsx(data,"D:/Code/Example.xlsx")

#Calling library readxl
library(readxl)

Example=read_excel("D:/Code/Example.xlsx")

plot(r1,Example$`A Efficiency`,xlab=expression('r'[(1)]), ylab=expression(A-Efficiency['d'[0]]),type="l",lty="solid",lwd=3)

plot(r1,Example$`D Efficiency`,xlab=expression('r'[(1)]), ylab=expression(D-Efficiency['d'[0]]),type="l",lty="solid",lwd=3)

plot(r1,Example$`E Efficiency`,xlab=expression('r'[(1)]), ylab=expression(E-Efficiency['d'[0]]),type="l",lty="solid",lwd=3)

