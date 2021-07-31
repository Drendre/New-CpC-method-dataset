%% Test against exact solution in 1D
%% loop for the timesteps 
clear; fprintf('\n')
fact=1.08;
tic
Mh=18; %How many different timestep we want 
D=1; xf=1; TIME=10; %Diffusion coefficient, lenght of the space and time interval, x parameter
dp=0.5; p0=0.2; Np=6; pf=p0+dp*Np; %Np=ceil((pf-p0)/dp+0.1); %How many different parameters we want
N=101; c1=0; c2=0; lamb=3; k0=2; % initial number of cells, coefficients of sin modes, wavelength
MaxD=zeros(Mh,Np+1);Axhstep=zeros(Mh,1); rfun=zeros(Mh,1);%%Errors and axis for the timestep-sizes
h=0.0012; Dx=xf/(N-1); %%initial timestep and Space-step
axis=zeros(N,1); for i=1:N axis(i)=(i-1)*Dx; end  %x axis
M1=zeros(N,N); 
u0=zeros(N,1); uex=zeros(N,1); %%initial values and exact solution
uee=zeros(N,1);utr=zeros(N,1);umi=zeros(N,1);ura=zeros(N,1);udf=zeros(N,1);u2=zeros(N,1);uk1=zeros(N,1); du=zeros(N,1);
E=D/(Dx)^2; %%matrix elements
k=zeros(N,1); %%sources, now zero
for i=1:N 
    k(i)=k0*D*sin(pi*axis(i));
u0(i)=c1*sin(pi*axis(i))+c2*sin(lamb*pi*axis(i)); 
     uex(i)=(c1*exp(-pi^2*D*TIME)+k0/pi^2*(1-exp(-pi^2*D*TIME)))*sin(pi*axis(i)) +c2*sin(lamb*pi*axis(i))*exp(-lamb^2*pi^2*D*TIME);  %%Initial and exact final results
end
M1(1,1)=-E;M1(1,2)=E; M1(N,N)=-E;M1(N,N-1)=E;
for i=2:N-1     M1(i,i-1)=E; M1(i,i)=-2*E; M1(i,i+1)=E; end
Neb1=zeros(N,1);Neb2=zeros(N,1); 
for i=1:N
    if (i>1) Neb1(i)=M1(i,i-1); end  %bal szomszéd
    if (i<N) Neb2(i)= M1(i,i+1); end %jobb szomszéd
   end
%%%%%%%%%%%%
EG = eig(M1);Emax=0; E0=-1; Emin=-2; %%Eigenvalues
for i=1:N if EG(i)<Emax Emax=EG(i); end
    if EG(i)>E0 
       E0=EG(i); end 
end
for i=1:N if EG(i)>Emin
      if EG(i)<E0 
           Emin=EG(i); end 
    end
end
hMAX= -2/Emax  %% max timestep for Expl Euler
StiffRatio=Emax/Emin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ih=1:Mh  %% Big loop for the timestep, ih is the Serial Number of run, 1 for the largest timestep; 
T = uint32(TIME/h)+1; 
if(h>0.0000001 && h<0.001) rfun(ih)=10000000*h*h; end
UC=zeros(N,1);Ux=zeros(N,1); UTEMP=zeros(N,1); UTEMP2=zeros(N,1); a0=zeros(N,1); %%Final U and temporary results
Axhstep(ih)=h; %%Timestep axis
for i=1:N  b(i)=-M1(i,i); ee(i)=exp(-b(i)*h); end
UC(:)=u0(:); %%%% 1 Stage CN
tic
for t=1:T-1
    for i=2:N-1  %%%Neighbours are constans
        a0(i)=k(i);
     if(i>1) a0(i)=a0(i)+Neb1(i)*UC(i-1); end  %%left neighbour
     if(i<N) a0(i)=a0(i)+Neb2(i)*UC(i+1); end     %%rigth neighbour
       UTEMP(i)=UC(i)*ee(i)+ a0(i)/b(i)*(1-ee(i));
    end
    UTEMP(1)=0;UTEMP(N)=0;
 UC(:)=UTEMP(:); 
end
toc
%%Constant errors
Kul=zeros(N,1); MaxKul=0; 
for i=1:N
     Kul(i)=UC(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,1)=MaxKul; 
fprintf('Const h: %f \nMaxK: %f \n',h,MaxKul)
%%%%%%%%%%% 2nd order Constant-x starts here
for ip=0:Np-1  %% Big loop for the parameter
p=p0+ip*dp;
    h1=p*h; 
Ux(:)=u0(:);  
for i=1:N  
ee1(i)=exp(-b(i)*h1);   
end
%% Loop for time
tic
for t=1:T-1
  for i=2:N-1  %%%Neighbours are constans, first stage
         a0(i)=k(i);
     if(i>1) a0(i)=a0(i)+Neb1(i)*Ux(i-1); end  %%left neighbour
     if(i<N) a0(i)=a0(i)+Neb2(i)*Ux(i+1); end     %%rigth neighbour  
       UTEMP(i)=Ux(i)*ee1(i)+ a0(i)/b(i)*(1-ee1(i));
  end
     UTEMP(1)=0;UTEMP(N)=0; %%Boundary conditions
UTEMP(:)=(1-1/2/p)*Ux(:)+1/2/p*UTEMP(:);
 for i=2:N-1  %%%Neighbours are constans, second stage
     a0(i)=k(i);
     if(i>1) a0(i)=a0(i)+Neb1(i)*UTEMP(i-1); end  %%left neighbour
     if(i<N) a0(i)=a0(i)+Neb2(i)*UTEMP(i+1); end     %%rigth neighbour    
       UTEMP2(i)=Ux(i)*ee(i)+ a0(i)/b(i)*(1-ee(i));
  end
  Ux(:)=UTEMP2(:);   Ux(1)=0;Ux(N)=0;
end
toc
%% Const-x errors
MaxKul=0; 
for i=1:N
    Kul(i)=Ux(i)-uex(i);
    if MaxKul<abs(Kul(i)) MaxKul=abs(Kul(i)); end
end
MaxD(ih,ip+2)=MaxKul; 
end
%fprintf('Const-x h: %f \nMaxK: %f \n',h,MaxKul)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  remh=rem(ih,3);
%      if (remh==2) h=2*h/5;
%      else h=h/2; end
h=h/fact;
end
toc
% plot(axis,u0,axis,uex,axis,UC,'--o','LineWidth',2)
%plot(Axhstep(:),MaxD(:,1),Axhstep(:),MaxD(:,2),Axhstep(:),MaxD(:,3),Axhstep(:),MaxD(:,4),Axhstep(:),MaxD(:,5),Axhstep(:),MaxD(:,6),Axhstep(:),MaxD(:,7),'--o','LineWidth',2)
plot(Axhstep(:),MaxD(:,1),Axhstep(:),MaxD(:,2),Axhstep(:),MaxD(:,3),Axhstep(:),MaxD(:,4),Axhstep(:),MaxD(:,5),Axhstep(:),MaxD(:,7),Axhstep(:),rfun(:),'--','LineWidth',1)
set(gca,'xScale','log')
set(gca,'yScale','log')
grid on
