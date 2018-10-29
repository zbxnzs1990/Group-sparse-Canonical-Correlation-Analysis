%2013/08/13
%dlin5@tulane.edu
%Group sparse Canonical correlation analysis for genomic data integration
%CCA-sparse group algorithm
%
%
%%
% Function group sparse CCA
%      Covariance matrix K decomposition as Loss with the 
%           sparse group Lasso Regularization
%
%% Problem
%
%  min  || K - duv||^2 + lamda1_2 \|u\|_1 + lamda1_1 * sum_j w_j ||u_{G_j}||+lamda2_2 \|v\|_1 + lamda2_1 * sum_j w_j ||v_{G_j}||
%
%% Input parameters:
%
%  K-         Covariance Matrix of size m x n
%                n is the number of variables from data set 1
%                m is the number of variables from data set 2
%  u0-        initial values for loading vector u
%  v0-        initial value  for loading vector v
%  lamda1-    Penalization for group lasso
%  lamda2-    Penalization for lasso
%  opts-      Optional inputs (default value: opts=[]) with struct('groupP',P,'groupPK',Pk,'groupQ',Q,'groupQK',Qk)
%
%% Output parameters:
%
%  u,v-         Solutions

function [u,v]=SCCA_sGrpLR(K,u0,v0,lamda1,lamda2,opts)
Minerr=1e-2;
diffu=0.1*Minerr;
diffv=0.1*Minerr;
erru=1;
errv=1;
itersteps=1;
itermaxsteps=100;
%lamda1 for group penalization
%lamda2 for lasso penalization
lamda1u=lamda1(1);%vector of group penalty for u
lamda1v=lamda1(2);%vector of group penalty for v
lamda2u=lamda2(1);
lamda2v=lamda2(2);
%%% Group & Others 
P=opts.groupP;%the number of groups
Pk=opts.groupPK;%the number of variables in group K
Q=opts.groupQ;%the number of groups
Qk=opts.groupQK;%the number of variables in group K
if (length(Qk) ~=Q || length(Pk) ~=P)
    error('\n Check the number of groups!\n');
end

[p,q]=size(K);

unew=u0;
vnew=v0;

while(itersteps<=itermaxsteps && (erru>diffu || errv>diffv))
    uold=unew;
    vold=vnew;
   %%update singular vector u
   %%compute K*v
    Kv=K*vnew;
    Kvlength=norm(Kv);
    if(Kvlength==0)
        Kvlength=1;
    end
    Kv=Kv/Kvlength;
   %update u in each group
   gInd0=1;%the beginning index of group    
   for i=1:P
        Ind_Pi=[gInd0:gInd0+Pk(i)-1];% the index of i-th group
        Kvh=Kv(Ind_Pi);
       SKvh=SoftThreshold(2*Kvh,lamda2u);
       J=norm(SKvh);
       if(J<=lamda1u(i))
           unew(Ind_Pi)=0;
       else
           unew(Ind_Pi)=0.5*(SKvh-lamda1u(i)*(SKvh/J));
       end
       gInd0=gInd0+Pk(i);%update the beginning index of group 
   end
    %normalize unew
    ulength=norm(unew);
    if(ulength==0)
        ulength=1;
    end
    unew=unew/ulength;
    
    %%update v
    %%compute Kt*u
    KTu=K'*unew;
    KTulength=norm(KTu);
    if(KTulength==0)
        KTulength=1;
    end
    KTu=KTu/KTulength;
    
   %update v in each group
   gInd0=1;%the beginning index of group    
   for j=1:Q
        Ind_Qj=[gInd0:gInd0+Qk(j)-1];% the index of j-th group
        KTuh=KTu(Ind_Qj);
       SKTuh=SoftThreshold(2*KTuh,lamda2v);
       J=norm(SKTuh);
       if(J<=lamda1v(i))
           vnew(Ind_Qj)=0;
       else
           vnew(Ind_Qj)=0.5*(SKTuh-lamda1v(i)*(SKTuh/J));
       end
       gInd0=gInd0+Qk(j);%update the beginning index of group 
   end
    %normalize vnew
    vlength=norm(vnew);
    if(vlength==0)
        vlength=1;
    end
    vnew=vnew/vlength;
    
    erru=max(abs(uold-unew));
    errv=max(abs(vold-vnew));
    itersteps=itersteps+1;
end
u=unew;
v=vnew;
