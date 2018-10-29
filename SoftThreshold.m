function S=SoftThreshold(x,lamda)
for i=1:length(x)
if(abs(x(i))>lamda)
    S(i,1)=sign(x(i))*(abs(x(i))-lamda);
else
    S(i,1)=0;
end
end