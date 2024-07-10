function[bRGN,wRGN,mRGN,sRGN,rRGN,cRGN]=RGN_script
n=1;%Number of runs
PopSize=50;
Iterations=1000;
Function_name='F1'
[Lb,Ub,Dim,Fun] = Get_CEC2005_Functions_details(Function_name)
for i=1:n
    [RGNbest,RGNfmin,bb]=KMN(Lb,Ub,Dim,Fun,Iterations,PopSize);
    RGNbest(i,:)=RGNbest;
    rRGN(i,:)=RGNfmin;
    eRGN(i,:)=bb;
end
disp('RGN runs completed');
cRGN=min(eRGN);
bRGN=min(rRGN);
wRGN=max(rRGN);
mRGN=mean(rRGN);
sRGN=std(rRGN);
RGNbest=min(RGNbest);