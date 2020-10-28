clc 
clear
close all 
path=input('������� ���� �� ����� � ���� \n ''C:\\����\\��\\�����\\file.mat''\n���� � ����\n ''file.mat'', ���� ���� ����� � �������� ����� \n')
load(path,'AA');%��������� ���������� �����
variant=input('\n������� ���� �������\n');
A=[AA(:,variant*2-1) AA(:,variant*2)];
n=50;
X=A(:,1);Y=A(:,2);
N=([1:1:50])';%������
tablet=table(N,X,Y);%�������
disp(tablet)
plot(X,Y,'.')
hold on
grid on
xlabel( 'X' );
ylabel( 'Y' );
Xs = sort(X);
Ys=sort(Y);
As=[Xs,Ys];
fprintf('�������������� ������� \n');
disp('Xs')
disp(Xs')
disp('Ys')
disp(Ys')
fprintf('\n\n')
% pointts=gca; %�������� ��������� figure
Xint=get(gca,'XTick');%��������� �� ��
Xintl=length(Xint);
Yint=get(gca,'YTick');%��������� �� �Y
Yintl=length(Yint);
Sqplotn=zeros(Xintl-1,Yintl-1);%������� ����������
for i=1:1:(Xintl-1)
    xtemp=find(Xint(i)<A(:,1)&A(:,1)<=Xint(i+1));
    for j=1:1:(Yintl-1)
        Sqplotn(i,j)=length(find(Yint(j)<A(xtemp,2)&A(xtemp,2)<=Yint(j+1)));
    end
    clear xtemp;
end
show=num2str([1 (Yint(2:end)) 1;[([Xint(2:end) 1])' [Sqplotn;sum(Sqplotn)] ([sum(Sqplotn') sum(sum(Sqplotn))])']]);
show(1,1:5)=' X\Y ';show(1,(end-1):end)='Ni';
show(end,1:4)=' Nj ';
disp(show)

%%%%%%%%%%%%%%%%%%% �� ������
nx=sum(X)/n;%�� ������� �
ny=sum(Y)/n;%�� ������� �
nsx=(1/(n-1))*(sum(X.^2)-n*nx^2);%��� �� ��������� �
nsy=(1/(n-1))*(sum(Y.^2)-n*ny^2);%��� �� ��������� �
nkxy=(1/(n-1))*(sum(X.*Y)-n*nx*ny);%����� ���� �� �� �������
npxy=nkxy/sqrt(nsy*nsx);%���.���� �� �� �������
%%%%%%%%%%%%%%%%%% �. ������
hX=(Xint(2)-Xint(1))/2; Xcenter=(Xint(1)+hX):2*hX:(Xint(end-1)+hX);%������ �������� �
hY=(Yint(2)-Yint(1))/2; Ycenter=(Yint(1)+hY):2*hY:(Yint(end-1)+hY);%������ �������� �
nid=sum(Sqplotn');%��������� ��� � 
ndj=sum(Sqplotn);%��������� ��� �
gx=sum(Xcenter.*nid)/n;%� ������� �
gy=sum(Ycenter.*ndj)/n;%� ������� �
gsx=(sum((Xcenter.^2).*nid)-n*(gx^2))/(n-1);%��� � ��������� �
gsy=(sum((Ycenter.^2).*ndj)-n*(gy^2))/(n-1);%��� � ��������� y
fSum=0;
for i=1:1:length(nid)
    for j=1:1:length(ndj)
       fSum=fSum+Xcenter(i)*Ycenter(j)*Sqplotn(i,j);
    end
end
gkxy=(fSum-n*gx*gy)/(n-1);%����� ���� �� � �������
gpxy=gkxy/sqrt(gsy*gsx);%���.���� �� �� �������
fprintf('\n\n')
matx=mean(X);maty=mean(Y);
matsx=var(X);matsy=var(Y);
mkxy=cov(A);matkxy=mkxy(2);
matpxy=corr(X,Y);
tablet3=table([nx;gx;matx],[ny;gy;maty],[nsx;gsx;matsx],[nsy;gsy;matsy],[nkxy;gkxy;matkxy],[npxy;gpxy;matpxy]);
tablet3.Properties.RowNames=({'ungroupped','groupped','matlab'});
tablet3.Properties.VariableNames={'xsr','ysr','sx','sy','kxy','pxy'};
disp(tablet3);
fprintf('\n\n')
disp(['�������� �������� H0: ' 961 ' = 0'])

disp(['������� ' 945 ':'])
a=input([945 ' = ']);
fprintf('\n\n �������� ����������� Z\n\n')
stud=tinv((1-a/2),n-2);
fprintf('t(%g)[%g]=%g\n',(1-a/2),(n-2),stud)

Zv=npxy*sqrt(n-2)/sqrt(1-npxy^2)

if(abs(Zv)>stud)
    disp('�������� ����������� � ������ H1. ���������� �������')
else
    disp('�������� �������')
end

fprintf('\n\n �������� ����������� U\n\n')
U=norminv(1-a/2);
fprintf('u(%g)=%g\n',(1-a/2),U)
Uv=sqrt(n-3)*atanh(npxy)


if(abs(Uv)>U)
    disp('�������� ����������� � ������ H1. ���������� �������')
else
    disp('�������� �������')
end

fprintf('\n\n������������ ������:\n')
left=tanh(atanh(npxy)-U/sqrt(n-3)-npxy/(2*(n-1)));
right=tanh(atanh(npxy)+U/sqrt(n-3)-npxy/(2*(n-1)));
fprintf(['\n%.4f<' 961 '[x,y]<%.4f\n'],left,right)
if(left>0||right<0)
    disp('�������� �� �������� 0, �.�. � ������������� ������������')
    disp('1-a ���������� ���������� ����� X � Y � ����� ����� ��������� ���������')
    syms x y x0;
    nyx=ny+npxy*sqrt(nsy/nsx)*(x-nx);
    fprintf('\nungroupped y(x)=%s\n with coeffs: %g\t%g',char(nyx),double(coeffs(nyx)))
    nxy=nx+npxy*sqrt(nsx/nsy)*(y-ny);
    fprintf('\nungroupped x(y)=%s\n with coeffs: %g\t%g',char(nxy),double(coeffs(nxy)))
    gyx=gy+gpxy*sqrt(gsy/gsx)*(x-gx);
    fprintf('\ngroupped y(x)=%s\n with coeffs: %g\t%g',char(gyx),double(coeffs(gyx)))
    gxy=gx+gpxy*sqrt(gsx/gsy)*(y-gy);
    fprintf('\ngroupped x(y)=%s\n with coeffs: %g\t%g\n',char(gxy),double(coeffs(gxy)))
    plot(X,subs(nyx,X),'k')
    plot(subs(nxy,Y),Y,'b'))
%     plot(X,subs(gyx,X),'r')
%     plot(X,subs(gxy,X),'g')
    legend('��������� �������������','y(x) �������','x(y) �������')%,'y(x) �����','x(y) �����')
 conyx=double(coeffs(nyx));anyx=conyx(2);bnyx=conyx(1);
 Qy=sum(Y.^2)-n*ny^2;
 Qr=(n-1)*(nkxy^2)/nsx;
 Qe=Qy-Qr;
 sOst=Qe/(n-2);
 R=Qr/Qy;
 fprintf('\na=%f,b=%f\n',anyx,bnyx)
 fprintf('\nQy=%f, Qr=%f, Qe=%f, ������ s^2=%f, \nR^2=%f\n',Qy,Qr,Qe,sOst,R)
 fprintf('\n �������� %f = %f\n',npxy,sign(anyx)*sqrt(R))
 fprintf('\n ����� �������� �������������\n')
 ti=tinv(1-a/2,n-2)
 xi1=chi2inv(1-a/2,n-2)
 xi2=chi2inv(a/2,n-2)
 fprintf('\n\n')
 aleft=anyx-ti*sqrt(sOst/((n-1)*nsx));
 aright=anyx+ti*sqrt(sOst/((n-1)*nsx));
 bleft=bnyx-ti*sqrt(sOst*sum(X.^2)/(n*(n-1)*nsx));
 bright=bnyx+ti*sqrt(sOst*sum(X.^2)/(n*(n-1)*nsx));
 dleft=(n-2)*sOst/(xi1);
 dright=(n-2)*sOst/(xi2);
 fprintf(['%.4f < a < %.4f \t %.4f < b < %.4f \t %.4f < ' 963 ' < %.4f\n'],aleft,aright,bleft,bright,dleft,dright)
 disp('������� ������������� ���������� ��� �������� �������� Y ��� � = �0')
 fprintf(['\ny0' 177 '%.4f' 8730 '%.4f*(1/%g + (x0-%.4f)^2/%.4f)\n'],ti,sOst,n,nx,(n-1)*nsx)
 disp('�������� ���������� �������� ��������� Y �� x:')
 if(aleft>0||aright<0)
     fprintf('�������� H0: a = 0 ����������� �� ������ ���������� a = %g ,\n��� ��� ������������� �������� �� ��������� ���� � ������������� ������������ %g\n',a,1-a)
     disp('��������� �������� ���� ���������, ��������� ���������� F')
     fv=(n-2)*Qr/Qe
     fprintf(['��� �����, fv ' 8713 '(%.4f;%.4f)\n'],finv(a/2,1,48),finv(1-a/2,1,48))
     disp('����� �������, �������� ��������� Y �� x ������������� �������')
 else
     fprintf('�������� H1: a = 0 ����������� �� ������ ���������� a = %g',a)
 end 
else
    disp('�������� �������� 0')
end

    


