clear;clc
path=input('������� ���� �� ����� (� ��������)\n');
load(path,'AA');
variant=input('������� ��� �������\n��� �������: ');
ISH=(AA(:,variant))';
disp('��� ���:')
disp(sprintf('%.2f    ',ISH(1:16)))
disp(sprintf('%.2f    ',ISH(17:32)))
disp(sprintf('%.2f    ',ISH(33:48)))
disp(sprintf('%.2f    ',ISH(49:50)))
disp('����������� ���:')
A=sort(ISH);
disp(sprintf('%.2f    ',A(1:16)))
disp(sprintf('%.2f    ',A(17:32)))
disp(sprintf('%.2f    ',A(33:48)))
disp(sprintf('%.2f    ',A(49:50)))
m=7;
V=50;
disp(['����� ������(w) � ' 916])
w=range(A);
delta=w/m;
disp(sprintf('w = %.3f',w))
disp(sprintf([916 ' = %.3f'],delta))
[N,Z]=hist(A,m);
N=(N)';
Z=(Z)';
Notn=N./V;
Nhis=Notn./delta;
Nc=cumsum(N);
NcOtn=Nc./V;
Nj=N./delta;
number=(1:7)';
rz(1)=min(A);
for i=2:8
rz(i)=rz(i-1)+delta;
end

%%%%%%%%%%%

for i=1:(length(rz)-2)
    tabss(i,1:4)={'[' rz(i) rz(i+1) ')'};
end
for i=(length(rz)-1)
    tabss(i,1:4)={'[' rz(i) rz(i+1) ']'};
end
tablet=table(number,tabss,Z,N,Notn,Nhis,Nc,NcOtn,Nj);
disp(['    �                       ������                 ��������  �������  ���.   ������   ������ ���.���.   nj/' 916])
disp(['    �������                                         �������           ����   �������    ���   ����.      '])
disp(tablet)


%%%%%%%%%%%
figure(1);clf
hold on
grid on
bar(Z,Nj,1,'g')
plot(Z,Nj,'b',Z,Nj,'bx')
title('����������� ������ � �������')
ylabel(['nj/' 916])

%%%%%%%%%%%%%%
figure(2);clf;hold on
bar(Z,Nhis,1,'g')
grid on
ylabel(['nj/(n' 916 ')'])
title('����������� ������������� ������')


%%%%%%%%%%%%%%
figure(3);clf;hold on
stairs([(min(A)-delta) Z' max(A)+delta],[0 NcOtn' 1])
grid on
ylabel('F(x)')
%%%%%%%%%%%%%%
disp('  ')
disp('  ')
disp('������ �� ���������������� ������')
% disp('������ �������')
hx=(A(25)+A(26))/2;
% disp(sprintf('hx = %.4f',hx))
% disp('������ ��������')
xn=sum(A)/V;
% disp(sprintf('xn = %.4f',xn))
% disp('���������� ��������� ���������')
Dx=(sum(A.^2)-V*xn^2)/V;
% disp(sprintf('Dx = %.4f',Dx))
% disp('����������� ���������')
S2=(sum(A.^2)-V*xn^2)/(V-1);
% disp(sprintf('S2 = %.4f',S2))
disp('������     ������    ���������  �����������')
disp('�������    ��������  ���������  ���������')
table2=table(hx,xn,Dx,S2);
disp(table2)
disp('  ')
disp('  ')
disp('______________________________________')
disp('������ �� �������������� ������')
disp('______________________________________')
% disp('������ ��������')
xnG=sum(Z.*N)/V;
% disp(sprintf('xnG = %.4f',xnG))
% disp('���������� ��������� ���������')
DxG=(sum(Z.^2.*N)-V*xnG^2)/V;
% disp(sprintf('DxG = %.4f',DxG))
% disp('����������� ���������')
S2G=(DxG*V)/(V-1);
% disp(sprintf('S2G = %.4f',S2G))
% disp('���������� ����')

for i=1:(length(rz)-1)
plotn(i)=length(find(A<rz(i+1)&A>rz(i)));
end
plotn(1)=plotn(1)+1;
if(sum(plotn)~=50)
    plotn(7)=plotn(7)+1;
end
[nd,adnum]=max(plotn);
ad=rz(adnum);
if(adnum-1==0)
    ndmin=0;
else
    ndmin=plotn(adnum-1);
end
if(adnum+1==8)
    ndplus=0;
else
    ndplus=plotn(adnum+1);
end
dxG=ad+((nd-ndmin)/(2*nd-ndmin-ndplus))*delta;
% disp(sprintf('dxG = %.4f',dxG))

% disp('������ �������')

h=0;
aio=(A(25)+A(26))/2;

for i=1:7
    if(rz(i)<aio&&aio<rz(i+1))
        ah=rz(i);
        h=i;
        break
    end
end
nh=plotn(h);
nhm=sum(plotn(1:(h-1)));
hxG=ah+((V/2-nhm)/nh)*delta;
% disp(sprintf('hxG = %.4f',hxG))
disp('    ������ |������  |���������|����������|����������')
disp('    �������|��������|���������|��������� |����')
table3=table(hxG,xnG,DxG,S2G,dxG);
disp(table3)

text={'����������� ����� �� ���������(������� 3/4) � �������(������� 4/5).';'��� �������, ����� ��� �������������?'};
choice=menu(text,'�����������', '����������', '�������������');
switch choice 
    case 1
        ypl=unifpdf([(min(A)-delta) A max(A)+delta],min(A),max(A));
        y=unifcdf([(min(A)-delta) A max(A)+delta],min(A)-delta,max(A)+delta);
    case 2
        ypl=normpdf([(min(A)-delta) A max(A)+delta],mean(A),sqrt(var(A)));
        y=normcdf([(min(A)-delta) A max(A)+delta],mean(A),sqrt(var(A)));
    case 3
        ypl=exppdf([(min(A)-delta) A max(A)+delta],mean(A));
        y=expcdf([(min(A)-delta) A max(A)+delta],mean(A));
    otherwise
        errordlg('�� ���� �������? ��, ���� �� ���� �����','�� ��� ����-�� �� �����')
        ypl=normpdf([(min(A)-delta) A max(A)+delta],mean(A),sqrt(var(A)));
        y=normcdf([(min(A)-delta) A max(A)+delta],mean(A),sqrt(var(A)));
end
figure(2)
plot([(min(A)-delta) A max(A)+delta],ypl)
figure(3)
plot([(min(A)-delta) A max(A)+delta],y)
axis([(min(A)-delta) (max(A)+delta) -0.1 1.1])
choice=menu('������� ������ ����� �� �������?','���','��, �������');
switch choice
    case 1
        figure(1)
        set(gca,'XTick',unique(sort([rz Z'])))
        set(gca,'YTick',unique((sort(Nj))'))
        figure(2)
        set(gca,'XTick',unique(sort([rz Z'])))
        set(gca,'YTick',unique((sort([Nhis' max(ypl)]))'))
        figure(3)
        set(gca,'XTick',unique(sort([rz Z'])))
        set(gca,'YTick',unique((sort([0 NcOtn' 1]))'))
    case 2
        errordlg('�� ��','�� ��')
    otherwise
        errordlg('�� ���� �������? ��, ���� �� ���� �����','�� ��� ����-�� �� �����')
end

disp('_________________��� 3_________________')
disp('')
dovver=input('������� ������������� �����������\n');
a=1-dovver;
disp('')
disp('������������� ���������')
nlowintm=xn-sqrt(S2/V)*tinv((1-a/2),(V-1));
nupintm=xn+sqrt(S2/V)*tinv((1-a/2),(V-1));
glowintm=xn-sqrt(S2G/V)*tinv((1-a/2),(V-1));
gupintm=xn+sqrt(S2G/V)*tinv((1-a/2),(V-1));

nlowintp=(V-1)*S2/chi2inv((1-a/2),V-1);
nupintp=(V-1)*S2/chi2inv((a/2),V-1);
glowintp=(V-1)*S2G/chi2inv((1-a/2),V-1);
gupintp=(V-1)*S2G/chi2inv((a/2),V-1);

disp('��� ���������������� �������')
disp(sprintf(['%.3f<m<%.3f            %.3f<' 963 '^2<%.3f'],nlowintm,nupintm,nlowintp,nupintp))

disp('��� �������������� �������')
disp(sprintf(['%.3f<m<%.3f            %.3f<' 963 '^2<%.3f'],glowintm,gupintm,glowintp,gupintp))
disp('___________')
disp('�������� ���������� ���������� � ������� ���������� ������� normfit')
[matozh,sigma,m_int,sigma_int]=normfit(A,a);
disp(sprintf('m=%.4f ���� m=xn=%.4f \nsigma=%.4f ���� sigma=sqrt(S2)=%.4f \nm_int=(%.4f %.4f) \nsigma_int=(%.4f %.4f)\ndisp_int=(%.4f %.4f)',matozh,xn,sigma,sqrt(S2),m_int,sigma_int,sigma_int.^2))
disp('___________')
disp('')
disp('_________________�������� �������_________________')
disp('��������:')
disp('H01: mx = xn+0.5sqrt(S2)')
M0=xn+0.5*sqrt(S2);
disp(sprintf('M0 = %.3f',M0))
if(M0>nlowintm&&M0<nupintm)
    disp('�������� �������')
else
    disp('�������� �� �������')
end
disp('H02: Dx = 2*S2')
A0=2*S2;
disp(sprintf('A0 = %.3f',A0))
if(A0>nlowintp&&A0<nupintp)
    disp('�������� �������')
else
    disp('�������� �� �������')
end
disp('______________')
disp('H01: mx = xn+0.5sqrt(S2G)')
M0=xn+0.5*sqrt(S2G);
disp(sprintf('M0 = %.3f',M0))
if(M0>glowintm&&M0<gupintm)
    disp('�������� �������')
else
    disp('�������� �� �������')
end
disp('H02: Dx = 2*S2G')
A0=2*S2G;
disp(sprintf('A0 = %.3f',A0))
if(A0>glowintp&&A0<gupintp)
    disp('�������� �������')
else
    disp('�������� �� �������')
end

disp(['____________________���������� ���������� ' 967 '^2____________________'])
disp(' ')
disp(' ')
rzsz=rz;
for i=1:(length(rzsz)-1)
    chplotn(i)=length(find(A<rzsz(i+1)&A>rzsz(i)));
end
nums=find(chplotn<5);
for i=1:length(nums)
    if(nums(i)==1)
        nums(i)=nums(i)+1;
    end
end
rzsz(nums)=[];
clear chplotn;
for i=1:(length(rzsz)-1)
    chplotn(i)=length(find(A<rzsz(i+1)&A>rzsz(i)));
end
pnum=length(rzsz)-1;
rzsz(1)=-inf;rzsz(length(rzsz))=inf;
for i=1:pnum
    p(i)=normcdf(((rzsz(i+1)-xnG)/sqrt(S2G)),0,1)-normcdf(((rzsz(i)-xnG)/sqrt(S2G)),0,1);
end
for i=1:(length(rzsz)-1)
    chplotn(i)=length(find(A<rzsz(i+1)&A>rzsz(i)));
end
numer=(1:pnum)';
njshtrih=(p.*50)';
nj=chplotn';
njdif=(((nj-njshtrih).^2)./njshtrih);
for i=1:(length(rzsz)-2)
    tabss2(i,1:4)={'[' rzsz(i) rzsz(i+1) ')'};
end
for i=(length(rzsz)-1)
    tabss2(i,1:4)={'[' rzsz(i) rzsz(i+1) ']'};
end

table4=table(numer,tabss2,nj,njshtrih,njdif);
 disp('    �����  |                ��������               | nj  | njtheor   |(nj-njtheor)^2)/njtheor')
 disp(table4)
chism=sum(njdif);
disp(sprintf(['(' 967 '^2)� = %.4f'],chism))
kchi=length(number)-length(numer);
chichi=chi2inv(1-a,pnum-kchi-1);
disp(sprintf(['(' 967 '^2)[%g][%g] = %.4f'],(1-a),pnum-kchi-1,chichi))

if chism<chichi
    disp('�������� � ���������� ������������� �� �����������')
else
    disp('�������� � ���������� ������������� �����������')
end
