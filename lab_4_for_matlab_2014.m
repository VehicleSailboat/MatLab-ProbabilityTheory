clc 
clear
close all 
path=input('введите путь до файла в виде \n ''C:\\путь\\до\\файла\\file.mat''\nлибо в виде\n ''file.mat'', если файл лежит в основной папке \n')
load(path,'AA');%загружаем переменные рядов
variant=input('\nвведите свой вариант\n');
A=[AA(:,variant*2-1) AA(:,variant*2)];
n=50;
X=A(:,1);Y=A(:,2);
N=([1:1:50])';%номера
tablet=table(N,X,Y);%таблица
disp(tablet)
plot(X,Y,'.')
hold on
grid on
xlabel( 'X' );
ylabel( 'Y' );
Xs = sort(X);
Ys=sort(Y);
As=[Xs,Ys];
fprintf('Группированная выборка \n');
disp('Xs')
disp(Xs')
disp('Ys')
disp(Ys')
fprintf('\n\n')
% pointts=gca; %получаем свойситва figure
Xint=get(gca,'XTick');%интервалы по оХ
Xintl=length(Xint);
Yint=get(gca,'YTick');%интервалы по оY
Yintl=length(Yint);
Sqplotn=zeros(Xintl-1,Yintl-1);%таблица плотностей
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

%%%%%%%%%%%%%%%%%%% нг данные
nx=sum(X)/n;%нг среднее х
ny=sum(Y)/n;%нг среднее у
nsx=(1/(n-1))*(sum(X.^2)-n*nx^2);%нсм нг дисперсия х
nsy=(1/(n-1))*(sum(Y.^2)-n*ny^2);%нсм нг дисперсия у
nkxy=(1/(n-1))*(sum(X.*Y)-n*nx*ny);%коэфф корр по нг выборке
npxy=nkxy/sqrt(nsy*nsx);%выб.коэф по нг выборке
%%%%%%%%%%%%%%%%%% г. данные
hX=(Xint(2)-Xint(1))/2; Xcenter=(Xint(1)+hX):2*hX:(Xint(end-1)+hX);%центры отрезков х
hY=(Yint(2)-Yint(1))/2; Ycenter=(Yint(1)+hY):2*hY:(Yint(end-1)+hY);%центры отрезков у
nid=sum(Sqplotn');%плотность для х 
ndj=sum(Sqplotn);%плотность для у
gx=sum(Xcenter.*nid)/n;%г среднее х
gy=sum(Ycenter.*ndj)/n;%г среднее у
gsx=(sum((Xcenter.^2).*nid)-n*(gx^2))/(n-1);%нсм г дисперсия х
gsy=(sum((Ycenter.^2).*ndj)-n*(gy^2))/(n-1);%нсм г дисперсия y
fSum=0;
for i=1:1:length(nid)
    for j=1:1:length(ndj)
       fSum=fSum+Xcenter(i)*Ycenter(j)*Sqplotn(i,j);
    end
end
gkxy=(fSum-n*gx*gy)/(n-1);%коэфф корр по г выборке
gpxy=gkxy/sqrt(gsy*gsx);%выб.коэф по нг выборке
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
disp(['Проверка гипотезы H0: ' 961 ' = 0'])

disp(['введите ' 945 ':'])
a=input([945 ' = ']);
fprintf('\n\n Проверка статистикой Z\n\n')
stud=tinv((1-a/2),n-2);
fprintf('t(%g)[%g]=%g\n',(1-a/2),(n-2),stud)

Zv=npxy*sqrt(n-2)/sqrt(1-npxy^2)

if(abs(Zv)>stud)
    disp('гипотеза отклоняется в пользу H1. Корреляция значима')
else
    disp('гипотеза принята')
end

fprintf('\n\n Проверка статистикой U\n\n')
U=norminv(1-a/2);
fprintf('u(%g)=%g\n',(1-a/2),U)
Uv=sqrt(n-3)*atanh(npxy)


if(abs(Uv)>U)
    disp('гипотеза отклоняется в пользу H1. Корреляция значима')
else
    disp('гипотеза принята')
end

fprintf('\n\nИнтервальная оценка:\n')
left=tanh(atanh(npxy)-U/sqrt(n-3)-npxy/(2*(n-1)));
right=tanh(atanh(npxy)+U/sqrt(n-3)-npxy/(2*(n-1)));
fprintf(['\n%.4f<' 961 '[x,y]<%.4f\n'],left,right)
if(left>0||right<0)
    disp('Интервал не содержит 0, т.е. с доверительной вероятностью')
    disp('1-a существует корреляция между X и Y и имеет смысл уравнение регрессии')
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
    legend('диаграмма распределения','y(x) негрупп','x(y) негрупп')%,'y(x) групп','x(y) групп')
 conyx=double(coeffs(nyx));anyx=conyx(2);bnyx=conyx(1);
 Qy=sum(Y.^2)-n*ny^2;
 Qr=(n-1)*(nkxy^2)/nsx;
 Qe=Qy-Qr;
 sOst=Qe/(n-2);
 R=Qr/Qy;
 fprintf('\na=%f,b=%f\n',anyx,bnyx)
 fprintf('\nQy=%f, Qr=%f, Qe=%f, оценка s^2=%f, \nR^2=%f\n',Qy,Qr,Qe,sOst,R)
 fprintf('\n проверка %f = %f\n',npxy,sign(anyx)*sqrt(R))
 fprintf('\n найдём квантили распределения\n')
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
 disp('границы доверительных интервалов для среднего значения Y при х = х0')
 fprintf(['\ny0' 177 '%.4f' 8730 '%.4f*(1/%g + (x0-%.4f)^2/%.4f)\n'],ti,sOst,n,nx,(n-1)*nsx)
 disp('проверим значимость линейной регрессии Y на x:')
 if(aleft>0||aright<0)
     fprintf('Гипотеза H0: a = 0 отклоняется на уровне значимости a = %g ,\nтак как доверительный интервал не накрывает нуль с доверительной вероятностью %g\n',a,1-a)
     disp('попробуем получить этот результат, используя статистику F')
     fv=(n-2)*Qr/Qe
     fprintf(['как видим, fv ' 8713 '(%.4f;%.4f)\n'],finv(a/2,1,48),finv(1-a/2,1,48))
     disp('Таким образом, линейная регрессия Y на x статистически значима')
 else
     fprintf('Гипотеза H1: a = 0 принимается на уровне значимости a = %g',a)
 end 
else
    disp('интервал содержит 0')
end

    


