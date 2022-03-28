clc;
clear;
close all;
tic;

fs=4000;   % 采样频率为4000Hz
t=(0:1/fs:(2-1/fs))';
N=length(t);
load('x5');
load('nt2');
x=x5;
y=x+nt2;
% plot(y);title('original signal x(t)')

%% ICEEMDAN分解
Nstd = 0.2;                % 正负高斯白噪声标准表
NR = 100;                  % 加入噪声的次数
MaxIter = 500;   
SNRFlag = 2;
% 最大迭代次数
%[imf, its]=ceemdan(y,Nstd,NR,MaxIter);
[imf, its]=iceemdan(y,Nstd,NR,MaxIter,SNRFlag);
%[imf, its]=emd(y);
%[imf, its]=eemd(y,Nstd,NR,MaxIter);
[m, n]=size(imf);
CC=zeros(1,12);  % 相关系数
figure;
for i=1:12
    subplot(6,2,i);plot(imf(i,:));ylabel(['IMF',num2str(i)]);
    CC(i)=corr(imf(i,:)',y,'type','Pearson');   % 相关系数
end

%% 小波包降噪：对含伪迹分量进行小波包分解
Y_imf=zeros(N,9);
v=[1,2,3,4,8,9];
%v=[1,2,3,4,5,9];
%v=[1,2,3,8,9];
for i=v
    y_imf=imf(i,:)';
    wpt=wpdec(y_imf,3,'db7');  % 小波包分解，3代表分解3层，'db8'使用Daubechies小波
    nodes=get(wpt,'tn');  % 小波包分解系数
    N_cfs=length(nodes);  %小波包系数个数
    ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
    nodes_ord=nodes(ord); %重排后的小波系数
    rex3=zeros(N,N_cfs);
    C1=zeros(1,N_cfs);
    
    
    
    ysoft=zeros(N,N_cfs);
    for j=1:N_cfs
        rex3(:,j)=wprcoef(wpt,nodes_ord(j));  % 得到第3层各个节点对应的信号
        thr=thselect(rex3(:,j),'rigrsure');              % 阈值获取
        ysoft(:,j)=wthresh(rex3(:,j),'s',thr);   % 进行软阈值处理
    end
    
    % 每一个imf分量的小波包分解系数重构
    Y_imf(:,i)=sum(ysoft,2);
end

%% 将去噪完成的模态分量与保留的未经处理的模态分量第一次重构
imf=imf';
l=[5,6,7];
%l=[6,7,8];
%l=[4,5,6,7];
Y=sum(Y_imf,2)+sum(imf(:,l),2);
%REMD(Y,x,fs);
CC1=corr(Y,x,'type','Pearson');   % 相关系数
disp(CC1);

%% 画图
% 时域波形对比图
figure;
subplot(3,1,1);plot(t,x);xlabel('t/s');ylabel('Amplitude/mV');title('original signal x(t)');
subplot(3,1,2);plot(t,y);xlabel('t/s');ylabel('Amplitude/mV');title('noise-containing signal y(t)');
subplot(3,1,3);plot(t,Y);xlabel('t/s');ylabel('Amplitude/mV');title('denoised signal Y(t)');

% 频谱对比图
figure;
[f,A] = PinPu(x,fs);
subplot(3,1,1);plot(f,A);xlabel('frequency/Hz');ylabel('Amplitude/mV');title('original signal x(t)');
[f,A] = PinPu(y,fs);
subplot(3,1,2);plot(f,A);xlabel('frequency/Hz');ylabel('Amplitude/mV');title('noise-containing signal y(t)');
[f,A] = PinPu(Y,fs);
subplot(3,1,3);plot(f,A);xlabel('frequency/Hz');ylabel('Amplitude/mV');title('denoised signal Y(t)');

%% 降噪指标
SNR1(x,y);
SNR1(x,Y);
% 降噪前信噪比
p1=sum(abs(y).^2)/N;
p2=sum(abs(nt2).^2)/N;
SNR(1)=10*log10(p1/p2);
% 降噪后信噪比
p3=sum(abs(Y).^2)/N;
p4=sum(abs(Y-x).^2)/N;
SNR(2)=10*log10(p3/p4);
%均方根误差
RMSE=sqrt((mean((Y-x).^2))/N);

toc;