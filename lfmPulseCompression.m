%% 
%    线性调频信号的脉冲压缩
%    作者：Alex_Zikun
% 
%    介绍：对线性调频信号进行仿真，输出其时频域的相关信息，并模拟回波信号，
%    对其进行脉冲压缩和加窗处理。
% 
%    实验记录：
%     1.线性调频信号时域包络、相位；实部、虚部
%     2.线性调频信号频谱幅频、相频特性；实部、虚部
%     3.两个目标回波的时域和频域波形
%     4.信号通过匹配滤波器的输出结果（脉冲压缩）。
%     5.用Hamming窗抑制脉冲压缩结果副瓣
%%  基本参数 
    clc;clear all;close all;

    T = 10e-6; % LFM周期/脉宽 10us
    B = 60e6; % LFM带宽 60Mhz
    fs = 100e6; % 采样率 100MHz
    K = B/T;
%%  模拟发射信号
    n = round(15*T*fs);
    t = linspace(-10*T, 10*T,n);

    lfmT = rectpuls(t,T).*exp(1j*pi*K*t.^2);
    lfmF = fftshift(fft(fftshift(lfmT)));
    f = linspace(-fs,fs,n);

    %% 时域绘图
        figure();
        plot(diff(phase(lfmT)));
        title('LFM信号的时间-频率变化趋势图');
        xlabel('时间');
        ylabel('频率');
        xlim([7200,7800])
    % 包络
        figure();
        subplot(2,2,1);
        plot(t,abs(lfmT));
        title('LFM信号时域包络');
        xlabel('t/s');
        ylabel('幅度');
        xlim([-1e-5,1e-5])
        ylim([-0.5,1.5])
    % 相位
        subplot(2,2,2);
        plot(t,phase(lfmT));
        title('LFM信号时域相位');
        xlabel('t/s');
        ylabel('相位');
        xlim([-5e-6,5e-6])
    % 实部
        subplot(2,2,3);
        plot(t,real(lfmT));
        title('LFM信号时域实部');
        xlabel('t/s');
        ylabel('幅度');
        xlim([-1.5e-6,1.5e-6]);
        ylim([-1,1]);
    % 虚部
        subplot(2,2,4);
        plot(t,imag(lfmT));
        title('LFM信号时域虚部');
        xlabel('t/s');
        ylabel('幅度');
        xlim([-1.5e-6,1.5e-6]);
        ylim([-1,1]);
    %% 频域绘图
        figure();
        subplot(2,2,1);
        plot(f,abs(lfmF));
        title('LFM信号幅频特性');
        xlabel('Hz');
        ylabel('幅度');
        
        subplot(2,2,2);
        plot(unwrap(angle(lfmF)));
        title('LFM信号相频特性');
        xlabel('Hz');
        ylabel('相位');
        
        subplot(2,2,3);
        plot(f,real(lfmF));
        title('LFM信号频谱实部');
        xlabel('Hz');
        ylabel('幅度');    
        xlim([-3e7,3e7]);
        
        subplot(2,2,4);
        plot(f,imag(lfmF));
        title('LFM信号频谱虚部');
        xlabel('Hz');
        ylabel('幅度');    
        xlim([-3e7,3e7]);
        
        
    
%%  模拟回波信号
% 延迟为50us和51us的两个信号( 5000点和5100点,t从4000开始存储 )
    %延迟
    t1=50e-6;
    t2=51e-6;
    
    %由于实际信号时间从0开始,时间轴重新定义
    echo1=rectpuls((t-t1),T).*exp(1j*pi*K*(t-t1).^2);
    echo2=rectpuls((t-t2),T).*exp(1j*pi*K*(t-t2).^2);
    echo=echo1+echo2;

     figure();
        subplot(2,2,1)
        plot(t,abs(echo1),'r');
        hold on;
        plot(t,abs(echo2),'b');
        title('两个回波信号');
        xlabel('t/s');
        ylabel('幅度');
        xlim([3.5e-5,7e-5])
        ylim([0,1.5])
        
        subplot(2,2,2)
        plot(t,real(echo));
        title('回波信号时域实部');
        xlabel('t/s');
        ylabel('幅度');
        xlim([4.9e-5 5.2e-5])
        
        subplot(2,2,3)
        plot(f,fftshift(abs(fft(echo))));
        title('回波信号幅频特性');
        xlabel('Hz');
        ylabel('幅度');
        
        subplot(2,2,4)
        plot(t,imag(echo));
        title('回波信号时域虚部');
        xlabel('t/s');
        ylabel('幅度');
        xlim([4.9e-5 5.2e-5])
 
    lfmT = rectpuls(t,T).*exp(1j*pi*K*t.^2);
    lfmF = fftshift(fft(fftshift(lfmT)));
    f = linspace(-fs,fs,n);
    
%%  脉冲压缩(匹配滤波）
% 方法1：
% 把信号变到频域 （2048点的FFT)[-fs/2,fs/2]
% 频域相乘 H(w)  =   *S(w)
% 逆傅里叶变换
% 方法2：
% 驻留相位原理获得传递函数

% 这里采用方法1
    %回波信号的频域形式
    echo_F=fftshift(fft(echo));
    %回波信号的匹配滤波器
    Hf=fftshift(fft(conj(lfmT)));
    %脉压结果
    Pc_F=echo_F.*Hf;
    %脉压时域结果
    Pc_T=ifftshift(ifft(Pc_F));
    
%用Hamming窗抑制副瓣
    %制作汉明窗
    Hm = [zeros(1,2300) hamming(9900)' zeros(1,2800)];
    %频域加窗
    Pcw_F = Hm.*Pc_F;
    Pcw_T=ifftshift(ifft(Pcw_F));
    
    figure();
    subplot(2,2,1)
    plot(t,abs(Pc_T));
    title('脉冲压缩后的时域波形');
    xlim([4.5e-5 5.5e-5])
    
    subplot(2,2,2)
    plot(f,abs(Pc_F));
    title('脉冲压缩后的频谱');
    
    subplot(2,2,3)
    plot(f,abs(Pcw_T));
    title('脉压加窗后的时域波形');
     xlim([4.5e7 5.5e7])
    
    subplot(2,2,4)
    plot(f,abs(Pcw_F));
    title('脉压加窗后的频谱');
 