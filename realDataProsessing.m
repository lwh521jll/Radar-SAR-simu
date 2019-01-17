%% 
%    二维分离SAR成像——实测数据处理
%    作者：Alex_Zikun
% 
%    介绍：实测数据为二维SAR回波数据，需要分别对方位向和距离向上进行脉冲压缩，得出的图像经过
%    简单的数字图像处理，输出

%%  基本参数和配置 
    clc;clear all;close all;
    load('raw_data_4k4k.mat');
    v_c =3e+8;%光速
    Tp =44.172e-6;%发射脉冲时间
    Br=59.44e6;%距离向带宽
    lamda=0.055;%波长
    f0=v_c/lamda;%载频
    vx=7200;%雷达平台运动速度
    R0=820e3;%场景中心最短斜距
    Kr=Br/Tp;%调频斜率
    Nr=4096;%距离向采样点数,必须要大于Tp*Br
    Fr=66.72839509333333e6;%距离向采样频率
    deta_t=1/Fr;%距离向采样时间间隔
    tr=2*R0/v_c+((0:Nr-1)-Nr/2)*deta_t;%距离向采样时间轴
    fr=((-Nr/2):(Nr/2-1))/Nr*Fr;
    PRF=1925;%PRF
    Na=4096;%方位向采样点数
    ta=((0:(Na-1))-Na/2)/PRF;
    T_sar=Na/PRF;%方位向采样时间
    fa=((-Na/2):(Na/2-1))/Na*PRF;
    theata_syn=2*atan(T_sar*vx/2/R0);%
    Ka=2*vx^2/(lamda*R0);
%%
    Echo=circshift(fft2(circshift(raw_data,[-Nr/2,-Na/2])),[Nr/2,Na/2]);% 变换到频域
    fr_matrix=repmat(fr.',1,Na); % 距离向脉冲压缩
    fa_matrix=repmat(fa,Nr,1); % 方位向脉冲压缩
    ref=exp(1j*pi*fr_matrix.^2/Kr);
    COMP=ref.*Echo; % 进行完距离向脉冲压缩

    Haz=exp(-1j*pi.*(fa_matrix.^2)./Ka);
    SAC=COMP.*Haz; % 进行完方位向脉冲压缩
    sac=circshift(ifft2(circshift(SAC,[-Nr/2,-Na/2])),[Nr/2,Na/2]);
    picture=abs(sac)/max(max(abs(sac)));
%% 图像处理
%     imshow(picture);title('处理前')
    picture2=histeq(abs(picture)); % 直方图均衡化
    picture2= filter2(fspecial('average',3),picture2)-0.3; % 平均值滤波降噪处理
    imshow(picture2);title('处理后')
    imwrite(picture2,'SAR_image.jpg');