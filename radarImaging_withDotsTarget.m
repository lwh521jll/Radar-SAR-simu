%% 
%    二维分离SAR成像算法
%    作者：Alex_Zikun
% 
%     介绍：先根据点目标分布，计算出对应的延时，再根据表达式计算回波数据进行仿真。再对回波数据
%    进行距离向和方位向上的脉冲压缩，得出二维图像。
% 
%    实验要求记录：
%     1.二维回波信号幅度、相位
%     2.距离向脉冲压缩结果的二维等高线图
%     3.点目标成像结果；二位等高线图，距离和方位剖面
%     4.用Hamming窗抑制成像结果副瓣
%%  基本参数和配置 
clc;clear all;close all;
    v_c =3e+8;%光速
    T =10e-6;%发射脉冲时间
    Br=60e6;%距离向带宽
    lamda=0.03;%波长
    f0=v_c/lamda;%载频
    vx=150;%雷达平台运动速度
    R0=15e3;%场景中心最短斜距
    Kr=Br/T;%调频斜率
    Nr=2048;%距离向采样点数,必须要大于T*Br
    Fr=100e6;%距离向采样频率
    deta_t=1/Fr;%距离向采样时间间隔
    tr=2*R0/v_c+((0:Nr-1)-Nr/2)*deta_t;%距离向采样时间轴
    fr=((-Nr/2):(Nr/2-1))/Nr*Fr;
    PRF=100;%PRF
    Na=60;%方位向采样点数
    ta=((0:(Na-1))-Na/2)/PRF;
    T_sar=Na/PRF;%方位向采样时间
    fa=((-Na/2):(Na/2-1))/Na*PRF;
    Ka=2*vx^2/(lamda*R0);
     
%% 根据目标点计算二维回波信号
    %%计算放置点的参数
    dot_num_a=1; % 方位向点个数
    deta_a=100; % 方位向点间距
    dot_num_r=3; % 距离向点个数
    deta_r=600; % 距离向点间距
    dot_xy_cell=cell(1,dot_num_a);
    
    middle_point_r=ceil(dot_num_r/2);
    middle_point_a=ceil(dot_num_a/2);
    line_x=vx*ta;
    line_y=zeros(1,Na);
    
    for i_dot_num_a=1:dot_num_a
        dot_xy=zeros(dot_num_r,2);
     for i_dot_num_r=1:dot_num_r
      dot_xy(i_dot_num_r,2)=(i_dot_num_r-middle_point_r)*deta_r;
      dot_xy(i_dot_num_r,1)=(i_dot_num_a-middle_point_a)*deta_a;
     end
     dot_xy_cell{1,i_dot_num_a}=dot_xy;
    end

    slant_range_cell=cell(1,dot_num_a);
    %计算每个点在所有方位的斜距
    for i_dot_num_a=1:dot_num_a
        slant_range=zeros(dot_num_r,Na);%single
        dot_xy=dot_xy_cell{1,i_dot_num_a};
        for i_dot_num_r=1:dot_num_r
        slant_range(i_dot_num_r,:)=sqrt((line_y-(R0+dot_xy(i_dot_num_r,2))).^2+(line_x-dot_xy(i_dot_num_r,1)).^2);%???
        end
        slant_range_cell{1,i_dot_num_a}=slant_range;
    end
    %计算每个点在所有方位的时延
    t_delay_cell=cell(1,dot_num_a);
    for i_dot_num_a=1:dot_num_a
         slant_range=slant_range_cell{1,i_dot_num_a};
         t_delay=slant_range*2/v_c;%single
         t_delay_cell{1,i_dot_num_a}=t_delay;
    end
    echo_data=zeros(Nr,Na);
    tr_matrix=repmat(tr.',1,Na);%矩阵化处理
    for i_dot_num_a=1:dot_num_a
        t_delay= t_delay_cell{1,i_dot_num_a};
       for i_dot_num_r=1:dot_num_r
       t_delay_matrix=repmat(t_delay(i_dot_num_r,:),Nr,1);%矩阵化处理
       echo_data=echo_data+rect((tr_matrix-t_delay_matrix)/T).*exp(1j*pi*Kr*(tr_matrix-t_delay_matrix).^2).*exp(-1j*2*pi*f0*t_delay_matrix);
       end
    end
    echo_data_phase=unwrap(angle(echo_data));
% 输出图像
    figure;
    subplot(2,1,1)
    imagesc(abs(echo_data));title('二维信号的幅度');
    subplot(2,1,2)
    imagesc(echo_data_phase);title('二维信号的相位');
 %% 脉冲压缩
 %距离向脉冲压缩
    Echo=circshift(fft2(circshift(echo_data,[-Nr/2,-Na/2])),[Nr/2,Na/2]);
    ref=exp(1j*pi*fr.^2/Kr);
    ref_matrix=repmat(ref.',1,Na);% 距离向二维匹配滤波器
    COMP=Echo.*ref_matrix; % 距离向匹配滤波后频域

 %方位向脉冲压缩
    fr_matrix=repmat(fr.',1,Na);
    fa_matrix=repmat(fa,Nr,1);
    Haz=exp(-1j*pi.*(fa.^2)./Ka);
    Haz_matrix=repmat(Haz,Nr,1); % 方位向二维匹配滤波器
    SAC=COMP.*Haz_matrix; % 再经方位向脉冲压缩后
    
%     [a b]=find(SAC==max(max(SAC))); % 找出最大值
%     figure;plot(abs(SAC(:,31)));
%     figure;plot(abs(SAC(1599,:)));
    sac=circshift(ifft2(circshift(SAC,[0,0])),[Nr/2,Na/2]); % 时域图像
    figure;
    subplot(2,1,1)
    imagesc(abs(sac));title('点目标成像时域结果')
    subplot(2,1,2)
    imagesc(abs(SAC));title('点目标成像频域结果')
%% 生成Hamming窗
    window=hamming(1151)*hamming(37).';
    window_hamming=zeros(2048,60);
    xx=10;
    window_hamming(450+xx:1600+xx,13:49)=window;
    WINDOW_data=window_hamming.*SAC;
    window_data=circshift(ifft2(circshift(WINDOW_data,[0,0])),[Nr/2,Na/2]);
    figure;
    subplot(2,1,1)
    imagesc(abs(window_data));title('加Hamming窗后时域')
    subplot(2,1,2)
    imagesc(abs(WINDOW_data));title('加Hamming窗后频域')
%% 剖面图
    i_x=1; % 第i个点
    i_y=1;

    x=1025+300*(i_x-1);
    y=31+400*(i_y-1);
    pou=sac(x-10:x+9,y-10:y+9); % 截取
    pou_range=pou(1:20,11);
    pou_azimuth=pou(11,1:20);
    range_db=ifft([fft(pou_range);zeros(980,1)]);
    tr_db=linspace(tr(x-10),tr(x+10),1000);

    azimuth_db=ifft([fft((pou_azimuth).');zeros(980,1)]);
    ta_db=linspace(ta(y-10),ta(y+10),1000);

    POU=fft2(pou);
    POU1=zeros(1000,1000);
    POU1(1:20,1:20)=POU;
    pou1=ifft2(POU1);
    a=max(abs(pou1));
    b=max(a);
    pou1_db=20*log10(abs(pou1)/b);
    c=max(pou1_db);
    dgx1=[-3 -13 -20 -30];
    figure;
    subplot(2,1,1)
    contour(pou1_db,dgx1);
    title('-3 -13 -20 -30(db)不加窗时二维图像');
    xlabel('方位向(m)');ylabel('距离向(m)');
% 加窗
    pou2=window_data(x-10:x+9,y-10:y+9);

    pou_range2=pou2(1:20,11);
    pou_azimuth2=pou2(11,1:20);
    range_db2=ifft([fft(pou_range2);zeros(980,1)]);
    tr_db=linspace(tr(x-10),tr(x+10),1000);

    azimuth_db2=ifft([fft((pou_azimuth2).');zeros(980,1)]);
    ta_db=linspace(ta(y-10),ta(y+10),1000);

    POU=fft2(pou2);
    POU1=zeros(1000,1000);
    POU1(1:20,1:20)=POU;
    pou1=ifft2(POU1);
    a=max(abs(pou1));
    b=max(a);
    pou1_db=20*log10(abs(pou1)/b);
    c=max(pou1_db);
    dgx2=[-3 -13 -20 -30];
    subplot(2,1,2)
    contour(pou1_db,dgx2);
    title('-3 -13 -20 -30(db)加窗时二维图像');
    xlabel('方位向(m)');ylabel('距离向(m)');
    
% 加窗对比输出
    figure;
    subplot(2,2,1)
    plot(tr_db,20*log10(abs(range_db)/max(abs(range_db))));
    title('不加窗的距离像')
    subplot(2,2,2)
    plot(ta_db,20*log10(abs(azimuth_db)/max(abs(azimuth_db))));
    title('不加窗的方位像')
    subplot(2,2,3)
    plot(tr_db,20*log10(abs(range_db2)/max(abs(range_db2))));
    title('加窗的距离像')
    subplot(2,2,4)
    plot(ta_db,20*log10(abs(azimuth_db2)/max(abs(azimuth_db2))));
    title('加窗的方位像')
    
   






