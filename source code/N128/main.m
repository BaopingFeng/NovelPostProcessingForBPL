%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    see paper "A Novel Post-Processing Method for BPL decoding of Polar Codes"
%    submit to IEEE Communication Letter, 2021.1.5
%    fengbaoping@buaa.edu.cn
%
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
rng('shuffle');
n=7;
N=2^n;                        
R=1/2;                        
K=fix(N*R);       
max_iternum = 200;     
% CRC is from 3gpp 5G standard CRC24C.     
crc_m=6;                                   
crcgen = crc.generator('Polynomial' ,[1 1 0 0 0 0 1], ...      % CRC6 generator
                           'InitialState', [0 0 0 0 0 0], ...
                           'FinalXOR',  [0 0 0 0 0 0]) ;
crcdet = crc.detector('Polynomial' ,[1 1 0 0 0 0 1], ...      % CRC6 detector
                           'InitialState', [0 0 0 0 0 0], ...
                           'FinalXOR',[0 0 0 0 0 0]) ;
poly = [1 1 0 0 0 0 1]; 
info_k=K;                                         % crc is included in K.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=7;                                                            % simulate 4 results.
minEbn0=2;                                            
maxEbn0=minEbn0+(j-1)*0.5;
Ebn0 = linspace(minEbn0,maxEbn0,j);
Esn0=Ebn0+10*log10(R);                            % BPSK modulation
SNR=Esn0-10*log10(0.5);                            % real input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_error=100;                      
max_frame=1000000;       
min_frame=1000;             
be=zeros(1,j);      
ber=zeros(1,j);
fe=zeros(1,j);
fer=zeros(1,j);
far=zeros(1,j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[free_sets,frozen_sets]=PWR15(N,K+crc_m);           % Channel Reliability based on 5G standard.
frozen = zeros(1,N);
frozen(1,free_sets) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Load the permutation set     
load('FactorGraphSetN128L32.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     6¡¢Simulation
for step=1:j 
    sigma = sqrt(1/((10^(SNR(step)/10))));               % noise variance
    i=0;   
    while ((fe(step) < frame_error) && (i < max_frame))||(i<min_frame)    
        u=randsrc(info_k,1,[0 1]);
        data_in=generate(crcgen,u);     
        uu=zeros(N,1);
        uu(free_sets)=data_in;
        % Encoding
        x = encode(uu);
        %AWGN channel
        channel_noise = sigma*randn(N,1);
        y0 = (1-2*x) + channel_noise;
        % Decoding
        llr=2*y0/(sigma*sigma);
        isTerm = 0;
        diff = 0;
        row = 32;              %size of factor graph set  /16
        ED = zeros(1,row);
        vL = zeros(N,row);
        for irow = 1:row           
            %===============================
            index = FactorGraphIndex(irow,:);  
            llrPermute = llr(index);                       % permutation map
            froPermute = frozen(index);
            %===============================
            [v,flg] =  NBPL(llrPermute, froPermute, N, max_iternum,poly,crc_m,index,frozen);
            vL(:,irow) = v;
            if flg(1) == 0            
                vx = encode(v);
                yx = 1-2*vx;                        
                for vi = 1:N                               % Compute Euclidean Distance
                    ED(1,irow) = ED(1,irow) + (yx(vi)-y0(vi))*(yx(vi)-y0(vi));
                end
            else
                ED(1,irow)=10000;
            end
        end
        [A,B] = sort(ED);
        v = vL(:,B(1));
        diff = N-sum(v==uu);
        be(step)=be(step)+diff;        
        i=i+1;   
        if diff~=0
            fe(step)=fe(step)+1;
        end
        [i fe(step)]
    end
    ber(step)=be(step)./(K*i);
    fer(step)=fe(step)/(i); 
end
semilogy(Ebn0,fer,'r-x','linewidth',2);
grid on;
axis([1 5 0.000001 1]);
xlabel('E_b/N_0(dB)');
ylabel('FER');