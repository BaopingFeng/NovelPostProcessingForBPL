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
n=10;
N=2^n;                        
R=1/2;                        
K=fix(N*R);       
max_iternum = 200;     
% CRC is from 3gpp 5G standard CRC24C.     
crc_m=24;                                   
crcgen = crc.generator('Polynomial' ,[1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1], ...      % CRC24C generator
                           'InitialState', [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], ...
                           'FinalXOR',  [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]) ;
crcdet = crc.detector('Polynomial' ,[1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1], ...      % CRC24C detector
                           'InitialState', [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0], ...
                           'FinalXOR',[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]) ;
poly = [1 1 0 1 1 0 0 1 0 1 0 1 1 0 0 0 1 0 0 0 1 0 1 1 1]; 
info_k=K-crc_m;                                         % crc is included in K.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=4;                                                            % simulate 4 results.
minEbn0=1.5;                                            
maxEbn0=minEbn0+(j-1)*0.5;
Ebn0 = linspace(minEbn0,maxEbn0,j);
Esn0=Ebn0+10*log10(R);                            % BPSK modulation
SNR=Esn0-10*log10(0.5);                            % real input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_error=100;                      
max_frame=10000000;       
min_frame=1000;             
be=zeros(1,j);      
ber=zeros(1,j);
fe=zeros(1,j);
fer=zeros(1,j);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[free_sets,frozen_sets]=PWR15(N,K);           % Channel Reliability based on 5G standard.
frozen = zeros(1,N);
frozen(1,free_sets) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Load the permutation set     
load('FactorGraphSetL32.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     6¡¢Simulation
for step=1:j 
    sigma = sqrt(1/((10^(SNR(step)/10))));               % sigma
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
        for irow = 1:row           
            %===============================
            index = permutation(irow,:);  
            llrPermute = llr(index);                       % permutation map
            froPermute = frozen(index);
            %===============================
            [v,flg] =  NBPL(llrPermute, froPermute, N, max_iternum,poly,crc_m,index,frozen);
            if flg(1) == 0
                diff = 0;
                isTerm = 1;
                break;
            end
        end
        if isTerm == 0
            diff = N-sum(v==uu);
        end
        be(step)=be(step)+diff;        
        i=i+1;    
        if diff~=0
            fe(step)=fe(step)+1;
        end
    end
    ber(step)=be(step)./(K*i);
    fer(step)=fe(step)/(i);    
end
semilogy(Ebn0,fer,'r-x','linewidth',2);
grid on;
axis([1 3.5 0.000001 1]);
xlabel('E_b/N_0(dB)');
ylabel('FER');
save BPL32N1024.mat