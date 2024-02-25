clc
clear 
close all
fprintf('This program started at %s\n', datestr(now,'HH:MM:SS'))
SNR_dB = -5:5:35;
a = [0.1 0.2 0.7];
PL = 10.^([0 6 12]./10);
N = 3;
thresh = 0.5;
PrCase1 = (1-thresh).*(1-thresh).*(1-thresh);
PrCase2 = (1-thresh).*(1-thresh).*thresh;
PrCase3 = (1-thresh).*thresh.*(1-thresh);
PrCase4 = (1-thresh).*thresh.*thresh;
PrCase5 = thresh.*(1-thresh).*(1-thresh);
PrCase6 = thresh.*(1-thresh).*thresh;
PrCase7 = thresh.*(1-thresh).*(1-thresh);
PrCase8 = thresh.*thresh.*thresh;
PrCases = [PrCase1 PrCase2 PrCase3 PrCase4 ...
    PrCase5 PrCase6 PrCase7 PrCase8];
a_su = 1;
Ps_dBm = 30;
Ps = 10^((Ps_dBm - 30)/10);
SNR = 10.^(SNR_dB./10);
% PL(1) = 1;
% PL(2) = 3.9811;
% PL(3) = 15.8489;
Va1 = [0.15/2 a(1) 0.050 1/3];
Va2 = [0.15/2 a(2) 0.475 1/3];
Va3 = [0.8500 a(3) 0.475 1/3];
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
SVa3 = sqrt(Va3);
Iter = 50000;
Data_length = 128;

%Lookup table for C-NOMA
MLD_LUT_C = [(1*SVa1(2)+1*SVa2(2)+1*SVa3(2)) ...
    (1*SVa1(2)+1*SVa2(2)-1*SVa3(2)) ...
    (1*SVa1(2)-1*SVa2(2)+1*SVa3(2)) ...
    (1*SVa1(2)-1*SVa2(2)-1*SVa3(2)) ...
    (-1*SVa1(2)+1*SVa2(2)+1*SVa3(2)) ...
    (-1*SVa1(2)+1*SVa2(2)-1*SVa3(2)) ...
    (-1*SVa1(2)-1*SVa2(2)+1*SVa3(2)) ...
    (-1*SVa1(2)-1*SVa2(2)-1*SVa3(2))];
PrC_normFact = MLD_LUT_C.*MLD_LUT_C*PrCases';
MLD_LUT_C = MLD_LUT_C/sqrt(PrC_normFact);
%Lookup table for PANOMA
MLD_LUT_P = [(1*SVa1(4)+1*SVa2(4)+1*SVa3(4)) ...
    (1*SVa1(1)+1*SVa2(1)-1*SVa3(1)) ...
    (1*SVa1(2)-1*SVa2(2)+1*SVa3(2)) ...
    (1*SVa1(3)-1*SVa2(3)-1*SVa3(3)) ...
    (-1*SVa1(3)+1*SVa2(3)+1*SVa3(3)) ...
    (-1*SVa1(2)+1*SVa2(2)-1*SVa3(2)) ...
    (-1*SVa1(1)-1*SVa2(1)+1*SVa3(1)) ...
    (-1*SVa1(4)-1*SVa2(4)-1*SVa3(4))];
PrP_normFact = MLD_LUT_P.*MLD_LUT_P*PrCases';
MLD_LUT_P = MLD_LUT_P/sqrt(PrP_normFact);
data1_LUT = [1 1 1 1 0 0 0 0];
data2_LUT = [1 1 0 0 1 1 0 0];
data3_LUT = [1 0 1 0 1 0 1 0];

for i_snr = 1:size(SNR_dB,2)
    i_block = 0;
    bit_error1_MLD = 0;
    bit_error2_MLD = 0;
    bit_error3_MLD = 0;
    bit_error1_MLD_PANOMA = 0;
    bit_error2_MLD_PANOMA = 0;
    bit_error3_MLD_PANOMA = 0;
    bit_error_su = 0;
    
    while(i_block<=Iter)
        %data generation
        data1 = (rand(1, Data_length)>thresh);
        data2 = (rand(1, Data_length)>thresh);
        data3 = (rand(1, Data_length)>thresh);
        data_su = (rand(1, Data_length)>0.5);
        s1 = 2.*data1 - 1;
        s2 = 2.*data2 - 1;
        s3 = 2.*data3 - 1;
        s_su = 2.*data_su - 1;
        DataComparison = ((data1==data3)&(data1==data2)&(data2==data3));
        DataComparison = DataComparison + (((data1~=data3)&(data1~=data2)&(data2==data3))|((data1==data3)&(data1==data2)&(data2==data3)));
        DataComparison = DataComparison + (((data1==data3)&(data1~=data2)&(data2~=data3))|((data1~=data3)&(data1~=data2)&(data2==data3))|((data1==data3)&(data1==data2)&(data2==data3)));
        DataComparison = DataComparison + (((data1~=data3)&(data1==data2)&(data2~=data3))|((data1==data3)&(data1~=data2)&(data2~=data3))|((data1~=data3)&(data1~=data2)&(data2==data3))|((data1==data3)&(data1==data2)&(data2==data3)));
        
        AdaptivePowerIndex = DataComparison; %1: 1&2 are equal; PANOMA
        A1 = [];                             %2: all different; conventional NOMA
        A2 = [];                             %3: 2&3 are equal; PANOMA
        A3 = [];                             %4: all equal; PANOMA
        for i=1:1:length(AdaptivePowerIndex)
            A1(i) = Va1(AdaptivePowerIndex(i));
            A2(i) = Va2(AdaptivePowerIndex(i));
            A3(i) = Va3(AdaptivePowerIndex(i));
        end
        %superposition coding
        x = sqrt(a(1)*Ps)*s1 + sqrt(a(2)*Ps)*s2 + sqrt(a(3)*Ps)*s3;
        x = x/sqrt(PrC_normFact);
        x_su = sqrt(a_su*Ps)*s_su;
        x_PANOMA = sqrt(A1*Ps).*s1 + sqrt(A2*Ps).*s2 + sqrt(A3*Ps).*s3;
        x_PANOMA = x_PANOMA/sqrt(PrP_normFact);
        %Rayleigh fading channel
        h_real = randn(1,3);
        h_img = randn(1,3);
        h = 1/sqrt(2).*(h_real + 1j.*h_img);
        h1 = 1/sqrt(PL(1))*h(1,1);
        h2 = 1/sqrt(PL(2))*h(1,2);
        h3 = 1/sqrt(PL(3))*h(1,2);
        h_su = 1*h(1,1);
        
        %AWGN
        Sf1 = sqrt(1/(2*N*SNR(i_snr)));
        Sf2 = sqrt(1/(2*N*SNR(i_snr)));
        Sf3 = sqrt(1/(2*N*SNR(i_snr)));
        Sf_su = sqrt(a_su/(2*SNR(i_snr)));
        n1 = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf1;
        n2 = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf2;
        n3 = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf3;
        n_su = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf_su;
        
        %Received signal
        y1 = h1*x + n1;
        y2 = h2*x + n2;
        y3 = h3*x + n3;
        y_su = h_su*x_su + n_su;
        y1_PANOMA = h1*x_PANOMA + n1;
        y2_PANOMA = h2*x_PANOMA + n2;
        y3_PANOMA = h3*x_PANOMA + n3;
        
        %Symbol detection and mapping into bits
        for i=1:1:length(MLD_LUT_C)
            MSE1(i,:) = abs(y1/h1 - MLD_LUT_C(i));
            MSE2(i,:) = abs(y2/h2 - MLD_LUT_C(i));
            MSE3(i,:) = abs(y3/h3 - MLD_LUT_C(i));
            MSE1_PANOMA(i,:) = abs(y1_PANOMA/h1 - MLD_LUT_P(1,i));
            MSE2_PANOMA(i,:) = abs(y2_PANOMA/h2 - MLD_LUT_P(1,i));
            MSE3_PANOMA(i,:) = abs(y3_PANOMA/h3 - MLD_LUT_P(1,i));
        end
        [min_metric1 decis1] = min(MSE1);
        [min_metric2 decis2] = min(MSE2);
        [min_metric3 decis3] = min(MSE3);
        [min_metric1_PANOMA decis1_PANOMA] = min(MSE1_PANOMA);
        [min_metric2_PANOMA decis2_PANOMA] = min(MSE2_PANOMA);
        [min_metric3_PANOMA decis3_PANOMA] = min(MSE3_PANOMA);
        data1_detected_MLD(1,:) = data1_LUT(decis1(1,:));
        data2_detected_MLD(1,:) = data2_LUT(decis2(1,:));
        data3_detected_MLD(1,:) = data3_LUT(decis3(1,:));
        data_su_detected = (y_su/h_su)>0;
        data1_detected_MLD_PANOMA(1,:) = data1_LUT(decis1_PANOMA(1,:));
        data2_detected_MLD_PANOMA(1,:) = data2_LUT(decis2_PANOMA(1,:));
        data3_detected_MLD_PANOMA(1,:) = data3_LUT(decis3_PANOMA(1,:));
        %Error counting
        bit_error_su = bit_error_su+sum(data_su_detected~=data_su);
        bit_error1_MLD=bit_error1_MLD+sum(data1_detected_MLD~=data1);
        bit_error2_MLD=bit_error2_MLD+sum(data2_detected_MLD~=data2);
        bit_error3_MLD=bit_error3_MLD+sum(data3_detected_MLD~=data3);
        bit_error1_MLD_PANOMA=bit_error1_MLD_PANOMA+sum(data1_detected_MLD_PANOMA~=data1);
        bit_error2_MLD_PANOMA=bit_error2_MLD_PANOMA+sum(data2_detected_MLD_PANOMA~=data2);
        bit_error3_MLD_PANOMA=bit_error3_MLD_PANOMA+sum(data3_detected_MLD_PANOMA~=data3);
        i_block=i_block+1;
    end
    %BER calculation
    BER1(i_snr)=bit_error1_MLD/Data_length/i_block;
    BER2(i_snr)=bit_error2_MLD/Data_length/i_block;
    BER3(i_snr)=bit_error3_MLD/Data_length/i_block;
    BER1_PANOMA(i_snr)=bit_error1_MLD_PANOMA/Data_length/i_block;
    BER2_PANOMA(i_snr)=bit_error2_MLD_PANOMA/Data_length/i_block;
    BER3_PANOMA(i_snr)=bit_error3_MLD_PANOMA/Data_length/i_block;
    BER_su(i_snr) = bit_error_su/Data_length/i_block;
    BER_avg(i_snr) = (BER1(i_snr) + BER2(i_snr) + BER3(i_snr))/3;
    BER_PANOMA_avg(i_snr) = (BER1_PANOMA(i_snr) + BER2_PANOMA(i_snr) + BER3_PANOMA(i_snr))/3;
    cap1(i_snr)=-((1-BER1(i_snr))*log2(1-BER1(i_snr))) ...
        - ((BER1(i_snr))*log2(BER1(i_snr)));
    cap1(i_snr)= 1 - cap1(i_snr);
    cap2(i_snr)=-((1-BER2(i_snr))*log2(1-BER2(i_snr))) ...
        - ((BER2(i_snr))*log2(BER2(i_snr)));
    cap2(i_snr)= 1 - cap2(i_snr);
    cap3(i_snr)=-((1-BER3(i_snr))*log2(1-BER3(i_snr))) ...
        - ((BER3(i_snr))*log2(BER3(i_snr)));
    cap3(i_snr)= 1 - cap3(i_snr);
    sumRate(i_snr) = cap1(i_snr) + cap2(i_snr) + cap3(i_snr);
    cap1_PANOMA(i_snr)=-((1-BER1_PANOMA(i_snr))*log2(1-BER1_PANOMA(i_snr))) ...
        - ((BER1_PANOMA(i_snr))*log2(BER1_PANOMA(i_snr)));
    cap1_PANOMA(i_snr)= 1 - cap1_PANOMA(i_snr);
    cap2_PANOMA(i_snr)=-((1-BER2_PANOMA(i_snr))*log2(1-BER2_PANOMA(i_snr))) ...
        - ((BER2_PANOMA(i_snr))*log2(BER2_PANOMA(i_snr)));
    cap2_PANOMA(i_snr)= 1 - cap2_PANOMA(i_snr);
    cap3_PANOMA(i_snr)=-((1-BER3_PANOMA(i_snr))*log2(1-BER3_PANOMA(i_snr))) ...
        - ((BER3_PANOMA(i_snr))*log2(BER3_PANOMA(i_snr)));
    cap3_PANOMA(i_snr)= 1 - cap3_PANOMA(i_snr);
    sumRate_PANOMA(i_snr) = cap1_PANOMA(i_snr) + cap2_PANOMA(i_snr) + cap3_PANOMA(i_snr);
    ca_su(i_snr)=-((1-BER_su(i_snr))*log2(1-BER_su(i_snr))) ...
        - ((BER_su(i_snr))*log2(BER_su(i_snr)));
    ca_su(i_snr)= 1 - ca_su(i_snr);
    
    tber = 0.5*(1 - sqrt(SNR./(1 + SNR)));
end
