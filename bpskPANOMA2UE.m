clc 
clear
close all
fprintf('This program started at %s\n', datestr(now,'HH:MM:SS'))

a = [0.05 0.95];
PL = 10.^([0 6]/10);
SNR_dB = -5:5:25;
N = 2;
thresh = 0.5;
PrCI = (1-thresh).*(1-thresh) + (thresh).*(thresh);
PrCase1 = (1-thresh).*(1-thresh);
PrCase2 = (1-thresh).*thresh;
PrCase3 = (1-thresh).*thresh;
PrCase4 = thresh.*thresh;
PrCases = [PrCase1 PrCase2 PrCase3 PrCase4];
a_su = 1;
Ps_dBm = 30;
Ps = 10^((Ps_dBm - 30)/10);
SNR = 10.^(SNR_dB./10);
Va1 = [a(1) 0.5];
Va2 = [a(2) 0.5];
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
Iter = 300000;
Data_length = 128;

%Lookup table for PANOMA
MLD_LUT_PANOMA = sqrt(Ps)*[(1*SVa1(2)+1*SVa2(2)) ...
    (1*SVa1(1)-1*SVa2(1)) ...
    (-1*SVa1(1)+1*SVa2(1)) ...
    (-1*SVa1(2)-1*SVa2(2))];
PrP_normFact = MLD_LUT_PANOMA.*MLD_LUT_PANOMA*PrCases';
MLD_LUT_PANOMA = MLD_LUT_PANOMA/sqrt(PrP_normFact);

%Lookup table for C-NOMA
MLD_LUT_CNOMA = sqrt(Ps)*[(1*sqrt(a(1))+1*sqrt(a(2))) ...
    (1*sqrt(a(1))-1*sqrt(a(2))) ...
    (-1*sqrt(a(1))+1*sqrt(a(2))) ...
    (-1*sqrt(a(1))-1*sqrt(a(2)))];
PrC_normFact = MLD_LUT_CNOMA.*MLD_LUT_CNOMA*PrCases';
MLD_LUT_CNOMA = MLD_LUT_CNOMA/sqrt(PrC_normFact);
data1_LUT = [1 1 0 0];
data2_LUT = [1 0 1 0];

for i_snr = 1:size(SNR_dB,2)
    i_block = 0;
    bit_error1_MLD = 0;
    bit_error2_MLD = 0;
    bit_error1_MLD_PANOMA = 0;
    bit_error2_MLD_PANOMA = 0;
    c1 = 0;
    c2 = 0;
    c1_PANOMA = 0;
    c2_PANOMA = 0;
    bit_error_su = 0;
    
    while(i_block<=Iter)
        %data generation
        data1 = (rand(1, Data_length)>thresh);
        data2 = (rand(1, Data_length)>thresh);
        data_su = (rand(1, Data_length)>0.5);
        s1 = 2.*data1 - 1;
        s2 = 2.*data2 - 1;
        s_su = 2.*data_su - 1;
        DataComparison = (data1==data2);
        AdaptivePowerIndex = DataComparison + 1; %1: different power & 2: equal power
        A1 = [];
        A2 = [];
        for i=1:1:length(AdaptivePowerIndex)
            A1(i) = Va1(AdaptivePowerIndex(i));
            A2(i) = Va2(AdaptivePowerIndex(i));
        end
        
        %superposition coding
        x = sqrt(a(1)*Ps)*s1 + sqrt(a(2)*Ps)*s2;
        x = x/sqrt(PrC_normFact);
        x_su = sqrt(a_su*Ps)*s_su;
        x_PANOMA = (sqrt(A1*Ps).*s1 + sqrt(A2*Ps).*s2);
        x_PANOMA = x_PANOMA/sqrt(PrP_normFact);
        %Rayleigh fading channel
        h_real = randn(1,2);
        h_img = randn(1,2);
        h = 1/sqrt(2).*(h_real + 1j.*h_img);
        h1 = 1/sqrt(PL(1))*h(1,1);
        h2 = 1/sqrt(PL(2))*h(1,2);
        h_su = 1*h(1,1);
        
        %AWGN
        Sf1 = sqrt(1/(2*N*SNR(i_snr)));
        Sf2 = sqrt(1/(2*N*SNR(i_snr)));
        Sf_su = sqrt(a_su/(2*SNR(i_snr)));
        n1 = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf1;
        n2 = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf2;
        n_su = (randn(1, Data_length)+1j*randn(1, Data_length))*Sf_su;
        
        %Received signal
        y1 = h1*x + n1;
        y2 = h2*x + n2;
        y_su = h_su*x_su + n_su;
        y1_PANOMA = h1*x_PANOMA + n1;
        y2_PANOMA = h2*x_PANOMA + n2;
        
        %Symbol detection and Bit mapping
        for i=1:1:4
            MSE1(i,:) = abs(y1/h1 - MLD_LUT_CNOMA(i));
            MSE2(i,:) = abs(y2/h2 - MLD_LUT_CNOMA(i));
            MSE1_PANOMA(i,:) = abs(y1_PANOMA/h1 - MLD_LUT_PANOMA(i));
            MSE2_PANOMA(i,:) = abs(y2_PANOMA/h2 - MLD_LUT_PANOMA(i));
            
        end
        [min_metric1 decis1] = min(MSE1);
        [min_metric2 decis2] = min(MSE2);
        [min_metric1_PANOMA decis1_PANOMA] = min(MSE1_PANOMA);
        [min_metric2_PANOMA decis2_PANOMA] = min(MSE2_PANOMA);
        
        data_su_detected = (y_su/h_su)>0;
        data1_detected_MLD(1,:) = data1_LUT(decis1(1,:));
        data2_detected_MLD(1,:) = data2_LUT(decis2(1,:));
        data1_detected_MLD_PANOMA(1,:) = data1_LUT(decis1_PANOMA(1,:));
        data2_detected_MLD_PANOMA(1,:) = data2_LUT(decis2_PANOMA(1,:));
        
        %Error counting
        bit_error_su=bit_error_su+sum(data_su_detected~=data_su);
        bit_error1_MLD=bit_error1_MLD+sum(data1_detected_MLD~=data1);
        bit_error2_MLD=bit_error2_MLD+sum(data2_detected_MLD~=data2);
        bit_error1_MLD_PANOMA=bit_error1_MLD_PANOMA+sum(data1_detected_MLD_PANOMA~=data1);
        bit_error2_MLD_PANOMA=bit_error2_MLD_PANOMA+sum(data2_detected_MLD_PANOMA~=data2);
        i_block=i_block+1;
    end
    %BER calculation
    BER_su(i_snr) = bit_error_su/Data_length/i_block;
    BER1_MLD(i_snr)=bit_error1_MLD/Data_length/i_block;
    BER2_MLD(i_snr)=bit_error2_MLD/Data_length/i_block;
    BER_avg(i_snr) = (BER1_MLD(i_snr) + BER2_MLD(i_snr))/2;
    BER1_MLD_PANOMA(i_snr)=bit_error1_MLD_PANOMA/Data_length/i_block;
    BER2_MLD_PANOMA(i_snr)=bit_error2_MLD_PANOMA/Data_length/i_block;
    BER_PANOMA_avg(i_snr) = (BER1_MLD_PANOMA(i_snr) + BER2_MLD_PANOMA(i_snr))/2;
    cap1(i_snr)=-((1-BER1_MLD(i_snr))*log2(1-BER1_MLD(i_snr))) ...
        - ((BER1_MLD(i_snr))*log2(BER1_MLD(i_snr)));
    cap1(i_snr)= 1 - cap1(i_snr);
    cap2(i_snr)=-((1-BER2_MLD(i_snr))*log2(1-BER2_MLD(i_snr))) ...
        - ((BER2_MLD(i_snr))*log2(BER2_MLD(i_snr)));
    cap2(i_snr)= 1 - cap2(i_snr);
    sumRate(i_snr) = cap1(i_snr) + cap2(i_snr);
    cap1_PANOMA(i_snr)=-((1-BER1_MLD_PANOMA(i_snr))*log2(1-BER1_MLD_PANOMA(i_snr))) ...
        - ((BER1_MLD_PANOMA(i_snr))*log2(BER1_MLD_PANOMA(i_snr)));
    cap1_PANOMA(i_snr)= 1 - cap1_PANOMA(i_snr);
    cap2_PANOMA(i_snr)=-((1-BER2_MLD_PANOMA(i_snr))*log2(1-BER2_MLD_PANOMA(i_snr))) ...
        - ((BER2_MLD_PANOMA(i_snr))*log2(BER2_MLD_PANOMA(i_snr)));
    cap2_PANOMA(i_snr)= 1 - cap2_PANOMA(i_snr);
    sumRate_PANOMA(i_snr) = cap1_PANOMA(i_snr) + cap2_PANOMA(i_snr);
    ca_su(i_snr)=-((1-BER_su(i_snr))*log2(1-BER_su(i_snr))) ...
        - ((BER_su(i_snr))*log2(BER_su(i_snr)));
    ca_su(i_snr)= 1 - ca_su(i_snr);
    tber = 0.5*(1 - sqrt(SNR./(1 + SNR)));
end

