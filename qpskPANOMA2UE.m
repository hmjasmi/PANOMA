% 04.08.2020
% Copyrights to The Universuty of Manchester
% Author: Hamad Yahya
% Ver2 including pathloss, normalization by the number of users and
% constellation normalization.
% This function is used to simulate NOMA technique to find the BER
% C-NOMA and PANOMA are considered
% N = 2 ==> u1 is NU while u2 is FU
% The input is the power allocation as a vector [a1 a2]
% SNR needs to be defined as a range or single value
% The number of iterations needs to be defined
% The channel model is the SISO broadcast Rayleigh flat fading
% Pathloss needs to be known
% QPSK is considered for both users

clc 
clear
close all
fprintf('This program started at %s\n', datestr(now,'HH:MM:SS'))

a = [0.05 0.95];
PL = 10.^([0 6]/10);
SNR_dB = -5:5:25;
Iter = 300000;
QPSK_map = sqrt(0.5)*[1+j,-1+j,1-j,-1-j];
dataMSB_LUT = [0 0 1 1];
dataLSB_LUT = [0 1 0 1];
qpskSym_key = [0 1 2 3];
Data_length = 256;              
Frame_length = Data_length/2;
a_su = 1;
N = 2;
f_su = 1;
Va1 = [a(1) 0.5];               %NU power coefficients vector
Va2 = [a(2) 0.5];               %FU power coefficients vector
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
SNR = 10.^(SNR_dB./10);
Ps_dBm = 30;                    %BS Tx power
Ps = 10^((Ps_dBm - 30)/10);

%Lookup table for C-NOMA
[qpskMesh1 qpskMesh2] = meshgrid(QPSK_map,QPSK_map);
[qpskSym_LUT1 qpskSym_LUT2] = meshgrid(qpskSym_key,qpskSym_key);
Data1LUT = de2bi(qpskSym_LUT1,'left-msb').';
Data2LUT = de2bi(qpskSym_LUT2,'left-msb').';
MLD_LUT_C = SVa1(1)*qpskMesh1 + SVa2(1)*qpskMesh2;
MLD_LUT_C = MLD_LUT_C(:).';
compReal = (((real(qpskMesh1(:)))) == ((real(qpskMesh2(:)))))+1;
compImag = (((imag(qpskMesh1(:)))) == ((imag(qpskMesh2(:)))))+1;
symbol1_LUT = qpskSym_LUT1(:).'+1;
symbol2_LUT = qpskSym_LUT2(:).'+1;

MLD_LUT_P = SVa1(compReal).*real(qpskMesh1(:)).' + SVa2(compReal).*real(qpskMesh2(:).') ...
           + 1j*SVa1(compImag).*imag(qpskMesh1(:)).' + 1j*SVa2(compImag).*imag(qpskMesh2(:).');
normFactor_P = MLD_LUT_P*MLD_LUT_P'/(16);
MLD_LUT_P = MLD_LUT_P/sqrt(normFactor_P);

%loop for different SNR values
parfor i_snr = 1:size(SNR_dB,2)
    tic
    i_block = 0;
    bit_error_su1 = 0;
    bit_error_su2 = 0;
    bit_error1 = 0;                %NU counter for bits in error
    bit_error2 = 0;                %FU counter for bits in error
    bit_error1_P = 0;
    bit_error2_P = 0;
    %User1 (PANOMA)
    ser_1_PANOMA = 0;correct1_1_PANOMA = 0;correct2_1_PANOMA = 0;correct3_1_PANOMA = 0;correct4_1_PANOMA = 0;
    ser12_1_PANOMA = 0;ser13_1_PANOMA = 0;ser14_1_PANOMA = 0;ser21_1_PANOMA = 0;ser23_1_PANOMA = 0;ser24_1_PANOMA = 0;
    ser31_1_PANOMA = 0;ser32_1_PANOMA = 0;ser34_1_PANOMA = 0;ser41_1_PANOMA = 0;ser42_1_PANOMA = 0;ser43_1_PANOMA = 0;
    %User2 (PANOMA)
    ser_2_PANOMA = 0;correct1_2_PANOMA = 0;correct2_2_PANOMA = 0;correct3_2_PANOMA = 0;correct4_2_PANOMA = 0;
    ser12_2_PANOMA = 0;ser13_2_PANOMA = 0;ser14_2_PANOMA = 0;ser21_2_PANOMA = 0;ser23_2_PANOMA = 0;ser24_2_PANOMA = 0;
    ser31_2_PANOMA = 0;ser32_2_PANOMA = 0;ser34_2_PANOMA = 0;ser41_2_PANOMA = 0;ser42_2_PANOMA = 0;ser43_2_PANOMA = 0;    %Monte-Carlo Simulation
    while (i_block<=Iter)
        
        %Data generation and symbol mapping
        data1 = (rand(1, Data_length)>0.5);
        data2 = (rand(1, Data_length)>0.5);
        qpsk_code1 = reshape(data1, 2, Frame_length); %MSB is first row
        qpsk_code2 = reshape(data2, 2, Frame_length); %MSB is first row
        decimal1 = 2*qpsk_code1(1,:) + qpsk_code1(2,:) + 1;
        decimal2 = 2*qpsk_code2(1,:) + qpsk_code2(2,:) + 1;
        s1 = QPSK_map(2*qpsk_code1(1,:) + qpsk_code1(2,:) + 1);
        s2 = QPSK_map(2*qpsk_code2(1,:) + qpsk_code2(2,:) + 1);
        TxCompReal = (((real(s1))) == ((real(s2))))+1;
        TxCompImag = (((imag(s1))) == ((imag(s2))))+1;
        A1Real = Va1(TxCompReal);
        A2Real = Va2(TxCompReal);
        A1Imag = Va1(TxCompImag);
        A2Imag = Va2(TxCompImag);
        %Superposition Coding (SC)
        x = sqrt(a(1)*Ps)*s1 + sqrt(a(2)*Ps)*s2;
        x_PANOMA = sqrt(A1Real).*real(s1) + sqrt(A2Real).*real(s2) ...
                   + 1j*sqrt(A1Imag).*imag(s1) + 1j*sqrt(A2Imag).*imag(s2);
        x_PANOMA = x_PANOMA./sqrt(normFactor_P);       
        x_su1 = sqrt(a_su*Ps)*s1;
        x_su2 = sqrt(a_su*Ps)*s2;
        %Received signal at NU and FU
        h_real = randn(1,2);
        h_img = randn(1,2);
        h = 1/sqrt(2).*(h_real + 1j.*h_img);        
        h1 = 1/sqrt(PL(1))*h(1,1);
        h2 = 1/sqrt(PL(2))*h(1,2);

        Sf1 = sqrt(1/(4*N*SNR(i_snr)));
        Sf2 = sqrt(1/(4*N*SNR(i_snr)));
        Sf_su_BPSK = sqrt(1/(2*SNR(i_snr)));
        Sf_su_QPSK = sqrt(1/(2*SNR(i_snr)));
        n1 = (randn(1, Frame_length)+1j*randn(1, Frame_length))*Sf1;
        n2 = (randn(1, Frame_length)+1j*randn(1, Frame_length))*Sf2;
        y1 = h1*x + n1;
        y2 = h2*x + n2;
        y1_PANOMA = h1*x_PANOMA + n1;
        y2_PANOMA = h2*x_PANOMA + n2;
        y_su1 = h1*x_su1 + n1;
        y_su2 = h2*x_su2 + n2;
        
        %Symbol detection and mapping into bits
        MSE_su1 = abs(y_su1/h1 - QPSK_map.');
        MSE_su2 = abs(y_su2/h2 - QPSK_map.');
        MSE1 = abs(y1/h1 - MLD_LUT_C.');
        MSE2 = abs(y2/h2 - MLD_LUT_C.');
        MSE1_PANOMA = abs(y1_PANOMA/h1 - MLD_LUT_P.');
        MSE2_PANOMA = abs(y2_PANOMA/h2 - MLD_LUT_P.');
        [min_metric_su decis_su1] = min(MSE_su1);      
        [min_metric_su decis_su2] = min(MSE_su2);                 
        [min_metric1 decis1] = min(MSE1);                 
        [min_metric2 decis2] = min(MSE2);                 
        [min_metric1_P decis1_P] = min(MSE1_PANOMA);                 
        [min_metric2_P decis2_P] = min(MSE2_PANOMA);  
        %Symbol detection
        symbol1_detected_PANOMA = symbol1_LUT(decis1_P);
        symbol2_detected_PANOMA = symbol2_LUT(decis2_P);

        dataMSB_detected_su1 = dataMSB_LUT(decis_su1);
        dataLSB_detected_su1 = dataLSB_LUT(decis_su1);
        dataMSB_detected_su2 = dataMSB_LUT(decis_su2);
        dataLSB_detected_su2 = dataLSB_LUT(decis_su2);
        data1MSB_detected_MLD = Data1LUT(1,decis1);
        data1LSB_detected_MLD = Data1LUT(2,decis1);
        data2MSB_detected_MLD = Data2LUT(1,decis2);
        data2LSB_detected_MLD = Data2LUT(2,decis2);
        data1MSB_detected_P = Data1LUT(1,decis1_P);
        data1LSB_detected_P = Data1LUT(2,decis1_P);
        data2MSB_detected_P = Data2LUT(1,decis2_P);
        data2LSB_detected_P = Data2LUT(2,decis2_P);
        
        %Input & output comparison (symbol-based) ==> user1 (PANOMA)
        for i_pctLength = 1:Frame_length
            if (symbol1_detected_PANOMA(i_pctLength)==1)
                if(decimal1(i_pctLength)==1)
                correct1_1_PANOMA = correct1_1_PANOMA + 1;
                elseif(decimal1(i_pctLength)==2)
                    ser12_1_PANOMA = ser12_1_PANOMA + 1; %e1
                elseif(decimal1(i_pctLength)==3)
                    ser13_1_PANOMA = ser13_1_PANOMA + 1; %e2
                elseif(decimal1(i_pctLength)==4)
                    ser14_1_PANOMA = ser14_1_PANOMA + 1; %e4
                end
            elseif (symbol1_detected_PANOMA(i_pctLength)==2)
                if(decimal1(i_pctLength)==1)
                    ser21_1_PANOMA = ser21_1_PANOMA + 1; %e1
                elseif(decimal1(i_pctLength)==2)
                    correct2_1_PANOMA = correct2_1_PANOMA + 1;
                elseif(decimal1(i_pctLength)==3)
                    ser23_1_PANOMA = ser23_1_PANOMA + 1; %e3
                elseif(decimal1(i_pctLength)==4)
                    ser24_1_PANOMA = ser24_1_PANOMA + 1; %e2
                end
            elseif (symbol1_detected_PANOMA(i_pctLength)==3)
                if(decimal1(i_pctLength)==1)
                    ser31_1_PANOMA = ser31_1_PANOMA + 1; %e1
                elseif(decimal1(i_pctLength)==2)
                    ser32_1_PANOMA = ser32_1_PANOMA + 1; %e3
                elseif(decimal1(i_pctLength)==3)
                    correct3_1_PANOMA = correct3_1_PANOMA + 1;
                elseif(decimal1(i_pctLength)==4)
                    ser34_1_PANOMA = ser34_1_PANOMA + 1; %e2
                end
            elseif (symbol1_detected_PANOMA(i_pctLength)==4)
                if(decimal1(i_pctLength)==1)
                    ser41_1_PANOMA = ser41_1_PANOMA + 1; %e3
                elseif(decimal1(i_pctLength)==2)
                    ser42_1_PANOMA = ser42_1_PANOMA + 1; %e1
                elseif(decimal1(i_pctLength)==3)
                    ser43_1_PANOMA = ser43_1_PANOMA + 1; %e2
                elseif(decimal1(i_pctLength)==4)
                    correct4_1_PANOMA = correct4_1_PANOMA + 1;
                end
            end
        end
        %Input & output comparison (symbol-based) ==> user2 (PANOMA)
        for i_pctLength = 1:Frame_length
            if (symbol2_detected_PANOMA(i_pctLength)==1)
                if(decimal2(i_pctLength)==1)
                correct1_2_PANOMA = correct1_2_PANOMA + 1;
                elseif(decimal2(i_pctLength)==2)
                    ser12_2_PANOMA = ser12_2_PANOMA + 1; %e1
                elseif(decimal2(i_pctLength)==3)
                    ser13_2_PANOMA = ser13_2_PANOMA + 1; %e2
                elseif(decimal2(i_pctLength)==4)
                    ser14_2_PANOMA = ser14_2_PANOMA + 1; %e4
                end
            elseif (symbol2_detected_PANOMA(i_pctLength)==2)
                if(decimal2(i_pctLength)==1)
                    ser21_2_PANOMA = ser21_2_PANOMA + 1; %e1
                elseif(decimal2(i_pctLength)==2)
                    correct2_2_PANOMA = correct2_2_PANOMA + 1;
                elseif(decimal2(i_pctLength)==3)
                    ser23_2_PANOMA = ser23_2_PANOMA + 1; %e3
                elseif(decimal2(i_pctLength)==4)
                    ser24_2_PANOMA = ser24_2_PANOMA + 1; %e2
                end
            elseif (symbol2_detected_PANOMA(i_pctLength)==3)
                if(decimal2(i_pctLength)==1)
                    ser31_2_PANOMA = ser31_2_PANOMA + 1; %e1
                elseif(decimal2(i_pctLength)==2)
                    ser32_2_PANOMA = ser32_2_PANOMA + 1; %e3
                elseif(decimal2(i_pctLength)==3)
                    correct3_2_PANOMA = correct3_2_PANOMA + 1;
                elseif(decimal2(i_pctLength)==4)
                    ser34_2_PANOMA = ser34_2_PANOMA + 1; %e2
                end
            elseif (symbol2_detected_PANOMA(i_pctLength)==4)
                if(decimal2(i_pctLength)==1)
                    ser41_2_PANOMA = ser41_2_PANOMA + 1; %e3
                elseif(decimal2(i_pctLength)==2)
                    ser42_2_PANOMA = ser42_2_PANOMA + 1; %e1
                elseif(decimal2(i_pctLength)==3)
                    ser43_2_PANOMA = ser43_2_PANOMA + 1; %e2
                elseif(decimal2(i_pctLength)==4)
                    correct4_2_PANOMA = correct4_2_PANOMA + 1;
                end
            end
        end

        %Error counting
        bit_error_su1=bit_error_su1+sum(dataMSB_detected_su1~=qpsk_code1(1,:)) ...
            +sum(dataLSB_detected_su1~=qpsk_code1(2,:));
        bit_error_su2=bit_error_su2+sum(dataMSB_detected_su2~=qpsk_code2(1,:)) ...
            +sum(dataLSB_detected_su2~=qpsk_code2(2,:));
        bit_error1=bit_error1+sum(data1MSB_detected_MLD~=qpsk_code1(1,:)) ...
            +sum(data1LSB_detected_MLD~=qpsk_code1(2,:));
        bit_error1_P=bit_error1_P+sum(data1MSB_detected_P~=qpsk_code1(1,:)) ...
            +sum(data1LSB_detected_P~=qpsk_code1(2,:));
        bit_error2=bit_error2+sum(data2MSB_detected_MLD~=qpsk_code2(1,:)) ...
            +sum(data2LSB_detected_MLD~=qpsk_code2(2,:));
        bit_error2_P=bit_error2_P+sum(data2MSB_detected_P~=qpsk_code2(1,:)) ...
            +sum(data2LSB_detected_P~=qpsk_code2(2,:));
        i_block=i_block+1;
    end
    
    %BER calculation
    BER1_MLD(i_snr)=bit_error1/Data_length/Iter;
    BER2_MLD(i_snr)=bit_error2/Data_length/Iter;
    BER1_P(i_snr)=bit_error1_P/Data_length/Iter;
    BER2_P(i_snr)=bit_error2_P/Data_length/Iter;
    BER_su1(i_snr) = bit_error_su1/Data_length/Iter;
    BER_su2(i_snr) = bit_error_su2/Data_length/Iter;
    
    SER_1_PANOMA(i_snr) = ser_1_PANOMA/Iter/Frame_length;
    SER12_1_PANOMA(i_snr) = ser12_1_PANOMA/Iter/Frame_length;
    SER13_1_PANOMA(i_snr) = ser13_1_PANOMA/Iter/Frame_length;
    SER14_1_PANOMA(i_snr) = ser14_1_PANOMA/Iter/Frame_length;
    SER21_1_PANOMA(i_snr) = ser21_1_PANOMA/Iter/Frame_length;
    SER23_1_PANOMA(i_snr) = ser23_1_PANOMA/Iter/Frame_length;
    SER24_1_PANOMA(i_snr) = ser24_1_PANOMA/Iter/Frame_length;
    SER31_1_PANOMA(i_snr) = ser31_1_PANOMA/Iter/Frame_length;
    SER32_1_PANOMA(i_snr) = ser32_1_PANOMA/Iter/Frame_length;
    SER34_1_PANOMA(i_snr) = ser34_1_PANOMA/Iter/Frame_length;
    SER41_1_PANOMA(i_snr) = ser41_1_PANOMA/Iter/Frame_length;
    SER42_1_PANOMA(i_snr) = ser42_1_PANOMA/Iter/Frame_length;
    SER43_1_PANOMA(i_snr) = ser43_1_PANOMA/Iter/Frame_length;
    Correct1_1_PANOMA(i_snr) = correct1_1_PANOMA/Iter/Frame_length;
    Correct2_1_PANOMA(i_snr) = correct2_1_PANOMA/Iter/Frame_length;
    Correct3_1_PANOMA(i_snr) = correct3_1_PANOMA/Iter/Frame_length;
    Correct4_1_PANOMA(i_snr) = correct4_1_PANOMA/Iter/Frame_length;
    
    SER_2_PANOMA(i_snr) = ser_2_PANOMA/Iter/Frame_length;
    SER12_2_PANOMA(i_snr) = ser12_2_PANOMA/Iter/Frame_length;
    SER13_2_PANOMA(i_snr) = ser13_2_PANOMA/Iter/Frame_length;
    SER14_2_PANOMA(i_snr) = ser14_2_PANOMA/Iter/Frame_length;
    SER21_2_PANOMA(i_snr) = ser21_2_PANOMA/Iter/Frame_length;
    SER23_2_PANOMA(i_snr) = ser23_2_PANOMA/Iter/Frame_length;
    SER24_2_PANOMA(i_snr) = ser24_2_PANOMA/Iter/Frame_length;
    SER31_2_PANOMA(i_snr) = ser31_2_PANOMA/Iter/Frame_length;
    SER32_2_PANOMA(i_snr) = ser32_2_PANOMA/Iter/Frame_length;
    SER34_2_PANOMA(i_snr) = ser34_2_PANOMA/Iter/Frame_length;
    SER41_2_PANOMA(i_snr) = ser41_2_PANOMA/Iter/Frame_length;
    SER42_2_PANOMA(i_snr) = ser42_2_PANOMA/Iter/Frame_length;
    SER43_2_PANOMA(i_snr) = ser43_2_PANOMA/Iter/Frame_length;
    Correct1_2_PANOMA(i_snr) = correct1_2_PANOMA/Iter/Frame_length;
    Correct2_2_PANOMA(i_snr) = correct2_2_PANOMA/Iter/Frame_length;
    Correct3_2_PANOMA(i_snr) = correct3_2_PANOMA/Iter/Frame_length;
    Correct4_2_PANOMA(i_snr) = correct4_2_PANOMA/Iter/Frame_length;
    toc
end
% table(SNR_dB',BER1_MLD',BER1_P',BER2_MLD',BER2_P')
e1_1_PANOMA = SER12_1_PANOMA + SER21_1_PANOMA + SER31_1_PANOMA + SER42_1_PANOMA;
e2_1_PANOMA = SER13_1_PANOMA + SER24_1_PANOMA + SER34_1_PANOMA + SER43_1_PANOMA;
e3_1_PANOMA = SER14_1_PANOMA + SER23_1_PANOMA + SER32_1_PANOMA + SER41_1_PANOMA;
correct_1_PANOMA = Correct1_1_PANOMA + Correct2_1_PANOMA + Correct3_1_PANOMA + Correct4_1_PANOMA;
Capacity_1_PANOMA = -e1_1_PANOMA.*log2(e1_1_PANOMA) - e2_1_PANOMA.*log2(e2_1_PANOMA) - e3_1_PANOMA.*log2(e3_1_PANOMA) - (correct_1_PANOMA.*log2(correct_1_PANOMA));
Capacity_1_PANOMA = 2 - Capacity_1_PANOMA;
e1_2_PANOMA = SER12_2_PANOMA + SER21_2_PANOMA + SER31_2_PANOMA + SER42_2_PANOMA;
e2_2_PANOMA = SER13_2_PANOMA + SER24_2_PANOMA + SER34_2_PANOMA + SER43_2_PANOMA;
e3_2_PANOMA = SER14_2_PANOMA + SER23_2_PANOMA + SER32_2_PANOMA + SER41_2_PANOMA;
correct_2_PANOMA = Correct1_2_PANOMA + Correct2_2_PANOMA + Correct3_2_PANOMA + Correct4_2_PANOMA;
Capacity_2_PANOMA = -e1_2_PANOMA.*log2(e1_2_PANOMA) - e2_2_PANOMA.*log2(e2_2_PANOMA) - e3_2_PANOMA.*log2(e3_2_PANOMA) - (correct_2_PANOMA.*log2(correct_2_PANOMA));
Capacity_2_PANOMA = 2 - Capacity_2_PANOMA;
sumRate = Capacity_1_PANOMA + Capacity_2_PANOMA;

avgBER = 0.5*(BER1_MLD+BER2_MLD);
avgBER_P = 0.5*(BER1_P+BER2_P);

Sd1 = sqrt(PL(1)./(4*N*SNR));
Sd2 = sqrt(PL(2)./(4*N*SNR));
Sd_su2 = sqrt(PL(2)./(SNR));
Sd_su1 = sqrt(PL(1)./(SNR));

gamma_su1 = (1./Sd_su1).^2;
gamma_su2 = (1./Sd_su2).^2;
BERth_su1 = 0.5*(1 - sqrt(0.5*gamma_su1./(1 + 0.5*gamma_su1)));
BERth_su2 = 0.5*(1 - sqrt(0.5*gamma_su2./(1 + 0.5*gamma_su2)));
gamma1_u1 = (sqrt(a(1))./(sqrt(2)*Sd1)).^2;
gamma2_u1 = ((sqrt(a(1))+sqrt(a(2)))./(sqrt(2)*Sd1)).^2;
gamma3_u1 = ((sqrt(a(1))+2*sqrt(a(2)))./(sqrt(2)*Sd1)).^2;
gamma4_u1 = ((-sqrt(a(1))+sqrt(a(2)))./(sqrt(2)*Sd1)).^2;
gamma5_u1 = ((-sqrt(a(1))+2*sqrt(a(2)))./(sqrt(2)*Sd1)).^2;

gamma1_u2 = ((sqrt(a(1))+sqrt(a(2)))./(sqrt(2)*Sd2)).^2;
gamma2_u2 = ((-sqrt(a(1))+sqrt(a(2)))./(sqrt(2)*Sd2)).^2;
BER2th = 0.25*(1 - sqrt(0.5*gamma1_u2./(1 + 0.5*gamma1_u2))) ...
    + 0.25*(1 - sqrt(0.5*gamma2_u2./(1 + 0.5*gamma2_u2)));
BER1th = 0.5*(1 - sqrt(0.5*gamma1_u1./(1 + 0.5*gamma1_u1)))  ...
    - 0.25*(1 - sqrt(0.5*gamma2_u1./(1 + 0.5*gamma2_u1))) ...
    + 0.25*(1 - sqrt(0.5*gamma3_u1./(1 + 0.5*gamma3_u1))) ...
    + 0.25*(1 - sqrt(0.5*gamma4_u1./(1 + 0.5*gamma4_u1))) ...
    - 0.25*(1 - sqrt(0.5*gamma5_u1./(1 + 0.5*gamma5_u1)));
avgBERth = 0.5*(BER1th+BER2th);
A11 = real(MLD_LUT_P(1));
Ag11 = real(MLD_LUT_P(5));
xgamma1_u1 = ((A11 - Ag11)./(2*Sd1)).^2;
xgamma2_u1 = ((3*A11+Ag11)./(2*Sd1)).^2;
xgamma3_u1 = ((Ag11)./(Sd1)).^2;
xgamma4_u1 = ((A11)./(Sd1)).^2;
xgamma5_u1 = ((3*Ag11+A11)./(2*Sd1)).^2;
xgamma1_u2 = ((A11)./(Sd2)).^2;
xgamma2_u2 = ((Ag11)./(Sd2)).^2;

BER1th_PANOMA = 0.5*(1 - sqrt(0.5*xgamma1_u1./(1 + 0.5*xgamma1_u1))) ...
               +1/2*0.5*(1 - sqrt(0.5*xgamma2_u1./(1 + 0.5*xgamma2_u1))) ...
               +1/2*0.5*(1 - sqrt(0.5*xgamma3_u1./(1 + 0.5*xgamma3_u1))) ...
               -1/2*0.5*(1 - sqrt(0.5*xgamma4_u1./(1 + 0.5*xgamma4_u1))) ...
               -1/2*0.5*(1 - sqrt(0.5*xgamma5_u1./(1 + 0.5*xgamma5_u1)));
BER2th_PANOMA = 1/2*0.5*(1 - sqrt(0.5*xgamma1_u2./(1 + 0.5*xgamma1_u2))) ...
               +1/2*0.5*(1 - sqrt(0.5*xgamma2_u2./(1 + 0.5*xgamma2_u2)));
avgBERth_PANOMA = 0.5*(BER1th_PANOMA+BER2th_PANOMA);
fprintf('This program finished at %s\n', datestr(now,'HH:MM:SS'))
%% Plotting
figure
semilogy(0,0,'--k',0,0,'-k',0,0,'-.k',0,0,':k',0,0,'ok')
hold on
semilogy(SNR_dB,BER1th,'--k',SNR_dB,BER1_MLD,'ok',SNR_dB,BER2th,'-k',SNR_dB,BER2_MLD,'ok',...
    SNR_dB,BERthBPSK_su,'-.k',SNR_dB,BER_su1,'ok',SNR_dB,BERthQPSK_su,':k',SNR_dB,BER_su2,'ok')
set(gca,'FontSize',20,'LineWidth',1')
hold off
xlabel('SNR [dB]');
ylabel('BER');
legend('U_{1}','U_{2}','SU | U_{1}','SU | U_{2}','Simulations')
axis([min(SNR_dB) max(SNR_dB) 1e-4 1])
grid on