function [BER_avg,BERth_NU,BERth_FU,sumRate] = CNOMA2UEavgBER(a,PL,SNR_dB)
N = 2;
Va1 = [a(1) 0.5];               %NU power coefficients vector
Va2 = [a(2) 0.5];               %FU power coefficients vector
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
SNR = 10.^(SNR_dB./10);

%Lookup table for C-NOMA
MLD_LUT_C = [(1*SVa1(1)+1*SVa2(1)) ...
    (1*SVa1(1)-1*SVa2(1)) ...
    (-1*SVa1(1)+1*SVa2(1)) ...
    (-1*SVa1(1)-1*SVa2(1))];
% Theoretical BER calculations
A11 = MLD_LUT_C(1);
A01 = MLD_LUT_C(3);
Sd1 = sqrt(PL(1)./(2*N*SNR));
Sd2 = sqrt(PL(2)./(2*N*SNR));
gamma1_u1 = (((A11 - A01)/2)./(Sd1)).^2;
gamma2_u1 = ((A01)./(Sd1)).^2;
gamma1c_u1 = ((A11)./(Sd1)).^2;
gamma2c_u1 = (((3*A11 + A01)/2)./(Sd1)).^2;
gamma3c_u1 = (((A11 + 3*A01)/2)./(Sd1)).^2;

gamma1_u2 = ((A11)./(Sd2)).^2;
gamma2_u2 = ((A01)./(Sd2)).^2;

BERth_NU = 0.5*(1 - sqrt(0.5*gamma1_u1./(1 + 0.5*gamma1_u1))) + 0.25*(1 - sqrt(0.5*gamma2_u1./(1 + 0.5*gamma2_u1)))...
    - 0.25*(1 - sqrt(0.5*gamma1c_u1./(1 + 0.5*gamma1c_u1))) - 0.25*(1 - sqrt(0.5*gamma3c_u1./(1 + 0.5*gamma3c_u1)))...
    + 0.25*(1 - sqrt(0.5*gamma2c_u1./(1 + 0.5*gamma2c_u1)));
BERth_FU = 0.25*(1 - sqrt(0.5*gamma1_u2./(1 + 0.5*gamma1_u2))) + 0.25*(1 - sqrt(0.5*gamma2_u2./(1 + 0.5*gamma2_u2)));
BER_avg = 0.5*(BERth_NU+BERth_FU);
c1=-((1-BERth_NU).*log2(1-BERth_NU)) ...
        - ((BERth_NU).*log2(BERth_NU));
c1 = 1 - c1;
c2=-((1-BERth_FU).*log2(1-BERth_FU)) ...
        - ((BERth_FU).*log2(BERth_FU));
c2 = 1 - c2;
sumRate = c1+c2;