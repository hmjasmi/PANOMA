function [BER_avg,BERth_NU,BERth_MU,BERth_FU,sumRate] = CNOMA3UEavgBER(a,PL,SNR_dB)
N = 3;
Va1 = [0.15/2 a(1) 0.050 1/3];    %NU power coefficients vector
Va2 = [0.15/2 a(2) 0.475 1/3];    %MU power coefficients vector
Va3 = [0.8500 a(3) 0.475 1/3];    %FU power coefficients vector
SVa1 = sqrt(Va1);
SVa2 = sqrt(Va2);
SVa3 = sqrt(Va3);
SNR = 10.^(SNR_dB./10);

%Lookup table for C-NOMA
MLD_LUT_C = [(1*SVa1(2)+1*SVa2(2)+1*SVa3(2)) ...
             (1*SVa1(2)+1*SVa2(2)-1*SVa3(2)) ...
             (1*SVa1(2)-1*SVa2(2)+1*SVa3(2)) ...
             (1*SVa1(2)-1*SVa2(2)-1*SVa3(2)) ...
             (-1*SVa1(2)+1*SVa2(2)+1*SVa3(2)) ...
             (-1*SVa1(2)+1*SVa2(2)-1*SVa3(2)) ...
             (-1*SVa1(2)-1*SVa2(2)+1*SVa3(2)) ...
             (-1*SVa1(2)-1*SVa2(2)-1*SVa3(2))];

% Theoretical BER calculations
Sd1 = sqrt(PL(1)./(2*N*SNR));
Sd2 = sqrt(PL(2)./(2*N*SNR));
Sd3 = sqrt(PL(3)./(2*N*SNR));
A = MLD_LUT_C(1);
B = MLD_LUT_C(5);
C = MLD_LUT_C(3);
D = MLD_LUT_C(7);

gamma1_u1 = (((A - B)/2)./(Sd1)).^2;
gamma2_u1 = (((B - C)/2)./(Sd1)).^2;
gamma3_u1 = ((D)./(Sd1)).^2;
gamma1c_u1 = (((2*A - C - D)/2)./(Sd1)).^2;
gamma2c_u1 = (((2*A - B - C)/2)./(Sd1)).^2;
gamma3c_u1 = (((2*A + C + D)/2)./(Sd1)).^2;
gamma4c_u1 = ((A)./(Sd1)).^2;
gamma5c_u1 = (((3*A + B)/2)./(Sd1)).^2;
gamma6c_u1 = (((2*A + B + C)/2)./(Sd1)).^2;
gamma7c_u1 = ((B)./(Sd1)).^2;
gamma8c_u1 = (((2*B - C - D)/2)./(Sd1)).^2;
gamma9c_u1 = (((3*B + C)/2)./(Sd1)).^2;
gamma10c_u1 = (((2*B + C + D)/2)./(Sd1)).^2;
gamma11c_u1 = (((3*B + A)/2)./(Sd1)).^2;
gamma12c_u1 = (((3*C + D)/2)./(Sd1)).^2;
gamma13c_u1 = ((C)./(Sd1)).^2;
gamma14c_u1 = (((3*C + B)/2)./(Sd1)).^2;
gamma15c_u1 = (((2*D + B + C)/2)./(Sd1)).^2;
gamma16c_u1 = (((3*D + C)/2)./(Sd1)).^2;

gamma1_u2 = (((2*A - B - C)/2)./(Sd2)).^2;
gamma2_u2 = (((B - C)/2)./(Sd2)).^2;
gamma3_u2 = ((C)./(Sd2)).^2;
gamma4_u2 = ((D)./(Sd2)).^2;
gamma1c_u2 = (((2*A + B + C)/2)./(Sd2)).^2;
gamma2c_u2 = ((A)./(Sd2)).^2;
gamma3c_u2 = (((3*B + C)/2)./(Sd2)).^2;
gamma4c_u2 = ((B)./(Sd2)).^2;
gamma5c_u2 = (((3*C + B)/2)./(Sd2)).^2;
gamma6c_u2 = (((2*D + B + C)/2)./(Sd2)).^2;

gamma1_u3 = ((A)./Sd3).^2;
gamma2_u3 = ((B)./Sd3).^2;
gamma3_u3 = ((C)./Sd3).^2;
gamma4_u3 = ((D)./Sd3).^2;

BERth_NU = 0.5*(1 - sqrt(0.5*gamma1_u1./(1 + 0.5*gamma1_u1))) + 0.25*(1 - sqrt(0.5*gamma2_u1./(1 + 0.5*gamma2_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma3_u1./(1 + 0.5*gamma3_u1))) ...
    + 0.25*(1 - sqrt(0.5*gamma1c_u1./(1 + 0.5*gamma1c_u1))) - 0.25*(1 - sqrt(0.5*gamma2c_u1./(1 + 0.5*gamma2c_u1))) ...
    + 0.25*(1 - sqrt(0.5*gamma3c_u1./(1 + 0.5*gamma3c_u1))) - 0.125*(1 - sqrt(0.5*gamma4c_u1./(1 + 0.5*gamma4c_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma5c_u1./(1 + 0.5*gamma5c_u1))) - 0.125*(1 - sqrt(0.5*gamma6c_u1./(1 + 0.5*gamma6c_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma7c_u1./(1 + 0.5*gamma7c_u1))) - 0.25*(1 - sqrt(0.5*gamma8c_u1./(1 + 0.5*gamma8c_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma9c_u1./(1 + 0.5*gamma9c_u1))) - 0.25*(1 - sqrt(0.5*gamma10c_u1./(1 + 0.5*gamma10c_u1))) ...
    - 0.125*(1 - sqrt(0.5*gamma11c_u1./(1 + 0.5*gamma11c_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma12c_u1./(1 + 0.5*gamma12c_u1))) - 0.125*(1 - sqrt(0.5*gamma13c_u1./(1 + 0.5*gamma13c_u1))) ...
    - 0.125*(1 - sqrt(0.5*gamma14c_u1./(1 + 0.5*gamma14c_u1))) ...
    + 0.125*(1 - sqrt(0.5*gamma15c_u1./(1 + 0.5*gamma15c_u1))) - 0.125*(1 - sqrt(0.5*gamma16c_u1./(1 + 0.5*gamma16c_u1)));
BERth_MU = 0.25*(1 - sqrt(0.5*gamma1_u2./(1 + 0.5*gamma1_u2))) ...
    + 0.25*(1 - sqrt(0.5*gamma2_u2./(1 + 0.5*gamma2_u2))) + 0.125*(1 - sqrt(0.5*gamma3_u2./(1 + 0.5*gamma3_u2))) ...
    + 0.125*(1 - sqrt(0.5*gamma4_u2./(1 + 0.5*gamma4_u2))) ...
    + 0.125*(1 - sqrt(0.5*gamma1c_u2./(1 + 0.5*gamma1c_u2))) - 0.125*(1 - sqrt(0.5*gamma2c_u2./(1 + 0.5*gamma2c_u2))) ...
     + 0.125*(1 - sqrt(0.5*gamma3c_u2./(1 + 0.5*gamma3c_u2))) - 0.125*(1 - sqrt(0.5*gamma4c_u2./(1 + 0.5*gamma4c_u2)))...
     - 0.125*(1 - sqrt(0.5*gamma5c_u2./(1 + 0.5*gamma5c_u2)))...
     - 0.125*(1 - sqrt(0.5*gamma6c_u2./(1 + 0.5*gamma6c_u2)));
BERth_FU = 0.125*(1 - sqrt(0.5*gamma1_u3./(1 + 0.5*gamma1_u3))) ...
    + 0.125*(1 - sqrt(0.5*gamma2_u3./(1 + 0.5*gamma2_u3))) + 0.125*(1 - sqrt(0.5*gamma3_u3./(1 + 0.5*gamma3_u3))) ...
    + 0.125*(1 - sqrt(0.5*gamma4_u3./(1 + 0.5*gamma4_u3))); 

BER_avg = 1/3*(BERth_NU+BERth_MU+BERth_FU);
c1=-((1-BERth_NU).*log2(1-BERth_NU)) ...
        - ((BERth_NU).*log2(BERth_NU));
c1 = 1 - c1;
c2=-((1-BERth_MU).*log2(1-BERth_MU)) ...
        - ((BERth_MU).*log2(BERth_MU));
c2 = 1 - c2;
c3=-((1-BERth_FU).*log2(1-BERth_FU)) ...
        - ((BERth_FU).*log2(BERth_FU));
c3 = 1 - c3;
sumRate = c1+c2+c3;