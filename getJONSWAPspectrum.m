function [Tperiodi,Ejonswap]=getJONSWAPspectrum(Tp_swell,Ho,faci);
omegai=2*pi/Tp_swell*faci;Tperiodi=2*pi./omegai;
[Ejonswap, domg] = JONSWAP(omegai , Ho, Tp_swell);
for i=1:length(omegai); a=JONSWAP(omegai(i)+[-0.5:0.01:0.5]*domg(i), Ho, Tp_swell);Ejonswap(i)=mean(a);end
Ejonswap=Ejonswap.*domg;
Ejonswap=Ejonswap/sum(Ejonswap);