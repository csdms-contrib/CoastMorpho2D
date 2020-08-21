function [ S, domg] = JONSWAP( Ohm, Hs, Tp)%
%Ohm is the omega  2*pi/T
%JONSWAP - Calculates the wave spectrum values for a JONSWAP spectrum

wp = 2*pi/Tp;
Gamma = 3.3;
 for x = 1:length(Ohm)
     if Ohm(x)<wp
         Sigma = 0.07;
     else
         Sigma = 0.09;
     end
     A = exp(-((Ohm(x)/wp-1)/(Sigma*sqrt(2)))^2);
     S(x) = 320*Hs^2*Ohm(x)^-5/Tp^4*exp(-1950*Ohm(x)^-4/Tp^4)*Gamma^A;
 end
 S=S*16;
 
 domg=NaN;
if  length(Ohm)>1
O1=(Ohm(1:end-1)+Ohm(2:end))/2;
domg=[O1(2:end)-O1(1:end-1)];domg=[domg(1) domg domg(end)];
end
%domg
 % Determine the frequency step from the frequency vector. Note that the
 % highest frequency step is extrapolated.
%  domg = zeros( size(Ohm) );
%  domg(1:end-1) = diff( Ohm );
%  domg(end) = domg(end-1);

% Ohm=Ohm/2/pi;
% O1=(Ohm(1:end-1)+Ohm(2:end))/2;
% domg=[O1(2:end)-O1(1:end-1)];domg=[domg(1) domg domg(end)];

 % Determine the amplitudes from the spectral values
 %Amp = sqrt( 2 * S .* domg );
 %Amp =  S .* domg ;
 
 %S=S.*domg/2/pi*100;
 %S=S.*domg/2/pi*100/Hs^2;

end
