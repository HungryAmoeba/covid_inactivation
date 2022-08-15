function Et = TEFM(Pstress, q, m_r, Q, w0, r, f)
% Calculates theoretical threshold electric field 
% magnitude in accordance with equation 12 in 
% Yang 2015


f = f(:);

%f must be converted to angular frequency
f = f.*(2*pi);
m = size(f,1);


%make vectorized implementaiton based on 
% freq.
Et = zeros(m, 1);

Et = Pstress*pi*r^2*(m_r^2*(-f.^2+w0^2).^2+(w0*m_r/Q)^2*f.^2).^(1/2)/(3.45*q*m_r*w0^2);

Et = Et.*2*pi



end