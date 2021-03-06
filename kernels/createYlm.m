function Ylm = createYlm( lmax, colatitude, longitude )

Ylm = zeros( size(longitude,1), (lmax+2)*(lmax+1)/2 );
for l = 0:2:lmax
    Pm = legendre(l,cos(colatitude'))';
    lconstant = sqrt((2*l + 1)/(4*pi));

    center = (l+1)*(l+2)/2 - l;
    Ylm(:,center) = lconstant*Pm(:,1);
    for m=1:l
        precoeff = lconstant * sqrt(factorial(l - m)/factorial(l + m));

		 if mod(m,2) == 1
			 precoeff = -precoeff;
		 end
		Ylm(:, center + m) = sqrt(2)*precoeff*Pm(:,m+1).*cos(m*longitude);
		Ylm(:, center - m) = sqrt(2)*precoeff*Pm(:,m+1).*sin(m*longitude);
	end
end
