function f = f_fFf(F)


f11 = F(1,1);
f12 = F(1,2);
f13 = F(1,3);
f21 = F(2,1);
f22 = F(2,2);
f23 = F(2,3);
f31 = F(3,1);
f32 = F(3,2);
f33 = F(3,3);

coeff2 = f11^2*f33-2*f11*f13*f31+f12^2*f33-2*f12*f13*f32+f21^2*f33-2*f21*f23*f31+f22^2*f33-2*f22*f23*f32; 
coeff1 = -f13^2*f33-f23^2*f33-f31^2*f33-f32^2*f33;
coeff0 = -f33^3
                     
                     
mu = roots([coeff2, coeff1, coeff0]);

est_f(1) = sqrt(mu(1));
est_f(2) = sqrt(mu(2));

if abs(imag(est_f(1))) < 1e-8
    if abs(imag(est_f(2))) < 1e-8
        f = est_f;
    else
        f = est_f(1);
    end
else
    f = est_f(2);
end
        
end