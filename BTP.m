%Generation of Plane Wave
p=zeros (1000);
for x=1:1000,
    for y=1:1000,  
        if (x-501) ^2 + (y-501)^2 <= (30)^2,
        p(x,y) = 1;
        end
    end
end

%Generation of Vortex Wave
p1=zeros (1000);
for x=1:1000,
    for y=1:1000,
        phi = atan2((y-501), (x-501));
        if (x-501) ^2 + (y-501) ^2 <= (30)^2,
            p1(x,y) = exp(1i*1*phi);
        end
    end
end

%Creating Mash Grid
dx =1; %% in mm, SLM pixel
dy=dx;
lx=999; %% mm
ly=lx;
[m,n]=meshgrid(-0.5*lx:dx:0.5*lx,+0.5*ly:-dy:-0.5*ly);
[wmax, smax] =size(m);

%Generation of Scatter
Rho_m=0.040;
A_g2=exp (-(m.^2+n.^2)/Rho_m^2);
C_a2=(fft2(rand(wmax,smax)-0.5));
H2=(fft2(A_g2));
C_d2= (ifft2(C_a2. *real(H2)));
C2=exp(1i*20*pi*(real(C_d2))/max(real(C_d2(:))));
C2=C2(1:1000,1:1000);
H1=fftshift(fft2(p.*C2));
H2=fftshift(fft2(p1. *C2));

%Amplitude Distribution and Stokes Parameter generation
Ex=H1. /max(abs(H1(:)));
Ey=H2. /max(abs(H2(:)));
I_0=(Ex+Ey). *conj(Ex+Ey);
I_1=I_0(1:1000,1:1000);
I_d=0;
I_av=0;

l_u=0;
for n_r=-350:350,
    l_u=l_u+1;
    l_v=0;
    for m_r=-350:350,
        l_v=l_v+1;
        Ir_n1=I_1(501+n_r,501+m_r);
        I_av=I_av+Ir_n1;
        clear Ir_n1 Ir_n1;
    end
end

l_u=0;
for n_r=-350:350,
    l_u=l_u+1;
    l_v=0;
    for m_r=-350:350,
        l_v=l_v+1;
        I_r1=I_1(351+n_r:650+n_r,351+m_r:650+m_r);
        I_r2=I_1(501+n_r,501+m_r);
        I_cor=(I_r1-(I_av/200^2)).*(I_r2-(I_av/200^2));
        I_d=I_d+I_cor;
    end
end

clear I_r1 I_r2 I_cor Ir_n1 Ir_n2;
I_1=I_d./(I_av)^2;
I_1

