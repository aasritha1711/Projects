clc;
clear;

for snr=0:4
    error(snr+1)=0;
    snr;
    c=[100 100 1000 1000 10000 100000];
    for j=1:c(snr+1)        
        n=100;
        u=randint(1,n-2);
        u(n)=0;
        trellis=poly2trellis(3,[5 7]);
        v=convenc(u,trellis);
        for i=1:2*n
            s(i)=2*v(i)-1;
        end
        r = awgn(s,snr,'measured');
        %BPSK Modulation
        x111=-1;x112=-1;x121=1;x122=1;
        x231=-1;x232=1;x241=1;x242=-1;
        x311=1;x312=1;x321=-1;x322=-1;
        x431=1;x432=-1;x441=-1;x442=1;
        %gamma calculation
        k=0;
       Eb=2;N0=Eb/(10^(snr/10));sigma=N0/2;
       for i=1:2:2*n
            k=k+1;
            gamma11(k)=exp((-(r(i)-x111)^2-(r(i+1)-x112)^2)/2/sigma);
            gamma12(k)=exp((-(r(i)-x121)^2-(r(i+1)-x122)^2)/2/sigma);
            gamma23(k)=exp((-(r(i)-x231)^2-(r(i+1)-x232)^2)/2/sigma);
            gamma24(k)=exp((-(r(i)-x241)^2-(r(i+1)-x242)^2)/2/sigma);
            gamma31(k)=exp((-(r(i)-x311)^2-(r(i+1)-x312)^2)/2/sigma);
            gamma32(k)=exp((-(r(i)-x321)^2-(r(i+1)-x322)^2)/2/sigma);
            gamma43(k)=exp((-(r(i)-x431)^2-(r(i+1)-x432)^2)/2/sigma);
            gamma44(k)=exp((-(r(i)-x441)^2-(r(i+1)-x442)^2)/2/sigma);
        end
        %alpha calculation:
        alpha1(1)=gamma11(1);alpha2(1)=gamma12(1);alpha3(1)=0;alpha4(1)=0;
        for i=2:n
            alpha1(i)=alpha1(i-1)*gamma11(i)+alpha3(i-1)*gamma31(i);
            alpha2(i)=alpha1(i-1)*gamma12(i)+alpha3(i-1)*gamma32(i);
            alpha3(i)=alpha2(i-1)*gamma23(i)+alpha4(i-1)*gamma43(i);
            alpha4(i)=alpha2(i-1)*gamma24(i)+alpha4(i-1)*gamma44(i);
        end
        %Beta calculation
        beta1(n-1)=gamma11(n);beta3(n-1)=gamma31(n);beta2(n-1)=0;beta4(n-1)=0;
        for i=n-2:-1:1
            beta1(i)=beta1(i+1)*gamma11(i+1)+beta2(i+1)*gamma12(i+1);
            beta2(i)=beta3(i+1)*gamma23(i+1)+beta4(i+1)*gamma24(i+1);
            beta3(i)=beta1(i+1)*gamma31(i+1)+beta2(i+1)*gamma32(i+1);
            beta4(i)=beta3(i+1)*gamma43(i+1)+beta4(i+1)*gamma44(i+1);
        end
        %L(ak) calculation
        lamda=gamma12(1)*beta2(1);
                L_a(1)=log((gamma11(1)*beta1(1))/(lamda));
        for i=2:n-2
            p=(alpha1(i-1)*gamma11(i)*beta1(i)+alpha2(i-1)*gamma23(i)*beta3(i)+...
                 alpha3(i-1)*gamma31(i)*beta1(i)+alpha4(i-1)*gamma43(i)*beta3(i));
            lamda=(alpha1(i-1)*gamma12(i)*beta2(i)+alpha2(i-1)*gamma24(i)*beta4(i)+...
                   alpha3(i-1)*gamma32(i)*beta2(i)+alpha4(i-1)*gamma44(i)*beta4(i));
            L_a(i)=log(p/lamda);
        end 
        for i=1:n-2
            if L_a(i)>0
                a(i)=0;
            else
                a(i)=1;
            end
        end
        a(n)=0;
        error(snr+1)=error(snr+1)+sum(abs(u-a));
    end    
    error(snr+1);
    bit_error_rate(snr+1)=error(snr+1)/100/c(snr+1)
end
snr=0:4;
semilogy(snr,bit_error_rate),grid
title('MAP');xlabel('Eb/N0(dB)');ylabel('Bit Error Rate');
