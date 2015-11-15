clear;
w=1;
snr=0:0.5:4;
for z=1:numel(snr)
    c = 100;
    error=0;
    for j =2:c+12
            N = 100;
            ip = randint(1,N-2) ;
            ip(N)=0;
            d=poly2trellis(3,[5 7]);
            v=convenc(ip,d);
            for i=1:2*N
             s(i) = 2*v(i) - 1;
            end
            
            r = awgn(s,snr(w),'measured');
        
        x111=-1;x112=-1;x121=1;x122=1;
        x231=-1;x232=1;x241=1;x242=-1;
        x311=1;x312=1;x321=-1;x322=-1;
        x431=1;x432=-1;x441=-1;x442=1;
        %gamma calculation
         Eb=2;N0=Eb/(10^(snr(w)/10));sigma=N0/2;
        
        k=0;
        
       for i=1:2:2*N
            k=k+1;
            gama11(k)=exp(log(-(r(i)-x111)^2-log(r(i+1)-x112)^2)/2/sigma);
            gama12(k)=exp(log(-(r(i)-x121)^2-log(r(i+1)-x122)^2))/2/sigma;
            gama23(k)=exp(log(-(r(i)-x231)^2-log(r(i+1)-x232)^2)/2/sigma);
            gama24(k)=exp(log(-(r(i)-x241)^2-log(r(i+1)-x242)^2)/2/sigma);
            gama31(k)=exp(log(-(r(i)-x311)^2-log(r(i+1)-x312)^2)/2/sigma);
            gama32(k)=exp(log(-(r(i)-x321)^2-log(r(i+1)-x322)^2)/2/sigma);
            gama43(k)=exp(log(-(r(i)-x431)^2-log(r(i+1)-x432)^2)/2/sigma);
            gama44(k)=exp(log(-(r(i)-x441)^2-log(r(i+1)-x442)^2)/2/sigma);
            
            
           % gama11(k)=((-(r(i)-x111)^2-(r(i+1)-x112)^2)/2/sigma);
            %gama12(k)=((-(r(i)-x121)^2-(r(i+1)-x122)^2)/2/sigma);
            %gama23(k)=((-(r(i)-x231)^2-(r(i+1)-x232)^2)/2/sigma);
            %gama24(k)=((-(r(i)-x241)^2-(r(i+1)-x242)^2)/2/sigma);
            %gama31(k)=((-(r(i)-x311)^2-(r(i+1)-x312)^2)/2/sigma);
            %gama32(k)=((-(r(i)-x321)^2-(r(i+1)-x322)^2)/2/sigma);
            %gama43(k)=((-(r(i)-x431)^2-(r(i+1)-x432)^2)/2/sigma);
            %gama44(k)=((-(r(i)-x441)^2-(r(i+1)-x442)^2)/2/sigma);
            
            
       end
               
        %alfa calculation
        alfa1(1)=gama11(1);alfa2(1)=(gama12(1));alfa3(1)=1;alfa4(1)=1;%assigning initial values to alpha
        for i=2:N %starting the iteration
          % alfa1(i)=alfa1(i-1)*gama11(i)+alfa3(i-1)*gama31(i);
          % alfa2(i)=alfa1(i-1)*gama12(i)+alfa3(i-1)*gama32(i);
          % alfa3(i)=alfa2(i-1)*gama23(i)+alfa4(i-1)*gama43(i);
          % alfa4(i)=alfa2(i-1)*gama24(i)+alfa4(i-1)*gama44(i);
        
          
          alfa1(i)=max(alfa1(i-1)+gama11(i),alfa3(i-1)+gama31(i));
          alfa2(i)=max(alfa1(i-1)+gama12(i),alfa3(i-1)+gama32(i));
          alfa3(i)=max(alfa2(i-1)+gama23(i),alfa4(i-1)+gama43(i));
          alfa4(i)=max(alfa2(i-1)+gama24(i),alfa4(i-1)+gama44(i));
          
         
        end
        %Beta calculation
        beta1(N-1)=(gama11(N));beta3(N-1)=(gama31(N));beta2(N-1)=1;beta4(N-1)=1;%initially assigning the terminating values to beta
        for i=N-2:-1:1
           % beta1(i)=beta1(i+1)*(gama11(i+1))+beta2(i+1)*(gama12(i+1));
           % beta2(i)=beta3(i+1)*(gama23(i+1))+beta4(i+1)*(gama24(i+1));
           % beta3(i)=beta1(i+1)*(gama31(i+1))+beta2(i+1)*(gama32(i+1));
           % beta4(i)=beta3(i+1)*(gama43(i+1))+beta4(i+1)*(gama44(i+1));
            
            beta1(i)=max(beta1(i+1)+(gama11(i+1)),beta2(i+1)+(gama12(i+1)));
            beta2(i)=max(beta3(i+1)+(gama23(i+1)),beta4(i+1)+(gama24(i+1)));
            beta3(i)=max(beta1(i+1)+(gama31(i+1)),beta2(i+1)+(gama32(i+1)));
            beta4(i)=max(beta3(i+1)+(gama43(i+1)),beta4(i+1)+(gama44(i+1)));
            
        end
         %Lambda calculation
         
          Lambda=gama12(1)*beta2(1);
               L_a(1)=log((gama11(1)*beta1(1))/(Lambda)); % initializing log likelihood ratio
            for i=2:N-2
              p=log(exp(alfa1(i-1)*gama11(i)*beta1(i))+exp(alfa2(i-1)*gama23(i)*beta3(i))+...
                 exp(alfa3(i-1)*gama31(i)*beta1(i))+exp(alfa4(i-1)*gama43(i)*beta3(i)));
               Lambda=log(exp(alfa1(i-1)*gama12(i)*beta2(i))+exp(alfa2(i-1)*gama24(i)*beta4(i))+...
                  exp(alfa3(i-1)*gama32(i)*beta2(i))+exp(alfa4(i-1)*gama44(i)*beta4(i)));
               L_a(i)=(p - Lambda);
           end 
        %Lambda=(gama12(1))+beta2(1);
        %L_a(1)=(((gama11(1))+beta1(1)) - (Lambda));% initializing log likelihood ratio
        %for i=2:N-2
         %   p =log(alfa1(i-1)+(gama11(i))+beta1(i),alfa2(i-1)+(gama23(i))+beta3(i),...
          %       alfa3(i-1)+(gama31(i))+beta1(i),alfa4(i-1)+(gama43(i))+beta3(i));
           % Lambda = log(alfa1(i-1)+(gama12(i))+beta2(i),alfa2(i-1)+(gama24(i))+beta4(i),...
            %       alfa3(i-1)+(gama32(i))+beta2(i),alfa4(i-1)+(gama44(i))+beta4(i));
            %L_a(i)=(p - Lambda);
        %end
        %estimating the input matrix a
        for i=1:N-2
            if (L_a(i))>0
                a(i)=0;
            else
                a(i)=1;
            end
        end
        a(N)=0;
        %difference between the input given and the estimated input
        error =  error + sum(xor(ip,a));
       
    end
       ber(z) = error/(N*c);%overall variation in the difference 
       
       w=w+1;
       
end

semilogy(snr,ber),grid
title('MAP');xlabel('Eb/N0(dB)');ylabel('Bit Error Rate');

