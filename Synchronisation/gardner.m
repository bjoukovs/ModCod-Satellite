function output = gardner(kappa,symb_rx,upsampling,n0,nsampling)
    % OUTPUT out_eps=output sampled at sampling time taking eps_tild into
    % account
    
    eps_tild=0;
    n=n0;
    
    eps_ev = [eps_tild];
    ns = [];
    
    output = zeros(1,length(symb_rx)/upsampling);
    
    i=0;
    
    while  n0 + (i+1)*upsampling + eps_tild*upsampling <= length(symb_rx) && i < length(symb_rx)/upsampling
        
        n = n0 + i*upsampling + eps_tild*upsampling;
        
        sample0 = getSample(symb_rx, n);
        sample1 = getSample(symb_rx, n + upsampling);
        sample_middle = getSample(symb_rx, n + upsampling/2);

        eps_tild = eps_tild + 2*kappa*real(sample_middle*(conj(sample1)-conj(sample0)));


        %n=n + upsampling + eps_tild*upsampling;
        
        ns = [ns n];
        eps_ev = [eps_ev, eps_tild];
        
        if n<0
            disp('error');
        end
        
        output(i+1) = sample0;
        i = i+1;
        
    end
    
    %figure;plot(ns);
    %hold on;
    %figure;plot(eps_ev);
    %hold off;
    
end



%Interpolation function
function s = getSample(samples,n)

   n_low = floor(n);
   n_high = ceil(n);
   
   %case if we fall just right on a sample
   if n_low == n_high
       s = samples(n_low);
   else
       %we have to interpolate linearly
       s_low = samples(n_low);
       s_high = samples(n_high);

       s = s_low + (s_high-s_low) * (n-n_low);
   end

end