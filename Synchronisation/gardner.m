function out_eps = gardner(kappa,symb_rx,upsampling,T)
% OUTPUT out_eps=output sampled at sampling time taking eps_tild into
% account

eps_tild=0;
n=0;
while n<=length(symb_tx)
    n=n+upsampling+eps_tild*upsampling;
    eps_tild=eps_tild+2*kappa*real()