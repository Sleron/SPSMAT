function out=Cby_frac_4(N,derv,X,u_)
%
% 
% Overview
% This function returns a Shifted fractional Chebyshev (4th kind) matrix functions. 
%     
%out = Cby_frac_4(N,derv,X,u_)
%
%inputs: 
%----------------------------------------------------------------- 
%| N     : integer          : From fractional Chebyshev(4th kind)|
%|                            sentence 0 to fractional Chebyshev |
%|                            sentence number N                  |
%| derv  : double           : derivative order                   |
%|   X   : [1xm] double     : Inputs of u_(x) in                 | 
%|                            fractional Chebyshev (4th kind)    |
%|                            functions                          |
%|   u_  : symbolic function: Shifting parameter                 |          
%-----------------------------------------------------------------    
%
%Output:
%-----------------------------------------------------------------    
%| out   : [mx(N+1)] double : shifted fractional Chebyshev       |
%|                            (4th kind) functions               |       
%-----------------------------------------------------------------    
% 
alpha_=0.5;
beta_=-0.5;

multiplier=zeros(1,length(X));
for i=0:N
  multiplier(i+1)=2^(2*i)/nchoosek(2*i,i);
 end% for

 out=multiplier.*jacobi_frac(N,alpha_,beta_,derv,X,u_);




end