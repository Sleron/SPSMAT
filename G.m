function out=G(n,a_,b_,u_)

f = @(x) u_;

% these next lines take the Anonymous function into a symbolic formula

pkg load symbolic
syms x;
u = f(x);

% now calculate the derivative of the function

g=diff(u, x);
g=function_handle(g);
%u_ == ax+b
a=g(0);
b=u(0);

%[A,B]
A=(-b-1)/a;
B=(1-b)/a;
G=zeros(n+1,n+1);

for i=1:1:n+1
 for j=1:1:n+1
  if i==j
    G(i,j)=(b_^2-a_^2)/((a_+b_+2*i+1)*(a_+b_+2*i+2));
   elseif i==j-1
    G(i,j)=2*(i+a_)*(i+b_)/((a_+b_+2*i)*(a_+b_+2*i+1));
   elseif i==j+1
    G(i,j)=2*(i)*(i+b_+a_)/((a_+b_+2*i)*(a_+b_+2*i-1));
  endif
 end%for
end%for


out=2*G/(B-A);
end%function
