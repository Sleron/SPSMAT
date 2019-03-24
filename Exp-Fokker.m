
%nonlinear fokker planck equation with GLGGL collocation method
function Fokker_nonlinear

alpha_=0.5;


N=7
teta=0.5
delta_t=0.001
Maxstep=1/delta_t
 
output_precision(20);
 


 x=Gnbr_zeros(N-1,alpha_+1)';
 x=[0,(x+1)/2,1];

 A1=eye(1,N+1);% for the Bounday Condition y(0,t)=0

 A2=eye(N+1,N+1);% I
 A2=A2(2:N,:);
 A3=[zeros(1,N),1];%for the Bounday Condition  y(1,t)=exp(t)

 D_1=D_lgr_Gnbr(N+1,alpha_,x', @(x) 2*x-1,'gauss_lbt');
 D_2=(0+diag(ones(N+1,1)*2)*D_1)*(diag(ones(N+1,1)*2)^(-1))*D_1;

 A4=D_2(2:N,:); 
 A5=D_1(2:N,:);


 C(1:(N+1),1)=x'.^2; % Initial condition
 
 for step=2:Maxstep
  
    v1=(delta_t)*(ones(1,N-1)'.*(1/3)+...
       C(2:N,step-1).*(4./(x(2:N)'.^2))-...
       8./(x(2:N)').*A5*C(1:N+1,step-1)+...
       (2*A4*C(1:N+1,step-1)));
 
    v2=(delta_t)*(x(2:N)'/3+2*A5*C(1:N+1,step-1));
  
    m1=ones(1,N-1)'-teta*v1;
    m2=teta*v2;
  
    t=delta_t*step;

 
    m_hat1=ones(1,N-1)'+(1-teta)*v1;
    m_hat2=(1-teta)*v2;

    A= [A1;  diag(m1)*A2-diag(m2)*(A5);A3];
    b(:,step-1)=[0;  (diag(m_hat1)*A2+diag(m_hat2)*(A5))*C(:,step-1);exp(t)];

    C(:,step) = A \ b(:,step-1);

 end %for step

%------------------plotting----------

 
max_poi=25; 
poi=linspace(0,1,max_poi);


for ti=1:1:Maxstep/40
for xi=1:1:max_poi
  y(xi,ti)=lgr_Gnbr_(N+1,poi(xi),x, @(x) 2*x-1)*C(1:(N+1),ti*40)...
  -(poi(xi)^2*exp((ti*40)*delta_t));
  end
end

[px py]=ndgrid(poi,poi);

figure(1);
surf(py,px,abs(y));

end %function

