%-----------------Nuclear Data from the Model-------------------------
global tau lambda_p beta_p beta_p PNL rho LAMBDA_p;


lambda_p =0.0787;
beta_p = 0.00755;
rho = 0.004;
LAMBDA_p=0.003;

%--------------------Initial conditions-------------------------------
global n_0 C_0;
n_0=1;
C_0=n_0*beta_p/(LAMBDA_p*lambda_p);
%---------------------------------------------------------------------
%----------------Calling the function -------------------------------
Iterative_step(100,0.01,rho,0.9,0.80,10)


%---------Definition of the time discretization-------------------------
%Section 6.3.2 of the paper.
%The variable "final_time" is the target time,
%the "step" is the lenght of the time in which was
%divided the interval. "order" is the fractional order of the derivative
%and "approx" is the number of terms used in equations (45) and (50)

function D = Iterative_step(final_time,step,rho,alpha_1,alpha_2,approx)
global vect_sol;
global lambda_p  beta_p; 
global n_0 C_0; 

vect_sol = [];
malla = ceil(final_time/step);
paso=step;
for i=0:malla
    i
    
    if i==0
        n_f = solution_neutrons(0,n_0,C_0,alpha_1,alpha_2,approx);
        c_f = solution_precursors(0,n_0,C_0,alpha_1,alpha_2,approx,i,rho);
    else
        n_f = solution_neutrons(paso,n_0,C_0,alpha_1,alpha_2,approx);
        c_f = solution_precursors(paso,n_0,C_0,alpha_1,alpha_2,approx,i,rho);
    end
    %Updating initial conditions 
    n_0=n_f;
    C_0=c_f;
    
    vect_sol = [vect_sol; i*step n_f c_f];
    size(vect_sol);
    i
end
vect_sol

filename = 'Densities_output_ramp_results_f_1.xlsx';
xlswrite(filename,vect_sol);

end



%-------------Analytical solution of the neutron density-------------
%Input = time, initial conditions, order, approax
%Output = neutron density
%In this case the parameter alpha_p is related to the fractional order.

function C = solution_neutrons(time,N0,C0,alpha_1,alpha_2,approx)
global tau lambda_p beta_p beta_p PNL rho LAMBDA_p;
%Definition of the constants a's and b's provided in Table2
a1=1;
a2=lambda_p^alpha_2; 
a3=(beta_p-rho)/(LAMBDA_p^alpha_1);
a4=-(lambda_p^alpha_2)*rho/(LAMBDA_p^alpha_1);
b1=N0;
b2=(N0+C0)*lambda_p^alpha_2;

sol = 0;
for m=0:approx
    factor1 = (-1)^m;
    factor2 = (a4/a1)^m;
    partial_sum=0;
    for k=0:m
      s1 = nchoosek(m,k);
      %s1 = factorial(m)/(factorial(k)*factorial(m-k))
      s2 = (time^((alpha_1+alpha_2)*m-alpha_2*k))/factorial(m);
      s3 = (a3/a4)^k;

      %for m==0, we have the standard Mittag-Leffler function
      if m==0
          M1=ml(-1*((lambda_p*time)^alpha_2),alpha_2,1+alpha_1*m-alpha_2*k);
          M2=ml(-1*((lambda_p*time)^alpha_2),alpha_2,alpha_2+alpha_1*m+1-alpha_2*k);
          
      %for m>=1, we have derivatives of the Mittag-Leffler function
      %and the Eq. 60) is used.
      else
          M1=factorial(m)*ml(-1*((lambda_p*time)^alpha_2),alpha_2,alpha_2*m+1+alpha_1*m-alpha_2*k);
          M2=factorial(m)*ml(-1*((lambda_p*time)^alpha_2),alpha_2,alpha_2*m+alpha_2+alpha_1*m+1-alpha_2*k);
          
      end
      SUMA_G=b1*M1+b2*(time^alpha_2)*M2;
      partial_sum=partial_sum+s1*s2*s3*SUMA_G;
    end
    sol=sol+factor1*factor2*partial_sum;
end
solution=sol;
C=solution;
end

%----------------Analytical solution of the precursors of delayed
%neutrons
function P = solution_precursors(time,N0,C0,alpha_1,alpha_2,approx,i,rho)
global lambda_p  beta_p rho LAMBDA_p;

a1=1;
a2=lambda_p^alpha_2; 
a3=(beta_p-rho)/(LAMBDA_p^alpha_1);
a4=-(lambda_p^alpha_2)*rho/(LAMBDA_p^alpha_1);
b1=N0;
b2=(N0+C0)*lambda_p^alpha_2;

c1=a1;
c2=a2+lambda_p^alpha_2*a1;
c3=a3;
c4=a2*lambda_p^alpha_2;
c5=a4+a3*lambda_p^alpha_2;
c6=a4*lambda_p^alpha_2;
d1=b1;
d2=b2;
sol = 0;
i
for m=0:approx
    factor1 = (-1)^m;
    D = [ ]; 
    partial_sum=0;
    D = particiones_c(m);
    
    for u=1:size(D,1)
      s1 = 1/(factorial(D(u,1))*factorial(D(u,2))*factorial(D(u,3))*factorial(D(u,4)));
      s2 = ((c6)^(D(u,1)))*((c5)^(D(u,2)))*((c4)^(D(u,3)))*((c3)^(D(u,4)));
      s_k = (alpha_2+alpha_1)*D(u,1)+alpha_1*D(u,2)+alpha_2*D(u,3)+(alpha_1-alpha_2)*D(u,4);
      s3 = (time^(alpha_2*m+alpha_1+s_k));
      if m==0
          M1=ml(-c2*(time^alpha_2),alpha_2,alpha_1+1+s_k);
          M2=ml(-c2*(time^alpha_2),alpha_2,alpha_1+alpha_2+1+s_k);
          
      else
          M1=factorial(m)*ml(-c2*(time^alpha_2),alpha_2,alpha_2*m+alpha_1+1+s_k,m+1);
          M2=factorial(m)*ml(-c2*(time^alpha_2),alpha_2,alpha_2*m+alpha_1+alpha_2+1+s_k,m+1);
         
      end
      SUMA_G=d1*M1+d2*(time^alpha_2)*M2;
      partial_sum=partial_sum+s1*s2*s3*SUMA_G;
    end
    sol=sol+factor1*partial_sum;
end
l=(beta_p/LAMBDA_p)*sol;
l1=C0*ml(-(lambda_p^alpha_2)*time^alpha_2,alpha_2,1);
solution=(beta_p/LAMBDA_p)*sol/c1+C0*ml(-(lambda_p^alpha_2)*time^alpha_2,alpha_2,1);
P=solution
end
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

%---------------------Kronecker's delta defined in Eq. (56)-----------
function A = deltakronecker(n,m)
if n==m
    A=1;
else
    A=0;
end
end
%------------------------------------------------------------------------


%------ Partitions of integers for the Precursors of delayed neutrons-----
%Input: a natural number, n
%Output: all the partitios of the number n, considering 4 elements, i.e.
%x_1+x_2+x_3+x_4=n, The partitions are stored as arrays 
%[x1,x_2,x_3,x_4]
function B = particiones_c(n)
L=[];

for k_0=0:n
    for k_1=0:n
        for k_2=0:n
            for k_3=0:n          
                if deltakronecker(k_0+k_1+k_2+k_3,n)==1
                        L = [L;k_0 k_1 k_2 k_3];
                end
            end
        end
    end
end
B =L;
end
%------------------------------------------------------------------------
