function [I1,I2] = I1I2TG(U1,K1)                                  
                                                                  
%  COMPUTE THE I1(U1,K1) AND I2(U1,K1) FOR U1>0 AND U1<0 
%  We need this function for U1<0 as we get the answer 'NaN' when we
%  substitude in a negative value of u1 into the Desmarais function.
%  This function REQUIRES the Desmarais function                                 

%To understand this, see Max Blair's A COMPILATION OF THE MATHEMATICS LEADING
%TO THE DOUBLET-LATTICE METHOD written in 1994. See equation 282 and its'
%explanation. Essentially, the symmetry properties (evenness and oddness)
%of the left hand side of equation 274 can be used to formulate I1 for
%negative values of U1 as well. Similar logic applies for I2.
    if U1 >= 0.0                                      
        [I1,I2]=Desmarais(U1,K1);                                    
    
    else
        [J1,J2]= Desmarais(-U1,K1);                            
        [C1,C2]= Desmarais(0.0,K1);
        I1=2*real(C1)-real(J1)+1i*imag(J1); % equation 282
        I2=2*real(C2)-real(J2)+1i*imag(J2);
    end
   
end 
      

%Most modern method used