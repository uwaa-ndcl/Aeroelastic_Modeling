function [output]=Roger_function(k_lst,Q_lst,beta_bar_lst)

        %k_lst is the list of reduced frequencies.

        %Q_lst is just the aerodynamic forces in the complex domain
        %normalized with respect to dynamic pressure.
        
        %The output is [p0_bar,p1_bar,p2_bar,p3_bar,.......]
        %These values are scalars corresponding to the ith row, jth column
        %elements of P matrices.

        p0_bar=Q_lst(1);
        N=length(k_lst);
        nbeta_bar=length(beta_bar_lst);

        L=zeros(2*N,2+nbeta_bar);
        R=zeros(2*N,1);
        
        for l=1:N
            L(2*l-1,1:2)=[0 -1*(k_lst(l))^2];
            L(2*l,1:2)=[k_lst(l) 0 ];
            R(2*l-1,:) = real(Q_lst(l))-p0_bar;
            R(2*l,:) =  imag(Q_lst(l));
            
            if length(nbeta_bar)>0
                for root=1:nbeta_bar %Loop over lag roots
                    L(2*l-1,root+2)=k_lst(l)^2/(k_lst(l)^2+beta_bar_lst(root)^2);
                    L(2*l,root+2)=(k_lst(l)*beta_bar_lst(root))/(k_lst(l)^2+beta_bar_lst(root)^2);
                end
            end
        end
        x = (L'*L) \ (L'*R); %x is a column vector with [p1_bar, p2_bar, .....]'
 
        output=[p0_bar,x'];  
        
end