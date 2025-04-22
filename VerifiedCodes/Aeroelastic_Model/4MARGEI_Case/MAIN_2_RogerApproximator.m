%% Obtain necessary Data from User Input

%%%Read the rest of the data from json
filename = 'user_input.json';
fid = fopen(filename, 'r');
rawData = fread(fid, '*char')';  % Read the file as a string
fclose(fid);

% Decode the JSON string into a MATLAB structure
data = jsondecode(rawData);

n_modes=data.n_modes;
NC=data.n_controlmodes;
lag_roots_lst=data.lag_roots_lst;   % list of beta_bar s
gust_lag_roots_lst=data.gust_lag_roots_lst; %list of gust_beta_bar s
RedFreq_lst=data.RedFreq_lst;  %%defined as wb/U
DoIUseGust=data.DoIUseGust;


% Load the variables: Q_lst and Q_gust_lst if gusts are used. 
load('IntermediaryOutput1.mat');

% Add to path the folder in which the subfunctions are contained.
addpath(fullfile(pwd, '..', 'MAIN2_helperfunctions'));


%% %Perform RFA approximation of aerodynamic forces 

%If reduced frequency in the user input does not include 0, add it in.
if ~ismember(0, RedFreq_lst)
    RedFreq_lst = [0;RedFreq_lst]; 
end

RedFreq_lst=sort(RedFreq_lst); % Rearrange the reduced frequencies in ascending order.

nLag=length(lag_roots_lst);
Pbar=zeros(n_modes+NC,n_modes+NC,3+nLag);

for ii=1:n_modes+NC
    for jj=1:n_modes+NC
        Pbar(ii,jj,:)=Roger_function(RedFreq_lst,Q_lst(ii,jj,:),lag_roots_lst); 
    end
end
        
Pbarss=Pbar(1:n_modes,1:n_modes,:);
Pbarsc=Pbar(1:n_modes,n_modes+1:n_modes+NC,:);

save('IntermediaryOutput2.mat', 'Pbarss', 'Pbarsc');

%%%The gust part remains unverified.

if DoIUseGust=="Yes"
    nLagG=length(gust_lag_roots_lst);
    PbarG=zeros(n_modes+NC,1,3+nLagG);
        for ii=1:n_modes+NC
            PbarG(ii,1,:)=Roger_function(RedFreq_lst,Q_gust_lst(ii,1,:),gust_lag_roots_lst);
        end
    PbarGs=PbarG(1:n_modes,1,:);
    save('IntermediaryOutput2.mat', 'PbarGs', '-append');
end



%%%Verification of Pbar By testing the RFA against the Q_lst

k_fine_mesh=linspace(0,max(RedFreq_lst),1000);
beta1_bar=lag_roots_lst(1); %<--- here, I assumed that the userinput.json  has two lag roots.
beta2_bar=lag_roots_lst(2);
P0_bar=Pbar(:,:,1);
P1_bar=Pbar(:,:,2);
P2_bar=Pbar(:,:,3);
P3_bar=Pbar(:,:,4);
P4_bar=Pbar(:,:,5);

for i=1:n_modes+NC
    for j=1:n_modes+NC
        roger_approx = P0_bar(i,j) + 1i * k_fine_mesh .* P1_bar(i,j) - k_fine_mesh.^2 .* P2_bar(i,j) ...
        + (1i * k_fine_mesh) ./ (1i * k_fine_mesh + beta1_bar) .* P3_bar(i,j) ...
        + (1i * k_fine_mesh) ./ (1i * k_fine_mesh + beta2_bar) .* P4_bar(i,j);
        
        plottingQ=Q_lst(i,j,:);
        plottingQ=squeeze(plottingQ);
        
        if i==1 && j==1  %just to check one of the Roger fits and not to get a million graphs
            figure;
            hold on;
            plot(RedFreq_lst,real(plottingQ), 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
            plot(k_fine_mesh,real(roger_approx), 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5)
            title(sprintf('Real Part of Roger Approx. vs. Tabulated Data for element %d , %d', i, j));
            xlabel('Reduced Frequency');
            ylabel('Value of A');
            hold off;
    
            figure;
            hold on;
            plot(RedFreq_lst,imag(plottingQ), 'o', 'MarkerEdgeColor', 'r', 'MarkerSize', 10)
            plot(k_fine_mesh,imag(roger_approx), 'o', 'MarkerEdgeColor', 'b', 'MarkerSize', 5)
            title(sprintf('Imaginary Part of Roger Approx. vs. Tabulated Data for element %d , %d', i, j));
            xlabel('Reduced Frequency');
            ylabel('Value of A');
            hold off;
        end
    end
end


