function [E,S_out] = find_equivalences(enscells_in_1, enscells_in_2, metric)
% Gets equivalences between ensembles detected and ensembles generated.
% INPUT:
%   enscells_in_1: <logical> N-by-nens_1 matrix with the core cells of each ensemble generated.
%   enscells_in_2: <logical> N-by-nens_2 matrix with the core cells of each ensemble detected.
%   metric: <string> with the metric used to find equivalences.
% OUTPUT:
%   E: <integer> max(nens_1,nens_2)-by-2 matrix with equivalences between ensembles
%     detected and ensembles generated. A row with values (i,j) means that i-th ensemble
%     generated is equivalent to j-th ensemble detected.
%   S_out: <float> nens_1-by-nens_2 matrix with the similarity values between ensembles
%     detected and ensembles generated. The value in component (i,j)-th represents the
%     similarity value between i-th ensemble generated and j-th ensemble detected.

% Get number of ensembles
nens_1 = size(enscells_in_1,2);
nens_2 = size(enscells_in_2,2);

% Pre-allocate memory
E = zeros(max(nens_1,nens_2),2); % Equivalence matrix
S = zeros(nens_1, nens_2); %Similarity matrix

% Get similarity values
for i = 1:nens_1 %Rows are ensembles 1
    for j = 1:nens_2 %Columns are ensembles 2
        S(i,j) = 1-pdist2(double(enscells_in_1(:,i))',double(enscells_in_2(:,j))',metric);
    end
end

S_out = S; %Save similarity matrix because it will change

% Get equivalences
for i = 1:size(E,1)
    max_value = max(S(:));
    if max_value < 0
        diff = abs(nens_1-nens_2);
        switch size(E,1)
            case nens_1
                ens_not_found = find(~ismember((1:nens_1)',E(:,1)));
                E(end-diff+1:end,1) = ens_not_found;
            case nens_2
                ens_not_found = find(~ismember((1:nens_2)',E(:,2)));
                E(end-diff+1:end,2) = ens_not_found;    
        end
        break;
    end
    [row_max, col_max] = find(S == max_value);
    if(length(row_max)>1)
        row_max = row_max(1);
        col_max = col_max(1);
    end
    S(row_max,:) = -1;
    S(:,col_max) = -1;
    E(i,:) = [row_max, col_max];
end  

end