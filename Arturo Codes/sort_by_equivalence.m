function [ensmat_in_sort, enscells_in_sort, E, S] = sort_by_equivalence(ensmat_in_1,enscells_in_1,ensmat_in_2,enscells_in_2,metric)
% Sort ensembles detected by its equivalences with ensembles generated.
% INPUT:
%   enscells_in_1: <logical> N-by-nens_1 matrix with the core cells of each ensemble generated.
%   enscells_in_2: <logical> N-by-nens_2 matrix with the core cells of each ensemble detected.
%   ensmat_in: <logical> nens_1-by-T matrix with the activity of each ensemble generated.
%   ensmat_in: <logical> nens_2-by-T matrix with the activity of each ensemble detected.
%   metric: <string> with the metric used to find equivalences.
% OUTPUT:
%   ensmat_in_sort: <logical> nens-by-T-by-M matrix with the activity of each ensemble
%     detected sorted with respect to ensembles generated.
%   enscells_in_sort: <logical> N-by-nens_1 matrix with the core cells of each ensemble
%     detected sorted with respect to ensembles generated.
%   E: <integer> max(nens_1,nens_2)-by-2 matrix with equivalences between ensembles
%     detected and ensembles generated. A row with values (i,j) means that i-th ensemble
%     generated is equivalent to j-th ensemble detected.
%   S: <float> nens_1-by-nens_2 matrix with the similarity values between ensembles
%     detected and ensembles generated. The value in component (i,j)-th represents the
%     similarity value between i-th ensemble generated and j-th ensemble detected.

% Get number of neurons, frames and ensembles generated
T = size(ensmat_in_1,2);
N = size(enscells_in_1,1);
nens = size(enscells_in_1,2);

% Pre-allocate memory
enscells_in_sort = false(N,nens);
ensmat_in_sort = false(nens,T);

% Find equivalences
[E,S] = find_equivalences(enscells_in_1, enscells_in_2, metric);

% Sort by equivalences
for i = 1:size(E,1)
    if(~E(i,1) || ~E(i,2))
        break
    end
    enscells_in_sort(:,E(i,1)) = enscells_in_2(:,E(i,2));
    index = logical(ensmat_in_2(E(i,2),:));
    ensmat_in_sort(E(i,1),index) = true;
end

end