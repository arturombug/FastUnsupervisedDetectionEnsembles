function [tpr,fpr] = tpr_fpr_pks(t_pks,e_pks)
% [tpr,fpr] = tpr_fpr_pks(t_pks,e_pks)
% Computes the true posive and false positive rate using the correct peak
% times in t_pks and the estimated peak times e_pks


%tp = sum(ismember(e_pks,t_pks));
tp = sum(eq(e_pks,t_pks));
%fp = sum(~ismember(e_pks,t_pks));
fp = sum(~eq(e_pks,t_pks));
tpr = tp/length(t_pks);
fpr = fp/length(e_pks);