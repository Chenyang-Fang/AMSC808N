function ind_stat = stat_selection(M,y)
florida = find(y == -1);
indiana = find(y == 1);
PF = length(florida);
PI = length(indiana);

score_stat = abs(sum(M(florida,:))./PF - sum(M(indiana,:))./PI);
[~,ind_stat] = maxk(score_stat,10000);
ind_stat = sort(ind_stat);
end