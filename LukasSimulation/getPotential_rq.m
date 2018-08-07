function val = getPotential_rq(x, r, V)

if x>=min(r) & x<=max(r)
    if ismember(x, r)
        val = V(x==r);
    else
        upper_bound = min(find(x<r));
        lower_bound = max(find(x>r));
        val = (V(upper_bound)+V(lower_bound))/2;
    end
elseif x > max(r)
    val = V(end);
else
    val = 0;
end