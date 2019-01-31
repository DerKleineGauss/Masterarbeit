function val = getPotential_rq(x, r, V)

if x>=min(r) & x<=max(r)
    if ismember(x, r)
        val = V(x==r);
    else      
        upper_bound = min(find(x<r));
        lower_bound = max(find(x>r));
%         val = (V(upper_bound)+V(lower_bound))/2;
        upper_weight = (r(upper_bound)-x)/(r(upper_bound)-r(lower_bound));
        lower_weight = (x-r(lower_bound))/(r(upper_bound)-r(lower_bound));
        val = upper_weight*V(upper_bound) + lower_weight*V(lower_bound);
    end
elseif x > max(r)
    val = V(end);
else
    val = 0;
end