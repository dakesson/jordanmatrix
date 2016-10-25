function [rest] = calculate_rest(x)
    rest = abs(round(x)-x);
end