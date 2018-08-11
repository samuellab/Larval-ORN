function  y  = objfun( x, a, b, c, d )
% used for fit the non-linear part.
% a example is Or42a's data.

y = a ./ (1 + exp(b * (c-x))) + d;

end