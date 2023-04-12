function e = tanner(y, syndrome, m,H)

%Calculate error pattern, different name??
e = (sum(xor(syndrome,m).*H) + y) > sum(H)/2;
end