function map = idv_to_mid(n)

map = zeros(n + 1, 2 ^ n);
  for r = 1:size(map, 2)
     map(sum(dec2bin(r - 1) == '1') + 1, r) = 1;
  end

end