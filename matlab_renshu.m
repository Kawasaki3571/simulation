 ix = 100;
 jx = 100;
 % p = zeros(ix + 1, jx + 1);
 % for i = 1 : ix + 1
 %     for j = 1 : jx + 1
 %         p(i,j) = sin(3.14 * 0.01 * (i^2 + j^2)^0.5);
 %     end
 % end
 dec_p(ix,jx);
 p = dec_p(ix,jx);
 disp(p);
 function p = dec_p(ix,jx)
 for i = 1 : ix + 1
     for j = 1 : jx + 1
         p(i,j) = sin(3.14 * 0.01 * (i^2 + j^2)^0.5);
     end
 end
 end

