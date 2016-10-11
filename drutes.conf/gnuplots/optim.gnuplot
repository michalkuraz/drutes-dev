 A = 0.113173248618806
  k = 0.00000515892375196475
  s = 0.000855461248783132
  
  f(x) = (s/A*(1-exp(-A*sqrt(x*3600)))+k*x*3600)*100 
  
  set xrange[0:0.75]
  
  plot f(x) ,  "../../out/obspt_RE_matrix-1.out" u 1:5 w l
