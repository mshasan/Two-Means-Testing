libname abcd 'U:\Documents\My Research (UCF)\SAS 1st paper Final Results\Paper2';

*libname abcd 'c:\shakilpaper';

options nonotes nosource nosource2  nocenter  errors =0;

proc iml;
RESET STORAGE = "Paper";
reset autoname; start main;

s = 10000;  mu1 = 0; sigma1= 1; mu2 = .35; sigma2= 1.5; 


/*
do n1 = 30 to 30 by 20; 
do n2 = 30 to 30 by 20;
do n =  10 to 10 by 20;
do rho = 0.10 to 0.9 by 0.1;
*/

/*
do n1 = 15 to 15 by 20; 
do n2 = 15 to 15 by 20;
do n =  15 to 15 by 20;
do rho = 0.10 to 0.9 by 0.1;
*/

/*
do n1 = 10 to 10 by 20; 
do n2 = 10 to 10 by 20;
do n =  30 to 30 by 20;
do rho = 0.1 to 0.9 by 0.1;
*/

/*
do n1 = 45 to 45 by 20; 
do n2 = 45 to 45 by 20;
do n =  10 to 10 by 20;
do rho = 0.10 to 0.9 by 0.1;
*/


/*
do n1 = 10 to 10 by 20; 
do n2 = 10 to 10 by 20;
do n =  45 to 45 by 20;
do rho = 0.10 to 0.9 by 0.1;
*/


do n1 = 40 to 40 by 20; 
do n2 = 25 to 25 by 20;
do n =  10 to 10 by 20;
do rho = -0.9 to -0.10 by 0.1;


/*
do n1 = 10 to 10 by 20; 
do n2 = 15 to 15 by 20;
do n =  30 to 30 by 20;
do rho = 0.10 to 0.9 by 0.1;
*/





Z_U = j(s,1,0);  *store test scores  in these matrices;
z_c = j(s, 1,0);
wt=j(s, 1, 0);
remlt = j(s, 1, 0);
ind = j(s, 1, 0);
 


do simu = 1 to s;
 

*sample one;
  s1 = j(n1, 1, 0);
do i = 1 to n1;
s1[i,1] = rand('normal', mu1, sigma1);
end;



*sample two;
 s2 = j(n2, 1, 0);
do i = 1 to n2;
s2[i,1] = rand('normal', mu2, sigma2);
end;


*sample of n pairs;

 mean = mu1//mu2;   cov = (sigma1*sigma1 || rho*sigma1*sigma2)//(rho*sigma1*sigma2 || sigma2*sigma2);  
x1x2 = randnormal(n, mean, cov);
diff = x1x2[, 1] - x1x2[ , 2]; 




* Calculation REML ;

blk = j(n1+n2+n, 1, 0);
do i = 1 to n1+n2+n; blk[i, 1] = i; end;
reml1 = blk||j(n1+n2+n, 1, 1)||(s1[, 1]//j(n2, 1, .)//x1x2[,1]); 
reml2 = blk||j(n1+n2+n, 1, 2)||(j(n1, 1, .)//s2[, 1]//x1x2[,2]);
reml = reml1//reml2; 

call sort(reml, {1, 2});


create remldata  var{blk trt y}; 
append from reml;
close remldata;




submit;
*proc sort data = remldata; 
*by blk trt;
*run;

ods listing close;
ods output tests3= type3test; proc mixed data =remldata ;  
class blk trt; model y = trt; repeated / type = un subject = blk; run;
*ods listing;

endsubmit;

use type3test;
read all into remlcalc;
close type3test;


call delete(remldata);
call delete(type3test);



if  remlcalc[1, ncol(remlcalc)] <  0.05 then remlt[simu, 1] = 1;


* END of calculations for REML;



*calculations of our Z  formula call it Z_U;

sigma12= var(s1); sigma22= var(s2);   
 
vdiff = var(diff);  
ps = x1x2;   
s12m = cov(ps);
s12 = s12m[1,2]; 

L1num = vdiff/n;  
comsig = ((n1-1)*sigma12 + (n2-1)*sigma22)/(n1+n2-2);
L1denom = sigma12/n1 + sigma22/n2  + vdiff/n;    * unequal var; 
*L1denom = (1/n1 + 1/n2)*comsig + vdiff/n;   * equal variance;
L1 = L1num/L1denom;  

num =  L1*(sum(s1[,1])/n1- sum(s2[,1])/n2) + (1- L1)*(sum(diff))/n;
denom = sqrt(  l1*l1*(sigma12/n1 + sigma22/n2)  + (1-l1)*(1-l1)*vdiff/n  ); *unequal variances;

*denom = sqrt(  l1*l1*((1/n1 + 1/n2)*comsig)  + (1-l1)*(1-l1)*vdiff/n  ); *equal variance;

ZUscore = num/denom;

theta = sqrt((sigma12/n1 + sigma22/n2)/(vdiff/n)); *unequal variances;
*theta = sqrt((1/n1 + 1/n2)*comsig/(vdiff/n)); * equal variance;

alpha =  theta/(1 + theta);

UL =  alpha;      

*fdf = n1+n2-2; * for equal variance;

fdf = ( ( (sigma12)/n1 + (sigma22)/n2)**2 )/(  (sigma12*sigma12)/(n1*n1*(n1-1) )
                   + (sigma22*sigma22)/(n2*n2*(n2-1) ));  *for unequal variance;


UA = UL*UL*(n-1)/(n-3) + (1-UL)*(1-UL)*(fdf)/(fdf-2);
UB = (UL**4)*3*(n-1)*(n-1)/((n-3)*(n-5)) + ((1-UL)**4)*3*((fdf)**2)/((fdf-2)*(fdf-4))
           + 6*(UL**2)*((1-UL)**2)*(n-1)*(fdf)/((n-3)*(fdf-2));
Uh = sqrt( (2*UB-3*UA*UA)/(UA*UB) );
Uf = ((4*UB - 6*UA*UA)/(UB -  3*UA*UA)) ;  


finalz = zuscore*Uh*sqrt(1+theta*theta)/(1 + theta); 


*if finalz > 0.5*(tinv(0.975, uf) + tinv(0.975, uf+1)) | finalz <  - (0.5*(tinv(0.975, uf) + tinv(0.975, uf+1))
                 then  nzu[simu, 1] = 1;


if finalz > ((1-(uf- int(uf)))*tinv(0.975, int(uf)) + (uf-int(uf))*tinv(0.975, int(uf)+1)) |
                      finalz <  - ((1-(uf- int(uf)))*tinv(0.975, int(uf)) + (uf-int(uf))*tinv(0.975, int(uf)+1))
                 then   z_u[simu, 1] = 1;




*calculation of Z-corrected from Looney's paper call it Z_C;
newdiff = (n/(n1+n))*x1x2[,1] - (n/(n2+n))*x1x2[, 2];
allx = s1//x1x2[ , 1];  ally = s2//x1x2[ , 2]; 
z_cnum = sum(allx)/(n1+n)  - sum(ally)/(n2+n);
z_cdenom= sqrt( n1*sigma12/((n1+n)*(n1+n)) + n2*sigma22/((n2+n)*(n2+n)) + (var(newdiff)/n));  *unequal variance;

*z_cdenom= sqrt( (n1/((n1+n)*(n1+n)) + n2/((n2+n)*(n2+n)))*comsig   + (var(newdiff)/n));   * for equal variance;

zcscore = z_cnum/z_cdenom;
*if zcscore > tinv(0.975, n1+n2+n-3)   | zcscore <  - tinv(0.975, n1+n2+n-3)   then  z_c[simu,1]= 1;
if zcscore > 1.96  |  zcscore <  - 1.96   then  z_c[simu,1]= 1;





*calculations of weighted average of paired and unpaired t;


fdf = ( ( (sigma12)/n1 + (sigma22)/n2)**2 )/(  (sigma12*sigma12)/(n1*n1*(n1-1) ) + (sigma22*sigma22)/(n2*n2*(n2-1) ));

TL = ((fdf)*(n-3))/((n-1)*(fdf-2) + (n-3)*(fdf));
A = tl*tl*(n-1)/(n-3) + (1-tl)*(1-tl)*(fdf)/(fdf-2);
B = (tl**4)*3*(n-1)*(n-1)/((n-3)*(n-5)) + ((1-tL)**4)*3*((fdf)**2)/((fdf-2)*(fdf-4))
           + 6*(TL**2)*((1-TL)**2)*(n-1)*(fdf)/((n-3)*(fdf-2));
h = sqrt( (2*B-3*A*A)/(A*B) );
f = ((4*B - 6*A*A)/(B -  3*A*A)) ; 


wtscore =  tL*(sum(diff)/n)/(sqrt(vdiff/n)) + (1-TL)*(sum(s1)/n1 - sum(s2)/n2)/(sqrt(  (sigma12)/n1 + (sigma22)/n2) );
hT = h*wtscore;

if ht > ((1-(f- int(f)))*tinv(0.975, int(f)) + (f-int(f))*tinv(0.975, int(f)+1)) |
          ht  < - ((1-(f- int(f)))*tinv(0.975, int(f)) + (f-int(f))*tinv(0.975, int(f)+1))  then wt[simu,1]=1;

*end of weighted average of paired and unpaired t;
		  *NOTE: Use it for both TYPE I error  and POWER;




*calculations of independent  t;

ind_tnum = sum(allx)/(n1+n)  - sum(ally)/(n2+n);
ind_tdenom= sqrt( (1/(n1+n) + 1/(n2+n))*(((n1+n-1)*sigma12 + ( n2+n-1)*sigma22)/(n1+n +n2 +n -2))  );
ind_t = ind_tnum/ind_tdenom;
if ind_t > tinv(0.975, n1+n2+n+n-2)  | ind_t < - tinv(0.975, n1+n2+n+n-2)   then  ind[simu,1]= 1;


end;

*calculate type1 error rates;


remltype1 = sum(remlt)/s;
zutype1 = sum(z_u)/s;
zctype1 = sum(z_c)/s;
wttype1 = sum(wt)/s;
indtype1 =sum(ind)/s;

result = result//(n1||n2||n||rho ||  zutype1 ||zctype1 ||wttype1 || indtype1 || remltype1 ) ;

end;end;end;end;

store result;

print result;

finish main;    run main; 


/*Transform matrix into dataset*/
proc iml;
RESET STORAGE = "Paper";
	LOAD result;
create abcd.power6 var _All_;
append;
quit;

proc transpose data=abcd.power6 out=abcd.power_6;
run;




