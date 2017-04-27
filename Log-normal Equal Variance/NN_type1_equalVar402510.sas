libname abcd 'U:\Documents\My Research (UCF)\SAS 1st paper Final Results\Paper';

options nonotes nosource nosource2  nocenter  errors =0;


proc iml;
RESET STORAGE = "paper";
reset autoname; start main;

s = 5000;  mu1 = 0; sigma1= 1; mu2 = 0; sigma2= 1; 



/*
do n1 = 30 to 30 by 20; 
do n2 = 30 to 30 by 20;
do n =  10 to 10 by 20;
do rhoL = -0.9 to 0.9 by 0.1;
*/

/*
do n1 = 15 to 15 by 20; 
do n2 = 15 to 15 by 20;
do n =  15 to 15 by 20;
do rhoL = -0.9 to 0.9 by 0.1;
*/

/*
do n1 = 10 to 10 by 20; 
do n2 = 10 to 10 by 20;
do n =  30 to 30 by 20;
do rhoL = -0.9 to 0.9 by 0.1;
*/


/*
do n1 = 45 to 45 by 20; 
do n2 = 45 to 45 by 20;
do n =  10 to 10 by 20;
do rhoL = -0.9 to 0.9 by 0.1;
*/


/*
do n1 = 10 to 10 by 20; 
do n2 = 10 to 10 by 20;
do n =  45 to 45 by 20;
do rhoL = -0.9 to 0.9 by 0.1;
*/


do n1 = 40 to 40 by 20; 
do n2 = 25 to 25 by 20;
do n =  10 to 10 by 20;
do rhoL = -0.9 to 0.9 by 0.1;


/*
do n1 = 10 to 10 by 20; 
do n2 = 15 to 15 by 20;
do n =  30 to 30 by 20;
do rhoL = -0.9 to 0.9 by 0.1;
*/

if (rhoL > -1/(sqrt((exp(sigma1*sigma1) - 1)*(exp(sigma2*sigma2) -1)))) then do;
if  ( log( 1 + rhoL*(sqrt((exp(sigma1*sigma1) - 1)*(exp(sigma2*sigma2) -1)) ))/(sigma1*sigma2)  < 1) & 
 (log( 1 + rhoL*(sqrt((exp(sigma1*sigma1) - 1)*(exp(sigma2*sigma2) -1)) ))/(sigma1*sigma2)  > -1 )   then do; 


Z_U = j(s,1,0);  *store test scores  in these matrices;
Z_UU = j(s,1,0);
z_c = j(s, 1,0);
z_cU = j(s, 1,0);
wt=j(s, 1, 0);
wtU = j(s, 1, 0);
remlt = j(s, 1, 0);
ind = j(s, 1, 0);
indU = j(s, 1, 0);
x1x2 = j(n, 2, 0);
 


do simu = 1 to s;
 

*sample one;
  s1 = j(n1, 1, 0);
do i = 1 to n1;
s1[i,1] =  exp(rand('normal', mu1, sigma1));
end;



*sample two;
 s2 = j(n2, 1, 0);
do i = 1 to n2;
s2[i,1] = exp(rand('normal', mu2, sigma2));
end;


*sample of n pairs;


rho = log(1 + rhoL*sqrt(  (exp(sigma1*sigma1) - 1)*(exp(sigma2*sigma2) -1)   ))/(sigma1*sigma2); 


mean = mu1//mu2;   cov = (sigma1*sigma1 || rho*sigma1*sigma2)//(rho*sigma1*sigma2 || sigma2*sigma2);  
x1x2normal = randnormal(n, mean, cov);  * random sample from bivariate normal;
x1x2[, 1] = exp(x1x2normal[, 1]);
x1x2[, 2] = exp(x1x2normal[ , 2]);  * x1x2 is now a paired random sample from bivariate lognormal;
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

ods  listing close;
ods output tests3= type3test; proc mixed data =remldata ;  
class blk trt; model y = trt; repeated / type = ar(1) subject = blk; run;
*ods  listing;

endsubmit;

use type3test;
read all into remlcalc;
close type3test;


call delete(remldata);
call delete(type3test);



if  remlcalc[1, ncol(remlcalc)] <=  0.05 then remlt[simu, 1] = 1;


* END of calculations for REML;





*calculations of our Z  formula call it Z_U  Z_UU;

sigma12= var(s1); sigma22= var(s2);   
 
vdiff = var(diff);  
ps = x1x2;   
s12m = cov(ps);
s12 = s12m[1,2]; 


fdf = ( ( (sigma12)/n1 + (sigma22)/n2)**2 )/(  (sigma12*sigma12)/(n1*n1*(n1-1) ) + (sigma22*sigma22)/(n2*n2*(n2-1) ));


L1num = vdiff/n;  
comsig = ((n1-1)*sigma12 + (n2-1)*sigma22)/(n1+n2-2);
L1denomU = sigma12/n1 + sigma22/n2  + vdiff/n;    * unequal var; 
L1denom = (1/n1 + 1/n2)*comsig + vdiff/n;   * equal variance;
L1 = L1num/L1denom;  
L1U = L1num/L1denomU;
num =  L1*(sum(s1[,1])/n1- sum(s2[,1])/n2) + (1- L1)*(sum(diff))/n;
numU= L1U*(sum(s1[,1])/n1- sum(s2[,1])/n2) + (1- L1U)*(sum(diff))/n;
denom = sqrt(  l1*l1*((1/n1 + 1/n2)*comsig)  + (1-l1)*(1-l1)*vdiff/n  );
denomU = sqrt(  l1U*l1U*(sigma12/n1 + sigma22/n2)  + (1-l1U)*(1-l1U)*vdiff/n  );

ZUscore = num/denom;
ZUscoreU = numU/denomU;


theta = sqrt((1/n1 + 1/n2)*comsig/(vdiff/n));
thetaU = sqrt( (sigma12/n1 + sigma22/n2)/(vdiff/n));

alpha =  theta/(1 + theta);
alphaU =  thetaU/(1 + thetaU);

UL =  alpha;      
ULU =  alphaU;      

UA = UL*UL*(n-1)/(n-3) + (1-UL)*(1-UL)*(n1+n2-2)/(n1+n2-4);

UAU = ULU*ULU*(n-1)/(n-3) + (1-ULU)*(1-ULU)*(fdf)/(fdf-2);

UB = (UL**4)*3*(n-1)*(n-1)/((n-3)*(n-5)) + ((1-UL)**4)*3*((n1+n2-2)**2)/((n1+n2-4)*(n1+n2-6))
           + 6*(UL**2)*((1-UL)**2)*(n-1)*(n1+n2-2)/((n-3)*(n1+n2-4));

UBU = (ULU**4)*3*(n-1)*(n-1)/((n-3)*(n-5)) + ((1-ULU)**4)*3*((fdf)**2)/((fdf-2)*(fdf-4))
           + 6*(ULU**2)*((1-ULU)**2)*(n-1)*(fdf)/((n-3)*(fdf-2));

Uh = sqrt( (2*UB-3*UA*UA)/(UA*UB) );
UhU = sqrt( (2*UBU-3*UAU*UAU)/(UAU*UBU) );

Uf = ((4*UB - 6*UA*UA)/(UB -  3*UA*UA)) ;  

UfU = ((4*UBU - 6*UAU*UAU)/(UBU -  3*UAU*UAU)) ;  


finalz = zuscore*Uh*sqrt(1+theta*theta)/(1 + theta); 

finalzU = zuscoreU*UhU*sqrt(1+thetaU*thetaU)/(1 + thetaU); 


if finalz > ((1-(uf- int(uf)))*tinv(0.975, int(uf)) + (uf-int(uf))*tinv(0.975, int(uf)+1)) |
                      finalz <  - ((1-(uf- int(uf)))*tinv(0.975, int(uf)) + (uf-int(uf))*tinv(0.975, int(uf)+1))
                 then   z_u[simu, 1] = 1;


if finalzU > ((1-(ufU- int(ufU)))*tinv(0.975, int(ufU)) + (ufU-int(ufU))*tinv(0.975, int(ufU)+1)) |
                      finalzU <  - ((1-(ufU- int(ufU)))*tinv(0.975, int(ufU)) + (ufU-int(ufU))*tinv(0.975, int(ufU)+1))
                 then   z_uU[simu, 1] = 1;

*End of calculation of Z_U and Z_UU - equal and unequal variances;








*calculation of Z-corrected from Looney's paper call it Z_C - equal and unequal variances;
newdiff = (n/(n1+n))*x1x2[,1] - (n/(n2+n))*x1x2[, 2];
allx = s1//x1x2[ , 1];  ally = s2//x1x2[ , 2]; 
z_cnum = sum(allx)/(n1+n)  - sum(ally)/(n2+n);
z_cdenomU = sqrt( n1*sigma12/((n1+n)*(n1+n)) + n2*sigma22/((n2+n)*(n2+n)) + (var(newdiff)/n));  *unequal variance;

z_cdenom= sqrt( (n1/((n1+n)*(n1+n)) + n2/((n2+n)*(n2+n)))*comsig   + (var(newdiff)/n));   * for equal variance;

zcscore = z_cnum/z_cdenom;
zcscoreU = z_cnum/z_cdenomU;

*if zcscore > tinv(0.975, n1+n2+n-3)   | zcscore <  - tinv(0.975, n1+n2+n-3)   then  z_c[simu,1]= 1;
if zcscore > 1.96   | zcscore <  - 1.96   then  z_c[simu,1]= 1;

*if zcscoreU > tinv(0.975, n1+n2+n-3)   | zcscoreU <  - tinv(0.975, n1+n2+n-3)   then  z_cU[simu,1]= 1;

if zcscoreU > 1.96   | zcscoreU <  - 1.96   then  z_cU[simu,1]= 1;

*end of calculation of Z-corrected from Looney's paper call it Z_C - equal and unequal  variances assumed;





*calculations of weighted average of paired and unpaired t : equal variance assumed;

*tL = (n-1)/(n1+n2+n-3);
TL = ((n1+n2-2)*(n-3))/((n-1)*(n1+n2-4) + (n-3)*(n1+n2-2));
A = tl*tl*(n-1)/(n-3) + (1-tl)*(1-tl)*(n1+n2-2)/(n1+n2-4);
B = (tl**4)*3*(n-1)*(n-1)/((n-3)*(n-5)) + ((1-tL)**4)*3*((n1+n2-2)**2)/((n1+n2-4)*(n1+n2-6))
           + 6*(TL**2)*((1-TL)**2)*(n-1)*(n1+n2-2)/((n-3)*(n1+n2-4));
h = sqrt( (2*B-3*A*A)/(A*B) );
f = ((4*B - 6*A*A)/(B -  3*A*A)) ; 

wtscore =  tL*(sum(diff)/n)/(sqrt(vdiff/n)) + (1-TL)*(sum(s1)/n1 - sum(s2)/n2)/(sqrt(comsig*(1/n1 + 1/n2)));


hT = h*wtscore;

if ht > ((1-(f- int(f)))*tinv(0.975, int(f)) + (f-int(f))*tinv(0.975, int(f)+1)) |
          ht  < - ((1-(f- int(f)))*tinv(0.975, int(f)) + (f-int(f))*tinv(0.975, int(f)+1))  then wt[simu,1]=1;

*end of calculations of weighted average of paired and unpaired t : equal variance assumed;




*calculations of weighted average of paired and unpaired t- unequal variance assumed;

fdf = ( ( (sigma12)/n1 + (sigma22)/n2)**2 )/(  (sigma12*sigma12)/(n1*n1*(n1-1) ) + (sigma22*sigma22)/(n2*n2*(n2-1) ));

TL = ((fdf)*(n-3))/((n-1)*(fdf-2) + (n-3)*(fdf));
A = tl*tl*(n-1)/(n-3) + (1-tl)*(1-tl)*(fdf)/(fdf-2);
B = (tl**4)*3*(n-1)*(n-1)/((n-3)*(n-5)) + ((1-tL)**4)*3*((fdf)**2)/((fdf-2)*(fdf-4))
           + 6*(TL**2)*((1-TL)**2)*(n-1)*(fdf)/((n-3)*(fdf-2));
h = sqrt( (2*B-3*A*A)/(A*B) );
f = ((4*B - 6*A*A)/(B -  3*A*A)) ; 


wtscoreU =  tL*(sum(diff)/n)/(sqrt(vdiff/n)) + (1-TL)*(sum(s1)/n1 - sum(s2)/n2)/(sqrt(  (sigma12)/n1 + (sigma22)/n2) );
hT = h*wtscoreU;

if ht > ((1-(f- int(f)))*tinv(0.975, int(f)) + (f-int(f))*tinv(0.975, int(f)+1)) |
          ht  < - ((1-(f- int(f)))*tinv(0.975, int(f)) + (f-int(f))*tinv(0.975, int(f)+1))  then wtU[simu,1]=1;


* end of calculations of weighted average of paired and unpaired t- unequal variance assumed;


		 





*calculations of independent  t;

fdf =  ((var(allx)/(n1+n) + var(ally)/(n2+n))**2 )/(((var(allx)/(n1+n))**2)/(n1+n-1)   + ((var(ally)/(n2+n))**2)/(n2+n-1));

ind_tnum = sum(allx)/(n1+n)  - sum(ally)/(n2+n);
ind_tdenom= sqrt( (1/(n1+n) + 1/(n2+n))*(((n1+n-1)*sigma12 + ( n2+n-1)*sigma22)/(n1+n +n2 +n -2))  );
ind_tdenomU= sqrt( var(allx)/(n1+n) +  var(ally)/(n2+n) );


ind_t = ind_tnum/ind_tdenom;

ind_tU = ind_tnum/ind_tdenomU;

if ind_t > tinv(0.975, n1+n2+n+n-2)  | ind_t < - tinv(0.975, n1+n2+n+n-2)   then  ind[simu,1]= 1;

if ind_tU > tinv(0.975, int(fdf) )  | ind_tU < - tinv(0.975, int(fdf))   then  indU[simu,1]= 1;

end;

*calculate type1 error rates;


remltype1 = sum(remlt)/s;
zutype1 = sum(z_u)/s;
zutype1U = sum(z_uU)/s;
zctype1 = sum(z_c)/s;
zctype1U = sum(z_cU)/s;
wttype1 = sum(wt)/s;
wttype1U= sum(wtU)/s;
indtype1 =sum(ind)/s;
indtype1U =sum(indU)/s;

result = result//(n1||n2||n||rhoL ||  zutype1 ||zctype1 ||wttype1 || indtype1 || remltype1 || wttype1U || Zctype1U || Zutype1U || indtype1U) ;

end; end;
end;end;end;end;


store result;

print result;

* create sas data set from magrix result;
create resultdata from result[colname={"n1" "n2" "n" "rhoL" "zutype1" "zctype1" "wttype1" "indtype1"
"remltype1" "wttype1U" "Zctype1U" "Zutype1U" "indtype1U"}] ;
append from result;
close resultdata;


finish main;    run main;

*Sending SAS data file "Daily" to your harddisk as Excel file;
PROC EXPORT DATA= resultdata
OUTFILE= 'U:\Documents\My Research (UCF)\SAS 1st paper Final Results\Paper\typeI_402510.xlsx'
DBMS=xlsx
REPLACE;
RUN; *End of sending SAS data file to your harddisk as Excel file; 






