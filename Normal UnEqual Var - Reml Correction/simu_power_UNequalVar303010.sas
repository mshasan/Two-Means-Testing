libname abcd 'U:\Documents\My Research (UCF)\SAS 1st paper Final Results\Paper';

*libname abcd 'c:\shakilpaper';

options nonotes nosource nosource2  nocenter  errors =0;

proc iml;
RESET STORAGE = "Paper";
reset autoname; start main;

s = 10000;  mu1 = 0; sigma1= 1; mu2 = .35; sigma2= 2; 



do n1 = 30 to 30 by 20; 
do n2 = 30 to 30 by 20;
do n =  10 to 10 by 20;
do rho = -.9 to 0.9 by 0.2;

/*
do n1 = 10 to 10 by 20; 
do n2 = 10 to 10 by 20;
do n =  30 to 30 by 20;
do rho = -.9 to 0.9 by 0.2;
*/

/*
do n1 = 45 to 45 by 20; 
do n2 = 45 to 45 by 20;
do n =  10 to 10 by 20;
do rho = -.9 to 0.9 by 0.2;
*/


/*
do n1 = 10 to 10 by 20; 
do n2 = 10 to 10 by 20;
do n =  45 to 45 by 20;
do rho = -.9 to 0.9 by 0.2;
*/

/*
do n1 = 40 to 40 by 20; 
do n2 = 25 to 25 by 20;
do n =  10 to 10 by 20;
do rho = -.9 to 0.9 by 0.2;
*/

/*
do n1 = 10 to 10 by 20; 
do n2 = 15 to 15 by 20;
do n =  30 to 30 by 20;
do rho = -.9 to 0.9 by 0.2;
*/





*store test scores  in these matrices;
remlt = j(s, 1, 0);

 
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
class blk trt; model y = trt / ddfm = kenwardroger; repeated / type = un subject = blk; run;
*ods listing;


endsubmit;

use type3test;
read all into remlcalc;
close type3test;


call delete(remldata);
call delete(type3test);



if  remlcalc[1, ncol(remlcalc)] <  0.05 then remlt[simu, 1] = 1;


* END of calculations for REML;
end;

*calculate type1 error rates;


remltype1 = sum(remlt)/s;

result = result//(n1||n2||n||rho || remltype1 ) ;

end;end;end;end;

store result;

print result;


* create sas data set from magrix result;

create resultdata from result[colname={"n1"  "n2"  "n"   "rho"  "remltype1" }] ;
append from result;
close resultdata;

* end of sas data set;

finish main;    run main; 




*Sending SAS data file  to your harddisk as Excel file;
PROC EXPORT DATA= resultdata 
            OUTFILE= 'U:\Documents\My Research (UCF)\SAS 1st paper Final Results\Paper\Power_303010.xlsx' 
            DBMS=xlsx 
            REPLACE;
RUN;  *End of sending  SAS data file to your harddisk as Excel file;





