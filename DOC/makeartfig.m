function makeartfig

%replace one of these!
fn1 = 'artificial_shabbybox_BICselRMSECV_1149536064_3_olam.binoct';
fn2 = 'artificial_shabbybox_BICselRMSECV_1149536064_3_finmod.binoct';

ld1str = sprintf('load -binary \'%s\' lambda blam',fn1);
eval(ld1str);
ld2str = sprintf('load -binary \'%s\' lam beta wvls colnums',fn2);
eval(ld2str);


blam	= abs(blam);

semilogy(blam,'c;blam;',colnums,blam(colnums),'r+;;');
saveplot('agoodart.eps');

endfunction
