function makeclumped

%replace one of these!
fn1 = 'kalivas_1148659686_4_olam.binoct';
fn2 = 'kalivas_1148620124_4_olam.binoct';

fn1 = 'kalivas_shabbybox_RMSECVselRMSECV_1149775582_4_olam.binoct';
fn2 = 'kalivas_shabbybox_BICselBIC_1149700380_4_olam.binoct';

ld1str = sprintf('load -binary \'%s\' lambda regb blam',fn1);
eval(ld1str);
blam1 = blam;

ld2str = sprintf('load -binary \'%s\' lambda regb blam',fn2);
eval(ld2str);
blam2 = blam;

offs  = 1e12;

semilogy(abs(blam1),'r;blam1;',offs * abs(blam2), 'b;blam2;');

saveplot('blams.eps');

plotopen('blams2.eps');
multiplot(1,2);
subwindow(1,1);
semilogy(abs(blam1),'r;blam1;');
subwindow(1,2);
semilogy(abs(blam2),'r;blam2;');
oneplot();
plotclose();


endfunction
