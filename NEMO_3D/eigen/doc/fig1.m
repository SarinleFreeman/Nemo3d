Az = [
   13456    54
   53824   109
  215296   567
  861184  6507
 3444736 46942
];

Ld = [
   13456    38
   53824    68
  215296   246
  861184  2321
 3444736 39093
];

Lex1 = [
   13456    22
   53824    21
  215296    51
  861184   269
 3444736  1783
13778944 16772
55115776 149606
];

Lex2 = [
   13456	   20
   53824	   30
  215296	   38
  861184	  243
 3444736	 1001
13778944	 8317
];


p_Az = loglog(Az(:,1), Az(:,2), '-+');
hold on
p_Ld = loglog(Ld(:,1), Ld(:,2), '-s');
p_Lex1 = loglog(Lex1(:,1), Lex1(:,2), '-o');
p_Lex2 = loglog(Lex2(:,1), Lex2(:,2), '-v');

legend([p_Az, p_Ld, p_Lex1, p_Lex2],...
  'PARPACK Arnoldi pzsdrv1', ...
  'PARPACK Lanczos pdsdrv1', ...
  'JPLNEMO Lanczos pexample', ...
  'JPLNEMO Lanczos pexample (single call to dstebz)');

xlabel('matrix rank');
ylabel('time (s)');

% prettyPlot('fs=16; lgs=12; ms=8; sty=''notebook'';')
% print -deps fig1.eps


proc = [2 3 6 9 18 27 54];

% dim = 7.46e5, 0.8 GHz
SS = [
567
440
284
214
172
175
173
];

% dim = 3e6, 0.8 GHz
LS = [
  4463
  3378
  2184
  1726
  2931
  2221
  1587
];

% dim = 3e6, 2.2 GHz
LF = [
 	1511	
	1054	
	 749	
	 616	
	 543	
	 541	
	 528	
];


p_SS = loglog(proc, SS, '-o');
hold on;
p_LS = loglog(proc, LS, '-v');
p_LF = loglog(proc, LF, '-s');

legend([p_SS, p_LS, p_LF],...
  'size=7.5 \times 10^5 (0.8 GHz)', ...
  'size=3.0 \times 10^6 (0.8 GHz)', ...
  'size=3.0 \times 10^6 (2.2 GHz)');

xlabel('number of processors');
ylabel('time (s)');

% prettyPlot('fs=16; lgs=12; ms=8; sty=''notebook'';')
% print -deps fig2.eps

