
allfigs = findobj('Type','Figure')
d = 1
for i = 1:length(allfigs)
	strd = num2str(d) 
		hgexport(i, ['suppleFig' strd '.eps], hgexport('factorystyle'), 'Format', 'eps');
	d=d+1;
end;
