# Makefile to generate fig23.pdf, fig4.pdf
# Michael Tschannen, 2016


all: fig23.pdf fig4.pdf

fig23.pdf: ceskm-distobslen.dat ceskm-distsigman.dat ceskm-obslenp.dat ceskm-obslensigman.dat ceskmfp-distobslen.dat ceskmfp-distsigman.dat ceskmfp-obslenp.dat ceskmfp-obslensigman.dat cesnnpc-distobslen.dat cesnnpc-distsigman.dat cesnnpc-obslenp.dat cesnnpc-obslensigman.dat
	pdflatex -halt-on-error -interaction=batchmode fig23.tex

ceskm-distobslen.dat:
	matlab -nojvm -r "synthPhaseDiag; exit(0);"

ceskm-distsigman.dat:
	

ceskm-obslenp.dat:
	

ceskm-obslensigman.dat:
	

ceskmfp-distobslen.dat:
	

ceskmfp-distsigman.dat:
	

ceskmfp-obslenp.dat:
	

ceskmfp-obslensigman.dat:
	

cesnnpc-distobslen.dat:
	

cesnnpc-distsigman.dat:
	

cesnnpc-obslenp.dat:
	

cesnnpc-obslensigman.dat:
	


fig4.pdf: cesnnpc.dat ceskm.dat
	pdflatex -halt-on-error -interaction=batchmode fig4.tex

cesnnpc.dat:
	matlab -nojvm -r "testeeg; exit(0);"

ceskm.dat:
	

clean:
	rm *.dat *.log *.gz *.aux 
cleanall:
	rm *.dat *.log *.pdf *.gz *.aux

