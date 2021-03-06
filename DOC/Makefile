
######################
######################
# makefile generated 
# wsel makefile 
# created by s.e.pav 
# $Id: Makefile 89 2006-01-26 20:59:08Z spav $
######################
######################

############### FLAGS ###############

# these don't work at times
# in that case set these by hand?

LATEX       := $(shell which latex)
BIBTEX      := $(shell which bibtex)
HTLATEX     := $(shell which htlatex.sh)
MAKEINDEX   := $(shell which makeindex)
GHOSTVIEW   := $(shell which gv)
EDITOR      := $(shell which vim)
PDF_VIEWER  := $(shell which xpdf)
ACROREAD		:= $(shell which acroread)
KGV     		:= $(shell which kghostview)
TODAY       := $(shell date +'%Y%b%d%H%m')
TIME        := $(shell perl -e 'print time();')
PAGER   		:= $(shell which less)
TAR         := $(shell which tar)
ASPELL  		:= $(shell which aspell)

TEXINPADD    = .:$(HOME)/sys/etc/tex:$(HOME)/sys/etc/tex/SEPtex:$(HOME)/work/math/TEX
PRETEX       = TEXINPUTS=$(TEXINPADD):$$TEXINPUTS
PREBIB       = BSTINPUTS=$(TEXINPADD):$$BSTINPUTS \
               BIBINPUTS=$(TEXINPADD):$$BIBINPUTS 

PREIDX       = INDEXSTYLE=$(TEXINPADD):$$INDEXSTYLE

#undoes psfrag for pdf
UNPSFRAG		 = perl $(HOME)/sys/bin/unpsfrag.pl
#unroll commands
DETEXIFY		 = perl $(HOME)/sys/perl/detexify.pl

SCREEN_SIZE  = normal
include	$(HOME)/sys/etc/.Makefile.local

PROJECT      = wsel
TEX_SOURCE   = $(PROJECT).tex
BIB_SOURCE   = $(PROJECT).bib
DVI_TARGET   = $(PROJECT).dvi
PS_TARGET    = $(PROJECT).ps 
PDF_TARGET   = $(PROJECT).pdf
BBLS         = $(PROJECT).bbl

#sub figures
QPS_FILES    := $(shell ls figs/*.qps)
JPG_FILES    := $(shell ls pics/*.jpg)
EPS_FILES    := $(patsubst %.qps,%.eps,$(QPS_FILES)) 
EPS_FILES    += $(patsubst %.jpg,%.eps,$(JPG_FILES)) 
SUB_FIGS     = $(EPS_FILES)

#SAVE
# tracked projects
PROJECTS     = $(PROJECT) $(PROJECT)_dis
#UNSAVE
# add on dependencies (subchapters of wsel)
TEX_EXTRAS   = 
# nonlocal dependencies
STY_FILES    = 

HOME_HOST    = shabbychef.ucsd.edu
DIR_STUB     = work
DOC_DIR      = $(DIR_STUB)/$(PROJECT)/DOC

# more worry
OWNER        = $(shell whoami)
THIS_HOST    = $(shell hostname)
CWD          = $(shell pwd)

# tar junk
BACKUP       = $(PROJECT)_backup_$(TODAY)
TAR_FILE     = $(BACKUP).tar
TARBZ_FILE   = $(TAR_FILE).bz2
TAR_MAJORS   = Makefile $(TEX_SOURCE) $(BIB_SOURCE) $(TEX_EXTRAS) 
TAR_MINORS   = *.ps *.eps *.fig *.pro *.sty 

STUB         = $(HOME)/sys/STUB
EMAIL        = spav@ucsd.edu

# cvs stuff
CVS          = cvs
# for external cvs directory on bates
CVS_REP      = /usr0/glmiller/CVS-SANGRIA
CVS_HOST     = bates.pscico.cs.cmu.edu
CVS_ROOT_DIR = :ext:$(CVS_HOST):$(CVS_REP)
# for external cvs directory on scicomp
CVS_REP      = /usr/people/spav/repository/private
CVS_HOST     = scicomp.ucsd.edu
CVS_ROOT_DIR = :ext:$(CVS_HOST):$(CVS_REP)
# for cvs directory locally on scicomp.
CVS_ROOT_DIR = $(HOME)/repository/private

# what to put into cvs
CVS_MAJORS   = $(TAR_MAJORS) $(TAR_MINORS) $(SUB_FIGS)
CVS_STUB     = spav/$(DIR_STUB)

# X geometry
X_GEOM_TINY  = -geometry 975x760+15+10
X_GEOM_NORM  = -geometry 1066x920+15+10
X_GEOM_WIDE  = -geometry 1230x920+15+10
X_GEOM_HUGE  = -geometry 1800x1450+15+10
X_GEOM_HUGE  = -geometry 1514x1165+1920+0
X_GEOM_HUGE  = -geometry 1834x1107+1920+0

# xdvi
XDVI_NORM    = -s 5
XDVI_HUGE    = -s 3

# ghostview
GVIEW_NORM   = --scalebase=2 --scale=2 
GVIEW_HUGE   = --scalebase=2 --scale=3

# now select
ifeq ($(SCREEN_SIZE),big)
X_GEOM_PORT  = $(X_GEOM_HUGE)
X_GEOM_LAND  = $(X_GEOM_WIDE)
XDVI_SCALE   = $(XDVI_HUGE)
GVIEW_SCALE  = $(GVIEW_HUGE)
else
X_GEOM_PORT  = $(X_GEOM_NORM)
X_GEOM_LAND  = $(X_GEOM_WIDE)
XDVI_SCALE   = $(XDVI_NORM)
GVIEW_SCALE  = $(GVIEW_NORM)
endif

# dvips
DVIPS_FLAGS  = -R -t letter
DVIEPS_FLAGS = -R -t letter -Pcmz -Pamz
DVIPDF_FLAGS = -R -t letter -Pcmz -Pamz

# xdvi
XDVI_FLAGS   = -safer $(XDVI_SCALE) -expert -hush $(X_GEOM_PORT) -paper us \
								-keep -postscript 0

# ghostview 
GVIEW_COMMON = --safer --spartan --media=letter $(GVIEW_SCALE) --resize 
# ghostview landscape
GVIEW_LAND   = $(GVIEW_COMMON) $(X_GEOM_LAND) --orientation=landscape
# ghostview portrait
GVIEW_PORT   = $(GVIEW_COMMON) $(X_GEOM_PORT) --orientation=portrait
GVIEW_FLAGS  = $(GVIEW_PORT)

# pdfviewer
XPDF_FLAGS   = $(X_GEOM_PORT) -z width -bg grey
KPDF_FLAGS   = -bg black 
PDFSEE_FLAGS = -z width -q -bg black -fullscreen
PDFSEE_FLAGS = -z width -q
PDFSEE_FLAGS = $(XPDF_FLAGS)

#kghostview 
KGV_FLAGS    = $(X_GEOM_PORT)

#aspell
ASPELL_FLAGS = 

#dot
DOT_FLAGS    = 

############## DEFAULT ##############

default : all

############## MARKERS ##############

.PHONY   : 
.SUFFIXES: .tex .bib .dvi .ps .pdf .eps
.PRECIOUS: %.dvi %.ps %.pdf %.jpg %.gif 

############ BUILD RULES ############

# compile and convert
%.dvi : %.tex $(STY_FILES)
		$(PRETEX) $(LATEX) $<
		if grep Citation $*.log > /dev/null; then $(PREBIB) $(BIBTEX) $*; $(PRETEX) $(LATEX) $*; fi
		if grep Rerun $*.log > /dev/null; then $(PRETEX) $(LATEX) $*; fi

%.dep.dot : %.tex %.dvi $(STY_FILES)
		echo digraph G \{ > $@;
		$(PRETEX) $(LATEX) $< | grep 'DEP:' | perl -pe 's/DEP://;s/\\hbox {}//g;' >> $@
		echo \} >> $@;

%.ps : %.dvi
		dvips $(DVIPS_FLAGS) -o $@ $< 

%.pdf : %.dvi
		dvips $(DVIPDF_FLAGS) -o $*.ps $< 
		ps2pdf $*.ps
		-rm -f $*.ps

# converting to pdf the hard way

%.pdftex : %.tex
		$(UNPSFRAG) $<
		cp _$< $@

%.pdf : %.pdftex
		$(PDFLATEX) $<

%.mergepdf :
		gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@ $(*)*.pdf

# roll commands from sty files

%.detex : %.tex
		@-if ! [ -r $@ ]; then cp $< $@; fi
		$(DETEXIFY) $@ >> $@;

# crappy ps from pdf.
%.craps : %.pdf
		a2ps -X'ascii' -M'Letter' -R --columns=1 -o $*.ps $<

# crappy ps from edp
%.craps : %.edp
		a2ps -X'ascii' -M'Letter' -R --columns=2 -o $*.ps $<

# reverse?
%.pdf : %.eps
		epstopdf $<

# crappy ascii anyone?
%.txt : %.ps
		ps2ascii $< > $@

# crappy html anyone?
%.html : %.tex
		$(HTLATEX) $< 
		-rm -f tex4ht.fls $*.4ct $*.4tc $*.lg $*.idv $*.xref $*.tmp

# tex extras
%.bbl : %.bib
		$(PREBIB) $(BIBTEX) $*

%.bbl : %.aux
		$(PREBIB) $(BIBTEX) $*

%.ind : %.idx
		$(PREIDX) $(MAKEINDEX) $*

# qt3 kseg ps files have unicode text.
# this is nogood for psfrag texthooks.
# change them here.

%.pps : %.qps
		@echo funky perl func to convert $< to hookable $@;
		@perl -n -e 'BEGIN { @FE=();$$HASF=0; } \
		if (m{^.*ENC-.*\[}) { \
			print;$$esc = 0;$$HASF ||= 1; \
			while (($$line = <>) && not $$esc) \
			{ print $$line;chomp($$line); \
				while ($$line =~ s{/([^\/\]]+)}{}) {push (@FE,$$1);} \
				$$esc = ($$line =~ m{\]\s*def\s*$$});} \
			%mymaps = qw( one 1 two 2 three 3 four 4 five 5 \
			  six 6 seven 7 eight 8 nine 9 zero 0 \
        period . equal = plus + );$$mymaps{space}=q[ ]; \
			while (($$k,$$v)=each(%mymaps)) { map {s/^$$k$$/$$v/;} @FE; } \
			print qq[/myT{gsave PCol SC ty MT 1 index dup length exch stringwidth pop 1 -1 roll -1 1 scale 180 rotate exch sub exch div exch 0 exch ashow grestore}D\n]; \
		} elsif ((m{(\d{1,3}\s+Y)?<([0-9A-F]+)>(-?\d+)\s+(-?\d+)\s+A?T\s*$$}) && $$HASF) \
		{ print qq[/Times-Roman findfont 18 scalefont setfont\n]; \
			$$cod = $$2;$$str = q[]; \
			while ($$cod) { $$ch = substr($$cod,0,4,q[]); \
										  $$str .= (hex($$ch))? $$FE[hex($$ch)]:q[ ]; } \
			print qq[$$1($$str)$$3 $$4 myT\n]; \
		} else { print $$_; }' < $< > $@;

# define a rule to find a bounding box

define findbb
@echo find real bounding box for $<;
@echo -n '   was: ';
@if grep Bounding $<; then\
perl -p -e 'BEGIN { $$spit = 0;\
$$bb = qx[gs -q -dNOPAUSE -dSAFER -dDELAYSAFER -dBATCH -sDEVICE=bbox -sOutputFile=/dev/null $< 2>&1 | awk "/\%\%Bounding/ {print}"];} \
$$spit ||= (/^%%Bounding/ and $$_ = qq[$$bb]);' $< > $@;\
else\
perl -p -e 'BEGIN { $$spit = 0;\
$$bb = qx[gs -q -dNOPAUSE -dSAFER -dDELAYSAFER -dBATCH -sDEVICE=bbox -sOutputFile=/dev/null $< 2>&1 | awk "/\%\%Bounding/ {print}"];} \
$$spit ||= ($$_ = qq[$$bb$$_]);' $< > $@;\
fi
@echo -n 'is now: ';
@head $@ | grep Bounding
endef

# make correct bounding box
%.eps : %.pps 
		$(findbb)

%.eps : %.ps 
		$(findbb)

# define a rule to convert color kseg figs to black and white
define convert_to_bw
@echo converting color $< to black and white $@;
@perl -p -e 'BEGIN { sub togrey { my $$avg = int(($$_[0] + $$_[1] + $$_[2])/3.0);\
return qq[$$avg $$avg $$avg]; } }\
s/^(\S+)\s(\S+)\s(\S+)(\sP1)$$/&togrey($$1,$$2,$$3) . $$4/e;\
s/^(\S+\s)(\S+)\s(\S+)\s(\S+)(\sBR)$$/$$1 . &togrey($$2,$$3,$$4) . $$5/e;\
s/^(\S+\s\S+\s)(\S+)\s(\S+)\s(\S+)(\s\S+\s\S+\sPE)$$/$$1 . &togrey($$2,$$3,$$4) . $$5/e;' $< > $@;
endef

# make a color fig generated by kseg into a black and white...
%.bw.ps : %.ps
		$(convert_to_bw)

%.bw.eps : %.eps
		$(convert_to_bw)

# define a rule to kill first page from a ps file.
# something wrong with octave -> gnuplot for multiplots

define killpageone
@echo killing page one in $<, output to $@;
@perl -n -e 'BEGIN { $$INPONE = 0; }\
/%%Page: 1 1/ && ($$INPONE = 1);\
/%%Page: 2 2/ && ($$INPONE = 0);\
s/%%Page: 2 2/%%Page: 1 1/;\
s/%%Pages: 2/%%Pages: 1/;\
$$INPONE || print;' $< > $@;
endef

# make a color fig generated by kseg into a black and white...
%.pp1.ps : %.ps
		$(killpageone)

%.pp1.eps : %.eps
		$(killpageone)

#get all the pages from a postscript file
%.pages.ps : %.ps
		gs -dNOPAUSE -dBATCH -sDEVICE=pswrite -sOutputFile="$*_page%d.ps" $<;

%.pages.eps : %.eps
		gs -dNOPAUSE -dBATCH -sDEVICE=pswrite -sOutputFile="$*_page%d.eps" $<;


# image conversions
%.ps	: %.gif
		gif2ps $< > $@

%.pnm	: %.jpg
		jpegtopnm $< > $@

%.pnm	: %.jpeg
		jpegtopnm $< > $@

%.ps	: %.pnm
		pnmtops $< > $@

%.png : %.pnm
		pnmtopng $< > $@

# convert from xfig
%.eps : %.fig
		/bin/rm -f $@ 
		fig2dev -L ps $< > $@

# make a graph
%.dot : %.flc
		perl $(HOME)/sys/perl/flowcharter.pl $< > $@

%.ps : %.dot
		dot $(DOT_FLAGS) -Tps $< -o $@

%.fig : %.dot
		dot $(DOT_FLAGS) -Tfig $< -o $@

# viewing targets
%.xsee : %.dvi
		xdvi $(XDVI_FLAGS) $*.dvi
		
%.psee : %.ps
		$(GHOSTVIEW) $(GVIEW_FLAGS) $*.ps

%.pdfsee : %.pdf
		$(PDF_VIEWER) $(PDFSEE_FLAGS) $*.pdf

# check a document
%.chk : %.dup %.spell

# check spelling
%.spell : %.tex
		$(ASPELL) $(ASPELL_FLAGS) --dont-tex-check-comments -t -l < $< | sort | uniq | $(PAGER)

# check duplicate words
%.dup : %.tex
		perl -an -F/\\s+/ -e 'BEGIN { $$last = q[]; $$line = 0; $$prevline = q[];}\
		$$line++;$$first = 1;\
		foreach $$word (@F) {\
		if ($$word eq $$last) {\
		if ($$first) { print qq[duplicate $$word, lines ],($$line-1),qq[-$$line:\n$$prevline$$_]; }\
		else { print qq[duplicate $$word, line $$line:\n$$_]; } }\
		$$last = $$word; $$first = 0; } \
		$$prevline = $$_;' < $< | $(PAGER)

# for slide shows

%.thumb.jpg : %.jpg
		-@cp $< $(patsubst %.jpg,%.thumb.jpg,$<)
		mogrify -resize 100x100 $(patsubst %.jpg,%.thumb.jpg,$<)

%.rotate.jpg : %.jpg
		mogrify -rotate 180  $< 

%.sshow.jpg : %.jpg
		-@make $*.rotate.jpg;
		-@make $*.thumb.jpg;

# audio stuff

# thank you linux HOWTOs
%.wav : %.mp3
		mpg123 --rate 44100 --stereo --buffer 3072 --resync -w $@ $<

%.wav : %.au
		sox $< $@

%.ogg : %.wav
		oggenc -o $@ $<

%.ogg : %.mp3
		mpg123 --rate 44100 --stereo --buffer 3072 --resync -w - $< | \
		oggenc -o $@ -

%.wavs : %.m3u
		-mkdir $*WAVS;
		cat $< | perl -ne 'BEGIN { $$trk = 0; } \
		chomp;if (not /^#/) { $$trk++; $$blnk = $$_; \
		$$blnk =~ s/\(/\\(/g; $$blnk =~ s/\)/\\)/g; \
		$$tf = $$blnk; $$tf =~ s/.mp3/.wav/; \
		$$tf =~ s/.*\/([^\/]+)$$/$$1/; $$tf =~ s/ /_/g; \
		$$tnm = sprintf(q[%3d],$$trk); $$tnm =~ s/ /0/g; \
		$$tf = qq[trk_$${tnm}_$${tf}]; \
		$$tf = qq[trk_$${tnm}_$${tf}]; \
		system qq[mpg123 -r 44100 --stereo -b 3072 -y -w $*WAVS/$$tf $$blnk;]; }' 


# put into a tar file
%.tar : %.manif
		@ $(TAR) -cvf $@ -T $< $< Makefile

%.tgz : %.manif
		@ $(TAR) -czvf $@ -T $< $< Makefile


############# CLEAN UP ##############

# clean up
%.clean : 
		-rm -f $*.aux $*.log $*.dvi $*.bbl $*.blg $*.toc $*.ilg $*.ind
		-rm -f $*.out $*.idx $*.lot $*.lof $*.brf
%.realclean : %.clean
		-rm -f $*.ps $*.pdf

############### RULES ###############

# an easy target
all: $(SUB_FIGS) $(DVI_TARGET) $(BBLS)  

dvis: $(patsubst %,%.dvi,$(PROJECTS))
pss: $(patsubst %,%.ps,$(PROJECTS))
pdfs: $(patsubst %,%.pdf,$(PROJECTS))

# targets
$(DVI_TARGET): $(TEX_SOURCE) $(TEX_EXTRAS) $(STY_FILES)
#$(PDF_TARGET): $(SUB_FIGS) $(PS_TARGET)

# see them 
presee: $(PROJECT).xsee
xsee: $(PROJECT).xsee
see: $(PROJECT).psee
pdfsee: $(PROJECT).pdfsee

# check it

spell: $(PROJECT).spell 

# damn bibunits?
bibunits: bu1.bbl bu2.bbl bu3.bbl bu4.bbl

# print
print: pss
		qup $(PS_TARGET)

# clean up
clean: $(patsubst %,%.clean,$(PROJECTS))
realclean: $(patsubst %,%.realclean,$(PROJECTS))

mailfinal: pss
		mpack -s wsel $(PS_TARGET) $(EMAIL)

######################
######################

safe: realclean
		@$(TAR) -cvf $(TAR_FILE) $(TAR_MAJORS) $(SUB_FIGS)
		@ -$(TAR) -uvf $(TAR_FILE) $(DVI_TARGET)
		@ -$(TAR) -uvf $(TAR_FILE) $(TAR_MINORS) 
		@ -$(TAR) -tvf $(TAR_FILE)
		rm -f $(TARBZ_FILE)
		bzip2 $(TAR_FILE)

remotesafe:
		@$(TAR) -cvf $(TAR_FILE) $(TAR_MAJORS)
		@ -$(TAR) -uvf $(TAR_FILE) $(TAR_MINORS) 
		rm -f $(TARBZ_FILE)
		bzip2 $(TAR_FILE)

# start a project w/ cvs
startcvs:
		cd .. ; $(CVS) -d$(CVS_ROOT_DIR) import -m 'adding $(PROJECT)' \
			$(CVS_STUB)/$(PROJECT) tag_$(PROJECT) $(OWNER) zero;
		$(CVS) -d$(CVS_ROOT_DIR) checkout CVSROOT/modules;
		echo $(PROJECT)  $(CVS_STUB)/$(PROJECT) >> CVSROOT/modules;
		$(CVS) -d$(CVS_ROOT_DIR) commit -m 'adding $(PROJECT)' CVSROOT;
		$(CVS) -d$(CVS_ROOT_DIR) release -d CVSROOT;
		cd ../../ ; rm -rf $(PROJECT) ; $(CVS) -d$(CVS_ROOT_DIR) \
			checkout $(PROJECT);

# add via cvs
add:
		@rm -if $(PS_TARGET)
		$(CVS) -d$(CVS_ROOT_DIR) add $(CVS_MAJORS)

# commit via cvs
commit:
		@rm -if $(PS_TARGET)
		$(CVS) -d$(CVS_ROOT_DIR) commit -m 'updating $(PROJECT)' 

# release via cvs
release:
		cd ../../ ; $(CVS) -d$(CVS_ROOT_DIR) release -d $(PROJECT)

# uncommit via cvs
unCVS:
		mkdir scratch
		-mv * scratch
		mv scratch/CVS .
		$(CVS) -d$(CVS_ROOT_DIR) remove
		$(CVS) -d$(CVS_ROOT_DIR) commit -m 'committing to remove $(PROJECT)' 
		mv scratch/* .
		rm -rf scratch

######################
######################
# remake the makefile from the stub directory
makefile:
		cp Makefile .oldMakefile;
		perl -pe 's/\s+PROJ_NAME/ $(PROJECT)/;s{\s+CWD_NAME}{ $(DIR_STUB)};' < $(STUB)/Makefile > Makefile;
		-diff Makefile .oldMakefile;
		-rm -i .oldMakefile;
	
