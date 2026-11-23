#!/bin/sh
# Build main_highlighted.tex/pdf: additions vs the submitted version in red, deletions hidden.
# Methods is diffed against the submitted Methods moved to its new place (submitted/main_reordered.tex).
set -e
cd "$(dirname "$0")/.."
S=submitted/main.tex
{ sed -n 1,110p $S; sed -n 364,477p $S; sed -n 111,363p $S; sed -n '478,$p' $S; } > submitted/main_reordered.tex
perl tools/latexdiff-so --flatten --math-markup=coarse submitted/main_reordered.tex main.tex > main_highlighted.tex 2>/dev/null
f=main_highlighted.tex
perl -0pi -e '
 s/\\providecommand\{\\DIFaddtex\}\[1\]\{\{\\protect\\color\{blue\}\\uwave\{#1\}\}\}/\\providecommand{\\DIFaddtex}[1]{{\\protect\\color{red}#1}}/;
 s/\\providecommand\{\\DIFdeltex\}\[1\]\{\{\\protect\\color\{red\}\\sout\{#1\}\}\}/\\providecommand{\\DIFdeltex}[1]{}/;
 s/\\providecommand\{\\DIFaddincludegraphics\}\[2\]\[\]\{\{\\color\{blue\}\\fbox\{\\DIFOincludegraphics\[#1\]\{#2\}\}\}\}/\\providecommand{\\DIFaddincludegraphics}[2][]{\\DIFOincludegraphics[#1]{#2}}/;
 s/(\\providecommand\{\\DIFdelincludegraphics\}\[2\]\[\]\{%)/\\providecommand{\\DIFdelincludegraphics}[2][]{}\n\\providecommand{\\DIFdelincludegraphicsUNUSED}[2][]{%/;
 s/\\DIFaddbegin \\begin\{aligned\}/\\DIFaddbegin \\color{red}\\begin{aligned}/g;
 s/\\(sub)*(section|paragraph)\{\\DIFdel\{[^\n]*\}\}\n%DIFAUXCMD\n\\addtocounter\{(sub)*(section|paragraph)\}\{-1\}%DIFAUXCMD\n//g;
 s/\\begin\{(displaymath|align\*)\}((?:(?!\\end\{\1\}).)*?)\\end\{\1\}%DIFAUXCMD/ ($2 =~ m|DIFadd|) ? $& : "%DIF removed deleted math"/gse;
' $f
latexmk -pdf -interaction=nonstopmode $f > /dev/null 2>&1 || true
grep -E '^!|Output written' main_highlighted.log
