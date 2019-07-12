%Class for thesis writing

\NeedsTeXFormat{LaTeX2e}[1999/12/01]
\ProvidesClass{myclass}
\LoadClass[a4paper,oneside,12pt]{book}

\usepackage[utf8]{inputenc}
\RequirePackage{etoolbox}
\DeclareOption{a5paper}{\PassOptionsToClass{\CurrentOption}{book}}
\DeclareOption{twoside}{\OptionNotUsed}
\newtoggle{viet}
\DeclareOption{vietnamese}{\toggletrue{viet}}
\newtoggle{nolof}
\DeclareOption{nolof}{\toggletrue{nolof}}
\newtoggle{nolot}
\DeclareOption{nolot}{\toggletrue{nolot}}
\newtoggle{noacr}
\DeclareOption{noacr}{\toggletrue{noacr}}
\newtoggle{noack}
\DeclareOption{noack}{\toggletrue{noack}}
\DeclareOption*{\ClassWarning{myclass}{Unknown option '\CurrentOption'}}

\ProcessOptions\relax

\usepackage[left = 2 cm, right = 2 cm, bottom = 2.8 cm, top = 2.8 cm]{geometry}
\usepackage[final]{graphicx}
\usepackage[export]{adjustbox} %change image alignment
\usepackage{multirow} %multirow table
\usepackage{wrapfig} %Wrapping fig and table with text
\usepackage{float}
\usepackage{amsmath,amsfonts,amssymb}
\usepackage[nottoc]{tocbibind} %Add List of Figures, Bilbiography,... to the Table of Content
\usepackage{notoccite} %suppress citation from appearing in ToC
\usepackage{url} %typesetting url
\urlstyle{sf}

%packages for caption
\usepackage{caption}
\usepackage{subcaption}
\captionsetup{
	labelfont=bf,
	format=hang,
	width=0.8\linewidth
}

\usepackage{xcolor} %define new color
\definecolor{cstblue}{HTML}{3333FF}
\definecolor{cstred}{HTML}{FF1A1A}
\definecolor{cstgreen}{HTML}{00B300}
\definecolor{maroon}{HTML}{800000}
\definecolor{lstblue}{HTML}{00008B}
\definecolor{lstgray}{HTML}{404040}

\usepackage[colorlinks,linktocpage=true]{hyperref}
\hypersetup{%
	colorlinks=true,
	linkcolor=cstred,
	filecolor=black,
	urlcolor=cstblue,
	citecolor=cstgreen%
	}
	
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[R]{\slshape\nouppercase{\rightmark}}
\fancyhead[L]{\slshape\nouppercase{\leftmark}}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.5pt}

\usepackage{listings} %for source code input

\lstdefinelanguage{JavaScript}{
  keywords={typeof, new, true, false, catch, function, return, null, catch, switch, var, if, in, while, do, else, case, break},
  ndkeywords={class, export, boolean, throw, implements, import, this},
  ndkeywordstyle=\bfseries,
  sensitive=false,
  comment=[l]{//},
  morecomment=[s]{/*}{*/},
  morestring=[b]',
  morestring=[b]"
}

%these values can be overriden locally using \begin{lstlisting}[key=value]
\lstset{%
basicstyle=\ttfamily\small,
keywordstyle=\color{lstblue}\bfseries,
identifierstyle=,
commentstyle=\color{lstgray},
stringstyle=\color{maroon},
numbers=left,
breakatwhitespace=false,
breaklines=true,
captionpos=b,
keepspaces=true,
showspaces=false,
showtabs=false,
showstringspaces=false,
escapeinside={(*}{*)}%
}

\iftoggle{viet}{%
	\usepackage[vietnamese=nohyphenation]{hyphsubst} %Avoid hyphenation warning
	\usepackage[vietnamese]{babel}%
	\usepackage{indentfirst}
}{%
	\usepackage[english]{babel}%
}


\usepackage{csquotes}
\usepackage[backend=biber,
	    style=numeric,
	    sorting=ynt,
	    citestyle=numeric-comp
	    ]{biblatex}
	    
\defbibheading{bibliography}[\refname]{%
\iftoggle{viet}{\renewcommand{#1}{Tài liệu tham khảo}}{~}
\chapter*{\centering\LARGE #1}%
\markboth{#1}{#1}
\addcontentsline{toc}{chapter}{#1}}
	    
\usepackage{tocloft} %Configure Table of Content, List of Figures, List of Tables

\renewcommand{\cfttoctitlefont}{\hfill\LARGE\bfseries}
\renewcommand{\cftaftertoctitle}{\hfill}
\renewcommand{\cftloftitlefont}{\hfill\LARGE\bfseries}
\renewcommand{\cftafterloftitle}{\hfill}
\renewcommand{\cftlottitlefont}{\hfill\LARGE\bfseries}
\renewcommand{\cftafterlottitle}{\hfill}

\renewcommand\cftfigpresnum{\figurename~}
\newlength\mylength
\settowidth\mylength{\cftfigpresnum}
\addtolength\cftfignumwidth{\mylength}

\renewcommand\cfttabpresnum{\tablename~}
\settowidth\mylength{\cfttabpresnum}
\addtolength\cfttabnumwidth{\mylength}

\usepackage[acronym]{glossaries} %Construct List of Acronyms 


\iftoggle{viet}{%
	\def\theGroup{ĐẠI HỌC QUỐC GIA TP.HỒ CHÍ MINH}%
}{%
	\def\theGroup{Vietnam National University - Ho Chi Minh City}%
}

\iftoggle{viet}{%
	\def\theUniversity{TRƯỜNG ĐẠI HỌC KHOA HỌC TỰ NHIÊN}%
}{%
	\def\theUniversity{University of Science}%
}

\iftoggle{viet}{%
	\def\theFaculty{KHOA VẬT LÝ - VẬT LÝ KỸ THUẬT}%
}{%
	\def\theFaculty{Faculty of Physics and Engineering Physics}%
}

\iftoggle{viet}{%
	\def\theDepartment{CHUYÊN NGÀNH VẬT LÝ TIN HỌC}%
}{%
	\def\theDepartment{Department of Physics and Computer Science}%
}

\iftoggle{viet}{%
	\def\theReport{KHOÁ LUẬN TỐT NGHIỆP}%
}{%
	\def\theReport{BACHELOR THESIS}%
}

\def\thesupervisorName{}
\newcommand{\supervisorName}[1]{\def\thesupervisorName{#1}}
\newcommand{\printcoverpage}{%
	\thispagestyle{empty}
	\begin{center}
	{\bfseries\parskip=0pt
	
	\theGroup
	\vspace*{0.1cm}
	
	\theUniversity
	\vspace*{0.1cm}
	
	\theFaculty
	\vspace*{0.1cm}
	
	\theDepartment\\
	\vspace*{0.1cm}
	------------------oOo------------------
	}
	\vspace*{1cm}
	
	{\bfseries
	\large
	\theReport}
	
	\vspace*{3cm plus 1cm}
	
	\iftoggle{viet}
	{\begin{flushleft}
		\textsl{\underline{Đề tài:}}
	\end{flushleft}}
	{~}
	{\huge\bfseries
		\@title\par	
	}
	\end{center}
	
	\vspace*{3cm plus 1 cm}
	
	\hfill
	{\bfseries\large
	\iftoggle{viet}{%
	\begin{tabular}{r l}
	\underline{SVTH}: & \@author\\
	\underline{CBHD}: & \thesupervisorName\\
	\end{tabular}%
	}{%
	\begin{tabular}{r l}
	\underline{Student}: & \@author\\
	\underline{Supervisor}: & \thesupervisorName\\
	\end{tabular}%
	}}
	
	\vfill
	\begin{center}
	-----------------------------------------\\
	\bfseries
	\iftoggle{viet}%
	{TP. HỒ CHÍ MINH - \the\year}%
	{HO CHI MINH CITY - \the\year}%
	\end{center}
	\clearpage
}

\def\theacknowledgements{}
\newcommand{\acknowledgements}[1]{\def\theacknowledgements{#1}}

\newcommand{\printfrontmatter}{%
		\frontmatter
		\iftoggle{noack}{~}{%
		\iftoggle{viet}{%
		\chapter*{\centering\LARGE Lời cảm ơn}%
		}{%
		\chapter*{\centering\LARGE Acknowledgements}%
		}
		\label{ch:ack}
		\thispagestyle{plain}
		\theacknowledgements
		
		\newpage
		}
		\tableofcontents
		
		\iftoggle{noacr}{~}
		{%
		\newpage
		\iftoggle{viet}{%
		\printnoidxglossary[type=\acronymtype,title={\protect\centering\Large Các kí hiệu viết tắt}, style=listdotted]
		\addcontentsline{toc}{chapter}{Các kí hiệu viết tắt}
		}{%
		\printnoidxglossary[type=\acronymtype,title={\protect\centering\Large Acronym}, style=listdotted]
		\addcontentsline{toc}{chapter}{Acronym}
		}%
		}
		
		\iftoggle{nolof}{~}
		{%
		\newpage
		\listoffigures
		}
		
		\iftoggle{nolot}{~}
		{%
		\newpage
		\listoftables
		}
		
		\mainmatter
}

\setlength{\parskip}{0.5cm}
\setlength{\parindent}{0.5cm}

\newcommand{\startintroduction}{%
		\iftoggle{viet}{%
		\chapter*{\centering\LARGE Lời giới thiệu}
		\addcontentsline{toc}{chapter}{Lời giới thiệu}%
		}{%
		\chapter*{\centering\LARGE Introduction}
		\addcontentsline{toc}{chapter}{Introduction}%
		}\label{ch:intro}
}

%%
%% End of File 'myclass.cls'