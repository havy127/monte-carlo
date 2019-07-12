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
\newtoggle{noacr}
\DeclareOption{noacr}{\toggletrue{noacr}}
\newtoggle{print}
\DeclareOption{print}{\toggletrue{print}}
\DeclareOption*{\ClassWarning{myclass}{Unknown option '\CurrentOption'}}
\ProcessOptions\relax

\iftoggle{viet}{%
	\usepackage[vietnamese=nohyphenation]{hyphsubst} %Avoid hyphenation warning
	\usepackage[vietnamese]{babel}%
	\usepackage{indentfirst}
}{%
	\usepackage[english]{babel}%
}
\usepackage[left = 1.38 in, right = 0.79 in, bottom = 1.15 in, top = 1.15 in]{geometry}
\usepackage[final]{graphicx}
\usepackage{multirow} %multirow table
\usepackage{wrapfig} %Wrapping fig and table with text
\usepackage{array} %extra options for table spec argument
\usepackage{float}
\usepackage{amsmath,amsfonts,amssymb}
\usepackage[nottoc]{tocbibind} %Add List of Figures, List of Tables,... to the Table of Content
\usepackage{notoccite} %suppress citation from appearing in ToC
\usepackage{path} %for file names, directory path, email... typesetting
\usepackage[table,figure]{totalcount} 

\usepackage[obeyspaces]{url} %typesetting url
\urlstyle{sf}

\usepackage{titlesec} %customize chapter and sections'... title
\titleformat{name=\chapter,numberless}[display]
{\normalfont\filcenter\bfseries}{}{0pt}{\LARGE}

\titlespacing*{name=\chapter,numberless}{0pt}{-55pt}{40pt}

\titleformat{name=\chapter}[display]
{\normalfont\Large\filcenter\bfseries}
{\titlerule[1pt]%
\vspace{1pt}%
\titlerule
\vspace{1pc}%
\LARGE\MakeUppercase{\chaptertitlename} \thechapter}
{1pc}
{\titlerule
\vspace{1pc}%
\LARGE}

\titlespacing*{name=\chapter}{0pt}{-40pt}{40pt}

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
\definecolor{lstgray}{HTML}{787878}

\usepackage[linktocpage=true]{hyperref}
\hypersetup{%
	colorlinks=true,
	breaklinks=true,
	linkcolor=cstred,
	filecolor=black,
	urlcolor=cstblue,
	citecolor=cstgreen%
	}
\iftoggle{print}{\hypersetup{%
	colorlinks=false,
	breaklinks=true}}
	
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[R]{\slshape\nouppercase{\rightmark}}
\fancyhead[L]{\slshape\nouppercase{\chaptertitlename\ \thechapter}}
\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0.5pt}

\fancypagestyle{supplement}{%
    \fancyhf{}
    \fancyhead[L]{\slshape\nouppercase{\leftmark}}
    \fancyhead[R]{\thepage}
    \renewcommand{\headrulewidth}{0.5pt}%
    }
    
    
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
frame=,
commentstyle=\color{lstgray},
stringstyle=\color{maroon},
numbers=left,
breakatwhitespace=false,
breaklines=true,
captionpos=b,
keepspaces=true,
showspaces=false,
showtabs=false,
escapeinside=,
showstringspaces=false,
inputencoding=utf8%
}

\usepackage{csquotes}
\usepackage[backend=biber,
	    style=numeric,
	    sorting=ynt,
	    citestyle=numeric-comp
	    ]{biblatex}
	    
\defbibheading{bibliography}[\refname]{%
\iftoggle{viet}{\renewcommand{#1}{Tài liệu tham khảo}}{~}
\chapter*{#1}
\markboth{#1}{#1}
\addcontentsline{toc}{chapter}{#1}}
	    
\usepackage{tocloft} %Configure Table of Content, List of Figures, List of Tables

\renewcommand{\cfttoctitlefont}{\hfill\LARGE\bfseries}
\renewcommand{\cftaftertoctitle}{\hfill}
\renewcommand{\cftloftitlefont}{\hfill\LARGE\bfseries}
\renewcommand{\cftafterloftitle}{\hfill}
\renewcommand{\cftlottitlefont}{\hfill\LARGE\bfseries}
\renewcommand{\cftafterlottitle}{\hfill}

\setlength{\cftbeforeloftitleskip}{-15pt}
\setlength{\cftbeforelottitleskip}{-15pt}
\setlength{\cftbeforetoctitleskip}{-15pt}

\renewcommand\cftfigpresnum{\bfseries\figurename~}
\newlength\mylength
\settowidth\mylength{\cftfigpresnum}
\addtolength\cftfignumwidth{\mylength}

\renewcommand\cfttabpresnum{\bfseries\tablename~}
\settowidth\mylength{\cfttabpresnum}
\addtolength\cfttabnumwidth{\mylength}

\usepackage[acronym]{glossaries} %Construct List of Acronyms 

\setlength{\parskip}{0.2cm}
\setlength{\parindent}{0.39in}
\linespread{1.3}

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

\iftoggle{viet}{%
	\def\thePlace{TP. HỒ CHÍ MINH}%
}{%
	\def\thePlace{HO CHI MINH CITY}%
}

\def\theDate{\the\year}

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
	
	\vspace*{2cm plus 1cm minus 0.5cm}
	
	\iftoggle{viet}
	{\begin{flushleft}
		\textsl{\Large\underline{Đề tài:}}
	\end{flushleft}}
	{~}
	{\huge\bfseries
		\@title\par	
	}
	\end{center}
	
	\vspace*{2cm plus 1 cm minus 0.5cm}
	
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
	\thePlace\ - \theDate
	\end{center}
	\clearpage
}

\def\theacknowledgements{}
\newcommand{\acknowledgements}[1]{\def\theacknowledgements{#1}}

\iftoggle{viet}{%
	\def\acrtitle{Các kí hiệu viết tắt}%
}{%
	\def\acrtitle{Acronym}%
}

\newcommand\conditionalLoF{\iftotalfigures\listoffigures\fi}
\newcommand\conditionalLoT{\iftotaltables\listoftables\fi}

\newcommand{\printfrontmatter}{%
		\frontmatter
		\ifdefempty{\theacknowledgements}{~}{
		\pagenumbering{gobble}
		\iftoggle{viet}{%
		\chapter*{Lời cảm ơn}\markboth{}{Lời cảm ơn}%
		}{%
		\chapter*{Acknowledgements}\markboth{}{Acknowledgements}%
		}
		\label{ch:ack}
		\theacknowledgements
		\clearpage
		\pagenumbering{roman}
		}
		
		{
		\tableofcontents
		}
		
		\iftoggle{noacr}{~}
		{%
		\clearpage
		
		\printnoidxglossary[type=\acronymtype,title=\acrtitle, style=listdotted]
		\addcontentsline{toc}{chapter}{\acrtitle}
		}
		
		{%
		\clearpage
		\conditionalLoF
		}
		
		{%
		\clearpage
		\conditionalLoT
		}
		
		\mainmatter
}

\newcommand{\startintroduction}{%
		\iftoggle{viet}{
		\chapter*{Lời giới thiệu}
		\markboth{}{Lời giới thiệu}
		\addcontentsline{toc}{chapter}{Lời giới thiệu}%
		}{%
		\chapter*{Introduction}
		\markboth{}{Introduction}
		\addcontentsline{toc}{chapter}{Introduction}%
		}\label{ch:intro}
}


\newcommand{\thebackmatter}{%
	\titleformat{\chapter}[display]
	{\normalfont\filcenter\Large\bfseries}{\chaptertitlename\ \thechapter}{20pt}{\LARGE}
	\titlespacing*{\chapter}{0pt}{-50pt}{40pt}
	\appendix 
	\pagestyle{supplement}
}

\newcommand{\reference}{
		\backmatter
		\printbibliography
		\thispagestyle{plain}
}

%%
%% End of File 'myclass.cls'