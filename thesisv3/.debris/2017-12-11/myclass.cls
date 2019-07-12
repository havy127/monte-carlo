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
\DeclareOption*{\ClassWarning{myclass}{Unknown option '\CurrentOption'}}


\ProcessOptions\relax

\usepackage[left = 2 cm, right = 2 cm, bottom = 2.8 cm, top = 2.8 cm]{geometry}
\usepackage[final]{graphicx}
\usepackage{amsmath,amsfonts,amssymb}
\usepackage{url} %typesetting url
\urlstyle{sf}

%packages for caption
\usepackage{caption}
\usepackage{subcaption}

\usepackage{xcolor} %define new color
\definecolor{cstblue}{HTML}{3333FF}
\definecolor{cstred}{HTML}{FF1A1A}
\definecolor{cstgreen}{HTML}{00B300}
\definecolor{maroon}{HTML}{800000}
\definecolor{lstblue}{HTML}{00008B}
\definecolor{lstgray}{HTML}{404040}

\usepackage[colorlinks]{hyperref}
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


\lstset{%
basicstyle=\ttfamily\footnotesize,
keywordstyle=\color{lstblue}\bfseries,
identifierstyle=,
commentstyle=\color{lstgray},
stringstyle=\color{maroon},
numbers=left,
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

\def\theauthorName{}
\newcommand{\authorName}[1]{\def\theauthorName{#1}}
\def\thesupervisorName{}
\newcommand{\supervisorName}[1]{\def\thesupervisorName{#1}}
\newcommand{\printcoverpage}{%
	\thispagestyle{empty}
	\begin{center}
	{\bfseries
	
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
	\vspace*{2cm}
	
	\iftoggle{viet}%
	{\Large\bfseries KHOÁ LUẬN TỐT NGHIỆP}%
	{\Large\bfseries BACHELOR THESIS}
	
	\vspace*{4cm}
	
	{\huge\bfseries
		\@title	
	}
	\end{center}
	
	\vspace*{4cm}
	\hfill
	
	\iftoggle{viet}{%
	\begin{tabular}{r l}
	\bfseries\underline{SVTH}: & \bfseries\theauthorName\\
	\bfseries\underline{CBHD}: & \bfseries\thesupervisorName\\
	\end{tabular}%
	}{%
	\begin{tabular}{r l}
	\bfseries\underline{Student}: & \bfseries\theauthorName\\
	\bfseries\underline{Supervisor}: & \bfseries\thesupervisorName\\
	\end{tabular}%
	}
	
	\vfill
	\begin{center}
	-----------------------------------------\\
	\iftoggle{viet}%
	{\bfseries TP. HỒ CHÍ MINH - \the\year}%
	{\bfseries HO CHI MINH CITY - \the\year}%
	\end{center}
}