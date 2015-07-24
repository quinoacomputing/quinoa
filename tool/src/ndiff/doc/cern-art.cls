%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CERN LHC/MTA Generic Article Class
% Version 1.0 (19/01/1998)
% Created by Laurent Deniau, CERN (Laurent.Deniau@cern.ch)
%
% Last Modification: Tue Jan 20 15:10:45 CET 1998
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Provided fields                 Default
% ----------------------------------------------
% \Author{L. Deniau}               Nobody
% \Institute{CERN LHC/MTA}         empty
% \DocReference{MTA--IN--98--51}   Draft
% \Title{How to use this class}    Notitle
% \Date{January 19, 1998}          \today
% \RevisionDate{\today}            empty
% \Email{Laurent.Deniau@cern.ch}   empty
% \Keywords{\LaTeX class}          empty
% \Distribution{Me, You}           empty
% \Participants{Him, Them}         empty
% \PSLogo{/tex/logos/mta-logo.ps}  empty
% \LogoHeight{2cm}                 1.75cm
%
% To build the title:
% \Maketitle
%
% To write a summary
% \Summary{This is the summary}
%
% To write an object header
% \Object{Object of the paper}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------- Identification Part --------------

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cern-art}[1998/01/01 CERN article class)]

% ---- Declarations and Execution of Options -----

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

% ------- Package and Class Loading Part  -------

\LoadClass[12pt,a4paper]{article} 

\RequirePackage{calc}
\RequirePackage{float}
\RequirePackage{ifthen}
\RequirePackage{graphicx}
\RequirePackage[english]{babel}

% -------------- Main Code Part -----------------

% Add 10% to the text height
\textheight  1.1\textheight
\topmargin -0.05\textheight

% Create a counter with 2 digits
%\newcounter{Year}\setcounter{Year}{\year-1900}

% Fancy line
\newcommand{\HRule}{%
  \par\noindent
  \parbox{\textwidth}{
    \setlength{\baselineskip}{1mm}
    \rule{\textwidth}{0.25mm}
    \rule{\textwidth}{0.5mm}
    }}

% List with labelwidth ajusted to the width of the argument
\newenvironment{Ventry}[1]{%
  \begin{list}{}{\renewcommand{\makelabel}[1]{##1:\hfil}
      \settowidth{\labelwidth}{#1:}
      \setlength{\leftmargin}{\labelwidth+\labelsep}}}
  {\end{list}}

% Title definition
\newcommand{\Maketitle}{%
  \parskip 0.25cm
  \noindent
  \begin{minipage}[t][\@LogoHeight]{0.5\textwidth}
	\vspace{-5pt}
    \IfFileExists{\@PSLogo}{%
      \includegraphics[height=\@LogoHeight,keepaspectratio=true]{\@PSLogo}}{}
  \end{minipage}
  \hfill
  \begin{minipage}[t][\@LogoHeight]{0.5\textwidth}
	\vspace{0pt}
    \begin{flushright}
      \textbf{\@DocReference} \\
      \@Date \\
      \ifthenelse{\equal{\@RevisionDate}{}}{}{revised \@RevisionDate}
      \ifthenelse{\equal{\@Email}{}}{}{\@Email}
    \end{flushright}
  \end{minipage}
  \vspace{4ex}
  \par\noindent{\LARGE\textbf{\textsf{\@Title}}}
  \par\noindent{\large\@Author}
  \par\noindent{\@Institute}
  \vspace{1ex}
  { \ifthenelse{\equal{\@Keywords\@Participants\@Distribution}{}}{}{
    \footnotesize\begin{Ventry}{Distributionxxx}
    \ifthenelse{\equal{\@Keywords}    {}}{}{\item[Keywords]     \@Keywords}
    \ifthenelse{\equal{\@Participants}{}}{}{\item[Participants] \@Participants}
    \ifthenelse{\equal{\@Distribution}{}}{}{\item[Distribution] \@Distribution}
    \end{Ventry}}}
  \HRule}

% Summary definition
\newcommand{\Summary}[1]{%
  \parskip 0.25cm
  {\em\par\noindent\textbf{Summary}\par\noindent #1}
  \HRule}

% Object definition
\newcommand{\Object}[1]{%
  \newlength{\objectlen}
  \settowidth{\objectlen}{\textsc{Object:}\hspace*{0.5cm}}
  \parskip 0.25cm
  {\par\noindent\textsc{Object:}\hspace*{0.5cm}
    \parbox[t]{.95\linewidth-\objectlen}{#1}}
  \vspace*{0.5cm}}

% Section definitions (come from article.cls)

\renewcommand{\part}{\secdef\Part\sPart}

\newcommand{\Part}[2][?]{%
  \refstepcounter{part}%
  \addcontentsline{toc}{part}{\textsf{\thepart\hspace{1em}#1}}%
  {\flushleft\normalfont\Large\bfseries\textsf{\partname~\thepart}
    \nobreak
    \flushleft\normalfont\huge\bfseries\textsf{#2}
    \nobreak}%
  \sectionmark{#1}\vspace{\baselineskip}}

\newcommand{\sPart}[1]{%
  {\flushleft\normalfont\huge\bfseries\textsf{#1}
    \nobreak}%
  \vspace{\baselineskip}}

\renewcommand{\section}{\@startsection
  {section}
  {1}
  {\z@}
  {-3.5ex \@plus -1ex \@minus -.2ex}
  {2.3ex \@plus.2ex}
  {\normalfont\Large\bfseries\textsf}}

\renewcommand{\subsection}{\@startsection
  {subsection}
  {2}
  {\z@}
  {-3.25ex\@plus -1ex \@minus -.2ex}
  {1.5ex \@plus .2ex}
  {\normalfont\large\bfseries\textsf}}

\renewcommand{\subsubsection}{\@startsection
  {subsubsection}
  {3}
  {\z@}
  {-3.25ex\@plus -1ex \@minus -.2ex}
  {1.5ex \@plus .2ex}
  {\normalfont\normalsize\bfseries\textsf}}

\renewcommand{\paragraph}{\@startsection
  {paragraph}
  {4}
  {\z@}
  {3.25ex \@plus1ex \@minus.2ex}
  {-1em}
  {\normalfont\normalsize\bfseries\textsf}}

\renewcommand{\subparagraph}{\@startsection
  {subparagraph}
  {5}
  {\parindent}
  {3.25ex \@plus1ex \@minus .2ex}
  {-1em}
  {\normalfont\normalsize\bfseries\textsf}}

% Floats definitions

\setcounter{topnumber}{5}
\setcounter{bottomnumber}{5}
\setcounter{totalnumber}{10}
\renewcommand{\topfraction}{1}
\renewcommand{\bottomfraction}{1}
\renewcommand{\textfraction}{0}
\renewcommand{\floatpagefraction}{1}

% Figures
\newfloat{fig}{hbt}{lof}[section]
\floatname{fig}{\normalfont\textsc{Fig}}
% Tables
\newfloat{tab}{hbt}{lot}[section]
\floatname{tab}{\normalfont\textsc{Tab}}
% Algorithm
\floatstyle{boxed}
\newfloat{algo}{hbt}{loa}[section]
\floatname{algo}{\normalfont\textsc{Algorithm}}
% Program
\floatstyle{ruled}
\newfloat{prog}{hbt}{lop}[section]
\floatname{prog}{\normalfont\textsc{Program}}
% Script
\newfloat{script}{hbt}{lop}[section]
\floatname{script}{\normalfont\textsc{Script}}
% default floatstyle
\floatstyle{plain}

% Variables
\def\Author#1{\gdef\@Author{#1}}
\gdef\@Author{Nobody}
\def\DocReference#1{\gdef\@DocReference{#1}}
\gdef\@DocReference{DRAFT}
\def\Institute#1{\gdef\@Institute{#1}}
\gdef\@Institute{}

\def\Title#1{\gdef\@Title{#1}}
\gdef\@Title{Notitle}
\def\Date#1{\gdef\@Date{#1}}
\gdef\@Date{\today}
\def\RevisionDate#1{\gdef\@RevisionDate{#1}}
\gdef\@RevisionDate{}

\def\Email#1{\gdef\@Email{#1}}
\gdef\@Email{}

\def\Keywords#1{\gdef\@Keywords{#1}}
\gdef\@Keywords{}
\def\Participants#1{\gdef\@Participants{#1}}
\gdef\@Participants{}
\def\Distribution#1{\gdef\@Distribution{#1}}
\gdef\@Distribution{}

\def\PSLogo#1{\gdef\@PSLogo{#1}}
\gdef\@PSLogo{}
\def\LogoHeight#1{\gdef\@LogoHeight{#1}}
\gdef\@LogoHeight{1.75cm}

% PDF flag

%\newif\ifpdf
%\ifx\pdfoutput\undefined
%  \pdffalse
%\else
%  \pdfoutput=1
%  \pdftrue
%\fi

%%% Local Variables: 
%%% mode: latex
%%% TeX-master: t
%%% End: 
