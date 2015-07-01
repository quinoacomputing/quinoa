\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{isisarticle}[2003/05/14 Version 0.01 Oliver Broeker]
\RequirePackage{graphicx}
\RequirePackage{multicol}
\RequirePackage{url}
\RequirePackage{xspace}

\newcommand{\colab}{\textsf{colab}\xspace}

\newcommand{\COLABaff}{%
  Sandia National Laboratories \\
  PO BOX 5800 Albuquerque, NM (USA) \\
  %%{\sl and} \\
  %%ETH Z\"{u}rich Computational Laboratory\\
  %%Hirschengraben 84\\
% ETH Zentrum\\
  %%8092 Z\"{u}rich, Switzerland
}
\newcommand{\defaultaffiliation}{\COLABaff}

\newcommand{\COLAB@colofon}{
    {\small
        {\bfseries 
          Sandia National Laboratories \\
          PO BOX 5800 Albuquerque, 
          New Mexico, USA \\
          %%{\sl and} \\
          %%ETH Z\"{u}rich Computational Laboratory\\
          %%Hirschengraben 84\\
          %%8092 Z\"{u}rich\\ 
          %%Switzerland
        }
      }
    }

\newif\if@email
\@emailfalse
\def\email#1{\@emailtrue\gdef\@email{\url{#1}}}

\newif\if@telephone
\@telephonefalse
\def\telephone#1{\@telephonetrue\gdef\@telephone{#1}}

\newif\if@http
\@httpfalse
\def\http#1{\@httptrue\gdef\@http{\url{#1}}}

\newif\if@correspondingauthorset
\@correspondingauthorsetfalse
\def\correspondingauthor#1{\@correspondingauthorsettrue\gdef\@correspondingauthor{#1}}

\newcommand{\COLAB@correspondingauthor}{
%%jhw\small
%jhw\textbf{Corresponding editor:}\\
%% jhw modified
\if@http{\@http}\fi
%%jhw \textbf{}\\
%%jhw \@correspondingauthor\\
%%jhw \if@telephone tel:\,\@telephone\\ \fi
%%jhw \if@email{\@email}\\\fi\if@http{\@http}\fi
}


\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions
\LoadClass[11pt,a4paper,twoside]{book}

%
%   Defining the page size
%
\ifcase \@ptsize
    % mods for 10 pt
    \oddsidemargin  0.15 in     %   Left margin on odd-numbered pages.
    \evensidemargin 0.35 in     %   Left margin on even-numbered pages.
    \marginparwidth 1 in        %   Width of marginal notes.
    \oddsidemargin 0.25 in      %   Note that \oddsidemargin = \evensidemargin
    \evensidemargin 0.25 in
    \marginparwidth 0.75 in
    \textwidth 5.875 in % Width of text line.
\or % mods for 11 pt
    \oddsidemargin 0.1 in      %   Left margin on odd-numbered pages.
    \evensidemargin 0.15 in    %   Left margin on even-numbered pages.
    \marginparwidth 1 in       %   Width of marginal notes.
    \oddsidemargin 0.125 in    %   Note that \oddsidemargin = \evensidemargin
    \evensidemargin 0.125 in
    \marginparwidth 0.75 in
    \textwidth 6.125 in % Width of text line.
\or % mods for 12 pt
    \oddsidemargin -10 pt      %   Left margin on odd-numbered pages.
    \evensidemargin 10 pt      %   Left margin on even-numbered pages.
    \marginparwidth 1 in       %   Width of marginal notes.
    \oddsidemargin 0 in      %   Note that \oddsidemargin = \evensidemargin
    \evensidemargin 0 in
    \marginparwidth 0.75 in
    \textwidth 6.375 true in % Width of text line.
\fi

\renewcommand{\textfraction}{0.2}
\renewcommand{\bottomfraction}{0.8}
\renewcommand{\topfraction}{0.8}
\renewcommand{\floatpagefraction}{0.8}

\newif\ifthanks
\newbox\@kapbox
\newcount\@kapcount

\newtoks\t@a
\newtoks\t@b
\long\def\append@item#1\to#2{\t@a={\@k@p{#1}}\t@b=\expandafter{#2}%
  \xdef#2{\the\t@b\the\t@a}}
\def\get@left#1\to#2{\expandafter\g@l#1\g@l#1#2}
\long\def\g@l\@k@p#1#2\g@l#3#4{\def#4{#1}\def#3{#2}}


\def\@authors{}
\def\@allauthors{}
\def\@institutes{}
\def\authors@institute{}

\def\author@#1{\append@item#1\to\@authors
  \append@item#1\to\@allauthors}
\def\institute@#1{\append@item#1\to\@institutes
  \expandafter\append@item\expandafter{\@authors}\to\authors@institute
  \gdef\@authors{}}

\let\author=\author@
\let\oneauthor=\author@
\let\institute=\institute@

\def\@formatname#1#2{\@kapcount\c@footnote\begingroup
  \def\thanks##1{\global\thankstrue}\setbox\@kapbox\vbox{#1}\endgroup
  \c@footnote\@kapcount
  \begingroup
  \ifthanks
  \def\fn##1\thanks##2{{##1}\hbox
    to 0pt{#2\hss}\thanks{##2}}\expandafter\fn#1
  \else
  \def\fn##1{{##1}#2}\expandafter\fn\expandafter{#1}\fi
  \endgroup\global\thanksfalse}


\def\@formatauthors{{\small\begingroup
    \newcount\cnt@authors
    \cnt@authors=0\def\@k@p##1{\advance\cnt@authors by 1}\@curauths
    \def\@k@p##1{\advance\cnt@authors by -1
      \ifnum\cnt@authors>1 \hbox{\@formatname{##1}{,}\penalty-1000}\space\fi
      \ifnum\cnt@authors=1 \hbox{\@formatname{##1}{} and\penalty-1000}\space\fi
      \ifnum\cnt@authors<1 \hbox{\@formatname{##1}{}}\par\fi}\@curauths
    \endgroup}}


\def\@formatinstitute{{\@curinst\par}}

\newif\if@runningtitleset
\newif\if@runningauthorset

\def\runningauthor#1{\if@runningauthorset\else
  \def\@runningauthor{{#1}}\fi
  \global\@runningauthorsettrue}
\gdef\@runningauthor{}


\def\@formatauthorsandinstitutes{\begingroup
  \newcount\cnt@institutes
  \newcount\cnt@authors
  \newbox\@authbox
  \cnt@authors=0\def\@k@p##1{\advance\cnt@authors by 1}\@allauthors
  \xdef\@authnames{}%
  \begingroup
  \def\thanks##1{}%
  \def\@k@p##1{\advance\cnt@authors by -1
    \ifnum\cnt@authors=0 \xdef\@authnames{\@authnames ##1}\fi
    \ifnum\cnt@authors=1 \xdef\@authnames{##1\ and }\fi
    \ifnum\cnt@authors>1 \xdef\@authnames{##1\ et al.}\cnt@authors=0\fi}%
  \global\setbox\@authbox=\hbox{\@allauthors}%
  \expandafter\runningauthor\expandafter{\@authnames}
  \endgroup
  \cnt@institutes=0
  \def\@k@p##1{\advance\cnt@institutes by 1}\@institutes
  \loop\ifnum\cnt@institutes>0
  \get@left\authors@institute\to\@curauths
  {\normalsize\sffamily\bfseries\noindent\@formatauthors}
  \get@left\@institutes\to\@curinst
  {\normalsize\sffamily\noindent\@formatinstitute}
  \ifnum\cnt@institutes=2
  \vskip 0.4\baselineskip
  %%      and   <<<==== this is ugly
  \fi
  \ifnum\cnt@institutes>1
  \vskip 0.4\baselineskip
  \fi
  \advance\cnt@institutes by -1
  \repeat\endgroup}

\newif\if@subtitle
\@subtitlefalse
\def\subtitle#1{\@subtitletrue\gdef\@subtitle{#1}}

\newif\if@affiliation
\@affiliationfalse
\def\affiliation#1{\@affiliationtrue\gdef\@affiliation{#1}}

\newif\if@status
\@statusfalse
\def\status#1#2{\@statustrue\gdef\@statushead{#1}\gdef\@statustail{#2}}

\newif\if@acknowledge
\@acknowledgefalse
\newcommand{\acknowledge}[1]{\@acknowledgetrue\gdef\@acknowledge{#1}}

\newif\if@COLABreport
\@COLABreportfalse
\newcommand{\COLABreport}[1]{\@COLABreporttrue\gdef\@COLABreport{#1}}

\newif\if@COLABreportdate
\@COLABreportdatefalse
\newcommand{\COLABreportdate}[1]{\@COLABreportdatetrue\gdef\@COLABreportdate{#1}}


\AtEndDocument{\COLAB@backcoverpage{\@COLABreport}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Headers and Footers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\ps@colabhead{%
   \let\@oddfoot\@empty
   \let\@evenfoot\@empty
   \def\@evenhead{%
       \parbox{\textwidth}{\sffamily\thepage\hfill\leftmark\\
       \rule[1.5ex]{\textwidth}{.1pt}}}
   \def\@oddhead{%
       \parbox{\textwidth}{\sffamily\rightmark\hfill\thepage\\
       \rule[1.5ex]{\textwidth}{.1pt}}}
   \let\@mkboth\markboth
   \def\subsectionmark##1{}
    \def\sectionmark##1{%
       \markboth{\runningauthor}{%
       \ifnum \c@secnumdepth >\m@ne
%jhw       JHWJHWSection \thesection\hskip 1em\relax
       \fi
       \projectname}}
}

%jhw   \def\sectionmark##1{%
%jhw       \markboth{\runningauthor}{%
%jhw       \ifnum \c@secnumdepth >\m@ne
%jhw           Section \thesection\hskip 1em\relax
%jhw       \fi
%jhw       ##1}}}

\pagestyle{colabhead}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sans Serif font for sections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\def\@sectionfont{\bf\sffamily}
\renewcommand{\section}{\@startsection {section}{1}{\z@}%
  {-3.5ex \@plus -1ex \@minus -.2ex}%
  {2.3ex \@plus.2ex}%
  {\reset@font\Large\@sectionfont}}
\renewcommand{\subsection}{\@startsection{subsection}{2}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {1.5ex \@plus .2ex}%
  {\reset@font\large\@sectionfont}}
\renewcommand{\subsubsection}{\@startsection{subsubsection}{3}{\z@}%
  {-3.25ex\@plus -1ex \@minus -.2ex}%
  {1.5ex \@plus .2ex}%
  {\reset@font\normalsize\@sectionfont}}
\renewcommand{\paragraph}{\@startsection{paragraph}{4}{\z@}%
  {3.25ex \@plus1ex \@minus.2ex}%
  {-1em}%
  {\reset@font\normalsize\bfseries}}
\renewcommand{\subparagraph}{\@startsection{subparagraph}{5}{\parindent}%
  {3.25ex \@plus1ex \@minus .2ex}%
  {-1em}%
  {\reset@font\normalsize\bfseries}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chapter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Front cover page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\maketitle}{%
                                %
                                %   Outside front cover page
                                %
  \thispagestyle{empty}
  \enlargethispage{20cm}  % just to make sure LaTeX doesn't
                                % generate a new page
  \vspace*{-3.5cm}\hspace*{-2cm}
  \noindent
%jhw \fbox{
%jhw \begin{minipage}[t][27cm]{5cm}
%jhw   \includegraphics[width=\linewidth,clip]{back}\\
%jhw   \vfill
%jhw   \includegraphics[width=\linewidth]{ethlogo}
%jhw \end{minipage}
%jhw  }
\hfill
  \begin{minipage}[t][27cm]{15cm}
     %%\strut\\[-3mm]
     \if@COLABreport
     \textsf{COLAB report series, Vol.\ \@COLABreport}%
     \if@COLABreportdate%
     , \@COLABreportdate
     \else
     , {\large NO DATE GIVEN}
     \fi
     \\[3cm]
     \else
     \mbox{}\\[5cm]
     \fi
     \noindent\@maketitle
     \if@status
     \vfill
     {\bf\large\@statushead\hspace{2mm}} 
     {\normalsize\@statustail}
     \fi    
  \end{minipage}
                                %
                                %   Inside front cover page
                                %
  \newpage
\thispagestyle{empty}
       \let\COLAB@tempevenhead\@evenhead
%%jhw  \def\@evenhead{%
%%jhw    \parbox{\textwidth}{\sffamily\@title\hfill\\
%%jhw      \rule[1.5ex]{\textwidth}{0.1pt}}}
%%jhw  \enlargethispage{2cm}
                                %\vspace*{-\COLAB@vup}
\strut
  \vfill
  %% nu het colofon zetten
  \if@correspondingauthorset%
  \noindent\begin{tabular}{@{}lr@{}} \hline
    \begin{minipage}[t]{.55\textwidth}
      \vspace*{2mm}
      \COLAB@colofon
    \end{minipage}&
    \begin{minipage}[t]{.45\textwidth}
      \vspace*{2mm}
      \noindent\COLAB@correspondingauthor
    \end{minipage}\\
    \vspace*{0.5ex}\mbox{}&\\
  \end{tabular}
  \else
  \hline\\
  \begin{center}\COLAB@colofon\end{center}
  \fi
  \setcounter{page}\z@
  \newpage
\newpage
% \enlargethispage{2cm}
%\thispagestyle{empty}
%\thispagestyle{empty}
  \let\@evenhead\COLAB@tempevenhead
%\thispagestyle{empty}
\newpage
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make Title 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      \def\@titlefont{\sffamily}
      \def\@maketitle{%
        {\LARGE\@titlefont\noindent\@title\\[.5ex]}
        \if@subtitle{\Large \@titlefont\noindent\@subtitle}\fi
        \hspace{0pt}\\[1.5cm]
\strut\\\vfil
        \@formatauthorsandinstitutes
      }



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  {\sffamily\bfseries\small #1:}\hspace{.5em}#2
  \vskip\belowcaptionskip}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Back cover page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The back cover page contains a list of technical reports that are already
% published within the series.

\newcommand{\COLAB@backcoverpage}[1]{
    \if@acknowledge\noindent\vfill
      \section*{Acknowledgements}
      \@acknowledge
    \fi
    \newpage
    \enlargethispage{2cm}
    \let\@oddhead\COLAB@oddtemphead
    \thispagestyle{empty}
    \if@COLABreport
        \strut\par
        \vspace*{\fill}
    \fi
    \noindent
    %% FIXME\hfill\includegraphics{ethlogo}\\
    \strut    
}



\newcommand{\COLAB@reportsfromfile}{
    \IfFileExists{references.bbl}{
        \input{references.bbl}
    }%else
    {
        \fbox{The file with technical reports ``references.bbl'' is missing}
    }
}



\endinput

\jajaja % just to show that this code doesn't do anything

\newlength{\COLAB@vup}   \setlength{\COLAB@vup}{3cm}
\newlength{\COLAB@vdown} \setlength{\COLAB@vdown}{2cm}
\newlength{\COLAB@tph}   \setlength{\COLAB@tph}{\textheight}
\addtolength{\COLAB@tph}{\COLAB@vup}
\addtolength{\COLAB@tph}{\COLAB@vdown}
\addtolength{\COLAB@tph}{-1cm}
\newlength{\COLAB@hleft} \setlength{\COLAB@hleft}{2cm}
\newlength{\COLAB@tpw}   \setlength{\COLAB@tpw}{\textwidth}
\addtolength{\COLAB@tpw}{\COLAB@hleft}
\addtolength{\COLAB@tpw}{-.5cm}

\newlength{\COLAB@vbannerwidth}  \setlength{\COLAB@vbannerwidth}{3cm}
\newlength{\COLAB@tpew}  \setlength{\COLAB@tpew}{\COLAB@tpw}
\addtolength{\COLAB@tpew}{-\COLAB@vbannerwidth}
\addtolength{\COLAB@tpew}{-17mm}

\newsavebox{\COLAB@vbanner}
\sbox{\COLAB@vbanner}{\includegraphics[bb=0 0 3cm 28cm,width=\COLAB@vbannerwidth]{\COLAB@vbannerfile}}

\renewcommand{\maketitle}{%
    %
    %   Outside cover page
    %
    \thispagestyle{empty}
    \enlargethispage{\COLAB@vdown}
    \vspace*{-\COLAB@vup}%
    \hspace*{-\COLAB@hleft}\mbox{
        \begin{minipage}[t][\COLAB@tph][b]{\COLAB@vbannerwidth}
            \usebox{\COLAB@vbanner}
        \end{minipage}\hspace{.4cm}\rule[-\COLAB@tph]{.1mm}{\COLAB@tph}\hspace{1cm}\mbox{
        \begin{minipage}[t][\COLAB@tph][t]{\COLAB@tpew}
            \if@COLABreport
                \vspace*{1.5ex}\noindent
                \textsf{COLAB technical report series, Vol.\ }\@COLABreport\\[4.5cm]
            \else
                \mbox{}\\[5cm]
            \fi
            \noindent\@maketitle
            \if@status\vfill
                \noindent{\sffamily
                \newlength{\COLAB@shw}\newlength{\COLAB@stw}
                \settowidth{\COLAB@shw}{\large\@statushead\hspace{2mm}}
                \setlength{\COLAB@stw}{\COLAB@tpew}
                \addtolength{\COLAB@stw}{-\COLAB@shw}
                \
                        \begin{minipage}[t]{\COLAB@shw}
                            \large\@statushead\hspace{2mm}
                    \end{minipage}
                        \begin{minipage}[t]{\COLAB@stw}
                            \normalsize\@statustail
                    \end{minipage}}
            \fi
        \end{minipage}}
    }
    %
    %   Inside cover page
    %
    \newpage
    \let\COLAB@tempevenhead\@evenhead
    \def\@evenhead{%
        \parbox{\textwidth}{\sffamily\@title\hfill\leftmark\\
        \rule[1.5ex]{\textwidth}{1pt}}}
    \enlargethispage{\COLAB@vdown}
    %\vspace*{-\COLAB@vup}
    \tableofcontents
    \vspace*{\fill}
    %% nu het colofon zetten
    \if@correspondingauthorset%
        \noindent\begin{tabular}{@{}l|r@{}}
            \hline\\
                \begin{minipage}[t]{.5\textwidth}
                    \begin{center}\COLAB@colofon\end{center}
                \end{minipage}&
                \begin{minipage}[t]{.5\textwidth}
                \vfill
                    \begin{center}\COLAB@correspondingauthor\end{center}
                \vfill
                \end{minipage}\\
            \vspace*{0.5ex}\mbox{}&\\ \hline
        \end{tabular}
    \else
        \hline\\
        \begin{center}\COLAB@colofon\end{center}
    \fi
    \setcounter{page}\z@
    \newpage
    \let\@evenhead\COLAB@tempevenhead
}
