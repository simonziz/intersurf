\NeedsTeXFormat{LaTeX2e}

\ProvidesClass{main}[2017/03/13 TELECOM Nancy Rapport de projet industriel, v1.0]

\LoadClass[a4paper,twoside,12pt]{report}

% appel des packages
\RequirePackage{graphicx}
\RequirePackage[french]{babel}
\RequirePackage{amsfonts,amssymb,amsmath,mathtools}
\RequirePackage[utf8]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{array}
\RequirePackage{libertine}
\RequirePackage[libertine]{newtxmath}
\RequirePackage[svgnames]{xcolor}
\RequirePackage{color}
\RequirePackage{dirtytalk}
\RequirePackage{listings}
\RequirePackage{pdfpages}
\RequirePackage{comment}
\RequirePackage{csquotes}
\RequirePackage{graphicx}
%glossaire
\RequirePackage[toc]{glossaries}
%références
%\RequirePackage{biblatex}

%placement des figures
\RequirePackage[section]{placeins}
\RequirePackage{rotating}
\RequirePackage{tablefootnote}
\AddThinSpaceBeforeFootnotes
\FrenchFootnotes
\RequirePackage[a4paper]{geometry}
\geometry{top=2cm, bottom=2cm,right=2cm,left=2.5cm}
\RequirePackage[colorlinks,linkcolor=black,citecolor=black]{hyperref}
\RequirePackage{acronym}
\RequirePackage{etoolbox}
\RequirePackage{authblk}
\RequirePackage{draftwatermark}
\RequirePackage{xcolor}
\RequirePackage{titlesec}
%\RequirePackage[chapter]{tocbibind}
%remove toc from toc
\RequirePackage[nottoc]{tocbibind}
% dessins et autres symboles
% symbole
\RequirePackage{tipa}
\RequirePackage{siunitx}
% drawing
\RequirePackage{tikz}
\usetikzlibrary{arrows.meta}
\usetikzlibrary{positioning}
\usetikzlibrary{decorations.pathreplacing}
\usetikzlibrary{shapes,arrows}
% symbols
\RequirePackage{pifont}
\newcommand{\cmark}{\huge\color{green}\ding{51}}%
\newcommand{\xmark}{\huge\color{red}\ding{55}}%
% maths
\RequirePackage{dsfont}
%gestion des espaces
\RequirePackage{xspace}
%retour à la ligne dans un tableau
\RequirePackage{makecell}
%nested item
\RequirePackage{enumitem}


\usepackage{subcaption}





%%%% début de l'affichage texttt de Jean%%%%
% Le package times mais sans la fonte courier
% qui apparait avec \texttt.
\ProvidesPackage{times}[\filedate\space\fileversion\space Times PSNFSS2e package]
\renewcommand{\sfdefault}{phv}
\renewcommand{\rmdefault}{ptm}
\renewcommand{\ttdefault}{pcr}
\def\fm#1{\text{\texttt{#1}}}
%%%% FIN %%%%

%%% Définition des macros
\def\source{\fm{srce}}
\def\cible{\fm{cible}}
\def\cassource{\fm{cas-srce}}
\def\cascible{\fm{cas-cible}}
\def\basecas{\fm{base-de-cas}}
\def\sol(#1){\fm{sol(#1)}}
\def\aliae{\textsc{Aliae}\xspace}
\def\loria{\textsc{Loria}\xspace}
\def\MmeB{Madame B.\xspace}

\lstset{ %
  backgroundcolor=\color{white},   % choose the background color; you must add \usepackage{color} or \usepackage{xcolor}
  basicstyle=\footnotesize,        % the size of the fonts that are used for the code
  breakatwhitespace=false,         % sets if automatic breaks should only happen at whitespace
  breaklines=true,                 % sets automatic line breaking
  captionpos=b,                    % sets the caption-position to bottom
  commentstyle=\color{Orange},    % comment style
  extendedchars=false,              % lets you use non-ASCII characters; for 8-bits encodings only, does not work with UTF-8
  frame=single,                    % adds a frame around the code
  keepspaces=true,                 % keeps spaces in text, useful for keeping indentation of code (possibly needs columns=flexible)
  keywordstyle=\color{DarkOrchid},       % keyword style
  numbers=left,                    % where to put the line-numbers; possible values are (none, left, right)
  numbersep=5pt,                   % how far the line-numbers are from the code
  numberstyle=\tiny\color{DarkGray}, % the style that is used for the line-numbers
  rulecolor=\color{DarkGray},         % if not set, the frame-color may be changed on line-breaks within not-black text (e.g. comments (green here))
  showspaces=false,                % show spaces everywhere adding particular underscores; it overrides 'showstringspaces'
  showstringspaces=false,          % underline spaces within strings only
  showtabs=false,                  % show tabs within strings adding particular underscores
  stepnumber=1,                    % the step between two line-numbers. If it's 1, each line will be numbered
  stringstyle=\color{DarkGreen},     % string literal style
  tabsize=2,                     % sets default tabsize to 2 spaces
  %title=\lstname                   % show the filename of files included with \lstinputlisting; also try caption instead of title
}

\newcommand{\email}[1]{\href{mailto:#1}{\nolinkurl{#1}}}
\newcommand{\universityYear}[1]{\the\numexpr #1 - 1 \relax--\the#1}

\@sc@wm@stampfalse
\SetWatermarkText{Confidentiel}
\SetWatermarkScale{.8}
\DeclareOption{confidential}{\@sc@wm@stamptrue}
\ProcessOptions

\addto\captionsfrench{%
  \renewcommand\listfigurename{Liste des illustrations}}

\titleformat{\chapter}
  {\normalfont\LARGE\bfseries}{\thechapter}{1em}{}

%\setlength{\parindent}{0pt} % définir la taille de l'alinéa
\parskip 10pt % définir l'espace entre deux paragraphes

\def\includeCompanyLogo{%
  \@ifundefined{reportCompanyLogoPath}{}{\includegraphics[height=5em]{\reportCompanyLogoPath}}
}

\def\includeLabLogo{%
  \@ifundefined{reportLabLogoPath}{}{\includegraphics[height=5em]{\reportLabLogoPath}}
}

\def\includeSchoolLogo{%
  \@ifundefined{reportSchoolLogoPath}{}{\includegraphics[height=5em]{\reportSchoolLogoPath}}
}

\def\reportAuthor{Nicolas Warnet}
\def\reportAuthorAddress{26 rue de Mon Désert}
\def\reportAuthorCity{54000 Nancy} % Adresse de l'élève
\def\reportAuthorPhone{+33 (0)6 07 98 88 88} % Téléphone de l'élève
\def\reportAuthorEmail{\email{nicolas.warnet@telecomnancy.eu}} % Courriel de l'élève

\def\reportTitle{Personnalisation d'un épisode de planification de soins} % Titre du mémoire
\def\reportLongTitle{Personnalisation d'un épisode de planification de soins - Application à la prise en charge de l'insomnie} % Titre plus long du mémoire

\def\reportCompany{ALIAE} % Nom de l'entreprise d'accueil
\def\reportCompanyAddress{14 rue de la Maison Blanche} % Rue de l'entreprise
\def\reportCompanyCity{75013 Paris} % Ville de l'entreprise
\def\reportCompanyPhone{téléphone} % Téléphone de l'entreprise
\def\reportCompanyEmail{\email{mail@company.fr}}
\def\reportCompanyLogoPath{figures/blc.png} % Logo de l'entreprise

\def\reportLab{LORIA}
\def\reportLabAddress{615 rue du Jardin botanique}
\def\reportLabCity{54506 Vand\oe{}uvre-lès-Nancy}
\def\reportLabPhone{+33 (0)3 83 59 20 00}
\def\reportLabLogoPath{figures/logo_loria.jpg}

\def\reportSchool{TELECOM Nancy} % Nom de l'école
\def\reportSchoolAddress{193 avenue Paul Muller} % Rue de l'école
\def\reportSchoolCity{CS 90172 Villers-lès-Nancy} % Ville de l'école
\def\reportSchoolPhone{+33 (0)3 83 68 26 00} % Téléphone de l'école
\def\reportSchoolEmail{\email{contact@telecomnancy.eu}}
\def\reportSchoolLogoPath{figures/school-logo.pdf}

\def\place{Nancy} % Ville pour la signature pour l'engagement anti-plagiat
\def\date{\today} % Date pour la signature de l'engagement anti-plagiat

\def\reportType{Mémoire d'ingénieur} % le gros titre sur la première page
\def\reportSchoolDiplomaYear{3} % année diplomante (3A)

\def\reportIndustrialSupervisor{Philippe Jolivet} % Prénom Nom de l'encadrant industriel
\def\reportAcademicSupervisor{Bruno Pinçon} % Prénom Nom de l'encadrant académique
\def\reportLoriaSupervisorOne{Jean Lieber}
\def\reportLoriaSupervisorTwo{Emmanuel Nauer}

\def\reportYear{\universityYear{\year}} % Année universitaire

\makeatletter
\def\maketitle{%
\pagestyle{empty}
\begin{center}
\includegraphics[width=4.5cm]{\reportSchoolLogoPath} %logo de ton école
\hspace{9mm}
\includegraphics[width=4.5cm]{\reportCompanyLogoPath} % ne pas toucher, permet de créer l'espace entre les deux logos
\hspace{9mm}
\includegraphics[width=4.5cm]{\reportLabLogoPath} %logo de ton entreprise

\vfill

{\fontsize{40}{40}\selectfont \reportType{}}

\rule{\textwidth}{1pt}

\bigskip

{\fontsize{22}{22}\selectfont \reportTitle{}}

\medskip

{\fontsize{16}{16}\selectfont 4 septembre 2017} % date de rendu du rapport

\medskip

{\fontsize{16}{16}\selectfont \reportAuthor{}}

\end{center}

\vfill

\begin{tabular}{@{}lr@{}}
\parbox[b]{.7\linewidth}{
\reportSchool{}\\
\reportSchoolAddress{},\\
\reportSchoolCity{}\\
\reportSchoolPhone{}\\
\smallskip

\reportCompany{} \\
\reportCompanyAddress{} \\
\reportCompanyCity{} \\
\reportCompanyPhone{} \\
\smallskip

\reportLab{} \\
\reportLabAddress{} \\
\reportLabCity{} \\
\reportLabPhone{} \\
} &
\parbox[b]{.7\linewidth}{
\textbf{Encadrant universitaire} \\
\reportAcademicSupervisor{}\\
\smallskip

\textbf{Encadrant industriel} \\
\reportIndustrialSupervisor{} \\
\smallskip

\textbf{Encadrants de recherche} \\
\reportLoriaSupervisorOne{} \\
\reportLoriaSupervisorTwo{} \\
\smallskip
}

\end{tabular}
\cleardoublepage

\newpage
\setcounter{page}{1}
\pagestyle{plain}
}
\makeatother
\makeatletter
\newcommand\insertAntiPlagiarismAgreement[2]{%

  \begin{center}
    \Large \bf Déclaration sur l'honneur de non-plagiat
  \end{center}

  \medskip
  \medskip

  {\bf Je soussigné(e),}

  {\bf Nom, prénom : #1}

  {\bf \'Elève-ingénieur(e) régulièrement inscrit(e) en \reportSchoolDiplomaYear{}\up{e} année à \reportSchool{}}

  {\bf Numéro de carte de l'étudiant(e) :} #2

  {\bf Année universitaire :} \reportYear{}

  {\bf Auteur(e) du document, mémoire, rapport ou code informatique intitulé :}

  \begin{center}
    \Large \reportLongTitle{}
  \end{center}
  \medskip

  Par la présente, je déclare m'être informé(e) sur les différentes formes de
  plagiat existantes et sur les techniques et normes de citation et référence.

  Je déclare en outre que le travail rendu est un travail original, issu de ma
  réflexion personnelle, et qu'il a été rédigé entièrement par mes soins.
  J'affirme n'avoir ni contrefait, ni falsifié, ni copié tout ou partie de
  l'\oe{}uvre d'autrui, en particulier texte ou code informatique, dans le but
  de me l'accaparer.

  Je certifie donc que toutes formulations, idées, recherches, raisonnements,
  analyses, programmes, schémas ou autre créations, figurant dans le document
  et empruntés à un tiers, sont clairement signalés comme tels, selon les
  usages en vigueur.

  Je suis conscient(e) que le fait de ne pas citer une source ou de ne pas la
  citer clairement et complètement est constitutif de plagiat, que le plagiat
  est considéré comme une faute grave au sein de l'Université, et qu'en cas de
  manquement aux règles en la matière, j'encourrais des poursuites non
  seulement devant la commission de discipline de l'établissement mais
  également devant les tribunaux de la République Fran\c{c}aise.

  \medskip
  \medskip

  \hspace{5cm}{\bf Fait à \place, le  \date{}}

  \hspace{5cm}{\bf Signature :}
\cleardoublepage
}


\makeatletter
\def\makesecondtitle{%
\begin{center}

{\fontsize{40}{40}\selectfont \reportType{}}
\rule{\textwidth}{1pt}

\bigskip

{\fontsize{22}{22}\selectfont \reportTitle{}}

\medskip

{\fontsize{18}{18}\selectfont \bf \reportAuthor{}}

\medskip

{\fontsize{14}{14}\selectfont \bf Année \reportYear{}}

\vfill
Stage de fin d'études réalisé dans l'entreprise \reportCompany{}
en partenariat avec le \reportLab{}
en vue de l'obtention du diplôme d'ingénieur de \reportSchool{}
\vfill
\end{center}

\begin{tabular}{@{}lr@{}}
\parbox[b]{.7\linewidth}{
\reportAuthor{} \\
\reportAuthorAddress{} \\
\reportAuthorCity{}\\
\reportAuthorPhone{} \\
} \end{tabular}

\begin{tabular}{@{}lr@{}}
\parbox[b]{.7\linewidth}{
\reportSchool{} \\
\reportSchoolAddress{} \\
\reportSchoolCity{} \\
\reportSchoolPhone{} \\
} & \includeSchoolLogo{} \\
\end{tabular}

\begin{tabular}{@{}lr@{}}
\parbox[b]{.7\linewidth}{
\reportCompany{} \\
\reportCompanyAddress{} \\
\reportCompanyCity{} \\
\reportCompanyPhone{} \\
} & \includeCompanyLogo{} \\
\end{tabular}

\begin{tabular}{@{}lr@{}}
\parbox[b]{.7\linewidth}{
\reportLab{} \\
\reportLabAddress{} \\
\reportLabCity{} \\
\reportLabPhone{} \\
} & \includeLabLogo{} \\
\end{tabular}

\medskip

\cleardoublepage
}
\makeatother
