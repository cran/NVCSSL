\name{yeast}
\alias{yeast}

\title{Yeast cell cycle data set}
\description{This data set contains the yeast cell cycle data that is analyzed in Bai et al. (2020), which is a subset of the original yeast data analyzed by Spellman et al. (1998). There are \eqn{n=47} yeast genes that were determined to periodically expressed by De Lichtenberg et al. (2005). For each of these genes, mRNA levels were measured at \eqn{n_i=18} time points from time 0 to 119 minutes (7 minutes apart), for a total of \eqn{N=846} data points. For covariates, we have the binding information for 96 transcription factors (TFs) taken from the  chromatin immunoprecipitation data set of Lee et al. (2002).   
 }
\usage{data(yeast)}
\format{
A dataframe with the following columns. 
\describe{
    \item{\code{mRNA}:}{mRNA levels for 47 yeast genes at 18 different time points.}
    \item{\code{Time}:}{time points from 0-119 minutes.}
    \item{\code{Gene}:}{gene identifiers (1-47).}
    \item{\code{ABF1-ZMS1}:}{Columns 4-99 contain the binding information on 96 TFs that are possibly associated with the yeast genes.}
}
}

\references{ 
Bai, R., Boland, M. R., and Chen, Y. (2020). "Fast algorithms and theory for high-dimensional Bayesian varying coefficient models." \emph{arXiv preprint arXiv:1907.06477}.

De Lichtenberg, U., Jensen, L. J., Fausboll, A., Jensen T. S., Bork, P., and Brunak, S. (2005). "Comparison of computational methods for the identification of cell cycle-regulated genes." \emph{Bioinformatics}, \bold{21}: 1164-1171.

Lee, T.I., Rinaldi, N.J., Robert, F., Odom, D.T., Bar-Joseph, Z., Gerber, G.K., Hannett, N.M., Harbison, C.T., Thomson, C.M., Simon, I., Zeitlinger, J., Jennings, E.G., Murray, H.L., Gordon, D.B., Ren, B., Wyrick, J.J., Tagne, J.B., Volkert, T.L, Fraenkel, E, Gifford, D.K., and Young, R.A (2002). "Transcriptional regulatory networks in Saccharomyces cerevisiae." \emph{Science}, \bold{298}: 799-804.

Spellman, P.T., Sherlock, G., Zhang, M.Q., Iyer, V.R., Anders, K., Eisen, M.B., Brown, P.O., Botstein, D., and Futcher, B. (1998). "Comprehensive identification of cell cycle-regulated genes of the yeast Saccharomyces cerevisiae by microarray hydrization. \emph{Molecular Biology of the Cell}, \bold{9}: 3273-3279.
}