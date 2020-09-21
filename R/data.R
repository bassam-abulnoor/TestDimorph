#' Summary statistics of baboon data collection data frame A dataset
#' containing summary statistics for low density lipoprotein (LDL) and
#' apolipoprotein B (apo B) levels in 604 baboons measured on two different
#' diets: a basal diet 'chow' and a high cholesterol, saturated fat diet
#' 'pink' (HCSF). The baboons were classified into one of three subspecies
#' (Papio hamadryas anubis, P.h. cynocephalus, or anubistcynocephalus hybrid).
#' Each animal was measured on each of the two diets.
#'
#' @format A list of 7 matrices.
#' \describe{ \item{R.res}{pooled within group
#' correlation matrix}
#' \item{M.mu}{Means of lipoproteins in different species
#' for males}
#' \item{F.mu}{Means of lipoproteins in different species for
#' females}
#' \item{m}{Male sample sizes} \item{f}{Female sample sizes}
#' \item{M.sdev}{Standard deviations for males} \item{F.sdev}{Standard
#' deviations for females} }
#' @note The baboon data collection were supported by NIH grant HL28972 and
#' NIH contract HV53030 to the Southwest Foundation for Biomedical Research
#' (Now: Texas Biomedical Research Institute), and funds from the Southwest
#' Foundation for Biomedical Research
#'#' @references
#'
#' Konigsberg LW (1991). *An historical note on the t-test for differences in sexual dimorphism
#' between populations.*. American journal of physical anthropology, 84(1), 93–96.
#'
"baboon.parms_list"


#' Summary statistics of baboon data collection list A dataset containing
#' summary statistics for low density lipoprotein (LDL) and apolipoprotein B
#' (apo B) levels in 604 baboons measured on two different diets: a basal diet
#' 'chow' and a high cholesterol, saturated fat diet 'pink' (HCSF). The
#' baboons were classified into one of three subspecies (Papio hamadryas
#' anubis, P.h. cynocephalus, or anubistcynocephalus hybrid). Each animal was
#' measured on each of the two diets.
#'
#' @format A data frame with 12 rows and 8 variables \describe{
#' \item{Trait}{Type of apolipoprotein} \item{Sub}{Type of species}
#' \item{M.mu}{Means of lipoproteins in different species for males}
#' \item{F.mu}{Means of lipoproteins in different species for females}
#' \item{m}{Male sample sizes} \item{f}{Female sample sizes}
#' \item{M.sdev}{Standard deviations for males} \item{F.sdev}{Standard
#' deviations for females}
#'
#'  }
#'
#' @note The baboon data collection were supported by NIH grant HL28972 and
#' NIH contract HV53030 to the Southwest Foundation for Biomedical Research
#' (Now: Texas Biomedical Research Institute), and funds from the Southwest
#' Foundation for Biomedical Research
#' @references
#'
#' Konigsberg LW (1991). *An historical note on the t-test for differences in sexual dimorphism
#' between populations.*. American journal of physical anthropology, 84(1), 93–96.
#'
#'
"baboon.parms_df"

#' Pooled within group correlation matrix for baboon data
#' @format A 4*4 numerical matrix
"R"

#' The Howells' craniometric data
#'
#' A subset of a dataset that consists of 82 craniometric measurements taken
#' from approximately two thousands and half human crania from 28
#' geographically diverse populations.
#'
#' @format A data frame with 441 rows and 10 variables:
#'  \describe{
#' \item{Sex}{'M' for male and 'F' for female}
#' \item{Pop}{Populations' names}
#' \item{GOL}{Glabello occipital length}
#' \item{NOL}{Nasio occipital length}
#' \item{BNL}{Bastion nasion length}
#' \item{BBH}{Basion bregma height}
#' \item{XCB}{Maximum cranial breadth}
#' \item{XFB}{Max frontal breadth}
#' \item{ZYB}{Bizygomatic breadth}
#' \item{AUB}{Biauricular breadth}
#'   }
#' @references
#'
#' Howells WW (1989). *Skull shapes and the map: craniometric analyses in the dispersion of modern
#' Homo.*. Papers of the Peabody Museum of Archaeology and Ethnology, 79.
#'
"Howells"
