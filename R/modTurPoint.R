#' @title Modified Turning Point Method for ED50 Estimation
#' @description 50 percent effective dose (abbreviated as ED50) is one of the most concerns
#' in the field of Anesthesiology. There are many methods to estimate ED50 and its confidence
#' interval. Turning point method is one of them proposed by Choi in 1990 and here is the
#' modified function to realize the ED50 estimation.
#' @param doseSeries A numeric vector. It can be the whole dose sequence given in a group of
#' experiments in the order. And it can also be the turning point sequence calculated.
#' @param onlyTurPoint A logical value indicating whether the \code{doseSeries} is a turning
#' point sequence or not, the defaut is value is \code{FALSE}.
#' @param confidence A value ranging from 0 to 1 that represents the confidence level. The
#' default value is 0.95.
#' @export
#' @return A list:
#' \item{Estimate}{The point estimate of ED59 using modified turning point method.}
#' \item{Confidence}{Confidence value input.}
#' \item{CI.lower}{The lower bound of confidence interval under a certain confidence level.}
#' \item{CI.upper}{The upper bound of confidence interval under a certain confidence level.}
#' @references Choi SC Interval estimation of the LD50 based on an up-and-down experiment. Biometrics 1990; 46: 485-92.
#' @examples
#' library(modTurPoint)
#' dose1 <- c(3.1, 3.2, 3.3, 3.2, 3.1, 3.2, 3.3, 3.2, 3.3)
#' modTurPoint(doseSeries = dose1, confidence = 0.9)
#' dose2 <- c(3.25, 3.15, 3.25, 3.25)
#' modTurPoint(doseSeries = dose2, onlyTurPoint = TRUE, confidence = 0.9)

modTurPoint <- function(doseSeries, onlyTurPoint = FALSE, confidence = 0.95)
{
  # First check out the data input:
  # Whether the dose data is a numeric vector...
  if(!is.vector(doseSeries))
    return(warning('The dose data input is not a vector!'))
  if(!is.numeric(doseSeries))
    return(warning('The type of data input is not numeric!'))

  # If the data input consists of the whole dose sequence,
  # the turning points are to be found out.
  if(!onlyTurPoint)
  {
    # Compute the fixed dose step and check if there is a mistake.
    doseDiff <- round(diff(doseSeries), 10)
    doseStep <- unique(abs(doseDiff))
    if(length(doseStep) > 1)
      return(warning('The dose step is not fixed!'))
    if(length(doseStep) == 0)
      return(warning('There is only one dose sample input!'))
    if(length(doseStep) == 1 & doseStep == 0)
      return(warning('The dose step is set to be zero!'))

    # Compute all the turning point values.
    temp1 <- diff(doseDiff)
    temp2 <- which(abs(temp1) == (2 * doseStep))
    if(length(temp1) == 0 | length(temp2) == 0)
      return(warning('The dose sequence input has no turning points!'))
    doseTurPoint <- (doseSeries[temp2] + doseSeries[temp2 + 1]) / 2
  }

  # If the data input consists of just the turning point values,
  # directly use them to calculate the modified turning point estimate.
  if(onlyTurPoint)
    doseTurPoint <- doseSeries

  # Compute the number of turning points.
  nTurPoint <- length(doseTurPoint)

  # Calculate ED50 estimate.
  if(nTurPoint == 1)
  {
    cat('ED50 estimate:', doseTurPoint, '\n')
    return(message('Additional message: There is only one turning point value!'))
  }
  meanTurPoint <- mean(doseTurPoint[-1])

  # Centralizing the turning point values and
  # define a, b, c mentioned in our article.
  cenTurPoint <- doseTurPoint - meanTurPoint
  a <- cenTurPoint[2]^2 + cenTurPoint[nTurPoint]^2
  b <- sum(cenTurPoint[-c(1, nTurPoint)] * cenTurPoint[-c(1, 2)])
  c <- sum(cenTurPoint[-c(1, 2, nTurPoint)]^2)

  # Solve the systerm of equations in the article.
  rhoHatFunc <- function(rho)
  {
    (b - rho * c) * (1 - rho^2) - rho * (a - 2 * rho * b + (1 + rho^2) * c) / (nTurPoint - 1)
  }
  rhoHat <- tryCatch(uniroot(rhoHatFunc,
                             interval = c(-1 + 10^(-5), 1-10^(-5)),
                             tol = 1e-9)$root, error = function(e) 'error')
  if(rhoHat == 'error') return(warning('Problem occured in solving the equations!\nMaybe the number of turning points is not enough.'))
  sigHat <- ((a - 2 * rhoHat * b + (1 + rhoHat^2) * c) / (nTurPoint - 1))^(0.5)

  # Compute standard error of ED50 estimate.
  i <- 1:(nTurPoint - 2)
  sd <- ((sigHat^2 / ((nTurPoint - 1) * (1 - rhoHat^2))) *
        (1 + sum((2 * (nTurPoint - i - 1) * rhoHat^i) / (nTurPoint - 1))))^(0.5)

  # Summarise ED50 estimate and its confidence interval.
  zAlpha <- qnorm(0.5 + confidence / 2)
  estimate <- meanTurPoint
  ciLower <- meanTurPoint - zAlpha * sd
  ciUpper <- meanTurPoint + zAlpha * sd
  output <- list(Estimate = estimate, Confidence = confidence,
                 CI.lower = ciLower, CI.upper = ciUpper)
  return(output)
}
