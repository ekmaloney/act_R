#libraries
library(tidyverse)

#downloaded from UGA website - should be most recent US one
dictionary <- 
  read.csv("C:/Users/emily/Documents/Grad School/Non-class Research/ACT-occ/Data/FullSurveyorUSMeans.csv")

#load equation information from Interact
equation <- read_table("C:/Users/emily/Documents/Grad School/Non-class Research/ACT-occ/equation info.txt", col_names = FALSE)
equation<- as.data.frame(equation) #put in better format

#label according to what term the coefficients correspond to
colnames(equation) <- c("Term", "Aep", "App", "Aap", "Bep", "Bpp", "Bap", "Oep", "Opp", "Oap")
rownames(equation) <- c("Constant", "Aeb", "Apb", "Aab", "Beb", "Bpb", "Bab", "Oeb", "Opb", "Oab",
                        "AeBe", "AeOp", "ApBp", "AaBa", "BeOe", "BeOp", "BpOe", "BpOp", "AeBeOe",
                        "AeBeOp")

#get rid of the Z000 column
equation <- equation[,-1]

transientimp <- function(actor, beh, object) {
  
  #make empty matrix to store info
  store <- matrix(NA, 20, 1)
  
  #something that won't affect the constant
  store[1,] <- 1
  
  #setting actor EPA values
  store[2,] <- dictionary$E[dictionary$term == actor] #Ae
  store[3,] <- dictionary$P[dictionary$term == actor] #Ap
  store[4,] <- dictionary$A[dictionary$term == actor] #Aa
  
  #setting behavior EPA values
  store[5,] <- dictionary$E[dictionary$term == beh] #Be
  store[6,] <- dictionary$P[dictionary$term == beh] #Bp
  store[7,] <- dictionary$A[dictionary$term == beh] #Ba
  
  #setting object EPA values
  store[8,] <- dictionary$E[dictionary$term == object] #Oe
  store[9,] <- dictionary$P[dictionary$term == object] #Op
  store[10,] <- dictionary$A[dictionary$term == object] #Oa
  
  store[11,] <- store[2,] * store[5,] #Ae*Be
  store[12,] <- store[2,] * store[9,] #Ae*Op
  store[13,] <- store[3,] * store[6,] #Ap*Bp
  store[14,] <- store[4,] * store[7,] #Aa*Ba
  store[15,] <- store[5,] * store[8,] #Be*Oe
  store[16,] <- store[5,] * store[9,] #Be*Op
  store[17,] <- store[6,] * store[8,] #Bp*Oe
  store[18,] <- store[6,] * store[9,] #Bp*Op
  store[19,] <- store[2,] * store[5,] * store[8,] #Ae*Be*Oe
  store[20,] <- store[2,] * store[5,] * store[9,] #Ae*Be*Op
  
  #to store the EPA values for actor, behavior, and object after situation
  postepa <- matrix(NA, 20, 9)
  
  #apply function 
  i<- 1
  for(i in 1:9) {
    postepa[,i] <- equation[,i]*store
    i+1
  }
  
  #put in better data format
  postepa <- as.data.frame(postepa)
  postepa <- apply(postepa, 2, sum) #get the sum of the equation to get the actual EPA values
  preepa <- store[2:10]
  
  #terms that have behavior in them
  z <- c(1, 1, 1, preepa[4], preepa[5], preepa[6], 1, 1, 1, 
         1, 1, 1, 1, postepa[4], postepa[5], postepa[6], 1, 1, 1,
         postepa[4], 1, postepa[5], postepa[6],
         postepa[4], postepa[4], postepa[5], postepa[5], postepa[4], postepa[4])
  z <- as.vector(z)
  
  #terms without behavior in them
  i <- c(preepa[1], preepa[2], preepa[3], 1, 1, 1, preepa[7], preepa[8], preepa[9], 1,
         postepa[1], postepa[2], postepa[3], 1, 1, 1, postepa[7], postepa[8], postepa[9],
         postepa[1], postepa[1], postepa[2], postepa[3], postepa[7], postepa[8],
         postepa[7], postepa[8], postepa[1]*postepa[7], postepa[1]*postepa[8])
  
  #make into matrix with these on the diagonal
  mat_i <- matrix(0, 29, 29)
  diag(mat_i) <- i
  
  #selection matrix
  #row 1 = where Be terms are in z
  r1 <- matrix(0, 1, 29)
  r1[1, 4] <- 1
  r1[1, 14] <- 1
  r1[1, 20] <- 1
  r1[1, 24] <- 1
  r1[1, 25] <- 1
  r1[1, 28] <- 1
  r1[1, 29] <- 1
  
  
  #row 2 = where Bp terms are in z
  r2 <- matrix(0, 1, 29)
  r2[1, 5] <- 1
  r2[1, 15] <- 1
  r2[1, 22] <- 1
  r2[1, 26] <- 1
  r2[1, 27] <- 1

  #row 3 = where Ba terms are in z
  r3 <- matrix(0, 1, 29)
  r3[1, 6] <- 1
  r3[1, 16] <- 1
  r3[1, 23] <- 1
  
  s <- rbind(r1, r2, r3)
  s <- as.matrix(s)
  
  g <- c(1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  
  #identity matrix
  identity <- matrix(0, 9, 9)
  diag(identity) <- 1
  
  #coefficient matrix
  m <- as.matrix(equation)
  
  #making h matrix with identity and coefficient info
  h1 <- rbind(identity, -1*m)
  h2 <- cbind(identity, -1*t(m))
  h <- h1 %*% h2
  
  #doing the parts of the equation for the solution
  term1 <- s %*% mat_i %*% h %*% mat_i %*% t(s)
  term1 <- solve(term1)
  term1 <- -1*term1
  
  term2 <- s %*% mat_i %*% h %*% mat_i %*% g
  
  sol <- term1 %*% term2

  return(sol)
}

test <- transientimp("immigrant", "scrutinize", "infant")
test
