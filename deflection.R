#provides deflection
deflcalc <- function(actor, beh, object) {
  
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
  postepa <- cbind(store[2:10], postepa) #put the initial EPA and post EPA in df together
  
  #calculate deflection
  deflection <- apply(postepa, 1, diff) #get difference b/w pre and post EPA
  deflection <- as.data.frame(deflection)
  deflection <- sapply(deflection, function(x) x^2) #sq differences
  deflection <- apply(deflection, 2, sum) #sum squared differences :) 
  return(deflection)
}
