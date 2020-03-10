library(lpSolveAPI)

Ions.List <- read.csv(file="E:/XSu/Papers/Adducts/Adduct Formula/Ions.csv",check.names = FALSE)
Input.DeltaM <- read.csv(file="E:/XSu/Papers/Adducts/Adduct Formula/Input_DeltaM.csv",check.names = FALSE)

Ions.Total <- nrow(Ions.List)

#The mixed integer programming aims to find the ion combinations that is closest to the Target.Mass
#This is to find the minimum of b that satisfies -b<=sum(ion.mass)-Target.Mass<=b

Get.Formula <- function(Target.Mass,Target.Mass.D2=0) {

  # the coefficients are for each of the ion in the list, Target.Mass and b 
  lpmodel<-make.lp(0,Ions.Total+2)
  
  #The bounds are the following
  #1. Each ion has an upper limit
  #2. The coefficient of Target.Mass is set to 1
  #3. If the delta mass in 2H-Ac is available (non-zero), the number of Ac is calculated and used as bound
  
  set.bounds(lpmodel, upper = Ions.List$Max, columns = c(1:Ions.Total))
  set.bounds(lpmodel, upper = 1, lower = 1, columns = Ions.Total+1)
  if(Target.Mass.D2!=0) {
    Ac.Number = round((Target.Mass.D2-Target.Mass)/3.0188)
    set.bounds(lpmodel, upper = Ac.Number, lower = Ac.Number, columns = 1)
  }
  
  #The constraints are the following
  #1. The "adduct" part is charge neutral. The X in adduct ion [M+X-H]- has charge of zero
  #2. sum(ion.mass)-Target.Mass<=b
  #3. sum(ion.mass)-Target.Mass>=-b
  
  add.constraint(lpmodel, c(Ions.List$Charge,0,0), "=", 0)
  add.constraint(lpmodel, c(Ions.List$Mass,-Target.Mass,-1), "<=", 0)
  add.constraint(lpmodel, c(Ions.List$Mass,-Target.Mass,1), ">=", 0)
  
  #All ions coefficients should be integers
  set.type(lpmodel,c(1:(Ions.Total)),type = "integer")
  #Minimize b
  set.objfn(lpmodel, c(rep(0,Ions.Total),0,1))
  
  lp.control(lpmodel,sense='min')
  
  solve(lpmodel)
  
  return(get.variables(lpmodel))
}

Formulae <- matrix(0,ncol=Ions.Total+2,nrow=nrow(Input.DeltaM))
for (i in 1:nrow(Input.DeltaM)) {
  Formulae[i,1] <- Input.DeltaM[i,1]
  Formulae[i,c(2:(Ions.Total+2))] <- Get.Formula(Input.DeltaM[i,1],Input.DeltaM[i,2])[c(1:Ions.Total,Ions.Total+2)]
}

colnames(Formulae) <- c("Delta_m/z",as.character(Ions.List$Ions),"Mass error")
write.csv(Formulae,"E:/XSu/Papers/Adducts/Adduct Formula/Output_DeltaM.csv")
