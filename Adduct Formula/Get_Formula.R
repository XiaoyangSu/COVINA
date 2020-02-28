library(lpSolveAPI)

Ions.List <- read.csv(file="E:/XSu/Papers/Adducts/Adduct Formula/Ions.csv",check.names = FALSE)
Input.DeltaM <- read.csv(file="E:/XSu/Papers/Adducts/Adduct Formula/Input_DeltaM.csv",check.names = FALSE)

Ions.Total <- nrow(Ions.List)

Get.Formula <- function(Target.Mass) {

  lpmodel<-make.lp(0,Ions.Total+2)
  
  set.bounds(lpmodel, upper = Ions.List$Max, columns = c(1:Ions.Total))
  set.bounds(lpmodel, upper = 1, lower = 1, columns = Ions.Total+1)
  
  add.constraint(lpmodel, c(Ions.List$Charge,0,0), "=", 0)
  add.constraint(lpmodel, c(Ions.List$Mass,-Target.Mass,-1), "<=", 0)
  add.constraint(lpmodel, c(Ions.List$Mass,-Target.Mass,1), ">=", 0)
  
  set.type(lpmodel,c(1:(Ions.Total)),type = "integer")
  set.objfn(lpmodel, c(rep(0,Ions.Total),-Target.Mass,1))
  
  lp.control(lpmodel,sense='min')
  
  solve(lpmodel)
  
  return(get.variables(lpmodel))
}

Formulae <- matrix(0,ncol=Ions.Total+2,nrow=nrow(Input.DeltaM))
for (i in 1:nrow(Input.DeltaM)) {
  Formulae[i,1] <- Input.DeltaM[i,1]
  Formulae[i,c(2:(Ions.Total+2))] <- Get.Formula(Input.DeltaM[i,1])[c(1:Ions.Total,Ions.Total+2)]
}

colnames(Formulae) <- c("Delta_m/z",as.character(Ions.List$Ions),"Mass error")
write.csv(Formulae,"E:/XSu/Papers/Adducts/Adduct Formula/Output_DeltaM.csv")
