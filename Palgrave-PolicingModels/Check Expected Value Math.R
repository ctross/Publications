
# Simulation to check analytic version
      Tb<-0.01
      Tw<-0.005

      Pb<-0.015
      Pw<-0.005

      Nw<-100000
      Nb<-50000

      Res<-rep(NA,1000000)
      for(i in 1:1000000)
      Res[i] <-(rbinom(1,Nb,Pb*Tb)/Nb)/((rbinom(1,Nw,Pw*Tw)+1)/Nw)

      mean(Res)                                # Simulation Expected Value
     (Tb/Tw)*(Pb/Pw)*(1-(1-Tw*Pw)^Nw)          # Analytic Expected Value





















