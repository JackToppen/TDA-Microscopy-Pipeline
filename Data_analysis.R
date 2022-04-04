
##### Retrieve data
Gata6=list()
HA=list()
concentration=c(0,5,15,25)
for(signal in 1:2)
{
  
  for(dapi in 1:4)
  {
    if(signal==1){Gata6[[dapi]]=list()}else{HA[[dapi]]=list()}
    for(well in 1:3)
    {
      if(signal==1){Gata6[[dapi]][[well]]=list()}else{HA[[dapi]][[well]]=list()}
      for(location in 1:5)
      {
        
        if(signal==1){Sig="Nanog_Gata6"}else{Sig="Nanog_HA"}
        
        Conc=concentration[dapi]
        setwd(paste0("~/Research/TDA_Data/xlsx_normalized/",Sig,"/concentration_",Conc))
        file=paste0(Sig,"_",Conc,"_",well," (raw tiles) ",location,".csv")
        
        if(signal==1)
        {
          Gata6[[dapi]][[well]][[location]]=read.csv(file)
        }else{
          HA[[dapi]][[well]][[location]]=read.csv(file)
        }
        
      }
    }
  }
}

#### Process GATA6 Data

# Set Percentiles
up_percentileR=.75
up_percentileG=.75


uP_RFP=list()
pop_RFP=list()
uP_GFP=list()
pop_GFP=list()

Mgata=list()
for(dapi in 1:4)
{
  
  Mgata[[dapi]]=list()
  for(well in 1:3)
  {
    uP_RFP[[well]]=list()
    pop_RFP[[well]]=list()
    uP_GFP[[well]]=list()
    pop_GFP[[well]]=list()
    Mgata[[dapi]][[well]]=list()
    for(location in 1:5)
    {
      RFP1=Gata6[[1]][[well]][[location]][,7] 
      GFP1=Gata6[[1]][[well]][[location]][,6]
      
      RFP4=Gata6[[4]][[well]][[location]][,7] 
      GFP4=Gata6[[4]][[well]][[location]][,6]
      
      uP_RFP[[well]][[location]]=quantile(RFP1 ,probs = up_percentileR) 
      uP_GFP[[well]][[location]]=quantile(GFP4 ,probs = up_percentileG) 
      
      pop_RFP[[well]][[location]]=length(RFP1)
      pop_GFP[[well]][[location]]=length(GFP4)
      
      uP_RFPv=uP_RFP[[well]][[location]]
      uP_GFPv=uP_GFP[[well]][[location]] 
      
      RFP=Gata6[[dapi]][[well]][[location]][,7] 
      GFP=Gata6[[dapi]][[well]][[location]][,6]
      
      # pdf(paste0("RFP_dox",concentration[dapi],"_well",well,"_location",location,".pdf")) 
      # qqnorm(RFP);qqline(RFP,col='red')
      # dev.off() 
      # 
      # pdf(paste0("GFP_dox",concentration[dapi],"_well",well,"_location",location,".pdf")) 
      # qqnorm(GFP);qqline(GFP,col='red')
      # dev.off() 
      
      M00=sum(((RFP<uP_RFPv) + (GFP<uP_GFPv))==2)/length(RFP)
      M11=sum(((RFP>uP_RFPv) + (GFP>uP_GFPv))==2)/length(RFP)
      M01=sum(((RFP<uP_RFPv) + (GFP>uP_GFPv))==2)/length(RFP)
      M10=sum(((RFP>uP_RFPv) + (GFP<uP_GFPv))==2)/length(RFP)
      
      Mgata[[dapi]][[well]][[location]]=matrix(c(M10,M00,M11,M01),ncol=2)
    }
  }
}

thresh_RFP=list()
thresh_GFP=list()
temp_sumRFP=list()
temp_sumGFP=list()
temp_multRFP=list()
temp_multGFP=list()

for(location in 1:5)
{
  temp_sumRFP[[location]]=0
  temp_sumGFP[[location]]=0
  temp_multRFP[[location]]=0
  temp_multGFP[[location]]=0
  
  for(well in 1:3)
  {
    temp_sumRFP[[location]]=pop_RFP[[well]][[location]][[1]]+temp_sumRFP[[location]]
    temp_sumGFP[[location]]=pop_GFP[[well]][[location]][[1]]+temp_sumGFP[[location]]
    
    temp_multRFP[[location]]=(uP_RFP[[well]][[location]]*pop_RFP[[well]][[location]][[1]])+temp_multRFP[[location]]
    temp_multGFP[[location]]=(uP_GFP[[well]][[location]]*pop_GFP[[well]][[location]][[1]])+temp_multGFP[[location]]
  }
  
  thresh_RFP[[location]]=temp_multRFP[[location]]/temp_sumRFP[[location]]
  thresh_GFP[[location]]=temp_multGFP[[location]]/temp_sumGFP[[location]]
}

## Output threshold tables (based on location in well)
thresh_RFP
thresh_GFP

## Output population distrbution tables
Mgata0 =matrix(rep(0,4),ncol=2)
Mgata5 =matrix(rep(0,4),ncol=2)
Mgata15=matrix(rep(0,4),ncol=2)
Mgata25=matrix(rep(0,4),ncol=2)


for(well in 1:3)
{
  for(location in 1:5)
  {
    Mgata0 =Mgata[[1]][[well]][[location]] +Mgata0
    Mgata5 =Mgata[[2]][[well]][[location]] +Mgata5
    Mgata15=Mgata[[3]][[well]][[location]] +Mgata15
    Mgata25=Mgata[[4]][[well]][[location]] +Mgata25
  }
}

Mgata0/sum(Mgata0)
Mgata5/sum(Mgata5)
Mgata15/sum(Mgata15)
Mgata25/sum(Mgata25)



#### Process HA Data

# Set Percentiles
up_percentileR=.75
up_percentileG=.75


uP_RFP=list()
pop_RFP=list()
uP_GFP=list()
pop_GFP=list()

MHA=list()
for(dapi in 1:4)
{
  
  MHA[[dapi]]=list()
  for(well in 1:3)
  {
    uP_RFP[[well]]=list()
    pop_RFP[[well]]=list()
    uP_GFP[[well]]=list()
    pop_GFP[[well]]=list()
    MHA[[dapi]][[well]]=list()
    for(location in 1:5)
    {
      RFP1=HA[[1]][[well]][[location]][,7] 
      GFP1=HA[[1]][[well]][[location]][,6]
      
      RFP4=HA[[4]][[well]][[location]][,7] 
      GFP4=HA[[4]][[well]][[location]][,6]
      
      uP_RFP[[well]][[location]]=quantile(RFP1 ,probs = up_percentileR) 
      uP_GFP[[well]][[location]]=quantile(GFP4 ,probs = up_percentileG) 
      
      pop_RFP[[well]][[location]]=length(RFP1)
      pop_GFP[[well]][[location]]=length(GFP4) 
      
      uP_RFPv=uP_RFP[[well]][[location]]
      uP_GFPv=uP_GFP[[well]][[location]] 
      
      RFP=HA[[dapi]][[well]][[location]][,7] 
      GFP=HA[[dapi]][[well]][[location]][,6]
      
      # pdf(paste0("HA_RFP_dox",concentration[dapi],"_well",well,"_location",location,".pdf")) 
      # qqnorm(RFP);qqline(RFP,col='red')
      # dev.off() 
      # 
      # pdf(paste0("HA_GFP_dox",concentration[dapi],"_well",well,"_location",location,".pdf")) 
      # qqnorm(GFP);qqline(GFP,col='red')
      # dev.off() 
      
      M00=sum(((RFP<uP_RFPv) + (GFP<uP_GFPv))==2)/length(RFP)
      M11=sum(((RFP>uP_RFPv) + (GFP>uP_GFPv))==2)/length(RFP)
      M01=sum(((RFP<uP_RFPv) + (GFP>uP_GFPv))==2)/length(RFP)
      M10=sum(((RFP>uP_RFPv) + (GFP<uP_GFPv))==2)/length(RFP)
      
      MHA[[dapi]][[well]][[location]]=matrix(c(M10,M00,M11,M01),ncol=2)
    }
  }
}

thresh_RFP=list()
thresh_GFP=list()
temp_sumRFP=list()
temp_sumGFP=list()
temp_multRFP=list()
temp_multGFP=list()

for(location in 1:5)
{
  temp_sumRFP[[location]]=0
  temp_sumGFP[[location]]=0
  temp_multRFP[[location]]=0
  temp_multGFP[[location]]=0
  
  for(well in 1:3)
  {
    temp_sumRFP[[location]]=pop_RFP[[well]][[location]][[1]]+temp_sumRFP[[location]]
    temp_sumGFP[[location]]=pop_GFP[[well]][[location]][[1]]+temp_sumGFP[[location]]
    
    temp_multRFP[[location]]=(uP_RFP[[well]][[location]]*pop_RFP[[well]][[location]][[1]])+temp_multRFP[[location]]
    temp_multGFP[[location]]=(uP_GFP[[well]][[location]]*pop_GFP[[well]][[location]][[1]])+temp_multGFP[[location]]
  }
  
  thresh_RFP[[location]]=temp_multRFP[[location]]/temp_sumRFP[[location]]
  thresh_GFP[[location]]=temp_multGFP[[location]]/temp_sumGFP[[location]]
}

## Output threshold tables (based on location in well)
thresh_RFP
thresh_GFP

## Output population distrbution tables
MHA0 =matrix(rep(0,4),ncol=2)
MHA5 =matrix(rep(0,4),ncol=2)
MHA15=matrix(rep(0,4),ncol=2)
MHA25=matrix(rep(0,4),ncol=2)


for(well in 1:3)
{
  for(location in 1:5)
  {
    MHA0 =MHA[[1]][[well]][[location]] +MHA0
    MHA5 =MHA[[2]][[well]][[location]] +MHA5
    MHA15=MHA[[3]][[well]][[location]] +MHA15
    MHA25=MHA[[4]][[well]][[location]] +MHA25
  }
}

MHA0/sum(MHA0)
MHA5/sum(MHA5)
MHA15/sum(MHA15)
MHA25/sum(MHA25)