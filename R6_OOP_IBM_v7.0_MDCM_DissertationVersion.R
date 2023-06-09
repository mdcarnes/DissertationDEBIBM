#Dynamic Energy Budget Individual Based Model of a typical Ozark Arkansas Timber Rattlesnake
#Built by Max Carnes-Mason for his Doctoral Dissertation entitled: "Shedding in the Timber Rattlesnake: Natural Patterns, Endocrinological Underpinnings, Temporal and Energetic Effort, and Integration as a Reptilian Life History Trait"
#Copywrite Maxwell D. Carnes-Mason, 2023
#All Rights Reserved

#This iteration is the base model.  It can be adapted and output programmed to access variables of interest.  

#If you are using this model and need help or want to talk to the author, please email me at; carnesmdm@gmail.com



#clear the r brain
rm(list=ls())

#load libraries
library(R6) #use the R6 library (S4 can be used with the {methods} library)
library(reshape2)
library(ggplot2)
library(dplyr)
library(uuid)


#Define Relevant Simulation Parameters
#Simulation Details
#nInd  <- 100                #number of snakes
#years <- 1                  #number of years
#days  <- 365*years          #number days
#Foraging Routine Details
MFS   <- 0.085              #Mean Foraging Success
LAM   <- 500                #Largest Available Mouse
MMS   <- 0.05               #a proproation of body size that an animal will not take a meal smaller than
#Digestion and passage time Routine Details
#DT    <- NULL               #edit: set within function forage #digestion time, set at 12 via beaupre and zaidan, 2012
convEff <- 0.80             #efficieny with which mouse is turned into metabolizable energy.  Beaupre 2001 sets it at 0.80
#Growth Routine Details
growthEff <- 0.6            #the efficieny of converting energy to biomass (60% via Beaupre 2001)
#Repro details
maturationSize <- 70        #set the size of sexual maturation (70cm = 248g = 186g maxiumum meal for females)
eSize <- 25                 #the weight (in g) of newborn snakes
lifeSpan <- 25              #number of years alive  allowed in the population

#shed parameters
shedLength <- 28            #set the duration of sheds
propShed <- 0.5             #proportion of the length of shed during which animals are unaffected by the shed

#survival parameters                     #probabilities of survival in different age/sex classes.  currently informed from the table in Olson et al., 2015
neonateSurvivorship <- 0.65
JuvenileSurvivorship <- 0.92
MaleSurvivorship <- 0.95
FemaleSurvivorship <- 0.835




#Define a snake class
Snake <- R6Class("Snake",
                 public = list(                  #define the class and slots
                 ID = NULL,                      #each ind. will have a unique ID
                 Tb = NULL,                      #Body Temperature
                 Mass = NULL,                    #Wet mass of individual
                 SVL = NULL,                     #Snout Vent Length (numeric, cm)
                 sexMature = logical(),          #is the animal sexually mature yet? true is mature false is juvenile
                 LTS = NULL,                     #Long term storage reserve (numeric)
                 STS = NULL,                     #Short term storage
                 safetyStore = NULL,             #amount of lts energy needed in reserve
                 FinG = NULL,                    #Food in Gut? a wet mass of rodent (in g) given by function forage
                 digestTime = NULL,              #Digestion time remaining, 0=no food in gut
                 alloGrowth = NULL,              #allocation of resources to growth
                 alloRepro = NULL,               #allocation of resources to reproduction
                 probGrowth = NULL,              #the ability to grow based on season, 0=no, 1=yes
                 probRepro = NULL,               #the ability to repro based on season, 0=no, 1=yes
                 shedDay = NULL,                 #counter to count down days remaining of a shed
                 vitel = logical(),              #vitellogenic?
                 shedEnergy = NULL,              #counter to keep track of total energy of shedding
                 daysVitel = NULL,               #counter to keep track of vitellogenesis
                 Er = NULL,                      #energetic value of skin removal
                 Eb = NULL,                      #energetic value of skin synthesis
                 Es = NULL,                      #energetic value sequesterd in shed skin 
                 totalMetabLife = NULL,          #keep track of total metabolic cost through whole life
                 SDAtotal = NULL,                #store the total SDA energy, calculated at meal aquisition, to be paid as increase in Met during digetion time.
                 SDAdaily = NULL,                #save the daily SDA cost at foraging, reset when digestion is complete (in "digest")
                 nnProducedYearly = NULL,        #save the number of nn produced in this year
                 yearRepro = logical(),          #save T/F whether an animal had a litter this year.
                 offspring = tibble::lst(),      #hold the creation of offspring with the individual for addition to the population
                 birthYear = NULL,               #track the birth year of new individuals
                 Dead = logical(),               #true/false if an animal has died
                 totalEB = NULL,                 #counter of total energy gained during the lifetime
                 initialize=function(ID = NA){   #function to populate a new ind.
                   self$ID <- ID                 #ID must be supplied upon ind. creation
                   self$Tb <- NA                 #Set body temp based on season in model
                   self$Mass <- 25               #default natal mass ~25g avg. from lab born data set n=13
                   self$SVL <- 33.2836           #SVL is equal to 11.73 *W^0.324 where W is wet mass, from Beaupre, 2001
                   self$sexMature <- F           #animals start as juveniles (when under a set SVL)
                   self$LTS <- 135                #animals start with 84kJ reserve energy  *update with data*   **sim animals kept dying at 25kJ**  edit 6Feb2023: Steve thinks maybe his original model should've read "5g fat" storage rather than 5kJ, which is approx. 135kJ *edit 7 Feb: programmed 84 kJ LTS just below the winter safety store at birth, so they will not grow right away.
                   self$STS <- 0                 #animals start with no STS
                   self$safetyStore <- NA        #safety store (in days) will be set by function reserveLTS
                   self$FinG <- 0                #animals start with no food in gut
                   self$digestTime <- 0          #digest time starts at 0
                   self$alloGrowth <- 1          #proportion of resources to growth
                   self$alloRepro <- 0           #proportion of resources to repro
                   self$probGrowth <- 1          #ability to grow based on season, 0=no, 1=yes
                   self$probRepro <- 0           #ability to repro based on season, 0=no, 1=yes
                   self$shedDay <- 0             #animals do not start in shed
                   self$vitel <- F               #all animals start out as not vitellogenic
                   self$shedEnergy <- 0          #no initial energy of shed
                   self$daysVitel <- 0           #number of days into the vitellogenic cycle (for all sexes, but should only populate for females)
                   self$Er <- 0                  #cost of removing skin set at start of shed generation (does not alter as animal grows during shed)
                   self$Eb <- 0                  #cost of synthesizing skin set at start of shed generation (does not alter as animal grows during shed)
                   self$Es <- 0                  #cost lost in skin set at start of shed generation (does not alter as animal grows during shed)
                   self$totalMetabLife <- 0      #start at 0
                   self$SDAtotal <- 0            #SDA starts at 0
                   self$SDAdaily <- 0            #daily Met cost of SDA starts at 0 when gut is empty
                   self$nnProducedYearly <- 0    #no neonates produced in a year at first
                   self$yearRepro <- F           #T/F if an animal has reproduced this year, starts as false.
                   self$offspring <- NA          #hold the offpsring created within the loops.
                   self$birthYear <- 1           #initialize at 1, but alter at birth.
                   self$Dead <- F                #they start out alive.
                   self$totalEB <- 84            #they have the energy from thier stored fat at the begining
                   }))

#make an environment class to make an object that is the environmental conditions
Environment <- R6Class("Environment",                                      #create an environment class
                       public = list(                                      #its got public attributes
                         day = NULL,                                       #environment will track day and ust that to set season and temp
                         #timeofDay = NULL,                                 #time of day will be set by loops later                                      
                         season = c("Spring","Summer","Winter","Autumn"),  #seaon can be one of 4 options              
                         temperature = NULL,                               #temp will set Tb of animals later                   
                         initialize = function(startDay=1){            #creation of environment requires input of a day and time
                           self$day = startDay                             #day populated at start                             
                           #self$timeofDay = time                           #time populated at start                       
                           self$season =  NA                               #season set by setSeason function in loops using solstice and equionx cutoffs.                         
                           self$temperature = NA                           #temp set by setTemp fucntion loop later                    
                             }))




#create male and female subclasses, they will have different methods for some things
#male snakes
maleSnake <- R6Class("maleSnake",                       #create male subclass
                     inherit = Snake,                   #inherit all valuse from super class "snake"
                     public = list(            
                       sex = c("Male","Female"),        #include a new slot, with options male and female
                       initialize = function(ID = NA){  #create a new initailize method to run when the animal is generated
                         super$initialize(ID=NA)        #use all of the slots and intiialize method of the super class
                         self$sex <- "Male"             #Set sex to "Male"
                         }
                     ))
#female snakes
femaleSnake <- R6Class("femaleSnake",                     #create female subclass
                       inherit = Snake,                   #inherit all values form superclass Snake
                       public = list(
                         sex = c("Male", "Female"),       #include a new slot with options male and female
                         nnProd = NULL,                   #number of neonates produced
                         reproStore = NULL,               #storage for Reproduction
                         reproCount = NULL,               #number of reproductive events
                         minEmbryoE = NULL,               #save the minimum amount of repro energy required for repro to be worth it that year.
                         dailyMCV = NULL,                 #daily cost of vitellogenesis will be saved here whenever a female is reproductive
                         initialize = function(ID = NA){  #create a new initilaize method that
                           super$initialize(ID = NA)      #run the initialize method from the superclass
                           self$sex <- "Female"           #add that we wnat to populate the sex with female for these individuals upon creation
                           self$nnProd <- 0               #animals start with no offspring
                           self$reproStore <- 0           #animals start not ready to reproduce.
                           self$dailyMCV <- 0             #start at nothing.
                           self$minEmbryoE <- 0           #start at nothing.
                           self$reproCount <- 0           #no intial repro events.
                           }
                       ))



#Define methods

#############################
#Methods for Snakes
#############################

#Forage Protocol
forage <- function(object) {                               #function to determine foraging success
   if(environment$season == "Winter"){                     #snakes do not forage in winter
     NULL}else{                                            #in all other seasons;
    mouse <- runif(1,0,1);                    #a number between 0 and 1 is generated from a uniform distribution
    largestMealtaken <- ((object$SVL^2.3502)/(235.99))
    mouseSize <- runif(1,min=(object$Mass*MMS), max = largestMealtaken);                  #mouseSize is a number between 5% of body mass and the largest meal the animal can take (by Rulons data) from a uniform dist.
    if(mouseSize > LAM){mouseSize = LAM}           #if the mouse is larger than avaiable, set it to the habitat maximum.
    maxMeal <- largestMealtaken                        #the biggest meal an animal can take based on its SVL and rulons data.
    if(mouseSize <= maxMeal & mouseSize >= (object$Mass*MMS)){   #if the mouse is not too big to eat, but bigger than the minimum meal size (a % of body size which an animal will not go below)
    if(mouse <= MFS){                                         #if the randomly selected probability of catching a mouse is less than the MFS (set above)
      mealProp <- mouseSize/object$Mass                         #calculate the proportion of body size that the meal is
      if(mealProp <= 0.15){                                    #if the meal is less than 15% body mass;
        duration <- ((22.25*(log10(object$Mass))) + 60.95)/24     #gives in hours, divide by 24 to get days.
        duration <- ceiling(duration)                            #round up to nearest whole day
      }else{
      if(mealProp > 0.15 & mealProp <= 0.40){                    #if meal size is between 15-40%; use the 30% meal size scalar
        duration <- ((116.14*(log10(object$Mass))) + 6.88)/24     #gives in hours, divide by 24 to get days.
        duration <- ceiling(duration)                            #round up to nearest whole day
      }else{
      if(mealProp > 0.40){                                       #if meal size is >40%, use the 50% meal size scalar
        duration <- ((150.33*(log10(object$Mass))) - 34.10)/24     #gives in hours, divide by 24 to get days.
        duration <- ceiling(duration)                            #round up to nearest whole day
      }else{
#        print(mealProp)
 #       print("DT scalar not working")
        }}}                        #if this routine is not catching all cases, tell me.
      DT <- duration                                            #set DT equal to the calculated duration
      # DT <- ((1975 * (object$Mass ^ 0.65)*(environment$temperature^-1.98))/24)  #calculate the digest time from the formula from Beaupre 2001 (divided by 24 to make it daily)
      # DT <- ceiling(DT)                                    #round the DT up to the nearest whole day.
      if(DT < 3){DT <- 3}else{NULL}                        #if digest time is less than 3 days, set it to 3 days (Beaupre, 2001)
 #     print(paste0("the digest time for the food item is ", DT))
 #     print(paste0("the meal is ", mealProp, " percent of body mass"))
      object$digestTime <- DT;                             #edit:changed to a scaling DT#12 days to digest a meal from beaupre and zaidan 2012 (no effect of size, body temp, sex)
      object$FinG <- mouseSize                           #Food in Gut is changed to a mass of food in gut (g)
      #object$SDAtotal <- (623.2498 + (-0.8603*object$Mass) + (83.3628*mouseSize))  #from Funk and Ortega (honors thesis, not in press yet); volume of CO2 produced during digetsiton....This value is not right, they are too large when compared to zaidan and beaupre (below)...could be a cumtrapz issue-JO
      object$SDAtotal <- ((12.158)*(object$Mass^0.160)*(mouseSize^1.030))    #following zaidan and beupare, 2003
      object$SDAtotal <- ((object$SDAtotal/27.42)/1000)      #convert from mL CO2 to kJ total and store to SDA total variable
      object$SDAdaily <- (object$SDAtotal / object$digestTime)  #set the daily amount of SDA.  this is reset in "digest" whne DT reaches 0.  should only trigger when gut is empty then filled by foraging.
      }else{                                                 #if the foraging event is unsuccessful, 
      NULL                                                 #do nothing
    }} else {NULL}                                         #if the meal is bigger than they can take, do nothing, they do not feed.
     }
#  print("forage")
  }

#Digestion (simple, almost identical to 2001 model)
digest <- function(object){
  digestPerDiem <- object$FinG/object$digestTime;         #amount of prey digested per day
  rodentDryMassPerDiem <- (digestPerDiem*0.25);           #calculate dry mass of food
  rodentEC <- rodentDryMassPerDiem*23.311;                #convert dry mass to energy content 23.311 kJ per g dry mass
  MECdaily <- rodentEC*convEff                            #convert energy content to daily metabolized energy multiplying by a set conversion efficiency (set at start)
  object$totalEB <- object$totalEB + MECdaily             #add metabolizable energy to the counter of total lifetime energy
  object$FinG <- object$FinG - digestPerDiem;             #new amount of food at end of day
  object$digestTime <- object$digestTime - 1;             #reduce number of days remaining by 1
  if(object$digestTime == 0){                             #if digest time now equals zero (it is complete)
    object$SDAtotal <- 0
    object$SDAdaily <- 0                                  #reset the daily and total SDA values, no more digestion, no more SDA
  }
  object$STS <- object$STS + MECdaily                     #add daily MEC to the short term stored energy
#  print("digest")
}

#updated metabolism routine using FMR scalar equation to produce variable scope values which change with body size.
metab <- function(object){
  if(environment$season == "Winter"){       #winter temp cycles
    z <- (-1.0768*(environment$temperature[1])) + (0.1091*((environment$temperature[1])^2)) - (0.00317*((environment$temperature[1])^3))      #winter RMR from agugliero Diss., assuming 10C for the duration of the season
    dailyRMR <- 24 * ((0.7546) * (object$Mass^1.0457) * (10^z))
  }
  if(environment$season == "Spring"){                   #spring temp cycles
    dailyRMR <- ((4*((0.00107)*(object$Mass^0.825)*(10^(0.0569*(environment$temperature[1])))))+ #from 1500 to 1800 (4 hours)
                   (5*((0.00091)*(object$Mass^0.799)*(10^(0.0628*(environment$temperature[2])))))+ #from 1900 to 2300 (5 hours)
                   (4*((0.00095)*(object$Mass^0.741)*(10^(0.0680*(environment$temperature[3])))))+ #from 0000 to 0300 (4 hours)
                   (4*((0.00120)*(object$Mass^0.727)*(10^(0.0643*(environment$temperature[4])))))+#from 0400 to 0700 (4 hours)
                   (4*((0.00124)*(object$Mass^0.777)*(10^(0.0590*(environment$temperature[5])))))+ #from 0800 to 1100 (4 hours)
                   (3*((0.00128)*(object$Mass^0.787)*(10^(0.0650*(environment$temperature[6])))))) #from 1200 to 1400 (3 hours) #this one isnt' in the beaupre and zaidan table.  just making it up
  }
  if(environment$season == "Summer"){                  #summer temp cycles
    dailyRMR <- ((4*((0.00107)*(object$Mass^0.825)*(10^(0.0569*(environment$temperature[1])))))+ #from 1500 to 1800 (4 hours)
                   (5*((0.00091)*(object$Mass^0.799)*(10^(0.0628*(environment$temperature[2])))))+ #from 1900 to 2300 (5 hours)
                   (4*((0.00095)*(object$Mass^0.741)*(10^(0.0680*(environment$temperature[3])))))+ #from 0000 to 0300 (4 hours)
                   (4*((0.00120)*(object$Mass^0.727)*(10^(0.0643*(environment$temperature[4])))))+#from 0400 to 0700 (4 hours)
                   (4*((0.00124)*(object$Mass^0.777)*(10^(0.0590*(environment$temperature[5])))))+ #from 0800 to 1100 (4 hours)
                   (3*((0.00128)*(object$Mass^0.787)*(10^(0.0650*(environment$temperature[6])))))) #from 1200 to 1400 (3 hours) #this one isnt' in the beaupre and zaidan table.  just making it up
  }
  if(environment$season == "Autumn"){                 #fall temp cycles
    dailyRMR <- ((4*((0.00107)*(object$Mass^0.825)*(10^(0.0569*(environment$temperature[1])))))+ #from 1500 to 1800 (4 hours)
                   (5*((0.00091)*(object$Mass^0.799)*(10^(0.0628*(environment$temperature[2])))))+ #from 1900 to 2300 (5 hours)
                   (4*((0.00095)*(object$Mass^0.741)*(10^(0.0680*(environment$temperature[3])))))+ #from 0000 to 0300 (4 hours)
                   (4*((0.00120)*(object$Mass^0.727)*(10^(0.0643*(environment$temperature[4])))))+#from 0400 to 0700 (4 hours)
                   (4*((0.00124)*(object$Mass^0.777)*(10^(0.0590*(environment$temperature[5])))))+ #from 0800 to 1100 (4 hours)
                   (3*((0.00128)*(object$Mass^0.787)*(10^(0.0650*(environment$temperature[6])))))) #from 1200 to 1400 (3 hours) #this one isnt' in the beaupre and zaidan table.  just making it up
  }
  dailyMC <- ((dailyRMR/1000)*27.42)       #convert RMR per day to kJ per day (mL C02 per day) -> 27.42 J per mL -> 1000 J/kJ
  # dailyMC <- dailyMC/5.9                #temporary scalar. steves 2nd model gives ~360J per day active met in a 29g snake in its first year.
  if(environment$season == "Summer"){
    prop <- 0.9073                              #percentage of total active season RMR that happens in the summer
    lengthDays <- 182                          #182 days of summer
  }
  if(environment$season == "Autumn"){
    prop <- 0.0459                              #percentage of total active season RMR that happens in the fall
    lengthDays <- 31                          #31 days of fall
  }
  if(environment$season == "Spring"){
    prop <- 0.0467                              #percentage of total active season RMR that happens in the spring
    lengthDays <- 31                          #31 days of fall
  }
  if(environment$season == "Winter"){
    prop <- 0                                 #no FMR during Winter, set prop to 0 and length to 1 to produce a scope of 0, preventing scope application
    lengthDays <- 1
    }                                        #in winter do nothing
  FMR <- (0.447*(object$Mass^1.13))          #daily FRM active season (mL CO2) from beuapre et al BORII, 2017
  dailyFMR <- (((FMR*244)*prop)/lengthDays)  #I *think* the FMR equation is daily values.  multiply by length of active season, then take the proportion for the season and divide by the number of days in that season
  dailyFMR <- ((dailyFMR/1000) * 27.42)          #convert from mL CO2 to kJ
  Scope <- dailyFMR/dailyMC
#  print(paste0("the daily kJ RMR is ", dailyMC))
#  print(paste0("scope for the day is ", Scope))
  if(object$digestTime > 0){                      #digestion increases cost of metabolism, on days when digesting,
    dailyMet <- dailyMC + (object$SDAdaily)                         #daily metabolic cost (kJ) increases by daily SDA on top of RMR - assumes activity is minimal at this time.
  }else{                                        #if the animal is digesting
    if(Scope < 1.000){                          #if scope is less than 1 (FMR equation is less than RMR daily)
      dailyMet <- dailyMC                       #no activity extra
    }else{                                      #if FMR exceeds RMR, multiply RMR by the calculated scope
      dailyMet <- dailyMC * Scope
    }
};
#  print(paste0("the daily met scaled is ", dailyMet))
  object$totalMetabLife <- object$totalMetabLife + dailyMet   #add a counter to track total metabolic cost of living.
  if(object$STS < dailyMet){                      #if STS is not sufficient to cover daily metabolism
    object$LTS <- object$LTS - (dailyMet - object$STS); #lower LTS by the difference between the 2,
    object$STS <- 0                              #and set STS to 0
  }else{                                         #if sts can cover daily metabolism
    object$STS <- object$STS - dailyMet}        #subtract daily Met from STS.
  # print("metab")
}

#Growth 
growth <- function(object){
  # EA <- object$LTS - object$safetyStore                           #the energy available for processes is equal to total LTS minus Safety net
  GE <- (object$STS  * object$alloGrowth * object$probGrowth) * growthEff   #energy avaialable for growth given by beaupre 2001
#  print(paste0("energy put towards growth is ", GE))
  #  if(EA > 0){                                                     #if the energy available for growth is greater than 0
    if(GE > 0){                                                   #do not grow if there is 0 allocated to growth
      snakeGrowthWet <- (GE/23.287)/(1-0.75)                      #available grow energy divided by density of tissue (dry) gives g of dry growth, divided by 1(total)-0.75 (percent water) gives an increase from dry to wet mass)
#      print(paste0("pregrowth SVL is", object$SVL))
      if(object$sex == "Male"){                                   #calculate maxiumum growth rates (in SVL per day) from Beaupre et al., 2017
        growMax <- (0.0507 - (0.00039*object$SVL)) * 100            #max cm/day in males (doubled from observed in our pop, assuming less food limiting)
        }else{
        growMax <- (0.0472 - (0.00043*object$SVL)) * 100            #max cm/day in females (doubled from observed in our pop, assuming less food limiting)
        }
      newSVL <- object$SVL + growMax                              #calculate a new SVL
      newMass <- ((newSVL/11.73)^3.08641975)                      #convert new svl (after max growth) to new mass using the inverse SVL/Mass relationship
      changeMass <- newMass - object$Mass                         #calculate a difference in mass (g wet growht)
      growMax <- changeMass                                       #update the max grow to give grams of wet tissue
      if(snakeGrowthWet > growMax){                               #if the snake is set to grow more than the daily max
        snakeGrowthWet <- growMax                                 #set the daily growth to the max for an animal of that size
        spentMax <- 27.287*(0.25*growMax)                         #convert from wetmass grown max to dry (25% dry tissue) to mass of tissue (27.287 kJ/gram dry)
        object$Mass <- object$Mass + snakeGrowthWet               #increase the snakes mass by the energy spent on growth.
        object$STS <- object$STS - (spentMax/growthEff)           #subtract the amount of kJ spent to grow (including that lost from inefficienty of matter conversion) from STS
#        print(paste0("the rate limited growth (in kJ from STS) is ", (spentMax/growthEff)))
        }else{                                                      #if the growth in length will not be greater than the daily maxiumum:
      object$Mass <- object$Mass + snakeGrowthWet                 #add wet mass to the snake, it will grow in length too.
      object$STS <- object$STS - (object$STS*object$alloGrowth*object$probGrowth) #save the amount of enregy used here to subtract from STS in model
      }
      if(object$SVL <= ((11.73)*object$Mass^0.324)){            #if the new mass of the snake gives a larger svl/
      object$SVL <- (11.73)*(object$Mass^0.324)                   #the svl will update and scale relative to mass
      }else{NULL}                                                 #if the new length would be shorter, don't update sVL, snakes cant shrink.
    }else{NULL}
#  print(paste0("growth: new mass is ", object$Mass))
#  print(paste0("growth: new SVL is ", object$SVL))
  }#else{NULL}}

#Repro Rules
reproduction <- function(object){
  EA <- object$reproStore                                #the energy available for processes is stored in repro store.
  RE <- EA * growthEff # * object$probRepro                              #repro energy is equal to availalbe energy * efficiency * probability of reproducting 
    if(object$sexMature == T){                      #if mature
      if(object$sex == "Male"){NULL}                    #if male     #this will be some amount of mate searching, so increased Metabolic cost, but need to fill in
      if(object$sex == "Female"){                       #if female
        availBody <- object$SVL * 0.55                           #only 55% of the body can be used for repro
        numEmbryo <- availBody / 4                               #each embryo takes up 4 cm
        minEmbryo <- numEmbryo - 3                               #minimum embryos for a viable pregnancy is possible - 3
        object$minEmbryoE <- (minEmbryo * (eSize*(1-0.75)*27.369))/growthEff        #the amount of energy for the minimum number of embryos (embryo size, dry mass, energy density) i is the day variable in the day loop divided by the growth eff scalar to account for energy loss.
  ###this chunk dictates vitellogenic status (primary initiation) and runs the daysVitel counter      
        if(i == 185 & object$vitel == F){                                       #on the first day of spring, make a decision about reproduction
          if(object$LTS >= (((object$safetyStore/30)*122) + object$minEmbryoE)){    #if the available energy from LTS is enough to produce a clutch and survive a winter;
            object$vitel <- T                                              #set the female as vitellogenic. 
            object$daysVitel <- object$daysVitel + 1                               #and set the counter for vitellogenesis to 0
          }else{object$vitel <- F}                                 #if they won't survive a pregnancy, do not reproduce this year.
          } else {                                                 #any other day but the first day of spring,
            if(object$vitel == T){                                 #if vitellogenic any other day the events of day 185, add one to the counter
              object$daysVitel <- object$daysVitel + 1             #add one to the day counter
          }else{NULL};                                        #if they are not vitellogenic, do nothing
        };                                            #do nothing
  ###this chunk dictates the cost of vitellogenesis if they fertilize at the start of spring 2 of the cycle  
        if(object$daysVitel == 366){                                       #if an animal entered vitellogeneisis last year, at the start of the next spring check if it is in "good shape" (ie can survive and make a litter with LTS), then start yolking.  the previous year, only water is invested more or less, so it is free, but dependent on bocy condition
          if(object$LTS >= (((object$safetyStore/30)*122) + object$minEmbryoE)){ #if the animal has enough stores to survive the winter (summer store of 10 days * 122 days for winter) they will continue the vitellogenesis
            energyYolk <- minEmbryo * 220                                  #the amount of energy needed for yolking is the minimum number of embryos times the enregy density of yolk (220kJ per ovum, Beaupre BORII)
            MCV <- ((1.087*(object$Mass^0.5324))*(energyYolk^0.3278))       #MCV in kJ is given by the scaling relationship from BORII 1.087*mass^0.5324 * Y (yolk energy)^0.3278
            object$dailyMCV <- (MCV + energyYolk) / 177                     #vitellogenesis has a daily cost equal to the cost of yolk plus the metabolic cost divided by the number of days (177 spring and summer)
            }else{                                              #if the animal does not have enough stores to survive the pregnancy, abort
             object$vitel <- F                                             #the animal aborts the pregnancy, not vitell
             object$daysVitel <- 0                                         #the counter resets.
             object$LTS <- object$LTS + object$reproStore                  #return repro energy (if any) to LTS
             object$reproStore <- 0                                        #reset reproStore
             }
          }
  ###this chunk allocates  the cost of vitellogenesis each day.  
        if(object$daysVitel > 366 & object$daysVitel < 546){                             #from the start of vitellogenisis to the end of the summer you need to save energy from your income or take it from LTS
          if(object$vitel == T){                                     #for vitellogenic females
          if(object$dailyMCV > (object$STS*object$alloRepro)){                    #if the cost for the day is greater than the available amount from the days income
            object$reproStore <- object$reproStore + object$dailyMCV                     #increase the reproStore by the daily needed
            object$STS <- object$STS - (object$STS*object$alloRepro)                     #decrease STS by the amount allowed to allocate
            object$LTS <- object$LTS - (object$dailyMCV - (object$STS*object$alloRepro)) #take the remainder necessary from LTS
          }
          if(object$dailyMCV <= (object$STS*object$alloRepro)){                   #if STS*allorepro is greater than or equal to the daily minimum amount;  ***importantly, this means that on "good days" a snake may save more than required, allowing her to get to a number of offspring greater than the minimum viable number
            object$reproStore <- object$reproStore + (object$STS*object$alloRepro)#increase repro stores by the avaialble STS fraction
            object$STS <- object$STS - (object$STS*object$alloRepro)              #decrease STS by the same.
          }}else{NULL}                                               #non vitellogenic females do nothing.
        }
  ###this chunk dictates parturition
        if(object$daysVitel == 546) {                                     #on the day of parturition;
          if(object$vitel == T){                                     #if vitellogenic;
            if(object$reproStore < object$minEmbryoE){                                   #if the amount of energy saved is not enough for minimum clutch;
              object$vitel <- F                                                   #the animal will abort the preganancy
              object$daysVitel <- 0                                               #days vitel resets
              object$LTS <- object$LTS + object$reproStore                        #stored repro energy will put back into LTS
              object$reproStore <- 0                                              #the store will be emptied
              }
            if(object$reproStore >= object$minEmbryoE){                                                         #if the stored repro energy is enough to make some embryos;
              availBody <- object$SVL * 0.55                           #only 55% of the body can be used for repro
              numEmbryo <- availBody / 4                               #each embryo takes up 4 cm
              minEmbryo <- numEmbryo - 3                               #minimum embryos for a viable pregnancy is possible - 3
              embryoMax <- (RE / (eSize * (1-0.75)*27.369))                       #use the reprostore accumulated during vitellogenesis to compute the maximum number of embryos that can be made.
                if(embryoMax > numEmbryo){embryoMax <- numEmbryo}                     #if the number of embryos she can produce is greater than she has space for, set the number equal to how many she can hold,
                else{NULL}                                                            #otherwise leave it alone.
                embryoMax <- floor(embryoMax)                                     #round number of embryos down to a whole number
                usedRE <- embryoMax * (eSize * (1-0.75)*27.369)                   #the cost of energy of producing that number of embryos
                object$reproStore <- object$reproStore - usedRE                   #reduce stored energy by the amound used on repro (usable repro energy scaled by percent allocated to repro)
                object$LTS <- object$LTS + object$reproStore                      #any leftover energy from vitellogenesis is put back into LTS
                object$reproStore <- 0                                            #empty repro store
                object$nnProducedYearly <- embryoMax                              #save the number of embryos produced this year to add to the population
                object$nnProd <- object$nnProd + embryoMax                        #add the number of embryos produced to the female
                object$reproCount <- object$reproCount + 1                       #add 1 to the number of reproductive events to keep track.
                object$daysVitel <- 0                                             #pregnancy complete, reset the vitel status
                object$vitel <- F                                                 #reset the vitellogensis status, post partum
                object$yearRepro <- T                                             #set repro do true for this year
                }else{NULL}                                      #if repro store and minembryoE arent right
          }
          else{NULL}                                      #if not vitel., do nothing
          }
          else{NULL}                                      #if not one of the appropriate dates, do nothing
          }
          else{NULL}}else{NULL} 
  # print("repro")
  }                           #if not male or female, or mature, do nothing
        
#maturation      #set animals to mature at over 70cm SVL
sexMaturation <- function(object){
  if(object$sexMature == F & object$SVL >= maturationSize){ #if SVL reaches 70 for a juvenile snake,
    object$sexMature <- T                       #change sex mature to true
  }else{NULL}                                   #if not change nothing.
# print("sexMatur")
}

#LTS rules
reserveLTS <- function(object){
  #RMR is the cost per hour over 24 hours; from beaupre 2001, the relationship between temp and RMR per hour changes over the day
  if(environment$season == "Winter"){       #winter temp cycles
    z <- (-1.0768*10) + (0.1091*(10^2)) - (0.00317*(10^3))      #winter RMR from agugliero Diss., assuming 10C for the duration of the season
    dailyRMR <- 24 * ((0.7546) * (object$Mass^1.0457) * (10^z))
  }
  if(environment$season == "Spring"){                   #spring temp cycles
    dailyRMR <- ((4*((0.00107)*(object$Mass^0.825)*(10^(0.0569*25))))+ #from 1500 to 1800 (4 hours)
                   (5*((0.00091)*(object$Mass^0.799)*(10^(0.0628*20))))+ #from 1900 to 2300 (5 hours)
                   (4*((0.00095)*(object$Mass^0.741)*(10^(0.0680*17))))+ #from 0000 to 0300 (4 hours)
                   (4*((0.00120)*(object$Mass^0.727)*(10^(0.0643*17))))+#from 0400 to 0700 (4 hours)
                   (4*((0.00124)*(object$Mass^0.777)*(10^(0.0590*22))))+ #from 0800 to 1100 (4 hours)
                   (3*((0.00128)*(object$Mass^0.787)*(10^(0.0650*23))))) #from 1200 to 1400 (3 hours) #this one isnt' in the beaupre and zaidan table.  just making it up
  }
  if(environment$season == "Summer"){                  #summer temp cycles
    dailyRMR <- ((4*((0.00107)*(object$Mass^0.825)*(10^(0.0569*32))))+ #from 1500 to 1800 (4 hours)
                   (5*((0.00091)*(object$Mass^0.799)*(10^(0.0628*30))))+ #from 1900 to 2300 (5 hours)
                   (4*((0.00095)*(object$Mass^0.741)*(10^(0.0680*28))))+ #from 0000 to 0300 (4 hours)
                   (4*((0.00120)*(object$Mass^0.727)*(10^(0.0643*25))))+#from 0400 to 0700 (4 hours)
                   (4*((0.00124)*(object$Mass^0.777)*(10^(0.0590*28))))+ #from 0800 to 1100 (4 hours)
                   (3*((0.00128)*(object$Mass^0.787)*(10^(0.0650*32))))) #from 1200 to 1400 (3 hours) #this one isnt' in the beaupre and zaidan table.  just making it up
  }
  if(environment$season == "Autumn"){                 #fall temp cycles
    dailyRMR <- ((4*((0.00107)*(object$Mass^0.825)*(10^(0.0569*25))))+ #from 1500 to 1800 (4 hours)
                   (5*((0.00091)*(object$Mass^0.799)*(10^(0.0628*22))))+ #from 1900 to 2300 (5 hours)
                   (4*((0.00095)*(object$Mass^0.741)*(10^(0.0680*18))))+ #from 0000 to 0300 (4 hours)
                   (4*((0.00120)*(object$Mass^0.727)*(10^(0.0643*16))))+#from 0400 to 0700 (4 hours)
                   (4*((0.00124)*(object$Mass^0.777)*(10^(0.0590*20))))+ #from 0800 to 1100 (4 hours)
                   (3*((0.00128)*(object$Mass^0.787)*(10^(0.0650*23))))) #from 1200 to 1400 (3 hours) #this one isnt' in the beaupre and zaidan table.  just making it up
  }#else{NULL}                           #if its none of these seasons do nothing
  dailyMC <- ((dailyRMR*27.42)/1000)       #convert RMR per day to kJ per day (mL C02 per day) -> L per day -> kJ per L CO2 (Ch. 3 diss.; nagy and gessaman)
  #seasonal safeties include daily Metabolic cost, scaled per season (active rates, Metabolic Scope; Beaupre 2001), times the seasonal safetly length                               
  if(environment$season == "Winter"){                  #in winter the long term reserve is non existent    
    object$safetyStore <- 0*dailyMC*2}
  if(environment$season == "Spring"){                  #in spring the snake must have 30 days reserve energy to grow
    object$safetyStore <- 30*dailyMC*3}                    
  if(environment$season == "Summer"){                  #in summer the snake must have 30 days reserve energy to grow
    object$safetyStore <- 30*dailyMC*3.6}
  if(environment$season == "Autumn"){                  #in autumn the snake must have 122 days reserve energy to grow
    object$safetyStore <- 122*dailyMC*3}#else{NULL}
  # print("reserveLTS")
  }

#ability to grow/repro by season and sex
canGrowRepro <- function(object){
  if(object$sexMature == F){              #if a juvenile
    if(environment$season == "Winter"){         #in winter the snake does not grow or repro    
      object$probGrowth <- 0
      object$probRepro <- 0}
    if(environment$season == "Spring"){         #in spring the snake grows but does not repro
      object$probGrowth <- 1
      object$probRepro <- 0}
    if(environment$season == "Summer"){         #in summer the snake repro but doesnt grow
      object$probGrowth <- 1
      object$probRepro <- 0}
    if(environment$season == "Autumn"){         #in autumn the snake grows but doesnt repro
      object$probGrowth <- 1
      object$probRepro <- 0}}else{        #if mature
      if(object$sex == "Male"){             #mature males
       if(environment$season == "Winter"){         #in winter the snake does not grow or repro    
        object$probGrowth <- 0
        object$probRepro <- 0}
       if(environment$season == "Spring"){         #in spring the snake grows but does not repro
        object$probGrowth <- 1
        object$probRepro <- 0}
       if(environment$season == "Summer"){         #in summer the male snakes dont grow (they reproduce but not for now)
        object$probGrowth <- 0
        object$probRepro <- 0}
       if(environment$season == "Autumn"){         #in autumn the snake grows but doesnt repro
        object$probGrowth <- 1
        object$probRepro <- 0}
    
     }else{                                #mature females
      if(environment$season == "Winter"){         #in winter the snake doesnt grow or repro   
        object$probGrowth <- 0
        object$probRepro <- 0}
      if(environment$season == "Spring"){         #in spring the snake devotes energy to growth and repro
        object$probGrowth <- 1
        object$probRepro <- 1}
      if(environment$season == "Summer"){         #in summer the snake devotes energy to growth and repro
        object$probGrowth <- 0
        object$probRepro <- 1}
      if(environment$season == "Autumn"){         #in autumn the snake devotes energy to growth and repro
        object$probGrowth <- 1
        object$probRepro <- 0}else{NULL}}
  }
# print("cangrorepro")
}

#growth allocation   #percentage of resources devoted to each activity.
alloRules <- function(object){              #function to set rules for allocation to growth and repro
  if(object$sexMature == F){      #if a juvenile
    if(environment$season == "Winter"){     #in the winter, just survive.  
      object$alloGrowth <- 0
      object$alloRepro <- 0
      }else{                                #other seasons juveniles only grow
      object$alloGrowth <- 1                # 100% grow
      object$alloRepro <- 0}}               # 0% Repro
  if(object$sexMature == T){                           #if mature
    if(object$sex == "Male"){         #mature males
      if(environment$season == "Winter"){   #in the winter, just survive.  
        object$alloGrowth <- 0
        object$alloRepro <- 0}
      else{                                #in other seasons males only grow
      object$alloGrowth <- 1                #100% growth
      object$alloRepro <- 0}}                #0% repro
    if(object$sex == "Female"){       #mature females
      if(environment$season == "Winter"){   #in the winter, just survive.  
        object$alloGrowth <- 0
        object$alloRepro <- 0
        }else{                              #other seasons        
        object$alloGrowth <- 0.015          #1.5 % growth
        object$alloRepro <- 0.985}          #98.5 % repro
    }
  }
# print("allorules")
}

#shedding rules for the mandatory early season shed
#I set this according to the time lines for 2 shedders.  I assumed mean date of first observation
#was about 20 days after the actual start of a shed (usually opaque basal and blue eyes occurs ~ 20 days)
#I assumed mean date of first shed obs. was 10 june (+/- 10 days) so 10 may - 1 june for start of shed

#edited 30 Jan 2023.  using an approach where enetering sehd is 1 fucntion, the cost of shed is another.
#this was implemented to avoid animals in shed being counted multiple times in each day interation.
# 

#1.15 control
shedRules <- function(object){
###this chunk determines if an animal qualifies for its neonatal shed
    if(i==1 & fi == 1){                     #if the animal has just been born, place it into shed (nn shed)
      object$shedDay <- shedLength                  #start a shed    #and calculate shed values
      object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
      object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
      object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
    }                                              #do nothing if not the first day of the first year
###this chunk determines if an animal will undergo a shed in the first window
    if(object$shedDay <= 0 & (fi >=2 | (fi ==1 & i > 90))){        #if an animal is not in shed currently and its either 1) in its 1st full summer (year 2) OR 2) past its winter of its first year (arbitrarily 90 days post parturition here);
       if(i == 279){                    #if an animal is not yet in shed by day 279 force it into shed.
          shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
          if(shedLiklihood < 0.70){  #if the generated number is less than a threshold, say 70 percent of anmials shed in this window (arbitrary for now)
          object$shedDay <- shedLength         #enter shed and calculate values
          object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
          object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
          object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
          }}

    if(i >= 259 & i < 279){        #if the date is between may 10 and june 1 (+/- 10 days around june 10th mean date of first shed obs. in 2 shedders minus 20 (days between shed start and observable cues), CM personal data)
          shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
          if(shedLiklihood < 0.50){  #if the generated number is less than a threshold, 80 percent of animals shed in this window (arbitrary for now)
          object$shedDay <- shedLength         #enter shed and calculate values
          object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
          object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
          object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
          }
       }
###this chunk determines if an animal will shed in the 2nd window
    if(i >= 322 & i < 342){                   #if the date is between July 11 and july 31 (+/- 10 days around August 10th mean date of first shed obs. in 2 shedders minus 20 (days between shed start and observable cues), CM personal data)
        #if its not in shed in the 2nd window, give animals a chance to shed
        shedChance <- runif(1,0,1)                      #draw a random number between 0 and 1
        if(shedChance <= 0.0077091277){                     #animals shed 1.15 times per year (males 1.13, femlaes 1.18), meaning average 15% of animals shed 2x per year, if this animal is one of those 15 percent (by random number), it will shed a 2nd time
                                                            #since probability stacks over the 20 days, I set the total chance (1-(1-p)^21) to 0.15, giving a chance, per draw of 0.007....  see math in yellow notebook or (https://www.quora.com/How-would-one-calculate-the-probability-of-the-occurrence-of-an-event-over-multiple-attempts-1-chance-over-10-tries)
          object$shedDay <- shedLength         #enter shed and calculate values
          object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
          object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
          object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
          #15% of animals will enter this 2nd shed
        }}

###this chunk determines if an animal will undergo its 2nd shed during the 2nd window

    }else{NULL}                     #if the animal is already in shed and is not a neonate; do nothing
  # print("shedrules")
  }
#1 time hard code
# shedRules <- function(object){
#   ###this chunk determines if an animal qualifies for its neonatal shed
#   if(i==1 & fi == 1){                     #if the animal has just been born, place it into shed (nn shed)
#     object$shedDay <- shedLength                  #start a shed    #and calculate shed values
#     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#   }                                              #do nothing if not the first day of the first year
#   ###this chunk determines if an animal will undergo a shed in the first window
#   if(object$shedDay <= 0 & (fi >=2 | (fi ==1 & i > 90))){        #if an animal is not in shed currently and its either 1) in its 1st full summer (year 2) OR 2) past its winter of its first year (arbitrarily 90 days post parturition here);
#     if(i == 279){                    #if an animal is not yet in shed by day 279 force it into shed.
#       shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#       if(shedLiklihood <= 1){  #if the generated number is less than a threshold, say 70 percent of anmials shed in this window (arbitrary for now)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#       }}
# 
#     # if(i >= 259 & i < 279){        #if the date is between may 10 and june 1 (+/- 10 days around june 10th mean date of first shed obs. in 2 shedders minus 20 (days between shed start and observable cues), CM personal data)
#     #   shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#     #   if(shedLiklihood < 0.50){  #if the generated number is less than a threshold, 80 percent of animals shed in this window (arbitrary for now)
#     #     object$shedDay <- shedLength         #enter shed and calculate values
#     #     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     #     object$Eb <- (12.386)*(object$Mass ^ 0.9099812) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #     object$Er <- (30.226)*(object$Mass ^ 0.2210836) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #   }
#     # }
#     ###this chunk determines if an animal will shed in the 2nd window
#     # if(i >= 322 & i < 342){                   #if the date is between July 11 and july 31 (+/- 10 days around August 10th mean date of first shed obs. in 2 shedders minus 20 (days between shed start and observable cues), CM personal data)
#     #   #if its not in shed in the 2nd window, give animals a chance to shed
#     #   shedChance <- runif(1,0,1)                      #draw a random number between 0 and 1
#     #   if(shedChance <= 0.0077091277){                     #animals shed 1.15 times per year (males 1.13, femlaes 1.18), meaning average 15% of animals shed 2x per year, if this animal is one of those 15 percent (by random number), it will shed a 2nd time
#     #     #since probability stacks over the 20 days, I set the total chance (1-(1-p)^21) to 0.15, giving a chance, per draw of 0.007....  see math in yellow notebook or (https://www.quora.com/How-would-one-calculate-the-probability-of-the-occurrence-of-an-event-over-multiple-attempts-1-chance-over-10-tries)
#     #     object$shedDay <- shedLength         #enter shed and calculate values
#     #     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     #     object$Eb <- (12.386)*(object$Mass ^ 0.9099812) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #     object$Er <- (30.226)*(object$Mass ^ 0.2210836) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #     #15% of animals will enter this 2nd shed
#     #   }}
#     #
#     ###this chunk determines if an animal will undergo its 2nd shed during the 2nd window
# 
#   }else{NULL}                     #if the animal is already in shed and is not a neonate; do nothing
#   # print("shedrules")
# }
#2 times hard code
# shedRules <- function(object){
#   ###this chunk determines if an animal qualifies for its neonatal shed
#   if(i==1 & fi == 1){                     #if the animal has just been born, place it into shed (nn shed)
#     object$shedDay <- shedLength                  #start a shed    #and calculate shed values
#     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#   }                                              #do nothing if not the first day of the first year
#   ###this chunk determines if an animal will undergo a shed in the first window
#   if(object$shedDay <= 0 & (fi >=2 | (fi ==1 & i > 90))){        #if an animal is not in shed currently and its either 1) in its 1st full summer (year 2) OR 2) past its winter of its first year (arbitrarily 90 days post parturition here);
#     # if(i == 279){                    #if an animal is not yet in shed by day 279 force it into shed.
#     #   shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#     #   if(shedLiklihood < 0){  #if the generated number is less than a threshold, say 70 percent of anmials shed in this window (arbitrary for now)
#     #     object$shedDay <- shedLength         #enter shed and calculate values
#     #     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     #     object$Eb <- (12.386)*(object$Mass ^ 0.9099812) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #     object$Er <- (30.226)*(object$Mass ^ 0.2210836) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #   }}
# 
#     if(i == 269){        #if the date is  may 20  (+/- 10 days around june 10th mean date of first shed obs. in 2 shedders minus 20 (days between shed start and observable cues), CM personal data)
#       shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#       if(shedLiklihood <= 1){  #if the generated number is less than a threshold, 80 percent of animals shed in this window (arbitrary for now)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#       }
#     }
#     ###this chunk determines if an animal will shed in the 2nd window
#     if(i == 332){                   #if the date is  July 21 ----- and july 31 (+/- 10 days around August 10th mean date of first shed obs. in 2 shedders minus 20 (days between shed start and observable cues), CM personal data)
#       #if its not in shed in the 2nd window, give animals a chance to shed
#       shedChance <- runif(1,0,1)                      #draw a random number between 0 and 1
#       if(shedChance <= 1){                     #animals shed 1.15 times per year (males 1.13, femlaes 1.18), meaning average 15% of animals shed 2x per year, if this animal is one of those 15 percent (by random number), it will shed a 2nd time
#         #since probability stacks over the 20 days, I set the total chance (1-(1-p)^21) to 0.15, giving a chance, per draw of 0.007....  see math in yellow notebook or (https://www.quora.com/How-would-one-calculate-the-probability-of-the-occurrence-of-an-event-over-multiple-attempts-1-chance-over-10-tries)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         #15% of animals will enter this 2nd shed
#       }}
# 
#     ###this chunk determines if an animal will undergo its 2nd shed during the 2nd window
# 
#   }else{NULL}                     #if the animal is already in shed and is not a neonate; do nothing
#   # print("shedrules")
# }
# #3 times hard code #sheds take 28 days, summer is 182 days, 182 - (3*28)=98.  98/3 = 32  so, periods of 32 days w/o shed then 28 day sheds in 3 cycles throughout the summer.  they will shed at the start of summer, then monthly after that.
# shedRules <- function(object){
#   ###this chunk determines if an animal qualifies for its neonatal shed
#   if(i==1 & fi == 1){                     #if the animal has just been born, place it into shed (nn shed)
#     object$shedDay <- shedLength                  #start a shed    #and calculate shed values
#     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#   }                                              #do nothing if not the first day of the first year
#   ###this chunk determines if an animal will undergo a shed in the first window
#   if(object$shedDay <= 0 & (fi >=2 | (fi ==1 & i > 90))){        #if an animal is not in shed currently and its either 1) in its 1st full summer (year 2) OR 2) past its winter of its first year (arbitrarily 90 days post parturition here);
#     # if(i == 279){                    #if an animal is not yet in shed by day 279 force it into shed.
#     #   shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#     #   if(shedLiklihood < 0.70){  #if the generated number is less than a threshold, say 70 percent of anmials shed in this window (arbitrary for now)
#     #     object$shedDay <- shedLength         #enter shed and calculate values
#     #     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     #     object$Eb <- (12.386)*(object$Mass ^ 0.9099812) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #     object$Er <- (30.226)*(object$Mass ^ 0.2210836) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     #   }}
#     #
#     if(i == 216){        #if the date is the first day of summer, shed
#       shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#       if(shedLiklihood <=1){  #if the generated number is less than a threshold, 80 percent of animals shed in this window (arbitrary for now)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#       }
#     }
#     if(i == 276){        #if the date is 32 days after the last shed, shed
#       shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#       if(shedLiklihood <=1){  #if the generated number is less than a threshold, 80 percent of animals shed in this window (arbitrary for now)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#       }
#     }
#     ###this chunk determines if an animal will shed in the 2nd window
#     if(i == 336){                   #if the date is 32 days after the lst shed, shed
#       #if its not in shed in the 2nd window, give animals a chance to shed
#       shedChance <- runif(1,0,1)                      #draw a random number between 0 and 1
#       if(shedChance <= 1){                     #animals shed 1.15 times per year (males 1.13, femlaes 1.18), meaning average 15% of animals shed 2x per year, if this animal is one of those 15 percent (by random number), it will shed a 2nd time
#         #since probability stacks over the 20 days, I set the total chance (1-(1-p)^21) to 0.15, giving a chance, per draw of 0.007....  see math in yellow notebook or (https://www.quora.com/How-would-one-calculate-the-probability-of-the-occurrence-of-an-event-over-multiple-attempts-1-chance-over-10-tries)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         #15% of animals will enter this 2nd shed
#       }}
# 
#     ###this chunk determines if an animal will undergo its 2nd shed during the 2nd window
# 
#   }else{NULL}                     #if the animal is already in shed and is not a neonate; do nothing
#   # print("shedrules")
# }
# #continuous hard code
# shedRules <- function(object){
#   ###this chunk determines if an animal qualifies for its neonatal shed
#   if(i==1 & fi == 1){                     #if the animal has just been born, place it into shed (nn shed)
#     object$shedDay <- shedLength                  #start a shed    #and calculate shed values
#     object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#     object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#     object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#   }                                              #do nothing if not the first day of the first year
#   ###this chunk determines if an animal will undergo a shed in the first window
#   if(object$shedDay <= 0 & (fi >=2 | (fi ==1 & i > 90)) & environment$season != "Winter"){        #if an animal is not in shed, but is out of its first summer of birth,  enter shed
#     shedLiklihood <- runif(1,0,1) #select a random number between 0 and 1
#       if(shedLiklihood <= 1){  #if the generated number is less than a threshold, say 70 percent of anmials shed in this window (arbitrary for now)
#         object$shedDay <- shedLength         #enter shed and calculate values
#         object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
#         object$Eb <- (15.21745)*(object$Mass ^ 0.8802506) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#         object$Er <- (36.40273)*(object$Mass ^ 0.1946641) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
#       }
#   }else{NULL}                     #if the animal is already in shed and is not a neonate; do nothing
#   # print("shedrules")
# }


shedCost <- function(object){
  if(object$shedDay > 0){                      #if an animal is in shed;
    shedEday <- (object$Eb + object$Es)/(shedLength)                               #shed cost per day includes the total cost of synth and energy lost depostited as skin divided by 28.  simplistic but semiaccurate.  likely is more normal distribution shaped
    storePerDay <- shedEday / growthEff                #the amount of energy required to make that amount of material is subject to the conversion from storage to tissue.
    if(object$STS <= storePerDay){                     #if the amount required for making skin exceeds STS,
      object$LTS <- object$LTS - (storePerDay - object$STS); #lower LTS by the difference between the 2,
      object$STS <- 0                              #and set STS to 0
      object$shedDay <- object$shedDay - 1              #reduce the shed counter
    }else{                                            #if its less, take it away from STS
      object$STS <- object$STS - storePerDay      
      object$shedDay <- object$shedDay - 1              #reduce the shed counter
    };
    if(object$shedDay == 0){                   #if this is the last day before shedding,
      if(object$STS <= object$Er){             #if sts is less than er
        object$LTS <- object$LTS - (object$Er - object$STS)      #subtract the differnce from LTS
        object$STS <- 0                                   #put STS at 0
      }else{                                        #if Er can all be covered by STS
      object$STS <- object$STS - object$Er                  #remove the physical cost of shedding from STS
      }
      object$shedEnergy <- object$shedEnergy + (object$Es + object$Eb + object$Er)
    }else{NULL}}                               #if its not the last day of shed, add no Er cost
    }


#constant shedding conditions
# shedRules <- function(object){
# if(object$shedDay == 0 & environment$season != "Winter"){                      #if the animal is not in shed and its the active season, put it in shed
# object$shedDay <- shedLength                  #start a shed    #and calculate shed values
# object$Es <- (0.258) * (object$Mass ^ 0.8795)                 #energetic content of skin (kJ) (includes scaling for body mass) from Blem and Zimmerman, 1986 - antilog transformed from the eq. presented in the paper
# object$Eb <- (12.386)*(object$Mass ^ 0.9099812) * 0.02742     #from ch3, mass scaling for biosynthesis, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
# object$Er <- (30.226)*(object$Mass ^ 0.2210836) * 0.02742     #from ch3, mass scaling for physical removal, converted from mL CO2 to KJ by 27.42 joules/mL CO2 following Gessaman and Nagy, 1988
# }
# }
# 




###############################
#Methods for Environment
###############################


#set season function for environment
setSeason <- function(object){                                        #set season will be used to set the season of the environment object
  #  if(object$day >= 1 & object$day <= 80){object$season <- "Winter"}                   #winter part A are jan1-march21
  #  if(object$day >= 81 & object$day <= 171){object$season <- "Spring"}                 #spring march 22-june20
  #  if(object$day >= 172 & object$day <= 266){object$season <- "Summer"}                #summers june21-sept23
  #  if(object$day >= 267 & object$day <= 355){object$season <- "Autumn"}                #falls sept 24 - dec 21
  #  if(object$day >= 356 & object$day <= 365){object$season <- "Winter"} else {NULL}    #winter part B dec 22-dec 31
  #}
  #alternate day scheme for season to start on 28 august based on seasons from Beaupre 2001.
  if(object$day >= 1 & object$day <= 33){object$season <- "Summer"}                   #winter part A are jan1-march21
  if(object$day >= 34 & object$day <= 64){object$season <- "Autumn"}                 #spring march 22-june20
  if(object$day >= 65 & object$day <= 184){object$season <- "Winter"}                #summers june21-sept23
  if(object$day >= 185 & object$day <= 215){object$season <- "Spring"}                #falls sept 24 - dec 21
  if(object$day >= 216 & object$day <= 365){object$season <- "Summer"} else {NULL}    #winter part B dec 22-dec 31
}


#alternate season set to change length of the active season
#here, winter summer extends and 3xtra month on each side, spring and fall same lenght, winter 2 months shorter
# setSeason <- function(object){
# if(object$day >= 1 & object$day <= 63){object$season <- "Summer"}                   #winter part A are jan1-march21
# if(object$day >= 64 & object$day <= 94){object$season <- "Autumn"}                 #spring march 22-june20
# if(object$day >= 95 & object$day <= 154){object$season <- "Winter"}                #summers june21-sept23
# if(object$day >= 155 & object$day <= 185){object$season <- "Spring"}                #falls sept 24 - dec 21
# if(object$day >= 216 & object$day <= 365){object$season <- "Summer"} else {NULL}    #winter part B dec 22-dec 31
# }
# 




#####
#this temp method will only be necessary if we loop by hours, right now they are hard set
#setTemp function encased below v
#####
#this current setTemp just sets the avergae daily per season for use in calculating passage time.
setTemp <- function(object){                #set the average seasonal temp to calculate passage time (DT) by season and mass.
  if(object$season == "Winter"){object$temperature <- c(10, 10, 10, 10, 10, 10)}
  if(object$season == "Spring"){object$temperature <- c(25, 20, 17, 17, 22, 23)}
  if(object$season == "Summer"){object$temperature <- c(32, 30, 28, 25, 28, 32)}
  if(object$season == "Autumn"){object$temperature <- c(25, 22, 18, 16, 20, 23)}
  #these average numbers are from multiplying the degree's by the length of time block, summing them, then dividing by 24 (beaupre 2001).
  
  #set temp function for the environment.  Temps is based on time of day and season.
  #setTemp <- function(object){                                        #set temperatures based on season and time of day from table in the back of beaupre 2001.
  #  if(object$season == "Winter"){                                    # Winters are 10 degrees all day, in hibernacula                   
  #   if(object$timeofDay >= 13 & object$timeofDay <20){
  #     object$temperature <- 10}
  #   if(object$timeofDay >= 20 & object$timeofDay <24){
  #     object$temperature <- 10}
  #   if(object$timeofDay >= 0 & object$timeofDay <4){
  #     object$temperature <- 10}
  #   if(object$timeofDay >= 4 & object$timeofDay <8){
  #     object$temperature <- 10}
  #   if(object$timeofDay >= 8 & object$timeofDay <13){
  #     object$temperature <- 10}else{print("time not set right")}}
  # if(object$season == "Spring"){                                    #in spring temps vary         
  #   if(object$timeofDay >= 13 & object$timeofDay <20){
  #     object$temperature <- 25}
  #   if(object$timeofDay >= 20 & object$timeofDay <24){
  #     object$temperature <- 20}
  #   if(object$timeofDay >= 0 & object$timeofDay <4){
  #     object$temperature <- 17}
  #   if(object$timeofDay >= 4 & object$timeofDay <8){
  #     object$temperature <- 17}
  #   if(object$timeofDay >= 8 & object$timeofDay <13){
  #     object$temperature <- 22}else{print("time not set right")}}
  # if(object$season == "Summer"){                                     #in summer more variable
  #   if(object$timeofDay >= 13 & object$timeofDay <20){
  #     object$temperature <- 32}
  #   if(object$timeofDay >= 20 & object$timeofDay <24){
  #     object$temperature <- 30}
  #   if(object$timeofDay >= 0 & object$timeofDay <4){
  #     object$temperature <- 28}
  #   if(object$timeofDay >= 4 & object$timeofDay <8){
  #     object$temperature <- 25}
  #   if(object$timeofDay >= 8 & object$timeofDay <13){
  #     object$temperature <- 28}else{print("time not set right")}}
  # if(object$season == "Autumn"){                                     #in autumn similar to spring
  #   if(object$timeofDay >= 13 & object$timeofDay <20){
  #     object$temperature <- 25}
  #   if(object$timeofDay >= 20 & object$timeofDay <24){
  #     object$temperature <- 22}
  #   if(object$timeofDay >= 0 & object$timeofDay <4){
  #     object$temperature <- 18}
  #   if(object$timeofDay >= 4 & object$timeofDay <8){
  #     object$temperature <- 16}
  #   if(object$timeofDay >= 8 & object$timeofDay <13){
  #     object$temperature <- 20}else{print("time not set right")}}
  # else{NULL}
  #}
}


#add these as seaosnal vectors that feed in to the metabolism time slots (6 slots)
#as setTemp(environment) :: temp <- vector of temps :: if summer vector(25, 28, 33, 31, 28, 22) :: metabolism 0.333 * W ^3.5 * (3.44((temps[1])))




##################################
#Methods for Population growth and reproduction
##################################

#add neonates to the population.
populationGrowth <- function(object){
  if(object$yearRepro == T & i == 365){                   #if the animal reprod this year and its the day of parturition; add animals to the years pool
    neonates <- as.integer(object$nnProducedYearly)        #make a new variable alled neonates that is the length of the number of neonanes
   # print(class(neonates))
    for(bar in 1:neonates){
      chanceSex <- runif(1,0,1)         #draw a random number to determine male or female
      if(chanceSex > 0.5){              #make a male snake
        newbb <- maleSnake$new(1)
        newSnakes[[UUIDgenerate()]] <- c(newSnakes[[UUIDgenerate()]], newbb)  #save the animal as a unique animal in the population
      }else{                            #make a female
        newbb <- femaleSnake$new(1)
        newSnakes[[UUIDgenerate()]] <- c(newSnakes[[UUIDgenerate()]], newbb)  #save the animal as a unique animal in the population
      }
    }
    object$offspring <- newSnakes          #save the output to the animal
   }else{NULL}                            #do nothing if the animal didnt reproduce or it is not the last day of the year
}

#death model - using vital rate data to shape the model    :: from Olson et al., 2015 for now...
populationDeath<- function(object){
  individualAge <- (fi - object$birthYear)+1        #as long as fi is year, current year-birth year = age. must add 1, value cannot be 0
  survivalProb <- runif(1,0,1)                  #draw a random number to test for this years survival
  if(individualAge == 1){                         #if the individual is in its first year, its a neonate
    if(survivalProb >= neonateSurvivorship){                     #65% survival in neonates (olson et al, table)
      object$Dead <- T
      object$Mass <- 0
  #    print("die by the stats")
    }
}
  if(object$sexMature == F & individualAge >= 2){                      #immature animals that are not neonates - so everyone who isnt a baby or an adult
    if(survivalProb >= JuvenileSurvivorship){
      object$Dead <- T
      object$Mass <- 0
   #   print("die by the stats")
    }
    }
  if(object$sexMature == T & object$sex == "Male"){                   #adult male survival
    if(survivalProb >= MaleSurvivorship){                                           #median value from bill brown via olson table
      object$Dead <- T
      object$Mass <- 0
   #   print("die by the stats")
    }
  }
  if(object$sexMature == T & object$sex == "Female"){                 #adult female survival
    if(survivalProb >= FemaleSurvivorship){                                          #median value form bill brown via olson table
      object$Dead <- T
      object$Mass <- 0
      object$LTS <- 0
   #   print("die by the stats")
    }
  }
  }





#########################
#Model Test Build
#########################


#test
#create new environment
environment <- Environment$new(1)  #starts day 1 

#create a new s1, brand new male
# s1 <- maleSnake$new(1)
#set time limit
days <- 365                #set number of "days"
years <- 15
numInd <- 5000                #set highest number of indidivuals
MFS <- 0.0875

# #set foraging parameters
# forageMean <- 0.065              #set the mean value for the draw of foraging success from a normal dist.
# forageSD <- 0.02                #set standard deviation for the above.

#build the population
numMales   <- 5000
numFemales <- 0

#create blank lists to allow the addition of animals to thhe population
population <- tibble::lst()                   #generate an empty list to hold the initial population
newSnakes <- tibble::lst()                    #create an empty list to hold the new animals.
newSnakes3 <- tibble::lst()                   #creat an empty list
cemetary <- tibble::lst()                     #empty list to hold the dead animals and remove them from the population

#initialize the population
for(animal in 1:(numMales)){
  newanimal <- maleSnake$new(1)
  population[[UUIDgenerate()]] <- c(population[[UUIDgenerate()]], newanimal)  #this works, saves the new animal under a key value in the population
  }         #generate males and add them to the population

# for(animal in 1:numFemales){
#   newanimal <- femaleSnake$new(1)
#   population[[UUIDgenerate()]] <- c(population[[UUIDgenerate()]], newanimal)
#   }       #generate females and add them to the population


#preload programming - start some of the animals off at a certain LTS, mass, length, age etc
# minPreload <- 250              #the minimum body size of animals at initialization
# maxPreload <- 250              #the maximum body size of animals at initialization
# 
# for(yi in population){                     #for each list of 1 in the population
#   for(yoo in unlist(yi)){                  #unlist that list and act upon each individual object
#     preloadMass <- runif(1, minPreload, maxPreload)        #draw from a uniform distribtuion to determine the mass at start
#     yoo$Mass <- preloadMass                 #set an initial mass
#     yoo$SVL <- (11.73)*(yoo$Mass^0.324)     #set the SVL given average mass/svl relationships
#     yoo$LTS <- 377                          #37.7 kJ per g of fat, start them with 10g of fat.
#     yoo$birthYear <- -8                     #set the age, negative means they were born before the start of the simulation
#     if(yoo$SVL >= maturationSize){
#     yoo$sexMature <- T                      #set sex maturity if SVL exceeds the set sexmaturation size at intialization
#     }
#     }
# }



#create empty matrices to receive output (years as rows, days as columns)
foodingut <- array(0, dim=c(years, days, numInd))   #how much food is in the gut (proportion)
LTstore <- array(0, dim=c(years, days, numInd))     #Long term stores
svl <- array(0, dim=c(years, days, numInd))         #Snout Vent Length
mass <- array(0, dim=c(years, days, numInd))        #Mass
shedcount <- array(0, dim=c(years, days, numInd))
totalMetabLife <- array(0, dim=c(years, days, numInd))  #number of neonates produced
totalEB <- array(0, dim=c(years, days, numInd))  #number of neonates produced
shedEnergy <- array(0, dim=c(years, days, numInd))  #kJ devoted to ecdysis


#create blank population tables
#vector of total population
popAlive <- rep(0, years)
popstatistics <- matrix(0,years, 5) %>% as.data.frame(.) 
colnames(popstatistics) <-  (c("Year", "nJuv", "nMale", "nFemale", "DeadthisYear"))




############
#Start model
############
for(fi in 1:years){                  #repeat the day series for x number of years (defined above)
          for(i in 1:days){                    #for each day;
          environment$day <- i               #set the day in the environment by the loop: this means that the day will always reset to 1 at start!
          setSeason(environment)             #use the function setSeason to set the season
          setTemp(environment)               #set the seasonal average temp to calculate passage time (DT)
    for(foo in population){            #loop over individuals in the population
          for(individual in unlist(foo)){            #since individuals are stored as key-value's, popualtion is actually a list of names with characteristics, this 4th loop is necessary to get down to the individual R6 object.  
  reserveLTS(individual);                    #set the seasonal value for reserve energy.  
  alloRules(individual);                     #determine time/energy allocation based on sex and maturity
  canGrowRepro(individual);                  #set proportion of energy to grow vs repro depending on season and sex/maturity
  if(i == 1){individual$nnProducedYearly <- 0} #reset the yearly production of neonates to 0 on day 1
#remove old animals from the population
if(lifeSpan < (fi - individual$birthYear)){      #if total lifespan (in years) is less than the animals' current age in years, 
 individual$LTS <- 0                          #"Kill" the animal
 individual$Mass <- 0
 individual$Dead <- T                         #make the individual Dead
  next           #skip to the next individual
  };
  
# if(i==1 & individual$Dead == F){                                    #on the first day of each year, for live animals
#   populationDeath(individual)                  #run the probabliity of yearly survival program to model population vital rates
# }else{NULL};
  
if(individual$Dead == T){                     #if this indiviual is dead,
  next                                          #skip to the next individual
}else{NULL};
  
  
#flow of events 
  
if(individual$LTS >= 0){                     #if the snake has energy;
 if(individual$vitel == T & individual$daysVitel > 395){NULL}#print(paste0("the module run is A"))}       #if the animal is 2ndary vitellogenic and in the summer (starts day 215, fertilized at 185, so 395 days after start of primary vitell), do noting
  else{                                  #if the animal isn't vitellogenic or it is but its before the fertilizaiotn date, behave as usual
  if(individual$shedDay >= (ceiling(shedLength*propShed)) | individual$shedDay == 0){               #if the animal is not in shed, or in the first two weeks (counter runs down from 28), it can eat and digest
    if(individual$FinG == 0){                  #if the animal has no individuald,
      forage(individual)                       #forage #forage sets the digest time too.
      }else{digest(individual)
        };               #if FinG > 0, it ate;#digest the individuald
  }else{
    if(individual$FinG > 0){                 #if in the last two weeks of shed, continue to digest but don't forage
      digest(individual)
      }
    else{NULL}                        #if no individuald in gut dont digest or forage.
  }                          
 }
  metab(individual);                         #pay daily metabolic costs

 shedRules(individual);                     #assess whether the animal is in shed (nn, 1st, 2nd)
 shedCost(individual);                      #pay any relevent costs for shedding

  
  if(individual$LTS > individual$safetyStore & individual$daysVitel <= 366){       #if the snake has enough energy in reserve to survive based on season;
    growth(individual)                       #then the animal will use energy to grow
    reproduction(individual)                 #allow for allocation to reproduction
    }else{NULL};                     #if the safety store is not met, do not grow.
  if(individual$daysVitel > 366){                 #if the animal is 2ndary vitellogenic
    #growth(individual)                      #gravid females will not grow
    reproduction(individual)                 #allow for allocation to reproduction
    }else{NULL};                     #close out the if statement

  individual$LTS <- individual$LTS + individual$STS;         #add any remaining STS to LTS
  individual$STS <- 0;                       #reset STS to 0 for the day.
  sexMaturation(individual)                 #check if the animal grows above 70cm; if so change to mature
  #run the neonate program to add animals to the population
 populationGrowth(individual)
  #print(paste0("the length of offspring after we reenter the loop is ", length(individual$offspring)))
  #implement population growth by adding individuals on the parturition day, and appending them to the population.
  if(individual$nnProducedYearly >= 1 & i==365){                                 #if there are snakes that have been born
    offspringtally <- tibble::lst(individual$offspring)    #save the animal's offspring to a new list
    newSnakes2 <- unlist(offspringtally)                   #rip the nested list out into single animals
    for(baby in newSnakes2){                               #for each baby born in the new list; 
      baby$birthYear <- fi                                 #log the birth year.
      newSnakes3[[UUIDgenerate()]] <- c(newSnakes3[[UUIDgenerate()]], baby)   #populate the newsnakes list with these new animals
    }
    population <- c(population, newSnakes3)   #add the new snakes to the population (this should happen at the end of day 365, new snakes should be born on day 1 of a year)
    individual$offspring <- tibble::lst()             #make the animals offpsring blank again
    newSnakes <- tibble::lst()                        #clear out the new snakes list
    newSnakes2 <- NULL                                #clear out newsnakes2 just in case.
    newSnakes3 <- tibble::lst()                       #clear out the newsnakes3 list too.
    individual$nnProducedYearly <- 0          #reset repro number to 0 AFTER individuals are added to the pop.
    individual$yearRepro <- F                 #reset repro to false at the end of the year
   }
  }else{                               #if the snake has no energy, it is dead.
    individual$Mass <- 0                     #set svl and mass to 0
    individual$Dead <- T                     #make individual Dead.
    next                            #move on to the next individual, individuals should leave 0's for all fields after they die, rather than a prolonged plateau
    };                       
#save outputs at each iteration at the position in the matrix column= day, row = year

            
  #           
  # foodingut[fi, i, (which(sapply(population, identical, foo)))] <- individual$FinG; #individuald in the gut, the complicated sapply term allows me to index the position of the list of individuals(individual) within the list of individuals (foo)
  # LTstore[fi, i, (which(sapply(population, identical, foo)))] <- individual$LTS;    #long term stores
  # svl[fi, i, (which(sapply(population, identical, foo)))] <- individual$SVL;        #svl
  # mass[fi, i, (which(sapply(population, identical, foo)))] <- individual$Mass       #mass
  # shedcount[fi, i, (which(sapply(population, identical, foo)))] <- individual$shedDay
  # if(individual$sex == "Female"){nnProduced[fi, i, (which(sapply(population, identical, foo)))] <- individual$nnProd}else{NULL} #neonates produced
  # if(individual$sex == "Female"){nnProducedYearly[fi, i, (which(sapply(population, identical, foo)))] <- individual$nnProducedYearly}else{NULL} #neonates produced
  shedEnergy[fi, i, (which(sapply(population, identical, foo)))] <- individual$shedEnergy
  totalEB[fi, i, (which(sapply(population, identical, foo)))] <- individual$totalEB
  totalMetabLife[fi, i, (which(sapply(population, identical, foo)))] <- individual$totalMetabLife
  


  #if animals died during the day, move them to the dead list instead of the population
  #havn't figured this out yet.
  }                                            #close indiviudals
  }                                           #close foo in population
#add a day to the environment, so we loop through seasons
environment$day <- environment$day + 1   #advance the day by 1
}                                             #close day loop
 popAlive[fi] <- length(population)             #
 popstatistics[fi,1] <- fi                      #populate the year of the simulation into the popstatistics data frame
 for(z in unlist(population)){      #for each individual in the population, lets count things and add them to thir corresponind cell in the data summary
                   if(z$LTS > 0){                  #count only life animals
                             if(z$sexMature ==F){            #for immature animals
                               popstatistics[fi,2] <- (popstatistics[fi,2]) + 1
                             }else{                          #for mature animals
                               if(z$sex == "Male"){             #count males
                                 popstatistics[fi,3] <- (popstatistics[fi,3]) + 1
                               }else{                           #count females
                                 popstatistics[fi,4] <- (popstatistics[fi,4]) + 1
                               }
                             }
                           }else{      #anyone with no LTS is dead.
                             popstatistics[fi,5] <- (popstatistics[fi,5]) + 1
                           }
 }              #save the number of individuals in each repro class each year.
 
             }                                             #close year loop
################
#End Test Model
################


popstatistics


























