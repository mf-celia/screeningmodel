##############################################################################
## A Simplified Model of the Cost-Effectiveness of Screening:               ##
## An Open-Source Teaching and Research Tool Coded in R                     ##
##############################################################################

# Authors: 
# James O'Mahony ,  PhD (1) 
# Yi-Shu Lin,  MSc,  MBA (1)	

# 1 Trinity College Dublin, Ireland

##################################################################################
# Please cite our publications when using this code
# O'Mahony JF. Simplified Model of the Cost-Effectiveness of Screening in R: a Teaching and Research Tool. 17th Biennial European Meeting of the Society for Medical Decision Making Leiden, the Netherlands, June 10-12, 2018. (2018). Medical Decision Making, 38(6), E372-E610. https://doi.org/10.1177/0272989X18793413

##################################################################################
# Reference for ICER calculation:
# Jalal H, et al. An Overview of R in Health Decision Sciences. Med. Decis. Making. 2017; 37(3): 735-746. 
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. Med. Decis. Making. 2018; (in press). 

##################################################################################
# See GitHub for code updates
# https://github.com/yishu-lin/screeningmodel.git

##################################### Background Setting #########################################
rm(list = ls())          # Delete everything that is in R's memory
time <- Sys.time()       # save system time 
set.seed(1)              # set a seed to be able to reproduce the same results

# Set the working file directory path. Paste your working directory within the quotes
# Ensure that the backslashes \ are changed to forward slashes /
# You can put as many file paths in as you like here,  R will take the last one that worked
# The try function simply masks error reporting: it will work if one of the directories attempted is correct
try(setwd("Your working directory 1"),  	silent = TRUE)        # Drive file path
try(setwd("Your working directory 2"),  	silent = TRUE)        # Drive file path
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))     # Drive file path

##################################### Model input #########################################
# Import the required data from the data sheets created by the Excel file
LifeTable			       <- data.frame(read.table(file = "InputFiles/LifeTable.txt",     		header = T))  # Interval_LifeTable; Survival
Incidence			       <- data.frame(read.table(file = "InputFiles/Incidence.txt",       	header = T))  # Age_Interval; Annual_Incidence; Cumulative_Failure
TestPerformance		   <- data.frame(read.table(file = "InputFiles/TestPerformance.txt",  header = T))  # Test; Sn; Sp
StageDuration	       <- data.frame(read.table(file = "InputFiles/StageDuration.txt",    header = T))  # Stage; Type; Scale; Shape
ScreenSchedule		   <- data.frame(read.table(file = "InputFiles/ScreenSchedule.txt",   header = T))  # ScheduleNumber; StartAge; StopAge; TestApplied1; Interval1
TreatmentSuccess	   <- data.frame(read.table(file = "InputFiles/TreatmentSuccess.txt", header = T))  # PreClinicalProbability; ClinicalProbability
DiscountRates		     <- data.frame(read.table(file = "InputFiles/DiscountRates.txt", 		header = T))  # Costs; Effects; DiscountYear
Costs				         <- data.frame(read.table(file = "InputFiles/Costs.txt", 				    header = T))  # PrimaryScreen; FollowUp; TreatmentScreenDetected; TreatmentClinical
DefineStages				 <- data.frame(read.table(file = "InputFiles/Stages.txt", 				  header = T)) 
source(file <- "InputFiles/Misc.txt")                                                                 # SampleSize; Threshold
source(file <- "InputFiles/Options.txt")                                                             

# Prepare the input table:
# the performance of five tests (including sensitivity and specificity): 10 parameters
# the discount rate for costs and effects,  and the discount year: 3 parameters
# the probability of treatment success: 2 parameters
# the costs of primary screening,  follow-up,  and treatment (for patients either detected by screens or not): 4 parameters
# the annual incidence rate by five year age group: 1 parameter but 21 variables to represent different groups
# Calculate the total number of parameters needed to run the simulation excluding the annual incidence rate
NumberofParameters <- 19
InputTable <- array(NA, dim=c(1, (NumberofParameters + nrow(Incidence) + (nrow(StageDuration) * 3))))
colnames(InputTable) <- c(paste("TestSensitivity", 1:length(TestPerformance[, "Sn"]), sep = ""), 
                          paste("TestSpecificity", 1:length(TestPerformance[, "Sp"]), sep = ""), 
                          "DiscountRateCosts", "DiscountRateEffects", "DiscountYear", 
                          "PreClinicalProbability", "ClinicalProbability", 
                          "CostofPrimaryScreen", "CostofFollowUp", "CostofTreatScreenDetected", "CostofTreatmentClinical", 
                          paste("Incidence", Incidence$Age_Interval, sep = "_"), 
                          paste("StageDuration", paste("stage", 1:length(StageDuration[, "Stage"]), sep = ""), "Type", sep = "_"), 
                          paste("StageDuration", paste("stage", 1:length(StageDuration[, "Stage"]), sep = ""), "Scale", sep = "_"), 
                          paste("StageDuration", paste("stage", 1:length(StageDuration[, "Stage"]), sep = ""), "Shape", sep = "_"))  # Set column names
# Import the data and list all the parameters
InputTable[1, ] <- c(as.numeric(TestPerformance[, "Sn"]), 
                     as.numeric(TestPerformance[, "Sp"]), 
                     as.numeric(DiscountRates), 
                     as.numeric(TreatmentSuccess), 
                     as.numeric(Costs), 
                     as.numeric(Incidence[, "Annual_Incidence"]), 
                     as.numeric(StageDuration[, "Type"]), 
                     as.numeric(StageDuration[, "Scale"]), 
                     as.numeric(StageDuration[, "Shape"])) # Combine all the parameters from different sources into one table

CurrentRun <- 1

##################################### Model structure #########################################
# Define number of health stages with the stage arrival 
Stages			     <- nrow(DefineStages)   # This model only features five states: (1)Disease Free; (2)Preclinical Disease; (3)Clinical Disease; (4)Cause-Specific Death; (5) Other-Cause Death.

##################################### Simulation: Natural History of Disease ###########################################
Outcomes <- array(NA, dim = c(SampleSize, 7)) # Create an array of the length of the sample size
# The intervening columns correspond to the arrival of the intermediate disease states
# This model studies the cohort with the same age, so here does not need to simulate the "Stage1Arrival"
colnames(Outcomes) <- c("PersonNumber", paste(DefineStages[2: nrow(DefineStages),"Name"]), "AllCauseDeath", "CureStatus") # Set column names
Outcomes[, "PersonNumber"] <- c(1:SampleSize)   # Set the first column to be the unique person-number for each individual

# Define the generic onset function which applies an age-specific probability of entering a specific stage
OnsetFunction <- function(x) {
  unlist(approx(probability, age, x, ties = max)[2], use.names = F) # Ties = max is required because of the possibility of multiple zero probabilities of disease at younger ages
}

CumulativeFailure	        <- 1
IncidenceTable	          <- data.frame(Incidence$Age_Interval,  as.numeric(InputTable[CurrentRun, paste("Incidence", Incidence$Age_Interval, sep = "_")]),  CumulativeFailure)
colnames(IncidenceTable)	<- c("AgeInterval", "AnnualIncidence", "CumulativeFailure")

for (c in c(2:nrow(IncidenceTable))) {
  IncidenceTable[c, "CumulativeFailure"]	<- (IncidenceTable[c-1, "CumulativeFailure"] - (IncidenceTable[c, "AnnualIncidence"] * (IncidenceTable[c, "AgeInterval"] - IncidenceTable[c-1, "AgeInterval"])))
}

# Define the specific values of probability and age for disease onset
probability	<-	IncidenceTable[, "CumulativeFailure"]
age			    <- 	IncidenceTable[, "AgeInterval"]
# The probability distribution need to be curtailed at the top and bottom for non-unique probability for the approx function to work as intended
age			    <- age[max(which(probability == max(probability))):length(probability)]            # Remove the lower bound values
probability	<- probability[max(which(probability == max(probability))):length(probability)]    # Remove the lower bound values
age		     	<- age[1:min(which(probability == min(probability)))]                              # Remove the upper bound values
probability	<- probability[1:min(which(probability == min(probability)))]                      # Remove the upper bound values
x	<- runif(SampleSize)                            # Create a vector of random values
Outcomes[, "OnsetAge"] <- OnsetFunction(x)        # Find the age of disease onset by applying the general function

# Define the specific values of probability and age for other cause death
# The intervals in the age-specific incidence are defined by those used in the life table
probability	<- LifeTable[, "Survival"]
age		    	<- LifeTable[, "Interval_LifeTable"]
# The probability distribution need to be curtailed at the top and bottom for non-unique probability for the approx function to work as intended
age			    <- age[max(which(probability == max(probability))):length(probability)]            # Remove the lower bound values
probability	<- probability[max(which(probability == max(probability))):length(probability)]    # Remove the lower bound values
age		     	<- age[1:min(which(probability == min(probability)))]                              # Remove the upper bound values
probability	<- probability[1:min(which(probability == min(probability)))]                      # Remove the upper bound values
x	<- runif(SampleSize)                            # Create a vector of random values
Outcomes[, "OtherCauseDeath"]	<- OnsetFunction(x) # Find the age at other-cause death by applying the generic onset function. 

# The model needs a loop here to go through the disease stages
for (Stage in 2:(Stages-2)) { # Apply the sojourn time to the stages, excluding "Disease Free", "Preclinical Disease" and "Other-Cause Death"
  # Retrieve the distribution type, scale and shape
  DurationType	<-	InputTable[CurrentRun, paste("StageDuration", paste("stage", Stage, sep = ""), "Type",  sep = "_")]
  DurationScale	<-	InputTable[CurrentRun, paste("StageDuration", paste("stage", Stage, sep = ""), "Scale", sep = "_")]
  DurationShape	<-	InputTable[CurrentRun, paste("StageDuration", paste("stage", Stage, sep = ""), "Shape", sep = "_")]
  
  if(DurationType == 1) {
    Duration <- rep(DurationScale, SampleSize) 																	     # set the preclinical distribution to be constant
  }
  if(DurationType == 2) {
    Duration <- -(log(1 - runif(SampleSize))) * DurationScale										     # set the preclinical distribution to be exponentially distributed
  }
  if(DurationType == 3) {
    Duration <- rweibull(SampleSize, shape = DurationShape, scale = DurationScale)   # set the preclinical distribution to be Weibull distributed    #Or you can use the code: Duration <- ((-log(1-runif(SampleSize)))^(1/DurationShape))*DurationScale
  }
  
  Outcomes[, Stage+1] <- Outcomes[, Stage] + Duration  # Find the end of the preclinical period by adding the onset to the duration
}

# This study assumes the maximum age of cause-specific death is 100. 
Outcomes[which(Outcomes[,"CauseSpecificDeath"]>100), "CauseSpecificDeath"] <- 100

# Find the all-cause death
Outcomes[, "AllCauseDeath"] <- pmin(Outcomes[, "CauseSpecificDeath"], Outcomes[, "OtherCauseDeath"], na.rm = TRUE)

# Censor the disease free, ie those who never develop disease within their lifetime
DiseaseFree   <- which((Outcomes[, "OnsetAge"] >= Outcomes[, "OtherCauseDeath"]) | (is.na(Outcomes[, "OnsetAge"]))) # Censor onset ages greater than death; first identify those who never develop disease
Outcomes_Sick <- Outcomes[-DiseaseFree, "PersonNumber"]

# Identify individuals presenting with clinical disease
ClinicalCases	<- Outcomes_Sick[which(Outcomes[Outcomes_Sick, "Clinical"] <= Outcomes[Outcomes_Sick, "OtherCauseDeath"])]

# Identify individuals cured successfully
Outcomes[, "CureStatus"]              <- runif(SampleSize)                              # Set the cure column equal to a random for these people
ClinicalTreatmentSuccess              <- InputTable[CurrentRun, "ClinicalProbability"]  # Retrieve the probability of successful treatment 
CuredCases                            <- ClinicalCases[which(Outcomes[ClinicalCases, "CureStatus"] <= ClinicalTreatmentSuccess)]  # Find treatment to be successful if probabilities less than probability of successful treatment
Outcomes[CuredCases, "AllCauseDeath"] <- pmax(Outcomes[CuredCases, "CauseSpecificDeath"], Outcomes[CuredCases, "OtherCauseDeath"], na.rm = TRUE)

##################################### Screening Strategy Generation ################################################################
# Retrieve the list of screening strategies
ScheduleNumbers	  	<- ScreenSchedule[, "ScheduleNumber"]
StartAges	      		<- ScreenSchedule[, "StartAge"]
StopAges		      	<- ScreenSchedule[, "StopAge"]
IntervalSwitchAges1	<- ScreenSchedule[, "IntervalSwitchAge1"]
IntervalSwitchAges2	<- ScreenSchedule[, "IntervalSwitchAge2"]
IntervalSwitchAges3	<- ScreenSchedule[, "IntervalSwitchAge3"]
Intervals1		    	<- ScreenSchedule[, "Interval1"]
Intervals2		     	<- ScreenSchedule[, "Interval2"]
Intervals3		    	<- ScreenSchedule[, "Interval3"]
Intervals4		    	<- ScreenSchedule[, "Interval4"]
TestSwitchAges	  	<- ScreenSchedule[, "ScreenSwitchAge"]
ScreenTests1	    	<- ScreenSchedule[, "TestApplied1"]
ScreenTests2        <- ScreenSchedule[, "TestApplied2"]

# Find number of strategies in the list of strategies to simulate
NumberOfStrategies  <- nrow(ScreenSchedule)

# Create Array to hold outcomes with column names: we may update the StrategyName to make something descriptive of the strategies
CostEffectivenessResults           <- array(NA, dim = c(0, 6))
colnames(CostEffectivenessResults) <- c("StrategyName", "UD_Effects", "UD_Costs", "D_Effects", "D_Costs", "ICER")
StrategyNames                      <- NULL   # Create blank object for strategy names

# Loop over the multiple strategies
for (Strategy in 1:NumberOfStrategies) {
  print(Strategy)
  print(NumberOfStrategies)
  # Retrieve the strategy-specific values
  ScheduleNumber		 <- ScheduleNumbers[Strategy]
  StartAge				   <- StartAges[Strategy]
  StopAge				     <- StopAges[Strategy]
  IntervalSwitchAge1 <- IntervalSwitchAges1[Strategy]
  IntervalSwitchAge2 <- IntervalSwitchAges2[Strategy]
  IntervalSwitchAge3 <- IntervalSwitchAges3[Strategy]
  Interval1				   <- Intervals1[Strategy]
  Interval2				   <- Intervals2[Strategy]
  Interval3				   <- Intervals3[Strategy]
  Interval4			   	 <- Intervals4[Strategy]
  TestSwitchAge	 	 	 <- TestSwitchAges[Strategy]
  ScreenTest1		   	 <- ScreenTests1[Strategy]
  ScreenTest2			   <- ScreenTests2[Strategy]
  
  # Paste together a descriptor of the screen strategy from features of the screen schedule
  StrategyName   <- paste(StartAge, "_", Interval1, "_", IntervalSwitchAge1, "_", Interval2, "_", 
                          IntervalSwitchAge2, "_", Interval3, "_", IntervalSwitchAge3, "_", Interval4, "_",
                          StopAge, "_", ScreenTest1, "_", TestSwitchAge, "_", ScreenTest2, sep = "")
  
  # Amend the screening schedule
  Intervals         <- 1         # Determine the number of intervals, and the default is one interval
  if (!is.na(IntervalSwitchAge1)) {Intervals = 2} else {}
  if (!is.na(IntervalSwitchAge2)) {Intervals = 3} else {}
  if (!is.na(IntervalSwitchAge3)) {Intervals = 4} else {}
  IntervalVector    <- c(Interval1, Interval2, Interval3, Interval4)  # Create a vector of the screening intervals
  ChangeAges        <-	sort(na.omit(c(StopAge, IntervalSwitchAge1, IntervalSwitchAge2, IntervalSwitchAge3)))  # Create a vector of ages at which screening changes
  Screens           <- NULL       # Create an empty vector to bind the sequence of screens
  for (Interval in 1:Intervals) { # Loops through the Intervals
    StopAge         <- ChangeAges[Interval]                                          # Redefine the StopAge as the Change age in this interval
    NumberOfScreens <- floor((StopAge - StartAge)/IntervalVector[Interval]) + 1      # Find the number of screens given rounding 
    StopAge         <- StartAge + (NumberOfScreens-1) * IntervalVector[Interval]     # Redefine the actual StopAge following rounding
    Screens         <- c(Screens, seq(StartAge, StopAge, IntervalVector[Interval]))  # Create the sequence of screens from the start and stop ages
    StartAge        <- max(Screens)                                                  # Update the final screen age as the Start Age for cases in which the interval changes
  }
  Screens           <- unique(Screens)  # Remove duplicates from the screen schedule
  NumberOfScreens   <- length(Screens)  # Update the number of screens
  
  # Create the screen schedule -an array with screening details  
  ScreenTable                   <- array(NA, dim = c(NumberOfScreens, 3))  # Create an array from the screen schedule
  colnames(ScreenTable)         <- c("ScreenNumber", "TestApplied", "ScreenAge")  # Name the columns
  ScreenTable[, "ScreenNumber"] <- c(1:NumberOfScreens) # Write in a list equivalent to the number of screens and the per round screening age
  ScreenTable[, "ScreenAge"]	  <- Screens
  ScreenTable[, "TestApplied"]  <- ScreenTest1          # Insert the appropriate test number
  ScreenTable[which(ScreenTable[, "ScreenAge"] >= TestSwitchAge), "TestApplied"] <- ScreenTest2
  
  # Now run screening
  # Make reference to screening protocol
  ScreenNumbers	<-	ScreenTable[, "ScreenNumber"]
  ScreenTests		<-	ScreenTable[, "TestApplied" ]
  ScreenAges		<-	ScreenTable[, "ScreenAge"   ]
  
  # Create the screen-adjusted outcomes from the natural history outcomes
  ScreenedOutcomes <- Outcomes

  # Create an array of the screen counts: This is an array recording the number of screens and consequent results
  ScreenCounts 	            	<- array(NA, dim=c(NumberOfScreens, 5))
  colnames(ScreenCounts)      <- c("ScreenAge", "N_Screens", "AllPostives", "TruePositives", "FalsePositives")  # Column Names for the screen counts
  ScreenCounts[, "ScreenAge"] <- ScreenAges  # Write in the screen ages
  
  for (ScreenNumber in c(ScreenNumbers)) {
    # Define the screen age
    ScreenAge <- ScreenAges[ScreenNumber]
    # Retrieve the test sensitivity and specificity 
    ScreenSn  <- switch(ScreenTests[ScreenNumber], InputTable[CurrentRun, "TestSensitivity1"], 
                                                   InputTable[CurrentRun, "TestSensitivity2"], 
                                                   InputTable[CurrentRun, "TestSensitivity3"], 
                                                   InputTable[CurrentRun, "TestSensitivity4"], 
                                                   InputTable[CurrentRun, "TestSensitivity5"])
    ScreenSp  <- switch(ScreenTests[ScreenNumber], InputTable[CurrentRun, "TestSpecificity1"], 
                                                   InputTable[CurrentRun, "TestSpecificity2"], 
                                                   InputTable[CurrentRun, "TestSpecificity3"], 
                                                   InputTable[CurrentRun, "TestSpecificity4"], 
                                                   InputTable[CurrentRun, "TestSpecificity5"])
    
    Alive          <- ScreenedOutcomes[, "AllCauseDeath"] >= ScreenAge
    NotDiagnosed   <- ScreenedOutcomes[, "Clinical"] >= ScreenAge
    NotDiagnosed[is.na(NotDiagnosed)] <- TRUE
    ScreenEligible <- Alive*NotDiagnosed
 
    # Identify those in the preclinical stage at the screen age: we don't need to exclude the screen eligible here
    AllPositives <- (which(ScreenedOutcomes[, "OnsetAge"] <= ScreenAge))[(which(ScreenedOutcomes[, "OnsetAge"] <= ScreenAge) %in% which(ScreenedOutcomes[, "Clinical"] >= ScreenAge))]
    
    # Identify the negatives as the complement of the positives from within the ScreenEligible set
    AllNegatives <- (which(ScreenEligible==1)[!(which(ScreenEligible==1) %in% AllPositives)])
    
    
    # Find the true positives by sampling without replacement over all positives in proportion to the test sensitivity
    TruePositives  <- sample(AllPositives, round(ScreenSn * length(AllPositives), 0), FALSE)
    # Find the false positives by sampling without replacement over the negatives in proportion to the test specificity
    FalsePositives <- sample(AllNegatives, round((1-ScreenSp) * length(AllNegatives), 0), FALSE)
    
    # It's useful to have a per screening round count of the number of positives,  true positives and false positives 
    ScreenCounts[ScreenNumber, 2:5] <- cbind(length(which(ScreenEligible==1)), length(AllPositives), length(TruePositives), length(FalsePositives))
    
    # Censor these successfully treated individuals and update all-cause death
    PreClinicalTreatmentSuccess <- InputTable[CurrentRun, "PreClinicalProbability"]  # Retrieve the probability of successful treatment following screen detection
    ScreenedCuredCases          <- TruePositives[which(ScreenedOutcomes[TruePositives, "CureStatus"] <= PreClinicalTreatmentSuccess)]
    
    # Now update the model to record the screen-adjusted outcomes
    ScreenedOutcomes[TruePositives, "Clinical"]           <- ScreenAge
    ScreenedOutcomes[ScreenedCuredCases, "AllCauseDeath"] <- pmax(ScreenedOutcomes[ScreenedCuredCases, "CauseSpecificDeath"], ScreenedOutcomes[ScreenedCuredCases, "OtherCauseDeath"], na.rm = TRUE)
  }
  
  ##################################### Effects and Costs of Screening Strategy ################################################################
  # Assess differences in effects
  # Retrieve discount rates for effects,  costs and the discount year
  DR_Effects   <- InputTable[CurrentRun, "DiscountRateEffects"]
  DR_Costs     <- InputTable[CurrentRun, "DiscountRateCosts"]
  DiscountYear <- InputTable[CurrentRun, "DiscountYear"]
  
  # Create a reference series of discount factors for 100 years
  DF_Effects <- (1 + DR_Effects) ^- (c(0:100))
  DF_Costs	 <- (1 + DR_Costs  ) ^- (c(0:100))
  
  # Now find the overall difference in LYG by comparing the two sets of outcomes
  # ScreenAdjusted <- which(!Outcomes[, "AllCauseDeath"] == ScreenedOutcomes[, "AllCauseDeath"])  # Identify those with lives extended by screening
  
  # Create an array for the difference in undiscounted and discounted health effects
  EffectsArray	         <- array(NA, dim=c(SampleSize, 4))
  colnames(EffectsArray) <- c("UD_NoScreen", "UD_WithScreen", "D_NoScreen", "D_WithScreen")  # Name the columns
  
  # Write in the undiscounted life years without and with screening respectively
  EffectsArray[, "UD_NoScreen"]	  <- Outcomes[, "AllCauseDeath"]
  EffectsArray[, "UD_WithScreen"]	<- ScreenedOutcomes[, "AllCauseDeath"]
  
  # Define the present value function
  PresentValue <- function(x) {
    A <- trunc(x, 0)
    B <- x-A
    sum(DF_Effects[1:A]) + B * DF_Effects[A + 1]
  }
  
  # Calculate the discounted life years with screening
  PresentValueFunction <-1
  for (i in c(2:100)) {
    PresentValueFunction <- c(PresentValueFunction, sum(DF_Effects[1:i]))
  }
  
  EffectsArray[, "D_WithScreen"] <- EffectsArray[, "UD_WithScreen"] - DiscountYear  # First adjust the discounted values for the discount year
  IntegerLYs <- trunc(EffectsArray[, "D_WithScreen"],0)
  DecimalLYs <- EffectsArray[, "D_WithScreen"]-trunc(EffectsArray[, "D_WithScreen"],0)

  # Sum up the undiscounted and discounted LYG
  UnDiscountedLYG <- sum(EffectsArray[, "UD_WithScreen"])
  DiscountedLYG   <- sum(PresentValueFunction[IntegerLYs])+sum(DecimalLYs*DF_Effects[IntegerLYs+1])
  
  # Assess differences in costs
  # Retrieve the costs parameters
  ScreenCost  		   <- InputTable[CurrentRun, "CostofPrimaryScreen"]
  TriageCost	  		 <- InputTable[CurrentRun, "CostofFollowUp"]
  EarlyTreatmentCost <- InputTable[CurrentRun, "CostofTreatScreenDetected"]
  LateTreatmentCost	 <- InputTable[CurrentRun, "CostofTreatmentClinical"]
  
  # Count the costs of treating clinical disease under screening
  # Count the treatment costs with screening
  ClinicalCases	  	<- which(ScreenedOutcomes[, "Clinical"] <= ScreenedOutcomes[, "AllCauseDeath"])  # Update the clinical cases post-screening
  ClinicalOnsetAges	<- round(ScreenedOutcomes[ClinicalCases, "Clinical"] - DiscountYear, 0)  # Identify the clinical ages and adjust for the discount year and round to integers
  DiscountWeight		<- DF_Costs[ClinicalOnsetAges]  # Find discount weight
  ScreenedUnDiscountTreatmentCosts <- sum(length(ClinicalCases) * LateTreatmentCost)
  ScreenedDiscountTreatmentCosts	 <- sum(DiscountWeight * LateTreatmentCost)  # Find the discounted treatment costs
  
  # First calculate the undiscounted amounts: the costs of primary screening, triage and treating preclinical disease under screening
  UD_ScreenCosts    <- sum(ScreenCounts[, "N_Screens"] * ScreenCost)  # The primary test cost applies to all screens
  UD_TriageCosts    <- sum((ScreenCounts[, "TruePositives"] + ScreenCounts[, "FalsePositives"]) * TriageCost)  # The cost of triage applies to all screen true and false positives
  UD_TreatmentCosts <- sum(ScreenCounts[, "TruePositives"] * EarlyTreatmentCost)  # The early treatment costs applies to all true positives (assuming a 100% specific triage test)
  
  # Now calculate the discounted values
  ScreenCostArray               <- array(NA, dim = c(length(ScreenNumbers), 4))  # Create a blank array for the discounted costs
  colnames(ScreenCostArray)     <- c("D_Ages", "D_Screen", "D_Triage", "D_EarlyTreatment")  # Name the columns in this array
  ScreenCostArray[, "D_Ages"]   <- ScreenCounts[, "ScreenAge"] - DiscountYear  # Adjust the screen age for the discount year
  ScreenCostArray[, "D_Ages"]   <- DF_Costs[ScreenCostArray[, "D_Ages"]]  # Apply the appropriate discount factors according to the screen ages 
  ScreenCostArray[, "D_Screen"] <- ScreenCostArray[, "D_Ages"] * ScreenCounts[, "N_Screens"] * ScreenCost  # Multiply the screen, triage and early treatment numbers by the discount factor
  ScreenCostArray[, "D_Triage"] <- ScreenCostArray[, "D_Ages"] * (ScreenCounts[, "TruePositives"] + ScreenCounts[, "FalsePositives"]) * TriageCost  
  ScreenCostArray[, "D_EarlyTreatment"] <- ScreenCostArray[, "D_Ages"] * ScreenCounts[, "TruePositives"] * EarlyTreatmentCost  
  
  # Totals for discounted and undiscounted values
  UndiscountedScreeningCost      <- UD_ScreenCosts + UD_TriageCosts  # Undiscounted total screening costs including triage costs
  UndiscountedTotalTreatmentCost <- UD_TreatmentCosts + ScreenedUnDiscountTreatmentCosts  # Total with screening undiscounted treatment costs including late treatment costs
  DiscountedScreeningCost        <- sum(ScreenCostArray[, c("D_Screen", "D_Triage")])  # Discounted total screening costs including triage costs
  DiscountedTotalTreatmentCost   <- sum(ScreenCostArray[, "D_EarlyTreatment"]) + ScreenedDiscountTreatmentCosts  # Total with screening discounted treatment costs including late treatment costs
  
  # Total costs under screening
  TotalUndiscountedWithScreenCosts <- UndiscountedScreeningCost + UndiscountedTotalTreatmentCost
  TotalDiscountedWithScreenCosts   <- DiscountedScreeningCost   + DiscountedTotalTreatmentCost
  
  #  Sum up the undiscounted and discounted costs
  CostEffectivenessResults <- rbind(CostEffectivenessResults, c(NA, UnDiscountedLYG, TotalUndiscountedWithScreenCosts, DiscountedLYG, TotalDiscountedWithScreenCosts, NA))
  
  # Create the vector of strategy names to be inserted later
  StrategyNames						 <- c(StrategyNames, StrategyName)

}# Close the loop over all the strategies


# Set the CostEffectivenessResults array as a data frame
CostEffectivenessResults <- as.data.frame(CostEffectivenessResults)
CostEffectivenessResults[, "StrategyName"] <- StrategyNames

# Calculate the discounted life years without screening
EffectsArray[, "D_NoScreen"] <- EffectsArray[, "UD_NoScreen"] - DiscountYear  # First adjust the discounted values for the discount year
EffectsArray[, "D_NoScreen"] <- sapply(EffectsArray[, "D_NoScreen"], PresentValue)  # Apply the discount function

IntegerLYs_NoScreen <- trunc(EffectsArray[, "D_NoScreen"],0)
DecimalLYs_NoScreen <- EffectsArray[, "D_NoScreen"]-trunc(EffectsArray[, "D_NoScreen"],0)

# Sum up the discounted LYG
DiscountedLYG_NoScreen <- sum(PresentValueFunction[IntegerLYs_NoScreen])+sum(DecimalLYs_NoScreen*DF_Effects[IntegerLYs_NoScreen+1])

# Count the costs of treating clinical disease under no screening
ClinicalCases     <- which(Outcomes[, "Clinical"] <= Outcomes[, "AllCauseDeath"]) # Identify the cases that present clinically
ClinicalOnsetAges	<- round(Outcomes[ClinicalCases, "Clinical"] - DiscountYear, 0)  # Identify the clinical ages and adjust for the discount year and round to integers
DiscountWeight		<- DF_Costs[ClinicalOnsetAges]  # Find discount weight
NoScreenUnDiscountTreatmentCosts <- sum(length(ClinicalCases) * LateTreatmentCost)
NoScreenDiscountTreatmentCosts   <- sum(DiscountWeight * LateTreatmentCost)  # Find the discounted treatment costs

CostEffectivenessResults <- rbind(CostEffectivenessResults, c("NoScreening", sum(EffectsArray[, "UD_NoScreen"]), NoScreenUnDiscountTreatmentCosts, DiscountedLYG_NoScreen, NoScreenDiscountTreatmentCosts, NA))

##################################### Calculation of ICERs ################################################################
maxWTP <- Inf
CEmat  <- cbind(Strategy = c(1:nrow(CostEffectivenessResults)), 
                Cost = as.numeric(CostEffectivenessResults[, "D_Costs"]), 
                Effectiveness = as.numeric(CostEffectivenessResults[, "D_Effects"]))

# Initialize some variables
costsCol <- 2
qalyCol  <- 3
numStrat <- nrow(CEmat)

# Find WTP levels to test so that all strategies on frontier will be captured
# This means testing on either side of all NMB intersections, which are just all the pairwise ICERs
ICERmat <- matrix(1, numStrat, numStrat)
for (i in 1:numStrat) {
  indexStrat <- matrix(1, numStrat, 3)
  indexStrat[, costsCol] <- indexStrat[, costsCol] * CEmat[i, costsCol]
  indexStrat[, qalyCol] <- indexStrat[, qalyCol] * CEmat[i, qalyCol]
  delCostQalys <- CEmat - indexStrat
  ICERmat[, i] <- delCostQalys[, costsCol] / delCostQalys[, qalyCol]
}  
intersections <- sort(unique(c(ICERmat)))
intersections <- intersections[is.finite(intersections)]
WTPtestPoints <- c(0, intersections [intersections >= 0 & intersections <= maxWTP], maxWTP)

# Find the strategy with the max NMB at each of the WTP test points
indiciesOfMax <- vector()
NMBmat <- matrix(0, numStrat, length(WTPtestPoints))
for (i in 1:length(WTPtestPoints) ) {
  NMBmat[, i] <- (WTPtestPoints[i] * CEmat[, qalyCol]) - CEmat[, costsCol]
}
if (is.infinite(maxWTP)) {   #WTP of infinity means costs are not considered
  NMBmat[, length(WTPtestPoints)] = CEmat[, qalyCol] - (0 * CEmat[, costsCol]); 
}
maxVals <- apply(NMBmat, 2, max)      # Find strategy that maximizes NMB at each WTP
for (i in 1:length(WTPtestPoints)) {  # Find all strategies that match max at each WTP
  indiciesOfMax <- c(indiciesOfMax, which( NMBmat[, i] == maxVals[i]))
}
frontier <- unique(indiciesOfMax)  # Find strategy that maximizes NMB at each WTP

BoundarySet <- CostEffectivenessResults[frontier,] # Define the boundary set as the non dominated strategies
for (i in 2:nrow(BoundarySet)) {
  BoundarySet[i, "ICER"] <- round((as.numeric(BoundarySet[i, "D_Costs"]) - as.numeric(BoundarySet[i-1, "D_Costs"])) / (as.numeric(BoundarySet[i, "D_Effects"]) - as.numeric(BoundarySet[i-1, "D_Effects"])), 3)
}

CostEffectivenessResults[frontier,"ICER"]  <- BoundarySet[,"ICER"]
CostEffectivenessResults[-frontier,"ICER"] <- "Dominated"

write.table(CostEffectivenessResults, paste("OutputFiles/LatestResults/OutputTable", Sys.Date(), ".txt", sep = ""), row.names = F, col.names = T, sep = '\t')

##################################### Simulation Ends ###########################################
runtime <- Sys.time() - time  # stop the clock
runtime





