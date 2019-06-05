# for testing the package, Shiny
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

#install.packages("shiny")
#install.packages('rsconnect')
#install.packages('Cairo')

library(shiny)
library(rsconnect)
library(Cairo)




#Import the required data from the data sheets created by the Excel file
LifeTable			=data.frame(read.table(file="InputFiles/LifeTable.txt",    			header=T))
Incidence			=data.frame(read.table(file="InputFiles/Incidence.txt",      		header=T))
TestPerformance		=data.frame(read.table(file="InputFiles/TestPerformance.txt",		header=T))
PreclinicalDuration	=data.frame(read.table(file="InputFiles/PreclinicalDuration.txt",	header=T))
ScreenSchedule		=data.frame(read.table(file="InputFiles/ScreenSchedule.txt",		header=T))
TreatmentSuccess	=data.frame(read.table(file="InputFiles/TreatmentSuccess.txt",		header=T))
DiscountRates		=data.frame(read.table(file="InputFiles/DiscountRates.txt",			header=T))
Costs				=data.frame(read.table(file="InputFiles/Costs.txt",					header=T))
source(file="InputFiles/Misc.txt")
source(file="InputFiles/Options.txt")
#Import the results of OWSA to the ICERTable
ICERTable			=data.frame(read.table(file="OutputFiles/LatestResults/ICERTable2019-05-08.txt",    			header=T))

NumberofSegments=19
strategies <- paste(ScreenSchedule$StartAge,ScreenSchedule$StopAge,ScreenSchedule$Interval1,ScreenSchedule$TestApplied1,sep="_")
NumberOfStrategies=length(strategies)


ui<-fluidPage(
  title="Cost-effectiveness Analysis",
  fluidRow(
    column(3,
           h2("Cost-effectiveness Analysis"),  
           radioButtons("AnalysisType", label = "Analysis Type", choices = list(
             "Deterministic Analysis" = "DeterministicAnalysis",
             "One-Way Sensitivity Analysis" = "Owsa"
           )),
           hr(), 
           h3("One-Way Sensitivity Analysis"),
           radioButtons("FigureRange", label = "Set the axis range, according to...", choices = list(
             "All the Parameters" = "Large",
             "The Chosen Parameter" = "Small"
           )),
           radioButtons("parameters", label = "Parameters", choices = list(
             "Test Sensitivity" =  "TestSensitivity1" ,                 
             "Test Specificity" = "TestSpecificity1" ,                
             "Discount Rate - Costs" = "DiscountRateCosts" ,    
             "Discount Rate - Effects" = "DiscountRateEffects" ,                
             "Discount Year" = "DiscountYear" ,            
             "PreClinical Probability" = "PreClinicalProbability" , 
             "Clinical Probability" = "ClinicalProbability" ,     
             "Cost of Primary Screen" = "CostofPrimaryScreen" ,          
             "Cost of FollowUp "= "CostofFollowUp" ,                
             "Cost of Treatment (Detected by Screens)"= "CostofTreatScreenDetected" ,
             "Cost of Treatment Clinical"= "CostofTreatmentClinical",
             "Incidence Rate" = "Incidence",
             "Preclinical Duration" = "PreclinicalDuration")
            )#,           
          # hr(), 
          
    ),
    column(8,
           plotOutput("distPlot", click = "plot_click"),
           hr()
    ),
    column(3,offset = 1,
           sliderInput("slider1", label = "Parameter Value", min = 1, max = (NumberofSegments+1), value=0, step = 1),
           br(),
           sliderInput("slider2", label = "Threshold", min = 10000, max = 200000, value=100000, step = 5000)
    ),
    column(3,offset = 2,
           h4("The Strategy under Oberservation"),
           sliderInput("startage", label = "Screen Startage", min = 20, max = 40, value=30, step = 5),
           sliderInput("stopage", label = "Screen Stopage", min = 50, max = 70, value=70, step = 5),
           sliderInput("interval", label = "Screen Interval", min = 1, max = 10, value=1, step = 1)
           #hr(),
    )
  )
)


server<-shinyServer(function(input, output) {
  output$distPlot <- renderPlot({
  
   #**************************************************ONE-WAY SENSITIVITY ANALYSES (OWSA)  END************************************************************************************************

      #Define the scenario we are going to plot on the CE plane:as.numeric(input$slider1)
      #Define the threshold
      if (input$AnalysisType == "DeterministicAnalysis") {Scenario<- 1} else{
      Scenario<- which(rownames(ICERTable)==paste(input$parameters,input$slider1,sep="_"))}
      Threshold<-input$slider2
      ObsearvationPoint<- paste(input$startage,"_",input$stopage,"_",input$interval,"_",1,sep="")
    
      
      PlotData<-array(0,dim=c(NumberOfStrategies,3))
      rownames(PlotData)<- strategies
      colnames(PlotData)<- c("D_Costs","D_Effects","ICER")
      
      #Define the costs and effects; costs here expressed in millions
      PlotData[,"D_Costs"]<-as.numeric(ICERTable[Scenario,][colnames(ICERTable) %in% paste("D_Costs",strategies,sep="_")])/10^6
      PlotData[,"D_Effects"]<-as.numeric(ICERTable[Scenario,][colnames(ICERTable) %in% paste("D_Effects",strategies,sep="_")])
      PlotData[,"ICER"]<-as.character(unlist(ICERTable[Scenario,][colnames(ICERTable) %in% paste("ICER",strategies,sep="_")]))
   
      if (input$AnalysisType == "DeterministicAnalysis") {
        MaxCost=round(max(as.numeric(ICERTable[Scenario,][colnames(ICERTable) %in% paste("D_Costs",strategies,sep="_")])/10^6),-3)+500
        MaxEffect=t=round(max(as.numeric(ICERTable[Scenario,][colnames(ICERTable) %in% paste("D_Effects",strategies,sep="_")])),-3)+1000
      } else if (input$FigureRange == "Small"){
      SameParameterD_Costs<-ICERTable[c(rownames(ICERTable) %in% paste((strsplit(paste(rownames(ICERTable)),split="_")[[Scenario]][1]),1:(NumberofSegments+1),sep="_")),][colnames(ICERTable) %in% paste("D_Costs",strategies,sep="_")]
      SameParameterD_Effects<-ICERTable[c(rownames(ICERTable) %in% paste((strsplit(paste(rownames(ICERTable)),split="_")[[Scenario]][1]),1:(NumberofSegments+1),sep="_")),][colnames(ICERTable) %in% paste("D_Effects",strategies,sep="_")]
      
      if (round(max(SameParameterD_Costs),-3)-max(SameParameterD_Costs)>=0)
      {MaxCost=round(max(SameParameterD_Costs)/10^6,-3)} else {MaxCost=round(max(SameParameterD_Costs)/10^6,-3)+500}
      
      if (round(max(SameParameterD_Effects),-3)-max(SameParameterD_Effects)>=0)
      {MaxEffect=round(max(SameParameterD_Effects),-3)} else {MaxEffect=round(max(SameParameterD_Effects),-3)+1000}
   
      }else {
        MaxCost=round(max(ICERTable[,colnames(ICERTable) %in% paste("D_Costs",strategies,sep="_")])/10^6,-3)+500
        MaxEffect=round(max(ICERTable[,colnames(ICERTable) %in% paste("D_Effects",strategies,sep="_")]),-3)+1000
      }    
      
      #*****************************************  To Plot the CE plane*********************************************************************************************************
      #Plot the CE plane
      Effects				=PlotData[,"D_Effects"]
      Costs				=PlotData[,"D_Costs"]
      
      
      #BoundarySet is the strategies whose ICER is not "SD" or "ED"
      BoundarySet<- PlotData[which(!(PlotData[,"ICER"]%in%c("SD","ED"))),]
      if (length(which(!(PlotData[,"ICER"]%in%c("SD","ED"))))==1){}else{
        
        #Do the same for the efficient frontier
        {
          EfficientEffects	=BoundarySet[,"D_Effects"]
          EfficientCosts		=BoundarySet[,"D_Costs"]}
        
        #OptimalStrategy:
        #Identify the highest of the cost-effective ICERs given the threshold
        BoundarySetWithoutRef<-BoundarySet[which(BoundarySet[,"ICER"]!=c("reference")),]
        
        if (min(as.numeric(BoundarySetWithoutRef[,"ICER"]))<=Threshold){
          OptimalStrategy<-rownames(BoundarySetWithoutRef)[which(BoundarySetWithoutRef[,"ICER"]==(max(as.numeric(BoundarySetWithoutRef[,"ICER"])[as.numeric(BoundarySetWithoutRef[,"ICER"])<=Threshold],na.rm=TRUE)))]
          CostEffectiveICER=BoundarySetWithoutRef[OptimalStrategy,"ICER"]
          #Find the corresponding costs and effect
          CostEffectiveEffects=BoundarySetWithoutRef[OptimalStrategy,"D_Effects"]
          CostEffectiveCosts	=BoundarySetWithoutRef[OptimalStrategy,"D_Costs"]
        }else{
          CostEffectiveICER=c("reference")
          CostEffectiveEffects=BoundarySet[,"D_Effects"][which(BoundarySet[,"ICER"]=="reference")]
          CostEffectiveCosts=BoundarySet[,"D_Costs"][which(BoundarySet[,"ICER"]=="reference")]
          
        }
        
        
        
        
        
        
        #Set the tick marks to incude zero and any big number(here are 10,000 and 2,000, separately)
        TicksEffects		=pretty(c(0,MaxEffect	))
        TicksCosts			=pretty(c(0,MaxCost	))
        
        
        
        #Set the ranges for both values
        xRange=range(TicksEffects)
        yRange=range(TicksCosts	 )
        
        #Plot the plane
        #Set the margin parameters
        par(mar=c(5, 6, 4, 1))
        plot(Effects,Costs,yaxt="n",xaxt="n", xlab="",ylab="", xlim=xRange,ylim=yRange, cex=1)
        #Add the axes
        axis(1, at = TicksEffects,	labels = formatC(TicksEffects,  big.mark = ",", format = "d"),las = 1)
        axis(2, at = TicksCosts	 , 	labels = formatC(TicksCosts,	big.mark = ",", format = "d"),las = 1)
        #Add in the X and Y axis labels with spacing
        title(paste("Scenario",strsplit(paste(rownames(ICERTable)[Scenario]),split="_")[[1]][2],":",strsplit(paste(rownames(ICERTable)[Scenario]),split="_")[[1]][1],"=",sprintf("%.3f",as.numeric(ICERTable[Scenario,paste(strsplit(paste(rownames(ICERTable)[Scenario]),split="_")[[1]][1])]),sep=" ")),ylab="Costs (?M)", xlab="Effects, LYG", mgp=c(3.75,1.75,0), cex.lab=1.2)

        
        #Add the efficient frontier
        sequence=order(as.numeric(EfficientEffects))
        lines (EfficientEffects[sequence],EfficientCosts[sequence],lwd=2)
        points(EfficientEffects,EfficientCosts,pch=19,lwd=2)
        #Add the optimal strategy
        points(CostEffectiveEffects,CostEffectiveCosts,pch=19,lwd=3,cex=1,col="forestgreen")
        #Add a legend
        #legend("topright",legend=c(paste("The Optimal Strategy",paste(ifelse(OptimalStrategy %in% "NA","NA",OptimalStrategy)))),pch=19,bty="o",pt.cex =1,col="forestgreen", cex=0.8)
        #Add a line to show the threshold
        abline(a=as.numeric(CostEffectiveCosts)-as.numeric(CostEffectiveEffects)*Threshold/10^6,b=Threshold/10^6,         col="red", lwd=3, lty=2)
        #Add the specific strategy that we want to observe
        points(PlotData[ObsearvationPoint,"D_Effects"],PlotData[ObsearvationPoint,"D_Costs"],pch=19,lwd=3,cex=1,col="red")
        #Add a legend
        #legend("topright",legend=c(paste("The Obsearved Strategy",paste(ifelse(ObsearvationPoint %in% "NA","NA",ObsearvationPoint)))),pch=19,bty="o",pt.cex =1,col="red", cex=0.8)
        
        #Create ICER labels: only when the strategy's "D_Effects" is higher than or equal to 4,000
        ICERLabels=gsub("NA","",noquote(format(c(BoundarySet[which(BoundarySet[,"D_Effects"]>=4000),"ICER"],CostEffectiveICER),big.mark=",")))
        
        #Add the ICERs to the frontier
        #Generate a label offset to place the values at an offset
        LabelYOffset=max(as.numeric(EfficientCosts))/40
        #Apply the text
        text(c(EfficientEffects[which(BoundarySet[,"D_Effects"]>=4000)],CostEffectiveEffects),c(as.numeric(EfficientCosts)[which(BoundarySet[,"D_Effects"]>=4000)],as.numeric(CostEffectiveCosts))-LabelYOffset,labels=prettyNum(ICERLabels,big.mark=",",scientific=FALSE), pos=4)
        
        
       
      }
    
  })
})
  
    
  

shinyApp(ui=ui,server=server)


