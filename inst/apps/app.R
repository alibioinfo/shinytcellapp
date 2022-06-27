

library(rtracklayer)
library(phylobase)
#
library(bsseq)
library(dmrseq)

library(BiocManager)
options(repos = BiocManager::repositories())

load("annoTrackmm10.RData")
regions <- read.csv("targeted.regions.csv",sep = ",",skip = 0, col.names = c("symbol","chr","start","end"))

shinyApp(
  ui = fluidPage(
    
    inputPanel(
      
      
      selectizeInput("n5", label = "Select gene of interest:", choices = c(regions$symbol), selected = "Foxp3",
                     options = NULL),
      selectInput("n1", label = "Cell type 1",
                  choices =  c("DMK", "Treg" ,"tTreg" , "Th0" ,"Th1" ,"Th17" ,"Th2" , "Th17_noTGFb" , "Th17_TGFb.KO","Th17_TGFb.WT"
                  ), selected = "DMK"),
      selectInput("n2", label = "Cell type 2",
                  choices = c("DMK", "Treg" ,"tTreg" , "Th0" ,"Th1" ,"Th17" ,"Th2" , "Th17_noTGFb" , "Th17_TGFb.KO","Th17_TGFb.WT")
                  , selected = "Treg", multiple = TRUE),
      
      selectInput("n3", label = "Select modification (5mC/5hmC)",
                  choices = c("5mC","5hmC"), selected = "5mC")
      
      
    ),
    
    mainPanel(
      fluidRow(
        plotOutput("plot1", width = "850px", height = "600px"),
        downloadButton(outputId = "down", label = "Download the plot")
      )
    )
  ),
  server = function(input, output) {
    
    observeEvent(input$n3,{
      if (input$n3 == "5mC") {
        
        n1 <- reactive(input$n1)
        n2 <- reactive(input$n2)
        n3 <- reactive(input$n3)
        n5 <- reactive(input$n5)
        load("BSobj_5mCpG_integrated.RData")
        
        
        
        
        pheno <- pData(BSobj)
        
        
        pheno$group[19] <- "Th17_noTGFb"
        pheno$group[20] <- "Th17_noTGFb"
        
        
        
        library(pals)
        
        group.colors <- rainbow(10)
        group.colors <- c(rainbow(10), alphabet(10))
        library(scales)
        
        
        pheno$col <- pheno$group
        pheno$col <- gsub("DMK", group.colors[1], pheno$col)
        pheno$col <- gsub("Th0", group.colors[2], pheno$col)
        pheno$col <- gsub("Th17_noTGFb", group.colors[3], pheno$col)
        pheno$col <- gsub("Th17_TGFb.KO", group.colors[6], pheno$col)
        pheno$col <- gsub("Th17_TGFb.WT", group.colors[8], pheno$col)
        pheno$col <- gsub("Th2", group.colors[9], pheno$col)
        pheno$col <- gsub("tTreg", group.colors[15], pheno$col)
        pheno$col <- gsub("Treg", group.colors[13], pheno$col)
        pheno$col <- gsub("Th17", group.colors[19], pheno$col)
        pheno$col <- gsub("Th1", group.colors[20], pheno$col)
        
        pData(BSobj) <- pheno
        
        output$plot1 <- renderPlot({
          
          BSobj3<- BSobj[ ,BSobj$group == as.character(n1()) ]
          BSobj4<- BSobj[ ,BSobj$group == as.character(n2()) ]  
          
          
          if (length(n2())>1){
            for ( i in 1:length(n2())) { 
              
              BSobj4<- cbind(BSobj4,BSobj[ ,BSobj$group == as.character(n2()[i]) ])
              
            }
          }
          
          
          BSobj5<-cbind(BSobj3,BSobj4)
          n52 <- as.character(n5())
          
          myPlot2 <- function(){
            plotDMRs(
              BSobj5, regions[ regions$symbol == as.character(n52), ], testCovariate="group", addRegions = NULL, 
              main = paste0(as.character(n5()), " Targeted Region_5mC"),
              extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)}
          
          myPlot2()
          
          output$down <- downloadHandler(
            filename =  function() {
              paste("5mC",as.character(n52),as.character(n1()),as.character(n2()), sep=".")
            },
            
            content = function(file) {
              png(file)
              print(myPlot2())
              dev.off()
              
            } 
          )
          
        })
        
        
      }
      
      if (input$n3 == "5hmC") {
        load("BSobj_5hmCpG_integrated.RData")
        n1 <- reactive(input$n1)
        n2 <- reactive(input$n2)
        n3 <- reactive(input$n3)
        n5 <- reactive(input$n5)
        pheno <- pData(BSobj)
        pheno$group[19] <- "Th17_noTGFb"
        pheno$group[20] <- "Th17_noTGFb"
        
        
        
        
        group.colors <- rainbow(20)
        
        
        pheno$col <- pheno$group
        pheno$col <- gsub("DMK", group.colors[1], pheno$col)
        pheno$col <- gsub("Th0", group.colors[2], pheno$col)
        pheno$col <- gsub("Th17_noTGFb", group.colors[3], pheno$col)
        pheno$col <- gsub("Th17_TGFb.KO", group.colors[6], pheno$col)
        pheno$col <- gsub("Th17_TGFb.WT", group.colors[8], pheno$col)
        pheno$col <- gsub("Th2", group.colors[9], pheno$col)
        pheno$col <- gsub("tTreg", group.colors[15], pheno$col)
        pheno$col <- gsub("Treg", group.colors[13], pheno$col)
        pheno$col <- gsub("Th17", group.colors[19], pheno$col)
        pheno$col <- gsub("Th1", group.colors[20], pheno$col)
        
        pData(BSobj) <- pheno
        
        
        
        
        output$plot1 <- renderPlot({
          
          BSobj3<- BSobj[ ,BSobj$group == as.character(n1()) ]
          BSobj4<- BSobj[ ,BSobj$group == as.character(n2()) ]  
          
          
          if (length(n2())>1){
            for ( i in 1:length(n2())) { 
              
              BSobj4<- cbind(BSobj4,BSobj[ ,BSobj$group == as.character(n2()[i]) ])
              
            }
          }
          
          
          BSobj5<-cbind(BSobj3,BSobj4)
          n52 <- as.character(n5())
          
          myPlot1 <- function(){
            plotDMRs(
              BSobj5, regions[ regions$symbol == as.character(n52), ], testCovariate="group", addRegions = NULL, 
              main = paste0(as.character(n5()), " Targeted Region_5hmC"),
              extend = 0, annoTrack = annoTrack, highlightMain = F, qval = F, stat = F, horizLegend = T)}
          myPlot1()
          output$down <- downloadHandler(
            filename =  function() {
              paste("5hmC", as.character(n52),as.character(n1()),as.character(n2()), sep=".")
            },
            
            content = function(file) {
              png(file)
              print(myPlot1())
              dev.off()
              
            } 
          )
        })
        
        
      }
    })
  }
)