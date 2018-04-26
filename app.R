#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a heatmap
ui <- fluidPage(
  titlePanel("Pregnancy associated genes of a live bearing cockroach"),
  sidebarLayout(
    sidebarPanel(
      img(src = "cockroach.png",height=120.3,width=254.65 ),
      selectInput("comp", label = "Groups compared:",
                  choices = c("Not Pregnant vs First Trimester"="NPvFT", "Not Pregnant vs Second Trimester"="NPvST", "Not Pregnant vs Third Trimester"="NPvTT", "Not Pregnant vs Male"="NPvM", "First Trimester vs Second Trimester"="FTvST","First Trimester vs Third Trimester"="FTvTT","First Trimester vs Male"="FTvM","SecondTrimester vs Third Trimester"="STvTT","Second Trimester vs Male"="STvM","Third Trimester vs Male"="TTvM"), selected = "Not Pregnant vs First Trimester"),
      
      sliderInput("pval", label = "Maximum P-value",
                  min = 0, max = 1, value = .05, step = 0.001),
      sliderInput("minfc",label="Minimum Log2 Fold Change from reference:", min=-1000,max=1000,value=1,step=5),
      sliderInput("maxfc",label="Maximum Log2 Fold Change from reference:", min=-1000,max=1000,value=10,step=5)
      
    ),
    # Show a plot of the generated distribution
    mainPanel(
      
      textOutput("rowcnt1"),
     plotlyOutput("selecteddata")
      
      #dataTableOutput("table"),
      #tableOutput("colnames")
      
    )
    
  )
)




# Define server logic required to draw a histogram
server <- function(input, output, session) {
  raw<-countframe
  output$selecteddata=renderPlotly({
    if(input$comp=="NPvFT") {
      filt=quote(NPvFT.padj<input$pval & NPvFT.FC<=input$maxfc & NPvFT.FC>=input$minfc)
      
    }
    else if (input$comp=="NPvST"){
      filt=quote(NPvST.padj<=input$pval & NPvST.FC<=input$maxfc & NPvST.FC>=input$minfc)
    }
    else if (input$comp=="NPvTT"){
      filt=quote(NPvTT.padj<=input$pval & NPvTT.FC<=input$maxfc & NPvTT.FC>=input$minfc)
    }
    else if (input$comp=="NPvM"){
      filt=quote(NPvM.padj<=input$pval & NPvM.FC<=input$maxfc & NPvM.FC>=input$minfc)
    }
    else if (input$comp=="FTvST"){
      filt=quote(FTvST.padj<=input$pval & FTvST.FC<=input$maxfc & FTvST.FC>=input$minfc)
    }
    else if (input$comp=="FTvTT"){
      filt=quote(FTvTT.padj<=input$pval & FTvTT.FC<=input$maxfc & FTvTT.FC>=input$minfc)
    }
    else if (input$comp=="FTvM"){
      filt=quote(FTvM.padj<=input$pval & FTvM.FC<=input$maxfc & FTvM.FC>=input$minfc)
    }
    else if (input$comp=="STvTT"){
      filt=quote(STvTT.padj<=input$pval & STvTT.FC<=input$maxfc & STvTT.FC>=input$minfc)
    }
    else if (input$comp=="STvM"){
      filt=quote(STvM.padj<=input$pval & STvM.FC<=input$maxfc & STvM.FC>=input$minfc)
    }
    else {
      filt=quote(TTvM.padj<=input$pval & TTvM.FC<=input$maxfc & TTvM.FC>=input$minfc)
    }
    
    dropit=function(x)(x[,2:16])
    dropit1=function(x)(x[,1])
    collabs<-c("Not.Pregnant1","Not.Pregnant2","Not.Pregnant3",
                                        "First.Trimester1","First.Trimester2","First.Trimester3",
                                        "Second.Trimester1","Second.Trimester2","Second.Trimester3",
                                        "Third.Trimester1","Third.Trimester2","Third.Trimester3",
                                        "Male1","Male2","Male3")
  
    m <- list(
      l = 110,
      t=150
    )  
  names1=as.list(dropit1(raw %>% filter_(filt)))
   t=as.matrix(dropit(raw %>% filter_(filt)))
   rownames(t)=names1
   colnames(t)=collabs
   plot_ly(x=colnames(t),y=rownames(t),z=t,type="heatmap")%>%
     layout(margin = m,xaxis=list(side="top",title="Pregnancy Stage"),autosize=F,width=700,height=3000)
      
    
  })
  
  
  rowcnt=reactive({
    if(input$comp=="NPvFT") {
      filt=quote(NPvFT.padj<input$pval & NPvFT.FC<=input$maxfc & NPvFT.FC>=input$minfc)
      
    }
    else if (input$comp=="NPvST"){
      filt=quote(NPvST.padj<=input$pval & NPvST.FC<=input$maxfc & NPvST.FC>=input$minfc)
    }
    else if (input$comp=="NPvTT"){
      filt=quote(NPvTT.padj<=input$pval & NPvTT.FC<=input$maxfc & NPvTT.FC>=input$minfc)
    }
    else if (input$comp=="NPvM"){
      filt=quote(NPvM.padj<=input$pval & NPvM.FC<=input$maxfc & NPvM.FC>=input$minfc)
    }
    else if (input$comp=="FTvST"){
      filt=quote(FTvST.padj<=input$pval & FTvST.FC<=input$maxfc & FTvST.FC>=input$minfc)
    }
    else if (input$comp=="FTvTT"){
      filt=quote(FTvTT.padj<=input$pval & FTvTT.FC<=input$maxfc & FTvTT.FC>=input$minfc)
    }
    else if (input$comp=="FTvM"){
      filt=quote(FTvM.padj<=input$pval & FTvM.FC<=input$maxfc & FTvM.FC>=input$minfc)
    }
    else if (input$comp=="STvTT"){
      filt=quote(STvTT.padj<=input$pval & STvTT.FC<=input$maxfc & STvTT.FC>=input$minfc)
    }
    else if (input$comp=="STvM"){
      filt=quote(STvM.padj<=input$pval & STvM.FC<=input$maxfc & STvM.FC>=input$minfc)
    }
    else {
      filt=quote(TTvM.padj<=input$pval & TTvM.FC<=input$maxfc & TTvM.FC>=input$minfc)
    }
    
    dropit=function(x)(x[,2:16])
    addnames=function(x)(colnames(x)<-c("Not.Pregnant1","Not.Pregnant2","Not.Pregnant3",
                                        "First.Trimester1","First.Trimester2","First.Trimester3",
                                        "Second.Trimester1","Second.Trimester2","Second.Trimester3",
                                        "Third.Trimester1","Third.Trimester2","Third.Trimester3",
                                        "Male1","Male2","Male3"))
    nrow(dropit(raw %>%
                  filter_(filt)))
    
    
  })
  
  contignames=reactive({
    if(input$comp=="NPvFT") {
      filt=quote(NPvFT.padj<input$pval & NPvFT.FC<=input$maxfc & NPvFT.FC>=input$minfc)
      
    }
    else if (input$comp=="NPvST"){
      filt=quote(NPvST.padj<=input$pval & NPvST.FC<=input$maxfc & NPvST.FC>=input$minfc)
    }
    else if (input$comp=="NPvTT"){
      filt=quote(NPvTT.padj<=input$pval & NPvTT.FC<=input$maxfc & NPvTT.FC>=input$minfc)
    }
    else if (input$comp=="NPvM"){
      filt=quote(NPvM.padj<=input$pval & NPvM.FC<=input$maxfc & NPvM.FC>=input$minfc)
    }
    else if (input$comp=="FTvST"){
      filt=quote(FTvST.padj<=input$pval & FTvST.FC<=input$maxfc & FTvST.FC>=input$minfc)
    }
    else if (input$comp=="FTvTT"){
      filt=quote(FTvTT.padj<=input$pval & FTvTT.FC<=input$maxfc & FTvTT.FC>=input$minfc)
    }
    else if (input$comp=="FTvM"){
      filt=quote(FTvM.padj<=input$pval & FTvM.FC<=input$maxfc & FTvM.FC>=input$minfc)
    }
    else if (input$comp=="STvTT"){
      filt=quote(STvTT.padj<=input$pval & STvTT.FC<=input$maxfc & STvTT.FC>=input$minfc)
    }
    else if (input$comp=="STvM"){
      filt=quote(STvM.padj<=input$pval & STvM.FC<=input$maxfc & STvM.FC>=input$minfc)
    }
    else {
      filt=quote(TTvM.padj<=input$pval & TTvM.FC<=input$maxfc & TTvM.FC>=input$minfc)
    }
    
    dropit1=function(x)(x[,1])
    addnames=function(x)(colnames(x)<-c("Not.Pregnant1","Not.Pregnant2","Not.Pregnant3",
                                        "First.Trimester1","First.Trimester2","First.Trimester3",
                                        "Second.Trimester1","Second.Trimester2","Second.Trimester3",
                                        "Third.Trimester1","Third.Trimester2","Third.Trimester3",
                                        "Male1","Male2","Male3"))
    dropit1(raw %>%
                  filter_(filt))
    
    
  })
  
  output$contigs=renderTable(contignames())
  
  output$rowcnt1=renderText(paste(rowcnt(), "genes match your criteria"))


  table1=renderTable(selecteddata(),rownames=TRUE)
  
 
output$heatmap=renderPlot({heatmaply(selecteddata)})
  
}
# Run the application 
shinyApp(ui = ui, server = server)

