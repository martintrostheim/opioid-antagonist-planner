library(shiny)
library(ggplot2)
library(plotly)

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      h4("Dose info"),
      numericInput("Dose", "PO naltrexone dose (mg)", value=50, step=1, min=0, max=NA),
      h4("Study info"),
      numericInput("n_participants", "Number of participants", value=1, min=1, max=NA),
      numericInput("n_sessions", "Number of naltrexone sessions", value=1, min=1, max=NA),
      numericInput("price", "Price per 50 mg tablet", value=NA, min=0, max=NA),
      htmlOutput("text_naltrexone_amount"),
      htmlOutput("text_cost"),
      h4("Recommendations"),
      htmlOutput("text_recommendations"),
      hr(),
      h4("Plot settings"),
      sliderInput("resolution", "Curve resolution (dots per minute)", value=10, step=1, min=1, max=100),
      checkboxInput("showDOR", "Show DOR blockade", value=FALSE),
      checkboxInput("showKOR", "Show KOR blockade", value=FALSE),
      htmlOutput("text_warning"),
      hr(),
      htmlOutput("text_citation"),
      hr(),
      textOutput("text_credits")
    ),
    
    mainPanel(
      plotOutput("Plot"),
      hr(),
      downloadButton("download_curve", "Download dose-blockade curve"),
      hr(),
      plotOutput("Plot_bar")
    )
  )
  
)

server <- function(input, output, session) {
  options(scipen = 999)
  
  output$text_warning <- renderText({
    note <- 'Note: '
    text <- 'DOR and KOR blockade are approximated from MOR blockade based on the affinity of naloxone to these receptors relative to MOR.'
    if(input$showDOR == TRUE | input$showKOR == TRUE){
      text <- paste0('<font color=\"#FF0000\">Warning! ', text, '</font>')
    } else{
      text <- paste0(note, text)
    }
    text
  })
  
  
  # Total naltrexone amount needed
  output$text_naltrexone_amount <- renderText({
    text_mg <- "Total mg needed: "
    text_tablets <- "50 mg tablets needed: "
    if(isTruthy(input$Dose) & isTruthy(input$n_participants) & isTruthy(input$n_sessions)){
      mg <- (input$Dose*input$n_participants*input$n_sessions)
      tablets <- ceiling(mg/50)
      text_mg <- paste0(text_mg, mg)
      text_tablets <- paste0(text_tablets, tablets)
    }
    paste0(text_mg, "<br>", text_tablets)
  })
  
  # Cost
  output$text_cost <- renderText({
    text_cost_total <- "Total cost: "
    if(isTruthy(input$Dose) & isTruthy(input$n_participants) & isTruthy(input$n_sessions) & isTruthy(input$price)){
      mg <- (input$Dose*input$n_participants*input$n_sessions)
      tablets <- ceiling(mg/50)
      cost_total <- tablets*input$price
      text_cost_total <- paste0(text_cost_total, cost_total)
    }
    paste0(text_cost_total)
  })
  
  # Recommendations
  output$text_recommendations <- renderText({
    paste0('<b>Delay</b><br>',
           'Recommended minimum delay between administration and outcome assessment: 1-2 hours<br>',
           '<br><b>Intersession interval</b><br>',
           'Recommended minimum intersession interval to avoid residual blockade: 15 days')
  })
  
  # Dose-blockade curves
  MOR.ED50 <- 5.59 # ED50 from Rabiner et al. (2011)
  DOR.ED50 <- MOR.ED50*79
  KOR.ED50 <- MOR.ED50*2
  
  f.MOR <- function(dose){return(dose/(dose+MOR.ED50)*100)}
  f.DOR <- function(dose){return(dose/(dose+DOR.ED50)*100)}
  f.KOR <- function(dose){return(dose/(dose+KOR.ED50)*100)}
  
  dose <- reactive({10^seq(-2, 3, input$resolution^-1)})
  
  data <- reactive({
    Receptors <- data.frame(Dose=dose(), Receptor=rep("MOR", length(dose())), Blockade=f.MOR(dose=dose()))
    if(input$showDOR == TRUE){
      DOR <- data.frame(Dose=dose(), Receptor=rep("DOR", length(dose())), Blockade=f.DOR(dose=dose()))
      Receptors <- rbind(Receptors, DOR)
    }
    if(input$showKOR == TRUE){
      KOR <- data.frame(Dose=dose(), Receptor=rep("KOR", length(dose())), Blockade=f.KOR(dose=dose()))
      Receptors <- rbind(Receptors, KOR)
    }
    Receptors
  })
  
  output$Plot <- renderPlot({
    p <- ggplot()+
      #geom_hline(yintercept=90, linetype="dashed", color = "black")+
      geom_line(data=data(), aes(x=log10(Dose), y=Blockade, color=Receptor))+
      ylim(0,100)+
      scale_x_continuous(breaks=-2:3, labels=10^(-2:3), sec.axis=dup_axis(labels=-2:3, name="log10(Dose [mg])"))+
      labs(x="Dose (mg)", y="Blockade (%)", title="Opioid Receptor Blockade with PO Naltrexone (t < 8 hours)")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))
    
    if(input$showDOR == TRUE & input$showKOR == FALSE){
      p <- p+
        scale_color_manual(name="Receptor", values=c("MOR"="black", "DOR"="red"), breaks=c("MOR", "DOR"))
    } else if(input$showDOR == FALSE & input$showKOR == TRUE){
      p <- p+
        scale_color_manual(name="Receptor", values=c("MOR"="black", "KOR"="blue"), breaks=c("MOR", "KOR"))
    } else if(input$showDOR == TRUE & input$showKOR == TRUE){
      p <- p+
        scale_color_manual(name="Receptor", values=c("MOR"="black", "DOR"="red", "KOR"="blue"), breaks=c("MOR", "DOR", "KOR"))
    } else if(input$showDOR == FALSE & input$showKOR == FALSE){
      p <- p+
        scale_color_manual(name="Receptor", values=c("MOR"="black"), breaks=c("MOR"))
    }
    p
  })
  
  # Download dose-blockade curve
  output$download_curve <- downloadHandler(
    filename <- function() {
      paste0("data-curve_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv")
    },
    content <- function(file) {
      write.csv(data(), file, row.names=FALSE)
    })
  
  # Blockade estimates
  res <- reactive({
    if(!isTruthy(input$Dose)){
      Dose <- 0
    } else(
      Dose <- input$Dose
    )
    
    Receptors <- data.frame(Dose=Dose, Receptor="MOR", Blockade=round(f.MOR(dose=Dose)))
    if(input$showDOR == TRUE){
      DOR <- data.frame(Dose=Dose, Receptor="DOR", Blockade=round(f.DOR(dose=Dose)))
      Receptors <- rbind(Receptors, DOR)
    }
    if(input$showKOR == TRUE){
      KOR <- data.frame(Dose=Dose, Receptor="KOR", Blockade=round(f.KOR(dose=Dose)))
      Receptors <- rbind(Receptors, KOR)
    }
    Receptors
  })
  
  output$Plot_bar <- renderPlot({
    title <- ifelse(isTruthy(input$Dose),
                    paste0("Opioid Receptor Blockade with ", input$Dose, " mg PO Naltrexone (t < 8 hours)"),
                    "Opioid Receptor Blockade with PO Naltrexone (t < 8 hours)")
    
    p <- ggplot(data=res(), aes(x=factor(Receptor, level=c("MOR", "DOR", "KOR")), y=Blockade, fill=Receptor))+
      labs(x="Receptor", y="Blockade (%)", title=title)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))+
      #scale_fill_manual(name="Receptor", values=c("MOR"="black", "DOR"="red", "KOR"="blue"), breaks=c("MOR", "DOR", "KOR"))+
      scale_y_continuous(breaks=seq(0,100,25), limits=c(0,110))
    
    if(input$showDOR + input$showKOR == 0){
      p <- p+
        geom_bar(stat="identity", width=0.5/3)
    } else if(input$showDOR + input$showKOR == 1){
      p <- p+
        geom_bar(stat="identity", width=0.5/2)
    } else if(input$showDOR + input$showKOR == 2){
      p <- p+
        geom_bar(stat="identity", width=0.5)
    }
    
    if(isTruthy(input$Dose)){
      p <- p+
        geom_text(aes(label=paste0(Blockade, "%")), vjust=-1)
    }

    if(input$showDOR == TRUE & input$showKOR == FALSE){
      p <- p+
        scale_fill_manual(name="Receptor", values=c("MOR"="black", "DOR"="red"), breaks=c("MOR", "DOR"))
    } else if(input$showDOR == FALSE & input$showKOR == TRUE){
      p <- p+
        scale_fill_manual(name="Receptor", values=c("MOR"="black", "KOR"="blue"), breaks=c("MOR", "KOR"))
    } else if(input$showDOR == TRUE & input$showKOR == TRUE){
      p <- p+
        scale_fill_manual(name="Receptor", values=c("MOR"="black", "DOR"="red", "KOR"="blue"), breaks=c("MOR", "DOR", "KOR"))
    } else if(input$showDOR == FALSE & input$showKOR == FALSE){
      p <- p+
        scale_fill_manual(name="Receptor", values=c("MOR"="black"), breaks=c("MOR"))
    }
    
    p
  })
  
  # Citation
  output$text_citation <- renderText({
    paste0("<b>Citation</b><br>",
           "[Reference to preprint/paper]<br>",
           "<br><b>Source</b><br>",
           '<a href="https://github.com/martintrostheim/opioid-antagonist-planner">GitHub</a><br>',
           "<br><b>Links</b><br>",
           '<a href="https://sirileknes.com/">Leknes Affective Brain Lab</a>')
  })
  
  # Credits and affiliations
  output$text_credits <- renderText({
    "Developed by Martin Tr\u00f8stheim"
  })
  
}

shinyApp(ui, server)