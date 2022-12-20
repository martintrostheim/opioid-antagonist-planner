library(shiny)
library(rhandsontable)
library(linpk)
library(ggplot2)
library(plotly)
library(pracma) # Contains Lambert W function needed to find ka
#library(shinythemes)
# https://benjaminrich.github.io/linpk/vignettes/linpk-intro.html

ui <- fluidPage(
  
  sidebarLayout(
    sidebarPanel(
      #shinythemes::themeSelector(),
      h4("Participant info"),
      numericInput("kg","Weight (kg)",value=70, step=1, min=0, max=NA),
      h4("Session info"),
      numericInput("Session","Session length (h)",value=1, step=1, min=0, max=NA),
      numericInput("Washout","Inter-session interval (h)",value=24, step=1, min=0, max=NA),
      textOutput("text_washout"),
      textOutput("text_recommended"),
      h4("Study info"),
      numericInput("n_participants", "Number of participants", value=1, min=1, max=NA),
      numericInput("n_sessions", "Number of naloxone sessions", value=1, min=1, max=NA),
      numericInput("price", "Price per 0.4 mg/ml vial or ampoule", value=NA, min=0, max=NA),
      htmlOutput("text_naloxone_amount"),
      htmlOutput("text_cost"),
      h4("Plot settings"),
      sliderInput("resolution", "Curve resolution (dots per minute)", value=10, step=1, min=1, max=100),
      checkboxInput("showDOR", "Simulate DOR blockade", value=FALSE),
      checkboxInput("showKOR", "Simulate KOR blockade", value=FALSE),
      checkboxInput("truncate", "Limit blockade to 0-100%", value=TRUE),
      hr(),
      htmlOutput("text_citation"),
      hr(),
      textOutput("text_credits")
    ),
    
    mainPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Main",
                           plotlyOutput("Plot"),
                           hr(),
                           downloadButton("download_curve", "Download blockade-time profile"),
                           hr(),
                           
                           fluidRow(
                             column(6,
                                    h4("Drug administration"),
                                    rHandsontableOutput("table_drug"),
                                    htmlOutput("text_drug"),
                                    downloadButton("download_dose", "Download drug administration info")
                             ),
                             column(6,
                                    h4("Outcome assessment"),
                                    rHandsontableOutput("table_measure"),
                                    htmlOutput("text_measure"),
                                    downloadButton("download_outcome", "Download outcome assessment info")
                             )
                           )),
                  tabPanel("Model parameters",
                           h4("Model parameters"),
                           numericInput("MOR.ED50", "ED50 MOR blockade (mg/kg)", value=0.002290531, min=0, max=NA),
                           numericInput("T.max", "Blockade tmax (minutes)", value=25, min=0, max=NA),
                           numericInput("Halflife", "Blockade half-life (minutes)", value=110, min=0, max=NA),
                           numericInput("DOR.ratio", "MOR:DOR affinity ratio (1:X)", value=41, min=0, max=NA),
                           numericInput("KOR.ratio", "MOR:KOR affinity ratio (1:X)", value=8, min=0, max=NA))
                  )
    )
  )
)  
 
server <- function(input, output, session) {
  options(scipen = 999)
  
  # Model from paper
  T.measurement <- 55
  #T.max <- 25
  #Halflife <- 110
  
  # MOR
  #MOR.ED50 <- 0.002290758
  
  # DOR
  #DOR.ratio <- 41
  DOR.ED50 <- reactive({input$MOR.ED50*input$DOR.ratio})
  
  # KOR
  #KOR.ratio <- 8
  KOR.ED50 <- reactive({input$MOR.ED50*input$KOR.ratio})
  

  k <- reactive({log(2)/input$Halflife})
  ka <- reactive({-1/input$T.max * lambertWn(-input$T.max*k()*exp(-input$T.max*k()))})
  vc <- 1
  cl <- reactive({k()*vc})
  
  # Blockade curves
  MOR.Blockade.0 <- function(Dose, MOR.ED50, k, T.measurement){
    #output <- (MOR.loglogit.c +(MOR.loglogit.d-MOR.loglogit.c)/(1+(Dose/MOR.loglogit.e)^MOR.loglogit.b)) * exp(-k*(0-T.measurement))
    output <- ((Dose*100)/(Dose+MOR.ED50)) * exp(-k*(0-T.measurement))
    return(output)
  }
  
  DOR.Blockade.0 <- function(Dose, DOR.ED50, k, T.measurement){
    #output <- (DOR.loglogit.c +(DOR.loglogit.d-DOR.loglogit.c)/(1+(Dose/DOR.loglogit.e)^DOR.loglogit.b)) * exp(-k*(0-T.measurement))
    output <- ((Dose*100)/(Dose+DOR.ED50)) * exp(-k*(0-T.measurement))
    return(output)
  }
  
  KOR.Blockade.0 <- function(Dose, KOR.ED50, k, T.measurement){
    #output <- (KOR.loglogit.c +(KOR.loglogit.d-KOR.loglogit.c)/(1+(Dose/KOR.loglogit.e)^KOR.loglogit.b)) * exp(-k*(0-T.measurement))
    output <- ((Dose*100)/(Dose+KOR.ED50)) * exp(-k*(0-T.measurement))
    return(output)
  }
  
  # Drug administration
  df_drug <- data.frame(Time = c(0, rep(NA_integer_, 14)),
                   Dose = c(0.4, rep(NA_integer_, 14)),
                   Unit = factor(c("mg", rep(NA_character_, 14)), levels=c("mg", "mg/kg", "mg/min", "mg/kg/min", "mg/h", "mg/kg/h",
                                                                           "mcg", "mcg/kg", "mcg/min", "mcg/kg/min", "mcg/h", "mcg/kg/h")), 
                   Duration = c(0, rep(NA_integer_, 14)),
                   stringsAsFactors = FALSE)
  
  output$table_drug <- renderRHandsontable({rhandsontable(df_drug,
                                                          stretchH='all',
                                                          colHeaders=c('<p title="Time point of IV naloxone administration (minutes)">Time<br>(minutes)</p>',
                                                                       '<p title="IV naloxone dose">Dose</p>',
                                                                       '<p title="IV naloxone dose unit (mg, mg/kg, mg/min, mg/kg/min, mg/h, mg/kg/h, mcg, mcg/kg, mcg/min, mcg/kg/min, mcg/h or mcg/kg/h)">Dose unit</p>',
                                                                       '<p title="Infusion duration (minutes). Set to 0 for bolus">Duration<br>(minutes)</p>')) %>%
      hot_col(2, format="0.00000")})
  
  info_drug <- reactive({
    DF <- hot_to_r(input$table_drug)
    DF <- na.omit(DF)
    })
  
  info_drug_Unit <- reactive({info_drug()$Unit})
  info_drug_Dose <- reactive({info_drug()$Dose})
  info_drug_Time <- reactive({info_drug()$Time})
  info_drug_Duration <- reactive({info_drug()$Duration})
  
  # Convert all doses to mg/kg
  info_drug_Dose_adj <- reactive({
    DF <- rep(NA, length(info_drug_Dose()))
    DF[info_drug_Unit()=="mg"] <- info_drug_Dose()[info_drug_Unit()=="mg"] / input$kg
    DF[info_drug_Unit()=="mg/kg"] <- info_drug_Dose()[info_drug_Unit()=="mg/kg"]
    DF[info_drug_Unit()=="mg/min"] <- info_drug_Dose()[info_drug_Unit()=="mg/min"] / input$kg * info_drug_Duration()[info_drug_Unit()=="mg/min"]
    DF[info_drug_Unit()=="mg/kg/min"] <- info_drug_Dose()[info_drug_Unit()=="mg/kg/min"] * info_drug_Duration()[info_drug_Unit()=="mg/kg/min"]
    DF[info_drug_Unit()=="mg/h"] <- info_drug_Dose()[info_drug_Unit()=="mg/h"] / input$kg / 60 * info_drug_Duration()[info_drug_Unit()=="mg/h"]
    DF[info_drug_Unit()=="mg/kg/h"] <- info_drug_Dose()[info_drug_Unit()=="mg/kg/h"] / 60 * info_drug_Duration()[info_drug_Unit()=="mg/kg/h"]
    
    DF[info_drug_Unit()=="mcg"] <- info_drug_Dose()[info_drug_Unit()=="mcg"] / input$kg / 1000
    DF[info_drug_Unit()=="mcg/kg"] <- info_drug_Dose()[info_drug_Unit()=="mcg/kg"] / 1000
    DF[info_drug_Unit()=="mcg/min"] <- info_drug_Dose()[info_drug_Unit()=="mcg/min"] / input$kg / 1000 * info_drug_Duration()[info_drug_Unit()=="mcg/min"]
    DF[info_drug_Unit()=="mcg/kg/min"] <- info_drug_Dose()[info_drug_Unit()=="mcg/kg/min"] / 1000 * info_drug_Duration()[info_drug_Unit()=="mcg/kg/min"]
    DF[info_drug_Unit()=="mcg/h"] <- info_drug_Dose()[info_drug_Unit()=="mcg/h"] / input$kg / 60 / 1000 * info_drug_Duration()[info_drug_Unit()=="mcg/h"]
    DF[info_drug_Unit()=="mcg/kg/h"] <- info_drug_Dose()[info_drug_Unit()=="mcg/kg/h"] / 60 / 1000 * info_drug_Duration()[info_drug_Unit()=="mcg/kg/h"]
    DF
  })
  
  # Convert dose to blockade at administration time point, assuming no absorption
  info_drug_Blockade_MOR <- reactive({MOR.Blockade.0(Dose=info_drug_Dose_adj(), MOR.ED50=input$MOR.ED50, k=k(), T.measurement=T.measurement)})
  info_drug_Blockade_DOR <- reactive({DOR.Blockade.0(Dose=info_drug_Dose_adj(), DOR.ED50=DOR.ED50(), k=k(), T.measurement=T.measurement)})
  info_drug_Blockade_KOR <- reactive({KOR.Blockade.0(Dose=info_drug_Dose_adj(), KOR.ED50=KOR.ED50(), k=k(), T.measurement=T.measurement)})
  
  
  output$text_drug <- renderText({
    paste0("Time = IV naloxone administration time point (minutes)<br>",
           "Dose = IV naloxone dose<br>",
           "Unit = IV naloxone dose unit<br>",
           "Duration = Infusion duration (minutes)"
    )
  })
  
  input_check <- reactive({
    check <- sum((isTruthy(input$kg) & input$kg > 0),
                 (isTruthy(input$Session) & input$Session > 0), 
                 isTruthy(info_drug_Dose()), 
                 isTruthy(info_drug_Time()),
                 isTruthy(info_drug_Duration()),
                 isTruthy(info_drug_Unit()),
                 isTruthy(input$MOR.ED50),
                 isTruthy(input$T.max),
                 isTruthy(input$Halflife),
                 isTruthy(input$DOR.ratio),
                 isTruthy(input$KOR.ratio))
  })  

    x <- reactive({
      if(input_check() == 11){
        seq(0 ,input$Session*60, input$resolution^-1)
        }
      })
    
    y.MOR <- reactive({
      if(input_check() == 11){
        y <- as.data.frame(pkprofile(x(), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_MOR(), dur=info_drug_Duration())))$conc
        if(input$truncate == TRUE){
          y <- ifelse(y > 100, 100, y)
          }
        }
      y
      })
    
    y.DOR <- reactive({
      if(input_check() == 11){
        y <- as.data.frame(pkprofile(x(), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_DOR(), dur=info_drug_Duration())))$conc
        if(input$truncate == TRUE){
          y <- ifelse(y > 100, 100, y)
        }
      }
      y
    })
    
    y.KOR <- reactive({
      if(input_check() == 11){
        y <- as.data.frame(pkprofile(x(), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_KOR(), dur=info_drug_Duration())))$conc
        if(input$truncate == TRUE){
          y <- ifelse(y > 100, 100, y)
        }
      }
      y
    })
    
    data <- reactive({
      if(input_check() == 11){
        data.frame(minutes=x(), MOR=y.MOR(), DOR=y.DOR(), KOR=y.KOR())
      } else{
        data.frame(minutes=NA, MOR=NA, DOR=NA, KOR=NA)
      }
    })
    
    ymax <- reactive({
      if(input$truncate == TRUE){
        100
        } else{
          if(!is.na(max(data()$MOR)) & max(data()$MOR) > 150){
            max(data()$MOR)
            } else{
              150
            }
          }
      })
    
output$Plot <- renderPlotly({
#output$Plot <- renderPlot({
    
  p <- ggplot()+
    labs(title="Opioid Receptor Blockade with IV Naloxone", y="Blockade (%)", x="Time (minutes after administration)")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), text=element_text(size=12))+
    ylim(0, ymax())
  
  if(isTruthy(input$Session) & input$Session > 0){
    p <- p+
      xlim(0, input$Session*60)
  }
  
  if(input_check() == 11){
    if(nrow(rect_measure()) > 0){
      p <- p+
        geom_rect(aes(xmin=rect_measure()$Start, xmax=rect_measure()$End, ymin=rect_measure()$ymin, ymax=rect_measure()$ymax), fill=alpha("blue", 0.1), color=alpha("blue", 0.3))+
        geom_text(aes(x=rect_measure()$Start+(rect_measure()$End-rect_measure()$Start)/2, y=rect_measure()$ymax/4, label=rect_measure()$Measure), color=alpha("blue", 0.3), size=10)
    }
    p <- p+
      geom_line(data=data(), aes(x=minutes, y=MOR, color="MOR"), lwd=1)+
      scale_color_manual(name="Receptor", values=c("MOR" = "black", "KOR" = "blue", "DOR" = "red")) # https://aosmith.rbind.io/2018/07/19/manual-legends-ggplot2/
    if(input$showDOR == TRUE){
      p <- p+
        geom_line(data=data(), aes(x=minutes, y=DOR, color="DOR"), lwd=1)
    }
    if(input$showKOR == TRUE){
      p <- p+
        geom_line(data=data(), aes(x=minutes, y=KOR, color="KOR"), lwd=1)
    }
  }
   ggplotly(p) %>%
     config(displayModeBar = FALSE) %>%
     layout(xaxis = list(fixedrange = TRUE), yaxis = list(fixedrange = TRUE), legend = list(itemclick = FALSE, itemdoubleclick = FALSE))
})
  
  # Outcome assessment
  df_measure <- data.frame(Measure = rep(NA_character_, 15),
                           Start = rep(NA_integer_, 15),
                           End = rep(NA_integer_, 15),
                           MOR = rep(NA_integer_, 15),
                           DOR = rep(NA_integer_, 15),
                           KOR = rep(NA_integer_, 15),
                           stringsAsFactors = FALSE)
  
  df_measure_previous <- reactive({df_measure})
  
  df_measure_changes <- reactive({
    if(is.null(input$table_measure)){return(df_measure_previous())}
    else if(!identical(df_measure_previous(),input$table_measure)){
      # hot.to.df function will convert your updated table into the dataframe
      DF <- hot_to_r(input$table_measure)
      # here the second column is a function of the first and it will be multipled by 100 given the values in the first column
      for(i in 1:nrow(DF)){
        if(!is.na(DF$Start[i]) & !is.na(DF$End[i]) & DF$Start[i] <= DF$End[i] & input_check() == 11){
          MOR <- as.data.frame(pkprofile(seq(DF$Start[i], DF$End[i], input$resolution^-1), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_MOR(), dur=info_drug_Duration())))$conc
          DOR <- as.data.frame(pkprofile(seq(DF$Start[i], DF$End[i], input$resolution^-1), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_DOR(), dur=info_drug_Duration())))$conc
          KOR <- as.data.frame(pkprofile(seq(DF$Start[i], DF$End[i], input$resolution^-1), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_KOR(), dur=info_drug_Duration())))$conc
          if(input$truncate == TRUE){
            MOR <- ifelse(MOR > 100, 100, MOR)
            DOR <- ifelse(DOR > 100, 100, DOR)
            KOR <- ifelse(KOR > 100, 100, KOR)
          }
          DF$MOR[i] <- mean(MOR)
          DF$DOR[i] <- mean(DOR)
          DF$KOR[i] <- mean(KOR)
          } else{
            DF$MOR[i] <- NA_integer_
            DF$DOR[i] <- NA_integer_
            DF$KOR[i] <- NA_integer_
          }
        }
      DF
    }
  })
  
  output$table_measure <- renderRHandsontable({rhandsontable(df_measure_changes(),
                                                             stretchH='all',
                                                             colHeaders=c('<p title="Name of outcome measure">Measure</p>',
                                                                          '<p title="Start time of outcome assessment (minutes)">Start</p>',
                                                                          '<p title="End time of outcome assessment (minutes)">End</p>',
                                                                          '<p title="Average MOR blockade (%)">MOR<br>blockade</p>',
                                                                          '<p title="Average DOR blockade (%)">DOR<br>blockade</p>',
                                                                          '<p title="Average KOR blockade (%)">KOR<br>blockade</p>'))})
  
  output$text_measure <- renderText({
    paste0("Measure = Measurement name<br>",
           "Start = Assessment start time (minutes)<br>",
           "End = Assessment end time (minutes)<br>",
           "MOR = Average MOR blockade (%)<br>",
           "DOR = Average DOR blockade (%)<br>",
           "KOR = Average KOR blockade (%)")
  })
  
  rect_measure <- reactive({
    DF <- df_measure_changes()
    DF$ymin <- 0
    DF$ymax <- ymax()
    DF <- DF[!is.na(DF$MOR),]
    DF$Measure[is.na(DF$Measure)] <- " "
    DF
  })
  
  info_measure <- reactive({
    DF <- df_measure_changes()
    DF <- DF[!is.na(DF$MOR),]
    DF$Measure[is.na(DF$Measure)] <- " "
    DF
  })
  
  # Download
  output$download_curve <- downloadHandler(
    filename <- function() {
      paste0("planoxone_data-curve_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv")
    },
    content <- function(file) {
      write.csv(data(), file, row.names=FALSE)
    })

  output$download_dose <- downloadHandler(
    filename <- function() {
      paste0("planoxone_data-dose_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv")
    },
    content <- function(file) {
      write.csv(info_drug(), file, row.names=FALSE)
    })

  output$download_outcome <- downloadHandler(
    filename <- function() {
      paste0("planoxone_data-outcome_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".csv")
    },
    content <- function(file) {
      write.csv(info_measure(), file, row.names=FALSE)
    })
  
  # Inter-session interval
  output$text_washout <- renderText({
    text <- "Residual MOR blockade at next session: "
    if(input_check() == 11 & isTruthy(input$Washout)){
      Blockade.wash <- mean(as.data.frame(pkprofile(c(input$Session*60+input$Washout*60), cl=cl(), vc=vc, ka=ka(), dose=list(t.dose=info_drug_Time(), amt=info_drug_Blockade_MOR(), dur=info_drug_Duration())))$conc)
      if(input$truncate == TRUE & Blockade.wash > 100){
        Blockade.wash <- 100
      }
      text <- paste0(text, round(Blockade.wash), "%")
    }
    text
  })
  
  output$text_recommended <- renderText({
    Washout.recommended <- c(5, 10)*input$Halflife/60
    paste0("Recommended minimum intersession interval: ", round(min(Washout.recommended),1), "-", round(max(Washout.recommended),1), " hours")
  })
  
  output$text_naloxone_amount <- renderText({
    text_amount <- "Naloxone dose per participant per session: "
    text_amount_total <- "Total amount of naloxone needed for this study: "
    if(input_check()==11 & isTruthy(input$n_participants) & isTruthy(input$n_sessions)){
      amount <- sum(info_drug_Dose_adj())
      amount_total <- amount*input$n_participants*input$n_sessions*input$kg
      text_amount <- paste0(text_amount, round(amount,5), " mg/kg (", round(amount*input$kg,5), " mg)")
      text_amount_total <- paste0(text_amount_total, round(amount_total,5), " mg")
    }
    paste0(text_amount, "<br>", text_amount_total)
  })
  
  output$text_cost <- renderText({
    text_cost_total <- "Total cost for this study: "
    if(input_check()==11 & isTruthy(input$n_participants) & isTruthy(input$n_sessions) & isTruthy(input$price)){
      price <- input$price/0.4
      cost <- sum(info_drug_Dose_adj()) * price * input$kg
      cost_total <- cost*input$n_participants*input$n_sessions
      text_cost_total <- paste0(text_cost_total, round(cost_total,5))
    }
    paste0(text_cost_total)
  })
  
  # Citation
  output$text_citation <- renderText({
    paste0("<b>Article (Open Access)</b><br>",
           #'Tr\u00f8stheim, M., Eikemo, M., Haaker, J., Frost, J. J., and Leknes, S. (2022). Opioid Antagonism in Humans: A Primer on Optimal Dose and Timing for Central Mu-Opioid Receptor Blockade. <i>bioRxiv</i>, 2022.02.25.481943. <a href="https://doi.org/10.1101/2022.02.25.481943">https://doi.org/10.1101/2022.02.25.481943</a><br>',
           '<a href="https://doi.org/10.1038/s41386-022-01416-z">Tr\u00f8stheim et al. (2022, <i>Neuropsychopharmacology</i>)</a></br>',
           "<br><b>Source code</b><br>",
           '<a href="https://github.com/martintrostheim/opioid-antagonist-planner">GitHub</a><br>',
           "<br><b>Links</b><br>",
           '<a href="https://sirileknes.com/">Leknes Affective Brain Lab</a><br>',
           '<a href="https://martintrostheim.shinyapps.io/plantrexone/">Plantrexone</a>')
  })
  
  # Credits and affiliations
  output$text_credits <- renderText({
    "Developed by Martin Tr\u00f8stheim"
  })
  
}

shinyApp(ui, server)