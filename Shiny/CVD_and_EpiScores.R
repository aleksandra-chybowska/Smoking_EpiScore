################################################################################

### Shiny app - CVD EpiScores (GenScot)

################################################################################

# Developed by Aleksandra Chybowska (https://github.com/Aleksandra-Chybowska), 
# Based on proteins, diseases and UKB application by Danni A Gadd (https://github.com/DanniGadd)
# Emails - a.d.chybowska@sms.ed.ac.uk, danni.gadd@ed.ac.uk

# Load packages
library(shiny)
library(readxl)
library(thematic) 
library(ggplot2) 
library(gridExtra)
library(patchwork)
library(forcats)
library(tidytext)
library(dplyr)
library(survminer)
library(survival)
library(tidyr)

##################################################

# Join in the correct naming for the plots as generated above 
full_annots = readRDS("/Volumes/marioni-lab/Ola/Lab/EpiScores/Annotations/EpiAnnots_prepped.RDS")

path = '/Volumes/marioni-lab/Ola/Lab/Cox/basic_dataset_no_assumptions/'

f1 = 'HRs_assign_epi.csv'
f2 = 'HRs_assign_epi_cTnI_corrected.csv'

mod = read.csv(paste0(path, f1))
mod = left_join(mod, full_annots, by=c("protein" = "ID")) # 109
mod$Model = "ASSIGN"
mod$Name = mod$Gene
mod = subset(mod, !duplicated(Name)) #103
mod = subset(mod, local > 0.05 & global > 0.05) #101
number_of_protein_episcores = nrow(mod)
threshold = 0.05/number_of_protein_episcores
bonf_corrected = subset(mod, p<threshold) #36
mod = subset(mod, p<0.05) # 67

troponin = read.csv(paste0(path, f2))
troponin = left_join(troponin, full_annots, by=c("protein" = "ID")) # 109
troponin$Model = "ASSIGN + cTnI"
troponin$Name = troponin$Gene
troponin = subset(troponin, !duplicated(Name)) #103
troponin = subset(troponin, local > 0.05 & global > 0.05) #101
number_of_protein_episcores = nrow(troponin)
threshold = 0.05/number_of_protein_episcores
bonf_corrected = subset(troponin, p<threshold) #33
troponin = subset(troponin, p<0.05) # 65

increase = subset(troponin, hr > 1)
decrease = subset(troponin, hr < 1)

x = rbind(mod, troponin)
episcore_names = x$Name
x = x %>% mutate(Name = reorder(Name, hr))
x$Outcome <- "Protein Episcores, p<0.05"
x = mutate(x, Name = reorder(Name, hr))

##################################################

### FORMAT INPUT DATA: PANEL ONE - COX PH 16 year 

# Load data
bind <- read.csv('/Volumes/marioni-lab/Ola/Lab/Cox/basic_dataset_no_assumptions/episcores_over_time_iterations_censor.csv')
protein_names <- sort(unique(bind$Prot_name))

#Set variable to dictate colour based on failures of local PH
bind$col_variable <- ifelse(bind$Shoenfeld == 'Local P < 0.05', "#D73526", "#4E84C4")

# Set variable to dictate shape based on whether P full < 0.05
bind$shape_variable <- ifelse(bind$Association == 'P < 0.05', 'circle', 'triangle')

# Set annnotations
Text1 <- paste(" Failures at the local protein level are shown in red.", "\n", "Circular points correspond to P values < 0.05, whereas triangles show no association.")

# Relabel for cases/controls cols
names(bind)[c(8,9)] <- c('Cases','Controls')

##################################################

### PANEL TWO - Survival plot

path = '/Volumes/marioni-lab/Ola/Lab/Cox/basic_dataset_no_assumptions/'
tb = read.csv(paste0(path, "episcores_and_assign_input_data.csv"), row.names = 1)
epi_names = colnames(tb[6:ncol(tb)])

##############################################################################################
### SHINY APP CODE 
### DEFINE UI

ui <- fluidPage(
  
  tags$head(tags$style(HTML(
    "
      .selectize-input, .selectize-dropdown {
        font-size: 16px;
      },
      h1 {
        font-size: 30px;
        margin-bottom: 30px;
      }
    
    "))),   
  
  theme = bslib::bs_theme(bootswatch = "superhero"), # sandstone
  
  headerPanel("Protein EpiScores and CVD risk"),
  
  tabsetPanel(               
    
    tabPanel("FOREST PLOT",
      sidebarLayout(
        sidebarPanel(width = 4, 
                      selectInput("p", "EpiScore", choices = episcore_names, multiple = TRUE),
                     
        ),
        mainPanel(
          fluidRow(
            column(12, 
              plotOutput("plot1")
            )
          )
        )
      )
    ),
    
    tabPanel("RISK OVER TIME",
      sidebarLayout(
        sidebarPanel(width = 4, selectInput("epi", "EpiScore", choices = protein_names, selected = "CRP")
        ),
        
        mainPanel(
          fluidRow(
            column(12, 
              plotOutput("plot2"),
              textOutput("dfStr")
            )
          )
        )
      )
    ),
    tabPanel("SURVIVAL PROBABILITY",
          sidebarLayout(
            sidebarPanel(width = 4, selectInput("scores", "EpiScore", choices = epi_names),
            numericInput("min", "Lower Quantile:", min = 1, max = 99, value = 25),
            numericInput("max", "Upper Quantile:", min = 1, max = 99, value = 75)
          ),
          mainPanel(
            fluidRow(
              column(12, 
                     plotOutput("plot3")

              )
            )
          )
        )
      )
    )
  )


### DEFINE SERVER LOGIC

server <- function(input, output, session) {
  
  ### TAB ONE - forest plot
  subset <- reactive({
    if (length(input$p) == 0) {
      return(x)
    }
    table = dplyr::filter(x, Name %in% input$p)
    return(table)
  })

  
  output$plot1 <- renderPlot({
    stacked = ggplot(subset(),aes(y=hr, x=Name, group=Model, colour=Model)) + 
      geom_point(size = 2, position = position_dodge(0.5))+
      geom_errorbar(aes(ymin = lci, ymax = uci),
                    position = position_dodge(0.5), width = 0.1)+
      ylab("Hazard Ratio per SD [95% Confidence Interval]")+ 
      xlab ("") +
      theme_bw() +
      geom_hline(yintercept = 1, linetype = "dotted")+
      theme(axis.text.x = element_text(size = 8, vjust = 0.5), 
            axis.text.y = element_text(size = 8), legend.position = "bottom",
            plot.title = element_text(size = 8))+ theme(legend.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = c("#4E84C4", "#D73526")) +
      coord_flip() 
    stacked
    
  }, res = 96, width="auto", height = 900)
  
  
  #####################################################################
  
  ### TAB TWO - year of follow up 
  
  # Subset to chosen protein and disease 
  filtered_data <- reactive({
    test <- dplyr::filter(bind, Prot_name == input$epi)
    return(test)
  })
  
  # Plot 
  output$plot2 <- renderPlot({
    prot_name  <- as.character(unique(filtered_data()$Prot_name))
    
    p1 <- ggplot(filtered_data(), aes(Iteration, Hazard.Ratio))+
      geom_point(aes(shape = Association), fill = filtered_data()$col_variable, colour = filtered_data()$col_variable, size = 5,
                 shape = filtered_data()$shape_variable)+ 
      geom_errorbar(aes(ymin = LCI, ymax = UCI), colour = filtered_data()$col_variable,
                    position = position_dodge(0.5), width = 0, linewidth = 0.75) +
      ylab("Hazard Ratio [95% CI]")+ xlab ("")+ theme_classic() +
      geom_hline(yintercept = 1, linetype = "dotted")+
      theme(title = element_text(colour = "aquamarine3"), axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "bottom",
            axis.text.y = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
      theme_classic() + ggtitle(paste0(prot_name, ' - CVD')) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    plot_data <- filtered_data() %>% 
      select('Iteration', 'Cases', 'Controls') %>%
      mutate(Year = Iteration,
             Cases = Cases,
             Controls = Controls) %>% 
      pivot_longer(c(Year, Cases, Controls), names_to = "layer", values_to = "label") 
    
    plot_data$layer <- factor(plot_data$layer, levels = c("Controls", "Cases", "Year"))
    
    p2 <-  ggplot(plot_data, aes(x = Iteration)) +
      geom_text(aes(y = factor(layer, c("Controls", "Cases", "Year")), label = label)) +
      labs(y = "", x = NULL) +
      theme_minimal() +
      theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
            panel.grid = element_blank(), strip.text = element_blank())
    
    p1 / p2  +  plot_layout(heights = c(4, 1))
    
  }, res = 96, width="auto")
  
  # Print coption for annotations
  output$dfStr <- renderText({
    print(Text1)
  })
  
  ###########################################################
  
  ### TAB THREE - survival plot
  
  findEpi <- reactive({
    selected = tb %>% select("age", "sex", "assign", "tte", "event")
    colname = input$scores
    idx = which(colnames(tb)==input$scores)
    col = tb[, idx]
    min = ifelse(input$min/100 <= 0, 1, input$min/100)
    max = ifelse(input$max/100 >= 100, 99, input$max/100)
    top = quantile(col, probs = c(min,max))[1]
    bottom = quantile(col, probs = c(min,max))[2]
    selected$protein = ifelse(tb[idx] >= top, 2,
                               ifelse(tb[idx] < bottom, 0, NA))
    selected = na.omit(selected)
    return(selected)
  })
  
  modelPrep <- reactive({
    mod1 <- survfit(Surv(tte, event) ~ protein, data=findEpi())
    return(mod1)
  })
  
  output$plot3 <- renderPlot({
    
    prot_name  <- as.character(unique(filtered_data()$Prot_name))
    
    plot = ggsurvplot(modelPrep(),
               data = findEpi(),
               conf.int=TRUE, # add confidence intervals
               pval=TRUE, # show the p-value for the log-rank test
               risk.table=FALSE, # show a risk table below the plot
               legend.labs=c(paste("<", "Lower"), 
                             paste(">", "Upper")), # change group labels
               legend.title=paste("Quantile"),  # add legend title
               palette=c("#4E84C4", "#D73526"), # change colors of the groups
               ylab="CVD free survival",
               xlab="Time (years)",
               risk.table.height=.2)
    
    plot
    
  }, res = 96, width="auto", height = 450)
  
}

# Run the application 
shinyApp(ui = ui, server = server)


