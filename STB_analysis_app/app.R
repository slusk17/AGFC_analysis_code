#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=
# The purpose of this application is to analyze data 
# collected during standardized gill netting for various 
# Moronid species
#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=

library(shiny)
library(ggplot2) # ggplot()
library(FSA) # lencat(), psdVal(), psdci()
library(dplyr) # group_by(), summarize(), mutate(), tally(), left_join()
library(lubridate) # year()
library(extrafont) # element_text(family = "Times New Roman")

rm(list=ls(all=TRUE)) # Clear environment

# Define UI for application that draws a histogram
ui = fluidPage(
    
    titlePanel("Striped Bass Analysis"),
    
    h5(HTML("The purpose of this application is to analyze data collected 
            during standardized gill netting for various Moronid species<br/><br/>")),
    h5(strong("Methods")),
    h5(HTML("Gillnetting was used to sample the *Moronid species* population following the Preferred 
            Method described in AGFC Striped Bass Species Management Plan 2nd Edition. Sampling was conducted 
            on *Dates* with *number of nets* experimental gillnets at *number of different locations* locations 
            for a total effort of *number of net nights* net nights (nn). All *Moronid species* were measured 
            (TL, nearest mm) and weighed (g). <br/><br/> Sampling statistics [Catch per unit effort (CPUE)] 
            and structural indices [Proportional Size Distribution-Quality (PSD-Q) and Proportional Size 
            Distribution-Memorable (PSD-M)] were calculated for *Moronid species* according to Neumann et al. 
            (2012). Relative standard error (RSE) was calculated for CPUE while standard error (SE) was 
            calculated for PSD-M.<br/><br/>")),
    h5(strong("Results")),
    h5(HTML("A total of *number of target Morone species individuals collected* *Morone species* 
            were collected via standardized sample for a mean catch per effort of *catch per net 
            night* fish/nn [*relative standard error*; *sample size* (RSE; n)]. Collected *Morone species* 
            ranged in size from *minimum size* to *maximum size* with an estimated PSD-Q of 
            *estimated PSD-Q* [*standard error* (SE)] and PSD-M of *estimated PSD-M* [*standard error* (SE)] 
            (*Figure X).<br/><br/>")),
    
    hr(),
    
    fluidRow(
        column(3, style = "background-color: rgba(214, 214, 214, 0.6)",
               fileInput("csvFile", "Select .CSV")),
        
        column(9,
               plotOutput("lengthfreq")
        )))

# Define server logic required to draw a histogram
server <- function(input, output) {

    dat <- eventReactive(input$csvFile, {
        read.csv(input$csvFile$datapath, na.strings = ".")
    })
    
    output$lengthfreq <- renderPlot({
        
        dat <- dat() #This sets the uploaded data "dat()" to object dat
        
        dat$date <- as.Date(dat$date, format = "%d/%m/%Y") # Convert "date" column to an actual date format
        
        num_nets <- dat %>% # Calculate the number of nets/night
            group_by(date) %>% # Group by date
            summarize(nets = n_distinct(net)) %>% # Count the number of nets/night
            mutate(year = year(date)) %>% # create a new column with the year
            group_by(year) %>% # Group by new year column
            summarize(num_nets = sum(nets)) # Sum the total number of nets used each year
        
        dat$lcat25 <- lencat(dat$length, w = 25) # Create a new column "lcat25" that assigns all fish to a 25 mm length group
        
        # Calculate the catch per effort and associated RSE#
        
        cpue <- dat %>% # Calculate the catch per effort by length group for plotting
            group_by(date, net, lcat25) %>% # Group_by year, net and lcat25
            tally() %>% # Tally the number of fish 
            mutate(year = year(date)) %>% # create a new column with the year
            left_join(num_nets) %>% # Join the "num_nets" dataframe with the "dat" dataframe 
            mutate(cpue = n/num_nets) # Create a cpue column
        
        estcpue <- dat %>% # Calculate the catch per effort by length group for plotting
            mutate(year = year(date)) %>% # create a new column with the year
            group_by(date, year, net) %>% # Group_by year, net and lcat25
            tally() %>% # Tally the number of fish per net
            group_by(year) %>% # Group by year
            summarize(meancpue = mean(n), # Calculate the mean cpue
                      SEcpue = sd(n, na.rm = T)/sqrt(sum(!is.na(n))), # Calculate the SE of CPUE
                      rsecpue = (SEcpue / meancpue) * 100) # Calculate the RSE of CPUE
        
        # Calculate PSD-Q/PSD-M and associated SE #
        
        psd <- dat %>%
            mutate(year = year(date)) %>% # create a new column with the year
            filter(length >= psdVal("Striped Bass (landlocked)")["stock"]) %>% # select only STB greater than or equal to "stock" size
            mutate(gcat = lencat(length, breaks = psdVal("Striped Bass (landlocked)"), # Create a column "gcat" and assign all individuals to a manegement category
                                 use.name = TRUE, drop.levels = TRUE)) %>%
            group_by(date, net, year) %>% # Group by date, net and year
            summarize(s = (prop.table(xtabs(~gcat))*100)["stock"], # Calculate the proportion of "stock" fish
                      q = (prop.table(xtabs(~gcat))*100)["quality"], # Calculate the proportion of "quality" fish
                      p = (prop.table(xtabs(~gcat))*100)["preferred"], # Calculate the proportion of "preferred" fish
                      m = (prop.table(xtabs(~gcat))*100)["memorable"], # Calculate the proportion of "memorable" fish
                      t = (prop.table(xtabs(~gcat))*100)["trophy"]) %>% # Calculate the proportion of "trophy" fish
            mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% # Convert "NA" values to 0 
            group_by(date, net, year) %>% # Group_by date, net and year
            summarize(psdq = (sum(q, p, m, t)/ sum(s, q, p, m, t)) * 100, # Calculate psdq
                      psdm = (sum(m, t)/ sum(s, q, p, m, t)) * 100) %>% # Caluclate psdm
            group_by(year) %>% # group by year
            summarize(meanpsdq = mean(psdq), # calculate mean psdq
                      SEpsdq = sd(psdq, na.rm = T)/sqrt(sum(!is.na(psdq))), # Calculate SE of psdq
                      meanpsdm = mean(psdm), # Calculate mean psdm
                      SEpsdm = sd(psdm, na.rm = T)/sqrt(sum(!is.na(psdm)))) # Calculate the SE of psdm
        
        #=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=
        # Create standard time lapse size structure figure
        #=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=#=~=
        
        windows()
        
        #tiff("CPUE_figure.jpeg", width = 15, height = 15, units = "cm", res = 300)
        
        ggplot(data = cpue, aes(x = lcat25, y = cpue, group = year)) + # Create basic plot of number of fish per length group 
            stat_summary(fun.y = sum, geom = "bar", fill = "grey", colour = "black") + # Create basic plot of number of fish per length group
            facet_wrap(~year, ncol = 1) + # Facet plots by year
            
            theme_classic() + # Select theme for plot
            theme (axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0), colour = "black"), #set the size, spacing, and color for the y-axis and x-axis titles
                   axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0), colour = "black"),
                   text = element_text(family = "Times New Roman"), #set the font type
                   plot.title = element_text(face = "bold", family = "Arial"), #modify plot title, the B in this case
                   legend.position = c(0.3,0.85), #position the legend on the figure
                   legend.text = element_text(size = 12), #adjust size of text for legend
                   plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"), #margin for the plot
                   axis.ticks.y = element_line(size = 0.5), #set size of the tick marks for y-axis
                   axis.ticks.x = element_line(size = 0.5), #set size of the tick marks for x-axis
                   axis.ticks.length = unit(0.2,"cm"), #adjust length of the tick marks
                   axis.text.y = element_text(colour = "black", size = 14, angle = 0, vjust = 0.5, hjust = 1, #set size and location of the tick labels for the y axis
                                              margin = margin(t = 0, r = 5, b = 0, l = 0)),
                   axis.text.x = element_text(colour = "black", size = 14, angle = 0, vjust = 0, hjust = 0.5, #set size and location of the tick labels for the x axis
                                              margin = margin(t = 5, r = 0, b = 0, l = 0)),
                   axis.line = element_line(colour = "black", size = 0.5, lineend = "square"), #set the axis size, color, and end shape
                   strip.text.y = element_blank(),
                   strip.text.x = element_blank()) + # Facet label size
            
            annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) + # Add x-axis tick line to all figures 
            
            geom_text(data = data.frame(label = as.factor(unique(cpue$year)), # Add year "labels" to each plot
                                        year = unique(cpue$year),
                                        x     = max(dat$lcat25, na.rm = T) - 50,
                                        y     = max(cpue$cpue, na.rm = T)),
                      mapping = aes(x = x, y = y, label = label), # add labels
                      family = "Times New Roman", # Font family
                      fontface = "bold", 
                      size = 4) + # Size of labels
            
            geom_text(data.frame(label = paste0("CPUE", " ", "=", " ",as.character(round(estcpue$meancpue)),"(",as.character(round(estcpue$rsecpue)),")"), # Add PSD text
                                 year  = unique(psd$year),
                                 x     = max(dat$lcat25, na.rm = T) - 50,
                                 y     = max(cpue$cpue, na.rm = T) * 0.75),
                      mapping = aes(x = x, y = y, label = label),
                      family = "Times New Roman",
                      fontface = "bold",
                      size = 3) +
            
            
            geom_text(data.frame(label = paste0("PSD-Q", " ", "=", " ",as.character(round(psd$meanpsdq)),"(",as.character(round(psd$SEpsdq)),")"), # Add PSD text
                                 year  = unique(psd$year),
                                 x     = max(dat$lcat25, na.rm = T) - 50,
                                 y     = max(cpue$cpue, na.rm = T) * 0.55),
                      mapping = aes(x = x, y = y, label = label),
                      family = "Times New Roman",
                      fontface = "bold",
                      size = 3) +
            
            geom_text(data.frame(label = paste0("PSD-M", " ", "=", " ",as.character(round(psd$meanpsdm)),"(",as.character(round(psd$SEpsdm)),")"), # Add PSD text
                                 year  = unique(psd$year),
                                 x     = max(dat$lcat25, na.rm = T) - 50,
                                 y     = max(cpue$cpue, na.rm = T) * 0.35),
                      mapping = aes(x = x, y = y, label = label),
                      family = "Times New Roman",
                      fontface = "bold",
                      size = 3) +
            
            scale_x_continuous(breaks = seq(min(dat$lcat25, na.rm = T)- 25, max(dat$lcat25, na.rm = T) + 25, 100), lim = c(min(dat$lcat25, na.rm = T) - 25, max(dat$lcat25, na.rm = T) + 25)) + # Set x axis tick intervals and labels 
            
            xlab("Length (mm)") + # Add x axis label
            ylab("Mean Striped Bass/net night") # Add y axis label
        
    })
}

shinyApp(ui = ui, server = server)
