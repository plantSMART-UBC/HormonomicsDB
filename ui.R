#hormonomicsDB v1.4.2
#July 21st 2021
#Authors: Ryland T. Giebelhaus, Lauren A.E. Erland, and Susan J. Murch
#PlantSMART research group at UBC Okanagan
#contact: Dr. Susan J. Murch. Email: susan.murch@ubc.ca
#Website: hormonomicsDB.com

#call shiny library in
library(shiny)
library(DT)

#####

ui <- fluidPage(

  titlePanel("HormonomicsDB"),

  sidebarLayout(
    sidebarPanel(
      strong("HormonomicsDB v1.4.2"),
      br(),
      p("Developed by: Ryland T. Giebelhaus, Lauren A.E. Erland, Susan J. Murch"),
      p("The University of British Columbia"),
      p("Contact: susan.murch@ubc.ca"),
      p("HormonomicsDB is a tool developed at UBC Okanagan by the plantSMART
        research team to allow users to process their untargeted metabolomics
        data to putativley identify plant hormones."),
      textOutput("timesRUN"), ##prints the output
    ),
    mainPanel(
      tabsetPanel(type = "tabs",


                  tabPanel("Instructions",
                           br(),

                           strong("Instructions"),
                           #edit this text to change the instructions
                           p("Use either the 'm/z screener' which searches against our hormonomics datasets or select the
                            'Custom Database Search' to upload your own dataset use our platform to perform your
                            own custom queries of your untargeted metabolomics data. View your output results in the tab
                            next to the tool you used then download your results as a .csv file."),
                           br(),

                           strong("Database Descriptions: "),
                           p(("PGR Monoisotopic and M+H: Only the monoisotopic mass and M+H adduct for the plant growth regulators in
                          ESI+ mode.")),
                           p(("PGR Adducts: Common adducts of plant growth regulators in ESI+ mode.")),
                           p(("PGR Biotransformations: Common predicted biotransformations of plant growth regulators in ESI+ mode.")),
                           p(("PGR Adduct and Biotransformations: Both adducts and predicted biotransformations for plant growth regulators
                          in ESI+ mode.")),
                           br(),

                           strong("Code Availability"),
                           p("All source code and files are available at https://github.com/plantSMART-UBC/HormonomicsDB"),
                           br(),
                           strong("Terms and Agreements"),
                           p("HormonomicsDB was developed for research use only and is not intended for use in diagnostic work.
                             Despite diligent validation and bug fixing, we are not responsible for any mistakes the application
                             makes in data processing. Considering this, please inform us immediately of any bugs that you encounter."),
                           p("We do not save any data that is uploaded to the server, it is
                             immediately deleted with every new session that you start."),
                           p("Please cite the github repository for this project when using HormonomicsDB in
                             any work."),
                           tags$hr(),),

                  tabPanel("M/Z Screener",

                           br(),

                           strong("Instructions: "),
                           p("Select which datasets to search from then select a search tolerance
                             and then upload your formatted data as a .csv and allow up to 3 minutes to perform
                             the search. After this is completed select how you want your data ordered and
                             view it in the 'Screener Output' tab."),
                           p("Note, 'ionization mode' only works for M+H and M-H adducts, as well as synthetic biotransformations
                             at the moment."),
                           
                           #ESI + or -
                           radioButtons("ionMode", "Select Ionization Mode: ",
                                        choices = list("ESI +" = 1,
                                                       "ESI -" = 2),
                                        selected = 1),

                           #gives checkboxes so user can select multiple datasets to search from at once
                           checkboxGroupInput("dataset", "Choose Dataset: ",
                                              choices = list("PGR Monoisotopic" = 1,
                                                             "PGR M+H or M-H" = 2,
                                                             "PGR Adducts (ESI + Only)" = 3,
                                                             "PGR Biotransformations" = 4),
                                              selected = 1),
                           
                           #controls for using +/- Da or +/- ppm
                           radioButtons("tolMode", "Select mass tolerance mode (Da or PPM): ",
                                        choices = list("+/- Da" = 1,
                                                       "+/- PPM" = 2),
                                        selected = 1),
                           
                  numericInput('tol', "Mass tolerance (+/- Da or PPM)", 0.02, min = 0, max = 500, step = 0.0001), #tolerance input
                  
                  fileInput('file1', 'Choose file to upload: ', #import csv button
                            accept = c(
                              'text/csv',
                              '.csv'
                            )),
           
       radioButtons("orderby", "Order By:",
                               choices = list("% RT Match" = 1,
                                              "Experimental RT" = 2,
                                              "Predicted RT" = 3,
                                              "Experimental m/z" = 4),
                               selected = 1),
                  downloadButton("downloadData", "Download Results"),
                  tags$hr(),),
                  tabPanel("Screener Output",
                           br(),
                           tableOutput('contents')),

                  tabPanel("Custom Database Search",
                           br(),
                           strong("Instructions: "),
                           p("Upload your custom dataset to search from then upload your experimental data and allow the tool to run, please allow
                             up to 3 minuites to run. Then view your results and download as a .csv"),
                           fileInput('dataset_input', 'Upload custom dataset: ', #import csv button
                                     accept = c(
                                       'text/csv',
                                       '.csv'
                                     )
                           ),
                           numericInput('tol2', "Mass tolerance (+/- Da)", 0.02, min = 0, max = 1, step = 0.0001), #slider bar input
                           fileInput('file2', 'Choose file to upload: ', #import csv button
                                     accept = c(
                                       'text/csv',
                                       '.csv'
                                     )
                           ),
                           radioButtons("orderby_shell", "Order By:",
                                        choices = list("Experimental RT" = 1,
                                                       "Predicted RT" = 2,
                                                       "% RT Match" = 3,
                                                       "Experimental m/z" = 4),
                                                        selected = 1),
                           downloadButton("downloadData_shell", "Download Results"),
                           tags$hr(),
                  ),

                  tabPanel("Custom Database Output",
                           br(),
                           tableOutput('contents_shell'),
                           ),
       
                  tabPanel("Complete DB",
                           br(),
                           DT::dataTableOutput("wholeDB")
                           ))
    )
  )
)


#####



