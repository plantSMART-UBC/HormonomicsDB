#hormonomicsDB v1.4.1
#July 21st 2021
#Authors: Ryland T. Giebelhaus, Lauren A.E. Erland, and Susan J. Murch
#PlantSMART research group at UBC Okanagan
#contact: Dr. Susan J. Murch. Email: susan.murch@ubc.ca
#Website: hormonomicsDB.com

#call shiny library in
library(shiny)

HORMONOMICSDBV11_adducts <- read.csv("adducts.csv") #read in the melted data
hormcsv_adducts <- data.frame(HORMONOMICSDBV11_adducts) #change name and format
colnames(hormcsv_adducts) <- c("Names", "Adduct", "mz") #gives all the column names the same names
v3_adducts <- hormcsv_adducts[,3]

HORMONOMICSDBV11_bts <- read.csv("biotransformations.csv") #read in the melted data
hormcsv_bts <- data.frame(HORMONOMICSDBV11_bts) #change name and format
colnames(hormcsv_bts) <- c("Names", "Adduct", "mz") #gives all the column names the same names
v3_bts <- hormcsv_bts[,3]

HORMONOMICSDBV11_MH <- read.csv("M_plus_H_adduct.csv") #read in the melted data
hormcsv_MH <- data.frame(HORMONOMICSDBV11_MH) #change name and format
colnames(hormcsv_MH) <- c("Names", "Adduct", "mz") #gives all the column names the same names
v3_MH <- hormcsv_MH[,3]

HORMONOMICSDBV11_monoiso <- read.csv("monoisotopic.csv") #read in the melted data
hormcsv_monoiso <- data.frame(HORMONOMICSDBV11_monoiso) #change name and format
colnames(hormcsv_monoiso) <- c("Names", "Adduct", "mz") #gives all the column names the same names
v3_monoiso <- hormcsv_monoiso[,3]

PRTs <- read.csv("RTs.csv")
PRTs <- data.frame(PRTs)

comp.classes <- read.csv("comp_classes.csv")
comp.classes <- data.frame(comp.classes)
colnames(comp.classes) <- c("Name", "Class")
#####
#json stuff here to count how many times the script has been run
library(rjson)

#read in json file from directory
jsonCOUNTER <- fromJSON(file = "counter.json")

#convert to dataframe
#reads in current count
#tool must be used (below) to make the counter increase
#actual data upload
jsonDF <- as.data.frame(jsonCOUNTER)
timesUsed <- jsonDF[1,1]

#####

ui <- fluidPage(

  titlePanel("HormonomicsDB"),

  sidebarLayout(
    sidebarPanel(
      strong("HormonomicsDB v1.4.0"),
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
                            'HormonomicsDB shell' to upload your own dataset use our platform to perform your
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
                             Dispite diligent validation and bug fixing, we are not responsible for any mistakes the application
                             makes in data processing. Considering this, please inform us immediately of any bugs that you encounter."),
                           p("We do not save any data that is uploaded to the server, it is
                             immediatley deleted with every new session that you start."),
                           p("Please acknowledge the aforementioned authors in any work where HormonomicsDB has been used."),
                           ),

                  tabPanel("M/Z Screener",

                           br(),

                           strong("Instructions: "),
                           p("Select which datasets to search from then select a search tolerance
                             and then upload your formatted data as a .csv and allow up to 3 minutes to perform
                             the search. After this is completed select how you want your data ordered and
                             view it in the 'Screener Output' tab."),

                           #gives checkboxes so user can select multiple datasets to search from at once
                           checkboxGroupInput("dataset", "Choose Dataset: ",
                                              choices = list("PGR Monoisotopic" = 1,
                                                             "PGR M+H" = 2,
                                                             "PGR Adducts" = 3,
                                                             "PGR Biotransformations" = 4),
                                              selected = 1),

                  numericInput('tol', "Mass tolerance (+/- Da)", 0.02, min = 0, max = 1, step = 0.0001), #tolerance input
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

                  tabPanel("HormonomicsDB Shell",
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

                  tabPanel("Shell Output",
                           br(),
                           tableOutput('contents_shell'),
                           ))
    )
  )
)


#####

server <- function(input, output) {

  #set max upload size to 10 Mb to prevent crashing
  options(shiny.maxRequestSize = 10*1024^2)

  #initalizes reactive (public) variables
  rvals <- reactiveValues(
    upload = NULL,
    expt.masses = NULL,
    results = NULL,
    results_appended = NULL,
    results_appended_1 = NULL,
    upload2.custom = NULL,
    custom.expt.masses = NULL,
    custom.results = NULL,
    custom.results.appended = NULL,
    custom.results.appended_1 = NULL,
    userdata = NULL,
    name_adduct_mass = NULL,
    v3 = NULL,
    name_and_PRT = NULL,
    custom_RTS = NULL,
    custom.col.names = NULL,
    download.custom.appended = NULL,
    download.noncustom = NULL,
    sample.names = NULL,
    input_vector = NULL,
    timesUsed = NULL
  )

  observeEvent(input$file1,{
    req(input$file1)

    #reads in the .csv file and saves to global envrionment. sets header to TRUE
    rvals$upload <- read.csv(input$file1$datapath, header = TRUE,)

    uploaded.data <- rvals$upload #takes user uploaded data and converts to local variable
    header.names <- colnames(uploaded.data) #saves the header names

    if (length(header.names > 2)) {
    sample.names <- header.names[3:length(header.names)]
    #saves sample names from 3rd index
    }
    else {
      sample.names <- c("N/A")
    }

    rvals$sample.names <- sample.names

    rvals$expt.masses <- uploaded.data[,1] #pulls experimental m/z's for use in search function

    #dataset selection algorithm
    #takes the users inputs (as a vector) and stores as a reactive var
    rvals$input_vector <- input$dataset

    #stores as a local var and converts to numeric (ints)
    user_databases <- as.numeric(rvals$input_vector)

    #empty list for sticking the databases into
    datalist_databases <- list()

    #for loop to operate a switch based
    #on the users selected inputs for databases
    for (i in 1:length(user_databases)) {

      #itterates over the selected options via a switch
      x <- switch(user_databases[i],
                  #add more databases here, they have to correspond
                  #to the order that they appear above in the
                  #user interface and must be EXACTLY
                  #the same format!
                  hormcsv_monoiso,
                  hormcsv_MH,
                  hormcsv_adducts,
                  hormcsv_bts)

      #stores as a datalist
      datalist_databases[[i]] <- x

    }

    #rbinds the selected databases into one dataframe
    data.base.full <- do.call(rbind, datalist_databases)

    #ensures that the column names are the same
    colnames(data.base.full) <- c("Names", "Adduct", "mz")

    #stores the 3rd column (mz) as its own vector
    data.base <- data.base.full[,3]

    ##search/query algorithm here
      marker.vect <- c() #create an empty vector for markers to go into, markers are where in the expt dataset a match occurs
      for (i in 1:length(rvals$expt.masses)){
        marker <- which(rvals$expt.masses[i] >= data.base-input$tol & rvals$expt.masses[i] <= data.base+input$tol)
        if (length(marker[i]) > 0){
          marker.vect[[i]] <- marker #converts to useable vector
        }
      } #loop to find markers (where in expt dataset match occurs)
      exp.match <- c() # empty vector for position of experimental matches to go
      for (j in 1:length(marker.vect)){
        if (length(unlist(marker.vect[j])) > 0){
          no.match <- length(unlist(marker.vect[j]))
          replicate.no.match <- (rep(j, length(unlist(marker.vect[j])))) #replicates for instances of j
          exp.match[[j]] <- replicate.no.match
        }
      }
      location.db <- as.vector(unlist(marker.vect)) # gives the matches from the data base as a vector
      database.hits <- data.base.full[location.db,] # prints out everything from db where there is a match
      position.experimental <- unlist(exp.match) # gives the postion of experimental hits from uploaded data
      results <- cbind(database.hits, rvals$upload[position.experimental,])
      colnames(results) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT")
      print_results <- results[,1:5]
      rvals$results <- print_results
      ##if ()
      rvals$results_appended <- results

      ##json stuff for counting times data is user upload
      ##through the front end of the main tool
      ##only counts when the tool has actually been used!!

      #adds single digit when run
      jsonDF[1,1] <- jsonDF[1,1] + 1
      rvals$timesUsed <- jsonDF[1,1]

      #convert back to JSON file
      jsonNewCount <- toJSON(jsonDF)

      #writes as JSON
      write(jsonNewCount, "counter.json")

  })

  output$contents <- renderTable({

    validate(
      need(input$file1 != "", label = "A .csv file with m/z and RT values")
    ) #eliminates the error message, so that the error message is much friendler


    results.rt <- merge(x = rvals$results_appended[,1:5], y = PRTs,
                        by.x = "Compound Name", by.y = "Compound.Name", all.x = TRUE) #merges the query results to predicted RTs
    colnames(results.rt) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT", "Predicted RT")
    delta.rt <- abs(results.rt[,5] - results.rt[,6]) #difference between expt and predicted RTs
    percent.delta.rt <- 100-(delta.rt/results.rt[,5])*100 #convert last line to a percentage
    results.for.download <- cbind(results.rt, percent.delta.rt) #add the percent in RT to the data
    colnames(results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                        "Experimental m/z", "RT", "Predicted RT", "Percent Match RT")
    results.for.download.mz.rt <- paste0(results.for.download[,4], "_",results.for.download[,5]) #concatanate the expt mz and rt together
    results.for.download <- cbind(results.for.download, results.for.download.mz.rt) #put the concat mz_rt into the data
    colnames(results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                        "RT", "Predicted RT", "Percent Match RT", "mzrt")
    results.for.download <- merge(x = results.for.download, y = comp.classes,
                        by.x = "Compound Name", by.y = "Name", all.x = TRUE) #merges the query results to predicted RTs
    colnames(results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                        "Experimental RT", "Predicted RT", "Percent Match RT", "mzrt", "Class")

    if (input$orderby == 1){
      results.for.display <- results.for.download[order(-results.for.download[,7]),]
    }
    else if (input$orderby == 2){
      results.for.display <- results.for.download[order(results.for.download[,5]),]
    }
    else if (input$orderby == 3){
      results.for.display <- results.for.download[order(results.for.download[,6]),]
    }
    else if (input$orderby == 4){
      results.for.display <- results.for.download[order(results.for.download[,4]),]
    } #loop to take UI input and sort the output data in the UI

    experimental.intensities <- rvals$results_appended #takes experimental intensities and makes a local variable
    experimental.intensities.mz.rt <- paste0(experimental.intensities[,4], "_", experimental.intensities[,5]) #concat to make mz_rt
    experimental.intensities <- cbind(experimental.intensities.mz.rt, experimental.intensities[,6:ncol(experimental.intensities)]) #add mz_rt to experimental intensities
    experimental.intensities <- unique(experimental.intensities) #drops duplicates (this is fine because isobars exist in rvals$results_appended)
    colnames(experimental.intensities) <- c("mzrt1") #change first column name to mzrt1

    download.results <- merge(x = results.for.download, y = comp.classes,
                              by.x = "Compound Name", by.y = "Name", all.x = TRUE) #pesudo-SQL left join for class in download file
    download.results <- merge(x = results.for.download, y = experimental.intensities,
                              by.x = "mzrt", by.y = "mzrt1", all.x = TRUE) #pseudo-SQL left join
    download.results <- as.matrix(download.results) #converts to matrix (to drop header)
    download.results <- matrix(download.results, ncol = ncol(download.results), dimnames = NULL) #drops header
    download.results <- data.frame(download.results) #convert back to dataframe (header = FALSE)
    sample.names <- rvals$sample.names #bring sample names back into local envrionment
    colnames(download.results) <- c("mz_rt", "Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                    "RT", "Predicted RT", "Percent Match RT", "Class", sample.names) #column names
    download.results <- download.results[,c(1,2,9,3,4,5,6,7,8,10:ncol(download.results))] #rearranges to bring compound class in
    download.results <- download.results[,2:ncol(download.results)] #drops mz_rt
    download.results <- data.frame(download.results) #convert this data to data.frame so it can be downloaded by the user

    rvals$download.noncustom <- download.results #assign the data for download to a global variable so it can be downloaded within the GUI envrionment by the user

    results.for.display <- results.for.display[, c(1, 9, 2, 3, 4, 5, 6, 7, 8)] #reorders to bring class into the mix
    results.for.display[,1:8] #displays the results we want in the GUI

  },
  digits = 4 #displays 4 decimal points
  )

  output$downloadData <- downloadHandler(
    filename = function(){
      paste("hormonomicsDB_results", Sys.Date(), ".csv", sep = "") #gives a unique name each time
    },
    content = function(file){
      write.csv(rvals$download.noncustom, file) #writes it to a .csv, reactive so it changes all the time
    }
  )

  observeEvent(input$dataset_input,{
    req(input$dataset_input)
    rvals$userdata <- read.csv(input$dataset_input$datapath, header = TRUE,)
    userdata <- rvals$userdata
    colnames(userdata) <- c("Compound.Name", "Adduct", "Mass", "RTP")

    rvals$name_adduct_mass <- userdata[,1:3]

    custom.unique.compounds <- length(unique(userdata[,1]))
    rvals$name_and_PRT <- userdata[1:custom.unique.compounds,]
    rvals$name_and_PRT <- rvals$name_and_PRT[,c(1,4)]


  })

  observeEvent(input$file2,{
    req(input$file2)
    rvals$upload2.custom <- read.csv(input$file2$datapath, header = TRUE,) #reads in the .csv file and saves to global envrionment, set header = TRUE

    uploaded.data <- rvals$upload2.custom
    custom.col.names <- colnames(uploaded.data[,3:ncol(uploaded.data)])
    rvals$custom.expt.masses <- uploaded.data[,1]

    name_adduct_mass <- rvals$name_adduct_mass
    v3 <- name_adduct_mass[,3]

    if (length(input$dataset_input) > 0){
      custom.marker.vect <- c() #create an empty vector for markers to go into
      for (i in 1:length(rvals$custom.expt.masses)){
        marker.custom <- which(rvals$custom.expt.masses[i] >= v3-input$tol2 & rvals$custom.expt.masses[i] <= v3+input$tol2)
        if (length(marker.custom[i]) > 0){
          custom.marker.vect[[i]] <- marker.custom
        }
      }
      custom.exp.match <- c() # empty vector for position of experimental matches to go
      for (j in 1:length(custom.marker.vect)){
        if (length(unlist(custom.marker.vect[j])) > 0){
          custom.no.match <- length(unlist(custom.marker.vect[j]))
          custom.vvv <- (rep(j, length(unlist(custom.marker.vect[j]))))
          custom.exp.match[[j]] <- custom.vvv
        }
      }
      custom.location.db <- as.vector(unlist(custom.marker.vect)) # gives the matches from the data base as a vector
      custom.db.hits <- name_adduct_mass[custom.location.db,] # prints out everything from db where there is a match
      custom.pos.exp <- unlist(custom.exp.match) # gives the postion of experimental hits from uploaded data
      custom.results <- cbind(custom.db.hits, rvals$upload2.custom[custom.pos.exp,])
      colnames(custom.results) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT")
      custom.print_results <- custom.results[,1:5]
      rvals$custom.results <- custom.print_results
      rvals$custom.results.appended <- custom.results
    }

    ##json stuff for counting times data is user upload
    ##through the front end of the custon queue tool
    ##only counts when the tool has actually been used!!

    #adds single digit when run
    jsonDF[1,1] <- jsonDF[1,1] + 1
    rvals$timesUsed <- jsonDF[1,1]

    #convert back to JSON file
    jsonNewCount <- toJSON(jsonDF)

    #writes as JSON
    write(jsonNewCount, "counter.json")

  })

  output$contents_shell <- renderTable({

    validate(
      need(input$file2 != "", label = "A .csv file with m/z and RT values")
    ) #eliminates the error message, so that the error message is much friendler


    custom.result.rt <- merge(x = rvals$custom.results.appended[,1:5], y = rvals$name_and_PRT,
                              by.x = "Compound Name", by.y = "Compound.Name", all.x = TRUE) #psuedo SQL left join
    colnames(custom.result.rt) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT", "Predicted RT") #rename cols
    custom.delta.rt <- abs(custom.result.rt[,5] - custom.result.rt[,6]) #computes difference between expt rt and predicted rt
    custom.percent.delta.rt <- 100-(custom.delta.rt/custom.result.rt[,5])*100 #converts to percentage
    custom.results.rt.delta.final <- cbind(custom.result.rt, custom.percent.delta.rt) #binds the percent diff to the results
    colnames(custom.results.rt.delta.final) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                                 "Experimental m/z", "RT", "Predicted RT",
                                                 "Percent Match RT", rvals$custom.col.names) #changes column names
    custom.results.rt.delta.final.mzrt <- paste0(custom.results.rt.delta.final[,4], "_", custom.results.rt.delta.final[,5]) #creates mz_rt for each hit
    custom.results.for.download <- cbind(custom.results.rt.delta.final, custom.results.rt.delta.final.mzrt) #binds mz_rt to rest of results
    colnames(custom.results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                               "Experimental m/z", "RT", "Predicted RT",
                                               "Percent Match RT", "mzrt") #column names

    if (input$orderby_shell == 1){
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,5]),]
    }else if (input$orderby_shell == 2){
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,6]),]
    }else if (input$orderby_shell == 3){
      shell.results.sorted <- custom.results.rt.delta.final[order(-custom.results.rt.delta.final[,7]),]
    }else if (input$orderby_shell == 4){
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,4]),]
    } #logic to order results output to GUI based on users input (order by)

    custom.experimental.intensities <- rvals$custom.results.appended #takes intensities and brings into local envrionment
    custom.experimental.intensities.mzrt <- paste0(custom.experimental.intensities[,4],
                                                   "_", custom.experimental.intensities[,5]) #creates mz_rt for each hit
    custom.experimental.intensities <- cbind(custom.experimental.intensities.mzrt,
                                             custom.experimental.intensities[,6:ncol(custom.experimental.intensities)]) #binds back together
    custom.experimental.intensities <- unique(custom.experimental.intensities) #removes duplicates
    colnames(custom.experimental.intensities) <- c("mzrt1") #renaming the mz_rt column

    custom.download.results <- merge(x = custom.results.for.download, y = custom.experimental.intensities,
                                     by.x = "mzrt", by.y = "mzrt1", all.x = TRUE) #inner join
    custom.download.results <- data.frame(custom.download.results) #df for downloading
    rvals$download.custom.appended <- custom.download.results #send to global envrionment so it can be downloaded from GUI

    shell.results.sorted[,1:7] #output final results into GUI for user to see

  },
  digits = 4 #always displays 4 decimal points
  )

  output$downloadData_shell <- downloadHandler(
    filename = function(){
      paste("hormonomicsDB_shell_results_", Sys.Date(), ".csv", sep = "") #gives a unique name each time
    },
    content = function(file){
      write.csv(rvals$download.custom.appended, file) #writes it to a .csv, reactive so it changes all the time
    }
  )


  output$timesRUN <- renderText({

    #prints the number of times the app has been used succesfully

    paste("The app has been run ", timesUsed, " times.")

  })

}

#make the server talk to the ui
shinyApp(ui = ui, server = server)
