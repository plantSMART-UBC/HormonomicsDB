##before server stuff
##Runs dependencies required to run the application

library(shiny)
library(DT)

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

#read in the entire database for displaying 
allCompsDB = read.csv("all_comps_db.csv")
colnames(allCompsDB) = c("Compound Name", "Class", "Formula", "PRT", "Monoisotopic", "M+H")

#####
#stuff here to count how many times the script has been run

# #read in csv file from directory
# dfm <- data.frame(read.csv("counterCSV.csv"))
# 
# #convert to dataframe
# #reads in current count
# #tool must be used (below) to make the counter increase
# #actual data upload
# timesUsed <- dfm[1,2]

##server

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
    
    #takes user uploaded data and converts to local variable
    uploaded.data <- rvals$upload
    #saves the header names
    header.names <- colnames(uploaded.data)
    
    if (length(header.names > 2)) {
      
      #saves sample names from 3rd index
      sample.names <- header.names[3:length(header.names)]
      
    } else {
      
      sample.names <- c("N/A")
      
    }
    
    rvals$sample.names <- sample.names
    
    #pulls experimetal m/z values for use in search function
    rvals$expt.masses <- uploaded.data[,1]
    
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
    
    if(input$ionMode == 1) {
      
      #if ESI+ is selected by the user keep the 
      #same masses since theyre ESI+
      data.base.full[,3] <- data.base.full[,3]
      
    } else if (input$ionMode == 2) {
      
      #if ESI- is selected by the user 
      #remove 2.01505 m/z which is mass 
      #of H less an extra electron
      data.base.full[,3] <- data.base.full[,3] - 2.01505
      
      if (input$dataset == 2) {
        
        #changes the adduct to M-H
        data.base.full[,2] <- "M-H"
        
      }
      
    }
    
    #stores the 3rd column (mz) as its own vector
    data.base <- data.base.full[,3]
    
    ##search/query algorithm here
    
    #creates empty vector for markers to go into.
    #markers are where in the expt dataset a match occurs
    #based on the query options
    marker.vect <- c()
    for (i in 1:length(rvals$expt.masses)){
      
      if (input$tolMode == 1){
      
        marker <- which(rvals$expt.masses[i] >= data.base-input$tol & rvals$expt.masses[i] <= data.base+input$tol)
      
        if (length(marker[i]) > 0){
        
        #converts to useable vector
        marker.vect[[i]] <- marker
        
        }
      }
      
      else if (input$tolMode == 2){
        
        #vector for absolute ppm error to go in (will all be positive)
        ppmAbsolute <- c()
        ppmTol <- c()
        
        #calculate the error for each entry in the database
        ppmAbsolute <- data.base*(1+input$tol/10^6) 
        ppmTol <-  ppmAbsolute - data.base
        
        marker <- which(rvals$expt.masses[i] >= data.base-ppmTol & rvals$expt.masses[i] <= data.base+ppmTol)
        
        if (length(marker[i]) > 0){
          
          #converts to useable vector
          marker.vect[[i]] <- marker
          
        }
        
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
    
    # gives the matches from the data base as a vector
    location.db <- as.vector(unlist(marker.vect))
    
    # prints out everything from db where there is a match
    database.hits <- data.base.full[location.db,]
    
    #gives position of experimental hits from uploaded data
    position.experimental <- unlist(exp.match)
    
    results <- cbind(database.hits, rvals$upload[position.experimental,])
    colnames(results) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT")
    print_results <- results[,1:5]
    rvals$results <- print_results
    
    rvals$results_appended <- results
    
    ##json stuff for counting times data is user upload
    ##through the front end of the main tool
    ##only counts when the tool has actually been used!!
    
    # #adds single digit when run
    # dfm[,2] <- dfm[,2] + 1
    # 
    # #saves as a variable that is read by the program
    # rvals$timesUsed <- dfm[1,2]
    # 
    # #saves the column with new number
    # counterDF <- dfm[,2]
    # 
    # #creating a blank matrix
    # emptyMatrix <- c()
    # 
    # #overwriting matrix with blank one
    # write.csv(emptyMatrix,
    #           "counterCSV.csv")
    # 
    # #saves counter number to the server as a .csv
    # write.csv(counterDF,
    #           "counterCSV.csv")
    
  })
  
  output$contents <- renderTable({
    
    validate(
      #error message here
      need(input$file1 != "", label = "A .csv file with m/z and RT values")
    )
    
    #merges the query results to predicted RTs
    results.rt <- merge(x = rvals$results_appended[,1:5], y = PRTs,
                        by.x = "Compound Name", by.y = "Compound.Name", all.x = TRUE)
    
    colnames(results.rt) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT", "Predicted RT")
    
    #difference between expt and predicted RTs
    delta.rt <- abs(results.rt[,5] - results.rt[,6])
    
    #convert last line to a percentage
    percent.delta.rt <- 100-(delta.rt/results.rt[,5])*100
    
    #add the percent in RT to the data
    results.for.download <- cbind(results.rt, percent.delta.rt)
    
    #results.for.download <- cbind(results.for.download, ppmError)
    
    colnames(results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                        "Experimental m/z", "RT", "Predicted RT", "Percent Match RT")
    
    #concats the expt m/z and rt together
    results.for.download.mz.rt <- paste0(results.for.download[,4], "_",results.for.download[,5])
    
    #put the concat ms_rt into the data
    results.for.download <- cbind(results.for.download, results.for.download.mz.rt)
    
    colnames(results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                        "RT", "Predicted RT", "Percent Match RT", "mzrt")
    
    #merges the query results to predicted RTs
    results.for.download <- merge(x = results.for.download, y = comp.classes,
                                  by.x = "Compound Name", by.y = "Name", all.x = TRUE)
    
    colnames(results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                        "Experimental RT", "Predicted RT", "Percent Match RT", "mzrt", "Class")
    
    #loop to take UI input and sort the output data in the UI
    if (input$orderby == 1){
      
      results.for.display <- results.for.download[order(-results.for.download[,7]),]
      
    } else if (input$orderby == 2){
      
      results.for.display <- results.for.download[order(results.for.download[,5]),]
      
    } else if (input$orderby == 3){
      
      results.for.display <- results.for.download[order(results.for.download[,6]),]
      
    } else if (input$orderby == 4){
      
      results.for.display <- results.for.download[order(results.for.download[,4]),]
      
    }
    
    #takes experimental intensities and makes a local variable
    experimental.intensities <- rvals$results_appended
    
    #concat to make mz_rt
    experimental.intensities.mz.rt <- paste0(experimental.intensities[,4], "_", experimental.intensities[,5])
    
    #add mz_rt to experimental intensities
    experimental.intensities <- cbind(experimental.intensities.mz.rt, experimental.intensities[,6:ncol(experimental.intensities)])
    
    #drops duplicates (this is fine because isobars exist in rvals$results_appended)
    experimental.intensities <- unique(experimental.intensities)
    
    #change first column name to mzrt1
    colnames(experimental.intensities) <- c("mzrt1")
    
    #pseudo-SQL left join for class in download file
    download.results <- merge(x = results.for.download, y = comp.classes,
                              by.x = "Compound Name", by.y = "Name", all.x = TRUE)
    
    #pseudo-SQL left join
    download.results <- merge(x = results.for.download, y = experimental.intensities,
                              by.x = "mzrt", by.y = "mzrt1", all.x = TRUE)
    
    #converts to matrix (to drop header)
    download.results <- as.matrix(download.results)
    
    #drops header
    download.results <- matrix(download.results, ncol = ncol(download.results), dimnames = NULL)
    
    #convert back to dataframe (header = FALSE)
    download.results <- data.frame(download.results)
    
    #bring sample names back into local envrionment
    sample.names <- rvals$sample.names
    
    #column names
    colnames(download.results) <- c("mz_rt", "Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                    "RT", "Predicted RT", "Percent Match RT", "Class", sample.names)
    
    #rearranges to bring compound class in
    download.results <- download.results[,c(1,2,9,3,4,5,6,7,8,10:ncol(download.results))]
    
    #drops mz_rt
    download.results <- download.results[,2:ncol(download.results)]
    
    #convert this data to data.frame so it can be downloaded by the user
    download.results <- data.frame(download.results)
    
    #assign the data for download to a global variable so it can be downloaded within the GUI envrionment by the user
    rvals$download.noncustom <- download.results
    
    #reorders to bring class into the mix
    results.for.display <- results.for.display[, c(1, 9, 2, 3, 4, 5, 6, 7, 8)]
    
    #compute ppm mass error (don't think people will want to sort by ppm?)
    ppmError <- ((results.for.display[,4] - results.for.display[,5])/results.for.display[,4])*10^6
    
    results.for.display <- cbind(results.for.display, ppmError)
    
    colnames(results.for.display) <- c("Compound Name", "Class", "Adduct/BT", "Actual m/z", "Experimental m/z",
                                       "RT", "Predicted RT", "Percent Match RT", "ppm", "ppm Mass Difference")
    
    #displays the results we want in the GUI
    results.for.display[,c(1:8,10)]
    
  },
  digits = 5 #displays 5 decimal points
  )
  
  output$downloadData <- downloadHandler(
    filename = function(){
      #gives a unique name each time
      paste("hormonomicsDB_results", Sys.Date(), ".csv", sep = "")
    },
    content = function(file){
      #writes it to a .csv, reactive so it changes all the time
      write.csv(rvals$download.noncustom, file)
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
    #reads in the .csv file and saves to global envrionment, set header = TRUE
    rvals$upload2.custom <- read.csv(input$file2$datapath, header = TRUE,)
    
    uploaded.data <- rvals$upload2.custom
    custom.col.names <- colnames(uploaded.data[,3:ncol(uploaded.data)])
    rvals$custom.expt.masses <- uploaded.data[,1]
    
    name_adduct_mass <- rvals$name_adduct_mass
    v3 <- name_adduct_mass[,3]
    
    if (length(input$dataset_input) > 0){
      
      #create an empty vector for markers to go into
      custom.marker.vect <- c()
      
      for (i in 1:length(rvals$custom.expt.masses)){
        
        marker.custom <- which(rvals$custom.expt.masses[i] >= v3-input$tol2 & rvals$custom.expt.masses[i] <= v3+input$tol2)
        
        if (length(marker.custom[i]) > 0){
          
          custom.marker.vect[[i]] <- marker.custom
          
        }
      }
      #empty vector for position of experimental matches to go
      custom.exp.match <- c()
      for (j in 1:length(custom.marker.vect)){
        
        if (length(unlist(custom.marker.vect[j])) > 0){
          
          custom.no.match <- length(unlist(custom.marker.vect[j]))
          custom.vvv <- (rep(j, length(unlist(custom.marker.vect[j]))))
          custom.exp.match[[j]] <- custom.vvv
          
        }
      }
      
      # gives the matches from the data base as a vector
      custom.location.db <- as.vector(unlist(custom.marker.vect))
      
      # prints out everything from db where there is a match
      custom.db.hits <- name_adduct_mass[custom.location.db,]
      
      # gives the postion of experimental hits from uploaded data
      custom.pos.exp <- unlist(custom.exp.match)
      custom.results <- cbind(custom.db.hits, rvals$upload2.custom[custom.pos.exp,])
      colnames(custom.results) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT")
      custom.print_results <- custom.results[,1:5]
      rvals$custom.results <- custom.print_results
      rvals$custom.results.appended <- custom.results
    }
    
    # #adds single digit when run
    # dfm[,2] <- dfm[,2] + 1
    # 
    # #saves as a variable that is read by the program
    # rvals$timesUsed <- dfm[1,2]
    # 
    # #saves the column with new number
    # counterDF <- dfm[,2]
    # 
    # #creating a blank matrix
    # emptyMatrix <- c()
    # 
    # #overwriting matrix with blank one
    # write.csv(emptyMatrix,
    #           "counterCSV.csv")
    # 
    # #saves counter number to the server as a .csv
    # write.csv(counterDF,
    #           "counterCSV.csv")
    
  })
  
  output$contents_shell <- renderTable({
    
    validate(
      #eliminates the error message, so that the error message is much friendler
      need(input$file2 != "", label = "A .csv file with m/z and RT values")
    )
    
    #psuedo SQL left join
    custom.result.rt <- merge(x = rvals$custom.results.appended[,1:5], y = rvals$name_and_PRT,
                              by.x = "Compound Name", by.y = "Compound.Name", all.x = TRUE)
    
    #rename cols
    colnames(custom.result.rt) <- c("Compound Name", "Adduct/BT", "Actual m/z", "Experimental m/z", "RT", "Predicted RT")
    
    #computes difference between expt rt and predicted rt
    custom.delta.rt <- abs(custom.result.rt[,5] - custom.result.rt[,6])
    
    #converts to %
    custom.percent.delta.rt <- 100-(custom.delta.rt/custom.result.rt[,5])*100
    
    #binds the percent diff to the results
    custom.results.rt.delta.final <- cbind(custom.result.rt, custom.percent.delta.rt)
    
    #changes column names
    colnames(custom.results.rt.delta.final) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                                 "Experimental m/z", "RT", "Predicted RT",
                                                 "Percent Match RT", rvals$custom.col.names)
    
    #creats msrt for each hit
    custom.results.rt.delta.final.mzrt <- paste0(custom.results.rt.delta.final[,4], "_", custom.results.rt.delta.final[,5])
    
    #binds mz_rt to rest of results
    custom.results.for.download <- cbind(custom.results.rt.delta.final, custom.results.rt.delta.final.mzrt)
    
    #column names
    colnames(custom.results.for.download) <- c("Compound Name", "Adduct/BT", "Actual m/z",
                                               "Experimental m/z", "RT", "Predicted RT",
                                               "Percent Match RT", "mzrt")
    
    #logic to order results output to GUI based on users input (order by)
    if (input$orderby_shell == 1){
      
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,5]),]
      
    } else if (input$orderby_shell == 2){
      
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,6]),]
      
    } else if (input$orderby_shell == 3){
      
      shell.results.sorted <- custom.results.rt.delta.final[order(-custom.results.rt.delta.final[,7]),]
      
    } else if (input$orderby_shell == 4){
      
      shell.results.sorted <- custom.results.rt.delta.final[order(custom.results.rt.delta.final[,4]),]
      
    }
    
    #takes intensities and brings into local envrionment
    custom.experimental.intensities <- rvals$custom.results.appended
    
    #creates mz_rt for each hit
    custom.experimental.intensities.mzrt <- paste0(custom.experimental.intensities[,4],
                                                   "_", custom.experimental.intensities[,5])
    
    #binds back together
    custom.experimental.intensities <- cbind(custom.experimental.intensities.mzrt,
                                             custom.experimental.intensities[,6:ncol(custom.experimental.intensities)])
    
    #removes duplicates
    custom.experimental.intensities <- unique(custom.experimental.intensities)
    
    #renaming the mz_rt column
    colnames(custom.experimental.intensities) <- c("mzrt1")
    
    #inner join
    custom.download.results <- merge(x = custom.results.for.download, y = custom.experimental.intensities,
                                     by.x = "mzrt", by.y = "mzrt1", all.x = TRUE)
    
    #DF for downloading
    custom.download.results <- data.frame(custom.download.results)
    
    #send to global envrio for download through GUI
    rvals$download.custom.appended <- custom.download.results
    
    #output final results into GUI for user to see
    shell.results.sorted[,1:7]
    
  },
  digits = 4 #always displays 4 decimal points
  )
  
  output$downloadData_shell <- downloadHandler(
    filename = function(){
      #gives a unique name each time
      paste("hormonomicsDB_shell_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file){
      #writes it to a .csv, reactive so it changes all the time
      write.csv(rvals$download.custom.appended, file)
    }
  )
  
  # 
  output$wholeDB = DT::renderDataTable({
    
    #display whole DB
    allCompsDB
    
  })
  
}