library(shiny)
library(Biostrings)
options(repos = BiocManager::repositories())
library(ape)
library(pvclust)

options(shiny.maxRequestSize = 30*1024^2)

shinyServer(function(input, output) {
    
    observeEvent(input$fnafile,{
        output$text1 <- renderText({ 
            paste("Distance: ", input$distance)
        })
        output$text2 <- renderText({ 
            paste("Method: ", input$method)
        })
        strainnames <- c()
        pentanucs <- data.frame(matrix(rep(NA, 512), nrow=1))[numeric(0), ]
        out5table <- c()
        nr <- length(input$fnafile$datapath)
        output$text3 <- renderText({
            paste(nr," files were uploaded.")
        })
        for (fna in input$fnafile$datapath){
            dnaseq <- readDNAStringSet(fna)
            seqname <- names(dnaseq)[1]
            seqname <- sub(' genome assembly, chromosome: 1','',seqname)
            seqname <- sub(' chromosome, complete genome','',seqname)
            seqname <- sub(', complete genome','',seqname)
            seqname <- sub(' genome assembly, chromosome: I','',seqname)
            seqname <- sub(' complete genome','',seqname)
            seqname <- sub(', complete sequence','',seqname)
            seqname <- sub(' genome assembly','',seqname)
            seqname <- sub(', chromosome','',seqname)
            seqname <- sub(' chromosome','',seqname)
            seqname <- sub(' whole genome','',seqname)
            seqname <- sub(' whole genome shotgun sequence','',seqname)
            seqname <- gsub(' ','_',seqname)
            seqname <- gsub(',','_',seqname)
            if (grepl("_NODE_",seqname)){
              sqnms <- strsplit(seqname, split="_NODE_")
              seqname <- sqnms[[1]][1]
            }
            output$fnaname <- renderText(seqname)
            out5 <- oligonucleotideFrequency(dnaseq, width=5)
            out5sum <- colSums(out5)
            out5table <- rbind(out5table,out5sum)
            tmp <- cbind(names(out5), out5)
            
            tmp.df <- data.frame(tmp)
            
            pentamers <- colnames(out5)
            rlist <- c()
            for (pntm in pentamers){
                nuc5 <- DNAString(pntm)
                cmp5 <- reverseComplement(nuc5)
                rpntm <- as.character(cmp5)
                if(!pntm %in% rlist){
                    tmp.df[,pntm] <- tmp.df[,pntm] + tmp.df[,rpntm]
                    tmp.df[,rpntm] <- NULL
                    rlist[[(length(rlist) + 1)]] <- rpntm
                }
            }
            sums5 <- colSums(tmp.df)
            strainnames <- c(strainnames,seqname)
            pentanucs <- rbind(pentanucs,sums5)
        }
        colnames(pentanucs) <- colnames(tmp.df)
        rownames(pentanucs) <- strainnames
        
        output$pentamertable <- downloadHandler(
            filename = "PentamerFrequencies.csv",
            content = function(file){
                write.csv(pentanucs,file,row.names = TRUE, quote = FALSE)
            }
        )
        
        pentanucs.m <- as.matrix(pentanucs)
        pentanucs.prop <- prop.table(pentanucs.m,1)
        
        output$proptable <- downloadHandler(
            filename = "Proportionals.csv",
            content = function(file){
                write.csv(pentanucs.prop,file,row.names = TRUE, quote = FALSE)
            }
        )
        
        if(nr<10){
            plot_height = 300
        }else{
            plot_height <- as.numeric(nr)*25
        }
        
        output$PhylogeneticTree <- renderPlot({
            dist <- dist(pentanucs.prop, method=input$distance)
            cluster <- hclust(dist, method=input$method)
            pentanucs.t <- t(pentanucs)
            btsrp <- pvclust(pentanucs.t, method.dist=input$distance,nboot=1000, 
                         method.hclust=input$method,r=seq(.5,1,by=.1))
            boot <- round((btsrp$edges$bp)*100)
            plot(as.phylo(cluster))
            nodelabels(boot,adj=c(1.2, -0.2), frame="n", cex=0.7)
        })
        
        output$ui_plot <- renderUI({
            plotOutput("PhylogeneticTree",height=plot_height)
        })
        
        output$image <- downloadHandler(filename = function(){"phylogeneteicTree.pdf"},
                                               content <- function(file){
                                                   pdf(file, height=(as.numeric(nr))/2.5, 
                                                       width=10, pointsize=10)
                                                   dist <- dist(pentanucs.prop, method=input$distance)
                                                   cluster <- hclust(dist, method=input$method)
                                                   pentanucs.t <- t(pentanucs)
                                                   btsrp <- pvclust(pentanucs.t, method.dist=input$distance,nboot=1000, 
                                                                    method.hclust=input$method,r=seq(.5,1,by=.1))
                                                   boot <- round((btsrp$edges$bp)*100)
                                                   plot(as.phylo(cluster))
                                                   nodelabels(boot,adj=c(1.2, -0.2), frame="n", cex=0.7)
                                                   dev.off()
                                               })
    })
}
)
