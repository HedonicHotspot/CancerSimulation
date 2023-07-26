set.seed(1)
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Simulation Conditions
# 1 means TRUE
cancer <- 1
chemo <- 0
fasted <- 0

# Initialize values
N <- 100                             # Initial number of cells
Rp <- 1000                           # Initial population resource amount
Rc <- rnorm(N, mean = 500, sd = 20)    # Initial cell resource amounts
A <- rnorm(N, mean = .98, sd = .1)      # Altruism/Cooperation scores
Ra <- rnorm(N, mean = 8, sd = 1)       # Resources acquired for a cell each round

# Create a data frame to store the cells' data
df <- data.frame(Rc = Rc, A = A, Ra = Ra, Evolver = FALSE)
df$A[df$A>1] <- 1    # Force A>1 to one
df$A[df$A<0] <- 0    # Force A<1 to zero
df$Ra[df$Ra<0] <- 0    # Force A<1 to zero

# Initialize counters and vectors to store data for plotting
round <- 0
round_limit <- 5000     # Limit on the number of rounds to prevent infinite loop

# Store initial data for plotting
cell_pop <- c(N)
resource_pool <- c(Rp)
mean_altruism <- c(mean(A))
mean_cell_resource <- c(mean(Rc))
median_cell_resource <- c(median(Rc))
mutPop <- c(0) # Percent of cancer cells

# Initialize Treatment Loop Variable
chemoing <- 0
fasting <- 0

# Simulation loop: continues until there are no cells or reaches round limit
while (N > 0 & round < round_limit) {
  # Cancer Detection
  if(round>0 & 
    (fasted==1 | chemo == 1) &
    (fasting == 0 & chemoing ==0) &
    (mean_altruism[length(mean_altruism)] < .94)
    ){
    fasting <- 50 # Number of fasts
    if(fasted==0){chemoing <- 35} # For chemo but no fast condition
  }
  
  # Regular feeding every 5 rounds
  if(round%%5==0){Rp <- Rp + N*40}
  
  
  # Extra feed for Chemo
  if(chemoing>0 & chemo==1){Rp <- Rp + N*5}
  
  # Fasting
  if(round%%5==0 & 
     fasted == 1 & 
     fasting>0){
      Rp <- Rp - N*30
      if(Rp<1){Rp<-0}
      fasting <- fasting -1
      if(fasting==0){chemoing <-35}
  }

  df$A[df$A>1] <- 1    # Force A>1 to one
  df$A[df$A<0] <- 0    # Force A<1 to zero
  df$Ra[df$Ra<0] <- 0    # Force A<1 to zero
  
  df <- df %>%
    mutate(
      # Cells can acquire more when they are altruistic
      Ra = (A^(2))*rnorm(nrow(df), mean = abs(rnorm(100, mean=10, sd=.1))[1], sd=.1),
      # Cell satiation state
      Sc = Rc / 1000,                       
      # Population satiation state
      Sp = 1 - 1/(10^(Rp / (N*10))),        
      # Amount given by each cell
      G = Ra * A^(2) * (Sc + (1 - Sp)^2) / 2,       
      # Amount kept by each cell
      K = Ra - G,                           
      # Amount taken by each cell
      Take = (10+ (1-A^2)*Rc/10) * (sqrt(Sp) + 1 - Sc) / 2  
    )

  
  # To limit taking when pool less than amount taken
  # Takes an amount relative to other takers
  if(sum(df$Take) > Rp){
    df$Take <- Rp*(df$Take/sum(df$Take))
  }
  
  
  # Update population resource pool
  Rp <- Rp + 3*sum(df$G) - sum(df$Take) # Giving is 3x more valuable than taking
  if(Rp < 0){ Rp <-0}
  
  # Update each cell's resource
  df <- df %>%
    mutate(
      Rc = Rc + K + Take
    )
  
  # Cells that don't have enough resources die
  df <- df %>%
    filter(Rc >= 10)
  
  # Cell metabolic consumption
  df <- df %>%
    mutate(Rc = Rc - 10)    # Cell consumes resource
  
  # Introduce the evolver cell at round 20
  if (cancer==1 & round == 1050) {
    df <- rbind(df, 
                data.frame(Rc = 500, A = .1, Ra = 8, Evolver = TRUE,
                           Sc = 500/1000, Sp = df$Sp[1], 
                           G = 8*.1*(.5+1+df$Sp[1])/2,
                           K = 8 - 8*.1*(.5+1+df$Sp[1])/2,
                           Take = (10+sqrt(500)) * (1-.1^2) * (df$Sp[1]^2 + 1 - .5) / 2
                )
    )
  }
  
  
  
  # Cell Division
  dividing_cells <- df %>% filter(Rc > 1000)
  if (nrow(dividing_cells) > 0) {
    # Generate new altruism values for daughter cells
    new_A_values <- dividing_cells$A + rnorm(nrow(dividing_cells), 
                                             mean = 0, sd = 0.0005) # Evolution component
    

    
    new_A_values <- ifelse(dividing_cells$Evolver,
                           # For Evolver True
                           ifelse(new_A_values < 0, 0, 
                                  ifelse(new_A_values > 1, 1, new_A_values)),
                           # For Evolver False
                           # Ensure scores are between 0 and 1 for evolvers
                           ifelse(new_A_values < 0, 0, 
                                  ifelse(new_A_values > 1, 1, new_A_values))) 
    

    
    # Create new cells with 500 resources and new altruism values
    dividing_cells <- dividing_cells %>%
      mutate(
        Rc = 500,
        A = new_A_values
      )
    dividing_cells <- dividing_cells[rep(row.names(dividing_cells), each = 2),]

    
    # Remove Parent Cells
    df <- df %>%
      filter(Rc <= 1000) %>%
      # Normal Cells die if low altruism
      filter( (Evolver ==1 ) |  
                (Evolver ==0 & A>.96) )
      
    # Chemotherapy
    # Random 95% of Dividing Cells Removed
    if((chemo==1) &
      (round !=0) &
      (chemoing>0)){
      
      if(chemoing>0){
        dividing_cells <-
          dividing_cells[sample(1:nrow(dividing_cells),
                                .05*nrow(dividing_cells)),]
      }
      chemoing <- chemoing -1
    }
    
    # Update cell data
      df <- rbind(df, dividing_cells)
  }

  
  # Cell Turnover
  
  if(round> 1200 & round%%2==0){
    # Oldest Undivided Cells Die First
    # keep <- floor(tanh(nrow(df)/(3*5000000))*nrow(df))
    # df <- df[keep:nrow(df),] 
    # Random Death
    keep <- 1-tanh(nrow(df)/(3*1000000))
    df <- df[sample(1:nrow(df),keep*nrow(df)),]
    # df <- df[sample(1:nrow(df),.99999*nrow(df)),]
  }
  
  # Update the number of cells
  N <- nrow(df)
  
  # System Status Check
  if(round != 0){
    if(round %% 200 == 0){
      print(summary(mean_altruism))
      print(paste("Round:", round))
      print(paste("Cell Population:",nrow(df)))
      print(summary(df$Evolver))
    }
  }
  
  
  # Store data for plotting
  if(N > 0 & round < round_limit){
    rownames(df) <- 1:N
    round <- round + 1
    cell_pop <- c(cell_pop, N)
    resource_pool <- c(resource_pool, Rp)
    mean_altruism <- c(mean_altruism, mean(df$A))
    mean_cell_resource <- c(mean_cell_resource, mean(df$Rc))
    median_cell_resource <- c(median_cell_resource, median(df$Rc))
    mutPop <- c(mutPop, nrow(df[df$Evolver,])/nrow(df))
  }
}

# Normalize the data by dividing each series by its maximum value
normalized_mean_altruism <- mean_altruism #/ max(mean_altruism)
normalized_resource_pool <- (resource_pool/cell_pop) / max(resource_pool/cell_pop)
normalized_cell_pop <- cell_pop / max(cell_pop)
normalized_mean_cell_resource <- mean_cell_resource / max(mean_cell_resource)

if(!cancer){plotTitle <- "Healthy System"}
if(cancer){plotTitle <- "No Treatment"}
if(fasted){plotTitle <- "Fasted"}
if(chemo){plotTitle <- "Chemotherapy"}
if(chemo & fasted){plotTitle <- "Chemotherapy and Fasted"}

# Plotting - Overlay of Normalized Mean Altruism, Population Resource Pool, Cell Population, and Mean Cell Resource
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(0:round, normalized_mean_altruism, type = "l", col = "blue", xlab = "Round", ylab = "", ylim=c(0, 1.2),
     main = plotTitle)
abline(h=.7, col="black", lty=2)
lines(0:round, mutPop, type = "l", col = "red")
lines(0:round, normalized_cell_pop, type = "l", col = "purple")
lines(0:round, normalized_mean_cell_resource, type = "l", col = alpha("orange", .8), lwd = .6)
# legend("topright", legend=c("Mean Altruism", "Population Resource Pool", "Cell Population", "Mean Cell Resource"),
#       col=c("blue", "green", "purple", "orange"), lty=1, cex=0.5)
