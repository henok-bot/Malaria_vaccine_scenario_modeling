
### Malaria model with out vaccine # Dera woreda -------------

#-------------------------------------------------
# Load required package
#-------------------------------------------------
library(readr) # For importing data
library(deSolve)  # For solving differential equations numerically
library(dplyr) # For data wrangling
library(ggplot2)
#-------------------------------------------------
# Load data
#-------------------------------------------------

Malaria <- read_csv("weekly_malaria_amh.csv", 
                    col_types = cols(date = col_date(format = "%m/%d/%Y")))


#-------------------------------------------------
# Data wrangling
#-------------------------------------------------
# Combine woredas with base package
Malaria_woredas <- aggregate(total_cases ~ woreda + year + epi_week, data = Malaria, sum) 

#make our data subset to fit
Dera <- with(Malaria_woredas, subset(cbind(total_cases, epi_week, year), woreda == 'Dera (AM)'))
Dera <- as.data.frame(Dera)
Dera <- arrange(Dera, year, epi_week)
Dera$C_cases <- c(Dera$total_cases[1],rep(NA,times=dim(Dera)[1]-1))
for (i in 1:(dim(Dera)[1]-1)){
  Dera$C_cases[i+1] <- Dera$total_cases[i+1]+Dera$C_cases[i]
}

plot(Dera$epi_week,Dera$total_cases)




#-------------------------------------------------
# Define the malaria model function
#-------------------------------------------------
mal_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S_H + I_H + R_H  # Total human population 
    M <- S_V + I_V        # Total adult vector population 
    
    # Force of infection: 
    lambda_h <- alpha_vh * beta_0*(1+psi*sin(2*pi*time/52)) * (I_V / M)   # mosquito->human
    lambda_v <- alpha_hv * beta_0*(1+psi*sin(2*pi*time/52)) * (I_H / N)   # human->mosquito
    
    #--------------------------------
    # Differential equations:
    #--------------------------------
    # Humans (SIR)
    dS_H <-  b*N - lambda_h * S_H - mu * S_H  + tau * R_H       # Susceptible individuals
    dI_H <- lambda_h* S_H - (mu + sigma + gamma) * I_H  # Infected individuals
    dR_H <- gamma * I_H - mu * R_H - tau * R_H                      # Recovered individuals
    dC <-  lambda_h* S_H# cumulative cases
    
    # Mosquitoes (SI)
    dS_V <- m * M - lambda_v * S_V - mu_v * S_V          # Susceptible individuals
    dI_V <- lambda_v * S_V - mu_v * I_V  # Infected individuals
    list(c(dS_H, dI_H, dR_H,  dC, dS_V, dI_V))  # Return the system of equations
    
  })
}

#--------------------------------------------------
# Epidemiological constants
#--------------------------------------------------
human_life_expectancy = 65*52 # Average human life expectancy(weeks)
mosquito_lifespan = 2 # (days converted to weeks) # changed
infectious_period = 9/7 # (days converted to weeks)

#---------------------------------------------------
# Model Parameters (Not estimated)
#---------------------------------------------------
b     = 0.000465  # Birth rate of humans "Crude birth rate is 24.2 and crude death rate is 5.7 per 1000 population per year"
m   = 0.254144  # Recruitment rate of vectors
mu    = 1/ (human_life_expectancy)    # Natural human death rate
mu_v  = 1/ (mosquito_lifespan)# Mortality rate of adult vectors
#sigma = 0.0000   # Disease-induced death rate (assuming insignificant death)
beta_0 = 26 # bites per week
#psi = 0.118557052 #amplitude
gamma = 1/ (infectious_period)      # Recovery rate
tau = 1/(126/7) # rate of waning immunity from recovered
#alpha_vh = 0.108575626 # Transmission probability from vectors to hosts
#alpha_hv = 0.118557052 # Transmission probability from hosts to vectors

#-------------------------------------------------
# Define model parameters (estimated)
#-------------------------------------------------
parameters <- c( 
  sigma = 0.001,   # Disease-induced death rate
  alpha_hv = 0.185, # 0.185
  alpha_vh = 0.092, # 0.092
  psi = 0.5 #amplitude # 0.5
)

#-------------------------------------------------
# Initial state values for each compartment at t = 0
#-------------------------------------------------
c0 <- 631
initial_conditions <- c(
  S_H = 305191,  # Initial number of susceptible humans
  I_H = 631,   # Initial number of infected humans
  R_H = 13081,    # Initial number of recovered humans
  C = c0, # cumlative cases
  S_V = 500000,  # Initial number of susceptible mosquitoes
  I_V = 443   # Initial number of infected mosquitoes
)
#A<-subset(Dera$total_cases, Dera$year == 2019 & Dera$epi_week >=34)
#sum(A)

#-------------------------------------------------
# Time points for simulation
#-------------------------------------------------
t <- seq(1, dim(Dera)[1]-52, by = 1)  # Simulate for 2020-2024 week time step

#-------------------------------------------------
# Solve the system of differential equations
#-------------------------------------------------
estimateMalaria <- function(initial_conditions, t, parameters, data) {
  model_output <- ode(y = initial_conditions, times = t,
                      func = mal_model, parms = parameters)
  df <- as.data.frame(model_output)
  inc <- c(c0,diff(df$C))
  #return(df)
  llik<-sum(dpois(c(data),inc,log=TRUE))
  #llik<-sum(dpois(c(data),df$C,log=TRUE))
  return(-llik)
}
#inc <- c(c0,diff(df$C))
#print(inc)

#estimateMalaria(initial_conditions, t, parameters, data=Dera$total_cases[53:313]) #see what the function does
#estimateMalaria(initial_conditions, t, parameters, data=Dera$C_cases[53:313])

#-------------------------------------------------
# Estimate the parameters
#-------------------------------------------------
NusianceFunction <-function(parameters){
  estimateMalaria(initial_conditions, t, parameters, data=Dera$total_cases[53:313])
  #estimateMalaria(initial_conditions, t, parameters, data=Dera$C_cases[53:313])
}

parametersEst<-optim(par=parameters,fn=NusianceFunction)
parametersEst$par

#-------------------------------------------------
# View the Result
#-------------------------------------------------
model_output <- ode(y = initial_conditions, times = t,
                    func = mal_model, parms = parametersEst$par)
#func = mal_model, parms = parameters)
df <- as.data.frame(model_output)
df$inc <- c(c0,diff(df$C))
## Compare Model and Data
matplot(
  df$time,
  cbind(df$inc, Dera$total_cases[53:313]),
  #cbind(df$C, Dera$C_cases[53:313]),
  #Dera$total_cases[53:313],
  type = "l",            # Line plot
  lty = 1,               # Solid lines
  lwd = 2,               # Line width
  col = c("blue", "green"),  # Colors for estimated outbreak and data
  #col = "green",
  xlab = "Time (weeks)", # xlabel
  ylab = "Number of Individuals", # ylabel
  main = "Malaria in Dera 2020-2024" # title
)

# Add legend to the plot
legend(
  "topleft",
  legend = c("Model Estimate", "Data"),
  col = c("blue", "green"),
  lty = 1,
  lwd = 2
)



## Compare Model and Data
#-------------------------------------------------
# Plot results using ggplot2 for clarity and aesthetics
#-------------------------------------------------
## Graph for Humans
ggplot(df, aes(x = t, y = c(I_H, R_H, V_H))) +
  geom_line(aes(y = S_H, color = "S_H"), lwd = 1) +
  #geom_line(aes(y = E_H, color = "E_H"), lwd = 1) +
  geom_line(aes(y = I_H, color = "I_H"), lwd = 1) +
  geom_line(aes(y = R_H, color = "R_H"), lwd = 1) +
  #geom_line(aes(y = V_H, color = "V_H"), lwd = 1)
  # Add manual color scales for better control and clarity
  scale_color_manual(
    values = c(
      "S_H" = "brown",
      "I_H"    = "red",
      "R_H"   = "green"
      #"V_H" = "blue",
      
    ),
    labels = c("Human Infectious", "Human Recovered", "Suceptible popn") # "Human Vaccinated") #, "Human Susceptible")
  ) +
  # Add title and axes labels
  labs(
    title = "Malaria Model Simulation with Juveniles",
    subtitle = "Population Dynamics over Time",
    x = "Time",
    y = "Population Size",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "right")

## Graph for Vectors
ggplot(df, aes(x = t, y = c(S_V, I_V))) +
  #geom_line(aes(y = A, color = "A"), lwd = 1) +
  geom_line(aes(y = S_V, color = "S_V"), lwd = 1) +
  #geom_line(aes(y = E_V, color = "E_V"), lwd = 1) +
  geom_line(aes(y = I_V, color = "I_V"), lwd = 1) +
  # Add manual color scales for better control and clarity
  scale_color_manual(
    values = c(
      "S_V" = "navy",
      "I_V" = "darkred"
    ),
    labels = c( "Vector Infected", 
               "Vector Susceptible")
  ) +
  # Add title and axes labels
  labs(
    title = "Malaria Model Simulation Vector population",
    subtitle = "Vector Population Dynamics over Time",
    x = "Time",
    y = "Population Size",
    color = NULL
  ) +
  theme_minimal() +
  theme(legend.position = "right")





# Trying to calculate the R square

# Time and cumulative cases
t_data <- seq(1, dim(Dera)[1]-52, by = 1) 
observed_cum_cases <- cumsum(Dera$total_cases[Dera$year>2019])

objective_func <- function(tspan, parameters) {
  model_output <- ode(
    y = initial_conditions, 
    times = t_data, 
    func = mal_model, 
    parms = parametersEst$par
  )
  cum_cases_model <- model_output[, "C"]
  return(cum_cases_model)
}


predicted_cum_cases <- objective_func(t_data, optimal_parameters)
SS_res <- sum((observed_cum_cases - predicted_cum_cases)^2)
SS_tot <- sum((observed_cum_cases - mean(observed_cum_cases))^2)
R_squared <- 1 - (SS_res / SS_tot)

# Print R-squared value
cat("R-squared:", R_squared, "\n")
