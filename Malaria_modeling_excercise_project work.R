
# Malaria with vaccine ----------------

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



#-------------------------------------------------
# Define the malaria model function
#-------------------------------------------------
mal_model_vac <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    N <- S_H + I_H + R_H + V_H # Total human population 
    M <- S_V + I_V        # Total adult vector population 
 
    # Force of infection: 
    lambda_h <- alpha_vh * beta_0*(1+psi*sin(2*pi*time/52)) * (I_V / M)   # mosquito->human
    lambda_v <- alpha_hv * beta_0*(1+psi*sin(2*pi*time/52)) * (I_H / N)   # human->mosquito
    
    #--------------------------------
    # Differential equations:
    #--------------------------------
    # Humans (SIR)
    dS_H <-  b*N - lambda_h * S_H - mu * S_H  + tau * R_H - (omega * S_H * tau_2)      # Susceptible individuals
    dI_H <- lambda_h* S_H - (mu + sigma + gamma) * I_H  # Infected individuals
    dR_H <- gamma * I_H - mu * R_H - tau * R_H                      # Recovered individuals
    dV_H <- (omega*S_H*tau_2) - mu * V_H 
    dC <-  lambda_h* S_H# cumulative cases
    
    # Mosquitoes (SI)
    dS_V <- m * M - lambda_v * S_V - mu_v * S_V          # Susceptible individuals
    dI_V <- lambda_v * S_V - mu_v * I_V  # Infected individuals
    list(c(dS_H, dI_H, dR_H,dV_H,  dC, dS_V, dI_V))  # Return the system of equations
   
  })
}

#--------------------------------------------------
# Epidemiological constants
#--------------------------------------------------
human_life_expectancy = 65*52 # Average human life expectancy(weeks)
mosquito_lifespan = 2 # (days converted to weeks) # changed
infectious_period = 9/7 # (days converted to weeks)

#---------------------------------------------------
# Model Parameters (Not estimated) - sourced from the vaccine model
#---------------------------------------------------
b     = 0.000465 # Birth rate of humans "Crude birth rate is 24.2 and crude death rate is 5.7 per 1000 population per year"
#phi = 12.8 # oviposition rate for eggs in water
#K = 1000000000 # carrying capacity of environment
m   = 0.254144  # Recruitment rate of vectors
#f = 0.5 #proportion of adults that are female
mu    = 1/ (human_life_expectancy)    # Natural human death rate
mu_v  = 1/ (mosquito_lifespan)# Mortality rate of adult vectors
#mu_a = 0.13# Mortality rate of aquatic stages of vectors
sigma = -0.641844026  # Disease-induced death rate ??
beta_0 = 26 # bites per week
psi = 0.336561225 #amplitude
gamma = 1/ (infectious_period)      # Recovery rate
tau = 1/(126/7) # rate of waning immunity from recovered
tau_2 = 3/4 # waning immunity from vaccine vaccine efficacy is 75% so 25% will be suseptible
alpha_vh = 0.008912758 # Transmission probability from vectors to hosts
alpha_hv = 0.007039342 # Transmission probability from hosts to vectors
omega = 0.10/52  # vaccination rate (4% of total population to be vaccinated annualy - divided by 52 to convert weekly)
#kappa = 0.061551*7 #Incubation rate in mosquitoes




# Bring the model here 

parameters <- c(b, mu, m, mu_v, sigma, beta_0,
                psi, gamma, tau, tau_2, alpha_hv,
                alpha_hv, omega)

c0 <- 631
initial_conditions <- c(
  S_H = 305191,  # Initial number of susceptible humans
  I_H = 631,   # Initial number of infected humans
  R_H = 13081,    # Initial number of recovered humans
  V_H = 0, 
  C = c0, # cumlative cases
  S_V = 500000,  # Initial number of susceptible mosquitoes
  I_V = 443   # Initial number of infected mosquitoes
)
t <- seq(1, dim(Dera)[1]-52, by = 1)  # Simulate for 2020-2024 week time step
model_output <- ode(y = initial_conditions, 
            times = t,
            func = mal_model_vac,
            parms = parameters)

df <- as.data.frame(model_output)
df <- as.data.frame(model_output)
inc <- c(c0,diff(df$C))
View(df)




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
  main = "Malaria in Dera 2020-2024(50% vaccination rate)" # title
)


# To be run with the different scenarios
df_50 <- as.data.frame(model_output)
df_50$inc <- c(c0,diff(df$C))

df_4 <- as.data.frame(model_output)
df_4$inc <- c(c0,diff(df$C))

df_10 <- as.data.frame(model_output)
df_10$inc <- c(c0,diff(df$C))


all_series <- cbind(
  df_4$inc,
  df_10$inc,
  df_50$inc,
  Dera$total_cases[53:313]  # observed data
)

matplot(
  df_4$time,
  all_series,
  type = "l",
  lty = c(1, 1, 1, 1),  # last one dashed for observed
  lwd = c(2, 2, 2, 1.5),
  col = c("blue", "green", "red", "grey50"),
  xlab = "Time (weeks)",
  ylab = "Incidence",
  main = "Malaria in Dera 2020–2024 under Different Vaccination Scenarios"
)

# Add legend
legend(
  "topleft",
  legend = c("4% vaccination", "10% vaccination", "50% vaccination", "Observed cases"),
  col = c("blue", "green", "red", "grey40"),
  lty = c(1, 1, 1, 1),
  lwd = 2
)

# Then for forecasting

forecast_period = 104 # Define timespan (weeks)
t_forc <- seq(min(t), max(t) + forecast_period, by = 1)



