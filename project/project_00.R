library(tidyverse)
rm(list = ls())
fd <- read_csv("project/gun_deaths full.csv")

# ---------------------#
# LOAD and CLEAN NAMES #
# -------------------- #
lower_names <- str_to_lower(names(fd))
wd <- fd |> `colnames<-`(lower_names)
names(wd)[str_detect(names(wd), "percent.+")] <- 
  str_replace(names(wd)[str_detect(names(wd), "percent.+")],"percent[s]{0,1}","pct")
names(wd) <- str_replace_all(names(wd),"\\.","_")
names(wd)[names(wd) == "paid_hunting_license_holders"]<-"hunt_lic"
names(wd)[names(wd) == "unemployment_rate"] <- "unem_rate"
names(wd)[names(wd) == "background_check_totals"] <- "gun_permits"


wd <- wd |> dplyr::select(state, year, 
                          total_firearm_deaths, population, 
                          hunt_lic, 
                          gun_permits, 
                          unem_rate) |>
      group_by(state, year) |>
      ungroup()

# even panel?
wd |> group_by(year) |> mutate(N = n()) %>% .$N |> summary()

# ~~~~~~~~~~~~#
# MISSINGNESS #
# ~~~~~~~~~~~~#
mi <- function(x) {
        apply(wd[,c("total_firearm_deaths","hunt_lic","unem_rate","gun_permits")],
        MARGIN = 2,
        FUN = function(x){ sum(is.na(x))} )
}
mi()
wd <- wd[complete.cases(wd),] #wd |> dplyr::filter(year > 1991 & year < 2015)
mi()


# ------------#
# TIME TRENDS #
# ------------#

# turn into ts object aggregated by year
gun_deaths <- wd |> 
              group_by(year) |>
              summarise(gun_deaths = sum(1e4 * total_firearm_deaths/population, 
                                         na.rm = TRUE)) %>% 
              .$gun_deaths |> ts()
xlabs <- unique(wd$year)
plot(gun_deaths, 
     ylab = "Gun deaths per 10,000", 
     xlab = "Year",
     xaxt = "n", 
     main = "Gun deaths\nin U.S per year")
axis(1, 
     at = seq(1, length(gun_deaths), by=3), 
     labels = xlabs[seq(1, length(gun_deaths), by=3)])


# gun purchases
gun_permits <- wd |> 
  group_by(year) |> 
  summarise(gun_permits = sum(1e4*gun_permits/population, 
                         na.rm = TRUE)) %>% 
  .$gun_permits |> ts()
xlabs <- unique(wd$year)
plot(gun_permits, 
     ylab = "Gun permits per 10,000", 
     xlab = "Year",
     xaxt = "n", 
     main = "Total number of guns purchased\nper year")
axis(1, at = seq(1, length(gun_permits), by=3), 
     labels = xlabs[seq(1, length(gun_permits), by=3)])


# hunting licenses
hunt_lic <- wd |> 
  group_by(year) |> 
  summarise(hunt_lic = sum(1e4*hunt_lic/population, 
                         na.rm = TRUE)) %>% 
  .$hunt_lic |> ts()

xlabs <- unique(wd$year)
plot(hunt_lic, 
     ylab = "Hunting licenses", 
     xlab = "Year",
     xaxt = "n", 
     main = "Total number hunting licenses issued\nper year")
axis(1, at = seq(1, length(hunt_lic), by=3), 
     labels = xlabs[seq(1, length(hunt_lic), by=3)])


# unem rate
unem_rate <- wd |> 
  group_by(year) |> 
  summarise(unem_rate = mean(unem_rate, 
                             na.rm = TRUE)) %>% 
  .$unem_rate |> ts()
xlabs <- unique(wd$year)
plot(unem_rate, 
     ylab = "Poverty rate", 
     xlab = "Year",
     xaxt = "n", 
     main = "U.S. unemployment rate from 1991 to 2016")
axis(1, at = seq(1, length(unem_rate), by=3), 
     labels = xlabs[seq(1, length(unem_rate), by=3)])

# all time plots
xlabs <- unique(wd$year)
plot(unique(wd$year), rethinking::standardize(gun_deaths), 
     type = "l", 
     col = "black", 
     xlab = "Year", 
     ylab ="standardized values", 
     ylim = c(-2.5, 2.5),
     main = "Standardized gun deaths, permits, and poverty rate \n in the U.S. from 1991 to 2017")
lines(unique(wd$year), 
      rethinking::standardize(gun_permits), type = "l", col = "blue")
lines(unique(wd$year), 
      rethinking::standardize(unem_rate), type = "l", col = "orange")
lines(unique(wd$year), 
      rethinking::standardize(hunt_lic), type = "l", col = "darkred")

legend("bottomleft", 
       c("gun deaths","gun purchases",
         "unem rate","hunt licenses"), 
       lty = 1, 
       col = c("black", "blue", "orange", "darkred"))


#----------#
# OUTLIERS #
#----------#
names(wd)[names(wd) == "total_firearm_deaths"] <- "gun_deaths"
wd <- wd |> mutate(log_gun_deaths = log(1e4 * gun_deaths/population))

show_outliers <- function(x, type = "gaus", tol = 1.5) {
  summ <- summary(x)
  q3 <- ifelse(type == "gaus", summ[5], summ[4] * -log(0.25))
  q1 <- ifelse(type == "gaus", summ[2], summ[4] * -log(0.75))
  iqr <- q3 - q1
  fence_upper <- q3 + (tol*iqr)
  fence_lower <- q1 - (tol*iqr)
  print(summ)
  outliers <- x[x<fence_lower | x>fence_upper]
  print(outliers)
  list("outliers" = outliers, "fu" = fence_upper, "fl" = fence_lower)
}

gun_deaths_os  <- show_outliers(wd$gun_deaths,  type = "exp", 3)
hunt_lic_os    <- show_outliers(wd$hunt_lic,    type = "exp", 3)
gun_permits_os <- show_outliers(wd$gun_permits, type = "exp", 3)

wd[wd$gun_permits >= gun_permits_os$fu, ]
wd[wd$gun_deaths >= gun_deaths_os$fu, ]

wd$gun_permits[wd$gun_permits > gun_permits_os$fu] <- gun_permits_os$fu


# -----------#
# INSTRUMENT #
# ---------- #

# correlation between hunting licenses and gun purchases
cor(wd$hunt_lic, wd$gun_permits)

cors <- vector(mode = "numeric", length(unique(wd$state)))
states <- vector(mode = "character", length(unique(wd$state)))
i = 0
for(x in unique(wd$state)){ 
  i = i + 1
  r <- cor(wd[wd$state == x,]$gun_permits,
           wd[wd$state == x,]$hunt_lic) |> round(2)
  cors[i] <- round(r, 2)
  states[i] <- x
  print(paste0(x, " =", r))
}
tibble(cors = cors, states = states) |> arrange(cors) |> View()


corsy <- vector(mode = "numeric", length(unique(wd$year)))
years <- vector(mode = "numeric", length(unique(wd$year)))
i = 0
for(x in unique(wd$year)){ 
  i = i + 1
  r <- cor(wd[wd$year == x,]$gun_permits,
           wd[wd$year == x,]$hunt_lic) |> round(2)
  corsy[i] <- round(r, 2)
  years[i] <- x
  print(paste0(x, " =", r))
}
tibble(corsy = corsy, years = years) |> arrange(corsy) |> View()
  
# -------------------- #
# OUTPUT ANALYSIS FILE #
# -------------------- #
an_data <- wd |> 
  ungroup() |> 
  dplyr::select(state, year, log_gun_deaths, gun_deaths, population,
                gun_permits, hunt_lic, unem_rate)
summary(an_data)
write_csv(an_data, "project/an_data.csv")



